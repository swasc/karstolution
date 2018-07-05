from . import constants, evaporation, cmodel_frac, O18EVA_MEAN, O18EVA
import numpy as np
import math

#main module for the ISOLUTION part of the Karstolution model, which deals with in-cave
#isotope fractionation. This is a translation of the matlab code from Deininger et al. (2012)
#preserving many of the comments. However, all components related to d13C from the original model
#have been removed (as only d18O is modelled in Karstolution)

try:
    from numba import jit
except ImportError:
    # do-nothing decorator
    jit = lambda x:x

@jit
def isotope_calcite(d, TC, pCO2, pCO2cave, h, V, phi, d18Oini, tt):
    """
    The isolution part of the model
    
    TODO:docstring
    
    Inputs
    ------
        - *d* 
            drip interval (s)
        - *TC*
            air temperature in degC
        - *pCO2*
            CO2 dissolved in drip water (ppmV -- ?)
        - *pCO2cave*
            CO2 in cave air (ppmV -- ?)
        - *h*
            humidity (percent, or 0--1 ?)
        - *V*
            ventilation (m/s ?)
        - *phi*
            mixing parameter for splash (1-no splash, 0-entire drop lost to splash)
        - *d18Oini*
            initial d18O of drip water (?)
        - *tt*
            timestep (months)

    
    Usage example:
    --------------
    isotope_calcite(drip_interval_ks1, cave_temp, drip_pco2, cave_pco2, h, v,
        phi,kststor118o,tt)
    
    """
    #boundary value constant parameters (equiv of BOUNDARY)
    R2smow = 0.00015575
    R18smow = 0.0020052
    R18vpdb = 0.0020672

    eva=evaporation.evaporation(TC, h, V)   #Evaporationrate (mol/l)
    e18_hco_caco, e18_hco_h2o, a18_m = cmodel_frac.cmodel_frac(TC)       #Fractionation factors
    TK = 273.15 + TC            #Absolute temperature (K)
    #Tau precipitation (s); t=d/a (according to Baker 98)
    Z = 0.0001/(1.188e-011 * TC**3 - 1.29e-011 * TC**2 + 7.875e-009 * TC + 4.844e-008)
    alpha = (1.188e-011 * TC**3 - 1.29e-011 * TC**2 + 7.875e-009 * TC + 4.844e-008)
    #Tau buffering, after Dreybrodt and Scholz (2010)
    T = 125715.87302 - 16243.30688*TC + 1005.61111*TC**2 - 32.71852*TC**3 + 0.43333*TC**4

    #Concentrations and mol mass
    #Concentrations of the spezies in the solution, with respect to soil pCO2
    outputsoil = constants.constants(TC, pCO2)
    #Concentrations of the spezies in the solution, with respect to cave pCO2
    outputcave = constants.constants(TC, pCO2cave)

    HCOSOIL = outputsoil[2][2]                   #HCO3- concentration, with respect to soil pCO2 (mol/l)
    HCOCAVE = outputcave[2][2]/np.sqrt(0.8)    #HCO3- concentration, with respect to cave pCO2 (mol/l)

    #Mol mass of the water, with respect to the volume of a single box (mol)
    h2o_ini = 0.1/18
    #Mol mass of the HCO3-(soil), with respect to the volume of a single box (mol)
    hco_ini = HCOSOIL*1e-4
    #Mol mass of the HCO3-(cave), with respect to the volume of a single box (mol)
    hco_eq = HCOCAVE*1e-4

    #Fractionation facors for oxygen isotope
    eps_m = a18_m - 1
    avl = (-7356./TK + 15.38)/1000. + 1
    abl = 1/(e18_hco_h2o + 1)
    a = 1/1.008
    f = 1/6.

    Rdrop18_w = ( (d18Oini / 1000.) + 1) * R18smow
    Rdrop18_b = Rdrop18_w / (e18_hco_h2o + 1)
    Rv18 = avl * Rdrop18_w

    # Definition of the input-parameter of the program "O18EVA" at the time
    # t=0s, e.g. at the beginning of the mix process

    hco_mix = hco_ini
    h2o_mix = h2o_ini
    HCOMIX = hco_mix / 1e-4

    r_hco18_mix = Rdrop18_b
    r_h2o18_mix = Rdrop18_w

    number = 0

    r18res = 0
    r18mix = 1

    # in the test case, this only takes one iteration
    while r18mix != r18res:

        number += 1

        r18_hco_res = r_hco18_mix

        temp = O18EVA.O18EVA(d, TC, pCO2, pCO2cave,  h, V, r_hco18_mix, r_h2o18_mix, Rv18, HCOMIX,
        h2o_mix,tt)

        hco_out = temp[3]                       #mol mass of hco
        h2o_out = temp[5]                       #mol mass of h2o

        r_hco18_out = temp[0]                   #oxygen isotopic ratio of hco
        r_h2o18_out = temp[1]                   #oxygen isotopic ratio of h2o

        #%%% 1) simple mixprocess
        hco_mix = phi*hco_ini + (1-phi)*hco_out                 #mixing of hco
        h2o_mix = phi*h2o_ini + (1-phi)*h2o_out                 #mixing of h2o

        phi_r_b = 1 / (1 + (1-phi) / phi*hco_out*hco_ini)       #isotopic mixing parameter of hco
        phi_r_w = 1 / (1 + (1-phi) / phi*h2o_out/h2o_ini)       #isotopic mixing parameter of h2o

        r_hco18_mix = phi_r_b*Rdrop18_b + (1-phi_r_b)*r_hco18_out        #new oxygen isotopic ratio of hco
        r_h2o18_mix = phi_r_w*Rdrop18_w + (1-phi_r_w)*r_h2o18_out        #new oxygen isotopic ratio of h2o

        H2O_mix = h2o_mix*18/1000.      #mol -> l (mol*(g/mol)/(g/kg)*(l/kg)
        HCOMIX = hco_mix / H2O_mix     #new hco condentration

        #end of the extended mixprocess

        r18mix = round(r_hco18_mix*10**13)
        r18res = round(r18_hco_res*10**13)
        
        # bail out of loop if we get invalid values (O18EVA can return NaN
        # if droplets evaporate completely)
        if np.isnan(r18mix) or np.isnan(r18res):
            break


    r_hco18, r_h2o18, hco, h2o, delta_1  = O18EVA_MEAN.O18EVA_MEAN(d,
                TC, pCO2, pCO2cave, h, V, r_hco18_mix, r_h2o18_mix, Rv18,
                HCOMIX, h2o_mix,tt)
    delta_0 = delta_1[0]
    delta_end= delta_1[-1]


    tmp = np.isnan(r_hco18)
    if tmp.all() != True:
        # use math.fsum to improve accuracy, compared to numpy.sum
        r_hco18_mean = math.fsum( (r_hco18[:-1] + np.diff(r_hco18,n=1)/2.) * (-np.diff(hco,n=1) / (hco[0]-hco[-1]) ))
        d18Ocalcite = (r_hco18_mean*(e18_hco_caco + 1)/R18vpdb - 1)*1000
        (r_hco18[[0,-1]]*(e18_hco_caco + 1)/R18vpdb-1)*1000

        r_h2o18_mean = math.fsum( (r_h2o18[:-1] + np.diff(r_h2o18,n=1)/2.) * (-np.diff(h2o,n=1) / (h2o[0]-h2o[-1]) ))
        d18Owater = (r_h2o18_mean/R18smow - 1)*1000

        d18Ovapor = (Rv18/R18smow - 1)*1000

        return d18Ocalcite

    else:
        return np.NaN

import math
from . import evaporation, cmodel_frac, constants
try:
    from numba import jit
except ImportError:
    # do-nothing decorator
    jit = lambda x:x
    
import numpy as np

@jit
def O18EVA(tmax, TC, pCO2, pCO2cave, h, v, R18_hco_ini, R18_h2o_ini, R18v, HCOMIX, h2o_new,tt):

    # Sourcecode to develope the evolution of the isotopic ratio of the oxygen
    # compostion of the oxygen isotopes 16O and 18O as a function of
    # temperature TC, supersaturation (pCO2), relative humidity (h) and wind
    # velocity (v). %(08.12.2010/m)

    eva = evaporation.evaporation(TC, h, v)
    e18_hco_caco, e18_hco_h2o, a18_m = cmodel_frac.cmodel_frac(TC)
    TK = 273.15 + TC
    #Tau precipitation (s); t=d/a (according to Baker 98)
    alpha_p = (1.188e-011 * TC**3 - 1.29e-011 * TC**2 + 7.875e-009 * TC + 4.844e-008)
    #Tau buffering, after Dreybrodt and Scholz (2010)
    T = 125715.87302 - 16243.30688*TC + 1005.61111*TC**2 - 32.71852*TC**3 + 0.43333*TC**4

    #Concentrations and mol mass
    outputcave = constants.constants(TC, pCO2cave)      #Concentrations of the spezies in the solution, with respect to cave pCO2
    HCOCAVE = outputcave[2][2]/np.sqrt(0.8)    #HCO3- concentration, with respect to cave pCO2 (mol/l)

    h2o_ini = h2o_new                          #Mol mass of the water, with respect to the volume of a single box (mol)
    H20_ini = h2o_ini*18/1000.
    hco_ini = HCOMIX*H20_ini                   #Mol mass of the HCO3-(soil), with respect to the volume of a single box (mol)
    #hco_eq = HCOCAVE*H20_ini                   #Mol mass of the HCO3-(cave), with respect to the volume of a single box (mol)

    #Fractionation facors for oxygen isotope
    eps_m = a18_m - 1
    avl = ((-7356./TK + 15.38)/1000. + 1)
    abl = 1/(e18_hco_h2o + 1)
    a = 1/1.008*1.003
    f = 1/6

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #Calculation of the 18R

    N_times = int(np.ceil(tmax + 1))
    
    #initialise arrays
    init = np.empty(N_times)
    init[:] = np.NaN
    r_hco18 = init.copy()
    r_h2o18 = init.copy()
    HCO = init.copy()
    hco = init.copy()
    H2O = init.copy()
    h2o = init.copy()

    r_hco18[0] = R18_hco_ini
    r_h2o18[0] = R18_h2o_ini


    if eva > 0 and tmax > np.floor(h2o_ini/eva):
        tmax = int(np.floor(h2o_ini/eva))
        #raise RuntimeError('DRIPINTERVALL IS TOO LONG, THE WATERLAYER EVAPORATES COMPLETLY FOR THE GIVEN d (tt={})'.format(tt))
        return (np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN)

    # adjust dt so that it's roughly 1 second, but divides evenly into tmax
    t = np.linspace(0, tmax, N_times)
    dt = t[1] - t[0]

    HCO[0] = HCOMIX                                  #Konzentration von HCO3-
    hco[0] = hco_ini                                 #Menge an HCO3-
    H2O[0] = H20_ini
    h2o[0] = h2o_new

    #"Restwassermenge" und daher konstant
    # "Residual water quantity" and therefore constant

    HCO_EQ = HCOCAVE                       #Equilibriumconcentration

    for ii in range(1, len(t)):

        delta = (H2O[ii-1]/1000.)/0.001

        # Update the state of aqueous solution to present timestep
        h2o[ii] = (h2o[ii-1] - eva*dt)              #Water (mol)
        d_h2o = -eva                              #Evaporationrate (mol/l)
        H2O[ii] = (h2o[ii]*18*1e-3)               #Water (l)

        #Evaporation
        HCO_temp = (HCO[ii-1] - HCO_EQ) * np.exp(-dt/(delta/alpha_p)) + HCO_EQ          #HCO3- concentration after timeintervall dt
        hco[ii] = (HCO_temp * H2O[ii-1])                                         #HCO3- mass (mol)

        HCO[ii] = (HCO_temp * (H2O[ii-1]/H2O[ii]))      #HCO3- concentration after timeintervall dt and the evaporation of water

        r_hco18[ii] = (r_hco18[ii-1] + ((eps_m*(hco[ii]-hco[ii-1])/hco[ii]-1/T) * r_hco18[ii-1] + abl/T*r_h2o18[ii-1]) * dt)
        if eva > 0:
            r_h2o18[ii] = (r_h2o18[ii-1] + 
                            ( ( hco[ii]/h2o[ii]/T 
                                - f/abl/h2o[ii]*(hco[ii]-hco[ii-1]) * r_hco18[ii] 
                                + ( d_h2o/h2o[ii]*(a*avl/(1-h)-1) - hco[ii]/h2o[ii]*abl/T) * r_h2o18[ii-1] 
                                - a*h/(1-h)*R18v/h2o[ii]*d_h2o
                              ) * dt
                            )
                           )
        else:
            # in the limit as eva -> 0, h=1, and d_h2o=0 
            # - so the commented term is zero
            r_h2o18[ii] = (r_h2o18[ii-1] + 
                            ( ( hco[ii]/h2o[ii]/T 
                                - f/abl/h2o[ii]*(hco[ii]-hco[ii-1]) * r_hco18[ii] 
                                + ( d_h2o/h2o[ii]*(a*avl/(1-h)-1) - hco[ii]/h2o[ii]*abl/T) * r_h2o18[ii-1] 
                                #- a*h/(1-h)*R18v/h2o[ii]*d_h2o
                              ) * dt
                            )
                           )

    #assert np.isnan( np.array([r_hco18[-1], r_h2o18[-1], HCO[-1], hco[-1], H2O[-1], h2o[-1]])  ).sum() == 0

    return (r_hco18[-1], r_h2o18[-1], HCO[-1], hco[-1], H2O[-1], h2o[-1])

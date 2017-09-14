from scipy import stats
import numpy as np
from isotope_calcite import isotope_calcite

def weibull_parameters_y(w,z, __cache=[None,None]):
    """
    Calculate y, derived from the Weibull distribution

    This is slow - so we remember the result for next time, as this
    tends to get called lots of times with the same values for w and z
    """
    if __cache[0] == (w,z):
        return __cache[1]
    else:
        x=np.linspace(0,2,12)
        v_1=stats.exponweib(w,z)
        y1=v_1.cdf(x)
        y=np.append([0],y1[1:]-y1[0:11])
        __cache[0] = (w,z)
        __cache[1] = y
    return y

def karst_process(tt,mm,evpt,prp,prpxp,tempp,d18o,d18oxp,dpdf,epdf,soilstorxp,soil18oxp,
epxstorxp,epx18oxp,kststor1xp,kststor118oxp,kststor2xp,kststor218oxp,config,calculate_drip,cave_temp,
calculate_isotope_calcite=True):
    """
    this function contains all the karst hydrological processes (based on the KarstFor code)
    followed by execution of the ISOLTUION code for in-cave processes (isotope_calcite module)
    """

    mf = config['monthly_forcing']
    #store size parameters  - soilstore, epikarst, ks1, ks2
    soilsize=config['soilstore']
    episize=config['epikarst']
    ks1size=config['ks1']
    ks2size=config['ks2']

    #making sure the init sizes don't exceed store capacity
    if soilstorxp>soilsize:
        soilstorxp=soilsize-1
    if epxstorxp>episize:
        epxstorxp=episize-1
    if kststor1xp>ks1size:
        kststor1xp=ks1size-1
    if kststor2xp>ks2size:
        kststor2xp=ks2size-1

    #overflow parameters
    epicap=config['epicap']
    ovcap=config['ovicap']
    #ensuring the overflow parameters are less than the store
    if epicap>=episize:
        epicap=episize-1
    if ovcap >= ks2size:
        ovcap=ks2size-1

    #average cave parameters for various months
    drip_interval = mf['drip_interval'][mm-1]
    drip_pco2=mf['drip_pco2'][mm-1]/1000000.0
    cave_pco2 = mf['cave_pco2'][mm-1]/1000000.0
    h = mf['rel_humidity'][mm-1]
    v = mf['ventilation'][mm-1]
    phi = config['mixing_parameter_phi']

    #making sure cave values don't become negative
    if v<0:
        v=0
    if drip_interval<0:
        drip_interval=0
    if drip_pco2<0:
        drip_pco2=0.0000000000000001
    if cave_pco2<0:
        cave_pco2=0.0000000000000001
    if h<0:
        h=0
    if phi<0:
        phi=0

    #making sure some cave values don't exceed one
    if h>=1:
        h =0.99
    if phi>1:
        phi=1


    #weibull parameters
    w=config['lambda_weibull'] #data_rest[6][0]
    z=config['k_weibull'] #data_rest[6][1]
    y = weibull_parameters_y(w,z)


    #parameterisable coefficients
    k_f1=config['f1'] #data_rest[0][0]        #f1 from soilstore to epikarst
    k_f3=config['f3'] #data_rest[0][1]        #f3 from epikarst to KS1
    k_f8=config['f8'] #data_rest[0][6]
    k_f5=config['f5'] #data_rest[0][2]        #f5 from KS1 to stal5
    k_f6=config['f6'] #data_rest[0][3]        #f6 from KS2 to stal1
    k_f7=config['f7'] #data_rest[0][4]        #f7 overflow from KS2 to KS1
    k_diffuse=config['k_diffuse'] #data_rest[0][5]   #diffuse flow from Epikarst to KS1
    k_e_evap=config['k_eevap'] #data_rest[2][0]    #epikarst evap (funct of ET for timestep) Used for both sources???
    k_evapf=config['k_d18o_soil'] #kdata_rest[2][1]     #soil evap d18o fractionation from somepaper????
    k_e_evapf=config['k_d18o_epi'] #data_rest[2][2]   #epikarst evap d18o fractionation ??? can use same value?
    i=config['i'] #data_rest[1][0]           #epikarst in bypass flow mixture to stal1, (<1 & i+j+k=1)
    j=config['j'] #data_rest[1][1]           #rain in bypass flow mixture to stal1, (<1 & i+j+k=1)
    k=config['k'] #data_rest[1][2]           #rain from last step in bypass flow mixture to stal1, (<1 & i+j+k=1)
    m=config['m'] #data_rest[1][3]           #epikarst in bypass flow mixture to stal2, (<1 & m+n=1)
    n=config['n'] #data_rest[1][4]           #rain in bypass flow mixture to stal2, (<1 & m+n=1)

    # these are placeholders if we're not running the ISOLTUION part of the model
    if not calculate_isotope_calcite:
        stal1d18o,stal2d18o,stal3d18o,stal4d18o,stal5d18o = (np.NaN, np.NaN,
                                                        np.NaN, np.NaN, np.NaN)

    #********************************************************************************************
    #starting going through the karst processes in a procedural manner (up-down)
    #making sure the soilstore does not become negative, whilst adding prp and removing evpt
    if soilstorxp + prp - evpt < 0:
        #evpt=0
        #^^^partially agreed can remove the above
        soilstor=0
    # if prp>=7:
        # soilstor=soilstorxp+prp-evpt
    else:
        soilstor=soilstorxp+prp-evpt

    #ensuring the soilstor does not exceed user-defined capacity
    if soilstor>soilsize:
        soilstor=soilsize

    #prevents any flux when surface is near-frozen. in this case, 0.0 degree c
    if tempp[0]>0.0:
        f1=soilstor*k_f1
    else:
        f1=0
    #updating the final soil store level (removing the F1 value)
    soilstor=soilstor-f1

    #increases epikarst store volume
    epxstor=epxstorxp+f1
    #draining from bottom first as gravity fed
    f3 = epxstor*k_f3
    #diffuse flow leaving epikarst and going to KS1
    #assuming diffuse flow follows a weibull distrubtion
    dpdf[0]=(epxstor-f3)*k_diffuse
    if epxstor-f3-dpdf[0]> epicap:
        f4=(epxstor-f3-dpdf[0]-epicap)
    else:
        f4=0

    #epikarst evaporation starts when soilstore is 10%
    #and increases with decreasing soil store
    #added the (1-4*...) term, can change the '4'
    if prp==0:
        e_evpt=k_e_evap*evpt
    elif soilstor<=0.1*soilsize:
        e_evpt=k_e_evap*evpt*(1-4*soilstor/soilsize)
    #elif condition is a second route for epikarst evaporation through
    #bypass flow which is the same route as used for stal 2 & 3
    else:
        e_evpt=0

    #calculating final epikarst value
    if epxstor-f3-f4-dpdf[0]-e_evpt<0:
        epxstor=0
    else:
        epxstor=epxstor-f3-f4-dpdf[0]-e_evpt

    #ensuring the epxstor does not exceed user-defined capacity
    if epxstor>episize:
        epxstor=episize

    #fluxes into and out of KS2
    kststor2=kststor2xp+f4
    if kststor2 > ovcap:
        f7=(kststor2-ovcap)*k_f7
    else:
        f7=0
    f6=(kststor2-f7)*k_f6
    kststor2=kststor2-f6-f7

    #ensuring the kststor2 does not exceed user-defined capacity
    if kststor2>ks2size:
        kststor2=ks2size

    #f8 bypass flow from surface rain to KS1
    # TODO: fix this magic number (rainfall threshold of 7)
    if prp>7:
        f8=prp*k_f8
    else:
        f8=0

    #fluxes into and out of KS1
    kststor1=kststor1xp+f3+sum(y*dpdf)+f7+f8
    f5=kststor1*k_f5
    kststor1=kststor1-f5

    #ensuring the kststor1 does not exceed user-defined capacity
    if kststor1>ks1size:
        kststor1=ks1size

    #mixing and fractionation of soil store d18o
    e=prp+soilstorxp
    if e<0.01:
        e=0.001
    f=soilstorxp/e
    g=prp/e
    # 0.03 term can be changed to enable evaporative fractionation in soil store
    h_1=d18o+(evpt*k_evapf)
    #mixing of soil d18o with prp and ???
    soil18o=(f*soil18oxp)+(g*h_1)

    #so if the soil value becomes positive it is reverted to original soild18o. Justified??
    if soil18o>0.0001:
        soil18o=soil18oxp

    #mixing and fractionation of epikarst store d18o
    a=f1
    b=a+epxstorxp
    #quick fix for divide-by-zero error when b is too small
    if b<=0.001:
        b=0.001
    c=(epxstorxp/b)*(epx18oxp+e_evpt*k_e_evapf)
    d=(a/b)*soil18o
    epx18o=c+d
    epdf[0]=epx18o

    #mixing of kststor2 d18o
    if f4<0.01:
        kststor218o=kststor218oxp
    else:
        b2=f4+kststor2xp
        c2=(kststor2xp/b2)*kststor218oxp
        d2=(f4/b2)*epx18o
        kststor218o=c2+d2

    #mixing of KS1 d18o
    b1=f3+kststor1xp+sum(y*dpdf)+f7+f8
    c1=(kststor1xp/b1)*kststor118oxp
    d1=(f3/b1)*epx18o
    e1=(sum(y*dpdf*epdf)/b1)
    g1=(f7/b1)*kststor218o
    h1=f8/b1*d18o
    kststor118o=c1+d1+e1+g1+h1

    #bypass flow (from epikarst and direct from rain)
    p=d18o
    r=d18oxp
    drip118o=(kststor118o*i)+(p*j)+(r*k)
    drip218o=(kststor118o*m)+(p*n)


    #if kststor2 is too low then a clearly outlying for the stal
    if kststor2<0.01:
        stal1d18o=-99.9
        drip_interval_ks2=9001
    else:
        if calculate_drip==True:
            #drip-interval: user inputted max drip-interval proportioned by store capacity
            drip_interval_ks2=int(drip_interval*ks2size/kststor2)
        else:
            drip_interval_ks2=int(drip_interval)
        if calculate_isotope_calcite:
            #running the ISOLUTION part of the model
            stal1d18o=isotope_calcite(drip_interval_ks2, cave_temp, drip_pco2, cave_pco2, h, v, phi,
            kststor218o,tt)


    #same drip interal calculation for epikarst store (stalagmite 4)
    if epxstor<0.01:
        stal4d18o=-99.9
        drip_interval_epi=9001
    else:
        if calculate_drip==True:
            drip_interval_epi=int(drip_interval*episize/epxstor)
        else:
            drip_interval_epi=int(drip_interval)
        if calculate_isotope_calcite:
            stal4d18o=isotope_calcite(drip_interval_epi, cave_temp, drip_pco2, cave_pco2, h, v, phi,
            epx18o,tt)

    #drip interval calculations for Karst Store 1, which includes the bypass stalagmites 2 and 3.
    #Drip interval for these are
    if kststor1<0.01:
        stal2d18o=-99.9
        stal3d18o=-99.9
        stal5d18o=-99.99
        drip_interval_ks1=9001
        drip_interval_stal3=9001
        drip_interval_stal2=9001
    else:
        if calculate_drip==True:
            drip_interval_ks1=int(drip_interval*ks1size/kststor1)
            ks1_temp3=kststor1+prp
            drip_interval_stal3=int(drip_interval*ks1size/ks1_temp3)
            ks1_temp2=ks1_temp3+prpxp
            drip_interval_stal2=int(drip_interval*ks1size/ks1_temp2)
        else:
            drip_interval_ks1=int(drip_interval)
            drip_interval_stal3=int(drip_interval)
            drip_interval_stal2=int(drip_interval)
        if calculate_isotope_calcite:
            stal5d18o=isotope_calcite(drip_interval_ks1, cave_temp, drip_pco2, cave_pco2, h, v,
            phi,kststor118o,tt)
            stal3d18o=isotope_calcite(drip_interval_stal3, cave_temp, drip_pco2, cave_pco2, h, v,
            phi,drip218o,tt)
            stal2d18o=isotope_calcite(drip_interval_stal2, cave_temp, drip_pco2, cave_pco2, h, v,
            phi,drip118o,tt)

    #returning the values to karstolution1.1 module to  be written to output
    return [tt,mm,f1,f3,f4,f5,f6,f7,soilstor,epxstor,kststor1,kststor2,soil18o,epx18o,kststor118o,
    kststor218o,dpdf[0],stal1d18o,stal2d18o,stal3d18o,stal4d18o,stal5d18o,drip_interval_ks2,
    drip_interval_epi,drip_interval_stal3,drip_interval_stal2,drip_interval_ks1,cave_temp]

#function which ensures a variable is not negative and if it is, assigna a small positive value
def checkzero(variable):
    if variable<=0:
        variable=0.001
    return variable

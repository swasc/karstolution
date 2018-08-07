from __future__ import division
import numpy as np
import pandas as pd
from . import karst_process

#this function unpacks and initialses some of the model parameters, reads the input file
#iterates each step (according to each entry of the input file)
#and writes the output into the defined file at each step

def karstolution(config,df_input,calculate_drip=True, calculate_isotope_calcite=True):
    # abbreviations...
    c = config
    ic = c['initial_conditions']
    mf = c['monthly_forcing']
    # number of months of history in weibull distribution, default of 12
    weibull_delay_months = config.get('weibull_delay_months', 12)
    weibull_delay_months = int(weibull_delay_months)

    # headers for the output
    output_columns = ['tt','mm','f1','f3','f4','f5','f6','f7','soilstor','epxstor',
    'kststor1','kststor2','soil18o','epx18o','kststor118o','kststor218o','dpdf[0]',
    'stal1d18o','stal2d18o','stal3d18o','stal4d18o','stal5d18o','drip_int_stal1',
    'drip_int_stal4','drip_int_stal3','drip_int_stal2','drip_int_stal5','cave_temp']

    #unpacking inital values taken from the configuration file (in list 'data')
    # soili, epikarsti, ks1i, ks2i, diffusei

    epx18oxp=ic['d18o_epikarst']         #inital d18O in epikarst
    epxstorxp=ic['epikarst']      #inital water level of epikarst
    soilstorxp=ic['soil']      #inital water level soil store
    soil18oxp=ic['d18o_soil']       #inital d18O in soil store
    kststor1xp=ic['ks1']      #inital water level of ks1 store
    kststor2xp=ic['ks2']      #inital water level of ks2 store
    kststor118oxp=ic['d18o_ks1']   #inital d18O in soil store
    kststor218oxp=ic['d18o_ks2']   #inital d18O in ks1 store
    d18oxp=ic['d18o_prevrain']
    prpxp=0                     #initial value for rain from 'previous' step
    #weibull distribution (diffuse flow) currently set to 12 months, can be increased
    dpdf=[ic['diffuse']]*weibull_delay_months       #inital water quantity in the weibull distribution (diffuse flow)
    epdf=[ic['d18o_diffuse']]*weibull_delay_months       #inital d18O in the weibull distribution (diffuse flow)


    #the rest of the data which we are not unpacking (but will do so later)
    #data_rest=[thing for thing in data[0:12]]


    #36 month surface temp list for purposes of coupling surface to cave; currently a dummy list
    tempp=[thing for thing in range(0,36)]
    #finding the average inputted cave value
    avr_cave=np.mean(mf['cave_temp'])
    #setting a dummy value that will be overwritten below
    difference=10

    output_rows = []
    num_timesteps = len(df_input)

    #reading the input file and using each row as one iteration of the model
    for index, row in df_input.iterrows():
        #using the headings of the columns
        tt=int(row['tt']) #step number (starts at 1 and ends at # of iterations)
        mm=int(row['mm']) #month: varies from 1-12
        evpt=float(row['evpt']) #value of evapotranspiration (mm)
        prp=float(row['prp'])  #value of precipitation (mm)
        #update the first value in the list with new input value
        tempp[0]=float(row['tempp']) #surface temperature (celcius)
        d18o=float(row['d18o']) #d18O value of rainfall (per mille)

        #for the first loop of the program filling tempp with the first input value
        if tt==1:
            tempp=[tempp[0] for thing in range(tt-1,36)]
            #the difference between the surface temp and cave temp
            #this value determines the set diff for rest of model
            difference=tempp[0]-avr_cave
        #seasonlity factor for that month based on GUI cave temp inputs
        seasonality = mf['cave_temp'][mm-1] - avr_cave

        #average surface temp of last 36 months
        avr_surfacet=sum(tempp)/len(tempp)
        #cave temp is the average surface temp - the set surface-cave temp difference from step 1
        #with an adjustment for seasonality
        cave_temp=avr_surfacet-difference+seasonality

        #passes each of the paramters through to the karst_process function...
        out=karst_process.karst_process(tt,mm,evpt,prp,prpxp,tempp,d18o,d18oxp,dpdf,epdf,soilstorxp,
        soil18oxp,epxstorxp,epx18oxp,kststor1xp,kststor118oxp,kststor2xp,kststor218oxp,config,
        calculate_drip, cave_temp,
        calculate_isotope_calcite=calculate_isotope_calcite)

        output_rows.append(out)

        # Note: this is what out contains:
        # 0  tt,mm,f1,f3,f4,
        # 5  f5,f6,f7,soilstor,epxstor,
        # 10 kststor1,kststor2,soil18o,epx18o,kststor118o,
        # 15 kststor218o,dpdf[0],stal1d18o,stal2d18o,stal3d18o,
        # 20 stal4d18o,stal5d18o,drip_interval_ks2, drip_interval_epi,drip_interval_stal3,
        # 25 drip_interval_stal2,drip_interval_ks1,cave_temp


        #update model terms for next iteration
        epx18oxp=out[13]
        epxstorxp=out[9]
        soilstorxp=out[8]
        soil18oxp=out[12]
        kststor1xp=out[10]
        kststor2xp=out[11]
        kststor118oxp=out[14]
        kststor218oxp=out[15]
        epdf[1:]=epdf[:-1]
        epdf[0]=out[13]
        dpdf[1:]=dpdf[:-1]
        dpdf[0]=out[9]
        d18oxp=d18o
        prpxp=prp
        tempp[1:36]=tempp[0:35]

    output_dataframe = pd.DataFrame.from_records(output_rows, columns=output_columns)
    return output_dataframe

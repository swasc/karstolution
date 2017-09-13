# Karstolution
The First Speleothem δ18O Model Integrating Karst Hydrological and In-Cave Fractionation Processes.  
Coupling of existing KarstFor (Bradley et al., 2010) and ISOLTUION (Deininger et al., 2012) models.   
A manuscript presenting Karstolution and with a case-study of Golgotha Cave is under review in GCA.  
A Windows executable GUI has also been created to allow easy use of the model. This is downloadable from www.connectedwaters.edu.au/karstolution. The wxpython code is also available in this GitHub despoitory. See below.  
This python code is to allow furhter customisation of the Karstolution model, improvements and allow more advanced model runs such as customisable plotting and batch runs.  

# Conceptual Figure
![Alt text](https://cloud.githubusercontent.com/assets/19680492/15954071/f8490ce4-2f15-11e6-822b-d1087f8248a9.png "Karstolution Conceptual Figure")

# Dependancies
Numpy  
Scipy  
(optional) matplotlib  
(optional) pandas  

# Installation

TODO

# Configuration File

The configuration is passed to the main model routine as a `dict`, i.e. a Python dictionary.
A convenient way of soring the configuration is as a [yaml](http://yaml.org/) formatted file.
See `example/example.py` for full details.

```yaml
f1 : 0.2
f3 : 0.008
f5 : 0.005
f6 : 0.002
f7 : 1.0
k_diffuse : 0.008
f8 : 0.001
i : 0.5
j : 0.25
k : 0.25
m : 0.75
n : 0.25
k_eevap : 0.0
k_d18o_soil : 0.03
k_d18o_epi : 0.0
soilstore : 200.0
epicap : 400.0
ovicap : 100.0
epikarst : 400.0
ks1 : 400.0
ks2 : 200.0
lambda_weibull : 1.5
k_weibull : 1.0
mixing_parameter_phi : 1.0
# these parameters are forced by a climatological monthly mean
# (so there needs to be a list of 12 values, January-December)
monthly_forcing : 
  drip_pco2 : [4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0,4000.0]
  cave_pco2 : [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0]
  rel_humidity : [0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95]
  ventilation : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  cave_temp : [10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0]
  drip_interval : [100.0,100.0,100.0,100.0,100.0,100.0,100.0,100.0,100.0,100.0,100.0,100.0]
initial_conditions :
  # initial level in each store
  soil : 50.0
  epikarst : 100.0
  ks1 : 230.0
  ks2 : 50.0
  diffuse : 30.0
  # initial oxygen-18 composition in each store (d18O, permille VSMOW)
  d18o_soil : -5.0
  d18o_epikarst : -4.0
  d18o_ks1 : -5.0
  d18o_ks2 : -4.0
  d18o_prevrain : -5.0
  d18o_diffuse : -4.0
```  
F1,F3,F5,F6,F7,k_diffuse,f8  
i,j,k,m,n  
k_eevap, k_d18o_soil, k_d18o_epi  
Cave temps (monthly: Jan-Dec)  
Drip intervals (monthly: Jan-Dec)  
Soilstore, epicap, ovicap, epikarst, KS1, KS2 (store sizes)  
Lambda, k (weibull)  
Drip pCO2 (monthly: Jan-Dec)  
Cave pCO2 (monthly: Jan-Dec)  
Rel Humidity (monthly: Jan-Dec)  
Ventilation (monthly: Jan-Dec)  
Mixing parameter (phi)  
Store initial values (mm): soil, epikarst, ks1, ks2, diffuse  
δ18O intial values (per mille):  soil, epikarst, ks1, ks2, prev rain, diffuse  

# Input file
The input file is a csv of climatic inputs, a similar format to that of KarstFor (example is provided).  
Note: the model steps are in months and the number of rows represents the number of model steps   
The columns are:  
tt: an id column with numbers from 1 to total number of model steps  
mm: representing the month of that model step (important for seasonality), values 1-12  
evpt: evapotranspiration (mm)  
prp: rainfall amount (mm)  
tempp: surface temperature (degree celcius)  
d18O: the δ18O of rainfall amount  


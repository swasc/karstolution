import wx
import csv
from functools import partial
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt

class Advanced(wx.Panel):
    def __init__(self,parent,data,name_config):
        wx.Panel.__init__(self,parent=parent)
        self.data1=data
        self.config=name_config
        self.font = wx.Font(14, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        posx=350
        #fluxes
        self.flabel=wx.StaticText(self,-1,'Fluxes',(30,posx),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.flabel=wx.StaticText(self,-1,'(1/month)',(110,posx+10),(200,-1))
        self.flabel=wx.StaticText(self,-1,'Bypass Composition',(190,posx),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        
        
        self.flabel=wx.StaticText(self,-1,'Evaporation Coeffs',(440,posx+30),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        
        posx+=35
        self.F1=wx.StaticText(self,-1,'F1: ',(10,posx))
        self.F1_v=wx.TextCtrl(self,-1,str(self.data1[0][0]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F1_v,id1=0,id2=0),self.F1_v)
        self.i=wx.StaticText(self,-1,'i (%, stal 1, epikarst): ',(190,posx))
        self.k_i=wx.TextCtrl(self,-1,str(self.data1[1][0]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_i,id1=1,id2=0),self.k_i)
        
        posx+=25
        self.F3=wx.StaticText(self,-1,'F3: ',(10,posx))
        self.F3_v=wx.TextCtrl(self,-1,str(self.data1[0][1]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F3_v,id1=0,id2=1),self.F3_v)
        self.j=wx.StaticText(self,-1,'j (%, stal 1, rain): ',(190,posx))
        self.k_j=wx.TextCtrl(self,-1,str(self.data1[1][1]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_j,id1=1,id2=1),self.k_j)
        epitext="k_eevap: epikarst remainder that"
        self.epit=wx.StaticText(self,-1,epitext,(440,posx))
        posx+=25
        self.F8=wx.StaticText(self,-1,'F8: ',(10,posx))
        self.F8_v=wx.TextCtrl(self,-1,str(self.data1[0][6]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F8_v,id1=0,id2=6),self.F8_v)
        self.k=wx.StaticText(self,-1,'k (%, stal 1, prev rain): ',(190,posx))
        self.k_k=wx.TextCtrl(self,-1,str(self.data1[1][2]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_k,id1=1,id2=2),self.k_k)
        self.e_epievap=wx.StaticText(self,-1,"evaporates (1/month): ",(440,posx))
        self.k_e_evap=wx.TextCtrl(self,-1,str(self.data1[2][0]),(590,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_e_evap,id1=2,id2=0),self.k_e_evap)
        posx+=25
        self.F5=wx.StaticText(self,-1,'F5: ',(10,posx))
        self.F5_v=wx.TextCtrl(self,-1,str(self.data1[0][2]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F5_v,id1=0,id2=2),self.F5_v)
        self.m=wx.StaticText(self,-1,'m (%, stal 2, epikarst): ',(190,posx))
        self.k_m=wx.TextCtrl(self,-1,str(self.data1[1][3]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_m,id1=1,id2=3),self.k_m)
        epitext2="d18O evap fract (per-mille*month/mm):"
        self.epit2=wx.StaticText(self,-1,epitext2,(440,posx))
        posx+=25
        self.F6=wx.StaticText(self,-1,'F6: ',(10,posx))
        self.F6_v=wx.TextCtrl(self,-1,str(self.data1[0][3]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F6_v,id1=0,id2=3),self.F6_v)
        self.n=wx.StaticText(self,-1,'n (%, stal 2, rain): ',(190,posx))
        self.k_n=wx.TextCtrl(self,-1,str(self.data1[1][4]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_n,id1=1,id2=4),self.k_n)
        self.e_evap=wx.StaticText(self,-1,"soil (k_d18o_soil): ",(440,posx))
        self.k_evapf=wx.TextCtrl(self,-1,str(self.data1[2][1]),(590,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_evapf,id1=2,id2=1),self.k_evapf)
        posx+=25
        self.F7=wx.StaticText(self,-1,'F7: ',(10,posx))
        self.F7_v=wx.TextCtrl(self,-1,str(self.data1[0][4]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F7_v,id1=0,id2=4),self.F7_v)
        self.n=wx.StaticText(self,-1,'*Note*: i + j + 1 = 1 and m + n = 1',(190,posx))
        self.e_evapf=wx.StaticText(self,-1,"epikarst (k_d18o_epi): ",(440,posx))
        self.k_e_evapf=wx.TextCtrl(self,-1,str(self.data1[2][2]),(590,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.k_e_evapf,id1=2,id2=2),self.k_e_evapf)
        pos=posx+25
        self.F8=wx.StaticText(self,-1,'k_Diffuse: ',(10,pos))
        self.F8_v=wx.TextCtrl(self,-1,str(self.data1[0][5]),(70,pos))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F8_v,id1=0,id2=5),self.F8_v)
        
        #cave temperature, drip int and humidity
        posx+=45
        self.label=wx.StaticText(self,-1,'Cave Temperature',(10,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(Celcius)',(230,posx+10))
        self.label=wx.StaticText(self,-1,'Min Drip Interval',(300,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(Seconds)',(505,posx+10))
        self.label=wx.StaticText(self,-1,'Rel Humidity',(600,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(%; 0<h<=1)',(760,posx+10))
        
        posx+=40
        self.Jan=wx.StaticText(self,-1,'Jan: ',(10,posx))
        self.Jan_temp=wx.TextCtrl(self,-1,str(self.data1[3][0]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jan_temp,id1=3,id2=0),self.Jan_temp)
        self.Jul=wx.StaticText(self,-1,'Jul: ',(150,posx))
        self.Jul_temp=wx.TextCtrl(self,-1,str(self.data1[3][6]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jul_temp,id1=3,id2=6),self.Jul_temp)
        self.Drip=wx.StaticText(self,-1,'Jan: ',(300,posx))
        self.JanDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][0]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JanDrip_interval,id1=4,id2=0),self.JanDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Jul: ',(440,posx))
        self.JulDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][6]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JulDrip_interval,id1=4,id2=6),self.JulDrip_interval)
        self.Hum=wx.StaticText(self,-1,'Jan: ',(600,posx))
        self.JanHumidity=wx.TextCtrl(self,-1,str(self.data1[9][0]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JanHumidity,id1=9,id2=0),self.JanHumidity)
        self.Hum=wx.StaticText(self,-1,'Jul: ',(740,posx))
        self.JulHumidity=wx.TextCtrl(self,-1,str(self.data1[9][6]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JulHumidity,id1=9,id2=6),self.JulHumidity)
        posx+=25
        self.Feb=wx.StaticText(self,-1,'Feb: ',(10,posx))
        self.Feb_temp=wx.TextCtrl(self,-1,str(self.data1[3][1]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Feb_temp,id1=3,id2=1),self.Feb_temp)
        self.Aug=wx.StaticText(self,-1,'Aug: ',(150,posx))
        self.Aug_temp=wx.TextCtrl(self,-1,str(self.data1[3][7]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Aug_temp,id1=3,id2=7),self.Aug_temp)
        self.Drip=wx.StaticText(self,-1,'Feb: ',(300,posx))
        self.FebDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][1]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.FebDrip_interval,id1=4,id2=1),self.FebDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Aug: ',(440,posx))
        self.JulDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][7]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JulDrip_interval,id1=4,id2=7),self.JulDrip_interval)
        self.Hum=wx.StaticText(self,-1,'Feb: ',(600,posx))
        self.FebHumidity=wx.TextCtrl(self,-1,str(self.data1[9][1]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.FebHumidity,id1=9,id2=1),self.FebHumidity)
        self.Hum=wx.StaticText(self,-1,'Aug: ',(740,posx))
        self.AugHumidity=wx.TextCtrl(self,-1,str(self.data1[9][7]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.AugHumidity,id1=9,id2=7),self.AugHumidity)
        posx+=25
        self.Mar=wx.StaticText(self,-1,'Mar: ',(10,posx))
        self.Mar_temp=wx.TextCtrl(self,-1,str(self.data1[3][2]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Mar_temp,id1=3,id2=2),self.Mar_temp)
        self.Sep=wx.StaticText(self,-1,'Sep: ',(150,posx))
        self.Sep_temp=wx.TextCtrl(self,-1,str(self.data1[3][8]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Sep_temp,id1=3,id2=8),self.Sep_temp)
        self.Cave=wx.StaticText(self,-1,'Mar: ',(300,posx))
        self.MarDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][2]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MarDrip_interval,id1=4,id2=2),self.MarDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Sep: ',(440,posx))
        self.SepDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][8]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.SepDrip_interval,id1=4,id2=8),self.SepDrip_interval)
        self.Hum=wx.StaticText(self,-1,'Mar: ',(600,posx))
        self.MarHumidity=wx.TextCtrl(self,-1,str(self.data1[9][2]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MarHumidity,id1=9,id2=2),self.MarHumidity)
        self.Hum=wx.StaticText(self,-1,'Sep: ',(740,posx))
        self.SepHumidity=wx.TextCtrl(self,-1,str(self.data1[9][8]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.SepHumidity,id1=9,id2=8),self.SepHumidity)
        posx+=25
        self.Apr=wx.StaticText(self,-1,'Apr: ',(10,posx))
        self.Apr_temp=wx.TextCtrl(self,-1,str(self.data1[3][3]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Apr_temp,id1=3,id2=3),self.Apr_temp)
        self.Oct=wx.StaticText(self,-1,'Oct: ',(150,posx))
        self.Oct_temp=wx.TextCtrl(self,-1,str(self.data1[3][9]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Oct_temp,id1=3,id2=9),self.Oct_temp)
        self.Hum=wx.StaticText(self,-1,'Apr: ',(300,posx))
        self.AprDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][3]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.AprDrip_interval,id1=4,id2=3),self.AprDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Oct: ',(440,posx))
        self.OctDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][9]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.OctDrip_interval,id1=4,id2=9),self.OctDrip_interval)
        self.Hum=wx.StaticText(self,-1,'Apr: ',(600,posx))
        self.AprHumidity=wx.TextCtrl(self,-1,str(self.data1[9][3]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.AprHumidity,id1=9,id2=3),self.AprHumidity)
        self.Hum=wx.StaticText(self,-1,'Oct: ',(740,posx))
        self.OctHumidity=wx.TextCtrl(self,-1,str(self.data1[9][9]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.OctHumidity,id1=9,id2=9),self.OctHumidity)
        posx+=25
        self.May=wx.StaticText(self,-1,'May: ',(10,posx))
        self.May_temp=wx.TextCtrl(self,-1,str(self.data1[3][4]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.May_temp,id1=3,id2=4),self.May_temp)
        self.Nov=wx.StaticText(self,-1,'Nov: ',(150,posx))
        self.Nov_temp=wx.TextCtrl(self,-1,str(self.data1[3][10]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Nov_temp,id1=3,id2=10),self.Nov_temp)
        self.Wind=wx.StaticText(self,-1,'May : ',(300,posx))
        self.MayDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][4]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MayDrip_interval,id1=4,id2=4),self.MayDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Nov: ',(440,posx))
        self.NovDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][10]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.NovDrip_interval,id1=4,id2=10),self.NovDrip_interval)
        self.Hum=wx.StaticText(self,-1,'May: ',(600,posx))
        self.MayHumidity=wx.TextCtrl(self,-1,str(self.data1[9][4]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MayHumidity,id1=9,id2=4),self.MayHumidity)
        self.Hum=wx.StaticText(self,-1,'Nov: ',(740,posx))
        self.NovHumidity=wx.TextCtrl(self,-1,str(self.data1[9][10]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.NovHumidity,id1=9,id2=10),self.NovHumidity)
        posx+=25
        self.Jun=wx.StaticText(self,-1,'Jun: ',(10,posx))
        self.Jun_temp=wx.TextCtrl(self,-1,str(self.data1[3][5]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jun_temp,id1=3,id2=5),self.Jun_temp)
        self.Dec=wx.StaticText(self,-1,'Dec: ',(150,posx))
        self.Dec_temp=wx.TextCtrl(self,-1,str(self.data1[3][11]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Dec_temp,id1=3,id2=11),self.Dec_temp)
        self.Wind=wx.StaticText(self,-1,'Jun : ',(300,posx))
        self.JunDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][5]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JunDrip_interval,id1=4,id2=5),self.JunDrip_interval)
        self.Drip=wx.StaticText(self,-1,'Dec: ',(440,posx))
        self.DecDrip_interval=wx.TextCtrl(self,-1,str(self.data1[4][11]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.DecDrip_interval,id1=4,id2=11),self.DecDrip_interval)
        self.Hum=wx.StaticText(self,-1,'Jun: ',(600,posx))
        self.JunHumidity=wx.TextCtrl(self,-1,str(self.data1[9][5]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JunHumidity,id1=9,id2=5),self.JunHumidity)
        self.Hum=wx.StaticText(self,-1,'Dec: ',(740,posx))
        self.DecHumidity=wx.TextCtrl(self,-1,str(self.data1[9][11]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.DecHumidity,id1=9,id2=11),self.DecHumidity)
        
        #Drip and Cave pCO2 and Wind
        posx+=40
        self.label=wx.StaticText(self,-1,'Drip pCO2',(80,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(ppmv)',(205,posx+10))
        self.label=wx.StaticText(self,-1,'Cave pCO2',(350,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(ppmv)',(480,posx+10))
        self.label=wx.StaticText(self,-1,'Ventilation',(620,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(m/s)',(760,posx+10))

        posx+=40
        self.Jan=wx.StaticText(self,-1,'Jan: ',(10,posx))
        self.Jan_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][0]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jan_Drip_pCO2,id1=7,id2=0),self.Jan_Drip_pCO2)
        self.Jul=wx.StaticText(self,-1,'Jul: ',(150,posx))
        self.Jul_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][6]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jul_Drip_pCO2,id1=7,id2=6),self.Jul_Drip_pCO2)
        self.Drip=wx.StaticText(self,-1,'Jan: ',(300,posx))
        self.JanCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][0]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JanCave_pCO2,id1=8,id2=0),self.JanCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Jul: ',(440,posx))
        self.JulCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][6]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JulCave_pCO2,id1=8,id2=6),self.JulCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'Jan: ',(600,posx))
        self.Janv=wx.TextCtrl(self,-1,str(self.data1[10][0]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Janv,id1=10,id2=0),self.Janv)
        self.Hum=wx.StaticText(self,-1,'Jul: ',(740,posx))
        self.Julv=wx.TextCtrl(self,-1,str(self.data1[10][6]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Julv,id1=10,id2=6),self.Julv)
        posx+=25
        self.Feb=wx.StaticText(self,-1,'Feb: ',(10,posx))
        self.Feb_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][1]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Feb_Drip_pCO2,id1=7,id2=1),self.Feb_Drip_pCO2)
        self.Aug=wx.StaticText(self,-1,'Aug: ',(150,posx))
        self.Aug_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][7]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Aug_Drip_pCO2,id1=7,id2=7),self.Aug_Drip_pCO2)
        self.Drip=wx.StaticText(self,-1,'Feb: ',(300,posx))
        self.FebCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][1]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.FebCave_pCO2,id1=8,id2=1),self.FebCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Aug: ',(440,posx))
        self.JulCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][7]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JulCave_pCO2,id1=8,id2=7),self.JulCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'Feb: ',(600,posx))
        self.Febv=wx.TextCtrl(self,-1,str(self.data1[10][1]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Febv,id1=10,id2=1),self.Febv)
        self.Hum=wx.StaticText(self,-1,'Aug: ',(740,posx))
        self.Augv=wx.TextCtrl(self,-1,str(self.data1[10][7]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Augv,id1=10,id2=7),self.Augv)
        posx+=25
        self.Mar=wx.StaticText(self,-1,'Mar: ',(10,posx))
        self.Mar_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][2]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Mar_Drip_pCO2,id1=7,id2=2),self.Mar_Drip_pCO2)
        self.Sep=wx.StaticText(self,-1,'Sep: ',(150,posx))
        self.Sep_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][8]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Sep_Drip_pCO2,id1=7,id2=8),self.Sep_Drip_pCO2)
        self.Cave=wx.StaticText(self,-1,'Mar: ',(300,posx))
        self.MarCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][2]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MarCave_pCO2,id1=8,id2=2),self.MarCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Sep: ',(440,posx))
        self.SepCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][8]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.SepCave_pCO2,id1=8,id2=8),self.SepCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'Mar: ',(600,posx))
        self.Marv=wx.TextCtrl(self,-1,str(self.data1[10][2]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Marv,id1=10,id2=2),self.Marv)
        self.Hum=wx.StaticText(self,-1,'Sep: ',(740,posx))
        self.Sepv=wx.TextCtrl(self,-1,str(self.data1[10][8]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Sepv,id1=10,id2=8),self.Sepv)
        posx+=25
        self.Apr=wx.StaticText(self,-1,'Apr: ',(10,posx))
        self.Apr_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][3]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Apr_Drip_pCO2,id1=7,id2=3),self.Apr_Drip_pCO2)
        self.Oct=wx.StaticText(self,-1,'Oct: ',(150,posx))
        self.Oct_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][9]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Oct_Drip_pCO2,id1=7,id2=9),self.Oct_Drip_pCO2)
        self.Hum=wx.StaticText(self,-1,'Apr: ',(300,posx))
        self.AprCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][3]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.AprCave_pCO2,id1=8,id2=3),self.AprCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Oct: ',(440,posx))
        self.OctCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][9]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.OctCave_pCO2,id1=8,id2=9),self.OctCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'Apr: ',(600,posx))
        self.Aprv=wx.TextCtrl(self,-1,str(self.data1[10][3]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Aprv,id1=10,id2=3),self.Aprv)
        self.Hum=wx.StaticText(self,-1,'Oct: ',(740,posx))
        self.Octv=wx.TextCtrl(self,-1,str(self.data1[10][9]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Octv,id1=10,id2=9),self.Octv)
        posx+=25
        self.May=wx.StaticText(self,-1,'May: ',(10,posx))
        self.May_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][4]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.May_Drip_pCO2,id1=7,id2=4),self.May_Drip_pCO2)
        self.Nov=wx.StaticText(self,-1,'Nov: ',(150,posx))
        self.Nov_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][10]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Nov_Drip_pCO2,id1=7,id2=10),self.Nov_Drip_pCO2)
        self.Wind=wx.StaticText(self,-1,'May : ',(300,posx))
        self.MayCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][4]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.MayCave_pCO2,id1=8,id2=4),self.MayCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Nov: ',(440,posx))
        self.NovCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][10]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.NovCave_pCO2,id1=8,id2=10),self.NovCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'May: ',(600,posx))
        self.Mayv=wx.TextCtrl(self,-1,str(self.data1[10][4]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Mayv,id1=10,id2=4),self.Mayv)
        self.Hum=wx.StaticText(self,-1,'Nov: ',(740,posx))
        self.Novv=wx.TextCtrl(self,-1,str(self.data1[10][10]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Novv,id1=10,id2=10),self.Novv)
        posx+=25
        self.Jun=wx.StaticText(self,-1,'Jun: ',(10,posx))
        self.Jun_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][5]),(40,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Jun_Drip_pCO2,id1=7,id2=5),self.Jun_Drip_pCO2)
        self.Dec=wx.StaticText(self,-1,'Dec: ',(150,posx))
        self.Dec_Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][11]),(180,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Dec_Drip_pCO2,id1=7,id2=11),self.Dec_Drip_pCO2)
        self.Wind=wx.StaticText(self,-1,'Jun : ',(300,posx))
        self.JunCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][5]),(330,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.JunCave_pCO2,id1=8,id2=5),self.JunCave_pCO2)
        self.Drip=wx.StaticText(self,-1,'Dec: ',(440,posx))
        self.DecCave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][11]),(470,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.DecCave_pCO2,id1=8,id2=11),self.DecCave_pCO2)
        self.Hum=wx.StaticText(self,-1,'Jun: ',(600,posx))
        self.Junv=wx.TextCtrl(self,-1,str(self.data1[10][5]),(630,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Junv,id1=10,id2=5),self.Junv)
        self.Hum=wx.StaticText(self,-1,'Dec: ',(740,posx))
        self.Decv=wx.TextCtrl(self,-1,str(self.data1[10][11]),(770,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.Decv,id1=10,id2=11),self.Decv)
        posx+=70
        self.Wind=wx.StaticText(self,-1,'',(300,posx))
        
        #store size
        self.flabel=wx.StaticText(self,-1,'Store Sizes',(490,10),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(630,20))
        self.soils=wx.StaticText(self,-1,'Soil Store: ',(475,40))
        self.soilcap=wx.TextCtrl(self,-1,str(self.data1[5][0]),(575,40))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.soilcap,id1=5,id2=0),self.soilcap)
        self.epics=wx.StaticText(self,-1,'Epikarst: ',(475,65))
        self.epistore=wx.TextCtrl(self,-1,str(self.data1[5][3]),(575,65))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.epistore,id1=5,id2=3),self.epistore)
        self.ks1=wx.StaticText(self,-1,'Karst Store 1: ',(475,90))
        self.ks1_size=wx.TextCtrl(self,-1,str(self.data1[5][4]),(575,90))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ks1_size,id1=5,id2=4),self.ks1_size)
        self.ks2=wx.StaticText(self,-1,'Karst Store 2: ',(475,115))
        self.ks2_size=wx.TextCtrl(self,-1,str(self.data1[5][5]),(575,115))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ks2_size,id1=5,id2=5),self.ks2_size)
        
        #Overflow limits
        self.flabel=wx.StaticText(self,-1,'Overflow Limits',(475,140),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(665,150))
        self.epic=wx.StaticText(self,-1,'Epicap: ',(475,175))
        self.epicap=wx.TextCtrl(self,-1,str(self.data1[5][1]),(530,175))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.epicap,id1=5,id2=1),self.epicap)
        self.ovc=wx.StaticText(self,-1,'Ovicap: ',(475,200))
        self.ovicap=wx.TextCtrl(self,-1,str(self.data1[5][2]),(530,200))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ovicap,id1=5,id2=2),self.ovicap)
        
        #weibull parameters
        self.flabel=wx.StaticText(self,-1,'Weibull Distribution',(460,225),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        weibt="For diffuse flow; more info on wiki"
        self.weibtext=wx.StaticText(self,-1,weibt,(475,260))
        self.lam=wx.StaticText(self,-1,'Lambda (scale): ',(475,285))
        self.lambd=wx.TextCtrl(self,-1,str(self.data1[6][0]),(570,285))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.lambd,id1=6,id2=0),self.lambd)
        self.kay=wx.StaticText(self,-1,'k (shape): ',(475,310))
        self.kayp=wx.TextCtrl(self,-1,str(self.data1[6][1]),(570,310))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kayp,id1=6,id2=1),self.kayp)
        weibt2="Click button to see distribution"
        self.weibtext2=wx.StaticText(self,-1,weibt2,(475,335))
        button1=wx.Button(self,-1,"Weibull Dist",(500,350))
        self.Bind(wx.EVT_BUTTON,self.weibull,button1)
        
        pos_init=100
        #Store Initial Values
        self.label=wx.StaticText(self,-1,'Store Initial Values',(700,pos_init),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(800,pos_init+30))
        pos_init+=50
        self.ss=wx.StaticText(self,-1,'Soil Store: ',(710,pos_init))
        self.sstore=wx.TextCtrl(self,-1,str(self.data1[12][0]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.sstore,id1=12,id2=0),self.sstore)
        pos_init+=25
        self.es=wx.StaticText(self,-1,'Epikarst: ',(710,pos_init))
        self.estore=wx.TextCtrl(self,-1,str(self.data1[12][1]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.estore,id1=12,id2=1),self.estore)
        pos_init+=25
        self.ks1=wx.StaticText(self,-1,'Karst Store 1: ',(710,pos_init))
        self.kstore1=wx.TextCtrl(self,-1,str(self.data1[12][2]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore1,id1=12,id2=2),self.kstore1)
        pos_init+=25
        self.ks2=wx.StaticText(self,-1,'Karst Store 2: ',(710,pos_init))
        self.kstore2=wx.TextCtrl(self,-1,str(self.data1[12][3]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore2,id1=12,id2=3),self.kstore2)
        pos_init+=25
        self.dif=wx.StaticText(self,-1,'Diffuse values: ',(710,pos_init))
        self.dpdf=wx.TextCtrl(self,-1,str(self.data1[12][4]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.dpdf,id1=12,id2=5),self.dpdf)
        pos_init+=25
        #d18o initial values
        self.label=wx.StaticText(self,-1,'d18O Initial Values',(700,pos_init),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(per-mille)',(785,pos_init+30))
        pos_init+=50
        self.ss=wx.StaticText(self,-1,'Soil Store: ',(710,pos_init))
        self.sstored18o=wx.TextCtrl(self,-1,str(self.data1[13][0]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.sstored18o,id1=13,id2=0),self.sstored18o)
        pos_init+=25
        self.es=wx.StaticText(self,-1,'Epikarst: ',(710,pos_init))
        self.estored18o=wx.TextCtrl(self,-1,str(self.data1[13][1]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.estored18o,id1=13,id2=1),self.estored18o)
        pos_init+=25
        self.ks1=wx.StaticText(self,-1,'Karst Store 1: ',(710,pos_init))
        self.kstore1d18o=wx.TextCtrl(self,-1,str(self.data1[13][2]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore1d18o,id1=13,id2=2),self.kstore1d18o)
        pos_init+=25
        self.ks2=wx.StaticText(self,-1,'Karst Store 2: ',(710,pos_init))
        self.kstore2d18o=wx.TextCtrl(self,-1,str(self.data1[13][3]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore2d18o,id1=13,id2=3),self.kstore2d18o)
        pos_init+=25
        self.ind18o=wx.StaticText(self,-1,'Rain prev d18O: ',(710,pos_init))
        self.rain_inid18o=wx.TextCtrl(self,-1,str(self.data1[13][4]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.rain_inid18o,id1=13,id2=4),self.rain_inid18o)
        pos_init+=25
        self.dif=wx.StaticText(self,-1,'Diffuse values: ',(710,pos_init))
        self.dpdfd18o=wx.TextCtrl(self,-1,str(self.data1[13][5]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.dpdfd18o,id1=13,id2=5),self.dpdfd18o)
        pos_init+=25
        #Cave Mixing Parameter
        self.label=wx.StaticText(self,-1,'Cave Mixing (phi)',(700,pos_init),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(%)',(910,pos_init+10))
        pos_init+=40
        self.Mixing=wx.StaticText(self,-1,'Mixing Parameter: ',(710,pos_init))
        self.phi=wx.TextCtrl(self,-1,str(self.data1[11][0]),(830,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.phi,id1=11,id2=0),self.phi)
        
        
        
    def assign(self,event,name,id1,id2):
        #self.SetStatusText(name.GetValue())
        self.data1[id1][id2]=float(name.GetValue())
        csvoutput=open(self.config,'w')
        writer=csv.writer(csvoutput)
        for list in self.data1:
            writer.writerow(list)
        
    def weibull(self,event):
        x = np.linspace(0,2.5,25)
        v=s.exponweib(self.data1[6][0],self.data1[6][1])
        q=v.pdf(x)
        plt.plot(x,q)
        plt.show()
# -*- coding: utf-8 -*-

"""
Tests

Run with pytest
"""
import os
import sys
import pytest

# try and make this script run from more than one directory
sys.path.append('.')
sys.path.append('..')

# disable numba for debugging purposes
os.environ['NUMBA_DISABLE_JIT'] = '1'

from Karstolution.isotope_calcite import isotope_calcite
from Karstolution.calcpco2 import calc_pco2

def test_zero_net_flux_case():
    # this case was blowing up (net calcite deposition --> zero)
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=1314.6646646646648, pCO2cave=1000,
          h=0.98, V=0.1, phi=1.0, d18Oini=0, tt=1) )
    print(ic,gr)

def test_golgotha_parameters():
    # these parameters are not behaving as expected
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=6000., pCO2cave=6000.,
          h=0.98, V=0.1, phi=1.0, d18Oini=0, tt=1) )
    print(ic)

def test_no_evaporation():
    # evaporation set to zero - model should still run because precipitation
    # is driven by gradient in CO2 between drip water and cave air
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=1.0, V=0.0, phi=1.0, d18Oini=0, tt=1) )
    print(ic)

def test_zero_ventilation():
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.95, V=0.0, phi=1.0, d18Oini=0, tt=1) )
    print(ic,gr)


def test_evaporation_tends_to_zero():
    # evaporation is zero
    ic1,gr1 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=1.0, V=0.0, phi=1.0, d18Oini=0, tt=1) )
    # evaporation close to zero
    ic2,gr2 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.9999, V=0.0001, phi=1.0, d18Oini=0, tt=1) )
    # evaporation closer to zero
    ic3,gr3 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.999999, V=0.000001, phi=1.0, d18Oini=0, tt=1) )
    # evaporation close to zero, but h not 1.0
    ic4,gr4 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.7, V=1e-6, phi=1.0, d18Oini=0, tt=1) )
    print(ic1,ic2,ic3,ic4)
    assert abs(ic1-ic3) < 1e-4

def test_humidity_sensitivity_when_evaporation_tends_to_zero():
    ic1,gr1 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.99, V=1e-7, phi=1.0, d18Oini=0, tt=1) )
    # evaporation close to zero
    ic2,gr2 = (isotope_calcite(d=500., TC=10., pCO2=16000., pCO2cave=6000.,
          h=0.7, V=1e-7, phi=1.0, d18Oini=0, tt=1) )
    print(ic1,ic2)

def test_no_gradient_no_evaporation():
    # evaporation is zero (rh 100%, no ventilation), pCO2 equals pCAVE, there should
    # be nothing happening (ic,gr undefined)
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=6000., pCO2cave=6000.,
          h=1.0, V=0.0, phi=1.0, d18Oini=0, tt=1) )
    print(ic,gr)

def test_no_gradient_weak_evaporation():
    # evaporation is close to zero, pCO2 equals pCAVE
    ic,gr = (isotope_calcite(d=500., TC=10., pCO2=6000., pCO2cave=6000.,
          h=0.9999, V=0.0001, phi=1.0, d18Oini=0, tt=1) )
    print(ic,gr)

def test_calc_pco2():
    pco2 = calc_pco2(1e-3, 21.)
    print(pco2)

def test_calc_pco2_fail():
    with pytest.raises(ValueError):
      pco2 = calc_pco2(0, 21.)
    

if __name__ == "__main__":
    #test_no_evaporation()
    #test_evaporation_tends_to_zero()
    #test_no_gradient_weak_evaporation()
    #test_zero_ventilation()
    #test_humidity_sensitivity_when_evaporation_tends_to_zero()
    test_calc_pco2()
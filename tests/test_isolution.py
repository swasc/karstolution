import os
import sys

sys.path.append('..')

# disable numba for debugging purposes
os.environ['NUMBA_DISABLE_JIT'] = '1'


def test_isolution():
    from Karstolution.isotope_calcite import isotope_calcite
    # this case was blowing up (net calcite deposition --> zero)
    ic = (isotope_calcite(d=500., TC=10., pCO2=1314.6646646646648, pCO2cave=1000, 
          h=0.98, V=0.1, phi=1.0, d18Oini=0, tt=1) )
    print(ic)
    # these parameters are not behaving as expected
    ic = (isotope_calcite(d=500., TC=10., pCO2=6000., pCO2cave=6000., 
          h=0.98, V=0.1, phi=1.0, d18Oini=0, tt=1) )
    print(ic)


if __name__ == "__main__":
    test_isolution()
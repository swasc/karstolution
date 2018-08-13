# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import numpy as np

from . import constants


def calc_pco2(ca, TC):
    """
    Calculate pCO2-equivalent at a given Ca2+ concentration and temperature


    Inputs
    ------
        *ca* - scalar
        Calcium 2+ ion concentration (mol/l)

        *TC* - scalar
        Temperature in deg C.

    Outputs
    -------
        *pCO2* - scalar
        CO2 equivalent volume mixing ratio (ppmV), with a relative error of 
        about 1e-6

    Algorithm
    ---------
    The following description is taken from the Isolution paper (Deininger and 
    Scholz, submitted 2018).

    Calcpco2 converts Ca 2+ concentrations (given in mol/l) in a pCO 2 -equivalent
    using the mass laws of the CO 2 -H 2 O-CaCO 3 -system assuming a chemical 
    equilibrium between all chemical species. CALCPCO2.m does not consider any 
    other ions occurring in natural cave drip waters, such as Mg 2+ . CALCPCO2.m 
    firstly calculates the Ca 2+ concentration for pCO 2 values ranging from 0 to 
    1,000,000ppmV subdivided into ten equidistant intervals (i.e., the Ca 2+ 
    concentration for 0, 100,000, 200,000ppmV, etc.). In a second step, the function
    finds the interval mirroring the real Ca 2+ concentration (e.g., the interval 
    from 0 to 100,000ppmV). Then step one and two are repeated until the real and 
    the calculated Ca 2+ concentration are similar, i.e., the pCO 2 interval 
    mirroring the real Ca 2+ concentration is again subdivided into ten equidistant
    intervals, and then step 2 is repeated.

    History
    -------
    13 August 2018: Translated from Deininger's Matlab function
    """
    # terminate when error is smaller than this
    epsilon = (ca * 1e-6) + 1e-20
    # number of steps to divide search interval into
    # (must be at least 2)
    n_steps = 10
    # pco2 interval to search
    pco2int = np.linspace(0, 1000000e-6, n_steps)
    run = True
    n_iter = 0
    ca_int = np.empty(n_steps)

    while run:
        pco2_step_size = pco2int[1] - pco2int[0]
        for ii, pco2ii in enumerate(pco2int):
            concentrations = constants.constants(TC, pco2ii)
            ca_int[ii] = concentrations[2][0]

        # the search region must always bracket the target value
        if not (ca_int[0] <= ca and ca_int[-1] >= ca):
            raise ValueError(
                "Unable to compute pCO2 for Ca concentration of {} mol/l.".
                format(ca))

        index_pco2 = np.argmin(np.abs(ca_int - ca))
        # decide where to search next
        if ca_int[index_pco2] < ca:
            # search the interval with higher pCO2 next
            pco2int = np.linspace(pco2int[index_pco2],
                                  pco2int[index_pco2] + pco2_step_size,
                                  n_steps)
        else:
            # search the interval with lower pCO2 next
            pco2int = np.linspace(pco2int[index_pco2] - pco2_step_size,
                                  pco2int[index_pco2], n_steps)

        if np.abs(ca_int[index_pco2] - ca) < epsilon:
            run = False

        n_iter += 1
        print(n_iter, ca, ca_int[index_pco2], pco2int[index_pco2] * 1e6)
        # no convergence after many iterations: something has gone wrong
        # tests converged in < 10 iterations.
        assert (n_iter < 10000)
    return pco2int[index_pco2] * 1e6


# original matlab code
#function output = CALCPCO2(ca, TC)
#
#%input variable is Ca2+ concentration of the drip water in mol/l and the
#%ave air temperature in deg C
#
#pco2int = [0:100000e-6:1000000e-6];
#
#
# run = 1;    %varialbe to stop the while loop when set to 0
#
# while run==1
#
#     ca_int = [];
#
#     for i = 1:length(pco2int)
#
#          concentrations = KONSTANTEN(TC, pco2int(i));
#          ca_int(end+1) = concentrations{3}(1);
#
#     end
#
#     [next_pco2, index_pco2] = min(abs(ca_int-ca));
#
#     if ca_int(index_pco2)<ca
#         pco2int = [pco2int(index_pco2):(pco2int(index_pco2+1)-pco2int(index_pco2))/10:pco2int(index_pco2+1)];
#     else
#         pco2int = [pco2int(index_pco2-1):(pco2int(index_pco2)-pco2int(index_pco2-1))/10:pco2int(index_pco2)];
#     end
#
#     if ca_int(index_pco2)==ca
#         run=0;
#     end
#
# end
#
# output = pco2int(index_pco2)*1e6;
#
# end
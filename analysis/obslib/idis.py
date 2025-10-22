#!/usr/bin/env python
import sys, os
import numpy as np
import copy
import pandas as pd
import scipy as sp
from scipy.interpolate import griddata

## matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preview']=True
import pylab as py
from matplotlib.ticker import MultipleLocator

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
#from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def get_xbins(data,kind):
    query = {}
    if kind == 'p':
        query[26] = data.query('X > 4.80e-3 and X <= 5.70e-3') # x = 0.005 (i = 26)
        query[25] = data.query('X > 7.20e-3 and X <= 9.30e-3') # x = 0.008 (i = 25)
        query[24] = data.query('X > 1.15e-2 and X <= 1.40e-2') # x = 0.013 (i = 24)
        query[23] = data.query('X > 1.65e-2 and X <= 1.90e-2') # x = 0.018 (i = 23)
        query[22] = data.query('X > 2.30e-2 and X <= 2.90e-2') # x = 0.026 (i = 22)
        query[21] = data.query('X > 3.40e-2 and X <= 3.80e-2') # x = 0.036 (i = 21)
        query[20] = data.query('X > 4.65e-2 and X <= 5.40e-2') # x = 0.050 (i = 20)
        query[19] = data.query('X > 6.50e-2 and X <= 7.30e-2') # x = 0.070 (i = 19)
        query[18] = data.query('X > 8.50e-2 and X <= 9.20e-2') # x = 0.090 (i = 18)
        query[17] = data.query('X > 9.80e-2 and X <= 10.3e-2') # x = 0.100 (i = 17)
        query[16] = data.query('X > 10.8e-2 and X <= 11.3e-2') # x = 0.110 (i = 16)
        query[15] = data.query('X > 13.6e-2 and X <= 14.6e-2') # x = 0.140 (i = 15)
        query[14] = data.query('X > 17.1e-2 and X <= 18.7e-2') # x = 0.180 (i = 14)
        query[13] = data.query('X > 21.7e-2 and X <= 23.7e-2') # x = 0.230 (i = 13)
        query[12] = data.query('X > 26.0e-2 and X <= 29.0e-2') # x = 0.280 (i = 12)
        query[11] = data.query('X > 33.0e-2 and X <= 36.0e-2') # x = 0.350 (i = 11)
        query[10] = data.query('X > 42.0e-2 and X <= 45.0e-2') # x = 0.430 (i = 10)
        query[9]  = data.query('X > 45.0e-2 and X <= 48.0e-2') # x = 0.460 (i = 9 )
        query[8]  = data.query('X > 48.0e-2 and X <= 52.0e-2') # x = 0.500 (i = 8 )
        query[7]  = data.query('X > 52.0e-2 and X <= 55.0e-2') # x = 0.530 (i = 7 )
        query[6]  = data.query('X > 55.0e-2 and X <= 58.0e-2') # x = 0.560 (i = 6 )
        query[5]  = data.query('X > 58.0e-2 and X <= 60.0e-2') # x = 0.590 (i = 5 )
        query[4]  = data.query('X > 60.0e-2 and X <= 63.0e-2') # x = 0.610 (i = 4 )
        query[3]  = data.query('X > 63.0e-2 and X <= 66.0e-2') # x = 0.650 (i = 3 )
        query[2]  = data.query('X > 73.0e-2 and X <= 76.0e-2') # x = 0.750 (i = 2 )
        query[1]  = data.query('X > 84.0e-2 and X <= 86.0e-2') # x = 0.850 (i = 1 )
        query[0]  = data.query('X > 86.0e-2 and X <= 90.0e-2') # x = 0.880 (i = 0 )

    if kind == 'd':
        query[18] = data.query('X > 6.50e-2 and X <= 7.30e-2') # x = 0.070 (i = 18)
        query[17] = data.query('X > 8.50e-2 and X <= 9.20e-2') # x = 0.090 (i = 17)
        query[16] = data.query('X > 9.80e-2 and X <= 10.3e-2') # x = 0.100 (i = 16)
        query[15] = data.query('X > 13.6e-2 and X <= 14.6e-2') # x = 0.140 (i = 15)
        query[14] = data.query('X > 17.1e-2 and X <= 18.7e-2') # x = 0.180 (i = 14)
        query[13] = data.query('X > 21.7e-2 and X <= 23.7e-2') # x = 0.230 (i = 13)
        query[12] = data.query('X > 26.0e-2 and X <= 29.0e-2') # x = 0.280 (i = 12)
        query[11] = data.query('X > 33.0e-2 and X <= 36.0e-2') # x = 0.350 (i = 11)
        query[10] = data.query('X > 42.0e-2 and X <= 45.0e-2') # x = 0.430 (i = 10)
        query[9]  = data.query('X > 45.0e-2 and X <= 48.0e-2') # x = 0.460 (i = 9 )
        query[8]  = data.query('X > 48.0e-2 and X <= 52.0e-2') # x = 0.500 (i = 8 )
        query[7]  = data.query('X > 52.0e-2 and X <= 55.0e-2') # x = 0.530 (i = 7 )
        query[6]  = data.query('X > 55.0e-2 and X <= 58.0e-2') # x = 0.560 (i = 6 )
        query[5]  = data.query('X > 58.0e-2 and X <= 60.0e-2') # x = 0.590 (i = 5 )
        query[4]  = data.query('X > 60.0e-2 and X <= 63.0e-2') # x = 0.610 (i = 4 )
        query[3]  = data.query('X > 63.0e-2 and X <= 66.0e-2') # x = 0.650 (i = 3 )
        query[2]  = data.query('X > 73.0e-2 and X <= 76.0e-2') # x = 0.750 (i = 2 )
        query[1]  = data.query('X > 84.0e-2 and X <= 86.0e-2') # x = 0.850 (i = 1 )
        query[0]  = data.query('X > 86.0e-2 and X <= 90.0e-2') # x = 0.880 (i = 0 )

    if kind == 'd/p':
        query[17] = data.query('X > 4.80e-3 and X <= 5.70e-3') # x = 0.005 (i = 17)
        query[16] = data.query('X > 7.20e-3 and X <= 9.30e-3') # x = 0.008 (i = 16)
        query[15] = data.query('X > 1.15e-2 and X <= 1.40e-2') # x = 0.013 (i = 15)
        query[14] = data.query('X > 1.65e-2 and X <= 1.90e-2') # x = 0.018 (i = 14)
        query[13] = data.query('X > 2.30e-2 and X <= 2.90e-2') # x = 0.025 (i = 13)
        query[12] = data.query('X > 2.90e-2 and X <= 3.80e-2') # x = 0.035 (i = 12)
        query[11] = data.query('X > 4.60e-2 and X <= 5.40e-2') # x = 0.050 (i = 11) 
        query[10] = data.query('X > 6.50e-2 and X <= 7.30e-2') # x = 0.070 (i = 10)
        query[9]  = data.query('X > 8.50e-2 and X <= 9.20e-2') # x = 0.090 (i = 9 )
        query[8]  = data.query('X > 10.7e-2 and X <= 11.3e-2') # x = 0.110 (i = 8 )
        query[7]  = data.query('X > 13.6e-2 and X <= 14.6e-2') # x = 0.140 (i = 7 )
        query[6]  = data.query('X > 14.6e-2 and X <= 18.7e-2') # x = 0.180 (i = 6 )
        query[5]  = data.query('X > 18.7e-2 and X <= 23.7e-2') # x = 0.230 (i = 5 )
        query[4]  = data.query('X > 23.7e-2 and X <= 29.0e-2') # x = 0.280 (i = 4 )
        query[3]  = data.query('X > 29.0e-2 and X <= 36.0e-2') # x = 0.350 (i = 3 )
        query[2]  = data.query('X > 36.0e-2 and X <= 47.6e-2') # x = 0.450 (i = 2 )
        query[1]  = data.query('X > 47.6e-2 and X <= 57.5e-2') # x = 0.550 (i = 1 )
        query[0]  = data.query('X > 57.5e-2 and X <= 70.5e-2') # x = 0.680 (i = 0 )

    if kind == 'HERA NC':
        query[24] = data.query('X < 3.50e-5')                 # x = 2.47e-05,2.928e-05,3.088e-05                                                            (i = 24)
        query[23] = data.query('X > 3.50e-5 and X <= 5.5e-5')  # x = 3.66e-05,4.06e-05,4.09e-05,4.323e-05,4.60e-05,5e-05,5.124e-05,5.22e-05,5.31e-05         (i = 23)
        query[22] = data.query('X > 5.50e-5 and X <= 9.5e-5')  # x = 5.92e-05,6.176e-05,6.83e-05,7.32e-05,7.54e-05,8.e-05,8.029e-05,8.18e-05,8.55e-05        (i = 22)
        query[21] = data.query('X > 9.50e-5 and X <= 1.5e-4')  # x = 9.515e-05,9.86e-05,1.0499e-04,1.118e-04,1.2443e-04,1.29e-04,1.30e-04,1.39e-04,1.392e-04 (i = 21)
        query[20] = data.query('X > 1.50e-4 and X <= 2.5e-4')  # x = 1.61e-04,1.741e-04,1.821e-04,2.0e-04,2.09e-04,2.276e-04,2.37e-04                        (i = 20)
        query[19] = data.query('X > 2.50e-4 and X <= 3.5e-4')  # x = 2.68e-04,2.9e-04,3.14e-04,3.2e-04,3.28e-04,3.45e-04                                     (i = 19)
        query[18] = data.query('X > 3.50e-4 and X <= 5.5e-4')  # x = 3.55e-04,3.88e-04,4.1e-04,4.603e-04,5e-04,5.31e-04                                      (i = 18)
        query[17] = data.query('X > 5.50e-4 and X <= 6.6e-4')  # x = 5.90e-04,5.92e-04,6.16e-04,6.341e-04,6.57e-04                                           (i = 17)
        query[16] = data.query('X > 7.50e-4 and X <= 8.5e-4')  # x = 7.0e-04,8.0e-04                                                                         (i = 16)
        query[15] = data.query('X > 8.50e-4 and X <= 0.98e-3') # x = 8.6e-04,9.2e-04,9.22e-04,9.4e-04                                                        (i = 15)
        query[14] = data.query('X > 0.98e-3 and X <= 1.9e-3')  # x = 0.001,0.0011,0.00124,0.0013,0.0015, 0.0016, 0.00172,0.00188                             (i = 14)
        query[13] = data.query('X > 0.0019  and X <= 0.0030')  # x = 0.002,0.00212,0.0025,0.0026,0.0027                                                      (i = 13)
        query[12] = data.query('X > 0.0030  and X <= 0.0040')  # x = 00.0032,0.033,0.0039                                                                    (i = 12)
        query[11] = data.query('X > 4.80e-3 and X <= 5.6e-3')  # x = 0.0050,0.053, 0.0066                                                                    (i = 11)
        query[10] = data.query('X > 7.90e-3 and X <= 8.6e-3')  # x = 0.0080,0.0085                                                                           (i = 10)
        query[9]  = data.query('X > 1.00e-2 and X <= 1.5e-2')  # x = 1.05,  0.0130, 0.014                                                                    (i = 9 )
        query[8]  = data.query('X > 1.90e-2 and X <= 2.2e-2')  # x = 0.020, 0.0219                                                                           (i = 8 )
        query[7]  = data.query('X > 3.00e-2 and X <= 3.4e-2')  # x = 0.032                                                                                   (i = 7 )
        query[6]  = data.query('X > 4.80e-2 and X <= 5.6e-2')  # x = 0.050, 0.0470                                                                           (i = 6 )
        query[5]  = data.query('X > 7.90e-2 and X <= 8.8e-2')  # x = 0.080, 0.0875                                                                           (i = 5 )
        query[4]  = data.query('X > 12.0e-2 and X <= 14e-2')   # x = 0.13                                                                                    (i = 4 )
        query[3]  = data.query('X > 17.0e-2 and X <= 19e-2')   # x = 0.18                                                                                    (i = 3 )
        query[2]  = data.query('X > 23.0e-2 and X <= 26e-2')   # x = 0.25                                                                                    (i = 2 )
        query[1]  = data.query('X > 38.0e-2 and X <= 42.0e-2') # x = 0.40                                                                                    (i = 1 )
        query[0]  = data.query('X > 62.0e-2 and X <= 66e-2')   # x = 0.65                                                                                    (i = 0 )

    if kind == 'HERA CC':
        query[7]  = data.query('X > 0.0079 and X <= 0.0081') # x = 0.008 (i = 7) 
        query[6]  = data.query('X > 0.010  and X <= 0.015')  # x = 0.013 (i = 6) 
        query[5]  = data.query('X > 0.030  and X <= 0.034')  # x = 0.032 (i = 5) 
        query[4]  = data.query('X > 0.079  and X <= 0.081')  # x = 0.080 (i = 4) 
        query[3]  = data.query('X > 0.120  and X <= 0.14')   # x = 0.13  (i = 3) 
        query[2]  = data.query('X > 0.230  and X <= 0.27')   # x = 0.25  (i = 2)
        query[1]  = data.query('X > 0.380  and X <= 0.42')   # x = 0.40  (i = 1)
        query[0]  = data.query('X > 0.630  and X <= 0.66')   # x = 0.65  (i = 0)

    if kind == 'HERA other':
        query[29] = data.query('X < 3.50e-5 ')                # x = 3.27e-05                                                                               (i = 29)
        query[28] = data.query('X > 3.50e-5 and X <= 9.5e-5')  # x = 4.09e-05,5e-05,5.73e-05,8e-05,8.18e-05                                                 (i = 28)
        query[27] = data.query('X > 9.50e-5 and X <= 2.5e-4')  # x = 9.86e-05,1.3e-04,1.39e-04,1.61e-04,2e-04,2.46e-04                                      (i = 27)
        query[26] = data.query('X > 2.50e-4 and X <= 3.5e-4')  # x = 2.68e-04,3.2e-04,3.28e-04,3.35e-04                                                     (i = 26)
        query[25] = data.query('X > 3.50e-4 and X <= 5.5e-4')  # x = 4.1e-04,5e-04                                                                          (i = 25)
        query[24] = data.query('X > 5.50e-4 and X <= 6.6e-4')  # x = 5.74e-04                                                                               (i = 24)
        query[23] = data.query('X > 6.90e-4 and X <= 8.5e-4')  # x = 8e-04                                                                                  (i = 23)
        query[22] = data.query('X > 8.50e-4 and X <= 0.98e-3') # x =  8.8e-04,9.1e-04,9.206e-04,9.344e-04,  9.545e-04                                       (i = 22)
        query[21] = data.query('X > 0.98e-3 and X <= 0.0014')  # x = 1.3e-03                                                                                (i = 21)
        query[20] = data.query('X > 1.3e-3  and X <= 1.9e-03') # 1 . 3660e-03,1.392e-03,1.397e-03,1.409-03,1.46e-03,1.479e-03,1.578e-03,1.585e-03,1.591e-03 (i = 20)
        query[19] = data.query('X > 0.0019  and X <= 0.0029')  # x = 2e-03,2.12e-03                                                                         (i = 19)
        query[18] = data.query('X > 0.0029  and X <= 0.0040')  # x = 3.2e-03                                                                                (i = 18)
        query[17] = data.query('X > 4.0e-3  and X <= 5.55e-3') # x = 5e-03                                                                                  (i = 17)
        query[16] = data.query('X > 5.55e-3 and X <= 7.95e-3') # x = 5.6e-03, 5.727e-03, 5.754e-03, 5.8e-03, 5.9e-03, 5.918e-03, 6e-03, 6.1e-03,6.2e-03,6.4e-03,6.6e-03,6.9e-03,7.3e-03,7.398e-03,7.4e-03,7.6e-03,7.9e-03                                                                                                                                (i = 16)
        query[15] = data.query('X > 7.95e-3 and X <= 8.9e-3')  # x = 8e-03,8.5e-03                                                                          (i = 15)
        query[14] = data.query('X > 9.00e-3 and X <= 9.9e-3')  # x = 9.1e-03,9.3e-03,9.864e-03                                                              (i = 14)
        query[13] = data.query('X > 9.90e-3 and X <= 1.25e-2') # x = 1e-02,1.04e-02,1.05e-02,1.09e-02,1.16e-02,1.21e-02                                     (i = 13)
        query[12] = data.query('X > 1.25e-2 and X <= 1.5e-2')  # x = 1.3e-02,1.4e-02                                                                        (i = 12)
        query[11] = data.query('X > 1.50e-2 and X <= 1.7e-2')  # x = 1.51e-02,1.52e-02,1.61e-02,1.660e-02                                                   (i = 11)
        query[10] = data.query('X > 1.70e-2 and X <= 1.99e-2') # x = 1.71e-02,1.85e-02,1.97e-02                                                             (i = 10)
        query[9] =  data.query('X > 1.99e-2 and X <= 2.2e-2')  # x = 2e-02                                                                                  (i = 9)
        query[8] =  data.query('X > 2.30e-2 and X <= 2.7e-2')  # x = 2.420e-02,2.61e-02                                                                     (i = 8)
        query[7] =  data.query('X > 3.00e-2 and X <= 3.4e-2')  # x = 3.20e-02                                                                               (i = 7)
        query[6] =  data.query('X > 4.80e-2 and X <= 5.6e-2')  # x = 0.050                                                                                  (i = 6)
        query[5] =  data.query('X > 7.90e-2 and X <= 8.8e-2')  # x = 0.080                                                                                  (i = 5)
        query[4] =  data.query('X > 12.0e-2 and X <= 14e-2')   # x = 0.13                                                                                   (i = 4)
        query[3] =  data.query('X > 17.0e-2 and X <= 19e-2')   # x = 0.18                                                                                   (i = 3)
        query[2] =  data.query('X > 23.0e-2 and X <= 26e-2')   # x = 0.25                                                                                   (i = 2)
        query[1] =  data.query('X > 30.0e-2 and X <= 45e-2')   # x = 0.4                                                                                    (i = 1)
        query[0] =  data.query('X > 60.0e-2 and X <= 70e-2')   # x = 0.65                                                                                   (i = 0)


    return query

def get_Q2bins(data,kind):
    query = {}
    if kind == 'BONuS':
        query[5]  = data.query('Q2 > 3.7')              # Q2 = 4.0 (i = 5)
        query[4]  = data.query('Q2 > 3.2 and Q2 <= 3.7') # Q2 = 3.4 (i = 4)
        query[3]  = data.query('Q2 > 2.7 and Q2 <= 3.2') # Q2 = 2.9 (i = 3)
        query[2]  = data.query('Q2 > 2.3 and Q2 <= 2.7') # Q2 = 2.4 (i = 2)
        query[1]  = data.query('Q2 > 1.9 and Q2 <= 2.3') # Q2 = 2.0 (i = 1)
        query[0]  = data.query('Q2 > 1.5 and Q2 <= 1.9') # Q2 = 1.7 (i = 0)
    if kind == 'MARATHON' or kind == 'JLab':
        query[0]  = data.query('Q2 > 1.5 and Q2 <= 20') # Q2 = 1.7 (i = 0)
    if kind == 'SLAC E140x':
        query[3]  = data.query('Q2 > 3.1')               # Q2 = 3.6 (i = 3)
        query[2]  = data.query('Q2 > 2.0 and Q2 <= 3.1') # Q2 = 3.0 (i = 2)
        query[1]  = data.query('Q2 > 0.9 and Q2 <= 2.0') # Q2 = 1.0 (i = 1)
        query[0]  = data.query('Q2 > 0.0 and Q2 <= 0.9') # Q2 = 0.5 (i = 0)
    if kind == 'JLab E06-009':
        query[2]  = data.query('Q2 > 3.5')               # Q2 = 4.0 (i = 2)
        query[1]  = data.query('Q2 > 2.5 and Q2 <= 3.5') # Q2 = 3.0 (i = 1)
        query[0]  = data.query('Q2 > 0.0 and Q2 <= 2.5') # Q2 = 2.0 (i = 0)

    if kind == 'clas6 p Q2 < 3' or kind == 'clas6 d Q2 < 3':
        query[6] = data.query('Q2 > 2.75 and Q2 <= 3.00') # (i = 0)
        query[5] = data.query('Q2 > 2.50 and Q2 <= 2.75') # (i = 0)
        
        query[4] = data.query('Q2 > 2.35 and Q2 <= 2.50') # (i = 0)
        query[3] = data.query('Q2 > 2.25 and Q2 <= 2.35') # (i = 0)
        
        query[2] = data.query('Q2 > 2.00 and Q2 <= 2.25') # (i = 0)
        query[1] = data.query('Q2 > 1.75 and Q2 <= 2.00') # (i = 0)
        query[0] = data.query('Q2 > 1.50 and Q2 <= 1.75') # (i = 0)
        #query[5]  = data.query('Q2 > 1.25 and Q2 <= 1.50') # (i = 0)
        #query[4]  = data.query('Q2 > 1.00 and Q2 <= 1.25') # (i = 0)
        #query[3]  = data.query('Q2 > 0.75 and Q2 <= 1.00') # (i = 0)
        #query[2]  = data.query('Q2 > 0.50 and Q2 <= 0.75') # (i = 0)
        #query[1]  = data.query('Q2 > 0.25 and Q2 <= 0.50') # (i = 0)
        #query[0]  = data.query('Q2 > 0.00 and Q2 <= 0.25') # (i = 0)

    if kind == 'clas6 p Q2 > 3' or kind == 'clas6 d Q2 > 3':
        query[6] = data.query('Q2 > 4.50 and Q2 <= 4.75') # (i = 0)
        query[5] = data.query('Q2 > 4.25 and Q2 <= 4.50') # (i = 0)
        query[4] = data.query('Q2 > 4.00 and Q2 <= 4.25') # (i = 0)
        query[3] = data.query('Q2 > 3.75 and Q2 <= 4.00') # (i = 0)
        query[2] = data.query('Q2 > 3.50 and Q2 <= 3.75') # (i = 0)
        query[1] = data.query('Q2 > 3.25 and Q2 <= 3.50') # (i = 0)
        query[0] = data.query('Q2 > 3.00 and Q2 <= 3.25') # (i = 0)

    if kind == 'BONuS n E4.223':
        query[4]  = data.query('Q2 > 3.2')               # Q2 = 3.4 (i = 4)
        query[3]  = data.query('Q2 > 2.7 and Q2 <= 3.2') # Q2 = 2.9 (i = 3)
        query[2]  = data.query('Q2 > 2.3 and Q2 <= 2.7') # Q2 = 2.4 (i = 2)
        query[1]  = data.query('Q2 > 1.9 and Q2 <= 2.3') # Q2 = 2.0 (i = 1)
        query[0]  = data.query('Q2 > 1.5 and Q2 <= 1.9') # Q2 = 1.7 (i = 0)

    if kind == 'BONuS n E5.262':
        query[6]  = data.query('Q2 > 4.5')               # Q2 = 4.8 (i = 6)
        query[5]  = data.query('Q2 > 3.7 and Q2 <= 4.5') # Q2 = 4.0 (i = 5)
        query[4]  = data.query('Q2 > 3.2 and Q2 <= 3.7') # Q2 = 3.4 (i = 4)
        query[3]  = data.query('Q2 > 2.7 and Q2 <= 3.2') # Q2 = 2.9 (i = 3)
        query[2]  = data.query('Q2 > 2.3 and Q2 <= 2.7') # Q2 = 2.4 (i = 2)
        query[1]  = data.query('Q2 > 1.9 and Q2 <= 2.3') # Q2 = 2.0 (i = 1)
        query[0]  = data.query('Q2 > 1.5 and Q2 <= 1.9') # Q2 = 1.7 (i = 0)

    if kind == 'E665 p' or kind == 'E665 d':
        query[5] = data.query('Q2 > 50')               
        query[4] = data.query('Q2 > 40 and Q2 <= 50') 
        query[3] = data.query('Q2 > 12 and Q2 <= 21') 
        query[2] = data.query('Q2 > 6.0 and Q2 <= 12') 
        query[1] = data.query('Q2 > 3.0 and Q2 <= 6.0') 
        query[0] = data.query('Q2 > 1.5 and Q2 <= 3.0') 

    if kind == 'JLab JLCEE96 p':
        query[8] = data.query('Q2 > 3.1') 
        query[7] = data.query('Q2 > 2.6 and Q2 <= 3.1') 
        query[6] = data.query('Q2 > 2.2 and Q2 <= 2.6') 
        query[5] = data.query('Q2 > 1.8 and Q2 <= 2.2') 
        query[4] = data.query('Q2 > 1.5 and Q2 <= 1.8') 
        query[3] = data.query('Q2 > 1.3 and Q2 <= 1.5') 
        query[2] = data.query('Q2 > 1.0 and Q2 <= 1.3') 
        query[1] = data.query('Q2 > 0.5 and Q2 <= 1.0') 
        query[0] = data.query('Q2 > 0.0 and Q2 <= 0.5') 

    if kind == 'JLab JLCEE96 d':
        query[7] = data.query('Q2 > 3.1') 
        query[6] = data.query('Q2 > 2.8 and Q2 <= 3.1') 
        query[5] = data.query('Q2 > 2.2 and Q2 <= 2.6') 
        query[4] = data.query('Q2 > 1.8 and Q2 <= 2.2') 
        query[3] = data.query('Q2 > 1.5 and Q2 <= 1.8') 
        query[2] = data.query('Q2 > 1.0 and Q2 <= 1.5') 
        query[1] = data.query('Q2 > 0.5 and Q2 <= 1.0') 
        query[0] = data.query('Q2 > 0.0 and Q2 <= 0.5') 


    return query

def get_thetabins(data,kind):
    query = {}
    if kind == 'JLab':
        query[5]  = data.query('theta > 65 and theta <= 75') # theta = 70 (i = 5)
        query[4]  = data.query('theta > 57 and theta <= 65') # theta = 60 (i = 4)
        query[3]  = data.query('theta > 50 and theta <= 57') # theta = 55 (i = 3)
        query[2]  = data.query('theta > 43 and theta <= 50') # theta = 45 (i = 2)
        query[1]  = data.query('theta > 40 and theta <= 43') # theta = 41 (i = 1)
        query[0]  = data.query('theta > 35 and theta <= 40') # theta = 38 (i = 0)
    if kind == 'HallC':
        query[4]  = data.query('thetac > 38 and thetac <= 40') # theta = 39 (i = 5)
        query[3]  = data.query('thetac > 32 and thetac <= 34') # theta = 33 (i = 4)
        query[2]  = data.query('thetac > 28 and thetac <= 30') # theta = 29 (i = 3)
        query[1]  = data.query('thetac > 24 and thetac <= 26') # theta = 25 (i = 2)
        query[0]  = data.query('thetac <= 26')                 # theta = 21 (i = 1)
    if kind == 'JLab E03-103 p':
        query[1]  = data.query('theta > 45')                 # theta = 50 (i = 1)
        query[0]  = data.query('theta > 38 and theta <= 45') # theta = 40 (i = 0)
    if kind=='JLab E03-103 d E5.011':
        query[3]  = data.query('theta > 40')                 # theta = 46 (i = 3)
        query[2]  = data.query('theta > 30 and theta <= 40') # theta = 36 (i = 2)
        query[1]  = data.query('theta > 25 and theta <= 30') # theta = 29 (i = 1)
        query[0]  = data.query('theta > 0  and theta <= 25') # theta = 24 (i = 0)
    if kind=='JLab E03-103 d E5.766':
        query[5]  = data.query('theta > 45')                 # theta = 50 (i = 5)
        query[4]  = data.query('theta > 38 and theta <= 45') # theta = 40 (i = 4)
        query[3]  = data.query('theta > 30 and theta <= 38') # theta = 32 (i = 3)
        query[2]  = data.query('theta > 24 and theta <= 30') # theta = 26 (i = 2)
        query[1]  = data.query('theta > 20 and theta <= 24') # theta = 22 (i = 1)
        query[0]  = data.query('theta > 0  and theta <= 20') # theta = 18 (i = 0)

    return query


def get_plot(query,cluster_i=0):

    #--generate dictionary with everything needed for plot

    plot = {_:{} for _ in ['theory','X','Q2','theta','value','alpha','std']}
    for key in query:
        theory = query[key]['thy-%d' % cluster_i]
        std    = query[key]['dthy-%d' % cluster_i]
        X      = query[key]['X']
        Q2     = query[key]['Q2']
        if 'theta' in query[key]: theta = query[key]['theta']
        value  = query[key]['value']
        alpha  = query[key]['alpha']
        #--sort by ascending Q2
        zx = sorted(zip(Q2,X))
        zt = sorted(zip(Q2,theory))
        zv = sorted(zip(Q2,value))
        za = sorted(zip(Q2,alpha))
        zs = sorted(zip(Q2,std))
        if 'theta' in query[key]: ztheta = sorted(zip(Q2,theta))
        plot['X'][key]      = np.array([zx[i][1] for i in range(len(zx))])
        plot['theory'][key] = np.array([zt[i][1] for i in range(len(zt))])
        plot['value'][key]  = np.array([zv[i][1] for i in range(len(zv))])
        plot['alpha'][key]  = np.array([za[i][1] for i in range(len(za))])
        plot['std'][key]    = np.array([zs[i][1] for i in range(len(zs))])
        plot['Q2'][key]     = np.array(sorted(Q2))
        if 'theta' in query[key]: 
            plot['theta'][key] = np.array([ztheta[i][1] for i in range(len(ztheta))])

    return plot

def get_theory(PLOT,nbins,loop=True,funcQ2=True,funcX=False,functheta=False):

    #--interpolate theoretical values across Q2 or X

    theory = {}

    if funcQ2:     svar = 'Q2'
    if funcX:      svar = 'X'
    if functheta:  svar = 'theta'

    theory = {_:{} for _ in [svar,'value','std']}
    for key in range(nbins):
        var = []
        thy = []
        std = []
        
        #--if plotting for multiple experiments, loop over and combine
        if loop:
            for exp in PLOT:
                var. extend(PLOT[exp][svar][key])
                thy.extend(PLOT[exp]['theory'][key])
                std.extend(PLOT[exp]['std'][key])
        else:
            var.extend(PLOT[svar][key])
            thy.extend(PLOT['theory'][key])
            std.extend(PLOT['std'][key])

        #--if nothing in bin, skip
        if len(var) == 0: continue

        vmin = np.min(var)
        vmax = np.max(var)
        theory[svar][key]  = np.geomspace(vmin,vmax,100)

        #--if more than one value, interpolate between them
        if len(var) > 1:
            theory['value'][key] = griddata(np.array(var),np.array(thy),theory[svar][key],method='linear')
            theory['std'][key]   = griddata(np.array(var),np.array(std),theory[svar][key],method='linear')
        else:
            theory['value'][key] = np.ones(100)*thy 
            theory['std'][key]   = np.ones(100)*std


    return theory

def get_details(exp):

    #--get details for plotting

    if exp=='NMC' :          color, marker, ms = 'firebrick', '*', 6
    if exp=='SLAC' :         color, marker, ms = 'darkgreen', '^', 6
    if exp=='BCDMS':         color, marker, ms = 'blue'     , 'o', 6

    if exp=='HERA 10026':    color, marker, ms = 'black'    , '.', 8 
    if exp=='HERA 10030':    color, marker, ms = 'firebrick', 'D', 8 
    if exp=='HERA 10031':    color, marker, ms = 'blue'     , '^', 6 
    if exp=='HERA 10032':    color, marker, ms = 'darkgreen', 's', 6 

    if exp=='HERA 10027':    color, marker, ms = 'darkgreen', 's', 6 
    if exp=='HERA 10028':    color, marker, ms = 'firebrick', '.', 8 
    if exp=='HERA 10029':    color, marker, ms = 'blue',      '^', 8 

    if exp=='JLab d':        color, marker, ms = 'firebrick', '*', 8 
    if exp=='JLab p':        color, marker, ms = 'darkgreen', '^', 8 
    if exp=='BONuS' :        color, marker, ms = 'blue'     , 'o', 6 

    if exp=='MARATHON dp':   color, marker, ms = 'blue', '*', 5 
    if exp=='MARATHON ht':   color, marker, ms = 'darkgreen', '^', 5 
    if exp=='MARATHON hd':   color, marker, ms = 'darkgreen', '^', 5 
    if exp=='MARATHON td':   color, marker, ms = 'darkgreen', '^', 5 
 
    if exp=='JLab hd':       color, marker, ms = 'firebrick', 'o', 5
    if exp=='JLab dp0':      color, marker, ms = 'firebrick', 'o', 5
    if exp=='JLab dp1':      color, marker, ms = 'firebrick', 'o', 5

    if exp=='SLAC E140x p':  color, marker, ms = 'firebrick', 'o', 5
    if exp=='SLAC E140x d':  color, marker, ms = 'darkgreen', '^', 5

    if exp=='JLab E06-009':  color, marker, ms = 'firebrick', 'o', 5

    if exp=='JLab E03-103 p':         color, marker, ms = 'firebrick', 'o', 5
    if exp=='JLab E03-103 d E5.011':  color, marker, ms = 'darkgreen', '^', 5
    if exp=='JLab E03-103 d E5.766':  color, marker, ms = 'darkgreen', '^', 5


    if exp=='clas6 p Q2 < 3':  color, marker, ms = 'firebrick', 'o', 5
    if exp=='clas6 d Q2 < 3':  color, marker, ms = 'firebrick', 'o', 5
    if exp=='clas6 p Q2 > 3':  color, marker, ms = 'firebrick', 'o', 5
    if exp=='clas6 d Q2 > 3':  color, marker, ms = 'firebrick', 'o', 5

    if exp=='BONuS n E4.223': color, marker, ms = 'firebrick', 'o', 5
    if exp=='BONuS n E5.262': color, marker, ms = 'firebrick', 'o', 5

    if exp=='E665 p':  color, marker, ms = 'firebrick', 'o', 5
    if exp=='E665 d':  color, marker, ms = 'firebrick', 'o', 5

    if exp=='JLab JLCEE96 p':  color, marker, ms = 'firebrick', 'o', 0 
    if exp=='JLab JLCEE96 d':  color, marker, ms = 'firebrick', 'o', 0 


    return color,marker,ms



def plot_proton(wdir, data):


    nrows, ncols = 2, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)
    ax21 = py.subplot(nrows, ncols, 3)
    ax22 = py.subplot(nrows, ncols, 4)

    #################################
    #--Plot F2p from BCDMS, SLAC, NMC
    #################################

    nbins = 27

    if 10020 and 10010 and 10016 in data:
        nmc   = data[10020]  ## NMC p
        slac  = data[10010]  ## SLAC p
        bcdms = data[10016]  ## BCDMS p

        DATA = {}
        DATA['NMC']   = pd.DataFrame(nmc)
        DATA['SLAC']  = pd.DataFrame(slac)
        DATA['BCDMS'] = pd.DataFrame(bcdms)

        PLOT = {}
        for exp in DATA:
            query = get_xbins(DATA[exp],'p')
            PLOT[exp] = get_plot(query)

        theory = get_theory(PLOT,nbins)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                Q2    = PLOT[exp]['Q2'][key]
                val   = PLOT[exp]['value'][key]*2.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax11.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory['value']:
            Q2   = theory['Q2'][key]
            mean = theory['value'][key]*2.0**float(key)
            std  = theory['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax11.plot(Q2,mean,linestyle='solid',color='black')
            thy_band  = ax11.fill_between(Q2,down,up,color='gold',alpha=1.0)

        #--plot labels
        ax11.text(2.4,   3.0e7,  r'$x=0.005\, (i=24)$', fontsize = 22)
        ax11.text(3.7,   1.2e7,  r'$x=0.008$',          fontsize = 22)
        ax11.text(5.8,   6.0e6,  r'$x=0.013$',          fontsize = 22)
        ax11.text(7.5,   3.0e6,  r'$x=0.018$',          fontsize = 22)
        ax11.text(12.0,  1.5e6,  r'$x=0.026$',          fontsize = 22)
        ax11.text(15.0,  8.0e5,  r'$x=0.036$',          fontsize = 22)
        ax11.text(21.0,  4.0e5,  r'$x=0.05$',           fontsize = 22)
        ax11.text(28.0,  2.0e5,  r'$x=0.07$',           fontsize = 22)
        ax11.text(28.0,  1.0e5,  r'$x=0.09$',           fontsize = 22)
        ax11.text(41.0,  5.0e4,  r'$x=0.10$',           fontsize = 22)
        ax11.text(38.0,  2.5e4,  r'$x=0.11$',           fontsize = 22)
        ax11.text(62.0,  1.0e4,  r'$x=0.14$',           fontsize = 22)
        ax11.text(70.0,  5.0e3,  r'$x=0.18$',           fontsize = 22)
        ax11.text(93.0,  2.5e3,  r'$x=0.23$',           fontsize = 22)
        ax11.text(125.0, 1.1e3,  r'$x=0.28$',           fontsize = 22)
        ax11.text(150.0, 4.0e2,  r'$x=0.35$',           fontsize = 22)
        ax11.text(200.0, 1.0e2,  r'$x=0.43$',           fontsize = 22)
        ax11.text(65.0,  60.0,   r'$x=0.46$',           fontsize = 22)
        ax11.text(15.0,  30.0,   r'$x=0.50$',           fontsize = 22)
        ax11.text(250.0, 8.0,    r'$x=0.53$',           fontsize = 22)
        ax11.text(22.0,  5.0,    r'$x=0.56$',           fontsize = 22)
        ax11.text(23.0,  2.5,    r'$x=0.59$',           fontsize = 22)
        ax11.text(24.0,  1.0,    r'$x=0.61$',           fontsize = 22)
        ax11.text(250.0, 0.2,    r'$x=0.65$',           fontsize = 22)
        ax11.text(250.0, 0.03,   r'$x=0.75$',           fontsize = 22)
        ax11.text(30.0,  0.008,  r'$x=0.85$',           fontsize = 22)
        ax11.text(30.0 , 0.002,  r'$x=0.88\, (i=0)$',   fontsize = 22)

        ax11.semilogy()
        ax11.semilogx()
        ax11.set_xlim(1.5, 5e2)
        ax11.set_ylim(0.0005, 1e8)

        ax11.tick_params(axis = 'both', labelsize = 30)

        ax11.yaxis.set_tick_params(which = 'major', length = 10)
        ax11.set_yticks([0.001,0.01, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])
        ax11.set_yticklabels([r'$10^{-3}$',r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$',r'$10^8$'])
        locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
        ax11.yaxis.set_minor_locator(locmin)
        ax11.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax11.yaxis.set_tick_params(which = 'minor', length = 5)

        ax11.xaxis.set_tick_params(which = 'major', length = 10)
        ax11.xaxis.set_tick_params(which = 'minor', length = 5)
        ax11.set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2])
        xtick_labels = [2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2]
        ax11.set_xticklabels(['$%0.0f$' % x for x in xtick_labels])

        ax11.text(0.05,0.05,r'\boldmath$F_2^p$',      transform = ax11.transAxes, size = 60)
        ax11.text(0.15,0.05,r'$(\times\, 2^{\, i})$', transform = ax11.transAxes, size = 40)

        ax11.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

        ax11.yaxis.set_ticks_position('both')
        ax11.xaxis.set_ticks_position('both')
        ax11.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
        ax11.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

        handles = [hand['BCDMS'],hand['NMC'],hand['SLAC'],(thy_band,thy_plot)]
        label1  = r'\textbf{\textrm{BCDMS}}'
        label2  = r'\textbf{\textrm{NMC}}'
        label3  = r'\textbf{\textrm{SLAC}}'
        label4  = r'\textbf{\textrm{JAM}}'
        labels  = [label1,label2,label3,label4]
        ax11.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #################################
    #--Plot F2p from Hall C
    #################################

    if 10003 in data:
        jlab_p = data[10003] ## JLab p 

        E  = 5.5 #--Elab
        Mp = 0.939

        theta = lambda X,Q2,M: np.arcsin(np.sqrt(Q2/(4*E**2 - 2*E*Q2/M/X)))*2*180/np.pi

        data[10003]['theta'] = theta(data[10003]['X'],data[10003]['Q2'],Mp)

        nbins = 6 

        DATA = {}
        DATA['JLab p'] = pd.DataFrame(jlab_p) 

        theory = {}
        PLOT = {}
        for exp in DATA:
            query = get_thetabins(DATA[exp],'JLab')
            PLOT[exp] = get_plot(query)

            theory[exp] = get_theory(PLOT[exp],nbins,loop=False)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                Q2    = PLOT[exp]['Q2'][key]
                val   = PLOT[exp]['value'][key]*2.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax12.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            for key in theory[exp]['value']:
                Q2   = theory[exp]['Q2'][key]
                mean = theory[exp]['value'][key]*2.0**float(key)
                std  = theory[exp]['std'][key]  *2.0**float(key)
                down = mean - std
                up   = mean + std
                thy_plot ,= ax12.plot(Q2,mean,linestyle='solid',color='black')
                thy_band  = ax12.fill_between(Q2,down,up,color='gold', alpha = 1.0)


        ax12.text(6.55,0.90,   r'$\theta=70^{\circ} \, (i=5)$',   fontsize = 30)
        ax12.text(6.2, 0.50,   r'$\theta=60^{\circ} $',           fontsize = 30)
        ax12.text(5.9, 0.28,   r'$\theta=55^{\circ} $',           fontsize = 30)
        ax12.text(5.2, 0.18,   r'$\theta=45^{\circ} $',           fontsize = 30)
        ax12.text(4.9, 0.10,   r'$\theta=41^{\circ} $',           fontsize = 30)
        ax12.text(4.6, 0.055,  r'$\theta=38^{\circ}  \, (i=0)$',  fontsize = 30)

        ax12.semilogy()

        ax12.set_xlim(3.5,   7.5)
        ax12.set_ylim(0.04,  3.5)

        ax12.set_xticks([4, 5, 6, 7])
        ax12.set_xticklabels([r'$4$', r'$5$', r'$6$', r'$7$'])
        lo12cmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
        ax12.yaxis.set_minor_locator(locmin)
        ax12.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax12.yaxis.set_tick_params(which = 'minor', length = 5)
        ax12.set_yticks([1e-1, 1e0])
        ax12.set_yticklabels([r'$10^{-1}$', r'$10^{0}$'])

        ax12.tick_params(axis = 'both', labelsize = 30)

        ax12.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

        ax12.text(0.05,  0.70,   r'\boldmath$F_2^p$',         transform = ax12.transAxes, size = 60)
        ax12.text(0.15,  0.70,   r'$(\, \times\, 2^{\, i})$', transform = ax12.transAxes, size = 40)

        ax12.xaxis.set_tick_params(which = 'major', length = 10)
        ax12.xaxis.set_tick_params(which = 'minor', length = 5)
        ax12.yaxis.set_tick_params(which = 'major', length = 10)
        ax12.yaxis.set_tick_params(which = 'minor', length = 5)

        ax12.yaxis.set_ticks_position('both')
        ax12.xaxis.set_ticks_position('both')
        ax12.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax12.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

        minorLocator = MultipleLocator(0.2)
        ax12.xaxis.set_minor_locator(minorLocator)

        handles = [hand['JLab p']]
        label1  = r'\textbf{\textrm{Hall C}}'
        labels  = [label1,label2]
        ax12.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #################################
    #--Plot NC cross-section from HERA at sqrt(s) = 318
    #################################

    if 10026 and 10030 in data:
        hera_10026 = data[10026] ## HERA e+p \sqrt(s)=318 NC
        hera_10030 = data[10030] ## HERA e-p \sqrt(s)=318 NC

        nbins = 25 

        DATA = {}
        DATA['HERA 10026'] = pd.DataFrame(hera_10026) 
        DATA['HERA 10030'] = pd.DataFrame(hera_10030) 

        theory = {}
        PLOT = {}
        for exp in DATA:
            query = get_xbins(DATA[exp],'HERA NC')
            PLOT[exp] = get_plot(query)

            theory[exp] = get_theory(PLOT[exp],nbins,loop=False)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                Q2    = PLOT[exp]['Q2'][key]
                val   = PLOT[exp]['value'][key]*2.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax21.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            for key in theory[exp]['value']:
                Q2   = theory[exp]['Q2'][key]
                mean = theory[exp]['value'][key]*2.0**float(key)
                std  = theory[exp]['std'][key]  *2.0**float(key)
                down = mean - std
                up   = mean + std
                thy_plot ,= ax21.plot(Q2,mean,linestyle='solid',color='black')
                thy_band  = ax21.fill_between(Q2,down,up,color='gold', alpha=1.0)


        ax21.text(3.4,1.65e7, r'$x=2.8\cdot 10^{-5}\, (i=24)$', fontsize = 22)
        ax21.text(5.5,8.65e6, r'$x=4.6\cdot 10^{-5}$',          fontsize = 22)
        ax21.text(8.0,4.6e6,  r'$x=7.3\cdot 10^{-5}$',          fontsize = 22)
        ax21.text(15.0,2.7e6, r'$x=1.2\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(22.0,1.4e6, r'$x=1.9\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(35.0,7.5e5, r'$x=3.1\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(45.0,3.9e5, r'$x=4.4\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(60.0,1.8e5, r'$x=6.2\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(80.0,8.5e4, r'$x=7.7\cdot 10^{-4}$',          fontsize = 22)
        ax21.text(86.0,4e4,   r'$x=9.1 \cdot 10^{-4}$',         fontsize = 22)
        ax21.text(160.0,2e4,  r'$x=0.0014$',                    fontsize = 22)
        ax21.text(260.0,1e4,  r'$x=0.0024$',                    fontsize = 22)
        ax21.text(440.0,4500, r'$x=0.0035$',                    fontsize = 22)
        ax21.text(740.0,2000, r'$x=0.0056$',                    fontsize = 22)
        ax21.text(1.5e3,970,  r'$x=0.0083$',                    fontsize = 22)
        ax21.text(1.7e3,410,  r'$x=0.013$',                     fontsize = 22)
        ax21.text(2.7e3,190,  r'$x=0.021$',                     fontsize = 22)
        ax21.text(3.9e3,77,   r'$x=0.032$',                     fontsize = 22)
        ax21.text(6.5e3,35,   r'$x=0.052$',                     fontsize = 22)
        ax21.text(1.15e4,19,  r'$x=0.084$',                     fontsize = 22)
        ax21.text(12.0,5.5,   r'$x=0.13$',                      fontsize = 22)
        ax21.text(18.0,2.25,  r'$x=0.18$',                      fontsize = 22)
        ax21.text(36.0,0.97,  r'$x=0.25$',                      fontsize = 22)
        ax21.text(43.0,0.25,  r'$x=0.40$',                      fontsize = 22)
        ax21.text(30.0,0.02,  r'$x=0.65 \, (i=0)$',             fontsize = 22)

        ax21.semilogy()
        ax21.semilogx()

        ax21.set_xlim(1.0,   2e5)
        ax21.set_ylim(0.004, 4e7)

        ax21.set_xticks([1.0, 10.0, 1e2, 1e3, 1e4, 1e5])
        ax21.set_xticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
        locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
        ax21.yaxis.set_minor_locator(locmin)
        ax21.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax21.yaxis.set_tick_params(which = 'minor', length = 5)
        ax21.set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
        ax21.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$'])

        ax21.tick_params(axis = 'both', labelsize = 30)

        ax21.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

        #ax21.text(1.2, 5e-2,  r'\boldmath$\sigma_r^{p,NC}$', size = 60)
        #ax21.text(13.0,7e-2, r'$(\, \times\, 2^{\, i})$',   size = 40)
        ax21.text(0.05, 0.10,  r'\boldmath$\sigma_r^{p \textrm{(NC)}}$', transform = ax21.transAxes, size = 60)
        ax21.text(0.05, 0.04, r'$(\, \times\, 2^{\, i})$',               transform = ax21.transAxes, size = 40)

        ax21.text(2.0e3, 3.0e4, r'$\sqrt{s}=318 \,\rm{GeV}$', fontsize = 40)

        ax21.xaxis.set_tick_params(which = 'major', length = 10)
        ax21.xaxis.set_tick_params(which = 'minor', length = 5)
        ax21.yaxis.set_tick_params(which = 'major', length = 10)
        ax21.yaxis.set_tick_params(which = 'minor', length = 5)

        ax21.yaxis.set_ticks_position('both')
        ax21.xaxis.set_ticks_position('both')
        ax21.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax21.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')


        handles = [hand['HERA 10026'],hand['HERA 10030']]
        label1  = r'\textbf{\textrm{HERA NC}} \boldmath$e^+p$'
        label2  = r'\textbf{\textrm{HERA NC}} \boldmath$e^-p$'
        labels  = [label1,label2]
        ax21.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #################################
    #--Plot NC cross-section from HERA at sqrt(s) other than 318
    #################################

    if 10027 and 10028 and 10029 in data:
        hera_10027 = data[10027] ## HERA e+p \sqrt(s)=300 NC
        hera_10028 = data[10028] ## HERA e+p \sqrt(s)=251 NC
        hera_10029 = data[10029] ## HERA e+p \sqrt(s)=225 NC

        nbins = 30

        DATA = {}
        DATA['HERA 10027'] = pd.DataFrame(hera_10027) 
        DATA['HERA 10028'] = pd.DataFrame(hera_10028) 
        DATA['HERA 10029'] = pd.DataFrame(hera_10029)
 
        theory = {}
        PLOT = {}
        for exp in DATA:
            query = get_xbins(DATA[exp],'HERA other')
            PLOT[exp] = get_plot(query)

            theory[exp] = get_theory(PLOT[exp],nbins,loop=False)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                Q2    = PLOT[exp]['Q2'][key]
                val   = PLOT[exp]['value'][key]*2.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax22.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            for key in theory[exp]['value']:
                Q2   = theory[exp]['Q2'][key]
                mean = theory[exp]['value'][key]*2.0**float(key)
                std  = theory[exp]['std'][key]  *2.0**float(key)
                down = mean - std
                up   = mean + std
                thy_plot ,= ax22.plot(Q2,mean,linestyle='solid',color=color)
                thy_band  = ax22.fill_between(Q2,down,up,color='gold', alpha = 1.0)

        ax22.text(2.5,    4.5e8,  r'$x=3.3 \cdot 10^{-5}\, (i=29)$', fontsize = 20)
        ax22.text(5.7,    2.5e8,  r'$x=5.9 \cdot 10^{-5}$',          fontsize = 20)
        ax22.text(19.0,   1.7e8,  r'$x=1.6 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(35.0,   8e7,    r'$x=3.1 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(35.0,   4.4e7,  r'$x=4.5 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(45.0,   2.1e7,  r'$x=6.2 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(57.0,   1.1e7,  r'$x=7.8 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(0.5e2,  5e6,    r'$x=9.1 \cdot 10^{-4}$',          fontsize = 20)
        ax22.text(0.9e2,  2.7e6,  r'$x=0.0011$',                     fontsize = 20)
        ax22.text(1e2,    1.15e6, r'$x=0.0015$',                     fontsize = 20)
        ax22.text(1.7e2,  6e5,    r'$x=0.0022$',                     fontsize = 20)
        ax22.text(2.5e2,  2.9e5,  r'$x=0.0034$',                     fontsize = 20)
        ax22.text(3.7e2,  1.2e5,  r'$x=0.0046$',                     fontsize = 20)
        ax22.text(0.45e3, 5.6e4,  r'$x=0.0064$',                     fontsize = 20)
        ax22.text(0.8e3,  2.8e4,  r'$x=0.0085$',                     fontsize = 20)
        ax22.text(0.6e3,  1.2e4,  r'$x=0.0094$',                     fontsize = 20)
        ax22.text(0.8e3,  6200.0, r'$x=0.011$',                      fontsize = 20)
        ax22.text(1.5e3,  3200.0, r'$x=0.0138$',                     fontsize = 20)
        ax22.text(1e3,    1400.0, r'$x=0.0156$',                     fontsize = 20)
        ax22.text(1e3,    700.0,  r'$x=0.0183$',                     fontsize = 20)
        ax22.text(1.9e3,  320.0,  r'$x=0.020$',                      fontsize = 20)
        ax22.text(1.1e3,  150.0,  r'$x=0.025$',                      fontsize = 20)
        ax22.text(2.5e3,  67.0,   r'$x=0.032$',                      fontsize = 20)
        ax22.text(60.0,   28.0,   r'$x=0.05$',                       fontsize = 20)
        ax22.text(60.0,   12.0,   r'$x=0.08$',                       fontsize = 20)
        ax22.text(60.0,   5.3,    r'$x=0.13$',                       fontsize = 20)
        ax22.text(60.0,   2.15,   r'$x=0.18$',                       fontsize = 20)
        ax22.text(160.0,  0.85,   r'$x=0.25$',                       fontsize = 20)
        ax22.text(60.0,   0.3,    r'$x=0.40$',                       fontsize = 20)
        ax22.text(70.0,   0.022,  r'$x=0.65 \, (i=0)$',              fontsize = 20)

        ax22.text(4e3,   28.0, r'$x=0.05$',                          fontsize = 20, color = 'darkgreen')
        ax22.text(6e3,   11.0, r'$x=0.08$',                          fontsize = 20, color = 'darkgreen')
        ax22.text(1.6e3, 4.6,  r'$x=0.13$',                          fontsize = 20, color = 'darkgreen')
        ax22.text(2.3e3, 2.0,  r'$x=0.18$',                          fontsize = 20, color = 'darkgreen')
        ax22.text(2.8e3, 0.8,  r'$x=0.25$',                          fontsize = 20, color = 'darkgreen')
        ax22.text(3.0e3, 0.2,  r'$x=0.40\, (i=1)$',                  fontsize = 20, color = 'darkgreen')

        ax22.semilogy()
        ax22.semilogx()
        ax22.set_xlim(1.0,  4e4)
        ax22.set_ylim(1e-2, 2e9)

        ax22.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

        #ax22.text(1.50 , 3.0e-2, r'\boldmath$\sigma_r^{p,NC}$', size = 60)
        #ax22.text(13.0 , 5.0e-2, r'$(\, \times\, 2^{\, i})$',   size = 40)
        ax22.text(0.05, 0.10,  r'\boldmath$\sigma_r^{p \textrm{(NC)}}$', transform = ax22.transAxes, size = 60)
        ax22.text(0.05, 0.04, r'$(\, \times\, 2^{\, i})$',               transform = ax22.transAxes, size = 40)

        ax22.set_xticks([1.0, 10.0, 1e2, 1e3, 1e4])
        ax22.set_xticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
        locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
        ax22.yaxis.set_minor_locator(locmin)
        ax22.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax22.yaxis.set_tick_params(which = 'minor', length = 3)
        ax22.set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
        ax22.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$', r'$10^8$', r'$10^9$'])

        ax22.yaxis.set_ticks_position('both')
        ax22.xaxis.set_ticks_position('both')
        ax22.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax22.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

        ax22.xaxis.set_tick_params(which = 'major', length = 10)
        ax22.xaxis.set_tick_params(which = 'minor', length = 5)
        ax22.yaxis.set_tick_params(which = 'major', length = 10)
        ax22.yaxis.set_tick_params(which = 'minor', length = 5)

        ax22.yaxis.set_ticks_position('both')
        ax22.xaxis.set_ticks_position('both')
        ax22.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax22.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

        ax22.tick_params(axis='both',labelsize=30)

        handles = [hand['HERA 10027'],hand['HERA 10028'],hand['HERA 10029']]
        label1  = r'\textbf{\textrm{HERA}}' + r'$\sqrt{s}=300$'# + ' ' +  r'\textrm{GeV}'
        label2  = r'\textbf{\textrm{HERA}}' + r'$\sqrt{s}=251$'# + ' ' +  r'\textrm{GeV}'
        label3  = r'\textbf{\textrm{HERA}}' + r'$\sqrt{s}=225$'# + ' ' +  r'\textrm{GeV}'
        labels  = [label1,label2,label3]
        ax22.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-proton.png' % (wdir))
    print('Saving figure to %s/gallery/dis-proton.png'%(wdir))
    py.close()

def plot_deuteron(wdir, data):

    if 10017 not in data: return
    if 10011 not in data: return
    if 10021 not in data: return

    nrows, ncols = 2, 2
    py.figure(figsize = (ncols * 12, nrows * 14))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)
    ax21 = py.subplot(nrows, ncols, 3)
    ax22 = py.subplot(nrows, ncols, 4)

    #################################
    #--Plot F2d from BCDMS and SLAC
    #################################

    nbins = 19

    bcdms = data[10017] ## BCDMS d
    slac  = data[10011] ## SLAC d

    DATA = {}
    DATA['SLAC']  = pd.DataFrame(slac)
    DATA['BCDMS'] = pd.DataFrame(bcdms)

    PLOT = {}
    for exp in DATA:
        query = get_xbins(DATA[exp],'d')
        PLOT[exp] = get_plot(query)

    theory = get_theory(PLOT,nbins)

    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            Q2    = PLOT[exp]['Q2'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax11.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

    #--plot theory interpolated between data points
    for key in theory['value']:
        Q2   = theory['Q2'][key]
        mean = theory['value'][key]*2.0**float(key)
        std  = theory['std'][key]  *2.0**float(key)
        down = mean - std
        up   = mean + std
        thy_plot ,= ax11.plot(Q2,mean,linestyle='solid',color='black')
        thy_band  = ax11.fill_between(Q2,down,up,color='gold',alpha=1.0)


    ax11.text(21.0,  9.0e4,  r'$x=0.07 \,(i=15)$', fontsize = 26)
    ax11.text(3.0,   4.0e4,  r'$x=0.09$',          fontsize = 26)
    ax11.text(45.0,  2.0e4,  r'$x=0.10$',          fontsize = 26)
    ax11.text(70.0,  1.0e4,  r'$x=0.14$',          fontsize = 26)
    ax11.text(75.0,  4.5e3,  r'$x=0.18$',          fontsize = 26)
    ax11.text(100.0, 2.0e3,  r'$x=0.23$',          fontsize = 26)
    ax11.text(140.0, 9.0e2,  r'$x=0.28$',          fontsize = 26)
    ax11.text(170.0, 3.0e2,  r'$x=0.35$',          fontsize = 26)
    ax11.text(190.0, 1.0e2,  r'$x=0.43$',          fontsize = 26)
    ax11.text(19.0,  58.0,   r'$x=0.46$',          fontsize = 26)
    ax11.text(17.0,  20.0,   r'$x=0.50$',          fontsize = 26)
    ax11.text(255.0, 4.00,   r'$x=0.53$',          fontsize = 26)
    ax11.text(22.0,  3.0,    r'$x=0.56$',          fontsize = 26)
    ax11.text(23.0,  1.50,   r'$x=0.59$',          fontsize = 26)
    ax11.text(25.0,  0.60,   r'$x=0.61$',          fontsize = 26)
    ax11.text(250.0, 0.14,   r'$x=0.65$',          fontsize = 26)
    ax11.text(250.0, 0.020,  r'$x=0.75$',          fontsize = 26)
    ax11.text(35.0,  0.005,  r'$x=0.85$',          fontsize = 26)
    ax11.text(35.0,  0.0015, r'$x=0.88 \,(i=0)$',  fontsize = 26)

    ax11.semilogy()
    ax11.semilogx()

    ax11.set_xlim(1.5, 1e3)
    ax11.set_ylim(0.001, 2e5)

    ax11.tick_params(axis = 'both', labelsize = 30)

    ax11.yaxis.set_tick_params(which = 'major', length = 10)
    ax11.yaxis.set_tick_params(which = 'minor', length = 5)
    ax11.set_yticks([0.01, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4,1e5])
    ax11.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$',r'$10^5$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax11.yaxis.set_minor_locator(locmin)
    ax11.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax11.yaxis.set_tick_params(which = 'minor', length = 5)

    ax11.xaxis.set_tick_params(which = 'major', length = 10)
    ax11.xaxis.set_tick_params(which = 'minor', length = 5)
    ax11.set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2])
    xticklabels=[2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2]
    ax11.set_xticklabels(['$%0.0f$'%x for x in xticklabels])

    ax11.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

    ax11.text(0.05,0.05,r'\boldmath$F_2^d$',      transform = ax11.transAxes, size = 60)
    ax11.text(0.15,0.05,r'$(\times\, 2^{\, i})$', transform = ax11.transAxes, size = 40)

    ax11.yaxis.set_ticks_position('both')
    ax11.xaxis.set_ticks_position('both')
    ax11.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
    ax11.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

    handles = [hand['BCDMS'],hand['SLAC'],(thy_band,thy_plot)]
    label1  = r'\textbf{\textrm{BCDMS}}'
    label2  = r'\textbf{\textrm{SLAC}}'
    label3  = r'\textbf{\textrm{JAM}}'
    labels  = [label1,label2,label3]
    ax11.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #################################
    #--Plot F2d from Hall C
    #################################

    if 10002 in data:
        jlab_d = data[10002] ## JLab d 


        E  = 5.5 #--Elab
        Mp = 0.939

        theta = lambda X,Q2,M: np.arcsin(np.sqrt(Q2/(4*E**2 - 2*E*Q2/M/X)))*2*180/np.pi

        data[10002]['theta'] = theta(data[10002]['X'],data[10002]['Q2'],Mp)

        nbins = 6 

        DATA = {}
        DATA['JLab d'] = pd.DataFrame(jlab_d) 

        theory = {}
        PLOT = {}
        for exp in DATA:
            query = get_thetabins(DATA[exp],'JLab')
            PLOT[exp] = get_plot(query)

            theory[exp] = get_theory(PLOT[exp],nbins,loop=False)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                Q2    = PLOT[exp]['Q2'][key]
                val   = PLOT[exp]['value'][key]*2.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax12.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            for key in theory[exp]['value']:
                Q2   = theory[exp]['Q2'][key]
                mean = theory[exp]['value'][key]*2.0**float(key)
                std  = theory[exp]['std'][key]  *2.0**float(key)
                down = mean - std
                up   = mean + std
                thy_plot ,= ax12.plot(Q2,mean,linestyle='solid',color='black')
                thy_band  = ax12.fill_between(Q2,down,up,color='gold', alpha = 1.0)


        ax12.text(6.55,0.80,   r'$\theta=70^{\circ} \, (i=5)$',   fontsize = 30)
        ax12.text(6.2, 0.45,   r'$\theta=60^{\circ} $',           fontsize = 30)
        ax12.text(5.9, 0.24,   r'$\theta=55^{\circ} $',           fontsize = 30)
        ax12.text(5.2, 0.16,   r'$\theta=45^{\circ} $',           fontsize = 30)
        ax12.text(4.9, 0.090,  r'$\theta=41^{\circ} $',           fontsize = 30)
        ax12.text(4.6, 0.050,  r'$\theta=38^{\circ}  \, (i=0)$',  fontsize = 30)

        ax12.semilogy()

        ax12.set_xlim(3.5,   7.5)
        ax12.set_ylim(0.04,  3.5)

        ax12.set_xticks([4, 5, 6, 7])
        ax12.set_xticklabels([r'$4$', r'$5$', r'$6$', r'$7$'])
        locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
        ax12.yaxis.set_minor_locator(locmin)
        ax12.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax12.yaxis.set_tick_params(which = 'minor', length = 5)
        ax12.set_yticks([1e-1, 1e0])
        ax12.set_yticklabels([r'$10^{-1}$', r'$10^{0}$'])

        ax12.tick_params(axis = 'both', labelsize = 30)

        ax12.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

        ax12.text(0.05,  0.70,   r'\boldmath$F_2^d$',         transform = ax12.transAxes, size = 60)
        ax12.text(0.15,  0.70,   r'$(\, \times\, 2^{\, i})$', transform = ax12.transAxes, size = 40)

        ax12.xaxis.set_tick_params(which = 'major', length = 10)
        ax12.xaxis.set_tick_params(which = 'minor', length = 5)
        ax12.yaxis.set_tick_params(which = 'major', length = 10)
        ax12.yaxis.set_tick_params(which = 'minor', length = 5)

        ax12.yaxis.set_ticks_position('both')
        ax12.xaxis.set_ticks_position('both')
        ax12.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax12.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

        minorLocator = MultipleLocator(0.2)
        ax12.xaxis.set_minor_locator(minorLocator)

        handles = [hand['JLab d']]
        label1  = r'\textbf{\textrm{Hall C}}'
        labels  = [label1]
        ax12.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #################################
    #--Plot F2d/F2p from NMC
    #################################

    nbins = 18
    nmc   = data[10021] ## NMC d/p

    DATA = {}
    DATA['NMC']   = pd.DataFrame(nmc)

    PLOT = {}
    for exp in DATA:
        query = get_xbins(DATA[exp],'d/p')
        PLOT[exp] = get_plot(query)

    theory = get_theory(PLOT,nbins)

    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            Q2    = PLOT[exp]['Q2'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax21.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

    #--plot theory interpolated between data points
    for key in theory['value']:
        Q2   = theory['Q2'][key]
        mean = theory['value'][key]*2.0**float(key)
        std  = theory['std'][key]  *2.0**float(key)
        down = mean - std
        up   = mean + std
        thy_plot ,= ax21.plot(Q2,mean,linestyle='solid',color='black')
        thy_band  = ax21.fill_between(Q2,down,up,color='gold',alpha=1.0)

    ax21.text(2.0,  1.2e5,  r'$x=0.005 \, (i=17)$', fontsize = 26)
    ax21.text(4.0,  0.63e5, r'$x=0.008$',           fontsize = 26)
    ax21.text(6.3,  3.1e4,  r'$x=0.013$',           fontsize = 26)
    ax21.text(8.2,  1.5e4,  r'$x=0.018$',           fontsize = 26)
    ax21.text(10.5, 7.8e3,  r'$x=0.025$',           fontsize = 26)
    ax21.text(16.5, 3.6e3,  r'$x=0.035$',           fontsize = 26)
    ax21.text(22.0, 1.9e3,  r'$x=0.05$',            fontsize = 26)
    ax21.text(30.0, 0.9e3,  r'$x=0.07$',            fontsize = 26)
    ax21.text(40,   4.7e2,  r'$x=0.09$',            fontsize = 26)
    ax21.text(52.0, 2.2e2,  r'$x=0.11$',            fontsize = 26)
    ax21.text(56.0, 1.1e2,  r'$x=0.14$',            fontsize = 26)
    ax21.text(74.0, 5.4e1,  r'$x=0.18$',            fontsize = 26)
    ax21.text(75.0, 2.5e1,  r'$x=0.23$',            fontsize = 26)
    ax21.text(108,  12.2,   r'$x=0.28$',            fontsize = 26)
    ax21.text(110,  6.0,    r'$x=0.35$',            fontsize = 26)
    ax21.text(110,  2.9,    r'$x=0.45$',            fontsize = 26)
    ax21.text(110,  1.35,   r'$x=0.55$',            fontsize = 26)
    ax21.text(110,  0.7,    r'$x=0.68 \,(i=0)$',    fontsize = 26)

    ax21.semilogy()
    ax21.semilogx()

    ax21.set_xlim(1.5, 3.5e2)
    ax21.set_ylim(0.5, 2e5)

    ax21.tick_params(axis = 'both', labelsize = 30)

    ax21.yaxis.set_tick_params(which = 'major', length = 10)
    ax21.yaxis.set_tick_params(which = 'minor', length = 5)
    ax21.set_yticks([1.0, 10.0, 1e2, 1e3, 1e4, 1e5])
    ax21.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax21.yaxis.set_minor_locator(locmin)
    ax21.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax21.yaxis.set_tick_params(which = 'minor', length = 5)

    ax21.xaxis.set_tick_params(which = 'major', length = 10)
    ax21.xaxis.set_tick_params(which = 'minor', length = 5)
    ax21.set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2])
    xticklabels=[2e0, 5e0, 1e1, 2e1, 5e1, 1e2]
    ax21.set_xticklabels(['$%0.0f$'%x for x in xticklabels])

    ax21.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

    ax21.text(0.60, 0.75,  r'\boldmath$F_2^d/F_2^p$',  transform = ax21.transAxes, size = 60)
    ax21.text(0.85, 0.75,  r'$(\, \times\, 2^{\, i})$',transform = ax21.transAxes, size = 40)

    ax21.yaxis.set_ticks_position('both')
    ax21.xaxis.set_ticks_position('both')
    ax21.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
    ax21.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

    handles = [hand['NMC']]
    label1  = r'\textbf{\textrm{NMC}}'
    labels  = [label1]
    ax21.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #################################
    #--Plot F2n/F2d from BONuS
    #################################

    if 10033 in data:
        bonus  = data[10033] ## BONuS n/d 

        nbins = 6 

        DATA = {}
        DATA['BONuS'] = pd.DataFrame(bonus) 

        PLOT   = {}
        theory = {}
        for exp in DATA:
            query = get_Q2bins(DATA[exp],'BONuS')
            PLOT[exp] = get_plot(query)

            theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)

        hand = {}
        #--plot data points
        for exp in PLOT:
            for key in range(nbins):
                X     = PLOT[exp]['X'][key]
                val   = PLOT[exp]['value'][key]*4.0**float(key)
                alpha = PLOT[exp]['alpha'][key]*4.0**float(key)
                color,marker,ms = get_details(exp)
                hand[exp] = ax22.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            for key in theory[exp]['value']:
                X    = theory[exp]['X'][key]
                mean = theory[exp]['value'][key]*4.0**float(key)
                std  = theory[exp]['std'][key]  *4.0**float(key)
                down = mean - std
                up   = mean + std
                thy_plot ,= ax22.plot(X,mean,linestyle='solid',color=color)
                thy_band  = ax22.fill_between(X,down,up,color='gold', alpha = 1.0)

        ax22.text(0.47, 0.05,   r'$Q^2 = 4.0 $'+'$\, (i=5)$',               transform = ax22.transAxes, fontsize = 30)
        ax22.text(0.55, 0.20,   r'$Q^2 = 3.4 $',                            transform = ax22.transAxes, fontsize = 30)
        ax22.text(0.62, 0.35,   r'$Q^2 = 2.9 $',                            transform = ax22.transAxes, fontsize = 30)
        ax22.text(0.74, 0.50,   r'$Q^2 = 2.4 $',                            transform = ax22.transAxes, fontsize = 30)
        ax22.text(0.78, 0.64,   r'$Q^2 = 2.0 $',                            transform = ax22.transAxes, fontsize = 30)
        ax22.text(0.30, 0.80,   r'$Q^2 = 1.7 \: \rm{GeV^2}$'+'$\, (i=0)$',  transform = ax22.transAxes, fontsize = 30)

        ax22.semilogy()

        ax22.set_xlim(0.2,   0.75)
        ax22.set_ylim(0.5,   4e3)

        ax22.set_xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
        ax22.set_xticklabels([r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$',r'$0.6$',r'$0.7$'])
        ax22.set_yticks([1e0, 1e1, 1e2, 1e3])
        ax22.set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$'])

        ax22.tick_params(axis = 'both', labelsize = 30)

        ax22.set_xlabel(r'\boldmath$x$',size=40)

        ax22.text(0.05, 0.90,   r'\boldmath$F_2^n/F_2^d$',   transform = ax22.transAxes, size = 60)
        ax22.text(0.30, 0.90,   r'$(\, \times\, 4^{\, i})$', transform = ax22.transAxes, size = 40)

        ax22.xaxis.set_tick_params(which = 'major', length = 10)
        ax22.xaxis.set_tick_params(which = 'minor', length = 5)
        ax22.yaxis.set_tick_params(which = 'major', length = 10)
        ax22.yaxis.set_tick_params(which = 'minor', length = 5)

        ax22.yaxis.set_ticks_position('both')
        ax22.xaxis.set_ticks_position('both')
        ax22.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
        ax22.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

        minorLocator = MultipleLocator(0.02)
        ax22.xaxis.set_minor_locator(minorLocator)

        handles = [hand['BONuS']]
        label1  = r'\textbf{\textrm{BONuS}}'
        labels  = [label1]
        #ax22.legend(handles,labels,loc=(0.02,0.6),  fontsize = 35, frameon = 0, handletextpad = 0.1, handlelength = 1.0)
        ax22.legend(handles,labels,loc='upper right',  fontsize = 35, frameon = 0, handletextpad = 0.1, handlelength = 1.0)

    py.tight_layout()
    py.savefig('%s/gallery/dis-deuteron.png' % (wdir))
    print('Saving figure to %s/gallery/dis-deuteron.png'%(wdir))
    py.close()

def plot_CC(wdir, data):

    if 10031 not in data: return
    if 10032 not in data: return

    nrows, ncols = 1, 1
    fig = py.figure(figsize = (ncols * 12, nrows * 14))
    ax = {}
    ax = fig.add_subplot(nrows, ncols, 1)

    hera_10031 = data[10031] ## HERA e+p \sqrt(s)=318 CC
    hera_10032 = data[10032] ## HERA e-p \sqrt(s)=318 CC

    nbins = 8 

    DATA = {}
    DATA['HERA 10031'] = pd.DataFrame(hera_10031) 
    DATA['HERA 10032'] = pd.DataFrame(hera_10032) 

    PLOT   = {}
    theory = {}
    for exp in DATA:
        query = get_xbins(DATA[exp],'HERA CC')
        PLOT[exp] = get_plot(query)

        theory[exp] = get_theory(PLOT[exp],nbins,loop=False)

    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            Q2    = PLOT[exp]['Q2'][key]
            val   = PLOT[exp]['value'][key]*5.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*5.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(Q2,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            Q2   = theory[exp]['Q2'][key]
            mean = theory[exp]['value'][key]*5.0**float(key)
            std  = theory[exp]['std'][key]  *5.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(Q2,mean,linestyle='solid',color='black')
            thy_band  = ax.fill_between(Q2,down,up,color='gold', alpha=1.0)


    ax.text(370.0,  1.4e5, r'$x=0.008 \, (i=7)$', fontsize = 28, color = 'darkgreen')
    ax.text(1200.0, 1.3e4, r'$x=0.013$',          fontsize = 28, color = 'darkgreen')
    ax.text(3700.0, 2e3,   r'$x=0.032$',          fontsize = 28, color = 'darkgreen')
    ax.text(6000.0, 400,   r'$x=0.08$',           fontsize = 28, color = 'darkgreen')
    ax.text(9800.0, 70,    r'$x=0.13$',           fontsize = 28, color = 'darkgreen')
    ax.text(400.0,  11.5,  r'$x=0.25$',           fontsize = 28, color = 'darkgreen')
    ax.text(600.0,  1.4,   r'$x=0.40$',           fontsize = 28, color = 'darkgreen')
    ax.text(1300.0, 0.035, r'$x=0.65 \, (i=0) $', fontsize = 28, color = 'darkgreen')

    ax.text(550.0,  0.6e5, r'$x=0.008 \, (i=7)$', fontsize = 28, color = 'blue')
    ax.text(1200.0, 0.8e4, r'$x=0.013$',          fontsize = 28, color = 'blue')
    ax.text(3700.0, 1000,  r'$x=0.032$',          fontsize = 28, color = 'blue')
    ax.text(4300.0, 180,   r'$x=0.08$',           fontsize = 28, color = 'blue')
    ax.text(0.8e4,  23,    r'$x=0.13$',           fontsize = 28, color = 'blue')
    ax.text(400.0,  5.2,   r'$x=0.25$',           fontsize = 28, color = 'blue')
    ax.text(370.0,  0.3,   r'$x=0.40 \, (i=1)$',  fontsize = 28, color = 'blue')

    ax.semilogy()
    ax.semilogx()

    ax.set_xlim(1.5e2,   4e4)
    ax.set_ylim(1.0e-2,  8e5)

    ax.set_xticks([2e2, 1e3, 1e4, 4e4])
    ax.set_xticklabels([r'$200$', r'$10^3$', r'$10^4$', r'$4\cdot10^4$'])
    ax.set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])
    ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])

    ax.tick_params(axis = 'both', labelsize = 30)

    ax.set_xlabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)

    #ax.text(1.8e2, 2.0e-2, r'\boldmath$\sigma_r^{p,CC}$', size = 60)
    #ax.text(5.6e2, 2.6e-2,  r'$(\, \times\, 5^{\, i})$',   size = 40)
    ax.text(2.0e2, 4.0e-2,  r'\boldmath$\sigma_r^{p \textrm{(CC)}}$', size = 60)
    ax.text(2.0e2, 2.0e-2,  r'$(\, \times\, 5^{\, i})$',   size = 40)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction = 'in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction = 'in')

    #ax.text(5.0e3,  2.0e4, r'$\sqrt{s}=318 \,\rm{GeV}$', fontsize = 40)

    handles = [hand['HERA 10031'],hand['HERA 10032'],(thy_band,thy_plot)]
    label1  = r'\textbf{\textrm{HERA CC}} \boldmath$e^+p$'
    label2  = r'\textbf{\textrm{HERA CC}} \boldmath$e^-p$'
    label3  = r'\textbf{\textrm{JAM}}'
    labels  = [label1,label2,label3]
    ax.legend(handles,labels,loc='upper right', fontsize = 35, frameon = 0, handletextpad = 0.3, handlelength = 1.0)
    py.tight_layout()

    py.savefig('%s/gallery/dis-CC.png' % (wdir))
    print('Saving figure to %s/gallery/dis-CC.png'%(wdir))
    py.close()

def plot_marathon(wdir, data):

    if 10050 not in data: return
    if 10051 not in data: return

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 5, nrows * 6))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #################################
    #--Plot F2d/F2p from MARATHON
    #################################

    dp   = data[10050] ## MARATHON d/p

    DATA = {}
    DATA['MARATHON dp']   = pd.DataFrame(dp)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'MARATHON')
        nbins = len(query)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            Q2    = PLOT[exp]['Q2'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax11.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax11.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax11.fill_between(X,down,up,color='gold',alpha=1.0)


    ax11.set_xlim(0.15, 0.45)
    ax11.set_ylim(0.76, 0.88)

    ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax11.yaxis.set_tick_params(which = 'major', length = 5)
    ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.xaxis.set_tick_params(which = 'major', length = 5)
    ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.set_xlabel(r'\boldmath$x$', size=30)
    ax11.xaxis.set_label_coords(0.95,0.00)

    ax11.text(0.05, 0.05, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 25)

    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.10)
    ax11.xaxis.set_minor_locator(minorLocator)
    ax11.xaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(0.01)
    majorLocator = MultipleLocator(0.05)
    ax11.yaxis.set_minor_locator(minorLocator)
    ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.2,0.4,0.6,0.8])

    handles, labels = [],[]
    handles.append(hand['MARATHON dp'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\boldmath$F_2^D/F_2^p$')
    labels.append(r'\textbf{\textrm{JAM}}')
    ax11.legend(handles,labels,loc='upper right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #################################
    #--Plot F2h/F2t from MARATHON
    #################################

    ht   = data[10051] ## MARATHON h/t

    DATA = {}
    DATA['MARATHON ht']  = pd.DataFrame(ht)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'MARATHON')
        nbins = len(query)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax12.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax12.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax12.fill_between(X,down,up,color='gold',alpha=1.0)


    ax12.set_xlim(0.10, 0.95)
    ax12.set_ylim(1.05, 1.32)

    ax12.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax12.yaxis.set_tick_params(which = 'major', length = 5)
    ax12.yaxis.set_tick_params(which = 'minor', length = 2.5)

    ax12.xaxis.set_tick_params(which = 'major', length = 5)
    ax12.xaxis.set_tick_params(which = 'minor', length = 2.5)

    ax12.set_xlabel(r'\boldmath$x$', size=30)
    ax12.xaxis.set_label_coords(0.95,0.00)

    ax12.text(0.50, 0.05, r'$Q^2 = 14 \cdot x ~ {\rm GeV^2}$', transform = ax12.transAxes, size = 20)

    minorLocator = MultipleLocator(0.04)
    majorLocator = MultipleLocator(0.2)
    ax12.xaxis.set_minor_locator(minorLocator)
    ax12.xaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.1)
    ax12.yaxis.set_minor_locator(minorLocator)
    ax12.yaxis.set_major_locator(majorLocator)
    ax12.set_xticks([0.2,0.4,0.6,0.8])

    handles, labels = [],[]
    handles.append(hand['MARATHON ht'])
    labels.append(r'\boldmath$F_2^{^3\rm{He}}/F_2^{^3\rm{H}}$')
    ax12.legend(handles,labels,loc='upper left', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.savefig('%s/gallery/dis-marathon.png' % (wdir))
    print('Saving figure to %s/gallery/dis-marathon.png'%(wdir))
    py.close()

def plot_marathon2(wdir, data):

    if 10052 not in data: return
    if 10053 not in data: return

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #################################
    #--Plot F2h/F2d from MARATHON
    #################################

    hd   = data[10052] ## MARATHON h/d

    DATA = {}
    DATA['MARATHON hd']   = pd.DataFrame(hd)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'MARATHON')
        nbins = len(query)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            Q2    = PLOT[exp]['Q2'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax11.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax11.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax11.fill_between(X,down,up,color='gold',alpha=1.0)


    ax11.set_xlim(0.10, 0.95)
    ax11.set_ylim(1.01, 1.16)

    ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax11.yaxis.set_tick_params(which = 'major', length = 5)
    ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.xaxis.set_tick_params(which = 'major', length = 5)
    ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.set_xlabel(r'\boldmath$x$', size=30)
    ax11.xaxis.set_label_coords(0.95,0.00)

    ax11.text(0.40, 0.05, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)

    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.10)
    ax11.xaxis.set_minor_locator(minorLocator)
    ax11.xaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(0.01)
    majorLocator = MultipleLocator(0.05)
    ax11.yaxis.set_minor_locator(minorLocator)
    ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.2,0.4,0.6,0.8])

    handles, labels = [],[]
    handles.append(hand['MARATHON hd'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\boldmath$F_2^{^3{\rm He}}/F_2^d$')
    labels.append(r'\textbf{\textrm{JAM}}')
    ax11.legend(handles,labels,loc='upper left', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #################################
    #--Plot F2t/F2d from MARATHON
    #################################

    td   = data[10053] ## MARATHON t/d

    DATA = {}
    DATA['MARATHON td']  = pd.DataFrame(td)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'MARATHON')
        nbins = len(query)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax12.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax12.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax12.fill_between(X,down,up,color='gold',alpha=1.0)


    ax12.set_xlim(0.10, 0.95)
    ax12.set_ylim(0.80, 1.00)

    ax12.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax12.yaxis.set_tick_params(which = 'major', length = 5)
    ax12.yaxis.set_tick_params(which = 'minor', length = 2.5)

    ax12.xaxis.set_tick_params(which = 'major', length = 5)
    ax12.xaxis.set_tick_params(which = 'minor', length = 2.5)

    ax12.set_xlabel(r'\boldmath$x$', size=30)
    ax12.xaxis.set_label_coords(0.95,0.00)

    #ax12.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax12.transAxes, size = 30)

    minorLocator = MultipleLocator(0.04)
    majorLocator = MultipleLocator(0.2)
    ax12.xaxis.set_minor_locator(minorLocator)
    ax12.xaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.1)
    ax12.yaxis.set_minor_locator(minorLocator)
    ax12.yaxis.set_major_locator(majorLocator)
    ax12.set_xticks([0.2,0.4,0.6,0.8])

    handles, labels = [],[]
    handles.append(hand['MARATHON td'])
    labels.append(r'\boldmath$F_2^{^3\rm{H}}/F_2^d$')
    ax12.legend(handles,labels,loc='upper right', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.savefig('%s/gallery/dis-marathon2.png' % (wdir))
    print('Saving figure to %s/gallery/dis-marathon2.png'%(wdir))
    py.close()

def plot_JLab_Hall_C_E03_103(wdir, data):

    if 10041 not in data: return

    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)

    #################################
    #--Plot sigma(h)/sigma(d) from JLab Hall C
    #################################

    nbins = 1

    hd   = data[10041] ## JLab Hall C h/d

    DATA = {}
    DATA['JLab hd']  = pd.DataFrame(hd)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'JLab')
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax11.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax11.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax11.fill_between(X,down,up,color='gold',alpha=1.0)


    ax11.set_xlim(0.28, 0.78)
    ax11.set_ylim(1.02, 1.09)

    ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax11.yaxis.set_tick_params(which = 'major', length = 5)
    ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.xaxis.set_tick_params(which = 'major', length = 5)
    ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    ax11.set_xlabel(r'\boldmath$x$', size=30)
    ax11.xaxis.set_label_coords(0.95,0.00)

    #ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.1)
    ax11.xaxis.set_minor_locator(minorLocator)
    ax11.xaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(0.004)
    majorLocator = MultipleLocator(0.02)
    ax11.yaxis.set_minor_locator(minorLocator)
    ax11.yaxis.set_major_locator(majorLocator)
    ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])
    ax11.set_yticks([1.04, 1.06, 1.08])

    handles, labels = [],[]
    handles.append(hand['JLab hd'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{JLab Hall C}}')
    labels.append(r'\textbf{\textrm{JAM}}')
    ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-JLab-Hall-C-E03-103.png' % (wdir))
    print('Saving figure to %s/gallery/dis-JLab-Hall-C-E03-103.png'%(wdir))
    py.close()

def plot_JLab_Hall_C_E12_10_002(wdir, data):

    if 10042 not in data: return

    nrows, ncols = 3, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax01 = py.subplot(nrows, ncols, 1)
    ax02 = py.subplot(nrows, ncols, 2)
    ax11 = py.subplot(nrows, ncols, 3)
    ax12 = py.subplot(nrows, ncols, 4)
    ax21 = py.subplot(nrows, ncols, 5)
    ax22 = py.subplot(nrows, ncols, 6)

    #################################
    #--Plot F2d/F2h from JLab Hall C
    #################################

    nbins = 5

    dp   = data[10042] ## JLab Hall C d/p

    DATA = {}
    DATA['JLab dp0']  = pd.DataFrame(dp).query("spec==0")
    DATA['JLab dp1']  = pd.DataFrame(dp).query("spec==1")

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_thetabins(DATA[exp],'HallC')
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if   exp=='JLab dp0' and key==0: ax=ax01
            elif exp=='JLab dp1' and key==0: ax=ax02
            elif exp=='JLab dp1' and key==1: ax=ax11
            elif exp=='JLab dp1' and key==2: ax=ax12
            elif exp=='JLab dp1' and key==3: ax=ax21
            elif exp=='JLab dp1' and key==4: ax=ax22
            else: continue
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]#*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]#*2.0**float(key)
            color,marker,ms = get_details(exp)
            Q2 = PLOT[exp]['Q2'][key]
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            if   exp=='JLab dp0' and key==0: ax=ax01
            elif exp=='JLab dp1' and key==0: ax=ax02
            elif exp=='JLab dp1' and key==1: ax=ax11
            elif exp=='JLab dp1' and key==2: ax=ax12
            elif exp=='JLab dp1' and key==3: ax=ax21
            elif exp=='JLab dp1' and key==4: ax=ax22
            else: continue
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]#*2.0**float(key)
            std  = theory[exp]['std'][key]  #*2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax.fill_between(X,down,up,color='gold',alpha=1.0)



    for ax in [ax01,ax02,ax11,ax12,ax21,ax22]:
        ax.set_xlim(0.20, 0.89)
        ax.set_ylim(0.69, 0.91)

        ax.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 30)

        ax.tick_params(axis='both', which = 'major', length = 8)
        ax.tick_params(axis='both', which = 'minor', length = 4)

        ax.set_xticks([0.4,0.6,0.8])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.set_yticks([0.7,0.8,0.9])
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    ax01.tick_params(labelbottom=False)
    ax02.tick_params(labelbottom=False)
    ax11.tick_params(labelbottom=False)
    ax12.tick_params(labelbottom=False)
    ax02.tick_params(labelleft=False)
    ax12.tick_params(labelleft=False)
    ax22.tick_params(labelleft=False)

    ax21.set_xlabel(r'\boldmath$x_{\rm bj}$', size=40)
    #ax21.xaxis.set_label_coords(0.95,0.00)
    ax22.set_xlabel(r'\boldmath$x_{\rm bj}$', size=40)
    #ax22.xaxis.set_label_coords(0.95,0.00)

    #ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    ax01.text(0.05, 0.08, r'\boldmath$F_2^{D}/F_2^{p}$', transform = ax01.transAxes, size = 60)

    ax01.text(0.05, 0.90, r'$4.5 < Q^2 < 7.3  ~ {\rm GeV}^2; \theta_c = 21^{\circ}$', transform = ax01.transAxes, size = 30)
    ax02.text(0.05, 0.90, r'$3.4 < Q^2 < 8.9  ~ {\rm GeV}^2; \theta_c = 21^{\circ}$', transform = ax02.transAxes, size = 30)
    ax11.text(0.05, 0.90, r'$4.6 < Q^2 < 8.9  ~ {\rm GeV}^2; \theta_c = 25^{\circ}$', transform = ax11.transAxes, size = 30)
    ax12.text(0.05, 0.90, r'$4.9 < Q^2 < 10.1 ~ {\rm GeV}^2; \theta_c = 29^{\circ}$', transform = ax12.transAxes, size = 30)
    ax21.text(0.05, 0.90, r'$5.2 < Q^2 < 11.1 ~ {\rm GeV}^2; \theta_c = 33^{\circ}$', transform = ax21.transAxes, size = 30)
    ax22.text(0.05, 0.90, r'$5.6 < Q^2 < 12.3 ~ {\rm GeV}^2; \theta_c = 39^{\circ}$', transform = ax22.transAxes, size = 30)


    handles, labels = [],[]
    handles.append(hand['JLab dp0'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{JLab Hall C}}')
    labels.append(r'\textbf{\textrm{JAM}}')
    ax21.legend(handles,labels,loc=(0.00,0.00), fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,left=0.06,right=0.99,top=0.99)
    py.savefig('%s/gallery/dis-JLab-Hall-C-E12-10-002.png' % (wdir))
    print('Saving figure to %s/gallery/dis-JLab-Hall-C-E12-10-002.png'%(wdir))
    py.close()



def plot_bonus_prelim(wdir, data, mode=1):

    if 1 not in data: return
    if 2 not in data: return
    if 3 not in data: return
    if 4 not in data: return
    if 5 not in data: return
    if 6 not in data: return

    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 7, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)

    dat_dp   = pd.read_excel('BONuS12dp.xlsx').to_dict(orient='list')
    dat_nd   = pd.read_excel('BONuS12nd.xlsx').to_dict(orient='list')
    dat_np   = pd.read_excel('BONuS12np.xlsx').to_dict(orient='list')

    dp_mean = np.array(dat_dp['JAM mean'])
    nd_mean = np.array(dat_nd['JAM mean'])
    np_mean = np.array(dat_np['JAM mean'])
    dp_std  = np.array(dat_dp['JAM std'])
    nd_std  = np.array(dat_nd['JAM std'])
    np_std  = np.array(dat_np['JAM std'])


    dp_replicas = np.array([dat_dp[_] for _ in list(dat_dp) if 'replica' in _])
    nd_replicas = np.array([dat_nd[_] for _ in list(dat_nd) if 'replica' in _])
    np_replicas = np.array([dat_np[_] for _ in list(dat_np) if 'replica' in _])

    nrep = len(dp_replicas)


    hand = {}
    X = dat_dp['X']
    Q2 = dat_dp['Q2']

    if mode == 0:
        ax = ax11
        for i in range(nrep):
            hand['dp'] ,= ax.plot(X,dp_replicas[i],color='firebrick',alpha=0.1,zorder=3)
            hand['nd'] ,= ax.plot(X,nd_replicas[i],color='darkgreen',alpha=0.1,zorder=2)
            hand['np'] ,= ax.plot(X,np_replicas[i],color='blue'     ,alpha=0.1,zorder=1)

    if mode == 1:
        ax = ax11
        #hand['dp'] = ax.fill_between(X,dp_mean-dp_std,dp_mean+dp_std,color='firebrick',alpha=0.8,zorder=3)
        #hand['nd'] = ax.fill_between(X,nd_mean-nd_std,nd_mean+nd_std,color='darkgreen',alpha=0.8,zorder=2)
        #hand['np'] = ax.fill_between(X,np_mean-np_std,np_mean+np_std,color='blue'     ,alpha=0.8,zorder=1)
        hand['dp'] = ax.errorbar(X,dp_mean,yerr=dp_std,color='firebrick',alpha=1.0,zorder=3,marker='o',ms=6,capsize=3.0,linestyle='none')
        hand['nd'] = ax.errorbar(X,nd_mean,yerr=nd_std,color='darkgreen',alpha=1.0,zorder=2,marker='^',ms=6,capsize=3.0,linestyle='none')
        hand['np'] = ax.errorbar(X,np_mean,yerr=np_std,color='blue'     ,alpha=1.0,zorder=1,marker='*',ms=6,capsize=3.0,linestyle='none')


    for ax in [ax11]:
        ax.set_xlim(0.10, 0.80)

        ax.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 30)

        ax.xaxis.set_tick_params(which = 'major', length = 8)
        ax.xaxis.set_tick_params(which = 'minor', length = 4)

        ax.set_xticks([0.2,0.4,0.6])
        ax.xaxis.set_minor_locator(MultipleLocator(0.10))

        ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=40)
        ax.xaxis.set_label_coords(0.95,0.00)

    ax11.set_ylim(0.30, 1.00)
    ax11.yaxis.set_tick_params(which = 'major', length = 8)
    ax11.yaxis.set_tick_params(which = 'minor', length = 4)
    ax11.set_yticks([0.4,0.6,0.8])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.10))


    ax11.text(0.28, 0.90, r'\textrm{\textbf{JAM Prediction}}', transform = ax11.transAxes, size = 40)
    ax11.text(0.45, 0.82, r'\textrm{BONuS Kinematics}', transform = ax11.transAxes, size = 30)

    handles, labels = [],[]
    handles.append(hand['dp'])
    handles.append(hand['nd'])
    handles.append(hand['np'])
    #handles.append((thy_band,thy_plot))
    labels.append(r'\boldmath$F_2^D/F_2^p$')
    labels.append(r'\boldmath$F_2^n/F_2^D$')
    labels.append(r'\boldmath$F_2^n/F_2^p$')
    ax11.legend(handles,labels,loc=(0.01,0.01), fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.10,left=None,right=0.99)
    filename = '%s/gallery/dis-bonus12-prelim'%wdir
    if mode==1: filename += '-bands'
    filename += '.png'
    py.savefig(filename)
    print('Saving figure to %s'%filename)
    py.close()


#--new plots
def plot_SLAC_E140x(wdir, data):

    if 10035 not in data: return
    if 10036 not in data: return

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #################################
    #--Plot F2p and F2D from SLAC E-140x
    #################################

    nbins = 4

    p   = data[10035]
    d   = data[10036]

    DATA = {}
    DATA['SLAC E140x p']  = pd.DataFrame(p)
    DATA['SLAC E140x d']  = pd.DataFrame(d)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'SLAC E140x')
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,functheta=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        if exp=='SLAC E140x p': ax = ax11
        if exp=='SLAC E140x d': ax = ax12
        for key in range(nbins):
            X     = PLOT[exp]['theta'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['theta'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color=color)
            thy_band  = ax.fill_between(X,down,up,color=color,alpha=0.5)

    for ax in [ax11,ax12]:
        ax.set_xlim(11,55)
        ax.set_ylim(0.82,1.48)
        ax.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 30)

        ax.set_xlabel(r'\boldmath$\theta$', size=35)
        ax.xaxis.set_label_coords(0.98,-0.02)

        ax.xaxis.set_tick_params(which = 'major', length = 8)
        ax.xaxis.set_tick_params(which = 'minor', length = 4)
        ax.yaxis.set_tick_params(which = 'major', length = 8)
        ax.yaxis.set_tick_params(which = 'minor', length = 4)

        ax.set_xticks([20,30,40,50])
        ax.set_yticks([0.9,1.0,1.1,1.2,1.3,1.4])

        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.xaxis.set_major_locator(MultipleLocator(10))

        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.10))

    ax12.tick_params(labelleft=False)

    ax12.text(0.05, 0.90, r'\textrm{\textbf{SLAC E140x}}', transform = ax12.transAxes, size = 40)
    ax11.text(0.85, 0.05, r'\boldmath$p$', transform = ax11.transAxes, size = 50)
    ax12.text(0.85, 0.05, r'\boldmath$D$', transform = ax12.transAxes, size = 50)

    #handles, labels = [],[]
    #handles.append((hand['SLAC E140x p'],thy_band,thy_plot))
    #labels.append(r'\textbf{\textrm{SLAC E140x}}')
    #ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.subplots_adjust(wspace=0.01,hspace=0.01)
    py.savefig('%s/gallery/dis-SLAC-E140x.png' % (wdir))
    print('Saving figure to %s/gallery/dis-SLAC-E140x.png'%(wdir))
    py.close()

def plot_JLab_E06_009(wdir, data):

    if 10077 not in data: return

    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)

    #################################
    #--Plot F2p and F2D from SLAC E-140x
    #################################

    nbins = 3

    d   = data[10077]

    DATA = {}
    DATA['JLab E06-009']  = pd.DataFrame(d)

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],'JLab E06-009')
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        ax = ax11
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*1.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*1.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*1.0**float(key)
            std  = theory[exp]['std'][key]  *1.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color=color)
            thy_band  = ax.fill_between(X,down,up,color=color,alpha=0.5)

    ax11.set_xlim(0.35, 0.68)
    ax11.set_xticks([0.40,0.50,0.60])
    ax11.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax11.xaxis.set_major_locator(MultipleLocator(0.10))

    ax11.set_ylim(0.05,0.23)
    ax11.set_yticks([0.10,0.15,0.20])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax11.yaxis.set_major_locator(MultipleLocator(0.05))

    ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    ax11.xaxis.set_tick_params(which = 'major', length = 8)
    ax11.xaxis.set_tick_params(which = 'minor', length = 4)
    ax11.yaxis.set_tick_params(which = 'major', length = 8)
    ax11.yaxis.set_tick_params(which = 'minor', length = 4)

    ax11.set_xlabel(r'\boldmath$x_{\rm bj}$', size=30)
    ax11.xaxis.set_label_coords(0.95,0.00)

    ax11.text(0.05, 0.05, r'\textrm{\textbf{JLab E06-009 }}'+r'\boldmath$D$', transform = ax11.transAxes, size = 25)
    ax11.text(0.42, 0.94, r'$Q^2 = 2~{\rm GeV}^2$', transform = ax11.transAxes, size = 20)
    ax11.text(0.25, 0.50, r'$Q^2 = 3~{\rm GeV}^2$', transform = ax11.transAxes, size = 20)
    ax11.text(0.55, 0.20, r'$Q^2 = 4~{\rm GeV}^2$', transform = ax11.transAxes, size = 20)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #minorLocator = MultipleLocator(0.02)
    #majorLocator = MultipleLocator(0.1)
    #ax11.xaxis.set_minor_locator(minorLocator)
    #ax11.xaxis.set_major_locator(majorLocator)
    #minorLocator = MultipleLocator(0.004)
    #majorLocator = MultipleLocator(0.02)
    #ax11.set_yticks([1.04, 1.06, 1.08])

    #handles, labels = [],[]
    #handles.append(hand['JLab E06-009'])
    #handles.append((thy_band,thy_plot))
    #labels.append(r'\textbf{\textrm{JLab E06-009}}')
    #ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-JLab_E06_009.png' % (wdir))
    print('Saving figure to %s/gallery/dis-JLab_E06_009.png'%(wdir))
    py.close()

def plot_JLab_E03_103(wdir, data):

    if 10045 not in data: return
    if 10046 not in data: return

    nrows, ncols = 1, 3
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)
    ax13 = py.subplot(nrows, ncols, 3)

    #################################
    #--Plot F2p and F2D from SLAC E-140x
    #################################


    p   = data[10045]
    d   = data[10046]

    DATA = {}
    DATA['JLab E03-103 p']  = pd.DataFrame(p)
    DATA['JLab E03-103 d E5.011']  = pd.DataFrame(d).query('Elab < 5.5')
    DATA['JLab E03-103 d E5.766']  = pd.DataFrame(d).query('Elab > 5.5')

    PLOT = {}
    theory = {}
    for exp in DATA:
        if exp=='JLab E03-103 p':        nbins = 2
        if exp=='JLab E03-103 d E5.011': nbins = 4
        if exp=='JLab E03-103 d E5.766': nbins = 6
        query = get_thetabins(DATA[exp],exp)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    k = 1.3
    #--plot data points
    for exp in PLOT:
        if exp=='JLab E03-103 p':        ax,nbins = ax11, 2
        if exp=='JLab E03-103 d E5.011': ax,nbins = ax12, 4
        if exp=='JLab E03-103 d E5.766': ax,nbins = ax13, 6
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*k**float(key)
            alpha = PLOT[exp]['alpha'][key]*k**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*k**float(key)
            std  = theory[exp]['std'][key]  *k**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color=color)
            thy_band  = ax.fill_between(X,down,up,color=color,alpha=0.5)


    for ax in [ax11,ax12,ax13]:
        ax.set_xlim(0.301,0.78)
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_major_locator(MultipleLocator(0.10))
        ax.set_xticks([0.40, 0.50, 0.60, 0.70])
        ax.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

        ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=30)
        ax.xaxis.set_label_coords(0.95,0.00)

        ax.xaxis.set_tick_params(which = 'major', length = 8)
        ax.xaxis.set_tick_params(which = 'minor', length = 4)
        ax.yaxis.set_tick_params(which = 'major', length = 8)
        ax.yaxis.set_tick_params(which = 'minor', length = 4)

    #ax11.set_ylim(1.02, 1.09)
    #ax11.set_yticks([1.04, 1.06, 1.08])



    ax11.text(0.05, 0.05, r'\textrm{\textbf{JLab E03-103 }}'+r'\boldmath$F_2$', transform = ax11.transAxes, size = 30)
    ax11.text(0.90, 0.90, r'\boldmath$p$', transform = ax11.transAxes, size = 30)
    ax12.text(0.90, 0.90, r'\boldmath$D$', transform = ax12.transAxes, size = 30)
    ax13.text(0.90, 0.90, r'\boldmath$D$', transform = ax13.transAxes, size = 30)

    ax12.text(0.60, 0.80, r'$E=5.011~{\rm GeV}$', transform = ax12.transAxes, size = 20)
    ax13.text(0.60, 0.80, r'$E=5.766~{\rm GeV}$', transform = ax13.transAxes, size = 20)

    ax11.text(0.05, 0.50, r'$\theta=40^\circ$', transform = ax11.transAxes, size = 20)
    ax11.text(0.25, 0.80, r'$\theta=50^\circ$', transform = ax11.transAxes, size = 20)

    ax12.text(0.20, 0.20, r'$\theta=24^\circ$', transform = ax12.transAxes, size = 20)
    ax12.text(0.30, 0.15, r'$29^\circ$',        transform = ax12.transAxes, size = 20)
    ax12.text(0.40, 0.10, r'$36^\circ$',        transform = ax12.transAxes, size = 20)
    ax12.text(0.50, 0.05, r'$46^\circ$',        transform = ax12.transAxes, size = 20)

    ax13.text(0.10, 0.05, r'$\theta=18^\circ$', transform = ax13.transAxes, size = 20)
    ax13.text(0.20, 0.05, r'$22^\circ$',        transform = ax13.transAxes, size = 20)
    ax13.text(0.30, 0.05, r'$26^\circ$',        transform = ax13.transAxes, size = 20)
    ax13.text(0.40, 0.05, r'$32^\circ$',        transform = ax13.transAxes, size = 20)
    ax13.text(0.50, 0.05, r'$40^\circ$',        transform = ax13.transAxes, size = 20)
    ax13.text(0.60, 0.05, r'$50^\circ$',        transform = ax13.transAxes, size = 20)


    #handles, labels = [],[]
    #handles.append(hand['JLab E03-103 p'])
    #handles.append((thy_band,thy_plot))
    #labels.append(r'\textbf{\textrm{JLab E03-103 p}}')
    #ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-JLab_E03_103.png' % (wdir))
    print('Saving figure to %s/gallery/dis-JLab_E03_103.png'%(wdir))
    py.close()


def plot_clas6_old(wdir, data):

    if 20005 not in data: return
    if 20006 not in data: return

    nrows, ncols = 2, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)
    ax21 = py.subplot(nrows, ncols, 3)
    ax22 = py.subplot(nrows, ncols, 4)

    #################################
    #--Plot F2p and F2D from SLAC E-140x
    #################################


    p   = data[20005]
    d   = data[20006]

    DATA = {}
    DATA['clas6 p Q2 < 3']  = pd.DataFrame(p).query("Q2 <= 3.0")
    DATA['clas6 d Q2 < 3']  = pd.DataFrame(d).query("Q2 <= 3.0")
    DATA['clas6 p Q2 > 3']  = pd.DataFrame(p).query("Q2 >  3.0")
    DATA['clas6 d Q2 > 3']  = pd.DataFrame(d).query("Q2 >  3.0")

    PLOT = {}
    theory = {}
    for exp in DATA:
        query = get_Q2bins(DATA[exp],exp)
        nbins = len(query)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        if exp=='clas6 p Q2 < 3': ax,nbins = ax11, 7
        if exp=='clas6 d Q2 < 3': ax,nbins = ax12, 7
        if exp=='clas6 p Q2 > 3': ax,nbins = ax21, 7
        if exp=='clas6 d Q2 > 3': ax,nbins = ax22, 7
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*2.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*2.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*2.0**float(key)
            std  = theory[exp]['std'][key]  *2.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax.fill_between(X,down,up,color='gold',alpha=1.0)


    #ax11.set_xlim(0.28, 0.78)
    #ax11.set_ylim(1.02, 1.09)

    #ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    #ax11.yaxis.set_tick_params(which = 'major', length = 5)
    #ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.xaxis.set_tick_params(which = 'major', length = 5)
    #ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.set_xlabel(r'\boldmath$x$', size=30)
    #ax11.xaxis.set_label_coords(0.95,0.00)

    ##ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #minorLocator = MultipleLocator(0.02)
    #majorLocator = MultipleLocator(0.1)
    #ax11.xaxis.set_minor_locator(minorLocator)
    #ax11.xaxis.set_major_locator(majorLocator)
    #minorLocator = MultipleLocator(0.004)
    #majorLocator = MultipleLocator(0.02)
    #ax11.yaxis.set_minor_locator(minorLocator)
    #ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])
    #ax11.set_yticks([1.04, 1.06, 1.08])

    handles, labels = [],[]
    handles.append(hand['clas6 p Q2 < 3'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{CLAS6}}')
    ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-clas6.png' % (wdir))
    print('Saving figure to %s/gallery/dis-clas6.png'%(wdir))
    py.close()

def plot_clas6_proton(wdir, data):

    if 10057 not in data: return

    data = data[10057]
    nrows, ncols = 7, 7
    py.figure(figsize = (ncols * 7, nrows * 7))

    _Q2 = np.unique(data['Q2'])
    nQ2 = len(_Q2)
    X     = data['X']
    val   = data['value']
    alpha = data['alpha']
    Q2    = data['Q2']

    mean  = data['thy-0']
    std   = data['dthy-0']
    up    = mean + std
    down  = mean - std

    hand = {}
    ax = {}
    for i in range(nQ2):
        ax[i] = py.subplot(nrows,ncols,i+1)
        color, marker, ms = 'firebrick', 'o', 6.0
        q2 = _Q2[i]
        hand['clas6'] = ax[i].errorbar(X[Q2==q2],val[Q2==q2],alpha[Q2==q2],marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')
        thy_plot ,= ax[i].plot(X[Q2==q2],mean[Q2==q2],linestyle='solid',color='black')
        thy_band  = ax[i].fill_between(X[Q2==q2],down[Q2==q2],up[Q2==q2],color='gold',alpha=1.0)
        



    #ax11.set_xlim(0.28, 0.78)
    #ax11.set_ylim(1.02, 1.09)

    #ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    #ax11.yaxis.set_tick_params(which = 'major', length = 5)
    #ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.xaxis.set_tick_params(which = 'major', length = 5)
    #ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.set_xlabel(r'\boldmath$x$', size=30)
    #ax11.xaxis.set_label_coords(0.95,0.00)

    ##ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #minorLocator = MultipleLocator(0.02)
    #majorLocator = MultipleLocator(0.1)
    #ax11.xaxis.set_minor_locator(minorLocator)
    #ax11.xaxis.set_major_locator(majorLocator)
    #minorLocator = MultipleLocator(0.004)
    #majorLocator = MultipleLocator(0.02)
    #ax11.yaxis.set_minor_locator(minorLocator)
    #ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])
    #ax11.set_yticks([1.04, 1.06, 1.08])

    handles, labels = [],[]
    handles.append(hand['clas6'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{CLAS6}}')
    ax[1].legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-clas6-proton.png' % (wdir))
    print('Saving figure to %s/gallery/dis-clas6-proton.png'%(wdir))
    py.close()

def plot_clas6_deuteron(wdir, data):

    if 10058 not in data: return

    data = data[10058]
    nrows, ncols = 10, 9
    py.figure(figsize = (ncols * 7, nrows * 7))

    _Q2 = np.unique(data['Q2'])
    nQ2 = len(_Q2)
    X     = data['X']
    val   = data['value']
    alpha = data['alpha']
    Q2    = data['Q2']

    mean  = data['thy-0']
    std   = data['dthy-0']
    up    = mean + std
    down  = mean - std

    hand = {}
    ax = {}
    for i in range(nQ2):
        ax[i] = py.subplot(nrows,ncols,i+1)
        color, marker, ms = 'firebrick', 'o', 6.0
        q2 = _Q2[i]
        hand['clas6'] = ax[i].errorbar(X[Q2==q2],val[Q2==q2],alpha[Q2==q2],marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')
        thy_plot ,= ax[i].plot(X[Q2==q2],mean[Q2==q2],linestyle='solid',color='black')
        thy_band  = ax[i].fill_between(X[Q2==q2],down[Q2==q2],up[Q2==q2],color='gold',alpha=1.0)
        



    #ax11.set_xlim(0.28, 0.78)
    #ax11.set_ylim(1.02, 1.09)

    #ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    #ax11.yaxis.set_tick_params(which = 'major', length = 5)
    #ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.xaxis.set_tick_params(which = 'major', length = 5)
    #ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.set_xlabel(r'\boldmath$x$', size=30)
    #ax11.xaxis.set_label_coords(0.95,0.00)

    ##ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #minorLocator = MultipleLocator(0.02)
    #majorLocator = MultipleLocator(0.1)
    #ax11.xaxis.set_minor_locator(minorLocator)
    #ax11.xaxis.set_major_locator(majorLocator)
    #minorLocator = MultipleLocator(0.004)
    #majorLocator = MultipleLocator(0.02)
    #ax11.yaxis.set_minor_locator(minorLocator)
    #ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])
    #ax11.set_yticks([1.04, 1.06, 1.08])

    handles, labels = [],[]
    handles.append(hand['clas6'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{CLAS6}}')
    ax[1].legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-clas6-deuteron.png' % (wdir))
    print('Saving figure to %s/gallery/dis-clas6-deuteron.png'%(wdir))
    py.close()


def plot_bonus_neutron(wdir, data):

    if 10061 not in data: return

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #################################
    #--Plot F2p and F2D from SLAC E-140x
    #################################


    n   = data[10061]

    DATA = {}
    DATA['BONuS n E4.223']  = pd.DataFrame(n).query('Elab < 5.0')
    DATA['BONuS n E5.262']  = pd.DataFrame(n).query('Elab > 5.0')

    PLOT = {}
    theory = {}
    for exp in DATA:
        if exp=='BONuS n E4.223': nbins = 5
        if exp=='BONuS n E5.262': nbins = 7
        query = get_Q2bins(DATA[exp],exp)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    #--plot data points
    for exp in PLOT:
        if exp=='BONuS n E4.223': ax,nbins = ax11, 5
        if exp=='BONuS n E5.262': ax,nbins = ax12, 7
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key]*3.0**float(key)
            alpha = PLOT[exp]['alpha'][key]*3.0**float(key)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')

        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key]*3.0**float(key)
            std  = theory[exp]['std'][key]  *3.0**float(key)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color='black')
            thy_band  = ax.fill_between(X,down,up,color='gold',alpha=1.0)


    #ax11.set_xlim(0.28, 0.78)
    #ax11.set_ylim(1.02, 1.09)

    #ax11.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

    #ax11.yaxis.set_tick_params(which = 'major', length = 5)
    #ax11.yaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.xaxis.set_tick_params(which = 'major', length = 5)
    #ax11.xaxis.set_tick_params(which = 'minor', length = 2.5)

    #ax11.set_xlabel(r'\boldmath$x$', size=30)
    #ax11.xaxis.set_label_coords(0.95,0.00)

    ##ax11.text(0.05, 0.90, r'\textrm{\textbf{MARATHON}}', transform = ax11.transAxes, size = 30)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #minorLocator = MultipleLocator(0.02)
    #majorLocator = MultipleLocator(0.1)
    #ax11.xaxis.set_minor_locator(minorLocator)
    #ax11.xaxis.set_major_locator(majorLocator)
    #minorLocator = MultipleLocator(0.004)
    #majorLocator = MultipleLocator(0.02)
    #ax11.yaxis.set_minor_locator(minorLocator)
    #ax11.yaxis.set_major_locator(majorLocator)
    #ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])
    #ax11.set_yticks([1.04, 1.06, 1.08])

    handles, labels = [],[]
    handles.append(hand['BONuS n E4.223'])
    handles.append((thy_band,thy_plot))
    labels.append(r'\textbf{\textrm{BONuS}}')
    ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-BONuS-neutron.png' % (wdir))
    print('Saving figure to %s/gallery/dis-BONuS-neutron.png'%(wdir))
    py.close()



def plot_E665(wdir, data):

    if 10062 not in data: return
    if 10063 not in data: return

    _data = data
    data = {}
    data['p'] = _data[10062]
    data['d'] = _data[10063]
    nrows, ncols = 3, 4
    py.figure(figsize = (ncols * 7, nrows * 7))

    thy_plot, thy_band = {}, {}
    hand = {}
    for tar in ['p','d']:
        _Q2 = np.unique(data[tar]['Q2'])
        nQ2 = len(_Q2)
        X     = data[tar]['X']
        val   = data[tar]['value']
        alpha = data[tar]['alpha']
        Q2    = data[tar]['Q2']

        mean  = data[tar]['thy-0']
        std   = data[tar]['dthy-0']
        up    = mean + std
        down  = mean - std

        if tar=='p': color, marker, ms = 'firebrick', 'o', 10.0
        if tar=='d': color, marker, ms = 'darkgreen', '^', 10.0

        if tar=='p': shift = 0.00
        if tar=='d': shift = 0.001
        ax = {}
        for i in range(nQ2):
            ax[i] = py.subplot(nrows,ncols,i+1)
            q2 = _Q2[i]
            hand['E665 %s'%tar] = ax[i].errorbar(X[Q2==q2]+shift,val[Q2==q2],alpha[Q2==q2],marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')
            thy_plot[tar] ,= ax[i].plot(X[Q2==q2]+shift,mean[Q2==q2],linestyle='solid',color=color)
            thy_band[tar]  = ax[i].fill_between(X[Q2==q2]+shift,down[Q2==q2],up[Q2==q2],color=color,alpha=0.5)

            if len(X[Q2==q2])==1:
                x = X[Q2==q2][0]
                x = np.array([x-0.1,x+0.1])
                ax[i].plot(x+shift,np.array([mean[Q2==q2],mean[Q2==q2]]),linestyle='solid',color=color)
                ax[i].fill_between(x+shift,down[Q2==q2],up[Q2==q2],color=color,alpha=0.5)

            if tar=='p':
                ax[i].text(0.04, 0.04, r'$Q^2 = %s~{\rm GeV}^2$'%np.round(q2,2), transform = ax[i].transAxes, size = 30)
             

    for i in range(12):
        ax[i].semilogx()
        ax[i].set_xlim(2e-3, 0.50)
        ax[i].tick_params(axis = 'x', which = 'both', direction='in', top = True, right = True, left=True, labelsize = 40)
        ax[i].tick_params(axis = 'y', which = 'both', direction='in', top = True, right = True, left=True, labelsize = 35)
        ax[i].xaxis.set_tick_params(which = 'major', length = 8)
        ax[i].xaxis.set_tick_params(which = 'minor', length = 4)
        ax[i].yaxis.set_tick_params(which = 'major', length = 8)
        ax[i].yaxis.set_tick_params(which = 'minor', length = 4)



    for i in range(8):
        ax[i].tick_params(labelbottom=False)

    for i in [1,2,3,5,6,7,9,10,11]:
        ax[i].tick_params(labelleft=False)
    

    for i in [8,9,10,11]:
        ax[i].set_xlabel(r'\boldmath$x$', size=60)
        ax[i].xaxis.set_label_coords(0.95,0.00)
        ax[i].tick_params(axis='x',pad=10)

    for i in [0,1,2,3]:
        ax[i].set_ylim(0.30,0.54)
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[i].set_yticks([0.40,0.50])

    for i in [4,5,6,7]:
        ax[i].set_ylim(0.25,0.56)
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[i].set_yticks([0.30,0.40,0.50])
    
    for i in [8,9,10,11]:
        ax[i].set_ylim(0.08,0.60)
        ax[i].yaxis.set_minor_locator(MultipleLocator(0.05))
        ax[i].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[i].set_yticks([0.10,0.20,0.30,0.40,0.50])



    ax[0].text(0.05, 0.83, r'\textrm{\textbf{E665 }}' + r'\boldmath$F_2$', transform = ax[0].transAxes, size = 70)
    #ax11.text(0.05, 0.90, r'\boldmath$F_2^{^3\rm{He}}/F_2^{d}$', transform = ax11.transAxes, size = 30)

    #ax11.set_xticks([0.30, 0.40, 0.50, 0.60, 0.70])

    handles, labels = [],[]
    handles.append((hand['E665 p'],thy_band['p'],thy_plot['p']))
    handles.append((hand['E665 d'],thy_band['d'],thy_plot['d']))
    labels.append(r'\textbf{\textrm{p}}')
    labels.append(r'\textbf{\textrm{D}}')
    ax[8].legend(handles,labels,loc=(0.01,0.60), fontsize = 50, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.subplots_adjust(wspace=0.01,hspace=0.01)
    py.savefig('%s/gallery/dis-E665.png' % (wdir))
    print('Saving figure to %s/gallery/dis-E665.png'%(wdir))
    py.close()

def plot_JLab_JLCEE96(wdir, data):

    if 10072 not in data: return
    if 10073 not in data: return

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 6, nrows * 7))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #################################
    #--Plot F2p and F2D
    #################################


    p   = data[10072]
    d   = data[10073]

    DATA = {}
    DATA['JLab JLCEE96 p']  = pd.DataFrame(p)
    DATA['JLab JLCEE96 d']  = pd.DataFrame(d)

    PLOT = {}
    theory = {}
    for exp in DATA:
        if exp=='JLab JLCEE96 p': nbins = 9
        if exp=='JLab JLCEE96 d': nbins = 8
        query = get_Q2bins(DATA[exp],exp)
        PLOT[exp] = get_plot(query)
        theory[exp] = get_theory(PLOT[exp],nbins,loop=False,funcX=True)


    hand = {}
    k = 2.0
    #--plot data points
    for exp in PLOT:
        i = 0
        if exp=='JLab JLCEE96 p': ax,nbins = ax11, 9
        if exp=='JLab JLCEE96 d': ax,nbins = ax12, 8
        for key in range(nbins):
            X     = PLOT[exp]['X'][key]
            if len(X)==0: continue
            val   = PLOT[exp]['value'][key]*k**float(i)
            alpha = PLOT[exp]['alpha'][key]*k**float(i)
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=2.5,linestyle='none')
            i+=1

        i = 0
        #--plot theory interpolated between data points
        for key in theory[exp]['value']:
            X    = theory[exp]['X'][key]
            if len(X)==0: continue
            mean = theory[exp]['value'][key]*k**float(i)
            std  = theory[exp]['std'][key]  *k**float(i)
            down = mean - std
            up   = mean + std
            thy_plot ,= ax.plot(X,mean,linestyle='solid',color=color)
            thy_band  = ax.fill_between(X,down,up,color=color,alpha=0.5)
            i+=1


    for ax in [ax11,ax12]:
        ax.set_xlim(0.301,0.78)
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_major_locator(MultipleLocator(0.10))
        ax.set_xticks([0.40, 0.50, 0.60, 0.70])
        ax.tick_params(axis = 'both', which = 'both', direction='in', top = True, right = True, labelsize = 20)

        ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=30)
        ax.xaxis.set_label_coords(0.95,0.00)

        ax.xaxis.set_tick_params(which = 'major', length = 8)
        ax.xaxis.set_tick_params(which = 'minor', length = 4)
        ax.yaxis.set_tick_params(which = 'major', length = 8)
        ax.yaxis.set_tick_params(which = 'minor', length = 4)

    #ax11.set_ylim(1.02, 1.09)
    #ax11.set_yticks([1.04, 1.06, 1.08])



    ax11.text(0.05, 0.05, r'\textrm{\textbf{JLab JLCEE96 }}'+r'\boldmath$F_2$', transform = ax11.transAxes, size = 30)
    #ax11.text(0.90, 0.90, r'\boldmath$p$', transform = ax11.transAxes, size = 30)
    #ax12.text(0.90, 0.90, r'\boldmath$D$', transform = ax12.transAxes, size = 30)

    #ax12.text(0.60, 0.80, r'$E=5.011~{\rm GeV}$', transform = ax12.transAxes, size = 20)

    #ax11.text(0.05, 0.50, r'$\theta=40^\circ$', transform = ax11.transAxes, size = 20)
    #ax11.text(0.25, 0.80, r'$\theta=50^\circ$', transform = ax11.transAxes, size = 20)

    #ax12.text(0.20, 0.20, r'$\theta=24^\circ$', transform = ax12.transAxes, size = 20)
    #ax12.text(0.30, 0.15, r'$29^\circ$',        transform = ax12.transAxes, size = 20)
    #ax12.text(0.40, 0.10, r'$36^\circ$',        transform = ax12.transAxes, size = 20)
    #ax12.text(0.50, 0.05, r'$46^\circ$',        transform = ax12.transAxes, size = 20)



    #handles, labels = [],[]
    #handles.append(hand['JLab E03-103 p'])
    #handles.append((thy_band,thy_plot))
    #labels.append(r'\textbf{\textrm{JLab E03-103 p}}')
    #ax11.legend(handles,labels,loc='lower right', fontsize = 25, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.savefig('%s/gallery/dis-JLab_JLCEE96.png' % (wdir))
    print('Saving figure to %s/gallery/dis-JLab_JLCEE96.png'%(wdir))
    py.close()



def plot_hera_hq(wdir,data):

    if 10037 not in data: return
    if 10038 not in data: return

    print('\ngenerating HERA heavy quark plot from %s'%(wdir))

    nrows,ncols=3,4
    fig = py.figure(figsize=(ncols*5,nrows*3))
    ax = {}
    for i in range(12):
        ax[i+1] = py.subplot(nrows,ncols,i+1)

    #######################
    #--plot absolute values
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in data:
        if   idx==10037: color,fmt = 'firebrick','o'
        elif idx==10038: color,fmt = 'darkgreen','^'
        else: continue        
   
        X = data[idx]['X']
        Q2 = data[idx]['Q2']
        nbins = len(np.unique(Q2))+1
        values = data[idx]['value']
        alpha = data[idx]['alpha']
        BIN = data[idx]['bin']
        _thy = data[idx]['thy-0']
        _std = data[idx]['dthy-0']
        for i in range(nbins):
            #lprint('generating bin: [%s/%s]'%(i+1,nbins))
            if (i+1) not in BIN: continue
            x,val,alp,thy,std = [],[],[],[],[]
            for j in range(len(values)):
                if BIN[j] != (i+1): continue
                val.append(values[j])
                alp.append(alpha[j])
                x.append(X[j])
                thy.append(_thy[j])
                std.append(_std[j])
                hand[idx] = ax[i+1].errorbar(x,val,yerr=alp,color=color,fmt='o',ms=2.0,capsize=3.0)
            ax[i+1].text(0.05,0.85,r'$Q^2=%d~{\rm GeV}^2$'%(np.unique(Q2)[i]),transform=ax[i+1].transAxes,size=22)
            #--plot mean and std of all replicas
            thy = np.array(thy)
            std = np.array(std)
            down = thy - std
            up   = thy + std
            thy_plot[idx] ,= ax[i+1].plot(x,thy,color=color)
            thy_band[idx]  = ax[i+1].fill_between(x,down,up,color=color,alpha=0.4)

            if len(x)==1:
                x = [x[0]*0.92,x[0]*1.08]
                down, up, thy = [down[0], down[0]], [up[0],up[0]], [thy[0],thy[0]]
                thy_plot[idx] ,= ax[i+1].plot(x,thy,color=color)
                thy_band[idx]  = ax[i+1].fill_between(x,down,up,color=color,alpha=0.4)


    for i in range(12):
        ax[i+1].semilogx()
        ax[i+1].semilogy()
        ax[i+1].tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        #minorLocator = MultipleLocator(0.1)
        #majorLocator = MultipleLocator(0.5)
        #ax[i+1].xaxis.set_minor_locator(minorLocator)
        #ax[i+1].xaxis.set_major_locator(majorLocator)
        ax[i+1].xaxis.set_tick_params(which='major',length=6)
        ax[i+1].xaxis.set_tick_params(which='minor',length=3)
        ax[i+1].yaxis.set_tick_params(which='major',length=6)
        ax[i+1].yaxis.set_tick_params(which='minor',length=3)
        #ax[i+1].set_xticks([0.5,1.0,1.5,2.0])

    for i in [1,5,9]:
        ax[i].set_xlim(2e-5,2e-2)
        ax[i].set_xticks([1e-4,1e-3,1e-2])

    for i in [2,6,10]:
        ax[i].set_xlim(6e-5,3e-2)
        ax[i].set_xticks([1e-4,1e-3,1e-2])

    for i in [3,7,11]:
        ax[i].set_xlim(1e-4,4e-2)
        ax[i].set_xticks([1e-3,1e-2])

    for i in [4,8,12]:
        ax[i].set_xlim(1e-4,6e-2)
        ax[i].set_xticks([1e-3,1e-2])

    for i in range(8):
        ax[i+1].tick_params(labelbottom=False)

    for i in [2,3,4,6,7,8,10,11,12]:
        ax[i].tick_params(labelleft=False)

    for i in [9,10,11,12]:
        ax[i].set_xlabel(r'\boldmath$x_{\rm bj}$',size=30)
        #ax[i].xaxis.set_label_coords(0.95,0.00)


    for i in [1,2,3,4]:
        ax[i].set_ylim(5e-4,0.9)

    for i in [5,6,7,8]:
        ax[i].set_ylim(5e-4,0.9)

    for i in [9,10,11,12]:
        ax[i].set_ylim(5e-4,0.9)


    #for i in [5,6,7,8]:
    #    ax[i].set_ylim(0.0,13)
    #    minorLocator = MultipleLocator(1)
    #    majorLocator = MultipleLocator(5)
    #    ax[i].yaxis.set_minor_locator(minorLocator)
    #    ax[i].yaxis.set_major_locator(majorLocator)
    #    ax[i].yaxis.set_tick_params(which='major',length=6)
    #    ax[i].yaxis.set_tick_params(which='minor',length=3)
    #    ax[i].set_yticks([5,10])

    #for i in [9,10,11,12]:
    #    ax[i].set_ylim(0.0,3.5)
    #    minorLocator = MultipleLocator(0.5)
    #    majorLocator = MultipleLocator(1)
    #    ax[i].yaxis.set_minor_locator(minorLocator)
    #    ax[i].yaxis.set_major_locator(majorLocator)
    #    ax[i].yaxis.set_tick_params(which='major',length=6)
    #    ax[i].yaxis.set_tick_params(which='minor',length=3)
    #    ax[i].set_yticks([1,2,3])

    #for i in [13,14,15,16]:
    #    ax[i].set_ylim(0.0,0.7)
    #    minorLocator = MultipleLocator(0.1)
    #    majorLocator = MultipleLocator(0.2)
    #    ax[i].yaxis.set_minor_locator(minorLocator)
    #    ax[i].yaxis.set_major_locator(majorLocator)
    #    ax[i].yaxis.set_tick_params(which='major',length=6)
    #    ax[i].yaxis.set_tick_params(which='minor',length=3)
    #    ax[i].set_yticks([0.0,0.2,0.4,0.6])


    ##ax[1] .text(0.05, 0.65, r'\boldmath${\rm d}^2 \sigma/{\rm d}z {\rm d}M_{h} ~ [{\rm nb}/{\rm GeV}]$',transform=ax[1].transAxes,size=25)
    #ax[1] .text(0.05, 0.75, r'\boldmath$\frac{{\rm d}^2 \sigma}{{\rm d}z {\rm d}M_{h}} ~ [{\rm nb}/{\rm GeV}]$',transform=ax[1].transAxes,size=25)
    #ax[1].text(0.05, 0.55, r'$\sqrt{s} = 10.58$'+' '+r'\textrm{GeV}', transform=ax[1].transAxes,size=20)

    #thy ,= ax[2].plot([],[],color='black')
    #thy_band = ax[2].fill_between([],[],[],color='gold',alpha=1.0)


    #ax[1].axis('off')

    handles,labels = [], []
    handles.append((thy_plot[10037],thy_band[10037],hand[10037]))
    handles.append((thy_plot[10038],thy_band[10038],hand[10038]))
    labels.append(r'\textbf{\textrm{HERA \boldmath$c$}}') 
    labels.append(r'\textbf{\textrm{HERA \boldmath$b$}}') 
    ax[9].legend(handles,labels,frameon=False,fontsize=22,loc='lower left',handletextpad = 0.2, handlelength = 1.0, labelspacing=0.7)

    py.tight_layout()
    py.subplots_adjust(hspace=0.02,wspace=0.02)


    checkdir('%s/gallery'%wdir)
    filename='%s/gallery/idis-hera-hq'%wdir
    filename+='.png'

    py.savefig(filename)
    print('Saving Belle plot to %s'%filename)


def plot_obs(wdir, kc):

    print('\nplotting dis data from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'idis' not in predictions['reactions']:
        print('inclusive DIS is not in data file')
        return
    #labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    #cluster = labels['cluster']

    data = predictions['reactions']['idis']

    for idx in data:
        predictions = copy.copy(data[idx]['prediction-rep'])
        shifts      = copy.copy(data[idx]['shift-rep'])
        #predictions = copy.copy(np.array(data[idx]['prediction-rep'])-np.array(data[idx]['shift-rep']))
        del data[idx]['prediction-rep']
        del data[idx]['residuals-rep']
        del data[idx]['shift-rep']
        del data[idx]['r-residuals']
        del data[idx]['n-residuals']
        if 'rres-rep' in data[idx]: del data[idx]['rres-rep']
        for ic in range(kc.nc[istep]):
            predictions_ic = [predictions[i] for i in range(len(predictions))]
            #predictions_ic = [predictions[i]-shifts[i] for i in range(len(predictions))]
            data[idx]['thy-%d' % ic] = np.mean(predictions_ic, axis = 0)
            data[idx]['dthy-%d' % ic] = np.std(predictions_ic, axis = 0)
            if 'X' in data[idx]: data[idx]['x'] = data[idx]['X']
            data[idx]['rQ2'] = np.around(data[idx]['Q2'], decimals = 0)
            data[idx]['rx'] = np.around(data[idx]['x'], decimals = 2)

    #plot_bonus_prelim(wdir,data,kc,istep)

    plot_hera_hq(wdir,data)

    plot_proton  (wdir, data)
    plot_deuteron(wdir, data)
    plot_CC      (wdir, data)
    plot_marathon(wdir, data)
    plot_marathon2(wdir, data)
    plot_JLab_Hall_C_E03_103(wdir, data)
    plot_JLab_Hall_C_E12_10_002(wdir, data)

    plot_SLAC_E140x(wdir, data)
    plot_JLab_E06_009(wdir, data)
    plot_JLab_E03_103(wdir, data)
    plot_clas6_proton(wdir, data)
    plot_clas6_deuteron(wdir, data)
    plot_bonus_neutron(wdir, data)
    plot_E665(wdir, data)
    plot_JLab_JLCEE96(wdir, data)
    
    return






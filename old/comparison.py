#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as  np

from core.coarse_decor import coarse_decor
from core.BLS_ms import BLS
from scipy.signal import lombscargle

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

P = 2.21857312
Tp = 2454037.612

with h5py.File('/data2/talens/2015Q2/LPE/tmp0_201504ALPE.hdf5', 'r') as f:
    lc1 = f['807144'].value
    
with h5py.File('/data2/talens/2015Q2/LPE/sigmas/tmp0_201504ALPE.hdf5', 'r') as f:
    lc2 = f['807144'].value

plt.subplot(211)
plt.plot(lc1['mag0'], '.')
plt.plot(lc2['mag0'], '.')
plt.subplot(212)
plt.plot(lc1['mag0'] - lc2['mag0'], '.')


plt.show()
    

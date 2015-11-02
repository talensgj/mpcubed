#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
        
with h5py.File('/data2/talens/3mEast/LBtests/huh.hdf5', 'r') as f:

    s = f['data/s'].value
    sigma2 = f['data/sigma2'].value

ind1, = np.where(np.all(np.isnan(sigma2), axis=0))
ind2, = np.where(np.all(np.isnan(sigma2), axis=1))

s = np.delete(s, ind1, axis=1)
s = np.delete(s, ind2, axis=0)
sigma2 = np.delete(sigma2, ind1, axis=1)
sigma2 = np.delete(sigma2, ind2, axis=0)

plt.subplot(211)
plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.subplot(212)
plt.imshow(sigma2.T, aspect='auto', cmap=viridis, vmin=0, vmax=.5)
plt.colorbar()

plt.show()

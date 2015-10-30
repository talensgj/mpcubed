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

with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter2.hdf5', 'r') as f:
    m1 = f['data/m'].value
    s1 = f['data/s'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter3.hdf5', 'r') as f:
    m2 = f['data/m'].value
    s2 = f['data/s'].value
    
plt.plot(m1, m1 - m2, '.')
plt.show()

delta = s1 - s2

ind1, = np.where(np.all(np.isnan(delta), axis=0))
ind2, = np.where(np.all(np.isnan(delta), axis=1))
s1 = np.delete(s1, ind1, axis=1)
s2 = np.delete(s2, ind1, axis=1)
delta = np.delete(delta, ind1, axis=1)
s1 = np.delete(s1, ind2, axis=0)
s2 = np.delete(s2, ind2, axis=0)
delta = np.delete(delta, ind2, axis=0)

plt.figure(figsize=(16,8))

ax = plt.subplot(221)
plt.imshow(s1.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()

plt.subplot(222, sharex=ax, sharey=ax)
plt.imshow(s2.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()

plt.subplot(223, sharex=ax, sharey=ax)
plt.imshow(delta.T, aspect='auto', cmap=viridis, vmin=-.05, vmax=.05)
plt.colorbar()

plt.show()

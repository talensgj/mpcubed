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

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter1.hdf5', 'r') as f:
    m1 = f['data/m'].value
    z1 = f['data/z'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter2.hdf5', 'r') as f:
    m2 = f['data/m'].value
    z2 = f['data/z'].value
    
plt.plot(m1, m1 - m2, '.')
plt.show()
   
z1 = z1.reshape((13502, 722))
z2 = z2.reshape((13502, 722)) 
delta = z1 - z2

plt.figure(figsize=(16,8))

ax = plt.subplot(221)
plt.imshow(z1.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()

plt.subplot(222, sharex=ax, sharey=ax)
plt.imshow(z2.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()

plt.subplot(223, sharex=ax, sharey=ax)
plt.imshow(delta.T, aspect='auto', cmap=viridis, vmin=-.05, vmax=.05)
plt.colorbar()

plt.show()

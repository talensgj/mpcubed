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

plt.figure(figsize=(16,8))

with h5py.File('/data2/talens/3mEast/fLC_20150605LPE.hdf5', 'r') as f:
    lc = f['data/807144'].value
    
ax = plt.subplot(311)
plt.title('05-06-2015')
plt.plot(lc['lst'], lc['flux0'], '.')
plt.ylabel('Flux')

with h5py.File('/data2/talens/3mEast/fLC_20150611LPE.hdf5', 'r') as f:
    lc = f['data/807144'].value
    
plt.subplot(312, sharex=ax, sharey=ax)
plt.title('11-06-2015')
plt.plot(lc['lst'], lc['flux0'], '.')
plt.ylabel('Flux')

with h5py.File('/data2/talens/3mEast/fLC_20150629LPE.hdf5', 'r') as f:
    lc = f['data/807144'].value
    
plt.subplot(313, sharex=ax, sharey=ax)
plt.title('29-06-2015')
plt.plot(lc['lst'], lc['flux0'], '.')
plt.ylabel('Flux')

plt.xlabel('LST')
    
plt.tight_layout()
plt.show()

exit()
with h5py.File('/data2/talens/3mEast/Stability/camip_20150629LPE_iter1.hdf5', 'r') as f:
    m1 = f['data/m'].value
    z1 = f['data/z'].value
    
with h5py.File('/data2/talens/3mEast/Stability/camip_20150611LPE_iter1.hdf5', 'r') as f:
    m2 = f['data/m'].value
    z2 = f['data/z'].value
   
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
plt.imshow(delta.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()

plt.show()

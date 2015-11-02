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

with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5_weights.hdf5', 'r') as f:
    mz = f['data/m'].value
    z = f['data/z'].value
    a = f['data/a'].value
    b = f['data/b'].value
    c = f['data/c'].value
    d = f['data/d'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter4_weights.hdf5', 'r') as f:
    ms = f['data/m'].value
    s = f['data/s'].value
    
fig = plt.figure()
ax = plt.subplot(221)
plt.plot(dec, mz-vmag, '.')
plt.subplot(222, sharex=ax, sharey=ax)
plt.plot(dec, ms-vmag, '.')
plt.subplot(223, sharex=ax)
plt.plot(dec, mz-ms, '.')
plt.ylim(-.1,.1)
plt.show()

Ax = np.sqrt(a**2 + b**2)
Ay = np.sqrt(c**2 + d**2)
phix = np.arctan2(b, a)
phiy = np.arctan2(d, c)
    
z = z.reshape((13502, 722))
Ax = Ax.reshape((272, 722))
Ay = Ay.reshape((272, 722))
phix = phix.reshape((272, 722))
phiy = phiy.reshape((272, 722))

plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.show()

plt.figure(figsize=(16,8))
ax = plt.subplot(121)
plt.imshow(Ax.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)
plt.colorbar()
plt.subplot(122, sharex=ax, sharey=ax)
plt.imshow(Ay.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)
plt.colorbar()
plt.tight_layout()
plt.show()

plt.figure(figsize=(16,8))
ax = plt.subplot(121)
plt.imshow(phix.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.colorbar()
plt.subplot(122, sharex=ax, sharey=ax)
plt.imshow(phiy.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.colorbar()
plt.tight_layout()
plt.show()

ind, = np.where(np.all(np.isnan(s), axis=0))
s = np.delete(s, ind, axis=1)
ind, = np.where(np.all(np.isnan(s), axis=1))
s = np.delete(s, ind, axis=0)

plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.show()

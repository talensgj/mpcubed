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

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day.hdf5', 'r') as f:
    m = f['data/m'].value
    z = f['data/z'].value
    
plt.plot(dec, m - vmag, '.')
plt.show()
    
z = z.reshape((13502, 722))

plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.show()

with h5py.File('/data2/talens/3mEast/LBtests/sky_15day.hdf5', 'r') as f:
    m = f['data/m'].value
    s = f['data/s'].value

plt.plot(dec, m - vmag, '.')
plt.show()

plt.scatter(ra, dec, c = m - vmag)
plt.show()

ind, = np.where(np.all(np.isnan(s), axis=0))
s = np.delete(s, ind, axis=1)

plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.show()

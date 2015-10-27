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

with h5py.File('/data2/talens/3mEast/LBtests/camip_test.hdf5', 'r') as f:
    z = f['data/z'].value
    
z = z.reshape((13502, 722))

plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar()
plt.show()

#with h5py.File('/data2/talens/3mEast/LBtests/sky_test.hdf5', 'r') as f:
    #s = f['data/s'].value

#plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
#plt.colorbar()
#plt.show()

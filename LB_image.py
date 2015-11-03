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
        
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter4_weights.hdf5', 'r') as f:
    
    sigma1 = f['data/sigma1'].value
    sigma2 = f['data/sigma2'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter5_weights.hdf5', 'r') as f:
    
    sigma1b = f['data/sigma1'].value
    sigma2b = f['data/sigma2'].value

plt.plot(vmag, sigma1, '.')
plt.plot(vmag, sigma1b, '.')
plt.show()

plt.plot(sigma1, sigma1 - sigma1b, '.')
plt.show()
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

#with h5py.File('/data2/talens/3mEast/LBtests/cam_15day.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/LBtests/reformat.hdf5') as g:
    
    #f.copy('header', g)
    
    #ascc = f['data/ascc'].value
    #m = f['data/m'].value
    #z = f['data/z'].value
    #idx, = np.where(~np.isnan(z))
    #values = z[idx]
    
    #grp = g.create_group('data')
    #grp.create_dataset('magnitudes/ascc', data = ascc)
    #grp.create_dataset('magnitudes/magnitude', data = m)
    
    #grp.create_dataset('camtrans/idx', data=idx)
    #grp.create_dataset('camtrans/value', data=values)
    #grp['camtrans'].attrs['grid'] = 'polar'
    #grp['camtrans'].attrs['nx'] = 13500
    #grp['camtrans'].attrs['ny'] = 720
    
with h5py.File('/data2/talens/3mEast/LBtests/sky_15day.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/LBtests/reformat.hdf5') as g:
    
    f.copy('header', g)
    
    ascc = f['data/ascc'].value
    m = f['data/m'].value
    s = f['data/s'].value
    idx, lstseq, = np.where(~np.isnan(s))
    values = s[idx, lstseq]
    
    tmax = s.shape[1]
    
    grp = g.create_group('data')
    grp.create_dataset('magnitudes/ascc', data = ascc)
    grp.create_dataset('magnitudes/magnitude', data = m)
    
    grp.create_dataset('skytrans/idx', data=idx)
    grp.create_dataset('skytrans/lstseq', data=lstseq)
    grp.create_dataset('skytrans/value', data=values)
    grp['skytrans'].attrs['grid'] = 'healpix'
    grp['skytrans'].attrs['nside'] = 8
    grp['skytrans'].attrs['lstrange'] = [0, tmax]

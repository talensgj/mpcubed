#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def delete_allnan(array):
    
    ind1, = np.where(np.all(np.isnan(array), axis=0))
    ind2, = np.where(np.all(np.isnan(array), axis=1))
    
    array = np.delete(array, ind1, axis=1)
    array = np.delete(array, ind2, axis=0)
    
    return array

with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value

hg = HealpixGrid(8)

with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter4.hdf5', 'r') as f:
    ms1 = f['data/magnitudes/m'].value
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s1 = f['data/skytrans/s'].value
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s1
s1 = tmp

with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter5.hdf5', 'r') as f:
    ms2 = f['data/magnitudes/m'].value
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s2 = f['data/skytrans/s'].value
   
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s2
s2 = tmp
    
# Magnitude calibration.
plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])
    
plt.subplot(gs[1])
plt.plot(vmag, ms1, '.', alpha=.2)
plt.plot(vmag, ms2, '.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.plot(ms1, ms1-ms2, '.', alpha=.2)
plt.ylim(-.5, .5)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()

# Sky profile.
s = s1 - s2
s = delete_allnan(s)

plt.figure(figsize=(13,12))

plt.suptitle('Sky', size='xx-large')

gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.1, vmax=.1)
plt.xlabel('Position')
plt.ylabel('Time')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('s')

plt.tight_layout()
plt.show()


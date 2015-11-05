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

with h5py.File('/data2/talens/3mEast/LBtests/June2.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value
    
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June2.hdf5', 'r') as f:
    ms = f['data/magnitudes/m'].value
    
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s = f['data/skytrans/s'].value

hg = HealpixGrid(8)
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp

skyidx = hg.find_gridpoint(ra, dec)

plt.plot(skyidx, ms-vmag, '.', alpha=.2)
plt.show()
    
# Magnitude calibration.
plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])
    
plt.subplot(gs[1])
plt.plot(vmag, ms, '.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.plot(dec, ms-vmag, '.', alpha=.2)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()

# Sky profile.
s = delete_allnan(s)

plt.figure(figsize=(13,12))

plt.suptitle('Sky', size='xx-large')

gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.xlabel('Position')
plt.ylabel('Time')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('s')

plt.tight_layout()
plt.show()

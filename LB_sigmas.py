#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from fLCfile import fLCfile
from core.coordinate_grids import PolarGrid, HealpixGrid

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
from viridis import viridis 

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
        
with h5py.File('/data2/talens/3mEast/LBtests/June1.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value
        
hg = HealpixGrid(8)

with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1.hdf5', 'r') as f:
    mz = f['data/magnitudes/m'].value
    sigma1 = f['data/magnitudes/sigma'].value
    
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s = f['data/skytrans/s'].value
    sigma2 = f['data/skytrans/sigma'].value

tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp

tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = sigma2
sigma2 = tmp

# Magnitude calibration.
plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])
    
plt.subplot(gs[1])
plt.errorbar(vmag, mz, yerr=sigma1, fmt='.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.errorbar(dec, mz-vmag, yerr=sigma1, fmt='.', alpha=.2)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()
    
# Sky profile.
s = delete_allnan(s)
sigma2 = delete_allnan(sigma2)

plt.figure(figsize=(13,12))

plt.suptitle('Sky', size='xx-large')

gs = gridspec.GridSpec(2, 4, height_ratios = [.5,10], width_ratios=[10,1,10,1])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow(s.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.xlabel('Position')
plt.ylabel('Time')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('s')

ax3 = plt.subplot(gs[1,2], xticks=[], yticks=[], sharex=ax1, sharey=ax1)
im = plt.imshow(sigma2.T, aspect='auto', cmap=viridis, vmin=-.1, vmax=.1)
plt.xlabel('Position')

ax4 = plt.subplot(gs[1,3])
plt.colorbar(im, cax=ax4).set_label('sigma')

plt.tight_layout()
plt.show()
    
    

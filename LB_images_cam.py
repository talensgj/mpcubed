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

with h5py.File('/data2/talens/3mEast/LBtests/June1.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/camip_June1_iter5.hdf5', 'r') as f:
    mz = f['data/magnitudes/m'].value
    
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value
    
pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, idx1, np.nan)
z = z - np.nanmean(z)

pg = PolarGrid(270, 720)
a = pg.put_values_on_grid(a, idx2, np.nan)
b = pg.put_values_on_grid(b, idx2, np.nan)
c = pg.put_values_on_grid(c, idx2, np.nan)
d = pg.put_values_on_grid(d, idx2, np.nan)

# Magnitude calibration.
plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])
    
plt.subplot(gs[1])
plt.plot(vmag, mz, '.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.plot(dec, mz-vmag, '.', alpha=.2)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()

# Transmission profile.
#z = z.reshape((13502, 722))
z = delete_allnan(z)

plt.figure(figsize=(13,12))

plt.suptitle('Transmission', size='xx-large')

gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=1.5)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('z')

plt.tight_layout()
plt.show()

# Intrapixel variations.
Ax = np.sqrt(a**2 + b**2)
Ay = np.sqrt(c**2 + d**2)
phix = np.arctan2(b, a)
phiy = np.arctan2(d, c)

Ax = delete_allnan(Ax)
Ay = delete_allnan(Ay)
phix = delete_allnan(phix)
phiy = delete_allnan(phiy)

plt.figure(figsize=(13,12))

plt.suptitle('Intrapixel Variations', size='xx-large')

gs = gridspec.GridSpec(3, 3, height_ratios = [.5,10,10], width_ratios=[10,10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
plt.title('x position', size='x-large')
plt.imshow(Ax.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)
plt.ylabel('Declination')

plt.subplot(gs[1,1], sharex=ax1, sharey=ax1)
plt.title('y position', size='x-large')
im = plt.imshow(Ay.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)

ax2 = plt.subplot(gs[1,2])
plt.colorbar(im, cax=ax2).set_label('Amplitude')

plt.subplot(gs[2,0], sharex=ax1, sharey=ax1)
plt.imshow(phix.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.subplot(gs[2,1], sharex=ax1, sharey=ax1)
im = plt.imshow(phiy.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.xlabel('Hour Angle')

ax3 = plt.subplot(gs[2,2])
plt.colorbar(im, cax=ax3).set_label('Phase')

plt.tight_layout()
plt.show()

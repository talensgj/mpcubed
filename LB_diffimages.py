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

with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value

pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter4.hdf5', 'r') as f:
    mz1 = f['data/magnitudes/m'].value

    ipx = f['data/camtrans/idx'].value
    z1 = f['data/camtrans/z'].value
    z1 = pg1.put_values_on_grid(z1, ipx, np.nan)

    ipx = f['data/intrapix/idx'].value
    a1 = f['data/intrapix/a'].value
    a1 = pg2.put_values_on_grid(a1, ipx, np.nan)
    b1 = f['data/intrapix/b'].value
    b1 = pg2.put_values_on_grid(b1, ipx, np.nan)
    c1 = f['data/intrapix/c'].value
    c1 = pg2.put_values_on_grid(c1, ipx, np.nan)
    d1 = f['data/intrapix/d'].value
    d1 = pg2.put_values_on_grid(d1, ipx, np.nan)

Ax1 = np.sqrt(a1**2 + b1**2)
Ay1 = np.sqrt(c1**2 + d1**2)
px1 = np.arctan2(b1, a1)
py1 = np.arctan2(d1, c1)

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5.hdf5', 'r') as f:
    mz2 = f['data/magnitudes/m'].value

    ipx = f['data/camtrans/idx'].value
    z2 = f['data/camtrans/z'].value
    z2 = pg1.put_values_on_grid(z2, ipx, np.nan)
    
    ipx = f['data/intrapix/idx'].value
    a2 = f['data/intrapix/a'].value
    a2 = pg2.put_values_on_grid(a2, ipx, np.nan)
    b2 = f['data/intrapix/b'].value
    b2 = pg2.put_values_on_grid(b2, ipx, np.nan)
    c2 = f['data/intrapix/c'].value
    c2 = pg2.put_values_on_grid(c2, ipx, np.nan)
    d2 = f['data/intrapix/d'].value
    d2 = pg2.put_values_on_grid(d2, ipx, np.nan)

Ax2 = np.sqrt(a2**2 + b2**2)
Ay2 = np.sqrt(c2**2 + d2**2)
px2 = np.arctan2(b2, a2)
py2 = np.arctan2(d2, c2)

plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])

plt.subplot(gs[1])
plt.plot(vmag, mz1, '.', alpha=.2)
plt.plot(vmag, mz2, '.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.plot(mz1, mz1-mz2, '.', alpha=.2)
plt.ylim(-.5, .5)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()

plt.figure(figsize=(13,12))

plt.suptitle('Transmission', size='xx-large')

gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow((z1-z2).T, aspect='auto', cmap=viridis, vmin=-.01, vmax=.01)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('z')

plt.tight_layout()
plt.show()

plt.figure(figsize=(13,12))

plt.suptitle('Intrapixel Variations', size='xx-large')

gs = gridspec.GridSpec(3, 3, height_ratios = [.5,10,10], width_ratios=[10,10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
plt.title('x position', size='x-large')
plt.imshow((Ax1-Ax2).T, aspect='auto', cmap=viridis, vmin=-.01, vmax=.01)
plt.ylabel('Declination')

plt.subplot(gs[1,1], sharex=ax1, sharey=ax1)
plt.title('y position', size='x-large')
im = plt.imshow((Ay1-Ay2).T, aspect='auto', cmap=viridis, vmin=-.01, vmax=.01)

ax2 = plt.subplot(gs[1,2])
plt.colorbar(im, cax=ax2).set_label('Amplitude')

plt.subplot(gs[2,0], sharex=ax1, sharey=ax1)
plt.imshow(np.sin(px1-px2).T, aspect='auto', cmap=viridis, vmin=-.1, vmax=.1)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.subplot(gs[2,1], sharex=ax1, sharey=ax1)
im = plt.imshow(np.sin(py1-py2).T, aspect='auto', cmap=viridis, vmin=-.1, vmax=.1)
plt.xlabel('Hour Angle')

ax3 = plt.subplot(gs[2,2])
plt.colorbar(im, cax=ax3).set_label('Phase')

plt.tight_layout()
plt.show()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

from coordinate_grids import PolarGrid

from matplotlib.ticker import FuncFormatter, MultipleLocator, NullLocator

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPS.hdf5') as f:
    
    ascc = f['header_table/ascc'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value
    
    pg = PolarGrid(13500, 720)
    decidx = pg.find_decidx(dec)
    
    here = decidx==360
    ascc = ascc[here]
    dec = dec[here]
    ra = ra[here]
    
    
    array = np.full((13502, len(ascc)), fill_value=np.nan)
    for i in range(len(ascc)):
        lc= f['data/'+ascc[i]]
        
        ha = np.mod(lc['lst']*15.-ra[i], 360.)
        haidx = pg.find_raidx(ha)
        
        array[haidx, i] = lc['flux0']
        
sort = np.argsort(dec)
array = array[:,sort]
array = array[1:-1]
array = np.roll(array, 13500/2, axis=0)
array = array/np.nanmedian(array, axis=0, keepdims=True)

xlim, ylim = np.where(np.isfinite(array))

xlim = pg.find_ra(xlim+1)/360*24-12

def periodic(x, pos):
    'The two args are the value and tick position'
    x = x%24
    x1 = np.floor(x)
    x2 = (x-np.floor(x))*60
    
    return '%02.fh%02.fm' % (x1, x2)

formatter = FuncFormatter(periodic)
majorLocator   = MultipleLocator(1./2)
minorLocator   = MultipleLocator(1./6)

plt.figure(figsize=(16,8))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_major_locator(NullLocator())
plt.imshow(array.T, aspect='auto', cmap=viridis, vmin=0, vmax=1.5, extent=(-12,12,0,1))
plt.xlabel('Hour Angle')
plt.xlim(np.amin(xlim)-12./13500, np.amax(xlim)+12./13500)
#plt.colorbar()
plt.tight_layout()
plt.show()

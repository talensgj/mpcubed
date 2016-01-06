#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
import healpy
from viridis import viridis

from core.index_functions import index_statistics
from core.coordinate_grids import HealpixGrid

with h5py.File('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 'r') as f:
    hdr = f['header_table']
    ra = hdr['ra'].value
    dec = hdr['dec'].value
    print np.amin(dec), np.amax(dec)
    print np.amin(ra), np.amax(ra)
    ax = plt.subplot(221)
    
    # Create the indices.
    hg = HealpixGrid(8)
    bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
    count = index_statistics(binnum, ra, statistic='count', keeplength=False)
    hmap = hg.put_values_on_grid(count, bins, fill_value=np.nan)
    print np.amax(count)
    healpy.mollview(hmap, min=1, max=200, sub=222, title='nside = 8', cmap=viridis, unit='# Stars')
    healpy.graticule(dpar = 10., dmer = 15.)
    ax.hist(count, bins = np.linspace(0.5, 200.5, 51), normed=True, histtype='step', label='nside = 8')
    
    hg = HealpixGrid(16)
    bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
    count = index_statistics(binnum, ra, statistic='count', keeplength=False)
    hmap = hg.put_values_on_grid(count, bins, fill_value=np.nan)
    
    healpy.mollview(hmap, min=1, max=150, sub=223, title='nside = 16', cmap=viridis, unit='# Stars')
    healpy.graticule(dpar = 10., dmer = 15.)
    ax.hist(count, bins = np.linspace(0.5, 150.5, 51), normed=True, histtype='step', label='nside = 16')
    
    hg = HealpixGrid(32)
    bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
    count = index_statistics(binnum, ra, statistic='count', keeplength=False)
    hmap = hg.put_values_on_grid(count, bins, fill_value=np.nan)
    
    healpy.mollview(hmap, min=1, max=150, sub=224, title='nside = 32', cmap=viridis, unit='# Stars')
    healpy.graticule(dpar = 10., dmer = 15.)
    ax.hist(count, bins = np.linspace(0.5, 150.5, 51), normed=True, histtype='step', label='nside = 32')
    
    plt.sca(ax)
    plt.title('2015-07-16 LPC')
    plt.legend()
    plt.xlim(.5,150.5)
    plt.xlabel('# Stars')
    
    plt.tight_layout()
    plt.show()

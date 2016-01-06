#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
import healpy
from viridis import viridis

from index_functions import index_statistics
from coordinate_grids import HealpixGrid

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', 'r') as f:
    hdr = f['table_header']
    dec = hdr['dec'].value
    
    plt.figure(figsize=(8,3))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    
    # Create the indices.
    bins = np.linspace(-90, 90, 361)
    cen = (bins[:-1]+bins[1:])/2
    bins = np.digitize(dec, bins)
    count = index_statistics(bins, dec, statistic='count', keeplength=False)
    
    ax2.plot(cen[np.unique(bins)], count, '.')
    ax1.hist(count, bins = np.linspace(0.5, 200.5, 41), normed=True, histtype='step', label='ny = 360')
    
    bins = np.linspace(-90, 90, 721)
    cen = (bins[:-1]+bins[1:])/2
    bins = np.digitize(dec, bins)
    count = index_statistics(bins, dec, statistic='count', keeplength=False)
    
    ax2.plot(cen[np.unique(bins)], count, '.')
    ax1.hist(count, bins = np.linspace(0.5, 200.5, 41), normed=True, histtype='step', label='ny = 720')
    
    bins = np.linspace(-90, 90, 1441)
    cen = (bins[:-1]+bins[1:])/2
    bins = np.digitize(dec, bins)
    count = index_statistics(bins, dec, statistic='count', keeplength=False)
    
    ax2.plot(cen[np.unique(bins)], count, '.')
    ax1.hist(count, bins = np.linspace(0.5, 200.5, 41), normed=True, histtype='step', label='ny = 1440')
    
    plt.sca(ax1)
    plt.title('2015-07-16 LPC')
    plt.legend()
    plt.xlim(.5,200.5)
    plt.xlabel('# Stars')
    
    plt.sca(ax2)
    plt.xlabel('Dec [deg]')
    plt.ylabel('# Stars')
    
    plt.tight_layout()
    plt.show()
    

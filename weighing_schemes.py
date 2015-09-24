#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

from index_functions import index_statistics

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150716LPC.hdf5', 'r') as g:
    
    ascc = f['header_table/ascc'].value
    vmag = f['header_table/vmag'].value
    nobs = f['header_table/nobs'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    plt.subplot(111, aspect='equal')
    for i in range(0,len(ascc),10):
        
        lc = f['data/'+ascc[i]]
        rc = g['data/'+ascc[i]]
        
        #scatter1 = index_statistics(lc['lstidx']//50, lc['eflux0'], statistic='mean')
        #scatter2 = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std')
        
        scatter1 = index_statistics(lc['lstidx']//50, lc['eflux0'], statistic='mean')
        scatter2 = index_statistics(lc['lstidx']//50, rc['ipcflux0'], statistic='std')
        
        plt.semilogy([vmag[i]]*len(scatter1), scatter2/scatter1, '.', c='k', alpha=.2)
        
    plt.axhline(1, c='r')
    plt.xlabel('V [mag]')
    plt.ylabel('mean(eflux0)/std(ipcflux0)')
    plt.xlim(2, 8.4)
    plt.ylim(1e-3,1e2)
    plt.show()
    plt.close()

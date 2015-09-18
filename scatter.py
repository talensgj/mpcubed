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
    
    #here = (dec>40)&(dec<46)
    here = (vmag<6.5)&(vmag>5.5)
    ascc = ascc[here]
    vmag = vmag[here]
    nobs = nobs[here]
    ra = ra[here]
    dec = dec[here]
    
    error = np.zeros(len(ascc))
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]].value
        rc1 = g['data/'+ascc[i]].value
        rc2 = g['data2/'+ascc[i]].value
        
        xbin = index_statistics(lc['lstidx']//50, lc['x'], statistic='mean')
        ybin = index_statistics(lc['lstidx']//50, lc['y'], statistic='mean')
        ebin = index_statistics(lc['lstidx']//50, rc1['cflux0'], statistic='std')/index_statistics(lc['lstidx']//50, rc1['cflux0'], statistic='mean')
        
        try:
            plt.scatter(xbin, ybin, c=np.log10(ebin), cmap=viridis)
        except: pass
        
        error[i] = np.nanstd(rc1['cflux0'])/np.nanmean(rc1['cflux0'])
        
plt.xlim(0,4008)
plt.ylim(0,2672)
plt.show()
        
ax = plt.subplot(221)
plt.semilogy(vmag, error, '.')
plt.subplot(222, sharey=ax)
plt.semilogy(nobs, error, '.')
plt.subplot(223, sharey=ax)
plt.semilogy(ra, error, '.')
plt.subplot(224, sharey=ax)
plt.semilogy(dec, error, '.')
plt.show()



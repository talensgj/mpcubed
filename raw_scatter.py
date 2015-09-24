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
    
    scatter1 = np.zeros(len(ascc))
    scatter2 = np.zeros(len(ascc))
    scatter3 = np.zeros(len(ascc))
    scatter4 = np.zeros(len(ascc))
    for i in range(len(ascc)):
    
        lc = f['data/'+ascc[i]].value
        rc1 = g['data/'+ascc[i]].value
        rc2 = g['data2/'+ascc[i]].value
        
        here = (lc['flag'] < 1) & (lc['flux0'] > 0) & (lc['sky'] > 0)# & (rc2['flags'] < 1)
        
        tmp = index_statistics(lc['lstidx'][here]//50, lc['flux0'][here], statistic='std')/index_statistics(lc['lstidx'][here]//50, lc['flux0'][here], statistic='mean')
        scatter1[i] = np.nanmedian(tmp)
        
        tmp = index_statistics(lc['lstidx'][here]//50, rc1['cflux0'][here], statistic='std')/index_statistics(lc['lstidx'][here]//50, rc1['cflux0'][here], statistic='mean')
        scatter2[i] = np.nanmedian(tmp)
        
        tmp = index_statistics(lc['lstidx'][here]//50, rc1['ipcflux0'][here], statistic='std')/index_statistics(lc['lstidx'][here]//50, rc1['ipcflux0'][here], statistic='mean')
        scatter3[i] = np.nanmedian(tmp)
        
        tmp = index_statistics(lc['lstidx'][here]//50, rc2['scflux0'][here], statistic='std')/index_statistics(lc['lstidx'][here]//50, rc2['scflux0'][here], statistic='mean')
        scatter4[i] = np.nanmedian(tmp)

plt.figure(figsize=(16,8))
ax = plt.subplot(221)
plt.semilogy(vmag, scatter1, '.', alpha=.2, label='flux0')
plt.semilogy(vmag, scatter4, '.', alpha=.2, label='scflux0')
plt.legend()
plt.xlabel('V [mag]')
plt.ylabel('std')
plt.ylim(1e-3, 1)
ax1 = plt.subplot(222, sharex=ax)
plt.semilogy(vmag, scatter2/scatter1, '.', alpha=.2)
plt.xlabel('V [mag]')
plt.ylabel('cflux0/flux0')
plt.subplot(223, sharex=ax1, sharey=ax1)
plt.xlabel('V [mag]')
plt.ylabel('ipcflux0/cflux0')
plt.semilogy(vmag, scatter3/scatter2, '.', alpha=.2)
plt.subplot(224, sharex=ax1, sharey=ax1)
plt.semilogy(vmag, scatter4/scatter3, '.', alpha=.2)
plt.ylim(1e-1,10)
plt.xlim(2, 8.4)
plt.xlabel('V [mag]')
plt.ylabel('scflux0/ipcflux0')
plt.tight_layout()
plt.show()



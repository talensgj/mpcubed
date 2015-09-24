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
    #here = (vmag<6.5)&(vmag>5.5)
    #ascc = ascc[here]
    #vmag = vmag[here]
    #nobs = nobs[here]
    #ra = ra[here]
    #dec = dec[here]
    
    v = np.array([])
    error = np.array([])
    
    for i in range(len(ascc)):
    
        lc = f['data/'+ascc[i]].value
        rc1 = g['data/'+ascc[i]].value
        rc2 = g['data2/'+ascc[i]].value
        
        if ascc[i] == '690978':
            plt.plot(lc['jdmid'], lc['flux0'], '.')
            plt.plot(lc['jdmid'], rc1['cflux0'], '.')
            plt.plot(lc['jdmid'], rc2['scflux0'], '.')
            plt.show()
        
        here = (lc['flag']<1)&(rc1['flags']<1)
        
        count = index_statistics(lc['lstidx'][here]//50, lc['x'][here], statistic='count')
        #xbin = index_statistics(lc['lstidx']//50, lc['x'], statistic='mean')
        #ybin = index_statistics(lc['lstidx']//50, lc['y'], statistic='mean')
        ebin = index_statistics(lc['lstidx'][here]//50, rc1['cflux0'][here], statistic='std')/index_statistics(lc['lstidx'][here]//50, rc1['cflux0'][here], statistic='mean')
        
        here = count == 50
        
        if np.any(ebin[here]>1e2): print ascc[i]
        
        v = np.append(v, [vmag[i]]*sum(here))
        error = np.append(error, ebin[here])
        
ax = plt.subplot(221)
plt.semilogy(v, error, '.')
#plt.subplot(222, sharey=ax)
#plt.semilogy(nobs, error, '.')
#plt.subplot(223, sharey=ax)
#plt.semilogy(ra, error, '.')
#plt.subplot(224, sharey=ax)
#plt.semilogy(dec, error, '.')
plt.show()



#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

from index_functions import index_statistics

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPE.hdf5') as f:
    
    ascc = f['table_header/ascc'].value
    
    for i in range(0, len(ascc), 50):
        lc = f['data/'+ascc[i]]
        
        
        idx = (lc['lstidx']//50).astype('int')
        sort = np.argsort(idx)
        isort = np.arange(len(idx))[sort]
        
        
        ym = index_statistics(idx, lc['y'], statistic='mean', keeplength=True)
        xm = index_statistics(idx, lc['x'], statistic='mean', keeplength=True)
        
        idx, uni = np.unique(lc['lstidx']//50, '.')

        here = np.zeros(len(lc['lstidx']), dtype='bool')
        for i in range(len(idx)):
            these = lc['lstidx']//50 == idx[i]
            
            if np.all(np.diff(lc['y'][these])%1 < 5e-2)&(np.sum(these)>25):
                here[these] = True

        plt.subplot(111, aspect='equal')
        plt.plot(lc['y'][here], lc['x'][here], '.', c='k')
        plt.ylim(4008,0)
        plt.xlim(2672,0)
        
        
plt.show()

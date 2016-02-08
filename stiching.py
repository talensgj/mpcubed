#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

import filters

data = ['/data2/talens/2015Q2/LPN/red0_201504ALPN.hdf5',
        '/data2/talens/2015Q2/LPE/red0_201504ALPE.hdf5',
        '/data2/talens/2015Q2/LPS/red0_201504ALPS.hdf5',
        '/data2/talens/2015Q2/LPW/red0_201504ALPW.hdf5',
        '/data2/talens/2015Q2/LPC/red0_201504ALPC.hdf5']

ASCC = np.array([])
NOBS = np.array([])

for filename in data:
    with h5py.File(filename, 'r') as f:
        grp = f['header_table']
        ascc = grp['ascc'].value
        nobs = grp['nobs'].value
        
        ASCC = np.append(ASCC, ascc)
        NOBS = np.append(NOBS, nobs)
        
ascc, idx = np.unique(ASCC, return_inverse=True)
nobs = np.bincount(idx, NOBS)

ascc = ascc[nobs > 1000]
for i in range(len(ascc)):
    
    for filename in data:
        with h5py.File(filename, 'r') as f:
        
            try:
                lc = f['data/'+ascc[i]].value
            except:
                continue
              
            jdmid = lc['jdmid']
            lst = lc['lst']
            mag0 = lc['mag0']
            emag0 = lc['emag0']
            nobs = lc['nobs']  
    
            emag0 = emag0/np.sqrt(nobs)
            
            select = (nobs == 50)
            jdmid = jdmid[select]
            lst = lst[select]
            mag0 = mag0[select]
            emag0 = emag0[select]
            
            weights = 1/emag0**2
            
            base = np.ptp(jdmid)
            chisq, pars, fit = filters.harmonic(jdmid, lst, mag0, weights, 2*base, 5)
            mag0 = mag0 - fit
                
            avg = np.mean(mag0)
            std = np.std(mag0)
            
            select = (np.abs(mag0 - avg) < 5*std)
                
        plt.errorbar(jdmid[select], mag0[select], yerr=emag0[select], fmt='.')
    plt.show()
        

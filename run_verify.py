#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

boxlstsq = '/data2/talens/2015Q2/boxlstsq/bls0_2015Q2_patch266.hdf5'
with h5py.File(boxlstsq, 'r') as f:
    
    grp = f['header']
    ascc = grp['ascc'].value
    chisq0 = grp['chisq0'].value
    
    grp = f['data']
    freq = grp['freq'].value
    dchisq = grp['dchisq'].value
    #hchisq = grp['hchisq'].value
    depth = grp['depth'].value
    
    grp = f['baseline']
    jdmid = grp['jdmid'].value
    mask = grp['mask'].value
    
nobs = np.sum(np.abs(mask - 1), axis=1)
flag = np.zeros(len(ascc), dtype='int')
    
# Best fit.
arg = np.argmax(dchisq, axis=0)
best_freq = freq[arg]
best_chisq = chisq0 - dchisq[arg,np.arange(len(ascc))]

# Good fit?
quality = best_chisq/nobs
args, = np.where(quality > 4)
flag[args] = flag[args] + 1
    
# Anti-transit ratio.
tmp = dchisq*np.sign(depth)
ratio = -np.amax(tmp, axis=0)/np.amin(tmp, axis=0)
args, = np.where(ratio < 1.5)
flag[args] = flag[args] + 2

# Data gaps.
q = (1.8/24)*best_freq**(2./3) # Fractional transit duration.
phase = np.outer(jdmid, best_freq)
phase = np.mod(phase, 1)
phase = np.sort(phase, axis=0)
gapsizes = np.diff(phase, axis=0)
gapsizes = np.vstack([gapsizes, 1 - np.ptp(phase, axis=0)])
gapsizes = np.amax(gapsizes, axis=0)/q
args, = np.where(gapsizes > 2.5)
flag[args] = flag[args] + 4 

args, = np.where(flag == 0)
for i in args:
    
    plt.title('ASCC {}'.format(ascc[i]))
    plt.plot(freq, dchisq[:,i])
    plt.axvline(best_freq[i], c='k')
    plt.show()
    plt.close()


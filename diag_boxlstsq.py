#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np 

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/2015Q2/bls0_2015Q2_patch266.hdf5', 'r') as f:
    
    grp = f['header']
    ascc = grp['ascc'].value
    
    grp = f['data']
    freq = grp['freq'].value
    dchisq = grp['dchisq'].value
    
nstars = len(ascc)

for i in range(nstars):
    
    plt.title('ASCC {}'.format(ascc[i]))
    plt.plot(freq, dchisq[:,i])
    plt.xlabel(r'Freq [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
    plt.show()
    plt.close()

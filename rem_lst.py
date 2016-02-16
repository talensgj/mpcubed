#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from scipy.interpolate import UnivariateSpline

import matplotlib.pyplot as plt

import filters

filename = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'

with h5py.File(filename, 'r') as f:
    grp = f['data/441706']
    jdmid = grp['jdmid'].value
    lst = grp['lst'].value
    mag = grp['mag0'].value
    emag = grp['emag0'].value
    nobs = grp['nobs'].value
    
emag = emag/np.sqrt(nobs)

select = (nobs == 50)
jdmid = jdmid[select]
lst = lst[select]
mag = mag[select]
emag = emag[select]

base = np.ptp(lst)
print 2*base

chisq, pars, fit1 = filters.harmonic(lst, mag, 1/emag**2, 24., 5)
print chisq

chisq, pars, fit2 = filters.harmonic(lst, mag, 1/emag**2, 2*base, 5)
print chisq
  
chisq, pars, fit3 = filters.harmonic(lst, mag, 1/emag**2, 3.33, 10)
print chisq
    
plt.subplot(211)
plt.errorbar(lst, mag, yerr=emag, fmt='.')
plt.plot(lst, fit3, '.')

plt.subplot(212)
plt.errorbar(lst, mag - fit3, yerr=emag, fmt='.')

plt.show()

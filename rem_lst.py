#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from scipy.interpolate import UnivariateSpline

import matplotlib.pyplot as plt

import filters
from fourierfuncs import fourier_mat

def median_filter(lst, mag, window):
    
    res = np.zeros(len(lst))
    for i in range(len(lst)):
        sel = np.abs(lst - lst[i]) < window/2
        res[i] = np.nanmedian(mag[sel])
    
    return res

filename = '/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5'

with h5py.File(filename, 'r') as f:
    grp = f['data/29188']
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

chisq, pars, fit1 = filters.harmonic(lst, mag, 1/emag**2, 24., 6)
print chisq

fit2 = median_filter(lst, mag, 4.)

plt.subplot(311)
plt.errorbar(lst, mag, yerr=emag, fmt='.')
plt.plot(lst, fit1, '.', label='Fourier')
plt.plot(lst, fit2, '.', label='Median')
plt.legend()

plt.subplot(312)
plt.errorbar(lst, mag - fit1, yerr=emag, fmt='.')
plt.errorbar(lst, mag - fit2, yerr=emag, fmt='.')

plt.subplot(313)
plt.plot(lst, fit1 - fit2, '.')
plt.ylim(-.005, .005)

plt.show()

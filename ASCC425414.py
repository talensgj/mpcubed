#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from scipy.signal import lombscargle

import matplotlib.pyplot as plt

import fourierfuncs

P = .5666

with h5py.File('/data2/talens/inj_signals/reference/red0_2015Q2LPE.hdf5', 'r') as f:
    grp = f['data/425414']
    jdmid = grp['jdmid'].value
    lst = grp['lst'].value
    mag0 = grp['mag0'].value
    emag0 = grp['emag0'].value
    nobs = grp['nobs'].value
    
emag0 = emag0/np.sqrt(nobs)
    
select = (nobs == 50) & (emag0 < .05)
jdmid = jdmid[select]
lst = lst[select]
mag0 = mag0[select]
emag0 = emag0[select]

#periods = np.logspace(1., -1, 5000)
#print periods[0], periods[-1]
#freqs = 1/periods
#print freqs[0], freqs[-1]
#pgram = lombscargle(jdmid, (mag0 - fit1).astype('float64'), 2*np.pi*freqs)

#plt.plot(freqs, np.sqrt(4*pgram/len(mag0)), '.')
#plt.axvline(1/P)
#plt.show()

#arg = np.argmax(pgram)
#P = 1/freqs[arg]

fit2 = np.zeros(len(jdmid))
for niter in range(50):
    
    pars1, fit1 = fourierfuncs.fourier_fit(jdmid, mag0 - fit2, 1/P, 20, 1/emag0**2)
    print pars1

    pars2, fit2 = fourierfuncs.fourier_fit(lst, mag0 - fit1, 1./6, 5, 1/emag0**2)
    print pars2
    
    plt.subplot(311)
    plt.plot(jdmid, mag0, '.')
    plt.plot(jdmid, fit1, '.')
    
    plt.subplot(312)
    plt.plot(jdmid, mag0 - fit1, '.')
    plt.plot(jdmid, fit2, '.')
    
    plt.subplot(313)
    plt.plot(jdmid, mag0 - fit1 - fit2, '.')
    
    plt.show()
    plt.close()


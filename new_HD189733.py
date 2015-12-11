#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as  np

from core.coarse_decor import coarse_decor
from core.BLS_ms import BLS
from scipy.signal import lombscargle

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

P = 2.21857312
Tp = 2454037.612

with h5py.File('/data2/talens/2015Q2/LPE/test.hdf5', 'r') as f:
    
    jdmid = f['data/807144/jdmid'].value
    lstseq = f['data/807144/lstseq'].value
    nobs = f['data/807144/nobs'].value
    mag0 = f['data/807144/mag0'].value
    emag0 = f['data/807144/emag0'].value
    clouds0 = f['data/807144/clouds0'].value

here = (nobs == 50)
jdmid = jdmid[here]
lstseq = lstseq[here]
mag0 = mag0[here]
emag0 = emag0[here]
clouds0 = clouds0[here]

emag0 = emag0/np.sqrt(50.)

plt.hist(emag0, bins = np.linspace(0, .2, 21))
plt.show()

here = (emag0 < .01)
jdmid = jdmid[here]
lstseq = lstseq[here]
mag0 = mag0[here]
emag0 = emag0[here]
clouds0 = clouds0[here]

plt.hist(clouds0, bins = np.linspace(-1, 1, 21))
plt.show()

here = (clouds0 < .05)
jdmid = jdmid[here]
lstseq = lstseq[here]
mag0 = mag0[here]
emag0 = emag0[here]
clouds0 = clouds0[here]

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1, idx2, mag0, emag0)

plt.plot(idx2, mag0, '.')
plt.plot(idx2, a[idx1] + b[idx2], '.')
plt.show()

mag0 = mag0 - (a[idx1] + b[idx2])

phase = (jdmid - Tp)/P
phase = np.mod(phase+.5, 1)-.5

plt.figure(figsize=(18, 5))
plt.plot(phase, mag0, '.')
plt.ylim(.1, -.1)
plt.xlim(-.5, .5)
plt.xlabel('Phase')
plt.ylabel(r'$\Delta m$')
plt.tight_layout()
plt.show()

freq, dchisq = BLS(jdmid, mag0, emag0)

plt.plot(freq, dchisq)
plt.axvline(1/P, c='k')
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()

plt.semilogy(freq, dchisq)
plt.axvline(1/P)
plt.show()

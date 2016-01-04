#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as  np

from core.coarse_decor import coarse_decor
from core.BLS_ms import BLS
from freq_fit import remove_freq

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

with h5py.File('/data2/talens/2015Q2/LPE/new_solver/red0_2015Q2LPE.hdf5', 'r') as f:
    
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
pars, fit = remove_freq(jdmid, mag0 - a[idx1], 1/.99727, 3, weights=1/emag0**2)

ax = plt.subplot(311)
plt.plot(idx2, mag0 - a[idx1], '.')
plt.plot(idx2, b[idx2], '.')

plt.subplot(312, sharex=ax, sharey=ax)
plt.plot(idx2, mag0 - a[idx1], '.')
plt.plot(idx2, fit, '.')

plt.subplot(313, sharex=ax)
plt.plot(idx2, b[idx2] - fit, '.')
plt.ylim(-.02, .02)

plt.show()

mag0 = mag0 - (a[idx1] + fit)

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

freq, dchisq, depth, hchisq = BLS(jdmid, mag0, emag0)

plt.figure(figsize=(18, 10))
plt.plot(freq, dchisq)
plt.axvline(1/P, c='k')
plt.ylim(0, 1200)
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()

plt.figure(figsize=(18, 10))
plt.plot(freq, depth)
plt.axvline(1/P, c='k')
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()

plt.figure(figsize=(18, 10))
plt.plot(freq, hchisq)
plt.axvline(1/P, c='k')
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()

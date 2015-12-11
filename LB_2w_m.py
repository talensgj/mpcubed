#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams

from core.coarse_decor import coarse_decor

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

P = 2.21857312
Tp = 2454037.612

with h5py.File('/data2/talens/2015Q2/LPE/nosigmas/tmp0_201506LPE.hdf5', 'r') as f:
    lc = f['807144'].value

with h5py.File('/data2/talens/2015Q2/LPE/nosigmas/tmp0_201506ALPE.hdf5', 'r') as f:
    lcA = f['807144'].value
    
with h5py.File('/data2/talens/2015Q2/LPE/nosigmas/tmp0_201506BLPE.hdf5', 'r') as f:
    lcB = f['807144'].value

lc = lc[lc['nobs']==50]
lcA = lcA[lcA['nobs']==50]
lcB = lcB[lcB['nobs']==50]

plt.subplot(211)
plt.plot(lc['jdmid'], lc['mag0'], '.')
plt.plot(lcA['jdmid'], lcA['mag0'], '.')
plt.plot(lcB['jdmid'], lcB['mag0'], '.')
plt.subplot(212)
plt.plot(lcA['jdmid'], lcA['mag0'] - lc['mag0'][:len(lcA['mag0'])], '.')
plt.plot(lcB['jdmid'], lcB['mag0'] - lc['mag0'][len(lcA['mag0']):], '.')

plt.show()



phase = (lc['jdmid'] - Tp)/P
phase = np.mod(phase+.5, 1)-.5

phaseA = (lcA['jdmid'] - Tp)/P
phaseA = np.mod(phaseA+.5, 1)-.5

phaseB = (lcB['jdmid'] - Tp)/P
phaseB = np.mod(phaseB+.5, 1)-.5

idx1 = lc['lstseq'] // 270
idx2 = lc['lstseq'] % 270
a, b, niter, chisq, npoints, npars = coarse_decor(idx1, idx2, lc['mag0'], lc['emag0'])
corr = a[idx1] + b[idx2]

idx1 = lcA['lstseq'] // 270
idx2 = lcA['lstseq'] % 270
a, b, niter, chisq, npoints, npars = coarse_decor(idx1, idx2, lcA['mag0'], lcA['emag0'])
corrA = a[idx1] + b[idx2]

idx1 = lcB['lstseq'] // 270
idx2 = lcB['lstseq'] % 270
a, b, niter, chisq, npoints, npars = coarse_decor(idx1, idx2, lcB['mag0'], lcB['emag0'])
corrB = a[idx1] + b[idx2]

plt.subplot(211)
plt.plot(phase, lc['mag0'] - corr, '.')
plt.plot(phaseA, lcA['mag0'] - corrA, '.')
plt.plot(phaseB, lcB['mag0'] - corrB, '.')

plt.subplot(212)
y = lc['mag0'] - corr
plt.plot(phaseA, lcA['mag0'] - corrA - y[:len(phaseA)], '.')
plt.plot(phaseB, lcB['mag0'] - corrB - y[len(phaseA):], '.')

plt.show()

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

x = np.array([])
y = np.array([])
err = np.array([])

plt.figure(figsize=(18, 5))

with h5py.File('/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5', 'r') as f:
    
    jdmid = f['data/807144/jdmid'].value
    lstseq = f['data/807144/lstseq'].value
    nobs = f['data/807144/nobs'].value
    mag0 = f['data/807144/mag0'].value
    emag0 = f['data/807144/emag0'].value
    clouds0 = f['data/807144/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)
phase = (jdmid - Tp)/P
phase = np.mod(phase+.5, 1)-.5
plt.errorbar(phase[here], mag0[here] - a[idx1[here]] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

x = np.append(x, jdmid[here])
y = np.append(y, mag0[here] - a[idx1[here]] - fit) 
err = np.append(err, emag0[here]/np.sqrt(nobs[here]))

with h5py.File('/data2/talens/2015Q2/LPE/new_solver/red0_2015Q2LPE.hdf5', 'r') as f:
    
    jdmid = f['data/807144/jdmid'].value
    lstseq = f['data/807144/lstseq'].value
    nobs = f['data/807144/nobs'].value
    mag0 = f['data/807144/mag0'].value
    emag0 = f['data/807144/emag0'].value
    clouds0 = f['data/807144/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)
phase = (jdmid - Tp)/P
phase = np.mod(phase+.5, 1)-.5
plt.errorbar(phase[here], mag0[here] - a[idx1[here]] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

x = np.append(x, jdmid[here])
y = np.append(y, mag0[here] - a[idx1[here]] - fit) 
err = np.append(err, emag0[here]/np.sqrt(nobs[here]))

with h5py.File('/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5', 'r') as f:
    
    jdmid = f['data/807144/jdmid'].value
    lstseq = f['data/807144/lstseq'].value
    nobs = f['data/807144/nobs'].value
    mag0 = f['data/807144/mag0'].value
    emag0 = f['data/807144/emag0'].value
    clouds0 = f['data/807144/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)
phase = (jdmid - Tp)/P
phase = np.mod(phase+.5, 1)-.5
plt.errorbar(phase[here], mag0[here] - a[idx1[here]] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

x = np.append(x, jdmid[here])
y = np.append(y, mag0[here] - a[idx1[here]] - fit) 
err = np.append(err, emag0[here]/np.sqrt(nobs[here]))

plt.ylim(.1, -.1)
plt.xlim(-.5, .5)
plt.xlabel('Phase')
plt.ylabel(r'$\Delta m$')
plt.tight_layout()
plt.show()

freq, dchisq, depth, hchisq = BLS(x, y, err)

plt.figure(figsize=(18, 10))
plt.plot(freq, dchisq)
plt.axvline(1/P, c='k')
plt.ylim(0, 1200)
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()


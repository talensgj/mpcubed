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

x = np.array([])
y = np.array([])
err = np.array([])

plt.figure(figsize=(18, 5))

with h5py.File('/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5', 'r') as f:
    
    jdmid = f['data/344324/jdmid'].value
    lstseq = f['data/344324/lstseq'].value
    nobs = f['data/344324/nobs'].value
    mag0 = f['data/344324/mag0'].value
    emag0 = f['data/344324/emag0'].value
    clouds0 = f['data/344324/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)
print pars
plt.errorbar(jdmid[here], mag0[here] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

with h5py.File('/data2/talens/2015Q2/LPE/new_solver/red0_2015Q2LPE.hdf5', 'r') as f:
    
    jdmid = f['data/344324/jdmid'].value
    lstseq = f['data/344324/lstseq'].value
    nobs = f['data/344324/nobs'].value
    mag0 = f['data/344324/mag0'].value
    emag0 = f['data/344324/emag0'].value
    clouds0 = f['data/344324/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)

plt.errorbar(jdmid[here], mag0[here] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

with h5py.File('/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5', 'r') as f:
    
    jdmid = f['data/344324/jdmid'].value
    lstseq = f['data/344324/lstseq'].value
    nobs = f['data/344324/nobs'].value
    mag0 = f['data/344324/mag0'].value
    emag0 = f['data/344324/emag0'].value
    clouds0 = f['data/344324/clouds0'].value

here = (nobs == 50) & (emag0 < .35) & (clouds0 <.05)

idx1 = lstseq // 270
idx2 = lstseq % 270

a, b, niter, chisq, npoint, npars = coarse_decor(idx1[here], idx2[here], mag0[here], emag0[here])
pars, fit = remove_freq(jdmid[here], mag0[here] - a[idx1[here]], 1/.99727, n=5)

plt.errorbar(jdmid[here], mag0[here] - fit, yerr = emag0[here]/np.sqrt(nobs[here]), fmt='.')

plt.xlabel('Time [JD]')
plt.ylabel('Magnitude')
plt.tight_layout()
plt.show()



#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from core import coarse_decor

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

P = 2.21857312
Tp = 2454037.612

filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
filelist = np.sort(filelist)

jdmid = np.array([])
lstidx = np.array([])
flux0 = np.array([])
eflux0 = np.array([])
x = np.array([])
y = np.array([])
flags = np.array([])

for filename in filelist:
    with h5py.File(filename, 'r') as f:
        try: lc = f['data/807144'].value
        except: pass
        else:
            jdmid = np.append(jdmid, lc['jdmid'])
            lstidx = np.append(lstidx, lc['lstidx'])
            flux0 = np.append(flux0, lc['flux0'])
            eflux0 = np.append(eflux0, lc['eflux0'])
            x = np.append(x, lc['x'])
            y = np.append(y, lc['y'])
            flags = np.append(flags, lc['flag'])
            
lstidx = lstidx.astype('int')
dayidx = np.floor(jdmid)
dayidx = dayidx - np.amin(dayidx)
dayidx = dayidx.astype('int')

here = (flux0 > 0) & (eflux0 > 0) & (flags < 1)
jdmid = jdmid[here]
dayidx = dayidx[here]
lstidx = lstidx[here]
flux0 = flux0[here]
eflux0 = eflux0[here]
x = x[here]
y = y[here]

mag0 = -2.5*np.log10(flux0)
emag0 = 2.5/np.log(10.)*eflux0/flux0

m, z, a, b, c, d, sigma, niter, chisq, npoints, npars = coarse_decor.coarse_positions(dayidx, lstidx, lstidx//50, x, y, mag0, emag0, verbose=True)

badday, = np.where(sigma > 0)
for i in badday:
    here = dayidx == i
    plt.plot(jdmid[here], mag0[here], '.')
    plt.show()

plt.subplot(311)
plt.plot(z, '.')
plt.subplot(312)
plt.plot(np.sqrt(a**2 + b**2), '.')
plt.plot(np.sqrt(c**2 + d**2), '.')
plt.subplot(313)
plt.plot(np.arctan2(b, a), '.')
plt.plot(np.arctan2(d, c), '.')
plt.show()

sol1 = m[dayidx] + z[lstidx]
sol2 = a[lstidx//50]*np.sin(2*np.pi*y) + b[lstidx//50]*np.cos(2*np.pi*y)
sol3 = c[lstidx//50]*np.sin(2*np.pi*x) + d[lstidx//50]*np.cos(2*np.pi*x)

plt.subplot(211)
plt.plot(lstidx, mag0, '.')
plt.plot(lstidx, sol1, '.')

plt.subplot(212)
plt.plot(lstidx, mag0 - sol1, '.')
plt.plot(lstidx, sol2, '.')
plt.plot(lstidx, sol3, '.')
plt.show()

mag0 = mag0 - sol1 - sol2 - sol3

binidx = np.ravel_multi_index((dayidx, lstidx//50), (31, 270))
count = np.bincount(binidx)
bin_jdmid = np.bincount(binidx, jdmid)/np.bincount(binidx)
bin_mag0 = np.bincount(binidx, mag0)/np.bincount(binidx)
bin_jdmid = bin_jdmid[count>49]
bin_mag0 = bin_mag0[count>49]

bin_phase = (bin_jdmid - Tp)/P
bin_phase = bin_phase - np.amin(np.floor(bin_phase))

bin_phase = bin_phase%1

plt.plot(bin_phase, bin_mag0, '.')
plt.show()


            
            
            
            
            

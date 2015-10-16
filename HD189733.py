#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

import TransitModels

from core import index_functions
from core import BLS_ms
from core import sysrem

P = 2.21857312
Tp = 2454037.612

filelist = glob.glob('/data2/talens/3mEast/fLC_201507??LPE.hdf5')
filelist = np.sort(filelist)

filelist2 = glob.glob('/data2/talens/3mEast/red_201507??LPE.hdf5')
filelist2 = np.sort(filelist2)

lst = np.array([])
lstidx = np.array([])
jdmid = np.array([])
mag = np.array([])
emag = np.array([])
flags = np.array([])

for filename, filename2 in zip(filelist, filelist2):
    print filename, filename2
    with h5py.File(filename, 'r') as f, h5py.File(filename2, 'r') as g:
    
        try:
            lc = f['data/807144'].value
        except:
            continue
            
        rc1 = g['data/807144'].value
        rc2 = g['data2/807144'].value
    
    lst = np.append(lst, lc['lst'])
    lstidx = np.append(lstidx, lc['lstidx'])
    jdmid = np.append(jdmid, lc['jdmid'])
    mag = np.append(mag, rc2['sipc_mag0'])
    emag = np.append(emag, rc1['emag0'])
    flags = np.append(flags, lc['flag']+rc1['flags']+rc2['flags'])

lstidx = lstidx.astype('int')

dayidx = np.floor(jdmid)
dayidx = dayidx-np.amin(dayidx)
dayidx = dayidx.astype('int')

binidx = np.ravel_multi_index((dayidx, lstidx//50), (120, 270))

here = (flags < 1)
lst = lst[here]
lstidx = lstidx[here]
jdmid = jdmid[here]
mag = mag[here]
emag = emag[here]
dayidx = dayidx[here]
binidx = binidx[here]

m = np.bincount(dayidx, mag/emag**2)/np.bincount(dayidx, 1./emag**2)
mag = mag - m[dayidx]

plt.plot(jdmid, mag, '.')
plt.show()

plt.plot(lst, mag, '.')
plt.show()

plt.plot(dayidx, mag, '.')
plt.show()

count = index_functions.index_statistics(binidx, binidx, statistic='count')
bin_jd = index_functions.index_statistics(binidx, jdmid, statistic='mean')
bin_lst = index_functions.index_statistics(binidx, lst, statistic='mean')
bin_mag = index_functions.index_statistics(binidx, mag, statistic='mean')
bin_emag = index_functions.index_statistics(binidx, mag, statistic='std')/np.sqrt(count)

plt.figure(figsize=(16,8))
plt.subplot(211)
plt.plot(bin_lst, bin_mag, '.', alpha=.1)
plt.ylim(-.1, .1)
plt.ylabel(r'$\Delta m$')

for i in range(5):
    a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lstidx, mag, emag)
    print chisq
    mag = mag - a1[dayidx]*a2[lstidx]

count = index_functions.index_statistics(binidx, binidx, statistic='count')
bin_jd = index_functions.index_statistics(binidx, jdmid, statistic='mean')
bin_lst = index_functions.index_statistics(binidx, lst, statistic='mean')
bin_mag = index_functions.index_statistics(binidx, mag, statistic='mean')
bin_emag = index_functions.index_statistics(binidx, mag, statistic='std')/np.sqrt(count)

plt.subplot(212)
plt.plot(bin_lst, bin_mag, '.', alpha=.1)
plt.ylim(-.1, .1)
plt.xlabel('LST [idx]')
plt.ylabel(r'$\Delta m$')
plt.show()

here = (count == 50)
bin_jd = bin_jd[here]
bin_lst = bin_lst[here]
bin_mag = bin_mag[here]
bin_emag = bin_emag[here]

plt.plot(bin_jd, bin_mag, '.')
plt.show()

phase = np.mod(bin_jd-Tp, P)/P
phase = np.mod(phase+.5, 1)

time = np.linspace(Tp, Tp+P, 600)
model = TransitModels.MA_Model(time, Tp, P, 0.0995*1.138/0.805, 215.1*0.03142/0.805, 0.0041, 85.51/180*np.pi, np.pi/2., 0., 0.)
mphase = np.mod((time-Tp)/P, 1)
model[np.abs(mphase)<.2] = 1.
model[np.abs(mphase-1)<.2] = 1.

plt.figure(figsize=(16,5))
plt.errorbar(phase-.5, bin_mag, yerr=bin_emag, fmt='.')
plt.plot(mphase-.5, model-1)
plt.xlim(-.5, .5)
plt.ylim(.1, -.1)
plt.xlabel('Phase')
plt.ylabel(r'$\Delta m$')
plt.show()

freq, chisq = BLS_ms.BLS(bin_jd, bin_mag, bin_emag)

plt.plot(freq, chisq, c='k')
plt.axvline(1/P, c='g')
plt.axvline(freq[np.argmax(chisq)], c='r')
plt.xlabel(r'Frequency [day$^{-1}$]')
plt.ylabel(r'$\Delta \chi^2$')
plt.show()

P = 1/freq[np.argmax(chisq)]
print chisq[np.argmax(chisq)]
print P
phase = np.mod(bin_jd/P, 1.)

plt.figure(figsize=(16,5))
plt.errorbar(phase, bin_mag, yerr=bin_emag, fmt='.')
plt.xlim(0, 1)
plt.ylim(.1, -.1)
plt.xlabel('Unreferenced Phase')
plt.ylabel(r'$\Delta m$')
plt.show()



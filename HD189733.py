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

filelist = glob.glob('/data2/talens/3mEast/fLC_20150???LPE.hdf5')
filelist = np.sort(filelist)

filelist2 = glob.glob('/data2/talens/3mEast/master_reduction/red_20150???LPE.hdf5')
filelist2 = np.sort(filelist2)

lst = np.array([])
jdmid = np.array([])
bjdmid = np.array([])
mag = np.array([])
emag = np.array([])

for filename, filename2 in zip(filelist, filelist2):
    print filename, filename2
    with h5py.File(filename, 'r') as f, h5py.File(filename2, 'r') as g:
    
        try:
            lc = f['data/807144'].value
        except:
            continue
        rc1 = g['data/807144'].value
        rc = g['data2/807144'].value
    
    here = (lc['flag'] < 1) & (rc['flags'] < 1) & np.isfinite(rc['sipc_mag0'])
    
    tmp_jd = lc['jdmid']
    tmp_lst = lc['lstidx'].astype('int')
    tmp_mag = np.nanmedian(rc['sipc_mag0'])-rc['sipc_mag0']
    tmp_emag = rc1['emag0']
    
    tmp_jd = tmp_jd[here]
    tmp_mag = tmp_mag[here]
    tmp_lst = tmp_lst[here]
    tmp_emag = tmp_emag[here]
    
    lst = np.append(lst, tmp_lst)
    jdmid = np.append(jdmid, tmp_jd)
    mag = np.append(mag, tmp_mag)
    emag = np.append(emag, tmp_emag)

dayidx = np.floor(jdmid)
dayidx = dayidx-np.amin(dayidx)
dayidx = dayidx.astype('int')
lst = lst.astype('int')

a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lst, mag, emag)

mag = mag - a1[dayidx]*a2[lst]

plt.plot(a2, '.')
plt.show()

a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lst, mag, emag)

mag = mag - a1[dayidx]*a2[lst]

plt.plot(a2, '.')
plt.show()

a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lst, mag, emag)

mag = mag - a1[dayidx]*a2[lst]

plt.plot(a2, '.')
plt.show()

a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lst, mag, emag)

mag = mag - a1[dayidx]*a2[lst]

plt.plot(a2, '.')
plt.show()

a1, a2, niter, chisq, npoints, npars = sysrem.sysrem(dayidx, lst, mag, emag)

mag = mag - a1[dayidx]*a2[lst]

plt.plot(a2, '.')
plt.show()

binidx = np.ravel_multi_index((dayidx, lst//50), (120, 270))

count = index_functions.index_statistics(binidx, binidx, statistic='count')
bin_jd = index_functions.index_statistics(binidx, jdmid, statistic='mean')
bin_mag = index_functions.index_statistics(binidx, mag, statistic='mean')
bin_emag = index_functions.index_statistics(binidx, mag, statistic='std')/np.sqrt(count)

here = count==50
bin_jd = bin_jd[here]
bin_mag = bin_mag[here]
bin_emag = bin_emag[here]

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
plt.ylim(-.1, .1)
plt.xlabel('Phase')
plt.ylabel(r'$\Delta m$')
plt.show()

freq, chisq = BLS_ms.BLS(bin_jd, bin_mag, bin_emag)

plt.plot(freq, chisq)
plt.axvline(1/P)
plt.show()

P = 1/freq[np.argmax(chisq)]

phase = np.mod(bin_jd, P)/P
phase = np.mod(phase+.5, 1)

plt.errorbar(phase, bin_mag, yerr=bin_emag, fmt='.')
plt.ylim(-.1, .1)
plt.show()



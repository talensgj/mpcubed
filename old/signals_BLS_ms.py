#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit
from fourierfuncs import fourier_fit
from BLS_ms import BLS_ms

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value

JDMID = np.array([])
MAG0 = np.array([])
EMAG0 = np.array([])
LSTSEQ = np.array([])
STARIDX = np.array([], dtype='int')

nstars = 500

flag = np.zeros(nstars)
for i in range(nstars):

    with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        jdmid = grp['jdmid'].value
        mag0 = grp['mag0'].value
        emag0 = grp['emag0'].value
        nobs = grp['nobs'].value
        lst = grp['lst'].value
        lstseq = grp['lstseq'].value
    
    emag0 = emag0/np.sqrt(nobs)

    BL = np.ptp(jdmid)
    if (P[i] > BL/9.):
        flag[i] = 1
        print 'Warning transit may not be detectable.'

    select = (nobs == 50) & (emag0 < .05)
    jdmid = jdmid[select]
    mag0 = mag0[select]
    emag0 = emag0[select]
    lst = lst[select]
    lstseq = lstseq[select]
    
    # Remove lst variations.
    pars, fit = fourier_fit(lst, mag0, 1/6., 5, 1/(emag0**2))
    mag0 = mag0 - fit
    
    JDMID = np.append(JDMID, jdmid)
    MAG0 = np.append(MAG0, mag0)
    EMAG0 = np.append(EMAG0, emag0)
    LSTSEQ = np.append(LSTSEQ, lstseq)
    STARIDX = np.append(STARIDX, [i]*len(mag0))
    
lstseq, args, idx = np.unique(LSTSEQ, return_index=True, return_inverse=True)
jdmid = JDMID[args]
ntimes = len(lstseq)

mag0 = np.full((nstars, ntimes), fill_value = np.nan)
mag0[STARIDX, idx] = MAG0

emag0 = np.full((nstars, ntimes), fill_value = np.nan)
emag0[STARIDX, idx] = EMAG0

freq, dchisq, depth, hchisq = BLS_ms(jdmid, mag0.T, emag0.T)

args = np.argmax(dchisq, axis=0)
Prec = 1/freq[args]

print P
print Prec

here = flag==0

plt.subplot(111, aspect='equal')
plt.plot(P[here], Prec[here], ls='.', c='g')
plt.plot(P[here], Prec[here], ls='.', c='r')
plt.plot(P, P, c='k')
plt.plot(P, .5*P, c='k')
plt.plot(P, 2*P, c='k')
plt.xlim(0, 15)
plt.ylim(0, 30)
plt.savefig('/home/talens/BLS_ms.png')

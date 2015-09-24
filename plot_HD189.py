#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np
import h5py

from sysrem import sysrem
from index_functions import index_statistics

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

from BLS_ms import BLS

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

period = 2.21857312
t0 = 2453988.80336

def box_model(jd, t0, period, q):
    
    phase = np.mod(jd-t0, period)/period
    phase = (phase-.5)%1
    phase = phase-np.floor(phase)
    
    here = np.abs(phase-.5) < .5*q 
    
    sol = np.ones(len(jd))
    sol[here] = .98
    
    return sol

ar1 = np.array([])
ar2 = np.array([])
ar3 = np.array([])

for cam in ['LPC', 'LPE', 'LPW']:

    filelist1 = glob.glob('/data2/talens/Jul2015/fLC_201507??%s.hdf5'%cam)
    filelist2 = glob.glob('/data2/talens/Jul2015/red_201507??%s.hdf5'%cam)
    filelist1 = np.sort(filelist1)
    filelist2 = np.sort(filelist2)

    jdmid = np.array([])
    lst = np.array([])
    lstidx = np.array([])
    flux = np.array([])
    eflux = np.array([])
    sky = np.array([])
    flags1 = np.array([])
    flags2 = np.array([])
    flags3 = np.array([])

    for i in range(len(filelist1)):
        
        with h5py.File(filelist1[i], 'r') as f, h5py.File(filelist2[i], 'r') as g:
            
            try:
                lc = f['data/807144'].value
                
            except:
                continue
            
            rc1 = g['data/807144'].value
            rc2 = g['data2/807144'].value
            
            jdmid = np.append(jdmid, lc['jdmid'])
            lst = np.append(lst, lc['lst'])
            lstidx = np.append(lstidx, lc['lstidx'])
            flux = np.append(flux, rc2['scflux0'])
            eflux = np.append(eflux, rc2['escflux0'])
            sky = np.append(sky, lc['flux1']/lc['flux0'])
            flags1 = np.append(flags1, lc['flag'])
            flags2 = np.append(flags2, rc1['flags'])
            flags3 = np.append(flags3, rc2['flags'])

    lstidx = lstidx.astype('int')
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - np.amin(dayidx)

    here = (flags1 < 1) & (flags2 < 1) & (flags3 < 1)
    
    jdmid = jdmid[here]
    lstidx = lstidx[here]
    dayidx = dayidx[here]
    flux = flux[here]
    eflux = eflux[here]

    sol = np.zeros(len(flux))
    for i in range(2):
        a1, a2 = sysrem(dayidx, lstidx, flux-sol, eflux)[:2]

        sol += a1[dayidx]*a2[lstidx]
    
    flux /= sol
    eflux /= sol

    binidx = np.ravel_multi_index((dayidx, lstidx//50), (31, 270))

    count = index_statistics(binidx, flux, statistic='count')
    bin_jd = index_statistics(binidx, jdmid, statistic='mean')
    bin_flux = index_statistics(binidx, flux, statistic='mean')
    bin_eflux = index_statistics(binidx, flux, statistic='std')/np.sqrt(count)

    here = count == 50

    ar1 = np.append(ar1, bin_jd[here])
    ar2 = np.append(ar2, bin_flux[here])
    ar3 = np.append(ar3, bin_eflux[here])

    plt.plot(bin_jd[here], bin_flux[here], '.')
    plt.plot(bin_jd[here], box_model(bin_jd[here], t0, period, .01))
    
plt.show()


freq, chisq = BLS(ar1, ar2, ar3)
arg = np.argmax(chisq)
print 1/freq[arg]

plt.subplot(211)
phase = np.mod(ar1-t0, period)/period
plt.errorbar(phase, ar2, yerr=ar3, fmt='.')
plt.subplot(212)
plt.plot(freq, chisq)
plt.axvline(1/period, c='k')
plt.axvline(2/period, c='k')
plt.axvline(3/period, c='k')
plt.axvline(1/(2*period), c='k')
plt.axvline(1/(3*period), c='k')
plt.axvline(1/(4*period), c='k')
plt.show()




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

def transit_model(t, t0, period, eta, delta, c=10.):
    
    t_p = period*np.sin(np.pi*(t-t0)/period)/(np.pi*eta)
    
    model = .5*delta*(2-np.tanh(c*(t_p+.5))+np.tanh(c*(t_p-.5)))+(1-delta)
    
    return model

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

period = 2.21857312
t0 = 2453988.80336

filelist1 = glob.glob('/data2/talens/Jul2015/fLC_201507??LPC.hdf5')
filelist2 = glob.glob('/data2/talens/Jul2015/red_201507??LPC.hdf5')
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
dayidx = dayidx[here]
lst = lst[here]
lstidx = lstidx[here]
flux = flux[here]
sky = sky[here]
eflux = eflux[here]

sol = np.zeros(len(flux))
for i in range(1):
    F, T = sysrem(dayidx, lstidx, flux-sol, eflux)[:2]
    sol += F[dayidx]*T[lstidx]
    
cflux = flux/sol
ecflux = eflux/sol

plt.subplot(311)
plt.plot(lstidx, flux, '.')
plt.subplot(312)
plt.plot(T, '.')
plt.subplot(313)
plt.plot(lstidx, cflux, '.')

plt.show()


binidx = np.ravel_multi_index((dayidx, lstidx//50), (np.amax(dayidx)+1, 270))

count = index_statistics(binidx, cflux, statistic='count')
bin_lst = index_statistics(binidx, lst, statistic='mean')
bin_jd = index_statistics(binidx, jdmid, statistic='mean')
bin_cflux = index_statistics(binidx, cflux, statistic='mean')
bin_ecflux = index_statistics(binidx, cflux, statistic='std')/np.sqrt(count)

ax = plt.subplot(211)
plt.plot(jdmid, cflux, '.', c='k', zorder=0)
plt.scatter(bin_jd, bin_cflux, zorder=1, c='r')
plt.plot(jdmid, sky, '.')
plt.subplot(212, sharex=ax)
plt.plot(jdmid, transit_model(jdmid, t0, period, .086, .02), '.')
plt.plot(bin_jd, transit_model(bin_jd, t0, period, .086, .02), '.')
plt.show()

here = count == 50

plt.errorbar(bin_jd[here], bin_cflux[here], yerr=bin_ecflux[here], fmt='o')
plt.plot(bin_jd, transit_model(bin_jd, t0, period, .086, .02), '.')
plt.ylim(.95, 1.05)
plt.show()
exit()
here = bin_ecflux > 0

freq, chisq = BLS(jdmid, cflux, ecflux)
plt.plot(freq, chisq, c='b')
freq, chisq = BLS(bin_jd[here], bin_cflux[here], bin_ecflux[here])
plt.plot(freq, chisq, c='r')
plt.axvline(1/period, c='k')
plt.show()




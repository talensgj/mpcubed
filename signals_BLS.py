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
from fourierfuncs import fourier_fit, lst_trend, fourier_fit2
from boxlstsq import boxlstsq

def filter1(lst, mag0, emag0):
    
    pars, fit = fourier_fit(lst, mag0, 1/6., 5, 1/(emag0**2))

    mag0 = mag0 - fit
    
    return mag0
        
def filter2(dayidx, lst, mag0, emag0):
    
    pars1, pars2, fit = lst_trend(dayidx, lst, mag0, 1/6., 5, 1/(emag0**2))

    mag0 = mag0 - fit
    
    return mag0

def filter3(jdmid, lst, mag0, emag0):
        
    fit2 = np.zeros(len(mag0))
    for i in range(10):
        pars1, fit1 = fourier_fit(lst, mag0 - fit2, 1/6., 5, 1/(emag0**2))
        pars2, fit2 = fourier_fit2(jdmid, mag0 - fit1, 1/(2*np.ptp(jdmid)), 10, 1/(emag0**2))
        
    mag0 = mag0 - fit1 - fit2

    return mag0

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value
 
nstars = 500
Prec = np.zeros(nstars)
flag = np.zeros(nstars)

for i in range(nstars):
    
    with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        jdmid = grp['jdmid'].value
        mag0 = grp['mag0'].value
        emag0 = grp['emag0'].value
        nobs = grp['nobs'].value
        lst = grp['lst'].value

    emag0 = emag0/np.sqrt(nobs)

    select = (nobs == 50) & (emag0 < .05)
    jdmid = jdmid[select]
    mag0 = mag0[select]
    emag0 = emag0[select]
    lst = lst[select]
    
    BL = np.ptp(jdmid)
    if (P[i] > BL/9.):
        flag[i] = 1
        print 'Warning transit may not be detectable.'
    
    tmp = np.floor(jdmid).astype('int')
    dayidx = tmp - np.amin(tmp)
    
    # Remove lst variations.
    mag0 = filter3(jdmid, lst, mag0, emag0)
    
    # Compute the BLS
    print 'ITERATION', i
    freq, dchisq, depth, hchisq = boxlstsq(jdmid, mag0, emag0)
    
    arg = np.nanargmax(dchisq)
    Prec[i] = 1/freq[arg]

    # Plot the result.
    phase1 = np.mod(jdmid*freq[arg], 1)
    phase2 = np.mod(jdmid/P[i], 1)
    
    plt.figure(figsize = (16,8))
    
    plt.subplot(311)
    plt.title(r'ASCC {}, $V = {:.1f}$, $\delta = {:.3f}$, $P = {:.3f}$'.format(ascc[i], vmag[i], delta[i], P[i]))
    plt.plot(freq, dchisq, c='k')
    plt.axvline(1/P[i], c='g')
    plt.axvline(1/Prec[i], c='r')
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
    
    plt.subplot(312)
    plt.plot(phase1, -mag0, '.')
    plt.ylim(-.1, .1)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    plt.subplot(313)
    plt.plot(phase2, -mag0, '.')
    plt.ylim(-.1, .1)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    plt.tight_layout()
    plt.savefig('/data2/talens/inj_signals/signals/filter3/ASCC{}.png'.format(ascc[i]))
    #plt.show()
    plt.close()

# Plot periods.
P = P[:nstars]
here = (flag == 0)
x = np.array([1,15])

plt.subplot(111)
plt.scatter(P[here], Prec[here], c='g')
plt.scatter(P[~here], Prec[~here], c='r')
plt.plot(x, x, c='k')
plt.plot(x, .5*x, c='k')
plt.plot(x, 2*x, c='k')
plt.xlim(0, 15)
plt.ylim(0, 30)
plt.xlabel('P [days]')
plt.ylabel(r'P$_{\rm{rec}}$ [days]')
plt.tight_layout()
plt.savefig('/data2/talens/inj_signals/signals/filter3.png')
plt.show()

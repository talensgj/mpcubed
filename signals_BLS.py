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
from BLS_ms import BLS

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value
 
Prec = np.zeros(len(P))
flag = np.zeros(len(P))

for i in [214]:
    
    with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        jdmid = grp['jdmid'].value
        mag0 = grp['mag0'].value
        emag0 = grp['emag0'].value
        nobs = grp['nobs'].value
        lst = grp['lst'].value

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
    
    # Remove lst variations.
    pars, fit = fourier_fit(lst, mag0, 1/6., 5, 1/(emag0**2))
    mag0 = mag0 - fit
    
    ## Model.
    #time = np.linspace(0, P[i], 1000)
    #model = transit.softened_box_model(time, P[i], Tp[i], delta[i], eta[i])
    #model = model*2.5/np.log(10.)
    #mphase = np.mod((time - Tp[i])/P[i], 1)
    #mphase = np.mod(mphase + .5, 1) - .5
    
    #sort = np.argsort(mphase)
    #mphase = mphase[sort]
    #model = model[sort]
    
    #phase = np.mod((jdmid - Tp[i])/P[i], 1)
    #phase = np.mod(phase + .5, 1) - .5
    
    ## Plot the result.
    #plt.figure(figsize = (18, 5))
    
    #plt.subplot(111)
    #plt.title(r'ASCC {}, $V = {:.1f}$, $\delta = {:.3f}$, $P = {:.3f}$'.format(ascc[i], vmag[i], delta[i], P[i]))
    #plt.plot(phase, -mag0, '.')
    #plt.plot(mphase, model)
    #plt.xlim(-.5, .5)
    #plt.ylim(-.1, .1)
    #plt.xlabel('Phase')
    #plt.ylabel(r'$\Delta m$')
    
    #plt.tight_layout()
    #plt.show()
    ##plt.savefig('/data2/talens/inj_signals/signals/fourier5/ASCC{}.png'.format(ascc[i]))
    #plt.close()

    # Compute the BLS
    print 'ITERATION', i
    freq, dchisq, depth, hchisq = BLS(jdmid, mag0, emag0)
    
    arg = np.nanargmax(dchisq)
    Prec[i] = freq[arg]

here = flag==0

plt.subplot(111, aspect='equal')
plt.scatter(P, Prec, c=flag)
plt.plot(P[here], P[here], c='g')
plt.plot(P[here], P[here], c='r')
plt.plot(P, .5*P, c='k')
plt.plot(P, 2*P, c='k')
plt.xlim(0, 15)
plt.ylim(0, 30)
plt.show()

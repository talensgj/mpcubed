#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit
import fourierfuncs

def harmonics_filter(jdmid, lst, mag0, emag0, Pjd, njd, Plst=24., nlst=5):
    
    weights = 1/emag0**2
    
    fmat = np.ones((len(mag0),1))
    
    if (njd > 0):
        mat_jd = fourierfuncs.fourier_mat(jdmid, 1/Pjd, njd)
        fmat = np.hstack([fmat, mat_jd])
        
    if (nlst > 0):
        mat_lst = fourierfuncs.fourier_mat(lst, 1/Plst, nlst)
        fmat = np.hstack([fmat, mat_lst])
    
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), mag0*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    
    chisq = np.sum(weights*(mag0 - fit)**2)/(len(mag0) - len(pars))
    
    return chisq, pars, fit

def avg_filter(jdmid, lstseq, mag0, emag0, Pjd, njd):
    
    weights = 1/emag0**2
    idx = lstseq%270
    
    fmat = np.ones((len(mag0),1))
    if (njd > 0):
        mat_jd = fourierfuncs.fourier_mat(jdmid, 1/Pjd, njd)
        fmat = np.hstack([fmat, mat_jd])
    
    fit2 = np.zeros(len(mag0))
    for i in range(10):
        
        pars1 = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), (mag0 - fit2)*np.sqrt(weights))[0]
        fit1 = np.dot(fmat, pars1)
        
        pars2 = np.bincount(idx, weights*(mag0 - fit1))/np.bincount(idx, weights)
        fit2 = pars2[idx]
        
    fit = fit1 + fit2
    chisq = np.sum(weights*(mag0 - fit)**2)/(len(mag0) - len(pars1) - len(np.unique(idx)))
        
    return chisq, fit

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value

nstars = len(ascc)
PARS = np.array([])
CHISQ = np.array([])

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

    select = (nobs == 50) & (emag0 < .05)
    jdmid = jdmid[select]
    mag0 = mag0[select]
    emag0 = emag0[select]
    lst = lst[select]
    lstseq = lstseq[select]
    
    if len(jdmid) < 50:
        print 'Few points.'
        continue
        
    #chisq, pars, fit = harmonics_filter(jdmid, lst, mag0, emag0, Pjd=2*np.ptp(jdmid), njd=5)
    #PARS = np.append(PARS, pars)
    #CHISQ = np.append(CHISQ, chisq)
    
    chisq, fit = avg_filter(jdmid, lstseq, mag0, emag0, Pjd=2*np.ptp(jdmid), njd=5)
    CHISQ = np.append(CHISQ, chisq)
    
    phase = np.mod((jdmid - Tp[i])/P[i], 1)
    phase = np.mod(phase + .5, 1) - .5
    
    time = np.linspace(0, P[i], 1000)
    model = transit.softened_box_model(time, P[i], Tp[i], delta[i], eta[i])
    model = model*2.5/np.log(10.)
    mphase = np.mod((time - Tp[i])/P[i], 1)
    mphase = np.mod(mphase + .5, 1) - .5
    
    sort = np.argsort(mphase)
    mphase = mphase[sort]
    model = model[sort]
    
    fig = plt.figure(figsize=(16,9))
    
    gs = gridspec.GridSpec(2, 1)
    
    ax1 = plt.subplot(gs[0])
    ax1.invert_yaxis()
    plt.title(r'ASCC {}, $V = {:.1f}$, $\delta = {:.3f}$, $P = {:.3f}$'.format(ascc[i], vmag[i], delta[i], P[i]))
    plt.errorbar(phase, mag0, yerr=emag0, fmt='.')
    plt.plot(phase, fit, '.')
    
    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.invert_yaxis()
    #plt.annotate(r'$\chi^2_\nu={:.2f}$'.format(chisq), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
    plt.errorbar(phase, mag0 - fit, yerr=emag0, fmt='.')
    plt.plot(mphase, -model)
    plt.xlim(-.5, .5)
    plt.ylim(.1, -.1)
    plt.xlabel('Phase')
    
    plt.tight_layout()
    #plt.show()
    plt.savefig('/data2/talens/inj_signals/signals/AF_njd5_2BL/ASCC{}.png'.format(ascc[i]))
    plt.close()
    
with h5py.File('AF_njd5_2BL.hdf5') as f:
    f.create_dataset('chisq', data=CHISQ)
    #f.create_dataset('pars', data=PARS)

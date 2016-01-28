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

import fourierfuncs

def filter1(lst, mag0, emag0):
    
    pars, fit = fourierfuncs.fourier_fit(lst, mag0, 1./24, 5, 1/emag0**2)

    return pars, fit
    
def filter2(jdmid, lst, mag0, emag0):
    
    n = np.floor(2*np.ptp(jdmid)/3.).astype('int')
    
    mat1 = fourierfuncs.fourier_mat(lst, 1/24., 5)
    mat2 = fourierfuncs.fourier_mat(jdmid, 1/(2*np.ptp(jdmid)), n)
    fmat = np.hstack([np.ones((len(mag0),1)), mat1, mat2])
    weights = 1/emag0**2
    
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), mag0*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    
    return pars, fit

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value

nstars = 500
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
    
    pars, fit = filter1(lst, mag0, emag0)
    chisq = (mag0 - fit)**2/emag0**2
    chisq = np.sum(chisq)/(len(mag0) - len(pars))
    
    fig = plt.figure(figsize=(16,9))
    
    gs = gridspec.GridSpec(2, 4)
    
    ax1 = plt.subplot(gs[0,0])
    plt.errorbar(lst, mag0, yerr=emag0, fmt='.')
    plt.plot(lst, fit, '.')
    
    ax2 = plt.subplot(gs[0,1], sharex=ax1)
    plt.plot(lst, mag0 - fit, '.')
    plt.annotate(r'$\chi^2_\nu = {:.2f}$'.format(chisq), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')

    pars, fit = filter2(jdmid, lst, mag0, emag0)
    chisq = (mag0 - fit)**2/emag0**2
    chisq = np.sum(chisq)/(len(mag0) - len(pars))
    
    plt.subplot(gs[0,2], sharex=ax1, sharey=ax1)
    plt.errorbar(lst, mag0, yerr=emag0, fmt='.')
    plt.plot(lst, fit, '.')
    
    plt.subplot(gs[0,3], sharex=ax1, sharey=ax2)
    plt.plot(lst, mag0 - fit, '.')
    plt.annotate(r'$\chi^2_\nu = {:.2f}$'.format(chisq), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
    
    plt.subplot(gs[1, :], sharey=ax1)
    plt.errorbar(jdmid, mag0, yerr=emag0, fmt='.')
    plt.plot(jdmid, fit, '.')
    
    plt.tight_layout()
    plt.show()

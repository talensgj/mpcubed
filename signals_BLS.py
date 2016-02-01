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
import fourierfuncs 
from boxlstsq import boxlstsq

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
Prec = np.zeros(nstars)
flag = np.zeros(nstars)
Dchisq = np.zeros(nstars)
Hchisq = np.zeros(nstars)
Depth = np.zeros(nstars)

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
        flag[i] = 1
        continue
    
    BL = np.ptp(jdmid)
    if (P[i] > BL/9.):
        flag[i] = 1
        print 'Warning transit may not be detectable.'
    
    # Remove lst variations.
    chisq, fit = avg_filter(jdmid, lstseq, mag0, emag0, Pjd=2*BL, njd=5)
    
    # Compute the BLS
    print 'ITERATION', i
    freq, dchisq, depth, hchisq = boxlstsq(jdmid, mag0 - fit, emag0)
    
    arg = np.nanargmax(dchisq)
    Prec[i] = 1/freq[arg]
    Dchisq[i] = dchisq[arg]
    Hchisq[i] = hchisq[arg]
    Depth[i] = depth[arg]

with h5py.File('AFperiods_njd5_2BL.hdf5') as f:
    f.create_dataset('Prec', data=Prec)
    f.create_dataset('flag', data=flag)
    f.create_dataset('dchisq', data=Dchisq)
    f.create_dataset('hchisq', data=Hchisq)
    f.create_dataset('depth', data=Depth)

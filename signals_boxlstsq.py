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

import filters
from boxlstsq import boxlstsq

# Read the list of stars with an injected signal.
with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    P = f['P'].value

# Process each star.
for i in range(len(ascc)):
    
    flag = 0
    
    with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        lstseq = grp['lstseq'].value
        lst = grp['lst'].value
        jdmid = grp['jdmid'].value
        mag0 = grp['mag0'].value
        emag0 = grp['emag0'].value
        nobs = grp['nobs'].value
        
    emag0 = emag0/np.sqrt(nobs)

    select = (nobs == 50)
    lstseq = lstseq[select]
    lst = lst[select]
    jdmid = jdmid[select]
    mag0 = mag0[select]
    emag0 = emag0[select]
    
    weights = 1/emag0**2
    
    nobs = len(jdmid)
    base = np.ptp(jdmid)
    
    # Flag the star.
    if (base/9 < P[i]):
        flag += 1
    
    if (nobs < 50):
        flag += 2
        continue
    
    # Remove lst variations.
    chisq, pars, fit = filters.harmonic(jdmid, lst, mag0, weights, 2*base, 5)
    mag0 = mag0 - fit
    
    # Compute the boxlstsq.
    freq, dchisq, depth, hchisq, chisq0 = boxlstsq(jdmid, mag0, weights)
    
    with h5py.File('signals_boxlstsq.hdf5') as f:
        grp = f.create_group(ascc)
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq)
        grp.create_dataset('hchisq', data=hchisq)
        grp.create_dataset('depth', data=depth)
        grp.attrs['nobs'] = nobs
        grp.attrs['base'] = base
        grp.attrs['flag'] = flag
        grp.attrs['chisq'] = chisq
        grp.attrs['pars'] = pars
        grp.attrs['chisq0'] = chisq0

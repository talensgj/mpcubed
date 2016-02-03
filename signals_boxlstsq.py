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

from boxlstsq import boxlstsq

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    vmag = f['vmag'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value

PARS = np.zeros((len(ascc), 13))
CHISQ = np.zeros(len(ascc))
for i in range(len(ascc)):
    
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
    
    if len(jdmid) < 50:
        continue
    
    # Remove lst variations.
    chisq, pars, fit = harmonics_filter(jdmid, lst, mag0, weights, 2*np.ptp(jdmid), 5)
    
    PARS[i] = pars
    CHISQ[i] = chisq
    
with h5py.File('HF_params.hdf5') as f:
    grp = f.create_group('nlst5min4')
    grp.create_dataset('chisq', data=CHISQ)
    grp.create_dataset('pars', data=PARS)

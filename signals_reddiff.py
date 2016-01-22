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

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value
    
for i in range(len(ascc)):
    
    with h5py.File('/data2/talens/inj_signals/reference/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        jdmid_ref = grp['jdmid'].value
        mag0_ref = grp['mag0'].value

    with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/' + ascc[i]]
        jdmid_inj = grp['jdmid'].value
        mag0_inj = grp['mag0'].value

    # Result.
    phase = np.mod((jdmid_ref - Tp[i])/P[i], 1)
    phase = np.mod(phase + .5, 1) - .5
    result = mag0_ref - mag0_inj
    
    # Model.
    time = np.linspace(0, P[i], 1000)
    model = transit.softened_box_model(time, P[i], Tp[i], delta[i], eta[i])
    model = model*2.5/np.log(10.)
    mphase = np.mod((time - Tp[i])/P[i], 1)
    mphase = np.mod(mphase + .5, 1) - .5
    
    sort = np.argsort(mphase)
    mphase = mphase[sort]
    model = model[sort]
    
    # Plot the result.
    plt.figure(figsize = (18, 5))
    
    plt.subplot(111)
    plt.title(r'ASCC {}, $\delta = {:.3f}$, $P = {:.3f}$'.format(ascc[i], delta[i], P[i]))
    plt.plot(phase, result, '.')
    plt.plot(mphase, model)
    plt.xlim(-.5, .5)
    plt.ylim(-.1, .1)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    plt.tight_layout()
    plt.savefig('/data2/talens/inj_signals/signals/red_diff/ASCC{}.png'.format(ascc[i]))
    plt.close()

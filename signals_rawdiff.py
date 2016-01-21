#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package.models import transit

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value
    
for i in range(len(ascc)):
    
    with h5py.File('/data2/talens/inj_signals/reference/fLC_201506ALPE.hdf5', 'r') as f:
        lc_ref = f['data/' + ascc[i]].value

    with h5py.File('/data2/talens/inj_signals/signals/fLC_201506ALPE.hdf5', 'r') as f:
        lc_inj = f['data/' + ascc[i]].value

    plt.plot(lc_ref['jdmid'], lc_inj['flux0']/lc_ref['flux0']-1, '.')
    plt.plot(lc_ref['jdmid'], transit.softened_box_model(lc_ref['jdmid'], P[i], Tp[i], delta[i], eta[i]))
    plt.show()

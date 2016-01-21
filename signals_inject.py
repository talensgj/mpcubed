#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import IO
from package.statistics import statistics
from package.models import transit

filelist = glob.glob('/data2/talens/inj_signals/signals/fLC_*.hdf5')

with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5', 'r') as f:
    ascc = f['ascc'].value
    delta = f['delta'].value
    P = f['P'].value
    Tp = f['Tp'].value
    eta = f['mu'].value
    
for filename in filelist:
    
    with h5py.File(filename) as f:
        grp = f['data']
    
        for i in range(len(ascc)):
            print i
            
            try:
                lc = grp[ascc[i]].value
            except:
                continue
            
            # Insert signal
            model = transit.softened_box_model(lc['jdmid'], P[i], Tp[i], delta[i], eta[i])
            model = model + 1
            
            lc['flux0'] = lc['flux0']*model
            lc['eflux0'] = lc['eflux0']*model
            
            # Save new lightcurve.
            del grp[ascc[i]]
            grp.create_dataset(ascc[i], data = lc)

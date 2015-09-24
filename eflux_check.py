#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5') as f:
    
    ascc = f['header_table/ascc'].value
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        
        error = np.sqrt(lc['flux0']+np.pi*(2.5)**2*lc['sky'])
        
        plt.plot(lc['eflux0'], error, '.')
        plt.show()

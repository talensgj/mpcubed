#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5', 'r') as f:
    
    ascc = f['data'].keys()
    
    for i in range(0,len(ascc),100):
        
        lc = f['data/'+ascc[i]]
        rc = f['data2/'+ascc[i]]
        
        plt.figure(figsize=(16,8))
        plt.subplot(311)
        plt.plot(lc['camtrans0'], '.')
        plt.plot(lc['intrapix0'], '.')
        plt.plot(rc['skytrans0'], '.')
        
        plt.subplot(312)
        plt.plot(lc['cflux0'], '.')
        plt.plot(lc['ipcflux0'], '.')
        plt.plot(rc['scflux0'], '.')
        
        plt.subplot(313)
        plt.plot(lc['flags'], '.')
        plt.plot(rc['flags'], '.')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        

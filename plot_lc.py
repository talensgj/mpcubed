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

with h5py.File('/data2/talens/Jul2015/coarsered_20150710LPC.hdf5', 'r') as f:
    
    ascc = f['data2'].keys()
    
    lc = f['data2/807144']
        
    plt.plot(lc['scflux0'], '.')
    plt.show()
    plt.close()
    
    for i in range(len(ascc)):
        
        lc = f['data2/'+ascc[i]]
        
        plt.plot(lc['scflux0'], '.')
        plt.show()
        plt.close()

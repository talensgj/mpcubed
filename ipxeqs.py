#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/3mEast/fLC_20150716LPE.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value

    for i in range(0, len(ascc), 50):
        
        lc = f['data/'+ascc[i]].value
        
        flux = lc['flux0']
        x = lc['x']
        y = lc['y']

        dx = x - (np.floor(x) + .5) 
        dy = y - (np.floor(y) + .5)
        
        r = np.sqrt(dx**2 + dy**2)

        ax = plt.subplot(411)
        plt.plot(flux, '.')
        plt.subplot(412, sharex=ax)
        plt.plot(np.sin(2*np.pi*dx), '.')
        plt.subplot(413, sharex=ax)
        plt.plot(np.sin(2*np.pi*dy), '.')
        plt.subplot(414, sharex=ax)
        plt.plot(np.sin(2*np.pi*r), '.')
        plt.show()
        plt.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'


with h5py.File('/data2/talens/3mEast/fLC_20150611LPE.hdf5') as f, h5py.File('/data2/talens/3mEast/red_std_w_20150611LPE.hdf5') as g:
    
    ascc = f['header_table/ascc'].value
    dec = f['header_table/dec'].value
    print dec[ascc=='807144']
    
    #here = (dec > 22.5) & (dec < 23)
    ascc = ['805488']
    
    for i in range(len(ascc)):
        
        try: lc = f['data/'+ascc[i]].value
        except: continue
        rc = g['data/'+ascc[i]].value
        
        plt.figure(figsize=(16,8))
        ax = plt.subplot(311)
        plt.title('ASCC '+ascc[i])
        plt.plot(rc['mag0'], '.')
        plt.subplot(312, sharex=ax)
        plt.plot(rc['camtrans0']+rc['intrapix0'], '.')
        plt.subplot(313, sharex=ax)
        plt.plot(rc['ipc_mag0']-np.nanmean(rc['ipc_mag0']), '.')
        plt.ylim(-.25, .25)
        plt.show()
        plt.close()

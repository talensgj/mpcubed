#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np
import h5py

from index_functions import index_statistics

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

period = 2.21857312
t0 = 2453988.80336

jdmid0 = np.array([])
flux0 = np.array([])

    
with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150716LPC.hdf5', 'r') as g:
    
    ascc = f['header_table/ascc'].value
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        rc1 = g['data/'+ascc[i]]
        rc2 = g['data2/'+ascc[i]]
        
        jdmid = lc['jdmid']
        flux0 = lc['flux0']
        camtrans0 = rc1['camtrans0']
        intrapix0 = rc1['intrapix0']
        skytrans0 = rc2['skytrans0']
        scflux0 = rc2['scflux0']
        day = np.floor(jdmid)[0]
        jdmid = jdmid-np.floor(jdmid)
        model = camtrans0*intrapix0*skytrans0
        
        plt.figure(figsize=(16,8))
        ax = plt.subplot(311)
        plt.title('ASCC '+ascc[i])
        plt.plot(jdmid, flux0, '.')
        plt.ylabel('flux0')
        plt.subplot(312, sharex=ax)
        plt.plot(jdmid, model, '.')
        plt.ylabel('model')
        plt.subplot(313, sharex=ax)
        plt.plot(jdmid, flux0/model, '.')
        plt.xlabel('Time [JD-%i]'%day)
        plt.ylabel('cflux0')
        plt.tight_layout()
        plt.show()
        plt.close()
        
        
        
        
    

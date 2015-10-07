#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from core import intrapix_weights

filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
filelist = np.sort(filelist)

filelist = filelist[22:]

for filename in filelist:
    ip = intrapix_weights.IntraPixel()
    ip.calculate(filename)

    cf = intrapix_weights.CameraFile(ip.camfile)
    cf.visualize()

#cf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5')

#with h5py.File('/data2/talens/Jul2015/verify_intrapix2.hdf5', 'r') as f:
    #ascc = f['data'].keys()
    
    #for i in range(0, len(ascc), 50):
        #lc = f['data/'+ascc[i]]
        
        #ax = plt.subplot(311)
        #plt.title(ascc[i])
        #plt.plot(lc['mag0'], '.')
        
        #plt.subplot(312, sharex=ax)
        #plt.plot(lc['camtrans0']+lc['intrapix0'], '.')

        #plt.subplot(313, sharex=ax)
        #plt.plot(lc['ipc_mag0'], '.')
        
        #plt.show()
        #plt.close()

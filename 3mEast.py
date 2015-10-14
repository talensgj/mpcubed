#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from core import intrapix_weights
from core import skytrans_weights

from time import time

filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
filelist = np.sort(filelist)
filelist = filelist[10:]

start = time()
for filename in filelist:
    ip = intrapix_weights.IntraPixel()
    ip.calculate(filename)

    cf = intrapix_weights.CameraFile(ip.camfile)
    cf.visualize()
print (time()-start)/60.
    
#filelist = glob.glob('/data2/talens/3mEast/camip_201508??LPE.hdf5')
#filelist = np.sort(filelist)

#for filename in filelist:
    #with h5py.File(filename, 'r') as f:
        #ntot = np.sum(f['header/npoints'].value)
        #chisq = np.sum(f['header/chisq'].value)/len(f['header/chisq'].value)
        #best = np.amax(f['header/chisq'].value)
        #worst = np.amin(f['header/chisq'].value)

        #print filename, ntot, chisq, best, worst

#filelist = glob.glob('/data2/talens/3mEast/fLC_201507??LPE.hdf5')
#filelist = np.sort(filelist)
#filelist = filelist[25:]
#print filelist

#cf = intrapix_weights.CameraFile('/data2/talens/3mEast/camip_20150818LPE.hdf5')
#st = skytrans_weights.SkyTransmission()
#for filename in filelist:
    #cf.correct(filename)
    #st.calculate(filename)
    #sf = skytrans_weights.SkyFile(st.skyfile)
    #sf.visualize()
    #sf.correct()
    
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

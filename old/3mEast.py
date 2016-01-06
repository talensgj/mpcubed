#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from core import intrapix_weights
from core import skytrans_weights

#filelist = glob.glob('/data2/talens/3mEast/fLC_201505??LPE.hdf5')
#filelist = np.sort(filelist)

#ip = intrapix_weights.IntraPixel()
#for filename in filelist:
    #ip.calculate(filename)

    #cf = intrapix_weights.CameraFile(ip.camfile)
    #cf.visualize()

#filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
#filelist = np.sort(filelist)

#cf = intrapix_weights.CameraFile('/data2/talens/3mEast/camip_20150611LPE.hdf5')
#st = skytrans_weights.SkyTransmission()
#for filename in filelist:
    #cf.correct(filename)
    #st.calculate(filename)
    #sf = skytrans_weights.SkyFile(st.skyfile)
    #sf.visualize()
    #sf.correct()
    
#filelist = glob.glob('/data2/talens/3mEast/fLC_201507??LPE.hdf5')
#filelist = np.sort(filelist)
    
#cf = intrapix_weights.CameraFile('/data2/talens/3mEast/camip_20150716LPE.hdf5')
#st = skytrans_weights.SkyTransmission()
#for filename in filelist:
    #cf.correct(filename)
    #st.calculate(filename)
    #sf = skytrans_weights.SkyFile(st.skyfile)
    #sf.visualize()
    #sf.correct()

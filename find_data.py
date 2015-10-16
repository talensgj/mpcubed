#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

filelist = glob.glob('/data2/mascara/LaPalma/2015????LPE/fLC/fLC_2015????LPE.hdf5')
filelist = np.sort(filelist)

for filename in filelist:
    with h5py.File(filename, 'r') as f:
        try: si = f['header/807144'].value
        except: print filename, 0
        else: print filename, si['nobs']

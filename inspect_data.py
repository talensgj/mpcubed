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

import TransitModels

from core import index_functions
from core import BLS_ms
from core import sysrem

P = 2.21857312
Tp = 2454037.612

filelist = glob.glob('/data2/mascara/LaPalma/20150609LPE/fLC/fLC_20150609LPE.hdf5')
filelist = np.sort(filelist)

lst = np.array([])
lstidx = np.array([])
jdmid = np.array([])
flux0 = np.array([])
eflux0 = np.array([])
flags = np.array([])

for filename in filelist:
    print filename
    with h5py.File(filename, 'r') as f:
    
        try:
            lc = f['data/807144'].value
        except:
            continue
    
        try:
            lst = np.append(lst, lc['lst'])
        except: pass
        else:
            lstidx = np.append(lstidx, lc['lstidx'])
            jdmid = np.append(jdmid, lc['jdmid'])
            flux0= np.append(flux0, lc['flux0'])
            eflux0 = np.append(eflux0, lc['eflux0'])
            flags = np.append(flags, lc['flag'])

phase = (jdmid - Tp)/P - 1400
phase = (phase + .5)%1 - .5

time = np.linspace(Tp, Tp+P, 600)
model = TransitModels.MA_Model(time, Tp, P, 0.0995*1.138/0.805, 215.1*0.03142/0.805, 0.0041, 85.51/180*np.pi, np.pi/2., 0., 0.)
mphase = (time - Tp)/P
mphase = (mphase + .5)%1 -.5

sort = np.argsort(mphase)
mphase = mphase[sort]
model = model[sort]

plt.figure(figsize=(16,8))
ax = plt.subplot(211)
plt.plot(phase, flux0, '.')
plt.ylim(0, 20000)
plt.subplot(212, sharex=ax)
plt.plot(mphase, model)
plt.xlim(-.1, .1)
plt.ylim(.95, 1.05)
plt.xlabel('Phase')
plt.ylabel('Flux')
plt.show()

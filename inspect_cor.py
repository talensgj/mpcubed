#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/2015Q2/LPE/cor_201506ALPE.hdf5', 'r') as f:

    lst = f['data/807144/lst'].value
    jdmid = f['data/807144/jdmid'].value
    lstseq = f['data/807144/lstseq'].value

    x = f['data/807144/x'].value
    y = f['data/807144/y'].value
    nobs = f['data/807144/nobs'].value

    plt.plot(jdmid, x, 'o')
    plt.show()

    plt.plot(jdmid, y, 'o')
    plt.show()
    
    plt.plot(jdmid, nobs, 'o')
    plt.show()

    camtrans = f['data/807144/camtrans'].value
    ecamtrans = f['data/807144/ecamtrans'].value
    
    plt.errorbar(jdmid, camtrans, yerr = ecamtrans, fmt='o')
    plt.show()

    skytrans = f['data/807144/skytrans'].value
    eskytrans = f['data/807144/eskytrans'].value
    
    plt.errorbar(jdmid, skytrans, yerr = eskytrans, fmt='o')
    plt.show()

    mag = f['data/807144/mag0'].value
    emag = f['data/807144/emag0'].value
    
    plt.errorbar(jdmid, mag, yerr = emag, fmt='o')
    plt.show()

    sky = f['data/807144/sky'].value
    esky = f['data/807144/esky'].value
    
    plt.errorbar(jdmid, sky, yerr = esky, fmt='o')
    plt.show()

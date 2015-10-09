#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/3mEast/fLC_20150618LPE.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    
    for i in range(0, len(ascc), 50):
        
        lc = f['data/'+ascc[i]]
        
        plt.figure(figsize=(16,8))
        ax = plt.subplot(211)
        plt.title('ASCC ' + ascc[i])
        plt.plot(lc['jdmid'], lc['flux0'], '.')
        plt.plot(lc['jdmid'], lc['flux1'], '.')
        plt.ylabel('Flux')
        plt.subplot(212, sharex=ax)
        plt.plot(lc['jdmid'], lc['sky'], '.')
        plt.xlabel('Time [JD]')
        plt.ylabel('Sky')
        plt.tight_layout()
        plt.show()
        plt.close()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data4/mascara/LaPalma/20150901LPS/fLC/fLC_20150901LPS.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    jdstart = f['header_table/jdstart'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    vmag = f['header_table/vmag'].value
    bmag = f['header_table/bmag'].value
    spectype = f['header_table/spectype'].value
    blend = f['header_table/blend'].value
    blendval = f['header_table/blendvalue'].value
    nobs = f['header_table/nobs'].value
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        print f['header/'+ascc[i]].value
        print jdstart[i], ra[i], dec[i], vmag[i], bmag[i], spectype[i], blend[i], blendval[i], nobs[i]
        ax = plt.subplot(421)
        plt.errorbar(lc['jdmid'], lc['flux0'], yerr=lc['eflux0'], fmt='.')
        plt.errorbar(lc['jdmid'], lc['flux1'], yerr=lc['eflux1'], fmt='.')
        plt.subplot(422, sharex=ax)
        plt.errorbar(lc['jdmid'], lc['sky'], yerr=lc['esky'], fmt='.')
        plt.plot(lc['jdmid'], lc['peak'])
        plt.subplot(423, sharex=ax)
        plt.plot(lc['jdmid'], lc['x'])
        plt.plot(lc['jdmid'], lc['y'])
        plt.subplot(424, sharex=ax)
        plt.plot(lc['jdmid'], lc['az'])
        plt.plot(lc['jdmid'], lc['alt'])
        plt.subplot(425, sharex=ax)
        plt.plot(lc['jdmid'], lc['ccdtemp'])
        plt.subplot(426, sharex=ax)
        plt.plot(lc['jdmid'], lc['exptime'])
        plt.subplot(427, sharex=ax)
        plt.plot(lc['jdmid'], lc['flag'])
        plt.subplot(428, sharex=ax)
        plt.plot(lc['jdmid'], lc['lst'])
        plt.plot(lc['jdmid'], lc['lstidx'])
        plt.plot(lc['jdmid'], lc['lstseq'])
        plt.show()
        plt.close()
        

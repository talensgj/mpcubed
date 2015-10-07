#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from index_functions import index_statistics

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5') as f:
    
    ascc = f['header_table/ascc'].value
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        
        error1 = np.sqrt(.67)*np.sqrt(lc['flux0']+np.pi*(2.5)**2*lc['sky'])
        error2 = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
        
        plt.plot(lc['eflux0'], error1, '.')
        plt.show()

        plt.plot(error1, error2, '.')
        plt.show()

        ax = plt.subplot(411)
        plt.errorbar(lc['jdmid'], lc['flux0'], yerr=lc['eflux0'], fmt='.')
        plt.subplot(412, sharex=ax, sharey=ax)
        plt.errorbar(lc['jdmid'], lc['flux0'], yerr=error1, fmt='.')
        plt.subplot(413, sharex=ax, sharey=ax)
        plt.errorbar(lc['jdmid'], lc['flux0'], yerr=error2, fmt='.')
        plt.subplot(414)
        plt.errorbar(lc['jdmid'], lc['sky'], yerr=lc['esky'], fmt='.')
        plt.show()

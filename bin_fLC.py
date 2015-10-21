#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from core import index_functions

def bin_fLC(lc):
    
    binidx = lc['lstidx'].astype('int')//50

    lstidx = np.unique(binidx)
    count = index_functions.index_statistics(binidx, None, statistic='count')
    flux0 = index_functions.index_statistics(binidx, lc['flux0'], statistic='mean')
    flux1 = index_functions.index_statistics(binidx, lc['flux1'], statistic='mean')
    eflux0 = index_functions.index_statistics(binidx, lc['flux0'], statistic='std')/np.sqrt(count)
    eflux1 = index_functions.index_statistics(binidx, lc['flux1'], statistic='std')/np.sqrt(count)
    jdmid = index_functions.index_statistics(binidx, lc['jdmid'], statistic='mean')
    lst = index_functions.index_statistics(binidx, lc['lst'], statistic='mean')
    
    #plt.plot(lc['lstidx'], lc['flux0'], '.')
    #plt.errorbar(lstidx*50, flux0, yerr=eflux0, fmt='o')
    #plt.show()
    
    arlist = [jdmid, lst, lstidx, flux0, eflux0, flux1, eflux1, count]
    names = ['jdmid', 'lst', 'lstidx', 'flux0', 'eflux0', 'flux1', 'eflux1', 'count']
    record = np.rec.fromarrays(arlist, names=names)
    
    return record
    
with h5py.File('/data2/talens/3mEast/fLC_20150601LPE.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/bfLC_20150601LPE.hdf5') as g:
    
    ascc = f['header_table/ascc'].value
    
    f.copy('header_table', g)
    
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]].value
        record = bin_fLC(lc)
        g.create_dataset('data/'+ascc[i], data=record)
        

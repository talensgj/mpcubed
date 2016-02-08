#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt


fLCfile = '/data2/talens/2015Q2/LPS/fLC_201506ALPS.hdf5'
with h5py.File(fLCfile, 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    
select = (dec < -29.75) & (dec > -30)
ascc = ascc[select]

for i in range(len(ascc)):
    
    with h5py.File(fLCfile, 'r') as f:
        lc = f['data/'+ascc[i]].value
    
    lstday = lc['lstseq']//13500
    lstidx = lc['lstseq']%13500
    
    if np.ptp(lstidx) < 2000: continue
        
    lstday = lstday - np.amin(lstday)
    
    tmp = np.full((np.amax(lstday) + 1, 13500), fill_value=np.nan)
    tmp[lstday, lstidx] = lc['flux0']    
        
    vmin = np.percentile(lc['flux0'], 1)
    vmax = np.percentile(lc['flux0'], 99)
        
    plt.title('ASCC {}'.format(ascc[i]))
    plt.imshow(tmp, aspect='auto', interpolation='None', vmin=vmin, vmax=vmax)
    plt.xlim(np.amin(lstidx)-.5, np.amax(lstidx)+.5)
    plt.colorbar()
    plt.show()
    plt.close()

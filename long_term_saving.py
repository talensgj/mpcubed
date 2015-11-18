#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from fLCfile import fLCfile

import matplotlib.pyplot as plt

from core.index_functions import index_statistics

filelist = glob.glob('/data2/mascara/LaPalma/20150[5-7]??LPE/fLC/fLC_*.hdf5')
filelist = np.sort(filelist)
print filelist
exit()
        
        
#first = True
#for filename in filelist:
    
    #if filename[-8:] == 'old.hdf5': print filename[-8:]
    
    #with h5py.File(filename, 'r') as f:
        
        #try: 
            #tmp = f['data/807144'].value
        #except:
            #pass
        #else:
            #if first:
                #lc = tmp
                #first = False
            #else:
                #lc = stack_arrays((lc, tmp), asrecarray=True)

## Write the lightcurve to file.
#with h5py.File('/data2/talens/ASCC807144_201505_201507.hdf5') as f:
    #f.create_dataset('LPE', data = lc)
    
#print lc.dtype.names
#('flag', 'flux0', 'eflux0', 'flux1', 'eflux1', 'sky', 'esky', 'peak', 'x', 'y', 'alt', 'az', 'ccdtemp', 'exptime', 'jdmid', 'lst', 'lstidx', 'lstseq') 

with h5py.File('/data2/talens/fLC_ASCC807144_201505_201507.hdf5') as f:
    lc = f['LPE'].value

tmag0 = 25 - 2.5*np.log10(lc['flux0'])
emag0 = 2.5/np.log(10.)*np.abs(lc['eflux0']/lc['flux0'])

tmag1 = 25 - 2.5*np.log10(lc['flux1'])
emag1 = 2.5/np.log(10.)*np.abs(lc['eflux1']/lc['flux1'])

binidx = lc['lstseq']//50
binidx = binidx - np.amin(binidx)

nobs = index_statistics(binidx, None, statistic='count')
flag = index_statistics(binidx, lc['flag'], statistic='sum')
flag = np.where(flag > 0, 1, 0)
mag0 = index_statistics(binidx, tmag0, statistic='mean')
emag0 = index_statistics(binidx, tmag0, statistic='std')/np.sqrt(nobs)
mag1 = index_statistics(binidx, tmag1, statistic='mean')
emag1 = index_statistics(binidx, tmag1, statistic='std')/np.sqrt(nobs)
jdmid = index_statistics(binidx, lc['jdmid'], statistic='mean')
lst = index_statistics(binidx, lc['lst'], statistic='mean')

#sky = index_statistics(binidx, lc['sky'], statistic='mean')
#esky = index_statistics(binidx, lc['esky'], statistic='std')
#peak = index_statistics(binidx, lc['peak'], statistic=np.nanmax)
#x = index_statistics(binidx, lc['x'], statistic='mean')
#y = index_statistics(binidx, lc['y'], statistic='mean')
#alt = index_statistics(binidx, lc['alt'], statistic='mean')
#az = index_statistics(binidx, lc['az'], statistic='mean')
#ccdtemp = index_statistics(binidx, lc['ccdtemp'], statistic='mean')
#eccdtemp = index_statistics(binidx, lc['ccdtemp'], statistic='std')
#exptime = index_statistics(binidx, lc['exptime'], statistic='mean')
#eexptime = index_statistics(binidx, lc['exptime'], statistic='std')

lstseq = np.unique(binidx)

arlist = [flag, mag0, emag0, mag1, emag1, nobs, jdmid, lst, lstseq]
names = ['flag', 'mag0', 'emag0', 'mag1', 'emag1', 'nobs', 'jdmid', 'lst', 'lstseq']
record = np.rec.fromarrays(arlist, names=names)

with h5py.File('/data2/talens/bfLC_ASCC807144_201505_201507.hdf5') as f:
    f.create_dataset('LPE', data=record)

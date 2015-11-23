#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

from numpy.lib.recfunctions import stack_arrays

from fLCfile import fLCfile
from core.index_functions import index_statistics

def long_term_example(filelist):
    
    # Create a list of observed stars.
    ascc = np.array([])
    for filename in filelist:
        obj = fLCfile(filename)
        tmp, = obj.read_header(['ascc'])
        ascc = np.append(ascc, tmp)
    ascc = np.unique(ascc)
        
    # Gather the full lightcurve of each star.
    for i in range(len(ascc)):
        print i, len(ascc)
        
        sID = ascc[i]
        
        first = True
        for filename in filelist:
            
            with h5py.File(filename, 'r') as f:
                
                try: 
                    tmp = f['data/' + sID].value
                except:
                    pass
                else:
                    if first:
                        lc = tmp
                        first = False
                    else:
                        lc = stack_arrays((lc, tmp), asrecarray=True)
        
        # Create a binned version.
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

        lstseq = np.unique(lc['lstseq']//50)

        arlist = [flag, mag0, emag0, mag1, emag1, nobs, jdmid, lst, lstseq]
        names = ['flag', 'mag0', 'emag0', 'mag1', 'emag1', 'nobs', 'jdmid', 'lst', 'lstseq']
        
        with h5py.File('/data2/talens/Saving/QuarterLPE.hdf5') as f:
            for name, array in zip(names, arlist):
                f.create_dataset('data/' + sID + '/' + name, data = array)
        
    return

filelist = glob.glob('/data2/mascara/LaPalma/20150[5-7]??LPE/fLC/fLC_*LPE.hdf5')
filelist = np.sort(filelist)
print filelist

long_term_example(filelist)

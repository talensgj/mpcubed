#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

from fLCfile import fLCfile

import matplotlib.pyplot as plt

from numpy.lib.recfunctions import stack_arrays

#filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
#filelist = np.sort(filelist)
#filelist = filelist[:13]

#print filelist

#ascc = np.array([])
#for filename in filelist:
    #obj = fLCfile(filename)
    #tmp, = obj.read_header(['ascc'])
    #ascc = np.append(ascc, tmp)

#ascc = np.unique(ascc)

#ra = np.zeros(len(ascc))
#dec = np.zeros(len(ascc))
#vmag = np.zeros(len(ascc))
#nobs = np.zeros(len(ascc))

#for i in range(len(ascc)):
    #sID = ascc[i]
    #print sID
    #first = True
    #for filename in filelist:
        
        #with h5py.File(filename, 'r') as f:
            #try: 
                #tmp = f['data/'+sID].value
            #except:
                #pass
            #else:
                #if first:
                    #lc = tmp
                    #tmp = f['header/'+sID].value
                    #ra[i] = tmp['ra']
                    #dec[i] = tmp['dec']
                    #vmag[i] = tmp['vmag']
                    #first=False
                #else:
                    #lc = stack_arrays((lc, tmp), asrecarray=True)
    #nobs[i] = lc.size
    #with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5') as f:
        #f.create_dataset('data/' + sID, data = lc)
            
#with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5') as f:
    #f.create_dataset('header_table/ascc', data = ascc)
    #f.create_dataset('header_table/ra', data = ra)
    #f.create_dataset('header_table/dec', data = dec)
    #f.create_dataset('header_table/nobs', data = nobs)
    #f.create_dataset('header_table/vmag', data = vmag)
        
with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5') as f:
    ascc = f['header_table/ascc'].value
    nobs = np.zeros(len(ascc))
    for i in range(len(ascc)):
        nobs[i] = f['data/'+ascc[i]].value.size
    f.create_dataset('header_table/nobs', data = nobs)
        
        
        
    
    
    
    

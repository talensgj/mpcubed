#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py 
import numpy as np

def fix_nobs(filename):

    with h5py.File(filename) as f:
        ascc = f['header_table/ascc'].value
        nobs = f['header_table/nobs'].value
        
        if (len(ascc) != len(f['header'].keys())) | (len(ascc) != len(f['data'].keys())):
            print 'The number of stars do not match.'
            print 'header: %i'%len(f['header'].keys())
            print 'header_table: %i'%len(ascc)
            print 'data: %i'%len(f['data'].keys())
        
        new_nobs = np.zeros(nobs.shape)
        for i in range(len(ascc)):
            
            new_nobs[i] = f['data/' + ascc[i]].size
            f['header/' + ascc[i]]['nobs'] = f['data/' + ascc[i]].size
            
        del f['header_table/nobs']
        f.create_dataset('header_table/nobs', data = new_nobs)
    
    return
    
if __name__ == '__main__':
    
    filelist = glob.glob('/data2/mascara/LaPalma/201506??LPE/fLC/201506??LPE.hdf5')
    
    for filename in filelist:
        fix_nobs(filename)

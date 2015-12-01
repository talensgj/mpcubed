#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py 
import numpy as np

with h5py.File('/data2/talens/3mEast/LBtests/test.hdf5') as f:
    ascc = f['header_table/ascc'].value
    nobs = f['header_table/nobs'].value
    
    new_nobs = np.zeros(nobs.shape)
    for i in range(len(ascc)):
        
        new_nobs[i] = f['data/' + ascc[i]].size
        
    del f['header_table/nobs']
    f.create_dataset('header_table/nobs', data = new_nobs)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py 

with h5py.File('/data2/mascara/LaPalma/20150631LPE/fLC/fLC_20150631LPE.hdf5', 'r') as f:
    ascc = f['header_table/ascc'].value
    nobs = f['header_table/nobs'].value
    
    for i in range(len(ascc)):
        
        if f['data/'+ascc[i]].value.size != nobs[i]:
            print ascc[i]

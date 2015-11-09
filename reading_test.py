#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/mascara/LaPalma/20150508LPE/fLC/fLC_20150508LPE.hdf5', 'r') as f:

    ascc = f['header_table/ascc'].value
    nobs = f['header_table/nobs'].value
    
    for i in range(len(ascc)):
        
        si = f['header/' + ascc[i]].value
        lc = f['data/' + ascc[i]].value
        
        if not np.allclose(lc['lstidx'], lc['lstseq']%13500):
            print 'Bad lstseq'
        
        if not (lc.size == nobs[i]):
            print 'nobs mismatch in table_header'
        
        if not (lc.size == si['nobs']):
            print 'nobs mismatch in header'

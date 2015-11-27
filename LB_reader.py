#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

def read_intermediate(filename):
    
    with h5py.File(filename, 'r') as f:
        ascc = f['header_table/ascc'].value
        nobs = f['header_table/nobs'].value
        
        camtransidx = np.zeros(np.sum(nobs))
        intrapixidx = np.zeros(np.sum(nobs))
        select = np.append(0, np.cumsum(nobs))
        
        for i in range(len(ascc)):
            camtransidx[select[i]:select[i+1]] = f['data/'+ascc[i]+'/camtransidx'].value
            intrapixidx[select[i]:select[i+1]] = f['data/'+ascc[i]+'/intrapixidx'].value
        
if __name__ == '__main__':
    read_intermediate('/data2/talens/3mEast/LBtests/intermediate_file.hdf5')
    

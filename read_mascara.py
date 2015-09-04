#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py

class ReadFile():
    
    def __init__(self, filename, filetype='fLC'):
        
        self.filename = filename
        self.filetype = filetype
        
        return
            
    def header(self, fields):
            
        with h5py.File(self.filename) as f:
            
            result = [f['table_header/'+field].value for field in fields] 

        return result
        
    def data(self, fields, ascc=None):
        
        if ascc is None:
            
            ascc, nobs = self.header(['ascc', 'nobs'])
            
        with h5py.File(self.filename) as f:
        
            result = np.empty((len(fields), np.sum(nobs)))
            select = np.append(0, np.cumsum(nobs))
            for i in range(len(ascc)):
                print ascc[i]
                for j in range(len(fields)):
                
                    result[j, select[i]:select[i+1]] = f['data/'+ascc[i]][fields[j]]
                
        return result
            
            

import matplotlib.pyplot as plt
import numpy as np

rf = ReadFile('/data2/talens/Jul2015/fLC_20150716LPC.hdf5')
ascc, ra, dec, vmag = rf.header(['ascc', 'ra', 'dec', 'vmag'])
result = rf.data(['jdmid', 'flux0'])

plt.plot(result[0], result[1], '.')
plt.show()

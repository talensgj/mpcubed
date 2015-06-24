#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py

class ReadMascaraFile():
    def __init__(filename):
        
        with h5py.File(filename) as f:
        
            ascc = f['header'].keys()
            
            ra = np.zeros(len(ascc))
            dec = np.zeros(len(ascc))
            Nobs = np.zeros(len(ascc))
            for i in range(len(ascc)):
                ra[i] = f['header'][ascc[i]]['ra']
                dec[i] = f['header'][ascc[i]]['dec']
                Nobs[i] = f['header'][ascc[i]]['nobs']

            self.ascc = ascc
            self.ra = ra
            self.dec = dec
            self.nobs = nobs
    
    def read_star(ascc, fields):
        
        with h5py.File(filename) as f:
            
            result = [f['data'][ascc][name] for name in fields]
        
        return result
                
                

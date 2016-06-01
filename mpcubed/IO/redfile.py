#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py

class redFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
        
        return
        
    def read_header(self, fields):
        
        hdr = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['header']
            
            for field in fields:
                hdr[field] = grp[field].value
        
        return hdr
        
    def read_lightcurve(self, ascc, fields):
        
        lc = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/'+ascc]
            
            for field in fields:
                lc[field] = grp[field].value
                
        return lc

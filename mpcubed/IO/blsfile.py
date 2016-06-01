#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py

class blsFile(object):
    
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
        
    def read_data(self, fields):
        
        data = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data']
            
            for field in fields:
                data[field] = grp[field].value
        
        return data

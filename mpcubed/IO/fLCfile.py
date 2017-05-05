#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

class Photfile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
    
        return
        
    def read_global(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            data = f['global'].attrs.items()
        
        return data
        
    def read_stars(self, fields=None):
        
        stars = dict()        
            
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['header_table']

            if fields is None:
                fields = grp.keys()
                
            for field in fields:
            
                if field in grp.keys():
                    stars[field] = grp[field].value
                else:
                    print 'Warning: skipping field {}, field not found.'.format(field)
                    
        return stars
    
    def read_lightcurves(self, ascc=None, fields=None, perstar=True):
        
        onestar = False        
        
        if ascc is None:
            stars = self.read_stars(['ascc'])
            ascc = stars['ascc']
            
        elif isinstance(ascc, basestring):
            onestar = True
            ascc = [ascc]
       
        nstars = len(ascc)
        curves = dict()
        nobs = np.zeros(nstars, dtype='int')
            
        # Read the data.
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data']
            ascc0 = set(grp.keys())            
            
            for i in range(nstars):
                
                if ascc[i] in ascc0:
                    curves[ascc[i]] = grp[ascc[i]].value
                    nobs[i] = len(curves[ascc[i]])
                else:
                    print 'Warning: skipping star {}, star not found.'.format(ascc[i])
                    
        if not curves:
            return curves
            
        # Select specified fields.
        if fields is not None:
            
            for i in range(nstars):
                curves[ascc[i]] = curves[ascc[i]][fields]
                    
        # Combine lightcurves.
        if not perstar:
            
            strides = np.append(0, np.cumsum(nobs))
            tmp = np.recarray(strides[-1], dtype=curves[ascc[0]].dtype)
            
            for i in range(nstars):
                tmp[strides[i]:strides[i+1]] = curves[ascc[i]]
                        
            curves = tmp                        
        
        if onestar:
            return curves[ascc[0]]             
            
        return curves

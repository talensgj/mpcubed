#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from fLCfile import fLCfile


def create_longbaseline(filelist, LBfile):
    """
        Merges a list of fLC files to a single large file.
    """    
    
    # Create a list of observed stars.
    ascc = np.array([])
    for filename in filelist:
        obj = fLCfile(filename)
        tmp, = obj.read_header(['ascc'])
        ascc = np.append(ascc, tmp)
    ascc = np.unique(ascc)
    
    # Gather the full lightcurve of each star.
    ra = np.zeros(len(ascc))
    dec = np.zeros(len(ascc))
    vmag = np.zeros(len(ascc))
    nobs = np.zeros(len(ascc))
    for i in range(len(ascc)):
        sID = ascc[i]
        
        first = True
        for filename in filelist:
            
            with h5py.File(filename, 'r') as f:
                
                try: 
                    tmp = f['data/' + sID].value
                except:
                    pass
                else:
                    if first:
                        lc = tmp
                        tmp = f['header/' + sID].value
                        ra[i] = tmp['ra']
                        dec[i] = tmp['dec']
                        vmag[i] = tmp['vmag']
                        first = False
                    else:
                        lc = stack_arrays((lc, tmp), asrecarray=True)
                        
        nobs[i] = lc.size
        
        # Write the lightcurve to file.
        with h5py.File(LBfile) as f:
            f.create_dataset('data/' + sID, data = lc)
    
    # Write the header to file.
    with h5py.File(LBfile) as f:
        f.create_dataset('header_table/ascc', data = ascc)
        f.create_dataset('header_table/ra', data = ra)
        f.create_dataset('header_table/dec', data = dec)
        f.create_dataset('header_table/nobs', data = nobs)
        f.create_dataset('header_table/vmag', data = vmag)
    
    return
    
        
if __name__ == '__main__':
    
        
    
    
    
    

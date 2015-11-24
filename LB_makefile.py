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
        
        ha = np.mod(lc['lst']*15. - ra[i], 360.)
        tmp = np.repeat(dec[i], nobs[i])
        
        camtransidx = camgrid.find_gridpoint(ha, tmp)
        intrapixidx = ipxgrid.find_gridpoint(ha, tmp)
        
        tmp1 = np.amin(lc['lstseq'])
        tmp2 = np.amax(lc['lstseq'])
        
        try:
            if tmp1 < lstmin: lstmin = tmp1
        except:
            lstmin = tmp1
        
        try:
            if tmp2 > lstmax: lstmax = tmp2
        except:
            lstmax = tmp2
        
        # Write the lightcurve to file.
        with h5py.File(LBfile) as f:
            f.create_dataset('data/' + sID, data = lc)
    
    staridx = np.arange(len(ascc))
    decidx = camgrid.find_decidx(dec)
    skyidx = skygrid.find_gridpoint(ra, dec)
    
    # Write the header to file.
    with h5py.File(LBfile) as f:
        f.create_dataset('header_table/ascc', data = ascc)
        f.create_dataset('header_table/ra', data = ra)
        f.create_dataset('header_table/dec', data = dec)
        f.create_dataset('header_table/nobs', data = nobs)
        f.create_dataset('header_table/vmag', data = vmag)
        f['data'].attrs['lstmin'] = lstmin
        f['data'].attrs['lstmax'] = lstmax
    
    return
    
        
if __name__ == '__main__':
    
        
    filelist = glob.glob('/data2/mascara/LaPalma/201506??LPE/fLC/fLC_*.hdf5')
    filelist = np.sort(filelist)
    
    #filelist = filelist[:5]
    filelist = filelist[15:]
    
    create_longbaseline(filelist, '/data2/talens/3mEast/LBtests/test.hdf5')
    

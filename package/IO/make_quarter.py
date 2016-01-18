#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

   
def _merge_global(filelist):
    
    with h5py.File(filelist[0]) as f:
        data = f['global'].attrs.items()
    
    return data
   
def _merge_headers(filelist):
    
    nfiles = len(filelist)
    
    # Create arrays.
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    vmag = np.array([])
    bmag = np.array([])
    spectype = np.array([])
    
    # Read data.
    for i in range(nfiles):
        with h5py.File(filelist[i], 'r') as f:
            grp = f['header_table']
            ascc = np.append(ascc, grp['ascc'].value)
            ra = np.append(ra, grp['ra'].value)
            dec = np.append(dec, grp['dec'].value)
            vmag = np.append(vmag, grp['vmag'].value)
            bmag = np.append(bmag, grp['bmag'].value)
            spectype = np.append(spectype, grp['spectype'].value)

    # Select unique entries.
    ascc, args = np.unique(ascc, return_index=True)

    # Put result in a dictionary.
    hdr = {}
    hdr['ascc'] = ascc
    hdr['ra'] = ra[args]
    hdr['dec'] = dec[args]
    hdr['vmag'] = vmag[args]
    hdr['bmag'] = bmag[args]
    hdr['spectype'] = spectype[args]
    
    return hdr
    
def _merge_data(filelist, ascc):
    
    nfiles = len(filelist)
    
    first = True
    for i in range(nfiles):
        
        with h5py.File(filelist[i], 'r') as f:
                
            # Try to read the star.
            try:
                tmp = f['data/' + ascc].value
            except:
                continue
            
        # Add the data to the lightcurve.
        if first:
            lc = tmp
            first = False
        else:
            lc = stack_arrays((lc, tmp), asrecarray=True)
                
    return lc

def make_quarter(filelist, outfile):
    """ Merge a list of reduced data files.
    
    Args:
        filelist (str): A list of reduced data files.
        outfile (str): The name of the merged file.
        
    """
    
    filelist = np.sort(filelist)
    
    # Write the global group.
    data = _merge_global(filelist)
    data = dict(data)
    with h5py.File(outfile) as f:
        grp = f.create_group('global')
        grp.attrs['station'] = data['station']
        grp.attrs['camera'] = data['camera']
    
    # Merge the headers.
    hdr = _merge_headers(filelist)
    with h5py.File(outfile) as f:
        grp = f.create_group('header')
        for key, value in hdr.iteritems():
            grp.create_dataset(key, data = value)
        
    # Merge the lightcurves.
    ascc = hdr['ascc']
    nstars = len(ascc)
    for i in range(nstars):
        
        # Read the data.
        lc = _merge_data(filelist, ascc[i])
        
        # Write the data.
        with h5py.File(outfile) as f:
            grp = f.create_group('data/' + ascc[i])
            for key in lc.dtype.names:
                grp.create_dataset(key, data = lc[key])

    return

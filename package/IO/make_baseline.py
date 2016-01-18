#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from fLCfile import fLCfile
from ..statistics import statistics

def verify_filelist(filelist):
    
    nfiles = len(filelist)
    
    args = []
    for i in range(nfiles):
        if os.path.isfile(filelist[i]):
            args.append(i)
        else:
            print filelist[i], 'does not exist.'
    
    filelist = filelist[args]
    
    if len(filelist) == 0:
        print 'No valid files.'
        print 'exiting...'
        exit()
    
    return filelist

def _lstrange(filename):
    
    # Read the lstseq.
    f = fLCfile(filename)
    lstseq = f.read_data(['lstseq'])
    
    # Find the minimum and the maximum.
    lstmin = np.amin(lstseq)
    lstmax = np.amax(lstseq)
    
    return lstmin, lstmax

def _merge_global(filelist):
    
    nfiles = len(filelist)
    
    # Create arrays.
    attrdict = {}
    
    aversion = []
    rversion = []
    cversion = []
    
    alt0 = np.zeros(nfiles)
    az0 = np.zeros(nfiles)
    th0 = np.zeros(nfiles)
    x0 = np.zeros(nfiles)
    y0 = np.zeros(nfiles)
    
    exptime = np.zeros(nfiles)
    ccdtemp = np.zeros(nfiles)
    
    # Read the global groups.
    for i in range(nfiles):
        with h5py.File(filelist[i], 'r') as f:
            grp = f['global']
            
            if (i < 1):
                
                try: attrdict['station'] = grp.attrs['station']
                except: attrdict['station'] = grp.attrs['STATION']
            
                try: attrdict['camera'] = grp.attrs['camera']
                except: attrdict['camera'] = grp.attrs['CAMERA']
            
                try: attrdict['naper'] = grp.attrs['naper']
                except: attrdict['naper'] = grp.attrs['NAPER']
                
                try: attrdict['aper0'] = grp.attrs['aper0']
                except: attrdict['aper0'] = grp.attrs['APER0']
                
                try: attrdict['aper1'] = grp.attrs['aper1']
                except: attrdict['aper1'] = grp.attrs['APER1']
                
                try: attrdict['nstaper'] = grp.attrs['nstaper']
                except: attrdict['nstaper'] = grp.attrs['NSTAPER']
                
                try: attrdict['skyrad0'] = grp.attrs['skyrad0']
                except: attrdict['skyrad0'] = grp.attrs['SKYRAD0']
                
                try: attrdict['skyrad1'] = grp.attrs['skyrad1']
                except: attrdict['skyrad1'] = grp.attrs['SKYRAD1']
            
            try: alt0[i] = grp.attrs['alt0']
            except: alt0[i] = grp.attrs['ALT0']
            
            try: az0[i] = grp.attrs['az0']
            except: az0[i] = grp.attrs['AZ0']
            
            try: th0[i] = grp.attrs['th0']
            except: th0[i] = grp.attrs['TH0']
            
            try: x0[i] = grp.attrs['x0']
            except: x0[i] = grp.attrs['X0']
            
            try: y0[i] = grp.attrs['y0']
            except: y0[i] = grp.attrs['Y0']
            
            try: aversion.append(grp.attrs['aversion'])
            except: aversion.append(grp.attrs['AVERSION'])
            
            try: rversion.append(grp.attrs['rversion'])
            except: rversion.append(grp.attrs['RVERSION'])
            
            try: cversion.append(grp.attrs['cversion'])
            except: cversion.append(grp.attrs['CVERSION'])
            
            try: exptime[i] = grp.attrs['exptime']
            except: exptime[i] = grp.attrs['EXPTIME']
            
            try: ccdtemp[i] = grp.attrs['ccdtemp']
            except: ccdtemp[i] = grp.attrs['CCDTEMP']
    
    # Put result in dictionary
    arrdict = {}
    arrdict['alt0'] = alt0
    arrdict['az0'] = az0
    arrdict['th0'] = th0
    arrdict['x0'] = x0
    arrdict['y0'] = y0
    arrdict['aversion'] = np.array(aversion)
    arrdict['rversion'] = np.array(rversion)
    arrdict['cversion'] = np.array(cversion)
    arrdict['ccdtemp'] = ccdtemp
    arrdict['exptime'] = exptime
    
    return attrdict, arrdict

def _merge_headers(filelist):
    
    nfiles = len(filelist)
    
    # Create arrays.
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    vmag = np.array([])
    bmag = np.array([])
    spectype = np.array([])
    blend = np.array([])
    blendvalue = np.array([])
    jdstart = np.array([])
    nobs = np.array([])
    
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
            blend = np.append(blend, grp['blend'].value)
            blendvalue = np.append(blendvalue, grp['blendvalue'].value)
            jdstart = np.append(jdstart, grp['jdstart'].value)
            nobs = np.append(nobs, grp['nobs'].value)
    
    # Select unique entries.
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    
    # Put result in dictionary.
    hdr = {}
    hdr['ascc'] = ascc
    hdr['ra'] = ra[args]
    hdr['dec'] = dec[args]
    hdr['vmag'] = vmag[args]
    hdr['bmag'] = bmag[args]
    hdr['spectype'] = spectype[args]
    hdr['blend'] = blend[args]
    hdr['blendvalue'] = blendvalue[args]
    hdr['jdstart'] = statistics.idxstats(idx, jdstart, statistic=np.amin)
    hdr['nobs'] = np.bincount(idx, nobs)

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

def make_baseline(filelist, outfile):
    """ Merge a list of fLC files.
    
    Args:
        filelist (str): A list of fLC files.
        outfile (str): The name of the merged file.
        
    """

    filelist = np.sort(filelist)
    filelist = verify_filelist(filelist)
    
    # Write the global group.
    lstmin, _ = _lstrange(filelist[0])
    _, lstmax = _lstrange(filelist[-1])
    attrdict, arrdict = _merge_global(filelist)
    with h5py.File(outfile) as f:
        grp = f.create_group('global')
        grp.attrs['lstmin'] = lstmin
        grp.attrs['lstmax'] = lstmax
        grp.create_dataset('filelist', data = filelist)
        
        for key, value in attrdict.iteritems():
            grp.attrs[key] = value
            
        for key, value in arrdict.iteritems():
            grp.create_dataset(key, data = value)
    
    # Merge the headers.
    hdr = _merge_headers(filelist)
    with h5py.File(outfile) as f:
        grp = f.create_group('header_table')
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
            f.create_dataset('data/' + ascc[i], data = lc)
        
    return
    

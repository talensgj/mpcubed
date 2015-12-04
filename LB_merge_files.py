#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from fLCfile import fLCfile
from core.index_functions import index_statistics

def verify_filelist(filelist):
    
    ndates = len(filelist)
    args = []
    
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            if 'global' in f.keys():
                args.append(i)
            else:
                print 'Removing', filelist[i]
                
    filelist = filelist[args]
    
    if len(filelist) == 0:
        print 'No valid files, exiting...'
        exit()
    
    return filelist

def lstrange(filename):
    
    obj = fLCfile(filename)
    lstseq = obj.read_data(['lstseq'])
    lstmin = np.amin(lstseq)
    lstmax = np.amax(lstseq)
    
    return lstmin, lstmax

def combine_global(filelist):
    
    ndates = len(filelist)
    
    attrdict = {}
    
    aversion = []
    rversion = []
    cversion = []
    
    alt0 = np.zeros(ndates)
    az0 = np.zeros(ndates)
    th0 = np.zeros(ndates)
    x0 = np.zeros(ndates)
    y0 = np.zeros(ndates)
    
    exptime = np.zeros(ndates)
    ccdtemp = np.zeros(ndates)
    
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            
            if (i < 1):
                
                try: attrdict['station'] = f['global'].attrs['station']
                except: attrdict['station'] = f['global'].attrs['STATION']
            
                try: attrdict['camera'] = f['global'].attrs['camera']
                except: attrdict['camera'] = f['global'].attrs['CAMERA']
            
                try: attrdict['naper'] = f['global'].attrs['naper']
                except: attrdict['naper'] = f['global'].attrs['NAPER']
                
                try: attrdict['aper0'] = f['global'].attrs['aper0']
                except: attrdict['aper0'] = f['global'].attrs['APER0']
                
                try: attrdict['aper1'] = f['global'].attrs['aper1']
                except: attrdict['aper1'] = f['global'].attrs['APER1']
                
                try: attrdict['nstaper'] = f['global'].attrs['nstaper']
                except: attrdict['nstaper'] = f['global'].attrs['NSTAPER']
                
                try: attrdict['skyrad0'] = f['global'].attrs['skyrad0']
                except: attrdict['skyrad0'] = f['global'].attrs['SKYRAD0']
                
                try: attrdict['skyrad1'] = f['global'].attrs['skyrad1']
                except: attrdict['skyrad1'] = f['global'].attrs['SKYRAD1']
            
            try: alt0[i] = f['global'].attrs['alt0']
            except: alt0[i] = f['global'].attrs['ALT0']
            
            try: az0[i] = f['global'].attrs['az0']
            except: az0[i] = f['global'].attrs['AZ0']
            
            try: th0[i] = f['global'].attrs['th0']
            except: th0[i] = f['global'].attrs['TH0']
            
            try: x0[i] = f['global'].attrs['x0']
            except: x0[i] = f['global'].attrs['X0']
            
            try: y0[i] = f['global'].attrs['y0']
            except: y0[i] = f['global'].attrs['Y0']
            
            try: aversion.append(f['global'].attrs['aversion'])
            except: aversion.append(f['global'].attrs['AVERSION'])
            
            try: rversion.append(f['global'].attrs['rversion'])
            except: rversion.append(f['global'].attrs['RVERSION'])
            
            try: cversion.append(f['global'].attrs['cversion'])
            except: cversion.append(f['global'].attrs['CVERSION'])
            
            try: exptime[i] = f['global'].attrs['exptime']
            except: exptime[i] = f['global'].attrs['EXPTIME']
            
            try: ccdtemp[i] = f['global'].attrs['ccdtemp']
            except: ccdtemp[i] = f['global'].attrs['CCDTEMP']
        
    arrdict = {}
    arrdict['alt0'] = alt0
    arrdict['az0'] = az0
    arrdict['th0'] = th0
    arrdict['x0'] = x0
    arrdict['y0'] = y0
    arrdict['aversion'] = aversion
    arrdict['rversion'] = rversion
    arrdict['cversion'] = cversion
    arrdict['ccdtemp'] = ccdtemp
    arrdict['exptime'] = exptime
        
    return attrdict, arrdict

def combine_header(filelist):
    
    ndates = len(filelist)
    
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
    
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            ascc = np.append(ascc, f['header_table/ascc'].value)
            ra = np.append(ra, f['header_table/ra'].value)
            dec = np.append(dec, f['header_table/dec'].value)
            vmag = np.append(vmag, f['header_table/vmag'].value)
            bmag = np.append(bmag, f['header_table/bmag'].value)
            spectype = np.append(dec, f['header_table/spectype'].value)
            blend = np.append(blend, f['header_table/blend'].value)
            blendvalue = np.append(blendvalue, f['header_table/blendvalue'].value)
            jdstart = np.append(jdstart, f['header_table/jdstart'].value)
            nobs = np.append(nobs, f['header_table/nobs'].value)
    
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    
    headerdict = {}
    
    headerdict['ascc'] = ascc
    headerdict['ra'] = ra[args]
    headerdict['dec'] = dec[args]
    headerdict['vmag'] = vmag[args]
    headerdict['bmag'] = bmag[args]
    headerdict['spectype'] = spectype[args]
    headerdict['blend'] = blend[args]
    headerdict['blendvalue'] = blendvalue[args]
    headerdict['jdstart'] = index_statistics(idx, jdstart, statistic=np.amin)
    headerdict['nobs'] = np.bincount(idx, nobs)
    
    return headerdict

def combine_data(filelist, ascc):
    
    ndates = len(filelist)
    nstars = len(ascc)
    
    stardict = {}
    first = np.ones(nstars, dtype='bool')
    for i in range(ndates):
        
        with h5py.File(filelist[i], 'r') as f:
            
            for j in range(nstars):
                
                try: 
                    tmp = f['data/' + ascc[j]].value
                except: 
                    pass
                else:
                    if first[j]:
                        stardict[ascc[j]] = tmp
                        first[j] = False
                    else:
                        stardict[ascc[j]] = stack_arrays((stardict[ascc[j]], tmp), asrecarray=True)

    return stardict
   
def LBfile(filelist, outfile):
        
    filelist = np.sort(filelist)

    filelist = verify_filelist(filelist)

    # Find the lstseq range.
    lstmin, _ = lstrange(filelist[0])
    _, lstmax = lstrange(filelist[-1])

    # Make the combined global group.
    attrdict, arrdict = combine_global(filelist)
    with h5py.File(outfile) as f:
        
        for key, value in arrdict.iteritems():
            f.create_dataset('global/' + key, data = value)
            
        f.create_dataset('global/filelist', data = filelist)
            
        for key, value in attrdict.iteritems():
            f['global'].attrs[key] = value
            
        f['global'].attrs['lstmin'] = lstmin
        f['global'].attrs['lstmax'] = lstmax
        
    # Make the combined header_table group.
    headerdict = combine_header(filelist)
    with h5py.File(outfile) as f:
        for key, value in headerdict.iteritems():
            f.create_dataset('header_table/' + key, data = value)
            
    ascc = headerdict['ascc']
    
    # Make the combined data group.
    nstars = len(ascc)
    for i in range(0, nstars, 50):
        
        stardict = combine_data(filelist, ascc[i:i+50])
    
        with h5py.File(outfile) as f:
            for key, value in stardict.iteritems():
                f.create_dataset('data/' + key, data = value)

    return
    
if __name__ == '__main__':
    
    filelist = glob.glob('/data2/mascara/LaPalma/201506??LPE/fLC/fLC_201506??LPE.hdf5')
    filelist = np.sort(filelist)
    #filelist = filelist[15:]
    print filelist
    
    LBfile(filelist, '/data2/talens/LongBaselines/fLC_201506LPE.hdf5')

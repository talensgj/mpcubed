# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:56:55 2017

@author: talens
"""

import os

import h5py
import numpy as np

import glob

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
        return filelist
    
    return filelist

def _index_files(filelist):
    
    nfiles = len(filelist)
    
    idx1 = np.array([], dtype='uint16')
    
    stars = dict()
    for i in range(nfiles):
        
        filename = filelist[i]        
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header_table']
            
            for key in grp.keys():
                
                ascc_ = grp['ascc'].value
                
                if key not in stars.keys():
                    stars[key] = grp[key].value
                else:
                    stars[key] = np.append(stars[key], grp[key].value)

            grp = f['data']

            dtype = grp[ascc_[0]].dtype

        idx1 = np.append(idx1, np.repeat(i, len(ascc_)))
    
    ascc, args, idx2 = np.unique(stars['ascc'], return_index=True, return_inverse=True)
    nstars = len(ascc)    
    nobs = np.zeros((nfiles, nstars), dtype='uint32')
    stars['nobs'] = stars['nobs'].astype('uint32')
    nobs[idx1, idx2] = stars['nobs']
    
    for key in stars.keys():
        stars[key] = stars[key][args]
    
    return stars, nobs, dtype

def _read_curves(filelist, ascc, nobs, dtype):
    
    nfiles = len(filelist)
    nstars = len(ascc)
    
    strides = np.row_stack([nstars*[0], np.cumsum(nobs, axis=0)]).astype('int')
    curves = {ascc[i]:np.recarray(strides[-1,i], dtype=dtype) for i in range(nstars)}
    
    for i in range(nfiles):
        
        filename = filelist[i]
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['data']
            
            for j in range(nstars):
                
                if (nobs[i,j] > 0):
                    
                    curves[ascc[j]][strides[i,j]:strides[i+1,j]] = grp[ascc[j]].value
                    
    return curves
    
def _read_global(filelist):
    
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

def combine_photometry(filename, filelist, nsteps=1000):
    
    filelist = np.sort(filelist)    
    filelist = verify_filelist(filelist)
    
    if len(filelist) == 0:
        return
        
    # Read the combined stars field and index the files.
    stars, nobs, dtype = _index_files(filelist)    
    
    stars['lstsqmin'] = np.zeros(len(stars['ascc']), dtype='uint32')
    stars['lstsqmax'] = np.zeros(len(stars['ascc']), dtype='uint32')
    
    nstars = len(stars['ascc'])
    for i in range(0, nstars, nsteps):
        
        # Read the combined lightcurves for a group of stars.
        curves = _read_curves(filelist, stars['ascc'][i:i+nsteps], nobs[:,i:i+nsteps], dtype)
             
        # Write the combined lightcurves for a group of stars.
        with h5py.File(filename) as f:
            
            for j in range(i, i+len(stars['ascc'][i:i+nsteps])):
             
                tmp = curves[stars['ascc'][j]]
                
                stars['nobs'][j] = len(tmp)
                stars['lstsqmin'][j] = tmp['lstseq'][0]
                stars['lstsqmax'][j] = tmp['lstseq'][-1]
                
                f.create_dataset('data/{}'.format(stars['ascc'][j]), data=tmp)    

    # Write the combined "header_table" field.
    with h5py.File(filename) as f:
        
        grp = f.create_group('header_table')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key])
            
    # Write the global group.
    lstmin = np.amin(stars['lstsqmin'])
    lstmax = np.amax(stars['lstsqmax'])
    
    attrdict, arrdict = _read_global(filelist)
    
    with h5py.File(filename) as f:
        
        grp = f.create_group('global')
        
        grp.attrs['lstmin'] = lstmin
        grp.attrs['lstmax'] = lstmax
        
        grp.create_dataset('filelist', data=filelist)
        
        for key, value in attrdict.iteritems():
            grp.attrs[key] = value
            
        for key, value in arrdict.iteritems():
            grp.create_dataset(key, data = value)

    return

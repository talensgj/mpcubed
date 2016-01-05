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
    """
    Given a list of fLC file it checks that they contain data.
    """
    
    ndates = len(filelist)
    args = []
    
    # Check if the global group exists, print a message if it does not.
    for i in range(ndates):
        
        if not os.path.isfile(filelist[i]):
            print filelist[i], 'does not exist.'
            continue
        
        with h5py.File(filelist[i], 'r') as f:
            if 'global' not in f.keys():
                print filelist[i], 'contains no global group.'
            else:
                args.append(i)
    
    # Remove the bad files.
    filelist = filelist[args]
    
    # If the list is now empty exit.
    if len(filelist) == 0:
        print 'No valid files.'
        print 'exiting...'
        exit()
    
    return filelist

def lstrange(filename):
    """ 
    Given an fLC file it reads all lstseq indices and return the minimum
    and maximum indices.
    """
    
    # Read the lstseq.
    f = fLCfile(filename)
    lstseq = f.read_data(['lstseq'])
    
    # Find the minimum and the maximum.
    lstmin = np.amin(lstseq)
    lstmax = np.amax(lstseq)
    
    return lstmin, lstmax

def combine_global(filelist):
    """
    Given a list of fLC files it reads the global and combines them into a
    single global.
    """
    
    ndates = len(filelist)
    
    # Create arrays.
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
    
    # Read the global groups.
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
    
    # Create dictionary.
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
    """
    Given a list of fLC files it reads the header and combines them into a
    single header.
    """
    
    ndates = len(filelist)
    
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
    
    # Read the headers.
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
    
    # Find the unique stars.
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    
    # Create the new header.
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
    """
    Given a list of fLC files and ascc numbers combines all the lightcurves
    for these stars.
    """
    
    ndates = len(filelist)
    nstars = len(ascc)
    
    # Create a dictionary to hold the lightcurves.
    stardict = {}
    first = np.ones(nstars, dtype='bool')
    
    # Read the data.
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            for j in range(nstars):
                
                # Try to read the star.
                try: 
                    tmp = f['data/' + ascc[j]].value
                except: 
                    continue
                
                # Add the data to the dictionary.
                if first[j]:
                    stardict[ascc[j]] = tmp
                    first[j] = False
                else:
                    stardict[ascc[j]] = stack_arrays((stardict[ascc[j]], tmp), asrecarray=True)

    return stardict
   
def merge_files(filelist, outfile):
    """ 
    Given a list of fLC files combines them to one file.
    """
        
    filelist = np.sort(filelist)

    # Check that the files exist and contain data.
    print 'Verifying the filelist.'
    filelist = verify_filelist(filelist)
    print 'Done.'

    if os.path.isfile(outfile):
        print 'Output file already exists:', outfile
        print 'exiting...'
        exit()

    # Find the lstseq range.
    lstmin, _ = lstrange(filelist[0])
    _, lstmax = lstrange(filelist[-1])

    # Make the combined global group.
    print 'Merging the global groups.'
    attrdict, arrdict = combine_global(filelist)
    with h5py.File(outfile) as f:
        
        f.create_dataset('global/filelist', data = filelist)
        f['global'].attrs['lstmin'] = lstmin
        f['global'].attrs['lstmax'] = lstmax
        
        for key, value in arrdict.iteritems():
            f.create_dataset('global/' + key, data = value)
            
        for key, value in attrdict.iteritems():
            f['global'].attrs[key] = value
    print 'Done.'
        
    # Make the combined header_table group.
    print 'Merging the header_table groups.'
    headerdict = combine_header(filelist)
    with h5py.File(outfile) as f:
        
        for key, value in headerdict.iteritems():
            f.create_dataset('header_table/' + key, data = value)
    print 'Done.'
            
    # Make the combined data group.
    print 'Merging the lightcurves.'
    ascc = headerdict['ascc']
    nstars = len(ascc)
    for i in range(0, nstars, 50):
        
        stardict = combine_data(filelist, ascc[i:i+50])
    
        with h5py.File(outfile) as f:
            
            for key, value in stardict.iteritems():
                f.create_dataset('data/' + key, data = value)
    print 'Done'

    return
    
if __name__ == '__main__':
    
    import argparse
    
    filepath = '/data2/mascara/LaPalma'
    outpath = '/data2/talens/LongBaselines'
    
    parser = argparse.ArgumentParser(description = 'Combine fLC files into one file.')
    parser.add_argument('date', help = 'YYYYMM, the year and month to combine files from.')
    parser.add_argument('camera', help = 'The camera to combine files from.')
    parser.add_argument('mode', type = int, choices = [0, 1, 2], help = 'The baseline to create.')
    parser.add_argument('--filepath', help = 'Path to the fLC files.', default = filepath)
    parser.add_argument('--outpath', help = 'Directory to write the combined file to.', default = outpath)
    args = parser.parse_args()
    
    if (args.mode == 0):
        part = 'A'
        dates = [args.date + '%.2i'%i + args.camera for i in range(1, 16)]
    if (args.mode == 1):
        part = 'B'
        dates = [args.date + '%.2i'%i + args.camera for i in range(16, 32)]
    if (args.mode == 2):
        part = ''
        dates = [args.date + '%.2i'%i + args.camera for i in range(1, 32)]
    
    filelist = [os.path.join(args.filepath, '%s/fLC/fLC_%s.hdf5'%(date, date)) for date in dates]
    outfile = os.path.join(args.outpath, 'fLC_%s%s%s.hdf5'%(args.date, part, args.camera))
    
    merge_files(filelist, outfile)

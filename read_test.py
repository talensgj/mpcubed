#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np
import h5py

def add_header_table(filename):
    
    print 'Processing file:', filename
    
    with h5py.File(filename, 'r+') as f:
        
        print f.keys()
        
        if 'header_table' in f.keys():
            print 'header_table exists.'
        
        elif 'table_header' in f.keys():
            print 'table_header exists, renaming...'
            f['header_table'] = f['table_header']
            del f['table_header']
    
        else:
            print 'Creating header_table.'
            
            ascc = f['header'].keys()
            ascc = np.array(ascc).astype('S7')
            
            jdstart = np.zeros(len(ascc))
            ra = np.zeros(len(ascc))
            dec = np.zeros(len(ascc))
            bmag = np.zeros(len(ascc))
            vmag = np.zeros(len(ascc))
            spectype = np.zeros(len(ascc), dtype='S9')
            blend = np.zeros(len(ascc), dtype='S1')
            blendvalue = np.zeros(len(ascc))
            nobs = np.zeros(len(ascc))
        
            for i in range(len(ascc)):
                
                starinfo = f['header/'+ascc[i]]
                
                try: jdstart[i] = starinfo['jdstart']
                except: jdstart[i] = starinfo['jdsart']
                ra[i] = starinfo['ra']
                dec[i] = starinfo['dec']
                bmag[i] = starinfo['bmag']
                vmag[i] = starinfo['vmag']
                spectype[i] = np.squeeze(starinfo['spectype'])
                blend[i] = np.squeeze(starinfo['blend'])
                blendvalue[i] = starinfo['blendvalue']
                nobs[i] = starinfo['nobs']
                
            grp = f.create_group('header_table')
            grp.create_dataset('ascc', data=ascc)
            grp.create_dataset('jdstart', data=jdstart)
            grp.create_dataset('ra', data=ra)
            grp.create_dataset('dec', data=dec)
            grp.create_dataset('bmag', data=bmag)
            grp.create_dataset('vmag', data=vmag)
            grp.create_dataset('spectype', data=spectype)
            grp.create_dataset('blend', data=blend)
            grp.create_dataset('blendvalue', data=blendvalue)
            grp.create_dataset('nobs', data=nobs)
        
        print
        
    return

def fix_lc(filename):
    
    print 'Processing file:', filename
    
    with h5py.File(filename, 'r+') as f:
        
        print f.keys()
        
        ascc = f['header_table/ascc'].value
    
        for i in range(len(ascc)):
            names = f['data/'+ascc[i-1]].value.dtype
            if names is not None: break
    
        print names
    
        for i in range(len(ascc)):
            lc = f['data/'+ascc[i]]
    
            try:
                f['data/'+ascc[i]]['flux0']
            except:
                print 'Encountered bad lightcurve, making recarray...'
                
                data = f['data/'+ascc[i]].value
                lc = np.rec.fromarrays([np.array(item) for item in data], dtype=names)
                
                del f['data/'+ascc[i]]
                f.create_dataset('data/'+ascc[i], data=lc)
    
        print
    
    return

def fix_nobs(filename):
    
    print 'Processing file:', filename
    
    with h5py.File(filename, 'r+') as f:
        
        print f.keys()
        
        ascc = f['header_table/ascc'].value
        
        nobs = np.zeros(len(ascc))
        for i in range(len(ascc)):
            lc = f['data/'+ascc[i]]
            
            if np.isscalar(lc['flux0']): nobs[i] = 1
            else: nobs[i] = len(lc['flux0'])
        
        try: del f['header_table/nobs']
        except : pass
        f.create_dataset('header_table/nobs', data=nobs)
    
    print
    
    return

filelist = glob.glob('/data2/talens/Jul2015/fLC_*')
filelist = np.sort(filelist)
#filelist = filelist[137:]
filelist = ['/data2/talens/Jul2015/'+obj for obj in ['fLC_20150708LPC.hdf5', 'fLC_20150709LPC.hdf5', 'fLC_20150728LPC.hdf5', 'fLC_20150729LPC.hdf5']]
for filename in filelist:
    add_header_table(filename)
    fix_lc(filename)
    fix_nobs(filename)

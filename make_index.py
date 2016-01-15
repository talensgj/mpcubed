#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from package import IO

def make_index(year, camera):
    
    filelist = glob.glob('/data2/mascara/LaPalma/'+year+'????'+camera+'/fLC/fLC_*.hdf5')
    filelist = np.sort(filelist)
    
    rversion = []
    cversion = []
    aversion = []
    alt0 = []
    az0 = []
    th0 = []
    x0 = []
    y0 = []
    
    for i in range(len(filelist)):
        if 'old' in filelist[i]:
            print filelist[i]
            
        f = IO.fLCfile(filelist[i])
        try:
            data = f.read_global()
        except:
            print 'no global', filelist[i]
            rversion.append('-')
            cversion.append('-')
            aversion.append('-')
            alt0.append(alt0[-1])
            az0.append(az0[-1])
            th0.append(th0[-1])
            x0.append(x0[-1])
            y0.append(y0[-1])
        else:
            data = dict(data)
            if len(data.keys()) == 6:
                print filelist[i]
                rversion.append('-')
                cversion.append('-')
                aversion.append('-')
                alt0.append(alt0[-1])
                az0.append(az0[-1])
                th0.append(th0[-1])
                x0.append(x0[-1])
                y0.append(y0[-1])
            else:
                try: rversion.append(data['rversion'])
                except: rversion.append(data['RVERSION'])
                try: cversion.append(data['cversion'])
                except: cversion.append(data['CVERSION'])
                try: aversion.append(data['aversion'])
                except: aversion.append(data['AVERSION'])
                
                try: alt0.append(data['alt0'])
                except: alt0.append(data['ALT0'])
                try: az0.append(data['az0'])
                except: az0.append(data['AZ0'])
                try: th0.append(data['th0'])
                except: th0.append(data['TH0'])
                try: x0.append(data['x0'])
                except: x0.append(data['X0'])
                try: y0.append(data['y0'])
                except: y0.append(data['Y0'])
            
    with h5py.File('/data2/talens/'+year+camera+'.hdf5') as f:
        f.create_dataset('filename', data = filelist)
        
        f.create_dataset('aversion', data = aversion)
        f.create_dataset('rversion', data = rversion)
        f.create_dataset('cversion', data = cversion)
        
        f.create_dataset('alt0', data = alt0)
        f.create_dataset('az0', data = az0)
        f.create_dataset('th0', data = th0)
        f.create_dataset('x0', data = x0)
        f.create_dataset('y0', data = y0)
        
    return filelist

make_index('2014', 'LPN')
make_index('2014', 'LPE')
make_index('2014', 'LPS')
make_index('2014', 'LPW')
make_index('2014', 'LPC')

make_index('2015', 'LPN')
make_index('2015', 'LPE')
make_index('2015', 'LPS')
make_index('2015', 'LPW')
make_index('2015', 'LPC')
    

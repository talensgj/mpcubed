#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from core.index_functions import index_statistics

import matplotlib.pyplot as plt



def pointing(filelist):
    
    ndates = len(filelist)
    
    Alt0 = np.zeros(ndates)
    Az0 = np.zeros(ndates)
    Theta0 = np.zeros(ndates)
    X0 = np.zeros(ndates)
    Y0 = np.zeros(ndates)
    
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            print f['global'].attrs.items()
            
            try: Alt0[i] = f['global'].attrs['ALT0']
            except: Alt0[i] = f['global'].attrs['alt0']
                
            try: Az0[i] = f['global'].attrs['AZ0']
            except: Az0[i] = f['global'].attrs['az0']
            
            try: Theta0[i] = f['global'].attrs['TH0']
            except: Theta0[i] = f['global'].attrs['th0']
            
            try: X0[i] = f['global'].attrs['X0']
            except: X0[i] = f['global'].attrs['x0']
            
            try: Y0[i] = f['global'].attrs['Y0']
            except: Y0[i] = f['global'].attrs['y0']
                
                
    plt.subplot(511)
    plt.plot((Az0-Az0[0])*60, '.')
    plt.ylim(-2,2)
    plt.subplot(512)
    plt.plot((Alt0-Alt0[0])*60, '.')
    plt.ylim(-2,2)
    plt.subplot(513)
    plt.plot((Theta0-Theta0[0])*60, '.')
    plt.ylim(-2,2)
    plt.subplot(514)
    plt.plot(X0-X0[0], '.')
    plt.ylim(-2,2)
    plt.subplot(515)
    plt.plot(Y0-Y0[0], '.')
    plt.ylim(-2,2)
    plt.show()
    plt.close()
    
    return
    
def coverage(filelist):
    
    ndates = len(filelist)
    
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    nobs = np.array([])
    for i in range(ndates):
        with h5py.File(filelist[i]) as f:
            ascc = np.append(ascc, f['header_table/ascc'].value)
            ra = np.append(ra, f['header_table/ra'].value)
            dec = np.append(dec, f['header_table/dec'].value)
            nobs = np.append(nobs, f['header_table/nobs'].value)
    
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    ra = ra[args]
    dec = dec[args]
    nobs = np.bincount(idx, nobs)
    
    nstars = len(ascc)
    ntot = np.sum(nobs)
    nmin = np.amin(nobs)
    nmax = np.amax(nobs)
    nmean = np.mean(nobs)
    
    print 'The baseline contains %i datapoints in %i stars.'%(ntot, nstars)
    print 'nmin = %i, nmax = %i, nmean = %i'%(nmin, nmax, nmean)
    
    plt.scatter(ra, dec, c=nobs)
    plt.xlabel('RA [deg]')
    plt.ylabel('Dec [deg]')
    plt.colorbar()
    plt.xlim(0, 360)
    plt.ylim(-90, 90)
    plt.show()
    plt.close()
    
def combined_header(filelist):
    
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
    ra = ra[args]
    dec = dec[args]
    vmag = vmag[args]
    bmag = bmag[args]
    spectype = spectype[args]
    blend = blend[args]
    blendvalue = blendvalue[args]
    jdstart = index_statistics(idx, jdstart, statistic=np.amin)
    nobs = np.bincount(idx, nobs)
    
    return ascc
    
def combined_lightcurve(filelist, ascc):
    
    ndates = len(filelist)
    nstars = len(ascc)
    
    first = np.ones(nstars, dtype='bool')
    new_dict = {}
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            for j in range(nstars):
            
                try:
                    tmp = f['data/' + ascc[j]].value
                except:
                    pass
                else:
                    if first[j]:
                        new_dict[ascc[j]] = tmp
                        first[j] = False
                    else:
                        new_dict[ascc[j]] = stack_arrays((new_dict[ascc[j]], tmp), asrecarray=True)
                    
    return new_dict
    
def file_generator(filelist):
    
    ascc = combined_header(filelist)
    
    nstars = len(ascc)
    chunk = 50
    for i in range(0, nstars, chunk):
        print i, nstars
        combined_lightcurve(filelist, ascc[i:i+chunk])
    
    return
    
def filelist_generator(fromdate, todate, camera):
    
    return
    
if __name__ == '__main__':
    
    filelist = glob.glob('/data2/mascara/LaPalma/201506??LPC/fLC/fLC_201506??LPC.hdf5')
    filelist = np.sort(filelist)
    filelist = filelist[16:]
    
    file_generator(filelist)

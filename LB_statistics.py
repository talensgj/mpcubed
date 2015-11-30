#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from core.index_functions import index_statistics

import matplotlib.pyplot as plt

from fLCfile import fLCfile

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
    
def combined_global(filelist):
    
    ndates = len(filelist)
    
    alt0 = np.zeros(ndates)
    az0 = np.zeros(ndates)
    th0 = np.zeros(ndates)
    x0 = np.zeros(ndates)
    y0 = np.zeros(ndates)
    
    aversion = []
    rversion = []
    cversion = []
    
    exptime = np.zeros(ndates)
    ccdtemp = np.zeros(ndates)
    
    for i in range(ndates):
        with h5py.File(filelist[i], 'r') as f:
            
            if (i < 1):
                
                try: station = f['global'].attrs['station']
                except: station = f['global'].attrs['STATION']
            
                try: camera = f['global'].attrs['camera']
                except: camera = f['global'].attrs['CAMERA']
            
                try: naper = f['global'].attrs['naper']
                except: naper = f['global'].attrs['NAPER']
                
                try: aper0 = f['global'].attrs['aper0']
                except: aper0 = f['global'].attrs['APER0']
                
                try: aper1 = f['global'].attrs['aper1']
                except: aper1 = f['global'].attrs['APER1']
                
                try: nstaper = f['global'].attrs['nstaper']
                except: nstaper = f['global'].attrs['NSTAPER']
                
                try: skyrad0 = f['global'].attrs['skyrad0']
                except: skyrad0 = f['global'].attrs['SKYRAD0']
                
                try: skyrad1 = f['global'].attrs['skyrad1']
                except: skyrad1 = f['global'].attrs['SKYRAD1']
            
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
            

    return station, camera, alt0, az0, th0, x0, y0, aversion, rversion, cversion, exptime, ccdtemp, naper, aper0, aper1, nstaper, skyrad0, skyrad1
   
def lstrange(filelist):
    
    obj = fLCfile(filelist[0])
    lstseq = obj.read_data(['lstseq'])
    lstmin = np.amin(lstseq)
    
    obj = fLCfile(filelist[-1])
    lstseq = obj.read_data(['lstseq'])
    lstmax = np.amax(lstseq)
    
    return lstmin, lstmax
    
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
    
    return ascc, ra, dec, vmag, bmag, spectype, blend, blendvalue, jdstart, nobs
    
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
    
    station, camera, alt0, az0, th0, x0, y0, aversion, rversion, cversion, exptime, ccdtemp, naper, aper0, aper1, nstaper, skyrad0, skyrad1 = combined_global(filelist)
    
    lstmin, lstmax = lstrange(filelist)
    
    with h5py.File('/data2/talens/3mEast/LBtests/test.hdf5') as f:
        grp = f.create_group('global')
        
        grp.attrs['station'] = station
        grp.attrs['camera'] = camera
        
        grp.attrs['naper'] = naper
        grp.attrs['aper0'] = aper0
        grp.attrs['aper1'] = aper1
        grp.attrs['nstaper'] = nstaper
        grp.attrs['skyrad0'] = skyrad0
        grp.attrs['skyrad1'] = skyrad1
        
        grp.create_dataset('filelist', data=filelist)
        
        grp.create_dataset('alt0', data=alt0)
        grp.create_dataset('az0', data=az0)
        grp.create_dataset('th0', data=th0)
        grp.create_dataset('x0', data=x0)
        grp.create_dataset('y0', data=y0)
        
        grp.create_dataset('rversion', data=rversion)
        grp.create_dataset('cversion', data=cversion)
        grp.create_dataset('aversion', data=aversion)
        
        grp.create_dataset('ccdtemp', data=ccdtemp)
        grp.create_dataset('exptime', data=exptime)
        
        grp.attrs['lstmin'] = lstmin
        grp.attrs['lstmax'] = lstmax
    exit()
    ascc, ra, dec, vmag, bmag, spectype, blend, blendvalue, jdstart, nobs = combined_header(filelist)
    
    with h5py.File('/data2/talens/test.hdf5') as f:
        f.create_dataset('header_table/ascc', data=ascc) 
        f.create_dataset('header_table/ra', data=ra) 
        f.create_dataset('header_table/dec', data=dec) 
        f.create_dataset('header_table/vmag', data=vmag) 
        f.create_dataset('header_table/bmag', data=bmag) 
        f.create_dataset('header_table/spectype', data=spectype) 
        f.create_dataset('header_table/blend', data=blend) 
        f.create_dataset('header_table/blendvalue', data=blendvalue) 
        f.create_dataset('header_table/jdstart', data=jdstart)
        f.create_dataset('header_table/nobs', data=nobs) 
    
    nstars = len(ascc)
    chunk = 50
    for i in range(0, nstars, chunk):
        print i, nstars
        new_dict = combined_lightcurve(filelist, ascc[i:i+chunk])
    
        with h5py.File('/data2/talens/test.hdf5') as f:
            for key, value in new_dict.iteritems():
                f.create_dataset('data/' + key, data=value)
    
    return
    
def filelist_generator(fromdate, todate, camera):
    
    return
    
if __name__ == '__main__':
    
    filelist = glob.glob('/data2/mascara/LaPalma/201506??LPE/fLC/fLC_201506??LPE.hdf5')
    filelist = np.sort(filelist)
    filelist = filelist[14:]
    print filelist
    
    file_generator(filelist)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from reduction_functions import make_array, make_transmission_map

import matplotlib.pyplot as plt
from time import time

import os
import glob

def make_idx(variable, range, nbins):
    
    bins = np.linspace(range[0], range[1], nbins+1)
    idx = np.searchsorted(bins, variable)
    
    if np.any(idx == 0) | np.any(idx == len(bins)):
        print 'Warning: there where values out of range.'
        exit()
        
    idx -= 1
    idx = idx.astype('int')
    
    offset = np.amin(idx)
    length = np.ptp(idx) + 1
    
    return idx, length, offset
    
def transmission(filename):

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lst_idx = np.zeros(np.sum(Nobs)).astype('int')
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst_idx[select[i]:select[i+1]] = data[ascc[i]]['lstidx']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux1']
            eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux1']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']

    # Discretize the ra and dec.
    ra_idx, length, offset = make_idx(ra, (0,360), 13500) # hardcoded...
    dec_idx, length, offset_dec = make_idx(dec, (-90,90), 1440) # hardcoded...

    # Create a unique ha idx.
    ha_idx = lst_idx - np.repeat(ra_idx, Nobs)
    ha_idx = np.mod(ha_idx, 13500) # hardcoded...

    offset_ha = np.amin(ha_idx)
    length = np.ptp(ha_idx)

    # Set bad data to NaN.
    here, = np.where((flags > 0))
    flux0[here] = np.nan
    eflux0[here] = np.nan

    # Cast the data and the errors to arrays suitable for sysrem.
    data = make_array(flux0, Nobs, ha_idx-offset_ha) # still neglecting periodicity...
    error = make_array(eflux0, Nobs, ha_idx-offset_ha)

    # Obtain the transmission map.
    transmission, flags, chi2, npoints, npars, niter = make_transmission_map(data, error, dec_idx-offset_dec)

    # NEED SOME QUALITY CHECKS.
    
    filename = 'T_' + os.path.basename(filename).split('_')[-1]
    print 'Result:', filename
    with h5py.File(filename) as f:
        aper = 'aper1'
        grp = f.create_group(aper)
        
        dset = grp.create_dataset('Transmission', data=transmission)
        dset.attrs['offset_HA'] = offset_ha
        dset.attrs['offset_Dec'] = offset_dec
        dset = grp.create_dataset('Flags', data=flags)
        dset.attrs['offset_HA'] = offset_ha
        dset.attrs['offset_Dec'] = offset_dec
        dset = grp.create_dataset('chi_sq', data=chi2)
        dset.attrs['eps'] = 1e-3 # hardcoded...
        dset = grp.create_dataset('niter', data=niter)
        dset.attrs['maxiter'] = 50 # hardcoded...
        grp.create_dataset('npoints', data=npoints)
        grp.create_dataset('npars', data=npars)

    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/talens/Apr2015LPE/fLC_*')
    
    for filename in filelist:
        print 'Data:', filename
        transmission(filename)

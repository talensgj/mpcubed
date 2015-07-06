#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import os
import glob

import matplotlib.pyplot as plt

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
def apply_transmission(filename, reference):
    
    # Read transmission map.
    with h5py.File(reference) as f:
        dset = f['aper0/Transmission']
        trans = dset.value
        offset_ha = dset.attrs['offset_HA']
        offset_dec = dset.attrs['offset_Dec']
        flags = f['aper0/Flags'].value
        
    with h5py.File(filename) as f:
        
        head, tail = os.path.split(filename)
        tail = 'Red_' + tail.split('_')[-1]
        head = '/data2/talens/Feb2015LPE'
        outfilename = os.path.join(head, tail)
        
        print 'Writing to', outfilename
       
        with h5py.File(outfilename) as g:
        
            hdr = f['table_header']
            ascc = hdr['ascc'].value
            ra = hdr['ra'].value
            dec = hdr['dec'].value
            
            ra_idx, length, offset = make_idx(ra, (0,360), 13500) # hardcoded...
            dec_idx, length, offset_dec = make_idx(dec, (-90,90), 1440) # hardcoded...
            
            data = f['data']
            for i in range(len(ascc)):
                
                if dec_idx[i] != 1206: continue # FOR TESTING MAKES IT DO ONLY ONE DEC BIN!
                
                lst_idx = data[ascc[i]]['lstidx'].astype('int')
                flux0 = data[ascc[i]]['flux0']
                eflux0 = data[ascc[i]]['eflux0']
                
                ha_idx = lst_idx - ra_idx[i]
                ha_idx = np.mod(ha_idx, 13500) # hardcoded...
                
                ind_dec = dec_idx[i]-offset_dec
                ind_ha = ha_idx-offset_ha
                
                if (ind_dec < 0) | (ind_dec >= trans.shape[0]):
                    print 'Warning: star outside transmission map.'
                    
                    tcurve = np.full(len(flux0), fill_value=np.nan)
                    tflags = np.full(len(flux0), fill_value=2)
                
                elif np.any(ind_ha < 0) | np.any(ind_ha >= trans.shape[1]):
                    print 'Warning: datapoints outside transmission map'
                
                    badpoints, = np.where((ind_ha < 0) | (ind_ha >= trans.shape[1]))
                    ind_ha[badpoints] = 0
                    
                    tcurve = trans[ind_dec, ind_ha]
                    tflags = flags[ind_dec, ind_ha]
                    
                    tcurve[badpoints] = np.nan
                    tflags[badpoints] = 2
                
                else:
                    
                    tcurve = trans[ind_dec, ind_ha]
                    tflags = flags[ind_dec, ind_ha]
                
                tflux0 = flux0/tcurve
                etflux0 = eflux0/tcurve
        
                g.create_dataset(ascc[i]+'/tflux0', data=tflux0)
                g.create_dataset(ascc[i]+'/etflux0', data=etflux0)
                g.create_dataset(ascc[i]+'/tflags', data=tflags)
    
    return 0

filelist = glob.glob('/data2/mascara/LaPalma/201502??LPE/fLC/fLC_*')
filelist = np.sort(filelist)

for i in range(len(filelist)):
    print 'Processing', filelist[i]
    apply_transmission(filelist[i], '/data2/talens/Feb2015LPE/T_20150203LPE.hdf5')

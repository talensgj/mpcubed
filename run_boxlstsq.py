#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package.coordinates import grids

import filters
from boxlstsq_ms import boxlstsq_ms


def header(data):
    
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    
    for filename in data:
        with h5py.File(filename, 'r') as f:
            grp = f['header']
            ascc_ = grp['ascc'].value
            ra_ = grp['ra'].value
            dec_ = grp['dec'].value
            
            ascc = np.append(ascc, ascc_)
            ra = np.append(ra, ra_)
            dec = np.append(dec, dec_)
            
    ascc, args = np.unique(ascc, return_index=True)
    ra = ra[args]
    dec = dec[args]
    
    return ascc, ra, dec

def data_as_array(filename, ascc):
    
    lstseq = np.array([])
    staridx = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    
    with h5py.File(filename, 'r') as f:
        
        for i in range(len(ascc)):
        
            try:
                grp = f['data/' + ascc[i]]
            except:
                continue
            
            lstseq_ = grp['lstseq'].value
            jdmid_ = grp['jdmid'].value
            lst_ = grp['lst'].value
            mag_ = grp['mag0'].value
            emag_ = grp['emag0'].value
            nobs_ = grp['nobs'].value
            
            emag_ = emag_/np.sqrt(nobs_)
            
            select = (nobs_ == 50)
            lstseq_ = lstseq_[select]
            jdmid_ = jdmid_[select]
            lst_ = lst_[select]
            mag_ = mag_[select]
            emag_ = emag_[select]
            
            lstseq = np.append(lstseq, lstseq_)
            staridx = np.append(staridx, [i]*len(lstseq_))
            jdmid = np.append(jdmid, jdmid_)
            lst = np.append(lst, lst_)
            mag = np.append(mag, mag_)
            emag = np.append(emag, emag_)
    
    staridx = staridx.astype('int')
    lstseq, args, idx = np.unique(lstseq, return_index=True, return_inverse=True)
    jdmid = jdmid[args]
    lst = lst[args]
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = mag
    mag = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = emag
    emag = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)), dtype='bool')
    tmp[staridx, idx] = True
    mask = ~tmp
    
    return jdmid, lst, mag, emag, mask

def trend_filter(jdmid, lst, mag, emag):
    
    nstars = mag.shape[0]
    weights = np.where(emag > 1e-3, 1./emag**2, 0.)
    
    for i in range(nstars):
        #chisq, pars, fit = filters.harmonic(lst, mag[i], weights[i], 24., 8)
        chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag[i], weights[i], 180., 20)
        mag[i] = mag[i] - fit
    
    return mag

def main():
    
    data = ['/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPS/red0_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5']
    
    # Read the combined header.
    ascc, ra, dec = header(data)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    
    for i in range(372, hg.npix):
        
        if i not in np.unique(skyidx): continue
        
        select = (skyidx == i)
        
        jdmid = np.array([])
        lst = np.array([])
        mag = np.array([[]]*len(ascc[select]))
        emag = np.array([[]]*len(ascc[select]))
        mask = np.array([[]]*len(ascc[select]))
        
        for filename in data:
            
            # Read the data.
            jdmid_, lst_, mag_, emag_, mask_ = data_as_array(filename, ascc[select])
            
            if len(jdmid_) == 0: continue
            
            # Correct for trends.
            mag_ = trend_filter(jdmid_, lst_, mag_, emag_)
            
            jdmid = np.append(jdmid, jdmid_)
            lst = np.append(lst, lst_)
            mag = np.append(mag, mag_, axis=1)
            emag = np.append(emag, emag_, axis=1)
            mask = np.append(mask, mask_, axis=1)
        
        if len(jdmid) == 0: continue
        
        # Run the box least-squares search.
        weights = np.where(mask, 0, 1/emag**2)
        freq, dchisq, depth, hchisq, chisq0 = boxlstsq_ms(jdmid, mag.T, weights.T)
        
        if freq is None: continue
        
        # Save the results to file.
        blsfile = '/data2/talens/2015Q2/bls0_2015Q2_patch{:03d}.hdf5'.format(i)
        with h5py.File(blsfile) as f:
            grp = f.create_group('header')
            grp.create_dataset('ascc', data=ascc[select])
            grp.create_dataset('chisq0', data=chisq0)
            
            grp = f.create_group('data')
            grp.create_dataset('freq', data=freq)
            grp.create_dataset('dchisq', data=dchisq)
            grp.create_dataset('hchisq', data=hchisq)
            grp.create_dataset('depth', data=depth)
            
            grp = f.create_group('baseline')
            grp.create_dataset('jdmid', data=jdmid)
            grp.create_dataset('mask', data=mask)

    return

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import multiprocessing as mp

import matplotlib.pyplot as plt

from package.coordinates import grids

import filters
from boxlstsq_ms_dev import boxlstsq_ms, phase_duration


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

def data_mc(data, ascc, detrend=True):
    
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([[]]*len(ascc))
    emag = np.array([[]]*len(ascc))
    mask = np.array([[]]*len(ascc), dtype='bool')
    trend = np.array([[]]*len(ascc))
    
    for filename in data:
        
        # Read the data.
        jdmid_, lst_, mag_, emag_, mask_ = data_as_array(filename, ascc)
        
        if len(jdmid_) == 0: continue
        
        # Correct for trends.
        if detrend:
            mag_, trend_ = trend_filter(jdmid_, lst_, mag_, emag_)
        
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        mag = np.append(mag, mag_, axis=1)
        emag = np.append(emag, emag_, axis=1)
        mask = np.append(mask, mask_, axis=1)
        trend = np.append(trend, trend_, axis=1)
        
    return jdmid, lst, mag, emag, mask, trend

def trend_filter(jdmid, lst, mag, emag):
    
    nstars = mag.shape[0]
    weights = np.where(emag > 1e-3, 1./emag**2, 0.)
    trend = np.zeros(mag.shape)
    
    for i in range(nstars):
        #chisq, pars, fit = filters.harmonic(lst, mag[i], weights[i], 24., 8)
        chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag[i], weights[i], 180., 20)
        trend[i] = fit
        mag[i] = mag[i] - fit
    
    return mag, trend

def bls_core(skyidx, ascc, jdmid, mag, emag, mask):
    
    print 'Computing boxlstsq for skypatch', skyidx
    
    # Run the box least-squares search.
    weights = np.where(mask, 0, 1/emag**2)
    freq, chisq0, dchisq, hchisq, depth, epoch, duration, nt = boxlstsq_ms(jdmid, mag.T, weights.T)
    
    if freq is None:
        print 'freq = None', skyidx
        return 
    
    # Find the peak in the periodogram.
    args1 = np.argmax(dchisq, axis=0)
    args2 = np.arange(mag.shape[0])
    
    best_freq = freq[args1]
    best_chisq = chisq0 - dchisq[args1, args2]
    best_depth = depth[args1, args2]
    best_epoch = epoch[args1, args2]
    best_duration = duration[args1, args2]
    best_nt = nt[args1, args2]
    
    # Create flags.
    flag = np.zeros(mag.shape[0], dtype='int')
    
    # Check the best fit chi-square.
    nobs = np.sum(~mask, axis=1)
    quality = best_chisq/nobs
    args, = np.where(quality > 4)
    flag[args] = flag[args] + 1
        
    # Check the anti-transit ratio.
    tmp = dchisq*np.sign(depth)
    ratio = -np.amax(tmp, axis=0)/np.amin(tmp, axis=0)
    args, = np.where(ratio < 1.5)
    flag[args] = flag[args] + 2
    
    # Check the phase coverage.
    q = phase_duration(best_freq, 1., 1.)
    phase = np.outer(jdmid, best_freq)
    phase = np.mod(phase, 1)
    phase = np.sort(phase, axis=0)
    gapsizes = np.diff(phase, axis=0)
    gapsizes = np.vstack([gapsizes, 1. - np.ptp(phase, axis=0)])
    gapsizes = np.amax(gapsizes, axis=0)/q
    args, = np.where(gapsizes > 2.5)
    flag[args] = flag[args] + 4 
    
    # Check the number of observed transits.
    ntransit = best_nt*(319.1262613/(24*3600))/best_duration
    args, = np.where(ntransit < 3.)
    flag[args] = flag[args] + 8
    
    # Save the results to file.
    blsfile = '/data2/talens/2015Q2_vmag/bls0_vmag_2015Q2_patch{:03d}.hdf5'.format(skyidx)
    #blsfile = '/home/talens/MASCARA/bls_test.hdf5'
    with h5py.File(blsfile) as f:
        grp = f.create_group('header')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('chisq0', data=chisq0)
        grp.create_dataset('period', data=1/best_freq)
        grp.create_dataset('depth', data=best_depth)
        grp.create_dataset('epoch', data=best_epoch)
        grp.create_dataset('duration', data=best_duration)
        grp.create_dataset('nt', data=best_nt)
        grp.create_dataset('flag', data=flag)
        
        grp = f.create_group('data')
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq)
        grp.create_dataset('hchisq', data=hchisq)
        grp.create_dataset('depth', data=depth)
        grp.create_dataset('epoch', data=epoch)
        grp.create_dataset('duration', data=duration)
        grp.create_dataset('nt', data=nt)

    return

def main():
    
    data = ['/data2/talens/2015Q2_vmag/LPN/red0_vmag_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2_vmag/LPS/red0_vmag_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2_vmag/LPW/red0_vmag_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2_vmag/LPC/red0_vmag_2015Q2LPC.hdf5']
    
    # Read the combined header.
    ascc, ra, dec = header(data)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    
    jobs = []
    
    for i in range(hg.npix):
        
        # Could loop over np.unique(skyidx), bu that is harder to restart.
        if i not in np.unique(skyidx): continue
        
        # Read the skybin for all cameras.
        select = (skyidx == i)
        jdmid, lst, mag, emag, mask = data_mc(data, ascc[select])
        
        # Make sure there was data.
        if len(jdmid) == 0: continue
        
        # Do not run if the baseline falls short of 60 days.
        if np.ptp(jdmid) < 60.: continue
        
        jobs.append((i, ascc[select], jdmid, mag, emag, mask))
    
    print 'Found {} sky-patches to perfrom the box least-squares on.'.format(len(jobs))
    print
        
    pool = mp.Pool(processes = 6)
    for i in range(len(jobs)):
        pool.apply_async(bls_core, args = jobs[i])
    pool.close()
    pool.join()
    
    return

if __name__ == '__main__':
    main()

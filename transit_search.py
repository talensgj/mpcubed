#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import h5py
import numpy as np
import multiprocessing as mp

from package.coordinates import grids

import filters
from boxlstsq_ms_dev import boxlstsq_ms, phase_duration

def read_header(filelist):
    
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    
    for filename in filelist:
        
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

def read_data_array(filename, ascc):
    
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
    
    # Get the array indices of the data.
    staridx = staridx.astype('int')
    lstseq, args, idx = np.unique(lstseq, return_index=True, return_inverse=True)
    
    # Get the time arrays.
    jdmid = jdmid[args]
    lst = lst[args]
    
    # Cast the data to arrays.
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = mag
    mag = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = emag
    emag = tmp
    
    # Create a mask.
    tmp = np.zeros((len(ascc), len(lstseq)), dtype='bool')
    tmp[staridx, idx] = True
    mask = ~tmp
    
    return jdmid, lst, mag, emag, mask

def detrend(jdmid, lst, mag, emag):
    
    nstars = mag.shape[0]
    weights = np.where(emag > 1e-3, 1./emag**2, 0.)
    trend = np.zeros(mag.shape)
    
    for i in range(nstars):
        
        chisq, pars, trend[i] = filters.masc_harmonic(jdmid, lst, mag[i], weights[i], 180., 20)
    
    return trend

def read_data(filelist, ascc):
    
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([[]]*len(ascc))
    emag = np.array([[]]*len(ascc))
    mask = np.array([[]]*len(ascc), dtype='bool')
    trend = np.array([[]]*len(ascc))
    
    for filename in filelist:
        
        # Read the data.
        jdmid_, lst_, mag_, emag_, mask_ = read_data_array(filename, ascc)
        
        if len(jdmid_) == 0: continue
        
        # Detrend the lightcurves.
        trend_ = detrend(jdmid_, lst_, mag_, emag_)
        mag_ = mag_ - trend_
        
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        mag = np.append(mag, mag_, axis=1)
        emag = np.append(emag, emag_, axis=1)
        mask = np.append(mask, mask_, axis=1)
        trend = np.append(trend, trend_, axis=1)
        
    return jdmid, lst, mag, emag, mask, trend

def search_skypatch(skyidx, ascc, jdmid, mag, emag, mask):
    
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

def transit_search(filelist):
    
    # Read the combined header of the files.
    ascc, ra, dec = read_header(filelist)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)

    jobs = []
    
    for i in range(hg.npix):
        
        if i not in np.unique(skyidx): continue
        
        # Read the skybin for all cameras.
        select = (skyidx == i)
        jdmid, lst, mag, emag, mask, trend = read_data(filelist, ascc[select])
        
        # Make sure there was data.
        if len(jdmid) == 0: continue
        
        # Do not run if the baseline falls short of 60 days.
        if np.ptp(jdmid) < 60.: continue
        
        jobs.append((i, ascc[select], jdmid, mag, emag, mask))
    
    print 'Found {} sky-patches to perfrom the box least-squares on.'.format(len(jobs))
    print
        
    pool = mp.Pool(processes = 6)
    for i in range(len(jobs)):
        pool.apply_async(search_skypatch, args = jobs[i])
    pool.close()
    pool.join()


def main():
    return

if __name__ == '__main__':
    main()

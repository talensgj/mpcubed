#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os
import glob

import h5py
import numpy as np
import multiprocessing as mp

from package.coordinates import grids
from package.statistics import statistics

import detrend
import boxlstsq

def read_header(filelist):
    """ Read the combined header given reduced lightcurves."""
    
    # Create arrays.
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    vmag = np.array([])
    sptype = np.array([])
    jdmin = np.array([])
    jdmax = np.array([])
    
    for filename in filelist:
        
        # Read the headers.
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            ascc_ = grp['ascc'].value
            ra_ = grp['ra'].value
            dec_ = grp['dec'].value
            vmag_ = grp['vmag'].value
            sptype_ = grp['spectype'].value
            jdmin_ = grp['jdmin'].value
            jdmax_ = grp['jdmax'].value
            
        ascc = np.append(ascc, ascc_)
        ra = np.append(ra, ra_)
        dec = np.append(dec, dec_)
        vmag = np.append(vmag, vmag_)
        sptype = np.append(sptype, sptype_)
        jdmin = np.append(jdmin, jdmin_)
        jdmax = np.append(jdmax, jdmax_)
    
    # Obtain the unique entries.
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    ra = ra[args]
    dec = dec[args]
    vmag = vmag[args]
    sptype = sptype[args]
    jdmin = statistics.idxstats(idx, jdmin, statistic=np.amin)
    jdmax = statistics.idxstats(idx, jdmax, statistic=np.amax)
    
    return ascc, ra, dec, vmag, sptype, jdmin, jdmax

def find_ns(lstseq):
    """ Number of sampling points in LST, takes wrapping into account."""
    
    lstidx = lstseq%270
    option1 = np.ptp(lstidx) + 1
    
    lstidx = np.mod(lstidx + 135, 270)
    option2 = np.ptp(lstidx) + 1
    
    if (option2 >= option1):
        return option1, False
    else:
        return option2, True
    
def read_data_array(filename, ascc):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    lstseq = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    staridx = np.array([])
    ns = 2*np.ones((len(ascc), 2), dtype='int')
    wrap = np.zeros(len(ascc), dtype='bool')
    
    # Read the lightcurves.
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
            
            if (len(jdmid_) == 0):
                continue
            
            lstseq = np.append(lstseq, lstseq_)
            jdmid = np.append(jdmid, jdmid_)
            lst = np.append(lst, lst_)
            mag = np.append(mag, mag_)
            emag = np.append(emag, emag_)
            staridx = np.append(staridx, [i]*len(lstseq_))
            ns[i, 0] = np.ptp(lstseq_) + 1
            ns[i, 1], wrap[i] = find_ns(lstseq_)

    ns = np.maximum(ns, 2)

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
    
    return lstseq, jdmid, lst, mag, emag, mask, ns, wrap

def fit_trend(jdmid, lst, mag, emag, mask, ns, wrap):
    """ Detrend an array of lightcurves."""
    
    nstars = mag.shape[0]
    weights = np.where(mask, 0., 1./emag**2)
    trend = np.zeros(mag.shape)
    
    for i in range(nstars):
        
        if wrap[i]:
            x2 = np.mod(lst+12., 24.)-12.
        else:
            x2 = np.copy(lst)
        
        pars, trend[i], chisq = detrend.new_harmonic(jdmid, x2, mag[i], weights[i], ns[i])

    return trend

def read_data(filelist, ascc):
    """ Read the data from a list of reduced lightcurve files."""
    
    # Create arrays.
    lstseq = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([[]]*len(ascc))
    emag = np.array([[]]*len(ascc))
    mask = np.array([[]]*len(ascc), dtype='bool')
    trend = np.array([[]]*len(ascc))
    cam = np.array([])
    
    camid = {'LPN':0, 'LPE':1, 'LPS':2, 'LPW':3, 'LPC':4}
    
    for filename in filelist:

        # Read the data.
        lstseq_, jdmid_, lst_, mag_, emag_, mask_, ns_, wrap_ = read_data_array(filename, ascc)
        
        if len(jdmid_) == 0:
            continue
        
        # Detrend the lightcurves.
        trend_ = fit_trend(jdmid_, lst_, mag_, emag_, mask_, ns_, wrap_)
        mag_ = mag_ - trend_
        
        lstseq = np.append(lstseq, lstseq_)
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        cam = np.append(cam, [camid[filename[-8:-5]]]*len(jdmid_))
        mag = np.append(mag, mag_, axis=1)
        emag = np.append(emag, emag_, axis=1)
        mask = np.append(mask, mask_, axis=1)
        trend = np.append(trend, trend_, axis=1)
        
    return lstseq, jdmid, lst, mag, emag, mask, trend, cam

def search_skypatch(skyidx, ascc, jdmid, mag, emag, mask, blsfile):
    """ Perform the box least-squares and flag non-detections."""
    
    print 'Computing boxlstsq for skypatch', skyidx
    
    # Run the box least-squares search.
    weights = np.where(mask, 0., 1./emag**2)
    freq, chisq0, dchisq, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag.T, weights.T)
    
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
    q = boxlstsq.phase_duration(best_freq, 1., 1.)
    phase = np.outer(jdmid, best_freq)
    phase = np.mod(phase, 1.)
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
    with h5py.File(blsfile) as f:
        
        grp = f.create_group('header')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('chisq0', data=chisq0, dtype='float32')
        grp.create_dataset('period', data=1./best_freq)
        grp.create_dataset('depth', data=best_depth, dtype='float32')
        grp.create_dataset('epoch', data=best_epoch)
        grp.create_dataset('duration', data=best_duration, dtype='float32')
        grp.create_dataset('nt', data=best_nt, dtype='uint32')
        grp.create_dataset('flag', data=flag, dtype='uint32')
        
        grp = f.create_group('data')
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq, dtype='float32')
        grp.create_dataset('hchisq', data=hchisq, dtype='float32')
        grp.create_dataset('depth', data=depth, dtype='float32')
        grp.create_dataset('epoch', data=epoch)
        grp.create_dataset('duration', data=duration, dtype='float32')
        grp.create_dataset('nt', data=nt, dtype='uint32')

    return

def search_skypatch_mp(jobs, nprocs):
    """ Use multiprocessing to perform the transit search."""
    
    pool = mp.Pool(processes = nprocs)
    for i in range(len(jobs)):
        pool.apply_async(search_skypatch, args = jobs[i])
    pool.close()
    pool.join()
    
    return
    
def transit_search(filelist, outdir, name, nprocs=6):
    """ Perform detrending and transit search given reduced lightcurves."""
    
    # Read the combined header of the files.
    ascc, ra, dec, vmag, sptype, jdmin, jdmax = read_header(filelist)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)

    # Determine the maximum baseline in each grid cell.
    blgrid = np.zeros(hg.npix)
    skyuni, idx = np.unique(skyidx, return_inverse=True)
    jdmax = statistics.idxstats(idx, jdmax, statistic=np.amax)
    jdmin = statistics.idxstats(idx, jdmin, statistic=np.amin)
    blgrid[skyuni] =  jdmax - jdmin
        
    jobs = []
    for i in [59]:#range(hg.npix):
        
        # Do not run if the baseline falls short of 60 days.
        if (blgrid[i] < 60.):
            print 'Skipping patch {}/{}, insufficient baseline.'.format(i, hg.npix) 
            continue
        
        # Read the stars in the skypatch.
        select = (skyidx == i)
        lstseq, jdmid, lst, mag, emag, mask, trend, cam = read_data(filelist, ascc[select])
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            continue
        
        # Filename for the output file. 
        blsfile = 'bls0_{}_patch{:03d}.hdf5'.format(name, i)
        blsfile = os.path.join(outdir, blsfile)
        
        # Check if multiprocessing was set.
        if (nprocs == 1):
            search_skypatch(i, ascc[select], jdmid, mag, emag, mask, blsfile)
        else:
            jobs.append((i, ascc[select], jdmid, mag, emag, mask, blsfile))
        
        # Compute the periodgrams in the joblist.
        if (len(jobs) == 5*nprocs):
            search_skypatch_mp(jobs, nprocs)
            jobs = []
    
    # Compute any remaining periodograms.
    if (len(jobs) != 0):
        search_skypatch_mp(jobs, nprocs)
    
    # Print a message indicating succes.
    print
    print 'Succesfully ran the boxlstsq on:'
    for filename in filelist:
        print ' ', filename
    print
    print 'The results were writen to:'
    print ' ', outdir
        
    return

def ensure_dir(path):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    return

def main():
    
    data = glob.glob('/data3/talens/2015Q?/LP?/red0_vmag_2015Q?LP?.hdf5')
    data = np.sort(data)

    print data
    
    outdir = '/data3/talens/boxlstsq/test'
    ensure_dir(outdir)
    
    transit_search(data, outdir, 'vmag_2015Q1234', nprocs=16)
    
    return

if __name__ == '__main__':
    main()

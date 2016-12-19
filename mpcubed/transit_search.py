#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os
import glob

import h5py
import numpy as np
import multiprocessing as mp

import misc
import boxlstsq
from coordinates import grids
from statistics import statistics
from systematics import detrend
import transit_statistics as stats

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
            
            try:
                mag_ = grp['mag0'].value
                emag_ = grp['emag0'].value
            except:
                mag_ = grp['mag1'].value
                emag_ = grp['emag1'].value
                
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
            ns[i, 1], wrap[i] = misc.find_ns(lstseq_)

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
        
        freq1, freq2, pars1, pars2, fit1, fit2, chisq = detrend.psf_variations(jdmid, x2, mag[i], weights[i], ns[i])
        trend[i] = fit1 + fit2

    return trend

def search_skypatch(skyidx, ascc, jdmid, mag, emag, mask, blsfile):
    """ Perform the box least-squares and flag non-detections."""
    
    print 'Computing boxlstsq for skypatch', skyidx
    
    # Run the box least-squares search.
    weights = np.where(mask, 0., 1./emag**2)
    freq, chisq0, dchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag.T, weights.T, (~mask).T)
    
    if freq is None:
        print 'freq = None', skyidx
        return 
    
    # Best fit parameters and statistics from the periodogram.
    freq0, dchisq0, epoch0, depth0, duration0, sde, atr = stats.bls_crit(freq, dchisq, epoch, depth, duration)
    
    # Create arrays.
    nstars = len(freq0)
    gap = np.zeros(nstars)
    sym = np.zeros(nstars)
    ntr = np.zeros(nstars)
    ntp = np.zeros(nstars)
    mst = np.zeros(nstars)
    eps = np.zeros(nstars)
    sne = np.zeros(nstars)
    sw = np.zeros(nstars)
    sr = np.zeros(nstars)
    snp = np.zeros(nstars)

    # Statistics from the lightcurve.
    for i in range(nstars):
        try:
            gap[i], sym[i], ntr[i], ntp[i], mst[i], eps[i], sne[i], sw[i], sr[i], snp[i] = stats.lc_crit(jdmid[~mask[i]], mag[i, ~mask[i]], emag[i, ~mask[i]], freq0[i], epoch0[i], depth0[i], duration0[i])
        except:
            pass
    
    # Save the results to file.
    with h5py.File(blsfile) as f:
        
        grp = f.create_group('header')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('chisq0', data=chisq0, dtype='float32')

        grp.create_dataset('period', data=1./freq0)
        grp.create_dataset('epoch', data=epoch0)
        grp.create_dataset('depth', data=depth0, dtype='float32')
        grp.create_dataset('duration', data=duration0, dtype='float32')
        
        grp.create_dataset('sde', data=sde, dtype='float32')
        grp.create_dataset('atr', data=atr, dtype='float32')
        grp.create_dataset('gap', data=gap, dtype='float32')
        grp.create_dataset('sym', data=sym, dtype='float32')
        grp.create_dataset('ntr', data=ntr, dtype='int32')
        grp.create_dataset('ntp', data=ntp, dtype='int32')
        grp.create_dataset('mst', data=mst, dtype='float32')
        grp.create_dataset('eps', data=eps, dtype='float32')
        grp.create_dataset('sne', data=sne, dtype='float32')
        grp.create_dataset('sw', data=sw, dtype='float32')
        grp.create_dataset('sr', data=sr, dtype='float32')
        grp.create_dataset('snp', data=snp, dtype='float32')

        grp = f.create_group('data')
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq, dtype='float32')

    return
    
def search_skypatch_mp(queue):
    """ Use multiprocessing to perform the transit search."""
    
    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:
            search_skypatch(*item)
    
    return

def transit_search(filelist, name, patches=None, outdir='/data3/talens/boxlstsq', nprocs=6):
    """ Perform detrending and transit search given reduced lightcurves."""
    
    print 'Trying to run the box least-squares on:' 
    for filename in filelist:
        print ' ', filename
    
    # Create directories for the output.
    outdir = os.path.join(outdir, name)
    blsdir = os.path.join(outdir, 'bls')
    misc.ensure_dir(blsdir)
    
    if (len(os.listdir(blsdir)) > 0):
        print 'Error: the output directory {} is not empty.'.format(outdir)
        exit()
    else:
        print 'Writing results to:', outdir
    
    # Write the a file containg the reduced data files used.
    fname = os.path.join(outdir, 'data.txt')
    np.savetxt(fname, filelist, fmt='%s')
    
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
    
    # Dtermine which skypatches to run.
    if patches is None:
        patches = range(hg.npix)
        
    if np.any(np.array(patches) >= hg.npix):
        print 'Error: patches greater than {} do not exist.'.format(hg.npix)
        exit()   
        
    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, search_skypatch_mp, (the_queue,))
        
    for i in patches:
        
        # Do not run if the baseline falls short of 60 days.
        if (blgrid[i] < 60.):
            print 'Skipping patch {}/{}, insufficient baseline.'.format(i, hg.npix) 
            continue
        
        # Read the stars in the skypatch.
        select = (skyidx == i)
        lstseq, jdmid, lst, mag, emag, mask, trend, cam = read_data(filelist, ascc[select])
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            print 'Skipping patch {}/{}, no good data found.'.format(i, hg.npix) 
            continue
        
        # Filename for the output file. 
        blsfile = 'bls0_{}_patch{:03d}.hdf5'.format(name, i)
        blsfile = os.path.join(blsdir, blsfile)
        
        the_queue.put((i, ascc[select], jdmid, mag, emag, mask, blsfile))
    
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
    
    return

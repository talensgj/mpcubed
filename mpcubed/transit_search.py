#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os
import glob

import h5py
import numpy as np
import multiprocessing as mp

import io, misc
import boxlstsq
from coordinates import grids
from systematics import detrend
import transit_statistics as stats

###############################################################################
### Functions for reading the reduced lightcurves.
###############################################################################

def fit_trend(jdmid, lst, mag, emag, ns, wrap):
    """ Detrend an array of lightcurves."""
    
    weights = 1./emag**2
        
    if wrap:
        x2 = np.mod(lst+12., 24.)-12.
    else:
        x2 = np.copy(lst)
    
    freq1, freq2, pars1, pars2, fit1, fit2, chisq = detrend.psf_variations(jdmid, x2, mag, weights, ns)

    return fit1 + fit2

def read_header(filelist):
    """ Read the combined header given reduced lightcurves."""
    
    # Create arrays.
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    vmag = np.array([])
    sptype = np.array([])
    
    for filename in filelist:
        
        f = io.PhotFile(filename)
        stars = f.read_stars(fields=['ascc', 'ra', 'dec', 'vmag', 'spectype'], grpname='header')
            
        ascc = np.append(ascc, stars['ascc'])
        ra = np.append(ra, stars['ra'])
        dec = np.append(dec, stars['dec'])
        vmag = np.append(vmag, stars['vmag'])
        sptype = np.append(sptype, stars['spectype'])
    
    # Obtain the unique entries.
    ascc, args = np.unique(ascc, return_index=True)
    ra = ra[args]
    dec = dec[args]
    vmag = vmag[args]
    sptype = sptype[args]
    
    return ascc, ra, dec, vmag, sptype
    
def detrended_lightcurves(filename, ascc, aper=0):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    staridx = np.array([])
    lstseq = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    trend = np.array([])
    
    # Read the lightcurves.
    f = io.PhotFile(filename)
    data = f.read_lightcurves(ascc=ascc, verbose=False)
        
    # Detrend and flatten the lightcurves.
    for i in range(len(ascc)):       
        
        try:
            lc = data[ascc[i]]
        except:
            continue
        
        mask = (lc['nobs'] == 50)
        lc = lc[mask]
        
        if (len(lc) < 1):
            continue
        
        # Detrend the lightcurves.
        ns = [0,0]
        ns[0] = np.ptp(lc['lstseq']) + 1
        ns[1], wrap = misc.find_ns(lc['lstseq'])
        ns = np.maximum(ns, 2)
        trend_ = fit_trend(lc['jdmid'], lc['lst'], lc['mag%i'%aper], lc['emag%i'%aper], ns, wrap)
        
        # Add the results to the arrays.
        staridx = np.append(staridx, [i]*len(lc))
        lstseq = np.append(lstseq, lc['lstseq'])
        jdmid = np.append(jdmid, lc['jdmid'])
        lst = np.append(lst, lc['lst'])
        mag = np.append(mag, lc['mag%i'%aper] - trend_)
        emag = np.append(emag, lc['emag%i'%aper])
        trend = np.append(trend, trend_)
    
    # Convert the lightcurves to 2D arrays.
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
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = trend
    trend = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)), dtype='bool')
    tmp[staridx, idx] = True
    mask = ~tmp
    
    return jdmid, lst, mag, emag, trend, mask

def read_data(filelist, ascc):
    """ Read the data from a list of reduced lightcurve files."""
    
    # Create arrays.
    jdmid = np.empty((0,))
    lst = np.empty((0,))
    mag = np.empty((len(ascc), 0))
    emag = np.empty((len(ascc), 0))
    trend = np.empty((len(ascc), 0))
    mask = np.empty((len(ascc), 0), dtype='bool')
    
    for filename in filelist:

        # Read the data.
        jdmid_, lst_, mag_, emag_, trend_, mask_ = detrended_lightcurves(filename, ascc)
        
        if (len(jdmid_) < 1):
            continue
        
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        mag = np.append(mag, mag_, axis=1)
        emag = np.append(emag, emag_, axis=1)
        trend = np.append(trend, trend_, axis=1)
        mask = np.append(mask, mask_, axis=1)
        
    return jdmid, lst, mag, emag, trend, mask

###############################################################################
### Functions for running the box least-squares algorithm.
###############################################################################

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
    ascc, ra, dec, vmag, sptype = read_header(filelist)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    
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
        
        # Read the stars in the skypatch.
        select = (skyidx == i)
        jdmid, lst, mag, emag, trend, mask = read_data(filelist, ascc[select])
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            print 'Skipping patch {}/{}, no good data found.'.format(i, hg.npix) 
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(jdmid) < 60.):
            print 'Skipping patch {}/{}, insufficient baseline.'.format(i, hg.npix) 
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

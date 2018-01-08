#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os
import glob

import h5py
import numpy as np
import multiprocessing as mp

from . import io, misc, boxlstsq
from .coordinates import grids
from .systematics import detrend
from . import transit_statistics as stats
from .statistics import statistics

###############################################################################
### Functions for reading the reduced lightcurves.
###############################################################################

def read_header(filelist):
    """ Read the combined header given reduced lightcurves."""
    
    # Create arrays.
    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    vmag = np.array([])
    
    for filename in filelist:
        
        f = io.PhotFile(filename)
        stars = f.read_stars(fields=['ascc', 'ra', 'dec', 'vmag'])
            
        ascc = np.append(ascc, stars['ascc'])
        ra = np.append(ra, stars['ra'])
        dec = np.append(dec, stars['dec'])
        vmag = np.append(vmag, stars['vmag'])
    
    # Obtain the unique entries.
    ascc, args = np.unique(ascc, return_index=True)
    ra = ra[args]
    dec = dec[args]
    vmag = vmag[args]
    
    return ascc, ra, dec, vmag
    
def detrended_lightcurves(filename, ascc, aper=0, method='legendre'):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    staridx = np.array([])
    lstseq = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    trend = np.array([])
    
    # Set the aperture.
    magstr = 'mag{}'.format(aper)
    emagstr = 'emag{}'.format(aper)
    
    # Read the lightcurves.
    f = io.PhotFile(filename)
    data = f.read_lightcurves(ascc=ascc, verbose=False)
        
    # Detrend and flatten the lightcurves.
    for i in range(len(ascc)):       
        
        # See if there was data.
        try:
            lc = data[ascc[i]]
        except:
            continue
        
        # Select data binned from 50 exposures.
        mask = (lc['nobs'] == 50)
        lc = lc[mask]
        
        # Check that there are at least 2 points.
        if (len(lc) < 2):
            continue

        # Get the julian date.
        try:
            jdmid_ = lc['jdmid'] # La Palma
        except:
            jdmid_ = lc['jd'] # bRing, La Silla
        
        # Detrend the lightcurves.
        if method is 'none':
            
            trend_ = np.zeros(len(jdmid_))
        
        elif method == 'legendre':        
        
            mat, fit0, fit1, fit2 = detrend.detrend_legendre(jdmid_, lc['lst'], lc['sky'], lc[magstr], lc[emagstr])        
            trend_ = fit0 + fit1 + fit2
            
        elif method == 'snellen':
            
            fit0, fit1, fit2 = detrend.detrend_snellen(jdmid_, lc['lstseq'], lc['sky'], lc[magstr], lc[emagstr])
            trend_ = fit0 + fit1 + fit2
            
        elif method == 'fourier':        
        
            ns = [0,0]
            ns[0] = np.ptp(lc['lstseq']) + 1
            ns[1], wrap = misc.find_ns(lc['lstseq'])
            ns = np.maximum(ns, 2)
            
            mat, fit0, fit1 = detrend.detrend_fourier(jdmid_, lc['lst'], lc[magstr], lc[emagstr], ns, wrap)
            trend_ = fit0 + fit1 
            
        else:
            raise ValueError('Unknown detrending method "{}"'.format(method))   
        
        # Add the results to the arrays.
        staridx = np.append(staridx, [i]*len(lc))
        lstseq = np.append(lstseq, lc['lstseq'])
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lc['lst'])
        mag = np.append(mag, lc[magstr] - trend_)
        emag = np.append(emag, lc[emagstr])
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

def read_data(filelist, ascc, aper=0, method='legendre'):
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
        jdmid_, lst_, mag_, emag_, trend_, mask_ = detrended_lightcurves(filename, ascc, aper, method)
        
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
    
    # Convert the uncertainties to weights.
    with np.errstate(divide='ignore'):
        weights = np.where(mask, 0., 1./emag**2)
        
    # Run the box least-squares search.    
    freq, chisq0, dchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag.T, weights.T, (~mask).T)
    
    # Create arrays.
    nstars = len(ascc)

    # Parameters.
    freq0 = np.ones(nstars)
    epoch0 = np.zeros(nstars)
    depth0 = np.zeros(nstars)
    duration0 = np.zeros(nstars)
    
    # Box-least squares statistics.
    sde = np.zeros(nstars)
    atr = np.zeros(nstars)    
    
    # Lightcurve statistics.
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

    for i in range(nstars):
        
        # Best-fit parameters.
        arg = np.argmax(dchisq[:,i])
        freq0[i] = freq[arg]
        epoch0[i] = epoch[arg,i]
        depth0[i] = depth[arg,i]
        duration0[i] = duration[arg,i]        

        if (sum(~mask[i]) > 1) & (dchisq[arg,i] > 0):
            # Statistics.
            sde[i], atr[i] = stats.boxlstsq_criteria(dchisq[:,i], depth[:,i])
            gap[i], sym[i], ntr[i], ntp[i], mst[i], eps[i], sne[i], sw[i], sr[i], snp[i] = stats.lightcurve_criteria(jdmid[~mask[i]], mag[i, ~mask[i]], emag[i, ~mask[i]], freq0[i], epoch0[i], depth0[i], duration0[i])
    
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

def transit_search(filelist, name, patches=None, aper=0, method='legendre', outdir='/data3/talens/boxlstsq', nprocs=6):
    """ Perform detrending and transit search given reduced lightcurves."""
    
    print 'Trying to run the box least-squares on aperture {} of:'.format(aper) 
    for filename in filelist:
        print ' ', filename
    
    # Create directories for the output.
    outdir = os.path.join(outdir, name + method)
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
    ascc, ra, dec, vmag = read_header(filelist)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    
    # Dtermine which skypatches to run.
    if patches is None:
        patches = np.unique(skyidx)
        
    if np.any(np.array(patches) >= hg.npix):
        print 'Error: patches greater than {} do not exist.'.format(hg.npix)
        exit()   
        
    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, search_skypatch_mp, (the_queue,))
        
    for idx in patches:
        
        # Read the stars in the skypatch.
        select = (skyidx == idx)
        jdmid, lst, mag, emag, trend, mask = read_data(filelist, ascc[select], aper, method)
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            print 'Skipping skypatch {}, no good data found.'.format(idx) 
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(jdmid) < 60.):
            print 'Skipping skypatch {}, insufficient baseline.'.format(idx) 
            continue        
        
        # Barycentric correction.
        ra, dec = hg.idx2radec(idx)
        jdbar = misc.barycentric_dates(jdmid, ra, dec)

        # Filename for the output file. 
        blsfile = 'bls{}_{}{}_patch{:03d}.hdf5'.format(aper, name, method, idx)
        blsfile = os.path.join(blsdir, blsfile)
        
        the_queue.put((idx, ascc[select], jdbar, mag, emag, mask, blsfile))
    
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
    
    return

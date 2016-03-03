#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os

import h5py
import numpy as np
import multiprocessing as mp

from package.coordinates import grids

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
    
    for filename in filelist:
        
        # Read the headers.
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            ascc_ = grp['ascc'].value
            ra_ = grp['ra'].value
            dec_ = grp['dec'].value
            vmag_ = grp['vmag'].value
            sptype_ = grp['spectype'].value
            
            ascc = np.append(ascc, ascc_)
            ra = np.append(ra, ra_)
            dec = np.append(dec, dec_)
            vmag = np.append(vmag, vmag_)
            sptype = np.append(sptype, sptype_)
    
    # Obtain the unique entries.
    ascc, args = np.unique(ascc, return_index=True)
    ra = ra[args]
    dec = dec[args]
    vmag = vmag[args]
    sptype = sptype[args]
    
    return ascc, ra, dec, vmag, sptype

def read_data_array(filename, ascc):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    lstseq = np.array([])
    staridx = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    
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

def fit_trend(jdmid, lst, mag, emag, mask):
    """ Detrend an array of lightcurves."""
    
    nstars = mag.shape[0]
    weights = np.where(mask, 0., 1./emag**2)
    trend = np.zeros(mag.shape)
    
    for i in range(nstars):
        
        pars, trend[i], chisq = detrend.masc_harmonic(jdmid, lst, mag[i], weights[i], 180., 21, 24., 6)
    
    return trend

def read_data(filelist, ascc):
    """ Read the data from a list of reduced lightcurve files."""
    
    # Create arrays.
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([[]]*len(ascc))
    emag = np.array([[]]*len(ascc))
    mask = np.array([[]]*len(ascc), dtype='bool')
    trend = np.array([[]]*len(ascc))
    
    for filename in filelist:
        
        # Read the data.
        jdmid_, lst_, mag_, emag_, mask_ = read_data_array(filename, ascc)
        
        if len(jdmid_) == 0:
            continue
        
        # Detrend the lightcurves.
        trend_ = fit_trend(jdmid_, lst_, mag_, emag_, mask_)
        mag_ = mag_ - trend_
        
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        mag = np.append(mag, mag_, axis=1)
        emag = np.append(emag, emag_, axis=1)
        mask = np.append(mask, mask_, axis=1)
        trend = np.append(trend, trend_, axis=1)
        
    return jdmid, lst, mag, emag, mask, trend

def search_skypatch(skyidx, ascc, jdmid, mag, emag, mask, blsfile):
    
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
        grp.create_dataset('depth', data=best_depth)
        grp.create_dataset('epoch', data=best_epoch)
        grp.create_dataset('duration', data=best_duration)
        grp.create_dataset('nt', data=best_nt, dtype='uint32')
        grp.create_dataset('flag', data=flag, dtype='uint32')
        
        grp = f.create_group('data')
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq, dtype='float32')
        grp.create_dataset('hchisq', data=hchisq, dtype='float32')
        grp.create_dataset('depth', data=depth)
        grp.create_dataset('epoch', data=epoch)
        grp.create_dataset('duration', data=duration)
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
    ascc, ra, dec, vmag, sptype = read_header(filelist)
    
    # Divide the stars in groups of neighbouring stars.
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)

    skypatches = np.unique(skyidx)
    jobs = []
    
    for i in range(hg.npix):
        
        # Check that the patch was observed.
        if i not in skypatches:
            continue
        
        # Read the stars in the skypatch.
        select = (skyidx == i)
        jdmid, lst, mag, emag, mask, trend = read_data(filelist, ascc[select])
        
        # Make sure there was data.
        if (len(jdmid) == 0):
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(jdmid) < 60.):
            continue
        
        blsfile = 'bls0_{}_patch{:03d}.hdf5'.format(name, i)
        blsfile = os.path.join(outdir, blsfile)
        
        # Check if multiprocessing was set.
        if (nprocs == 1):
            search_skypatch(i, ascc[select], jdmid, mag, emag, mask, blsfile)
        else:
            jobs.append((i, ascc[select], jdmid, mag, emag, mask, blsfile))
        
        # Compute the periodgrams in the joblist.
        if (len(jobs) == nprocs):
            search_skypatch_mp(jobs, nprocs)
            jobs = []
    
    # Compute any remaining periodograms.
    if (len(jobs) != 0):
        search_skypatch_mp(jobs, nprocs)
        
    return


def main():
    
    data = ['/data2/talens/2015Q2_vmag/LPN/red0_vmag_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2_vmag/LPS/red0_vmag_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2_vmag/LPW/red0_vmag_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2_vmag/LPC/red0_vmag_2015Q2LPC.hdf5']
    
    outdir = '/data2/talens/2015Q2_vmag'
    
    transit_search(data, outdir, 'vmag_2015Q2', nprocs=6)
    
    return

if __name__ == '__main__':
    main()

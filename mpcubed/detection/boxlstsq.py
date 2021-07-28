#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np
from numba import jit
import multiprocessing as mp

from .. import io, misc, statistics, models
from . import detrend, criteria
from ..calibration import grids

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m]  

np.random.seed(19910909)

###############################################################################
### The Box Least-Squares Algorithm (Kovacs+ 2002, Ofir 2014)
###############################################################################

def phase_duration(freq, M, R):
    """ Compute the expected transit duration."""

    freq = freq/SecInDay
    M = M*Msun
    R = R*Rsun

    q = (2.*np.pi)**(2./3)/np.pi * R/(G*M)**(1./3) * (freq)**(2./3)
    
    return q

def freq_max(M, R):
    """ Compute the maximum orbital frequency from the Roche limit."""
    
    M = M*Msun
    R = R*Rsun
    
    fmax = 1./(2.*np.pi) * np.sqrt(G*M/(3.*R)**3)
    
    fmax = fmax*SecInDay
    
    return fmax
    
def freqs(fmin, fmax, S, OS, M, R):
    """ Compute the ideal sample frequencies. """
    
    fmin = fmin/SecInDay
    fmax = fmax/SecInDay
    S = S*SecInDay
    M = M*Msun
    R = R*Rsun
    
    A = (2.*np.pi)**(2./3)/np.pi * R/(G*M)**(1./3) * 1/(S * OS)
    xmax = (fmax**(1./3) - fmin**(1./3) + A/3.)*3./A
    xmax = np.rint(xmax).astype('int')
    xvec = np.arange(1, xmax + 1)
    C = fmin**(1./3) - A/3.
    freq = (A/3.*xvec + C)**3.
    
    freq = freq*SecInDay
    
    return freq

@jit(nopython=True)
def mycounts(idx, weights, nbins):
    
    hist = np.zeros((nbins, weights.shape[1]))
    for i in range(len(idx)):
        hist[idx[i]] += weights[i]
        
    return hist

def cumsum_to_grid(time, bin_edges, wflux, weights, mask):
    
    ndim = wflux.ndim    
    
    if (ndim == 1):
        wflux = np.atleast_2d(wflux).T
        weights = np.atleast_2d(weights).T
        mask = np.atleast_2d(mask).T

    # Get the sizes of the arrays.
    nbins = len(bin_edges) 
    npoints, nstars = wflux.shape 

    # Find the bin indices.
    idx = np.searchsorted(bin_edges, time, 'right') 
    
    # Compute the sums in the bins.
    r_cum = mycounts(idx, weights, nbins)
    s_cum = mycounts(idx, wflux, nbins) 
    n_cum = mycounts(idx, mask, nbins)
    
    # Compute the cumulative sum.
    r_cum = np.cumsum(r_cum, axis=0)
    s_cum = np.cumsum(s_cum, axis=0)
    n_cum = np.cumsum(n_cum, axis=0)

    if (ndim == 1):
        r_cum = r_cum[:,0]
        s_cum = s_cum[:,0]
        n_cum = n_cum[:,0]

    return r_cum, s_cum, n_cum

def boxlstsq(time, flux, weights, mask, exp_time=320./86400., **options):
    """ Compute the box least-square periodogram for multiple stars."""
    
    time = np.asarray(time)
    flux = np.asarray(flux)
    weights = np.asarray(weights)
    
    if (len(time) != len(flux)):
        raise ValueError('time and flux must have the same first dimension.')
    
    if (np.ndim(flux) > 2):
        raise ValueError('flux must be a 1d or a 2d array.')
    
    if (np.shape(flux) != np.shape(weights)):
        raise ValueError('flux and weights must be the same shape.')
    
    # Get the options.
    M = options.pop('M', 1.)
    R = options.pop('R', 1.)
    ES = options.pop('ES', 3)
    Qmin = options.pop('Qmin', 3)
    Qmax = options.pop('Qmax', 3)
    NumSteps = Qmin*Qmax
    freq = options.pop('freq', None)
    fmin = options.pop('fmin', None)
    fmax = options.pop('fmax', None)
    OS = options.pop('OS', 3)
    
    # Determine the baseline of the data.
    S = np.ptp(time)
    
    # Frequency sampling.
    if freq is not None:
        print 'Using user specified frequencies, ignoring fmin, fmax, OS.'
    
    else:

        # Defaults for fmin and fmax.
        fmin_def = np.maximum(12./S, 1./30) # Reasonable limits on longest period.
        fmax_def = freq_max(M, R) # Smallest orbit determined by the Roche limit.
        
        if fmin is None:
            fmin = fmin_def
        elif (fmin < fmin_def):
            print 'Warning: minimum frequency contains <2 events or exceeds 30 days.'
            
        if fmax is None:
            fmax = fmax_def
        elif (fmax > fmax_def):
            print 'Warning: maximum frequency is inside the Hill radius.'
                
        if (fmax < fmin):
            raise ValueError('fmax must be larger than fmin.')
        
        # Compute optimal sampling frequencies in range fmin, fmax.
        freq = freqs(fmin, fmax, S, OS, M, R)
    
    # Use the transit duration to determine the sampling grid.
    q = phase_duration(freq, M, R)
    nbins = np.ceil((Qmin*ES)/q).astype('int')
    
    # Prepare the data for the loop.
    time0 = np.amin(time)
    time = time - time0
    mask = mask.astype('float64')
    t = np.nansum(weights, axis=0)  # Sum of weights.
    with np.errstate(invalid='ignore'):
        flux = flux - np.nansum(weights*flux, axis=0)/t  # Subtract average
    wflux = weights*flux
    chisq0 = np.nansum(weights*flux**2., axis=0)  # Best fit constant model.

    # Create arrays.
    nfreq = len(freq)
    
    if (np.ndim(flux) == 1):
        res_shape = (nfreq,)
    else:
        res_shape = (nfreq, flux.shape[1])
    
    dchisq_inc = np.zeros(res_shape)
    dchisq_dec = np.zeros(res_shape)
    
    depth = np.zeros(res_shape)
    epoch = np.zeros(res_shape)
    duration = np.zeros(res_shape)
    nt = np.zeros(res_shape)
    
    # Loop over frequency domain.
    for i in xrange(nfreq):
        
        # Phase fold the lightcurve. 
        phase = np.mod(freq[i]*time, 1.)

        # Sum the data on a regular grid.
        bins = np.linspace(0., 1., nbins[i] + 1)
        r_cum, s_cum, n_cum = cumsum_to_grid(phase, bins, wflux, weights, mask)
        nmax = n_cum[-1]

        # Extend the grid to account for all epochs.
        bins = np.append(bins, bins[1:NumSteps*ES] + bins[-1])
        n_cum = np.append(n_cum, n_cum[1:NumSteps*ES] + n_cum[-1], axis=0)
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES] + r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES] + s_cum[-1], axis=0)
        
        # Calculate the indices of the start and end of transit for all epochs
        # and durations.
        i1, i2 = np.indices((nbins[i], NumSteps))
        i2 = i1 + ES*(i2 + 1)
       
        i1 = i1.ravel()
        i2 = i2.ravel()
       
        # Compute arrays of n, r and s values.
        n = n_cum[i2] - n_cum[i1]
        r = r_cum[i2] - r_cum[i1]
        s = s_cum[i2] - s_cum[i1]

        # Compute the epoch and duration values.
        epoch_tmp = (bins[i1] + bins[i2])/(2*freq[i])
        duration_tmp = (bins[i2] - bins[i1])/freq[i]

        # Compute the depth and dchisq values.
        with np.errstate(invalid='ignore', divide='ignore'):
            depth_tmp = s*t/(r*(t - r))
            dchisq_tmp = s**2*t/(r*(t - r))

        # Set dchisq to zero if all data is out-of-transit or in-transit.
        dchisq_tmp[n < 1] = 0
        dchisq_tmp[(nmax - n) < 1] = 0

        depth_tmp[n < 1] = 0
        depth_tmp[(nmax - n) < 1] = 0

        # Set dchisq to zero if too few points are in-transit or too few out-of-transit.
        nmin = duration_tmp/exp_time
        if dchisq_dec.ndim > 1:

            depth_tmp[n < nmin[:, None]] = 0
            dchisq_tmp[n < nmin[:, None]] = 0
            depth_tmp[(nmax - n) < nmin[:, None]] = 0
            dchisq_tmp[(nmax - n) < nmin[:, None]] = 0

        else:

            depth_tmp[n < nmin] = 0
            dchisq_tmp[n < nmin] = 0
            depth_tmp[(nmax - n) < nmin] = 0
            dchisq_tmp[(nmax - n) < nmin] = 0

        # Select the best solution.
        args_inc_1d = np.nanargmax(dchisq_tmp*np.sign(depth_tmp), axis=0)
        args_dec_1d = np.nanargmin(dchisq_tmp*np.sign(depth_tmp), axis=0)
        
        if (np.ndim(flux) == 2):
            args_2d = np.arange(flux.shape[1])
            args_inc_2d = (args_inc_1d, args_2d)
            args_dec_2d = (args_dec_1d, args_2d)
        else:
            args_inc_2d = args_inc_1d
            args_dec_2d = args_dec_1d
        
        dchisq_inc[i] = dchisq_tmp[args_inc_2d]
        dchisq_dec[i] = dchisq_tmp[args_dec_2d]

        depth[i] = depth_tmp[args_dec_2d]
        epoch[i] = epoch_tmp[args_dec_1d]
        duration[i] = duration_tmp[args_dec_1d]
        nt[i] = n[args_dec_2d]
        
    epoch = epoch + time0
        
    return freq, chisq0, dchisq_inc, dchisq_dec, depth, epoch, duration, nt

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

def _read_data(filename, ascc, aper=0):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    lc = []
    staridx = []
    
    # Read the lightcurves.
    f = io.PhotFile(filename)
    data = f.read_lightcurves(ascc=ascc, verbose=False)
        
    for i in range(len(ascc)):       
        
        # See if there was data.
        try:
            tmp = data[ascc[i]]
        except:
            continue

        # Select data binned from >45 exposures.
        tmp = tmp[tmp['nobs'] > 45]
        
        # Check that there are at least 2 points.
        if (len(tmp) < 2):
            continue
                    
        # Add the results to the arrays.
        lc.append(tmp)
        staridx.append(i*np.ones_like(tmp['lstseq']))

    if (len(lc) == 0):
        return [], [], 0

    lc = np.concatenate(lc)
    staridx = np.concatenate(staridx)
    
    # Create a 1D array of time.
    lstseq, args, idx = np.unique(lc['lstseq'], return_index=True, return_inverse=True)

    dtype = [('lstseq', 'uint32'), ('jd', 'float64'), ('lst', 'float64')]

    time = np.recarray((len(args),), dtype=dtype)
    
    time['lstseq'] = lstseq
    try: time['jd'] = lc['jd'][args]
    except: time['jd'] = lc['jdmid'][args]
    time['lst'] = lc['lst'][args]

    # Cast the data to a 2D array.
    dtype = [('lstseq', 'uint32'), ('jd', 'float64'), ('lst', 'float64'), ('mag', 'float64'), ('emag', 'float64'), ('x', 'float64'), ('y', 'float64'), ('sky', 'float64'), ('trend', 'float64'), ('mask', 'bool')]

    lc2d = np.recarray((len(args), len(ascc)), dtype=dtype)
    lc2d[:,:] = 0
    
    lc2d['lstseq'][idx,staridx] = lc['lstseq']
    try: lc2d['jd'][idx,staridx] = lc['jd']
    except: lc2d['jd'][idx,staridx] = lc['jdmid']
    lc2d['lst'][idx,staridx] = lc['lst']
    lc2d['mag'][idx,staridx] = lc['mag{}'.format(aper)]
    lc2d['emag'][idx,staridx] = lc['emag{}'.format(aper)]
    lc2d['x'][idx,staridx] = lc['x']
    lc2d['y'][idx,staridx] = lc['y']
    lc2d['sky'][idx,staridx] = lc['sky']
    lc2d['mask'][idx,staridx] = True

    return time, lc2d, len(args)

def read_data(filelist, ascc, aper=0):
    """ Read the data from a list of reduced lightcurve files."""
    
    # Create arrays.
    time = []
    lc2d = []
    nobs = []
    
    for filename in filelist:

        # Read the data.
        time_, lc2d_, nobs_ = _read_data(filename, ascc, aper)
        
        if (nobs_ == 0):
            continue
        
        time.append(time_)
        lc2d.append(lc2d_)
        nobs.append(nobs_)

    if len(time) > 0:
        time = np.concatenate(time)
        lc2d = np.concatenate(lc2d)
        nobs = np.stack(nobs)

    return time, lc2d, nobs

def clip_outliers(lc, nsig=3., floor=0.05):
            
    mask = lc['mask']
    
    if np.any(mask):
    
        # Compute the median and MAD of the data.
        m0 = np.median(lc['mag'][mask] - lc['trend'][mask])
        m1 = statistics.mad(lc['mag'][mask] - lc['trend'][mask])
        
        # Remove outliers from the data.
        cutoff = np.maximum(nsig*m1, floor)
        lc['mask'] = lc['mask'] & (np.abs(lc['mag'] - lc['trend'] - m0) < cutoff)

    return lc

def transit_params(nsets, P_range=[1., 5.], p_choices=np.sqrt([0.005, 0.01, 0.02]), b_choices=[0.0, 0.5], rho_choices=[0.4, 0.9, 1.4]):
    
    names = ['T0', 'P', 'T14', 'p', 'b']
    formats = ['float64', 'float64', 'float32', 'float32', 'float32']
    pars = np.recarray((nsets,), names=names, formats=formats)
    
    pars['P'] = np.random.rand(nsets)*(P_range[1] - P_range[0]) + P_range[0]
    pars['T0'] = np.random.rand(nsets)*pars['P']
    pars['p'] = np.random.choice(p_choices, nsets)
    pars['b'] = np.random.choice(b_choices, nsets)

    rho = np.random.choice(rho_choices, nsets)
    axis = models.dens2axis(rho, pars['P'])
    pars['T14'] = models.axis2T14(axis, pars['P'], pars['p'], pars['b'])
    
    return pars

def inject_transit(lc, lc_pars):
    
    mask = lc['mask']
    
    if np.any(mask):
        
        # Compute the model lightcurve.
        phase, model = models.transit_model(lc['jd'][mask], lc_pars.tolist())
        
        # Add the model to the data.
        lc['mag'][mask] = lc['mag'][mask] + model
        
    return lc

###############################################################################
### Functions for running the box least-squares algorithm.
###############################################################################

def search_skypatch(item, ascc, time, lc2d, nobs, method, blsfile, inj_pars=None):
    """ Perform the box least-squares and flag non-detections."""
    
    print 'Computing boxlstsq for:', item
    
    for i in range(lc2d.shape[1]):
            
        # Inject transit signals in the lightcurves.
        if (inj_pars is not None) & (i > 0):
            lc2d[:,i] = inject_transit(lc2d[:,i], inj_pars[i])
            
        # Perform the secondary calibration.
        lc2d[:,i] = detrend.remove_trend(lc2d[:,i], nobs, method=method)
        
        # Mask outlying data points.
        lc2d[:,i] = clip_outliers(lc2d[:,i])
    
    # Apply the secondary calibration.
    lc2d['mag'] = lc2d['mag'] - lc2d['trend']
    
    # Convert the uncertainties to weights.
    with np.errstate(divide='ignore'):
        weights = np.where(lc2d['mask'], 1./lc2d['emag']**2, 0.)
        
    # Run the box least-squares search.    
    freq, chisq0, dchisq_inc, dchisq_dec, depth, epoch, duration, nt = boxlstsq(time['jd'], -lc2d['mag'], weights, lc2d['mask'])
    
    # Create arrays.
    nstars = len(ascc)

    names = ['epoch', 'period', 'duration', 'depth']
    formats = ['float64', 'float64', 'float64', 'float64']
    box_pars = np.recarray((nstars,), names=names, formats=formats)
    box_pars[:] = 0
    
    names = ['sde', 'atr', 'gap', 'sym', 'ntr', 'ntp', 'mst', 'eps', 'sne', 'sw', 'sr', 'snp']
    formats = ['float32', 'float32', 'float32', 'float32', 'int32', 'int32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32'] 
    blscrit = np.recarray((nstars,), names=names, formats=formats)
    blscrit[:] = 0

    # Find best-fit parameters and evaluate quality criteria. 
    for i in range(nstars):
        
        # Best-fit parameters.
        arg = np.argmax(dchisq_dec[:,i])
        box_pars[i] = epoch[arg,i], 1./freq[arg], duration[arg,i], depth[arg,i]
        
        lc = lc2d[:,i]
        mask = lc['mask']
        
        # Quality criteria.
        if np.any(mask) & (dchisq_dec[arg,i] > 0):
            
            sde, atr = criteria.boxlstsq_criteria(dchisq_dec[:,i], dchisq_inc[:,i])
            gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp = criteria.lightcurve_criteria(time['jd'][mask], -lc['mag'][mask], lc['emag'][mask], box_pars[i])
    
            blscrit[i] = sde, atr, gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp
            
    # Save the results to file.
    io.write_boxlstsq(blsfile, ascc, chisq0, box_pars, blscrit, freq, dchisq_inc, dchisq_dec, inj_pars)

    return
    
def search_skypatch_mp(queue):
    """ Use multiprocessing to perform the transit search."""
    
    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:

            # Unpack the arguments.
            item, method, tmpfile, blsfile, inj_pars = item

            # Read the temporary file.
            npzfile = np.load(tmpfile)
            ascc, time, lc2d, nobs = npzfile['ascc'], npzfile['time'], npzfile['lc2d'], npzfile['nobs']

            # Remove the temporary file.
            os.remove(tmpfile)

            # Call the search algorithm.
            search_skypatch(item, ascc, time, lc2d, nobs, method, blsfile, inj_pars)
    
    return

def run_boxlstsq(filelist, name, aper=0, method='legendre', declims=[-90.,90.], outdir='/data3/talens/boxlstsq', nprocs=6, injection=False, nstars=5000, nsignals=11):
    """ Perform detrending and transit search given reduced lightcurves."""

    print 'Trying to run the box least-squares on aperture {} of:'.format(aper) 
    for filename in filelist:
        print ' ', filename
    
    # Create directories for the output.
    if injection:
        name = 'inj_' + name
        
    outdir = os.path.join(outdir, name + method)
    blsdir = os.path.join(outdir, 'bls')
    io.ensure_dir(blsdir)
    
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
                
    if not injection:
    
        # Divide the stars in groups of neighbouring stars.
        hg = grids.HealpixGrid(8)
        skyidx = hg.radec2idx(ra, dec)
        
        items = np.unique(skyidx)
        
        ra_, dec_ = hg.idx2radec(items)
        mask = (dec_ >= declims[0]) & (dec_ <= declims[1])
        items = items[mask]
            
    else:
        
        # Select stars for injection.
        mask = (dec >= declims[0]) & (dec <= declims[1])
        p = np.where(mask, 1./vmag**2, 0.)
        p = p/np.sum(p)
        
        items = np.random.choice(len(ascc), nstars, replace=False, p=p)
    
    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, search_skypatch_mp, (the_queue,))
        
    for item in items:
        
        if not injection:
            ascc_ = ascc[skyidx == item]
        else:
            ascc_ = [ascc[item]]*nsignals
            
            # Generate model parameters.
            inj_pars = transit_params(nsignals)
        
        # Read the stars in the skypatch.
        time, lc2d, nobs = read_data(filelist, ascc_, aper=aper)
    
        # Make sure there was data.
        if (len(time) == 0):
            print 'Skipping {}, no good data found.'.format(item) 
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(time['jd']) < 60.):
            print 'Skipping {}, insufficient baseline.'.format(item) 
            continue    
            
        # Compute the barycentric julian date.
        if not injection:
            ra_, dec_ = hg.idx2radec(item)
        else:
            ra_, dec_ = ra[item], dec[item]
            
        time['jd'] = misc.barycentric_dates(time['jd'], ra_, dec_)
        
        # Filenames for the output file and temporary data file.
        if not injection:
            tmpfile = 'tmp_patch{:03d}.npz'.format(item)
            blsfile = 'bls{}_{}{}_patch{:03d}.hdf5'.format(aper, name, method, item)
        else:
            tmpfile = 'tmp_ascc{}.npz'.format(ascc[item])
            blsfile = 'bls{}_{}{}_ascc{}.hdf5'.format(aper, name, method, ascc[item])

        tmpfile = os.path.join(blsdir, tmpfile)
        blsfile = os.path.join(blsdir, blsfile)

        # Save the data arrays to a temporary file.
        # Data volumes are too large to put directly in the queue.
        np.savez(tmpfile, ascc=ascc_, time=time, lc2d=lc2d, nobs=nobs)

        if not injection:
            the_queue.put((item, method, tmpfile, blsfile, None))
        else:
            the_queue.put((ascc[item], method, tmpfile, blsfile, inj_pars))
    
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
    
    # Create table for the database.
    io.table_boxlstsq(outdir)
    
    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Run the box least-squares on a collection of reduced data.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to process')
    parser.add_argument('name', type=str,
                        help ='the name of this box least-squares run')
    parser.add_argument('-a', '--aper', type=int, default=0,
                        help ='the aperture to use', dest='aper')
    parser.add_argument('-m', '--method', type=str, default='legendre',
                        help ='the detrending method', dest='method')
    parser.add_argument('-d', '--declims', type=float, nargs=2, default=[-90.,90.],
                        help='the declination range to process', dest='declims')
    parser.add_argument('-o', '--outdir', type=str, default='/data3/talens/boxlstsq',
                        help ='the location to write the results', dest='outdir')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    parser.add_argument('--injection', action='store_true', 
                        help='when set perform a signal recovery test', dest='injection')
    parser.add_argument('--nstars', type=int, default=5000,
                        help='the number of stars to use for injection', dest='nstars')
    parser.add_argument('--nsignals', type=int, default=11,
                        help='the number of signals to inject per star', dest='nsignals')
    
    args = parser.parse_args()
    
    run_boxlstsq(args.files, args.name, args.aper, args.method, args.declims, args.outdir, args.nprocs, args.injection, args.nstars, args.nsignals)
    
    return
    
if __name__ == '__main__':
    
    main()
    
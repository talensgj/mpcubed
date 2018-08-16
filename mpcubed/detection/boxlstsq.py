#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np
import multiprocessing as mp

from .. import io, misc, statistics
from . import detrend, criteria
from ..calibration import grids

from exofit import transit

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

def cumsum_to_grid(time, bin_edges, flux, weights, mask):
    
    ndim = flux.ndim    
    
    if (ndim == 1):
        flux = np.atleast_2d(flux).T
        weights = np.atleast_2d(weights).T
        mask = np.atleast_2d(mask).T

    # Get the sizes of the arrays.
    nbins = len(bin_edges) 
    npoints, nstars = flux.shape 

    # Find the bin indices.
    idx = np.searchsorted(bin_edges, time, 'right') 
    
    # Create arrays.
    r_cum = np.zeros((nbins, nstars))
    s_cum = np.zeros((nbins, nstars))
    n_cum = np.zeros((nbins, nstars))
    
    # Compute the sums in the bins.
    for i in range(nstars):
        r_cum[:,i] = np.bincount(idx, weights[:,i], minlength=nbins)
        s_cum[:,i] = np.bincount(idx, weights[:,i]*flux[:,i], minlength=nbins) 
        n_cum[:,i] = np.bincount(idx, mask[:,i], minlength=nbins)
    
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
        fmin_def = 12./S # Three transits in temporal baseline of ground based survey.
        fmax_def = freq_max(M, R) # Smallest orbit determined by the Roche limit.
        
        if fmin is None:
            fmin = fmin_def
        elif (fmin < fmin_def):
            print 'Warning: minimum frequency contains <2 events.'
            
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
    t = np.nansum(weights, axis=0) # Sum of weights.
    with np.errstate(invalid='ignore'):
        flux = flux - np.nansum(weights*flux, axis=0)/t # Subtract average
    chisq0 = np.nansum(weights*flux**2., axis=0) # Best fit constant model.
    
    # Create arrays.
    nfreq = len(freq)
    
    if (np.ndim(flux) == 1):
        res_shape = (nfreq,)
    else:
        res_shape = (nfreq, flux.shape[1])
    
    dchisq = np.zeros(res_shape)
    
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
        r_cum, s_cum, n_cum = cumsum_to_grid(phase, bins, flux, weights, mask)
        
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
        
        epoch_tmp = (bins[i1] + bins[i2])/(2*freq[i])
        duration_tmp = (bins[i2] - bins[i1])/freq[i]
        
        # Find the best fit.
        with np.errstate(invalid='ignore', divide='ignore'):
            dchisq_tmp = s**2*t/(r*(t - r))
        
        # Set dchisq to zero if all data is out-of-transit or in-transit.
        nmax = np.sum(mask, axis=0)
        dchisq_tmp[n < 1] = 0
        dchisq_tmp[(nmax - n) < 1] = 0
        
        # Set dchisq to zero if too few points are in-transit.
        nmin = duration_tmp/exp_time
        if (dchisq.ndim > 1):
            dchisq_tmp[n < nmin[:,None]] = 0
        else:
            dchisq_tmp[n < nmin] = 0

        # Select the best solution.
        args_1d = np.nanargmax(dchisq_tmp, axis=0)
        
        if (np.ndim(flux) == 2):
            args_2d = np.arange(flux.shape[1])
            args_2d = (args_1d, args_2d)
        else:
            args_2d = args_1d

        n = n[args_2d]
        r = r[args_2d]
        s = s[args_2d]
        
        dchisq[i] = dchisq_tmp[args_2d]
        
        with np.errstate(invalid='ignore', divide='ignore'):
            depth[i] = s*t/(r*(t - r))
            
        epoch[i] = epoch_tmp[args_1d]
        duration[i] = duration_tmp[args_1d]
        nt[i] = n
        
    epoch = epoch + time0
        
    return freq, chisq0, dchisq, depth, epoch, duration, nt

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
  
class GrowArray(object):
    
    def __init__(self, capacity, growth_factor=2, dtype=None):
        
        self.data = np.zeros((capacity,), dtype=dtype)
        self.capacity = capacity
        self.dtype = dtype
        self.growth_factor = growth_factor
        self.size = 0
        
        return
    
    def _extend(self):
        
        self.capacity *= self.growth_factor
            
        newdata = np.zeros((self.capacity,), dtype=self.dtype)
        newdata[:self.size] = self.data[:self.size]
        self.data = newdata  
        
        return
    
    def append(self, x):
        
        nx = x.size
        
        while (self.size + nx) > self.capacity:
            self._extend()

        self.data[self.size:self.size+nx] = x
        self.size += nx
        
        return
    
    def finalize(self):
        
        return self.data[:self.size]

def cast_to_2d(i, j, values, shape):
    
    arr = np.zeros(shape, dtype=values.dtype)
    arr[i,j] = values
    
    return arr

def detrended_lightcurves(filename, ascc, vmag, aper=0, method='legendre', options={}, inj_pars=None):
    """ Read the lightcurves of a group of stars."""
    
    # Create arrays.
    staridx = GrowArray(5000, dtype='uint32')
    lstseq = GrowArray(5000, dtype='uint32')
    jdmid = GrowArray(5000)
    lst = GrowArray(5000)
    mag = GrowArray(5000)
    emag = GrowArray(5000)
    trend = GrowArray(5000)
    mask = GrowArray(5000, dtype='bool')
    
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

        # Select data binned from >45 exposures.
        lc = lc[lc['nobs'] > 45]
        
        # Check that there are at least 2 points.
        if (len(lc) < 2):
            continue

        if f.get_siteid() != 'LP':
            lc[magstr] = lc[magstr] - vmag[i]

        # Get the julian date.
        try:
            jdmid_ = lc['jdmid'] # La Palma
        except:
            jdmid_ = lc['jd'] # bRing, La Silla
        
        if (inj_pars is not None) & (i > 0):
            phase, model = transit.circular_model(jdmid_, inj_pars[i])
        else:
            model = 0.
        
        # Detrend the lightcurves.
        if method is 'none':
            
            trend_ = np.zeros(len(jdmid_))
            mask_ = np.ones_like(trend_, dtype='bool')
        
        elif method == 'legendre':        
        
            mat, fit0, fit1, fit2 = detrend.detrend_legendre(jdmid_, lc['lst'], lc['sky'], lc[magstr] + model, lc[emagstr], **options)        
            trend_ = fit0 + fit1 + fit2
            mask_ = np.ones_like(trend_, dtype='bool')
            
        elif method == 'snellen':
            
            fit0, fit1, fit2, mask_ = detrend.detrend_snellen(jdmid_, lc['lstseq'], lc['x'], lc['y'], lc['sky'], lc[magstr] + model, lc[emagstr], **options)
            trend_ = fit0 + fit1 + fit2
            
        elif method == 'fourier':        
        
            ns = [0,0]
            ns[0] = np.ptp(lc['lstseq']) + 1
            ns[1], wrap = misc.find_ns(lc['lstseq'])
            ns = np.maximum(ns, 2)
            
            mat, fit0, fit1 = detrend.detrend_fourier(jdmid_, lc['lst'], lc[magstr] + model, lc[emagstr], ns, wrap, **options)
            trend_ = fit0 + fit1 
            mask_ = np.ones_like(trend_, dtype='bool')
            
        else:
            raise ValueError('Unknown detrending method "{}"'.format(method))   
        
        # Add the results to the arrays.
        staridx.append(i*np.ones_like(lc['lstseq']))
        lstseq.append(lc['lstseq'])
        jdmid.append(jdmid_)
        lst.append(lc['lst'])
        mag.append(lc[magstr] + model - trend_)
        emag.append(lc[emagstr])
        trend.append(trend_)
        mask.append(mask_)
    
    staridx = staridx.finalize()
    lstseq = lstseq.finalize()
    jdmid = jdmid.finalize()
    lst = lst.finalize()
    mag = mag.finalize()
    emag = emag.finalize()
    trend = trend.finalize()
    mask = mask.finalize()

    # Convert the lightcurves to 2D arrays.
    lstseq, args, idx = np.unique(lstseq, return_index=True, return_inverse=True)

    jdmid = jdmid[args]
    lst = lst[args]
    
    shape = (len(ascc), len(lstseq))
    mag = cast_to_2d(staridx, idx, mag, shape)
    emag = cast_to_2d(staridx, idx, emag, shape)
    trend = cast_to_2d(staridx, idx, trend, shape)
    mask = cast_to_2d(staridx, idx, mask, shape)
    
    return jdmid, lst, mag, emag, trend, mask

def read_data(filelist, ascc, vmag, aper=0, method='legendre', options={}, inj_pars=None):
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
        jdmid_, lst_, mag_, emag_, trend_, mask_ = detrended_lightcurves(filename, ascc, vmag, aper, method, options, inj_pars)
        
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

def search_skypatch(skyidx, ascc, jdmid, mag, emag, mask, blsfile, inj_pars):
    """ Perform the box least-squares and flag non-detections."""
    
    print 'Computing boxlstsq for skypatch', skyidx
    
    # Convert the uncertainties to weights.
    with np.errstate(divide='ignore'):
        weights = np.where(mask, 1./emag**2, 0.)
        
    # Run the box least-squares search.    
    freq, chisq0, dchisq, depth, epoch, duration, nt = boxlstsq(jdmid, mag.T, weights.T, mask.T)
    
    # Create arrays.
    nstars = len(ascc)

    names = ['period', 'epoch', 'depth', 'duration']
    formats = ['float64', 'float64', 'float64', 'float64']
    boxpars = np.recarray((nstars,), names=names, formats=formats)
    boxpars[:] = 0
    
    names = ['sde', 'atr', 'gap', 'sym', 'ntr', 'ntp', 'mst', 'eps', 'sne', 'sw', 'sr', 'snp']
    formats = ['float32', 'float32', 'float32', 'float32', 'int32', 'int32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32'] 
    blscrit = np.recarray((nstars,), names=names, formats=formats)
    blscrit[:] = 0

    # Find best-fit parameters and evaluate quality criteria. 
    for i in range(nstars):
        
        # Best-fit parameters.
        arg = np.argmax(dchisq[:,i])
        boxpars[i] = 1./freq[arg], epoch[arg,i], depth[arg,i], duration[arg,i]
        
        # Quality criteria.
        if (sum(mask[i]) > 1) & (dchisq[arg,i] > 0):
            
            jdmid_, mag_, emag_ = jdmid[mask[i]], mag[i, mask[i]], emag[i, mask[i]]
            
            sde, atr = criteria.boxlstsq_criteria(dchisq[:,i], depth[:,i])
            gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp = criteria.lightcurve_criteria(jdmid_, mag_, emag_, boxpars[i])
    
            blscrit[i] = sde, atr, gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp
            
    # Save the results to file.
    io.write_boxlstsq(blsfile, ascc, chisq0, boxpars, blscrit, freq, dchisq, inj_pars)

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

def run_boxlstsq(filelist, name, patches=None, aper=0, method='legendre', outdir='/data3/talens/boxlstsq', nprocs=6):
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
        jdmid, lst, mag, emag, trend, mask = read_data(filelist, ascc[select], vmag[select], aper=aper, method=method)
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            print 'Skipping skypatch {}, no good data found.'.format(idx) 
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(jdmid) < 60.):
            print 'Skipping skypatch {}, insufficient baseline.'.format(idx) 
            continue    
        
        # Mask outliers.
        for i in range(mag.shape[0]):
            
            m0 = np.median(mag[i][mask[i]])
            m1 = statistics.mad(mag[i][mask[i]])
            
            cutoff = np.maximum(3*m1, 0.05)
            mask[i] = mask[i] & (np.abs(mag[i] - m0) < cutoff)

        # Barycentric correction.
        ra, dec = hg.idx2radec(idx)
        jdbar = misc.barycentric_dates(jdmid, ra, dec)

        # Filename for the output file. 
        blsfile = 'bls{}_{}{}_patch{:03d}.hdf5'.format(aper, name, method, idx)
        blsfile = os.path.join(blsdir, blsfile)
        
        the_queue.put((idx, ascc[select], jdbar, mag, emag, mask, blsfile, None))
    
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
    
    return

def transit_params(nsets, P_range=[1., 5.], p_choices=np.sqrt([0.005, 0.01, 0.02]), b_choices=[0.0, 0.5], rho_choices=[0.4, 0.9, 1.4]):
    
    names = ['T0', 'P', 'T14', 'p', 'b']
    formats = ['float64', 'float64', 'float32', 'float32', 'float32']
    pars = np.recarray((nsets,), names=names, formats=formats)
    
    pars['P'] = np.random.rand(nsets)*(P_range[1] - P_range[0]) + P_range[0]
    pars['T0'] = np.random.rand(nsets)*pars['P']
    pars['p'] = np.random.choice(p_choices, nsets)
    pars['b'] = np.random.choice(b_choices, nsets)

    rho = np.random.choice(rho_choices, nsets)
    axis = transit.dens2axis(rho, pars['P'])
    pars['T14'] = transit.axis2duration(axis, pars['P'], pars['p'], pars['b'])
    
    return pars

def run_injection(filelist, name, nobj=5000, ninj=11, aper=0, method='legendre', outdir='/data3/talens/boxlstsq', nprocs=6):
    """ Perform detrending and transit search given reduced lightcurves."""
    
    print 'Trying to run the box least-squares on aperture {} of:'.format(aper) 
    for filename in filelist:
        print ' ', filename
    
    # Create directories for the output.
    outdir = os.path.join(outdir, 'inj_' + name + method)
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
    
    # Select stars for injection.
    p = 1./vmag**2
    p = p/np.sum(p)
    indices = np.random.choice(len(ascc), nobj, replace=False, p=p)
    
    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, search_skypatch_mp, (the_queue,))
    
    for idx in indices:
        
        ascc_ = [ascc[idx]]*ninj
        vmag_ = [vmag[idx]]*ninj
        
        # Generate model parameters.
        inj_pars = transit_params(ninj)
        
        # Read the star with injected models.
        jdmid, lst, mag, emag, trend, mask = read_data(filelist, ascc_, vmag_, aper=aper, method=method, inj_pars=inj_pars)
    
        # Make sure there was data.
        if (len(jdmid) == 0):
            print 'Skipping star {}, no good data found.'.format(ascc[idx]) 
            continue
        
        # Do not run if the baseline falls short of 60 days.
        if (np.ptp(jdmid) < 60.):
            print 'Skipping star {}, insufficient baseline.'.format(ascc[idx]) 
            continue    
        
        # Mask outliers.
        for i in range(mag.shape[0]):
            
            m0 = np.median(mag[i][mask[i]])
            m1 = statistics.mad(mag[i][mask[i]])
            
            cutoff = np.maximum(3*m1, 0.05)
            mask[i] = mask[i] & (np.abs(mag[i] - m0) < cutoff)

        # Barycentric correction.
        jdbar = misc.barycentric_dates(jdmid, ra[idx], dec[idx])

        # Filename for the output file. 
        blsfile = 'bls{}_{}{}_ascc{}.hdf5'.format(aper, name, method, ascc[idx])
        blsfile = os.path.join(blsdir, blsfile)
        
        the_queue.put((ascc[idx], ascc_, jdbar, mag, emag, mask, blsfile, inj_pars))
    
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
    
    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Run the box least-squares on a collection of reduced data.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to process')
    parser.add_argument('name', type=str,
                        help ='the name of this box least-squares run')
    parser.add_argument('-p', '--patches', type=int, nargs='+', default=None,
                        help ='the sky patch(es) to process', dest='patches')
    parser.add_argument('-a', '--aper', type=int, default=0,
                        help ='the aperture to use', dest='aper')
    parser.add_argument('-m', '--method', type=str, default='legendre',
                        help ='detrending method', dest='method')
    parser.add_argument('-o', '--outdir', type=str, default='/data3/talens/boxlstsq',
                        help ='the location to write the results', dest='outdir')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()
    
    run_boxlstsq(args.files, args.name, args.patches, args.aper, args.method, args.outdir, args.nprocs)
    
    return
    
if __name__ == '__main__':
    
    main()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m]  

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

def boxlstsq(time, flux, weights, mask, **options):
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
    
    # Dtermine the baseline of the data.
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
            print 'Error: fmax must be larger than fmin.'
            return None, None, None, None, None, None, None, None
        
        # Compute optimal sampling frequencies in range fmin, fmax.
        freq = freqs(fmin, fmax, S, OS, M, R)
    
    # Use the transit duration to determine the sampling grid.
    q = phase_duration(freq, M, R)
    nbins = np.ceil((Qmin*ES)/q)
    
    # Prepare the data for the loop.
    time0 = np.amin(time)
    time = time - time0
    t = np.nansum(weights, axis=0) # Sum of weights.
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

        nmax = np.amax(n_cum, axis=0)
        
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
       
        # Compute arrays of n, r, s and q values.
        n = n_cum[i2] - n_cum[i1]
        r = r_cum[i2] - r_cum[i1]
        s = s_cum[i2] - s_cum[i1]
        
        epoch_tmp = (bins[i1] + bins[i2])/(2*freq[i])
        duration_tmp = (bins[i2] - bins[i1])/freq[i]
        
        # Find the best fit.
        dchisq_tmp = s**2*t/(r*(t - r))
        
        nmin = duration_tmp/(320./(24.*3600.))
        if (dchisq.ndim > 1):
            dchisq_tmp[n < nmin[:,None]] = 0
        else:
            dchisq_tmp[n < nmin] = 0

        dchisq_tmp[(nmax - n) < 1] = 0
        
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
        
        depth[i] = s*t/(r*(t - r))
        epoch[i] = epoch_tmp[args_1d]
        duration[i] = duration_tmp[args_1d]
        nt[i] = n
        
    epoch = epoch + time0
        
    return freq, chisq0, dchisq, depth, epoch, duration, nt

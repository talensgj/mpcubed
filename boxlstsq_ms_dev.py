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
    
def cumsum_to_grid(time, bin_edges, weights):

    nt = len(time)
    nedges = len(bin_edges)
    nw = len(weights)
    
    if (np.diff(bin_edges) < 0).any():
        raise AttributeError('bins must increase monotonically.')
    
    if (nt != nw):
        raise ValueError('weights should have the same first dimension as a.')
    
    if (weights.ndim > 1):
        nstars = np.shape(weights)[1]
        n = np.zeros((nedges,nstars), weights.dtype)
        zero = np.array([0]*nstars, weights.dtype)
    else:
        n = np.zeros((nedges,), weights.dtype)
        zero = np.array(0, weights.dtype)
    
    block = 65536
    for i in np.arange(0, len(time), block):
        sa = time[i:i+block]
        sw = weights[i:i+block]
        tmp = np.cumsum(sw, axis=0)
        cw = np.vstack(([zero,], tmp))
        bin_index = np.r_[sa.searchsorted(bin_edges[:-1], 'left'), sa.searchsorted(bin_edges[-1], 'right')]
        n += cw[bin_index]
    
    return n

def boxlstsq_ms(time, flux, weights, **options):
    """ Compute the box least-square periodogram for multiple stars."""
    
    time = np.asarray(time)
    flux = np.asarray(flux)
    weights = np.asarray(weights)
    
    if (len(time) != len(flux)):
        raise ValueError('time and flux must have the same first dimension.')
    
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
        elif (fmax < fmax_def):
            print 'Warning: maximum frequency is inside the Hill radius.'
                
        if (fmax < fmin):
            print 'Error: fmax must be larger than fmin.'
            return None, None, None, None, None
        
        # Compute optimal sampling frequencies in range fmin, fmax.
        freq = freqs(fmin, fmax, S, OS, M, R)
    
    print 'Using M=%.2f Msun and R=%.2f Rsun for the stellar parameters.'%(M, R)
    print 'Sampling frequencies between fmin=%.2e, fmax=%.2e per day with OS=%i'%(fmin, fmax, OS)
    print 'Sampling transit epochs at %.2f time(s) the smallest sampled transit.'%(1./ES)
    print 'Sampling transits between %.2f and %.2f time(s) the expected duration in %i steps.'%(1./Qmin, Qmax, NumSteps)
    
    # Use the transit duration to determine the sampling grid.
    q = phase_duration(freq, M, R)
    nbins = np.ceil((Qmin*ES)/q)
    
    # Prepare the data for the loop.
    time = time - np.amin(time)
    t = np.nansum(weights, axis=0) # Sum of weights.
    flux = flux - np.nansum(weights*flux, axis=0)/t # Subtract average
    chisq0 = np.nansum(weights*flux**2., axis=0) # Best fit constant model.
    
    # Create arrays.
    nfreq = len(freq)
    
    dchisq = np.zeros((nfreq, flux.shape[1]))
    hchisq = np.zeros((nfreq, flux.shape[1]))
    
    depth = np.zeros((nfreq, flux.shape[1]))
    epoch = np.zeros((nfreq, flux.shape[1]))
    duration = np.zeros((nfreq, flux.shape[1]))
    nt = np.zeros((nfreq, flux.shape[1]))
    
    # Loop over frequency domain.
    for i in xrange(nfreq):
        
        # Phase fold the lightcurve. 
        phase = np.mod(time*freq[i], 1.)
        sort = np.argsort(phase)
        
        phase_fold = np.take(phase, sort, axis=0)
        flux_fold = np.take(flux, sort, axis=0)
        weights_fold = np.take(weights, sort, axis=0)

        # Sum the data on a regular grid.
        bins = np.linspace(0., 1., nbins[i] + 1)
        n_cum = cumsum_to_grid(phase_fold, bins, np.where(weights_fold > 0, 1, 0))
        r_cum = cumsum_to_grid(phase_fold, bins, weights_fold)
        s_cum = cumsum_to_grid(phase_fold, bins, weights_fold*flux_fold)
        q_cum = cumsum_to_grid(phase_fold, bins, weights_fold*flux_fold**2)
        
        # Extend the grid to account for all epochs.
        bins = np.append(bins, bins[1:NumSteps*ES] + bins[-1])
        n_cum = np.append(n_cum, n_cum[1:NumSteps*ES] + n_cum[-1], axis=0)
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES] + r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES] + s_cum[-1], axis=0)
        q_cum = np.append(q_cum, q_cum[1:NumSteps*ES] + q_cum[-1], axis=0)
        
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
        q = q_cum[i2] - q_cum[i1]
        
        # Find the best fit.
        dchisq_tmp = s**2*t/(r*(t - r))
        dchisq_tmp[n < 1] = 0
        
        epoch_tmp = (bins[i1] + bins[i2])/(2*freq[i])
        duration_tmp = (bins[i2] - bins[i1])/freq[i]
        
        # Select the best solution.
        args1 = np.nanargmax(dchisq_tmp, axis=0)
        args2 = np.arange(flux.shape[1])
        
        n = n[args1, args2]
        r = r[args1, args2]
        s = s[args1, args2]
        q = q[args1, args2]
        
        dchisq[i] = dchisq_tmp[args1, args2]
        hchisq[i] = chisq0 - s**2/(t - r) - q
        
        depth[i] = s*t/(r*(t - r))
        epoch[i] = epoch_tmp[args1]
        duration[i] = duration_tmp[args1]
        nt[i] = n
        
    return freq, chisq0, dchisq, hchisq, depth, epoch, duration, nt

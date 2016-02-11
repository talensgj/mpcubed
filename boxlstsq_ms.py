#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m]  

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
        cw = np.concatenate(([zero,], sw.cumsum(axis=0)))
        bin_index = np.r_[sa.searchsorted(bin_edges[:-1], 'left'), sa.searchsorted(bin_edges[-1], 'right')]
        n += cw[bin_index]
    
    return n

def boxlstsq_ms(time, flux, weights, **options):
    
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
    
    M *= Msun # [kg]
    R *= Rsun # [m]
    time = time*SecInDay # [s]
    S = np.ptp(time) # [s]  
    
    # Frequency sampling.
    if freq is not None:
        nfreq = len(freq)
        print 'Using user specified frequencies, ignoring fmin, fmax, OS.'
    
    else:
        
        # Defaults for fmin and fmax.
        fmin_def = 9./S # Three transits in temporal baseline of ground based survey.
        fmax_def = 1./(2.*np.pi)*np.sqrt(G*M/(3.*R)**3.) # Smallest orbit determined by the Hill radius.
        
        if fmin is None:
            fmin = fmin_def
        elif (fmin < fmin_def):
            print 'Warning: minimum frequency contains <2 events.'
            
        if fmax is None:
            fmax = fmax_def
        elif (fmax < fmax_def):
            print 'Warning: maximum frequency is inside the Hill radius.'
                
        if (fmax < fmin):
            raise AttributeError('fmax must be larger than fmin.')
        
        # Compute sample frequencies in range fmin, fmax
        A = ((2.*np.pi)**(2./3.)*R)/(np.pi*(G*M)**(1./3.)*S*OS)
        xmax = int(np.rint((fmax**(1./3.) - fmin**(1./3.))*3./A + 1.))
        xvec = np.array(range(1, xmax))
        freq = (A/3.*(xvec - 1.) + fmin**(1./3.))**3.
        nfreq = len(freq)
    
    print 'Using M=%.2f Msun and R=%.2f Rsun for the stellar parameters.'%(M/Msun, R/Rsun)
    print 'Sampling frequencies between fmin=%.2e, fmax=%.2e per day with OS=%i'%(fmin*SecInDay, fmax*SecInDay, OS)
    print 'Sampling transit epochs at %.2f time(s) the smallest sampled transit.'%(1./ES)
    print 'Sampling transits between %.2f and %.2f time(s) the expected duration in %i steps.'%(1./Qmin, Qmax, NumSteps)
    
    # Calculate the expected transit durations.
    # Calculate the number of bins in the sampling grid.
    q = ((2.*np.pi)**(2./3)*R*freq**(2./3))/(np.pi*(G*M)**(1./3))
    nbins = np.ceil((Qmin*ES)/q)
    #nbins = np.minimum(nbins, len(time))
    
    # Prepare the data for the loop.
    time = time - time[0]
    t = np.nansum(weights, axis=0)
    flux = flux - np.nansum(weights*flux, axis=0)/t # Subtract average
    chisq0 = np.nansum(weights*flux**2., axis=0) # Best fit constant model.
    
    # Create arrays.
    nfreq = len(freq)
    dchisq = np.zeros((nfreq, flux.shape[1]))
    depth = np.zeros((nfreq, flux.shape[1]))
    hchisq = np.zeros((nfreq, flux.shape[1]))
    
    # Loop over frequency domain.
    for i in xrange(nfreq):
        
        # Phase fold the lightcurve. 
        phase = (time*freq[i])%1.
        phase_sort = np.argsort(phase)
        
        phase_fold = np.take(phase, phase_sort, axis=0)
        flux_fold = np.take(flux, phase_sort, axis=0)
        weights_fold = np.take(weights, phase_sort, axis=0)

        # Sum the data on a regular grid.
        bins = np.linspace(0., 1., nbins[i] + 1)
        r_cum = cumsum_to_grid(phase_fold, bins, weights_fold)
        s_cum = cumsum_to_grid(phase_fold, bins, weights_fold*flux_fold)
        q_cum = cumsum_to_grid(phase_fold, bins, weights_fold*flux_fold**2)
        
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES+1] + r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES+1] + s_cum[-1], axis=0)
        q_cum = np.append(q_cum, q_cum[1:NumSteps*ES+1] + q_cum[-1], axis=0)
        
        # Calculate the indices of the start and end of transit for all epochs
        # and durations.
        i1, i2 = np.indices((nbins[i], NumSteps))
        i1 = i1 + 1
        i2 = i1 + ES*(i2 + 1)
       
        i1 = i1.ravel()
        i2 = i2.ravel()
       
        # Compute arrays of r and s values.
        r = r_cum[i2] - r_cum[i1 - 1]
        s = s_cum[i2] - s_cum[i1 - 1]
        q = q_cum[i2] - q_cum[i1 - 1]
        
        # Find the quality of the best fit.
        depth_tmp = s*t/(r*(t - r))
        dchisq_tmp = s*depth_tmp
        hchisq_tmp = chisq0 - s**2/(t - r) - q
        
        args1, args2 = np.where(~np.isfinite(dchisq_tmp))
        dchisq_tmp[args1, args2] = 0
        depth_tmp[args1, args2] = 0
        hchisq_tmp[args1, args2] = chisq0[args2]
            
        args1 = np.nanargmax(dchisq_tmp, axis=0)
        args2 = np.arange(flux.shape[1])
        
        depth[i] = depth_tmp[args1, args2]
        dchisq[i] = dchisq_tmp[args1, args2]
        hchisq[i] = hchisq_tmp[args1, args2]
    
    return freq*SecInDay, dchisq, depth, hchisq, chisq0

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:42:00 2016

@author: talens
"""

import numpy as np

from models import transit
from mpcubed.statistics import statistics

def sigma_clip(array, axis=None, sigma=5., niter=5):
    """ Compute a robust mean and standard deviation."""  

    weights = np.ones(array.shape)
    for i in range(niter):
        
        m0 = np.sum(weights*array, axis=axis)/np.sum(weights, axis=axis)
        m1 = np.sum(weights*(array - m0)**2., axis=axis)/np.sum(weights, axis=axis)
        m1 = np.sqrt(m1)
    
        weights = np.where(np.abs(array - m0) < sigma*m1, 1., 0.)
        
    return m0, m1

def bls_crit(freq, dchisq, epoch, depth, duration):

    # Optimal parameters.
    args_1d = np.argmax(dchisq, axis=0)
    
    if (np.ndim(dchisq) == 2):
        args_2d = np.arange(dchisq.shape[1])
        args_2d = (args_1d, args_2d)
    else:
        args_2d = args_1d
        
    freq0 = freq[args_1d]
    dchisq0 = dchisq[args_2d]
    epoch0 = epoch[args_2d]
    depth0 = depth[args_2d]
    duration0 = duration[args_2d]

    # Signal Detection Efficiency (SDE).    
    mu, sigma = sigma_clip(dchisq, axis=0)
    sde = (dchisq0 - mu)/sigma    
    
    # Anti-transit ratio.
    tmp = dchisq*np.sign(depth)
    atr = -np.amax(tmp, axis=0)/np.amin(tmp, axis=0)

    return freq0, dchisq0, epoch0, depth0, duration0, sde, atr  

def phase_gap(jdmid, freq, epoch, duration):
    
    # Compute and sort the phases.
    phase = freq*(jdmid - epoch)
    phase = np.mod(phase, 1)
    phase = np.sort(phase)
    
    sym = np.sum(phase < .5, dtype='float')/len(phase)    
    
    # Find the largest gap in the coverage.
    gap = np.amax(np.diff(phase))
    gap = np.maximum(gap, 1 - np.ptp(phase))
    gap = gap/(freq*duration)
    
    return gap, sym    
    
def transits(jdmid, freq, epoch, duration):

    # Determine to which orbit and phase each point belongs.
    phase = freq*(jdmid - epoch)
    phase = phase - np.floor(np.amin(phase))
    orbit = np.floor(phase + .5).astype('int')
    phase = np.mod(phase + .5, 1.) - .5
    
    # Select in-transit data points.
    mask = (np.abs(phase) < .5*freq*duration)
    
    # Determine how many points in each transit.
    intransit = np.bincount(orbit[mask])
    intransit = intransit[intransit > 0]

    return len(intransit), intransit

def ellipsoidal_variations(jdmid, mag, emag, freq, epoch, depth, duration):
    
    weights = 1/emag**2
    
    # Evaluate the box-model. 
    boxfit = transit.box(jdmid, 1/freq, epoch, -depth, duration)
    m = np.sum(weights*(mag - boxfit))/np.sum(weights) # Out-of-transit level.
    boxfit = boxfit + m   
    
    # Find in-transit points.
    phase = freq*(jdmid - epoch)
    phase = np.mod(phase + .5, 1.) - .5
    mask = (np.abs(phase) < .5*freq*duration)
    
    # Create the ellipsoidal model.
    p = np.cos(4.*np.pi*phase)
    weights[mask] = 0. # Only fit out-of-transit data.    
    
    # Find the amplitude and S/N of the model.
    u = np.sum(weights*(mag - boxfit)*p)
    v = np.sum(weights*p**2.)
    eps = u/v
    sne = u/np.sqrt(v)     
        
    return eps, sne

def pink_noise(jdmid, mag, freq, epoch, depth, duration):
    
    # The best fit model.
    boxfit = transit.box(jdmid, 1/freq, epoch, -depth, duration) 
    residuals = mag - boxfit    
  
    # The white noise term.
    sigma_white = statistics.mad(residuals)     
    
    # Bin the data to the transit duration.
    nbins = np.ceil(np.ptp(jdmid)/duration)
    bins = duration*np.arange(nbins+1) + np.amin(jdmid)
   
    npbin, edges = np.histogram(jdmid, bins=bins)
    resbin, edges = np.histogram(jdmid, bins=bins, weights=residuals)
    resbin = resbin/npbin
    resbin = resbin[npbin > 0]    
    
    # The red noise term.
    sigma_red = statistics.mad(resbin)    
    
    return sigma_white, sigma_red

def lc_crit(jdmid, mag, emag, freq, epoch, depth, duration):
    
    # Phase coverage.
    gap, sym = phase_gap(jdmid, freq, epoch, duration)
    
    # Number of transits and number of in-transit points.
    ntr, ntp = transits(jdmid, freq, epoch, duration)
    ntp = np.sum(ntp)

    # M statistic.
    mst = np.mean(mag) - np.median(mag)

    # Ellipsoidal variations.
    eps, sne = ellipsoidal_variations(jdmid, mag, emag, freq, epoch, depth, duration)

    # Signal-to-noise.
    sw, sr = pink_noise(jdmid, mag, freq, epoch, depth, duration)
    snp = np.sqrt(depth**2/(sw**2/ntp + sr**2/ntr))

    return gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp
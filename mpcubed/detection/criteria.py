# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:42:00 2016

@author: talens
"""

import numpy as np

from .. import statistics
from ..models import transit

def boxlstsq_criteria(dchisq, depth):

    arg = np.argmax(dchisq)

    with np.errstate(invalid='ignore', divide='ignore'):

        # Signal Detection Efficiency (SDE).    
        mu, sigma = statistics.sigma_clip(dchisq)
        sde = (dchisq[arg] - mu)/sigma    
        
        # Anti-transit ratio.
        tmp = dchisq*np.sign(depth)
        atr = -np.amax(tmp)/np.amin(tmp)

    return sde, atr  

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

    mask = (npbin > 0)    
    resbin = resbin[mask]/npbin[mask]
    
    # The red noise term.
    sigma_red = statistics.mad(resbin)    
    
    return sigma_white, sigma_red

def lightcurve_criteria(jdmid, mag, emag, freq, epoch, depth, duration):
    
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
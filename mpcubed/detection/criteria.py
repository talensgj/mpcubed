# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:42:00 2016

@author: talens
"""

import numpy as np

from .. import statistics
from ..models import transit

def compute_sde(dchisq, width=2500, center=250):

    # Get the location of the peak and convert to signal residue.
    arg = np.argmax(dchisq)
    sr = np.sqrt(dchisq)

    # Select the region for computing the background and noise.
    x = np.abs(np.arange(len(sr)) - arg)
    mask = (x <= width) & (x > center)
    
    if np.all(sr[mask] == 0.): # Some periodograms are pathalogical.
        return 0.
    
    # Compute the signal detection efficiency.
    mu, std = statistics.sigma_clip(sr[mask])
    sde = (sr[arg] - mu)/std
    
    return sde

def boxlstsq_criteria(dchisq, depth):

    with np.errstate(invalid='ignore', divide='ignore'):

        # Signal Detection Efficiency (SDE).    
        sde = compute_sde(dchisq)   
        
        # Anti-transit ratio.
        tmp = dchisq*np.sign(depth)
        atr = -np.amax(tmp)/np.amin(tmp)

    return sde, atr  

def phase_gap(jdmid, boxpars):
    
    # Compute and sort the phases.
    phase = (jdmid - boxpars['epoch'])/boxpars['period']
    phase = np.mod(phase, 1)
    phase = np.sort(phase)
    
    sym = np.sum(phase < .5, dtype='float')/len(phase)    
    
    # Find the largest gap in the coverage.
    gap = np.amax(np.diff(phase))
    gap = np.maximum(gap, 1 - np.ptp(phase))
    gap = gap*boxpars['period']/boxpars['duration']
    
    return gap, sym    
    
def transits(jdmid, boxpars):

    # Determine to which orbit and phase each point belongs.
    phase = (jdmid - boxpars['epoch'])/boxpars['period']
    phase = phase - np.floor(np.amin(phase))
    orbit = np.floor(phase + .5).astype('int')
    phase = np.mod(phase + .5, 1.) - .5
    
    # Select in-transit data points.
    mask = (np.abs(phase) < .5*boxpars['duration']/boxpars['period'])
    
    # Determine how many points in each transit.
    intransit = np.bincount(orbit[mask])
    intransit = intransit[intransit > 0]

    return len(intransit), intransit

def ellipsoidal_variations(jdmid, mag, emag, boxpars):
    
    weights = 1/emag**2
    
    # Evaluate the box-model. 
    boxfit = transit.box(jdmid, boxpars['period'], boxpars['epoch'], -boxpars['depth'], boxpars['duration'])
    m = np.sum(weights*(mag - boxfit))/np.sum(weights) # Out-of-transit level.
    boxfit = boxfit + m   
    
    # Find in-transit points.
    phase = (jdmid - boxpars['epoch'])/boxpars['period']
    phase = np.mod(phase + .5, 1.) - .5
    mask = (np.abs(phase) < .5*boxpars['duration']/boxpars['period'])
    
    # Create the ellipsoidal model.
    p = np.cos(4.*np.pi*phase)
    weights[mask] = 0. # Only fit out-of-transit data.    
    
    # Find the amplitude and S/N of the model.
    u = np.sum(weights*(mag - boxfit)*p)
    v = np.sum(weights*p**2.)
    eps = u/v
    sne = u/np.sqrt(v)     
        
    return eps, sne

def pink_noise(jdmid, mag, boxpars):
    
    # The best fit model.
    boxfit = transit.box(jdmid, boxpars['period'], boxpars['epoch'], -boxpars['depth'], boxpars['duration']) 
    residuals = mag - boxfit    
  
    # The white noise term.
    sigma_white = statistics.mad(residuals)     
    
    # Bin the data to the transit duration.
    nbins = np.ceil(np.ptp(jdmid)/boxpars['duration'])
    bins = boxpars['duration']*np.arange(nbins+1) + np.amin(jdmid)
   
    npbin, edges = np.histogram(jdmid, bins=bins)
    resbin, edges = np.histogram(jdmid, bins=bins, weights=residuals)

    mask = (npbin > 0)    
    resbin = resbin[mask]/npbin[mask]
    
    # The red noise term.
    sigma_red = statistics.mad(resbin)    
    
    return sigma_white, sigma_red

def lightcurve_criteria(jdmid, mag, emag, boxpars):
    
    # Phase coverage.
    gap, sym = phase_gap(jdmid, boxpars)
    
    # Number of transits and number of in-transit points.
    ntr, ntp = transits(jdmid, boxpars)
    ntp = np.sum(ntp)

    # M statistic.
    mst = np.mean(mag) - np.median(mag)

    # Ellipsoidal variations.
    eps, sne = ellipsoidal_variations(jdmid, mag, emag, boxpars)

    # Signal-to-noise.
    sw, sr = pink_noise(jdmid, mag, boxpars)
    snp = np.sqrt(boxpars['depth']**2/(sw**2/ntp + sr**2/ntr))

    return gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp

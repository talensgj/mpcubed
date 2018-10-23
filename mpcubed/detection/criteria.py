# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:42:00 2016

@author: talens
"""

import numpy as np

from .. import statistics

###############################################################################
### Criteria from the box least-squares periodogram.
###############################################################################

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

###############################################################################
### Criteria from the lightcurve.
###############################################################################

def box_model(time, box_pars):
    """A simple box model of a transit."""
    
    if (len(box_pars) != 4):
        raise ValueError('box_pars must have length 4')
    
    T0, P, T14, delta = box_pars
    
    phase = (time - T0)/P
    phase = np.mod(phase+.5, 1.)-.5
    model = np.where(np.abs(phase) < .5*T14/P, delta, 0.)
    
    return phase, model

def phase_gap(jd, box_pars):
    
    # Compute and sort the phases.
    phase = (jd - box_pars['epoch'])/box_pars['period']
    phase = np.mod(phase, 1)
    phase = np.sort(phase)
    
    sym = np.sum(phase < .5, dtype='float')/len(phase)    
    
    # Find the largest gap in the coverage.
    gap = np.amax(np.diff(phase))
    gap = np.maximum(gap, 1 - np.ptp(phase))
    gap = gap*box_pars['period']/box_pars['duration']
    
    return gap, sym    
    
def transits(jd, box_pars):

    # Determine to which orbit and phase each point belongs.
    phase = (jd - box_pars['epoch'])/box_pars['period']
    phase = phase - np.floor(np.amin(phase))
    orbit = np.floor(phase + .5).astype('int')
    phase = np.mod(phase + .5, 1.) - .5
    
    # Select in-transit data points.
    mask = (np.abs(phase) < .5*box_pars['duration']/box_pars['period'])
    
    # Determine how many points in each transit.
    intransit = np.bincount(orbit[mask])
    intransit = intransit[intransit > 0]

    return len(intransit), intransit

def ellipsoidal_variations(jd, mag, emag, box_pars):
    
    weights = 1/emag**2
    
    # Evaluate the box-model. 
    phase, model = box_model(jd, box_pars)
    m = np.sum(weights*(mag - model))/np.sum(weights) # Out-of-transit level.
    model = model + m   
    
    # Find in-transit points.
    phase = (jd - box_pars['epoch'])/box_pars['period']
    phase = np.mod(phase + .5, 1.) - .5
    mask = (np.abs(phase) < .5*box_pars['duration']/box_pars['period'])
    
    # Create the ellipsoidal model.
    p = np.cos(4.*np.pi*phase)
    weights[mask] = 0. # Only fit out-of-transit data.    
    
    # Find the amplitude and S/N of the model.
    u = np.sum(weights*(mag - model)*p)
    v = np.sum(weights*p**2.)
    eps = u/v
    sne = u/np.sqrt(v)     
        
    return eps, sne

def pink_noise(jd, mag, box_pars):
    
    # The best fit model.
    phase, model = box_model(jd, box_pars) 
    residuals = mag - model   
  
    # The white noise term.
    sigma_white = statistics.mad(residuals)     
    
    # Bin the data to the transit duration.
    nbins = np.ceil(np.ptp(jd)/box_pars['duration'])
    bins = box_pars['duration']*np.arange(nbins+1) + np.amin(jd)
   
    npbin, edges = np.histogram(jd, bins=bins)
    resbin, edges = np.histogram(jd, bins=bins, weights=residuals)

    mask = (npbin > 0)    
    resbin = resbin[mask]/npbin[mask]
    
    # The red noise term.
    sigma_red = statistics.mad(resbin)    
    
    return sigma_white, sigma_red

def lightcurve_criteria(jd, mag, emag, box_pars):
    
    # Phase coverage.
    gap, sym = phase_gap(jd, box_pars)
    
    # Number of transits and number of in-transit points.
    ntr, ntp = transits(jd, box_pars)
    ntp = np.sum(ntp)

    # M statistic.
    mst = np.mean(mag) - np.median(mag)

    # Ellipsoidal variations.
    eps, sne = ellipsoidal_variations(jd, mag, emag, box_pars)

    # Signal-to-noise.
    sw, sr = pink_noise(jd, mag, box_pars)
    snp = np.sqrt(box_pars['depth']**2/(sw**2/ntp + sr**2/ntr))

    return gap, sym, ntr, ntp, mst, eps, sne, sw, sr, snp

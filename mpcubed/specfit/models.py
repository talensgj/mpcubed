#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import legendre

import core

def spec_model(pars, wave, spec, espec, wave_mod, spec_mod, deltav, order=3):
    """ Apply Doppler shift and rotational broadening and fit the baseline. """      
    
    # Apply Doppler shift and rotational broadening.
    wave_mod = core.dshift(wave_mod, pars[0])
    spec_mod = core.broaden(spec_mod, deltav, pars[1])
       
    # Interpolate to the data.
    spec_mod = np.interp(wave, wave_mod, spec_mod)
        
    # Fit the baseline using a polynomial.
    x = (wave - np.amin(wave))/np.ptp(wave)
    x = 2.*x - 1.
    y = spec/spec_mod
    yerr = espec/spec_mod    
    
    mat = legendre.legvander(x, order)
    coefs = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    baseline = np.dot(mat, coefs)
    
    return baseline*spec_mod, baseline
    
def spec_chisq(pars, wave, spec, espec, wave_mod, spec_mod, deltav, order=3): 
    """ Evaluate the chi-squared value of the model. """    
    
    # Evaluate the model.
    fit, coefs = spec_model(pars, wave, spec, espec, wave_mod, spec_mod, deltav, order)
    
    # Compute the chi-squared value.
    chisq = (spec - fit)**2/espec**2
    chisq = np.sum(chisq)
    
    return chisq
    
def spec_lnlike(pars, wave_data, spec_data, wave_mod, spec_mod, deltav, order=3):
    """ Evaluate the log-likelihood value of the model. """    
    
    # Evaluate the model.
    fit = spec_model(pars, wave_data, spec_data, wave_mod, spec_mod, deltav, order, method='L-BFGS-B')
    
    # Compute the log-likelihood value.
    lnlike = (spec_data - fit)**2/(spec_data) - np.log(2.*np.pi*spec_data)
    lnlike = -.5*np.sum(lnlike)
    
    return lnlike
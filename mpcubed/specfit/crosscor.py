# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 10:46:49 2016

@author: talens
"""

import numpy as np

import core

def crosscor(x, y):
    """ Compute the normalized cross-correlation of two arrays."""

    x = x - np.mean(x)
    y = y - np.mean(y)

    ccf = np.sum(x*y)/np.sqrt(np.sum(x**2)*np.sum(y**2))    
    
    return ccf    
    
def crosscor_spec(wave, spec, wave_m, spec_m, vcor):
    """ Cross-correlate two spectra at the specified velocities."""

    npoints = len(vcor)
    ccf = np.zeros(npoints)

    for i in range(npoints):
        
        # Shift and interpolate the model to the data.
        wave_tmp = core.dshift(wave_m, vcor[i])
        spec_tmp = np.interp(wave, wave_tmp, spec_m)      
        
        # Compute the cross-correlation.
        ccf[i] = crosscor(spec, spec_tmp)
    
    return ccf
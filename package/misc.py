#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def flux2mag(flux, eflux = None, m0 = 25):
    
    mag = m0 - 2.5*np.log10(flux)
    
    if eflux is not None:
        emag = 2.5/np.log(10.)*np.abs(eflux/flux)
        return mag, emag
    
    return mag

def phase(time, period, time_ref = 0., fold = True):
    
    phase = (time - time_ref)/period
    
    if fold:
        phase = np.mod(phase, 1.)
        
    return phase
    
def mad(data, K = 1.4826, axis = None):
    
    med = np.nanmedian(data, axis = axis, keepdims = True)
    mad = np.nanmedian(np.abs(data - med), axis = axis)
    
    return K*mad

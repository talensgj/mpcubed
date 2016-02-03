#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def flux2mag(flux, eflux=None, m0=25.):
    
    mag = m0 - 2.5*np.log10(flux)
    
    if eflux is not None:
        emag = 2.5/np.log(10.)*np.abs(eflux/flux)
        return mag, emag
    
    return mag

def phase(time, period, time_ref=0., fold=True):
    
    phase = (time - time_ref)/period
    
    if fold:
        phase = np.mod(phase, 1.)
        
    return phase

def zero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

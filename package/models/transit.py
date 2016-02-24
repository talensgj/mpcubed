#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def box_model(t, P, Tp, delta, eta):
    """ Model the transit as a square dip in the lightcurve.
    
    Args:
        t (float): The times at which to evaluatethe transit.
        P (float): The period of the orbit.
        Tp (float): The time of mid-transit.
        delta (float): The fractional depth of the transit.
        eta (float): The transit duration.
        
    Returns:
        model: The model value for each time t.
        
    """
    
    phase = np.mod((t - Tp)/P, 1.)
    phase = np.mod(phase + .5, 1.) - .5
    mask = np.abs(phase) < .5*eta/P
    model = np.where(mask, -delta, 0)
    
    return model

def softened_box_model(t, P, Tp, delta, eta, c):
    """ Model the transit as a softened square dip in the lightcurve.
    
    Args:
        t (float): The times at which to evaluatethe transit.
        P (float): The period of the orbit.
        Tp (float): The time of mid-transit.
        delta (float): The fractional depth of the transit.
        eta (float): The transit duration.
        c (float): Parameter that detrmines the softening. Larger values
            result in more box-like transits.
        
    Returns:
        model: The model value for each time t.
    
    """
    
    x = P*np.sin(np.pi*(t - Tp)/P)/(np.pi*eta)
    model = .5*delta*(np.tanh(c*(x - .5)) - np.tanh(c*(x + .5)))
    
    return model

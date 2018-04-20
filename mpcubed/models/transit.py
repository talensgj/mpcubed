#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def box(t, P, Tp, delta, eta, m=0.):
    """ Model the transit as a square dip in the lightcurve.
    
    Args:
        t (float): The times at which to evaluate the model.
        P (float): The period of the orbit.
        Tp (float): The time of mid-transit.
        delta (float): The fractional depth of the transit.
        eta (float): The transit duration.
        m (float): The out-of-transit level.
        
    Returns:
        model: The model value for each time t.
        
    """
    
    phase = np.mod((t - Tp)/P, 1.)
    phase = np.mod(phase + .5, 1.) - .5
    mask = np.abs(phase) < .5*eta/P
    model = np.where(mask, -delta, 0)
    
    return model + m

def box_chisq(pars, t, y, yerr):
    """ The chi-square value of the box model.
    
    Args:
        pars (float): The parameters of the box model.
        t (float): The times at which to evaluate the model.
        y (float): The data values corresponding to t.
        yerr (float): The measurement errors on y.
        
    """
    
    y0 = box(t, *pars)
    chisq = (y - y0)**2/yerr**2
    
    return np.sum(chisq)

def softbox(t, P, Tp, delta, eta, c=10., m=0.):
    """ Model the transit as a softened square dip in the lightcurve.
    
    Args:
        t (float): The times at which to evaluatethe transit.
        P (float): The period of the orbit.
        Tp (float): The time of mid-transit.
        delta (float): The fractional depth of the transit.
        eta (float): The transit duration.
        c (float): Parameter that detrmines the softening. Larger values
            result in more box-like transits.
        m (float): The out-of-transit level.
        
    Returns:
        model: The model value for each time t.
    
    """
    
    x = P*np.sin(np.pi*(t - Tp)/P)/(np.pi*eta)
    model = .5*delta*(np.tanh(c*(x - .5)) - np.tanh(c*(x + .5)))
    
    return model + m
    
def softbox_chisq(pars, t, y, yerr):
    """ The chi-square value of the softbox model.
    
    Args:
        pars (float): The parameters of the softbox model.
        t (float): The times at which to evaluate the model.
        y (float): The data values corresponding to t.
        yerr (float): The measurement errors on y.
        
    """
    
    y0 = softbox(t, *pars)
    chisq = (y - y0)**2/yerr**2
    chisq = np.sum(chisq)

    if not np.isfinite(chisq):
        chisq = 1e24
    
    return chisq

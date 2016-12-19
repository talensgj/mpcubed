#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import occultquad

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

def projected_seperation(t, t0, P, a, b):    
    """ Compute projected seperation for circular orbits. """        
    
    # Compute the phase and projected separation.
    phase = 2.*np.pi*(t - t0)/P
    z = a*np.sqrt(np.sin(phase)**2. + (b/a)**2.*np.cos(phase)**2.)

    return phase/(2.*np.pi), z

def circ_model(t, Tp, P, tt, b, p, u1, u2):
    """ Mandel & Agol model for circular transits. """    
    
    # Compute the semi-major axis using Seager 2003 eq 8.
    x = np.sin(tt*np.pi/P)**2.
    a = np.sqrt(((1. + p)**2. - (b**2.)*(1. - x))/x) 
    print a
    # Calculate the model lightcurve.
    phase, z = projected_seperation(t, Tp, P, a, b)
    mod = occultquad.occultquad(z, u1, u2, p)[0]
    mod = np.where(np.abs(phase%1. - .5) < .25, 1., mod) # No secondary eclipse.
        
    return phase, mod

def circ_chisq(pars, t, y, yerr, u1, u2, magnitudes=True):
    """ Compute the chi-squared for circular orbits. """    
    
    Tp, P, tt, b, p = pars
    weights = 1./yerr**2.    
    
    # Evaluate the model.
    phase, mod = circ_model(t, Tp, P, tt, b, p, u1, u2)

    # Evaluate the baseline.
    if magnitudes:
        mod = -2.5*np.log10(mod)
        m0 = np.sum(weights*(y - mod))/np.sum(weights)
        chisq = weights*(y - mod - m0)**2.
    else:
        f0 = np.sum(weights*y*mod)/np.sum(weights*mod**2.)
        chisq = weights*(y - f0*mod)
    
    # Compute the chi-squared.
    chisq = np.sum(chisq) 
   
    # Minimization routines do not work with infinities.
    if not np.isfinite(chisq):
        chisq = 1e24
    
    return chisq
    
def circ_lnlike(pars, t, y, yerr, u1, u2, magnitudes=True):
    """ Compute the log-likelihood for circular orbits. """     
    
    Tp, P, tt, b, p = pars
    weights = 1./yerr**2.    
    
    # Evaluate the model.
    phase, mod = circ_model(t, Tp, P, tt, b, p, u1, u2)
    
    # Evaluate the baseline.
    if magnitudes:
        mod = -2.5*np.log10(mod)
        m0 = np.sum(weights*(y - mod))/np.sum(weights)
        lnlike = weights*(y - mod - m0)**2. - np.log(2.*np.pi*yerr**2.)  
    else:
        f0 = np.sum(weights*y*mod)/np.sum(weights*mod**2.)
        lnlike = weights*(y - f0*mod)**2. - np.log(2.*np.pi*yerr**2.)  
    
    # Compute the log-likelihood.
    lnlike = -.5*np.sum(lnlike)    
    
    return lnlike
    
def circ_lnprior(pars):
    
    b = pars[3]
    p = pars[4]
    
    if not (0 < b < (1 + p)):
        return -np.inf
    
    return 0.
    
def circ_lnprob(pars, t, y, yerr, u1, u2, magnitudes=True):
    
    lp = circ_lnprior(pars)
    
    if not np.isfinite(lp):
        return -np.inf
        
    lnlike = circ_lnlike(pars, t, y, yerr, u1, u2, magnitudes)        
        
    return lp + lnlike

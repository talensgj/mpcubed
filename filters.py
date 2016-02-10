#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import fourierfuncs

def harmonic(x, y, weights, P, n, const=False):
    """ Filter long term and lst variations from a lightcurve using sinusoids."""
    
    # Fit n harmonics of the period.
    fmat = fourierfuncs.fourier_mat(x, 1./P, n)
        
    # Allow the offset to vary.
    if const:
        offset = np.ones((len(x),1))
        fmat = np.hstack([fmat, offset])
    
    # Calculate the best fit.
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    chisq = np.sum(weights*(y - fit)**2)
    
    return chisq, pars, fit

def masc_harmonic(jdmid, lst, mag0, weights, Pjd, njd, Plst=24., nlst=5):
    """ Filter long term and lst variations from a lightcurve using sinusoids."""
    
    # Fit with a constant.
    #fmat = np.ones((len(mag0),1))
    
    # Plus harmonics in the baseline.
    if (njd > 0):
        mat_jd = fourierfuncs.fourier_mat(jdmid, 1/Pjd, njd)
        #fmat = np.hstack([fmat, mat_jd])
        fmat = mat_jd
        
    # Plus harmonics in a sidereal day.
    if (nlst > 0):
        mat_lst = fourierfuncs.fourier_mat(lst, 1/Plst, nlst)
        fmat = np.hstack([fmat, mat_lst])
    
    # Calculate the best fit.
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), mag0*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    chisq = np.sum(weights*(mag0 - fit)**2)
    
    return chisq, pars, fit

def hybrid(jdmid, lstseq, mag0, emag0, Pjd, njd):
    
    weights = 1/emag0**2
    idx = lstseq%270
    
    fmat = np.ones((len(mag0),1))
    if (njd > 0):
        mat_jd = fourierfuncs.fourier_mat(jdmid, 1/Pjd, njd)
        fmat = np.hstack([fmat, mat_jd])
    
    fit2 = np.zeros(len(mag0))
    for i in range(10):
        
        pars1 = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), (mag0 - fit2)*np.sqrt(weights))[0]
        fit1 = np.dot(fmat, pars1)
        
        pars2 = np.bincount(idx, weights*(mag0 - fit1))/np.bincount(idx, weights)
        fit2 = pars2[idx]
        
    fit = fit1 + fit2
    chisq = np.sum(weights*(mag0 - fit)**2)/(len(mag0) - len(pars1) - len(np.unique(idx)))
        
    return chisq, fit

def sysrem(mag0, weights):

    par2 = np.ones((mag0.shape[0], 1))
    for i in range(5):
        
        par1 = np.nansum(weights*mag0*par2, axis=1, keepdims=True)/np.nansum(weights*par2**2, axis=1, keepdims=True)
        par2 = np.nansum(weights*mag0*par1, axis=0, keepdims=True)/np.nansum(weights*par1**2, axis=0, keepdims=True)
    
    fit = np.outer(par1, par2)
    
    return fit

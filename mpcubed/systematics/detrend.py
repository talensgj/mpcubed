#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import legendre

from ..models import fourier

def psf_variations(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.)):
    """ Fit long-term and psf-variations using sine and cosine waves.
    
    Args:
        jdmid (float): The julian dates of the observations in UT days.
        lst (float): The sidereal times of the observations in sidereal hours.
        value (float): The observations.
        weights (float): The weights corresponding to the observations.
        ns (int): Tuple containig the sample baseline of the observations in 
                  julian dates and sidereal times.
        step (float): Tuple containing the sample interval in UT days and sidereal hours.
    
    Returns:
        freq1 (float): Frequencies used to fit the long-term variations in UT days.
        freq2 (float): Frequencies used to fit the psf-variations in sidereal hours.
        pars1 (float): The best-fit amplitudes of the long-term variations.
        pars2 (float): The best-fit amplitudes of the psf-variations.
        fit1 (float): The best-fitting long-term variations.
        fit2 (float): The best fitting psf-variations.
        chisq (float): The best fit chi-square value.
        
    """
    
    # Compute the frequencies and matrix for the JD times.
    freq1 = fourier.fftfreq(step[0], ns[0])
    freq1 = np.append(np.amin(freq1)/2., freq1)
    freq1 = freq1[freq1 < 1/6.]
    
    if (len(freq1) != 0):
        mat1 = fourier.fourier_mat(jdmid, freq1)
    else:
        mat1 = np.array([[]]*len(jdmid))
        
    # Compute the frequencies and matrix for the LST times.
    freq2 = fourier.fftfreq(step[1], ns[1])
    freq2 = np.append(np.amin(freq2)/2., freq2)
    freq2 = freq2[freq2 < 1/.5]
    
    if (len(freq2) != 0):
        mat2 = fourier.fourier_mat(lst, freq2)
    else:
        mat2 = np.array([[]]*len(jdmid))
        
    # Compute the full matrix of basis functions.
    mat = np.hstack([mat1, mat2])
    
    if (mat.shape[1] == 0):
        return [], [], [], [], np.zeros(len(jdmid)), np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    # Calculate the best fit parameters.
    pars = fourier.fit_mat(value, weights, mat)
    
    # Evaluate the fit.
    n = 2*len(freq1)
    fit1 = np.sum(mat[:,:n]*pars[:n], axis=1)
    fit2 = np.sum(mat[:,n:]*pars[n:], axis=1)
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit1 - fit2)**2
    chisq = np.sum(chisq)

    return freq1, freq2, pars[:n], pars[n:], fit1, fit2, chisq

def wrap_lst(lst):
    
    sort = np.argsort(lst)
    gap = np.amax(np.diff(lst[sort])) 
    arg = np.argmax(np.diff(lst[sort]))
    gap0 = (np.amin(lst) + 24.) - np.amax(lst) # Gap across lst=0.
    
    if (gap > gap0):
        
        lst = np.mod(lst - lst[sort][arg+1], 24.)
    
    return lst
    
def scale(x):
    
    xmax = np.amax(x)
    xmin = np.amin(x)    
    
    x = (x - xmin)/(xmax - xmin)
    x = 2.*x - 1.    
    
    return x

def detrend_pol(jd, lst, mag, emag, scale0=3., scale1=.25):
    
    # Offset.
    mat = [np.ones(len(mag))]    
        
    # Long-term variations.
    deg0 = int(np.floor(np.ptp(jd)/scale0))
    x = scale(jd)
    mat0 = legendre.legvander(x, deg0)[:,1:]
    mat.append(mat0)        
        
    # PSF variations.
    lst = wrap_lst(lst)
    deg1 = int(np.floor(np.ptp(lst)/scale1))
    x = scale(lst)
    mat1 = legendre.legvander(x, deg1)[:,1:]
    mat.append(mat1)
    
    # Solve.
    mat = np.column_stack(mat)
    pars = np.linalg.lstsq(mat/emag[:,None], mag/emag)[0]

    # Evaluate.
    fit0 = np.sum(pars[:deg0+1]*mat[:,:deg0+1], axis=1)
    fit1 = np.sum(pars[deg0+1:]*mat[:,deg0+1:], axis=1)
    
    return fit0, fit1

#def masc_harmonic(jdmid, lst, value, weights, Pjd, njd, Plst, nlst, cnst=False):
    #""" Compute the best fit model for the long-term and LST trends. """
    
    ## Compute the frequencies and matrix for the JD times.
    #freq_jd = fourier.frequencies(Pjd, njd, False)
    #mat_jd = fourier.fourier_mat(jdmid, freq_jd)
    
    ## Compute the frequencies and matrix for the LST times.
    #freq_lst = fourier.frequencies(Plst, nlst, False)
    #mat_lst = fourier.fourier_mat(lst, freq_lst)
    
    ## Compute the full matrix of basis functions.
    #if cnst: 
        #ones = np.ones((len(value),1))
        #mat = np.hstack([ones, mat_jd, mat_lst])
    #else:
        #mat = np.hstack([mat_jd, mat_lst])
    
    ## Calculate the best fit.
    #pars = fourier.fit_mat(value, weights, mat)
    #fit = np.dot(mat, pars)
    
    ## Calculate the chi-square value of the fit.
    #chisq = weights*(value - fit)**2
    #chisq = np.sum(chisq)
    
    #return pars, fit, chisq
    
    
#def hybrid(jdmid, lstidx, mag, weights, step, ns, maxiter=10):
    
    #freqs = fourier.fftfreq(step, ns)
    #freqs = freqs[freqs < 1/6.]
    #freqs = np.append(np.amin(freqs)/2, freqs)
    #mat = fourier.fourier_mat(jdmid, freqs)
    
    ##x = jdmid[:,None]/np.amax(jdmid)
    ##mat = np.hstack([x, x**2, mat])
    
    #fit2 = 0
    #for niter in range(maxiter):
    
        #pars1 = fourier.fit_mat(mag - fit2, weights, mat)
        #fit1 = np.dot(mat, pars1)
        
        #x0 = np.bincount(lstidx, weights)
        #x1 = np.bincount(lstidx, weights*(mag - fit1))

        #pars2 = x1/x0
        #fit2 = pars2[lstidx]
        
    #chisq = weights*(mag - fit1 - fit2)**2.
    #chisq = np.sum(chisq)
    
    #pars = np.append(pars1, pars2[np.unique(lstidx)]) 
    
    #return pars, fit1 + fit2, chisq

#def new_harmonic(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.), cnst=False):
    
    ## Compute the frequencies and matrix for the JD times.
    #freq1 = fourier.fftfreq(step[0], ns[0], False)
    #freq1 = np.append(np.amin(freq1)/2, freq1)
    #freq1 = freq1[freq1 < 1/6.]
    #if (len(freq1) != 0):
        #mat1 = fourier.fourier_mat(jdmid, freq1)
    #else:
        #mat1 = np.array([[]]*len(jdmid))
        
    ## Compute the frequencies and matrix for the LST times.
    #freq2 = fourier.fftfreq(step[1], ns[1], False)
    #freq2 = np.append(np.amin(freq2)/2, freq2)
    #freq2 = freq2[freq2 < 1/.5]
    #if (len(freq2) != 0):
        #mat2 = fourier.fourier_mat(lst, freq2)
    #else:
        #mat2 = np.array([[]]*len(jdmid))
        
    ## Compute the full matrix of basis functions.
    #if cnst: 
        #ones = np.ones((len(value),1))
        #mat = np.hstack([ones, mat1, mat2])
    #else:
        #mat = np.hstack([mat1, mat2])
    
    #if (mat.shape[1] == 0):
        #return [], np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    ## Calculate the best fit.
    #pars = fourier.fit_mat(value, weights, mat)
    #fit = np.dot(mat, pars)
    
    ## Calculate the chi-square value of the fit.
    #chisq = weights*(value - fit)**2
    #chisq = np.sum(chisq)
    
    #return pars, fit, chisq

#def new_harmonic2(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.), cnst=False):
    
    ## Compute the frequencies and matrix for the JD times.
    #freq1 = fourier.fftfreq(step[0], ns[0], False)
    #freq1 = np.append(np.amin(freq1)/2, freq1)
    #freq1 = freq1[freq1 < 1/6.]
    #if (len(freq1) != 0):
        #mat1 = fourier.fourier_mat(jdmid, freq1)
    #else:
        #mat1 = np.array([[]]*len(jdmid))
        
    ## Compute the frequencies and matrix for the LST times.
    #freq2 = fourier.fftfreq(step[1], ns[1], False)
    #freq2 = np.append(np.amin(freq2)/2, freq2)
    #freq2 = freq2[freq2 < 1/.5]
    #if (len(freq2) != 0):
        #mat2 = fourier.fourier_mat(lst, freq2)
    #else:
        #mat2 = np.array([[]]*len(jdmid))
        
    ## Compute the full matrix of basis functions.
    #if cnst: 
        #ones = np.ones((len(value),1))
        #mat = np.hstack([ones, mat1, mat2])
    #else:
        #mat = np.hstack([mat1, mat2])
    
    #if (mat.shape[1] == 0):
        #return [], np.zeros(len(jdmid)), np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    ## Calculate the best fit.
    #pars = fourier.fit_mat(value, weights, mat)
    #fit = np.dot(mat, pars)
    
    #n = 2*len(freq1)
    #fit1 = np.dot(mat[:,:n], pars[:n])
    #fit2 = np.dot(mat[:,n:], pars[n:])
    
    ## Calculate the chi-square value of the fit.
    #chisq = weights*(value - fit)**2
    #chisq = np.sum(chisq)
    
    #return pars, fit1, fit2, chisq

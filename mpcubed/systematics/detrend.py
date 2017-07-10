#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import legendre

from ..models import fourier

def detrend_fourier(jdmid, lst, mag, emag, ns, wrap, step=(.003693591 ,320./3600.)):
    """ Fit long-term and psf-variations using sine and cosine waves.
    
    Args:
        jd (float): The julian dates of the observations in UT days.
        lst (float): The sidereal times of the observations in sidereal hours.
        mag (float): The data.
        emag (float): The uncertainties on the data.
        ns (int): Tuple containig the sample baseline of the observations in 
                  julian dates and sidereal times.
        wrap (bool): If True wrap the sidereal times by 12 hours.
        step (float): Tuple containing the sample interval in UT days and sidereal hours.
    
    Returns:
        mat (array): The matrix of babsis functions used in the fit.
        fit1 (float): The best-fitting long-term variations.
        fit2 (float): The best fitting psf-variations.
        
    """
    
    if wrap:
        lst = np.mod(lst+12., 24.)-12.
    
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
        return mat, np.zeros(len(jdmid)), np.zeros(len(jdmid))
    
    # Calculate the best fit parameters.
    pars = fourier.fit_mat(mag, emag, mat)
    
    # Evaluate the fit.
    n = 2*len(freq1)
    fit1 = np.sum(mat[:,:n]*pars[:n], axis=1)
    fit2 = np.sum(mat[:,n:]*pars[n:], axis=1)

    return mat, fit1, fit2

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

def detrend_legendre(jd, lst, mag, emag, scale0=3., scale1=.25):
    """ Fit long-term and psf-variations using legendre polynomials.
    
    Args:
        jd (float): The julian dates of the observations in UT days.
        lst (float): The sidereal times of the observations in sidereal hours.
        mag (float): The data.
        emag (float): The uncertainties on the data.
        scale0 (float): Smallest scale of the long-term variations on days.
        scale1 (float): Smallest scale of the psf-variations on days.
    
    Returns:
        mat (array): The matrix of babsis functions used in the fit.
        fit1 (float): The best-fitting long-term variations.
        fit2 (float): The best fitting psf-variations.
        
    """    
    
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
    
    return mat, fit0, fit1

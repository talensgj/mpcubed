#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import legendre

from .. import statistics
from ..models import fourier

###############################################################################
### Helper functions.
###############################################################################

def wrap_lst(lst):
    """ Wrap the LST. """
    
    sort = np.argsort(lst)
    gap = np.amax(np.diff(lst[sort])) 
    arg = np.argmax(np.diff(lst[sort]))
    gap0 = (np.amin(lst) + 24.) - np.amax(lst) # Gap across lst=0.
    
    if (gap > gap0):
        
        lst = np.mod(lst - lst[sort][arg+1], 24.)
    
    return lst
    
def scale(x):
    """ Scale a coordinate to run from -1 to 1. """
    
    xmax = np.amax(x)
    xmin = np.amin(x)    
    
    x = (x - xmin)/(xmax - xmin)
    x = 2.*x - 1.    
    
    return x

def moving_mean(x, y, yerr=None, window=3.):
    """ Compute a moving mean along the x-axis. """

    # Set the weights.
    if yerr is None:
        weights = np.ones_like(y)
    else:
        weights = 1/yerr**2

    # Sums for computing the mean.
    sum1 = np.append(0, np.cumsum(weights*y))
    sum2 = np.append(0, np.cumsum(weights))

    # Indices at the start and end of the window.
    i = np.searchsorted(x, x - window/2.)
    j = np.searchsorted(x, x + window/2.)

    # Compute the mean.
    mean = (sum1[j] - sum1[i])/(sum2[j] - sum2[i])
    
    return mean

###############################################################################
### Detrending methods.
###############################################################################

def detrend_legendre(jd, lst, sky, mag, emag, scale0=3., scale1=.25, sig=3., maxiter=1):
    """ Fit long-term and psf-variations using legendre polynomials.
    
    Args:
        jd (float): The julian dates of the observations in UT days.
        lst (float): The sidereal times of the observations in sidereal hours.
        sky (float): The sky values associated with the data.
        mag (float): The data.
        emag (float): The uncertainties on the data.
        scale0 (float): Smallest scale of the long-term variations on days.
        scale1 (float): Smallest scale of the psf-variations on days.
    
    Returns:
        mat (array): The matrix of babsis functions used in the fit.
        fit1 (float): The best-fitting long-term variations.
        fit2 (float): The best fitting psf-variations.
        
    """       
    
    mat = []    
        
    # Long-term variations.
    deg0 = int(np.floor(np.ptp(jd)/scale0))
    x = scale(jd)
    mat0 = legendre.legvander(x, deg0)[:,1:]
    mat.append(mat0)        
        
    # PSF variations.
    lst = wrap_lst(lst)
    deg1 = int(np.floor(np.ptp(lst)/scale1))
    x = scale(lst)
    mat1 = legendre.legvander(x, deg1)
    mat.append(mat1)

    # Sky dependency.    
    mat.append(mat1*sky[:,None]/np.amax(sky))    
    
    # Solve.
    mat = np.column_stack(mat)
    mask = np.ones_like(jd, dtype='bool')
    
    for niter in range(maxiter):
        
        # Compute the best fit.
        pars = np.linalg.lstsq(mat[mask]/emag[mask,None], mag[mask]/emag[mask])[0]
        fit = np.sum(pars*mat, axis=1)
        
        # Compute the residuals.
        res = mag - fit
        m0 = np.nanmedian(res)
        m1 = statistics.mad(res)

        # Update the mask.
        mask_old = np.copy(mask)
        mask = np.abs(res - m0) < sig*m1         
        
        # End the iterative process.
        if np.all(mask == mask_old):
            break
        
    # Evaluate the best fit.
    fit0 = np.sum(pars[:deg0]*mat[:,:deg0], axis=1)
    fit1 = np.sum(pars[deg0:deg0+deg1+1]*mat[:,deg0:deg0+deg1+1], axis=1)
    fit2 = np.sum(pars[deg0+deg1+1:]*mat[:,deg0+deg1+1:], axis=1)

    return mat, fit0, fit1, fit2

def psfsky(lstidx, sky, mag, emag):
    
    _, idx = np.unique(lstidx, return_inverse=True)
    
    nobs = np.bincount(idx)
    
    xbar = np.bincount(idx, sky/emag**2)/np.bincount(idx, 1/emag**2)
    ybar = np.bincount(idx, mag/emag**2)/np.bincount(idx, 1/emag**2)
    
    b = np.bincount(idx, (sky - xbar[idx])*(mag - ybar[idx])/emag**2)/np.bincount(idx, (sky - xbar[idx])**2/emag**2)
    b = np.where(nobs > 1, b, 0.)
    a = ybar - b*xbar

    fit1 = a[idx]
    fit2 = b[idx]*sky  
    
    return fit1, fit2   
    
def detrend_snellen(jd, lstseq, sky, mag, emag, window=3., maxiter=50, dtol=1e-3):
    
    lstidx = (lstseq % 270)
    fit0 = np.zeros(len(jd))
    fit1 = np.zeros(len(jd)) 
    
    fit = np.zeros_like(mag)
    for niter in range(maxiter):
            
        fit1, fit2 = psfsky(lstidx, sky, mag - fit0, emag)
        fit0 = moving_mean(jd, mag - fit1 - fit2, emag, window)
        
        if niter > 0:
            
            if np.all(np.abs(fit - fit0 - fit1 - fit2) < dtol):
                break
            
        fit = fit0 + fit1 + fit2
        
    return fit0, fit1, fit2

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

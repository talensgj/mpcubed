#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from .. import misc, statistics

###############################################################################
### Helper functions.
###############################################################################

def wrap_lst(lst):
    """ Wrap the LST. """
    
    sort = np.argsort(lst)
    gap = np.amax(np.diff(lst[sort])) 
    arg = np.argmax(np.diff(lst[sort]))
    gap0 = 24. - np.ptp(lst) # Gap across lst=0.
    
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

###############################################################################
### Detrending methods.
###############################################################################

def polynomial(jd, mag, emag, deg=0):
    
    # Create the matrix.
    x = scale(jd)
    mat = np.polynomial.legendre.legvander(x, deg)
    
    # Compute the best fit.
    pars = np.linalg.lstsq(mat/emag[:,np.newaxis], mag/emag, rcond=None)[0]
    
    # Evaluate the best fit.
    trend = np.sum(pars*mat, axis=1)
    chisq = np.sum(((mag - trend)/emag)**2)
    
    return trend, mat, pars, chisq

def legendre(jd, lst, sky, mag, emag, s_jd=5.0, s_lst=0.25, sig=10., maxiter=50):
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
    deg0 = int(np.floor(np.ptp(jd)/s_jd))
    x = scale(jd)
    mat0 = np.polynomial.legendre.legvander(x, deg0)[:,1:]
    mat.append(mat0)        
        
    # PSF variations.
    lst = wrap_lst(lst)
    deg1 = int(np.floor(np.ptp(lst)/s_lst))
    x = scale(lst)
    mat1 = np.polynomial.legendre.legvander(x, deg1)
    mat.append(mat1)

    # Sky dependency.    
    mat.append(mat1*sky[:,np.newaxis]/np.amax(sky))    
    
    # Solve.
    mat = np.column_stack(mat)
    mask = np.ones_like(jd, dtype='bool')
    
    for niter in range(maxiter):

        # Compute the best fit.
        pars = np.linalg.lstsq(mat[mask]/emag[mask,np.newaxis], mag[mask]/emag[mask], rcond=None)[0]
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
            
        if np.all(~mask):
            break
        
    # Evaluate the best fit.
    trend = np.sum(pars*mat, axis=1)
    chisq = np.sum(((mag - trend)/emag)**2)

    return trend, mat, pars, chisq

def linfit(lstidx, x, y, sky, mag, emag):
    
    sort = np.argsort(lstidx)
    invsort = np.argsort(sort)
    
    lstidx = lstidx[sort]
    x = x[sort]
    y = y[sort]
    sky = sky[sort]
    mag = mag[sort]
    emag = emag[sort]
    
    _, idx = np.unique(lstidx, return_inverse=True)
    
    nobs = np.bincount(idx)
    strides = np.append(0, np.cumsum(nobs))
    
    xbar = np.bincount(idx, x)/np.bincount(idx)
    ybar = np.bincount(idx, y)/np.bincount(idx)
    
    mat = np.column_stack([np.ones_like(mag), x-xbar[idx], y-ybar[idx], sky])
    
    pars = np.zeros((len(nobs), 4))
    pars[:,0] = np.bincount(idx, mag/emag**2)/np.bincount(idx, 1/emag**2)
    
    for i in range(len(nobs)):
        
        if nobs[i] < 5:
            continue
            
        i1 = strides[i]
        i2 = strides[i+1]
        
        pars[i] = np.linalg.lstsq(mat[i1:i2]/emag[i1:i2,np.newaxis], mag[i1:i2]/emag[i1:i2], rcond=None)[0]

    trend = np.sum(pars[idx]*mat, axis=1)

    return trend[invsort], (nobs > 4)[idx][invsort] 

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

def local_linear(jd, lstseq, x, y, sky, mag, emag, window=5., maxiter=50, dtol=1e-3):
    
    lstidx = (lstseq % 270)
    
    trend = np.zeros_like(mag)
    trend0 = np.zeros_like(mag)
    trend1 = np.zeros_like(mag)
    
    for niter in range(maxiter):
            
        trend0, mask = linfit(lstidx, x, y, sky, mag - trend1, emag)
        trend1 = moving_mean(jd, mag - trend0, emag, window)
        
        if niter > 0:
            
            if np.all(np.abs(trend - trend0 - trend1) < dtol):
                break
            
        trend = trend0 + trend1
            
    # Evaluate the best fit.  
    trend = trend0 + trend1
    chisq = np.sum(((mag - trend)/emag)**2)
        
    return trend, mask, chisq 

def fourier_matrix(time, freq):
    
    tmp = 2.*np.pi*np.outer(time, freq)
    mat = np.column_stack([np.sin(tmp), np.cos(tmp)])
    
    return mat

def fourier(jd, lst, mag, emag, ns, wrap, step=(.003693591 ,320./3600.)):
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
    freq1 = np.fft.rfftfreq(ns[0], step[0])[1:]
    freq1 = np.append(np.amin(freq1)/2., freq1)
    freq1 = freq1[freq1 < 1/6.]
    
    if (len(freq1) != 0):
        mat1 = fourier_matrix(jd, freq1)
    else:
        mat1 = np.empty((len(jd), 0))
        
    # Compute the frequencies and matrix for the LST times.
    freq2 = np.fft.rfftfreq(ns[1], step[1])[1:]
    freq2 = np.append(np.amin(freq2)/2., freq2)
    freq2 = freq2[freq2 < 1/.5]
    
    if (len(freq2) != 0):
        mat2 = fourier_matrix(lst, freq2)
    else:
        mat2 = np.empty((len(jd), 0))
        
    # Compute the full matrix of basis functions.
    mat = np.column_stack([mat1, mat2])
    
    if (mat.shape[1] == 0):
        return np.zeros_like(jd), mat, [], np.sum((mag/emag)**2)
    
    # Calculate the best fit parameters.
    pars = np.linalg.lstsq(mat/emag[:,np.newaxis], mag/emag, rcond=None)[0]
    
    # Evaluate the fit.
    trend = np.sum(mat*pars, axis=1)
    chisq = np.sum(((mag - trend)/emag)**2)

    return trend, mat, pars, chisq

###############################################################################
### Remove trend from light curve.
###############################################################################
    
def remove_trend(lc, nobs, model=None, method='legendre', options={}):

    strides = np.append(0, np.cumsum(nobs))

    # Loop over blocks of data.
    for i in range(nobs.size):
        
        i1 = strides[i]
        i2 = strides[i+1]
        
        # Check if there is data in this block.
        mask = lc['mask'][i1:i2]

        if np.all(~mask):
            continue
        
        # Extract the good data in the block.
        tmp = lc[i1:i2][mask]
        
        # Subtract the model.
        if model is None:
            mag = tmp['mag']
        else:
            mag = tmp['mag'] - model[i1:i2][mask]
        
        # Detrend the lightcurves.
        if method == 'none':
            pass
        
        elif method == 'polynomial':
            
            trend, mat, pars, chisq = polynomial(tmp['jd'], mag, tmp['emag'], **options)
        
            lc['trend'][i1:i2][mask] = trend
        
        elif method == 'legendre':        
        
            trend, mat, pars, chisq = legendre(tmp['jd'], tmp['lst'], tmp['sky'], mag, tmp['emag'], **options)

            lc['trend'][i1:i2][mask] = trend

        elif method == 'loclin':
            
            trend, flag, chisq = local_linear(tmp['jd'], tmp['lstseq'], tmp['x'], tmp['y'], tmp['sky'], mag, tmp['emag'], **options)
            
            lc['trend'][i1:i2][mask] = trend
            lc['mask'][i1:i2][mask] = flag
              
        elif method == 'fourier':        
    
            ns = [0,0]
            ns[0] = np.ptp(tmp['lstseq']) + 1
            ns[1], wrap = misc.find_ns(tmp['lstseq'])
            ns = np.maximum(ns, 2)
            
            trend, mat, pars, chisq = fourier(tmp['jd'], tmp['lst'], mag, tmp['emag'], ns, wrap, **options)
            
            lc['trend'][i1:i2][mask] = trend
            
        else:
            raise ValueError('Unknown detrending method "{}"'.format(method)) 
    
    return lc

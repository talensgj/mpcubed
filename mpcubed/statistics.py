#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def mad(data, K = 1.4826, axis = None):
    """ Compute the median absolute deviation of an array.
    
    Args:
        array: An array of data.
        K (float): A scale factot so that the mad of a gaussian distribution
            is equal to it standard deviation.
        axis (int): The axis along whic to compute the mad.
        
    Returns:
        The mad, scaled by K.
    
    """
    
    med = np.nanmedian(data, axis = axis, keepdims = True)
    mad = np.nanmedian(np.abs(data - med), axis = axis)
    
    return K*mad
    
def sigma_clip(array, axis=None, ndev=5., niter=5):
    """ Compute a robust mean and standard deviation."""  

    weights = np.ones(array.shape)
    for i in range(niter):
        
        mu = np.sum(weights*array, axis=axis)/np.sum(weights, axis=axis)
        sigma = np.sqrt(np.sum(weights*(array - mu)**2., axis=axis)/np.sum(weights, axis=axis))
    
        weights = np.where(np.abs(array - mu) < ndev*sigma, 1., 0.)
        
    return mu, sigma

def bin_data(x, y, bins, yerr=None):
    
    # Set the weights.
    if yerr is None:
        weights = np.ones_like(y)
    else:
        weights = 1./yerr**2

    # Compute quantities in the bins.
    x_bin = (bins[:-1] + bins[1:])/2.
    w_sum, bins = np.histogram(x, weights=weights, bins=bins)
    wy_sum, bins = np.histogram(x, weights=weights*y, bins=bins)
    wysq_sum, bins = np.histogram(x, weights=weights*y**2, bins=bins)
    
    # Remove empty bins.
    mask = (w_sum > 0)
    x_bin = x_bin[mask]
    w_sum = w_sum[mask]
    wy_sum = wy_sum[mask]
    wysq_sum = wysq_sum[mask]
    
    # Compute results.
    y_bin = wy_sum/w_sum

    if yerr is None:
        yerr_bin = np.sqrt(wysq_sum/w_sum - (wy_sum/w_sum)**2)/np.sqrt(w_sum)
    else:
        yerr_bin = 1./np.sqrt(w_sum)
        
    return x_bin, y_bin, yerr_bin

def idxstats(indices, values, statistic='mean', keeplength=False):
    """ Compute a statistic for all values with the same index.
    
    Args:
        indices (int): An array of indices.
        values (float): An array of values.
        statistic (string or function): The statistic to compute on values
            that have the same index may be any of 'mean', 'std', 'count',
            'sum', 'median' or a function. Default is 'mean'.
        keeplength (bool): If True the return will have the same shape as
            values, otherwise it will be the same length as the number of 
            unique indices. Default is False.
            
    Returns:
        result: The statistic computed on values with the same index.
    
    """
    
    # Check that the statistic is valid.
    known_stats = ['mean', 'std', 'count', 'sum', 'median']
    if not callable(statistic) and statistic not in known_stats:
        raise ValueError('invalid statistic {}'.format(statistic))
    
    # Make sure the indices are integers.
    indices = indices.astype('int') # It would be better to return an error if they are not integer, but how...
    
    # If only one scalar index is given bincount will not work.
    if np.isscalar(indices):
        
        if statistic == 'std':
            return 0.
        if statistic == 'count':
            return 1
        elif statistic in ['mean', 'sum', 'median']:
            return values
        else:
            return statistic(values)
    
    # Count the number of datapoints at each index.
    flatcount = np.bincount(indices)
    
    # Obtain the unique indices
    a = flatcount.nonzero()
    
    # Compute the desired statistic.
    if statistic == 'mean':
        flatsum = np.bincount(indices, values)
        result = flatsum / flatcount
        
    elif statistic == 'std':
        flatsum = np.bincount(indices, values)
        flatsum2 = np.bincount(indices, values ** 2)
        result = np.sqrt(flatsum2 / flatcount - (flatsum / flatcount) ** 2)

    elif statistic == 'count':
        result = flatcount
        
    elif statistic == 'sum':
        result = np.bincount(indices, values)
        
    elif statistic == 'median':
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = np.median(values[indices == i])
        
    else:
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = statistic(values[indices == i])
    
    # Return the statistic in the desired format.
    if not keeplength: 
        return result[a]
    else:
        return result[indices]

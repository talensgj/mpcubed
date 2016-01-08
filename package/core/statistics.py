#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def idxstats(indices, values, statistic='mean', keeplength=False):
    """Compute some statistic on the values that have the same index."""
    
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

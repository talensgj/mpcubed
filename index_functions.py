#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def index_statistics(indices, values, statistic='mean', keeplength=False):
    
    indices = indices.astype('int')
    
    known_stats = ['mean', 'std', 'count', 'sum', 'median']
    if not callable(statistic) and statistic not in known_stats:
        raise ValueError('invalid statistic %r' % (statistic,))
    
    flatcount = np.bincount(indices, None)
    a = flatcount.nonzero()
    
    if statistic == 'mean':
        flatsum = np.bincount(indices, values)
        result = flatsum / flatcount
        
    elif statistic == 'std':
        
        flatsum = np.bincount(indices, values)
        flatsum2 = np.bincount(indices, values ** 2)
        result = np.sqrt(flatsum2[a] / flatcount[a] - (flatsum[a] / flatcount[a]) ** 2)

    elif statistic == 'count':
        result = flatcount[a]
        
    elif statistic == 'sum':
        flatsum = np.bincount(indices, values)
        result = flatsum[a]
        
    elif statistic == 'median':
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = np.median(values[indices == i])
        result = result[a]
        
    else:
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = statistic(values[indices == i])
        result = result[a]
    
    if not keeplength:    
        return result[a]
    else:
        return result[indices]


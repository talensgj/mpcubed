# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 10:50:59 2016

@author: talens
"""

import numpy as np

def sigma_clip(x, sigma=5., niter=5):
    """ Compute a robust mean and median."""    
    
    mask = np.ones(len(x), dtype='bool')
    for i in range(niter):
        
        m0 = np.nanmean(x[mask])
        m1 = np.nanstd(x[mask])
        
        mask = (np.abs(x - m0) < sigma*m1)
        
    return m0, m1
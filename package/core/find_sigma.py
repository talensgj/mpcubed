#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
    
def find_sigma(ind1, res, err, maxiter=10, eps=1e-3):
    
    err1 = np.zeros(np.amax(ind1)+1)
    err2 = 2*np.ones(np.amax(ind1)+1)
    
    wgt = 1./(err*err + (err1[ind1])**2)
    term2 = np.bincount(ind1, res*res*wgt*wgt)
    term1 = np.bincount(ind1, wgt)
    diff1 = term2 - term1
    
    args1, = np.where(diff1 < 1e-10)
    
    wgt = 1./(err*err + (err2[ind1])**2)
    term2 = np.bincount(ind1, res*res*wgt*wgt)
    term1 = np.bincount(ind1, wgt)
    diff2 = term2 - term1
    
    args2, = np.where(diff2 > 1e-10)
    
    for niter in range(maxiter):
        err3 = (err2 + err1)/2.
        
        wgt = 1./(err*err + (err3[ind1])**2)
        term2 = np.bincount(ind1, res*res*wgt*wgt)
        term1 = np.bincount(ind1, wgt)
        diff3 = term2 - term1
        
        here = (diff3 > 1e-10 )
        err1[here] = err3[here]
        diff1[here] = diff3[here]
        err2[~here] = err3[~here]
        diff2[~here] = diff3[~here]
        
        if np.all((err2-err1) < eps):
            break
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
            
    return err3

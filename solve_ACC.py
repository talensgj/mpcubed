#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def find_err(res, err, maxiter=50):

    err1 = 0.
    err2 = 2.
        
    wgt = 1./(err*err + err1*err1)
    term2 = np.sum(res*res*wgt*wgt)
    term1 = np.sum(wgt)
    diff1 = term2 - term1
    
    if (diff1 < 1e-10):
        return 0.
        
    wgt = 1./(err*err + err2*err2)
    term2 = np.sum(res*res*wgt*wgt)
    term1 = np.sum(wgt)
    diff2 = term2 - term1
    
    if (diff2 > 1e-10):
        return 2.
        
    for niter in range(maxiter):
        err3 = ( err2 + err1 ) / 2.
        if ( err2-err1 < 3e-4 ):
            return err3
          
        wgt = 1./(err*err + err3*err3)
        term2 = np.sum(res*res*wgt*wgt)
        term1 = np.sum(wgt)
        diff3 = term2 - term1
          
        if (diff3 > 1e-10 ):
            err1 = err3
            diff1 = diff3
        else:
            err2 = err3
            diff2 = diff3

    return 0.

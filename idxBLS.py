#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def idxBLS(time, flux, flux_err):
    
    S = np.max(time)-np.min(time) # [s]  
    weights = 1./flux_err**2.
    
    t = np.sum(weights, axis=0)
    flux = flux - np.sum(flux*weights, axis=0)/t # Subtract average
    chisq0 = np.sum(flux**2.*weights, axis=0) # B

    print S/2

    freq = np.arange(50, np.floor(S/2)).astype('int')
    
    NumSteps = 9
    ES = 3
    
    dchisq = np.zeros(len(freq))
    for i in xrange(len(freq)):
        
        subidx = time%freq[i]
        
        r = np.bincount(subidx, weights, minlength=freq[i])
        s = np.bincount(subidx, flux*weights, minlength=freq[i])
        
        r_cum = np.cumsum(r)
        s_cum = np.cumsum(s)
        
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES+1]+r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES+1]+s_cum[-1], axis=0)
        
        i1, i2 = np.indices((freq[i]/5, NumSteps))
        i1 += 1
        i2 = i1+(i2+1)*ES
        
        r = r_cum[i2]-r_cum[i1-1]
        s = s_cum[i2]-s_cum[i1-1]
        if np.sum(r == t) > 0 : print 'warn'
        # Find the quality of the best fit.
        dchisq_tmp = (s**2.*t)/(r*(t-r))
        dchisq[i] = np.nanmax(dchisq_tmp, axis=(0,1))
        
    return dchisq

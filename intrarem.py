#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def trans_intrapixel(ind1, ind2, ind3, y, flux, eflux, maxiter=500, eps=1e-3, verbose=False):
    
    npoints = len(flux)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    weights = 1/eflux**2
    
    T = np.ones(np.amax(ind2)+1)
    a = np.zeros(np.amax(ind3)+1)
    b = np.zeros(np.amax(ind3)+1)
    
    for niter in range(maxiter):
        
        F = np.bincount(ind1, flux*T[ind2]*(a[ind3]*sn+b[ind3]*cs+1)*weights)/np.bincount(ind1, (T**2)[ind2]*(a[ind3]*sn+b[ind3]*cs+1)**2*weights)
        T = np.bincount(ind2, flux*F[ind1]*(a[ind3]*sn+b[ind3]*cs+1)*weights)/np.bincount(ind2, (F**2)[ind1]*(a[ind3]*sn+b[ind3]*cs+1)**2*weights)
        
        a = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(b[ind3]*cs+1))*F[ind1]*T[ind2]*sn*weights)/np.bincount(ind3, (F[ind1]*T[ind2]*sn)**2*weights)
        b = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(a[ind3]*sn+1))*F[ind1]*T[ind2]*cs*weights)/np.bincount(ind3, (F[ind1]*T[ind2]*cs)**2*weights)
    
        if (niter > 0):

            critF = np.nanmax(np.abs((F-F_old)/F_old))
            critT = np.nanmax(np.abs((T-T_old)/T_old))
            crita = np.nanmax(np.abs((a-a_old)/a_old))
            critb = np.nanmax(np.abs((b-b_old)/b_old))
            
            if (critF < eps) & (critT < eps) & (crita < eps) & (critb < eps):
                break
                
        if verbose:
            print 'niter = %i'%niter
            
        F_old = np.copy(F)
        T_old = np.copy(T)
        a_old = np.copy(a)
        b_old = np.copy(b)
        
    chi_tmp = (flux - F[ind1]*T[ind2]*(a[ind3]*sn+b[ind3]*cs+1))**2*weights
    chisq = np.sum(chi_tmp)/(npoints-npars)
            
    return F, T, a, b, niter, chisq, npoints, npars

def test(ind1, ind2, ind3, y, flux, eflux, maxiter=500, eps=1e-3, verbose=False):
    
    npoints = len(flux)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    weights = 1/eflux**2
    
    T = np.ones(np.amax(ind2)+1)
    a = np.zeros(np.amax(ind3)+1)
    b = np.zeros(np.amax(ind3)+1)
    
    for niter in range(maxiter):
        
        F = np.bincount(ind1, weights*(flux-a[ind3]*sn-b[ind3]*cs)*T[ind2])/np.bincount(ind1, weights*(T**2)[ind2])
        T = np.bincount(ind2, weights*(flux-a[ind3]*sn-b[ind3]*cs)*F[ind1])/np.bincount(ind2, weights*(F**2)[ind1])
        
        a = np.bincount(ind3, weights*(flux-F[ind1]*T[ind2]-b[ind3]*cs)*sn)/np.bincount(ind3, weights*sn**2)
        b = np.bincount(ind3, weights*(flux-F[ind1]*T[ind2]-a[ind3]*sn)*cs)/np.bincount(ind3, weights*cs**2)
    
        if (niter > 0):

            critF = np.nanmax(np.abs((F-F_old)/F_old))
            critT = np.nanmax(np.abs((T-T_old)/T_old))
            crita = np.nanmax(np.abs((a-a_old)/a_old))
            critb = np.nanmax(np.abs((b-b_old)/b_old))
            
            if (critF < eps) & (critT < eps) & (crita < eps) & (critb < eps):
                break
                
        if verbose:
            print 'niter = %i'%niter
                
            
        F_old = np.copy(F)
        T_old = np.copy(T)
        a_old = np.copy(a)
        b_old = np.copy(b)
        
    chi_tmp = weights*(flux - F[ind1]*T[ind2] - a[ind3]*sn - b[ind3]*cs)**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
            
    return F, T, a, b, niter, chisq, npoints, npars

#def bin_intrapixel(ind, y, flux, eflux, maxiter=500, eps=1e-3, verbose=False):
    
    #npoints = len(flux)
    #npars = len(np.unique(ind))
    
    #cs = np.cos(2*np.pi*y)
    #sn = np.sin(2*np.pi*y)
    #w = 1/eflux**2
    
    #a = np.zeros(np.amax(ind)+1)
    #b = np.zeros(np.amax(ind)+1)
    
    #niter = 0
    #while (niter < maxiter):
        
        #a = np.bincount(ind, (flux-(b[ind]*cs+1))*sn*w)/np.bincount(ind, (sn)**2*w)
        #b = np.bincount(ind, (flux-(a[ind]*sn+1))*cs*w)/np.bincount(ind, (cs)**2*w)
        
        #if (niter > 0):
            
            #crita = np.nanmax(np.abs((a-a_old)/a_old))
            #critb = np.nanmax(np.abs((b-b_old)/b_old))
            
            #if (crita < eps) & (critb < eps):
                #break
                
        #if verbose:
            #print 'niter = %i'%niter
            
        #a_old = np.copy(a)
        #b_old = np.copy(b)
    
        #niter += 1
        
    #chi_tmp = (flux - (a[ind]*sn+b[ind]*cs+1))**2*w
    #chisq = np.sum(chi_tmp)
    
    #return a, b, niter, chisq, npoints, npars

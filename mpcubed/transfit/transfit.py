# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:49:00 2017

@author: talens
"""

import numpy as np

import batman
import emcee
#from scipy import optimize

###############################################################################
### Parameter helper functions.
###############################################################################

def duration2axis(T14, P, p, b):    
    
    sinq = np.sin(np.pi*T14/P)
    a_sq = ((1 + p)**2 - b**2*(1 - sinq**2))/sinq**2

    return np.sqrt(a_sq)
    
def axis2duration(a, P, p, b):  
    
    sinq_sq = ((1 + p)**2 - b**2)/(a**2 - b**2)
    T14 = P/np.pi*np.arcsin(np.sqrt(sinq_sq))
    
    return T14

def impact2inc(b, a):
    
    return np.arccos(b/a)*180./np.pi

def inc2impact(i, a):
    
    return a*np.cos(i*np.pi/180.)
    
def axis2dens(a, P):
    
    G = 498.211    
    
    dens = (3.*np.pi*a**3)/(G*P**2)
    
    return dens

###############################################################################
### Transit model.
###############################################################################

def transit_circ(lc_pars, ld_pars, time=None):        
    
    F0, T0, P, T14, p, b = lc_pars    
    
    a = duration2axis(T14, P, p, b)    
    i = impact2inc(b, a)    
    
    params = batman.TransitParams()
    params.t0 = T0
    params.per = P
    params.rp = p      
    params.a = a                       
    params.inc = i
    params.ecc = 0.
    params.w = 90.
                      
    params.limb_dark = ld_pars[0]        
    params.u = ld_pars[1]
    
    if time is None:
        time = params.t0 + np.linspace(0, params.per, 500)
    
    m = batman.TransitModel(params, time)        
    
    phase = (time - params.t0)/params.per    
    model =  -2.5*np.log10(F0*m.light_curve(params))
    
    return phase, model

def baseline(mag, emag, model, mat):  
    
    pars = np.linalg.lstsq(mat/emag[:,None], (mag - model)/emag)[0]
    fit = np.sum(mat*pars, axis=1)    
    
    return pars, fit

def lnlike_circ(lc_pars, ld_pars, time, mag, emag, mat=None):

    # Priors on the parameters.
    if (lc_pars[0] < 0):
        return -np.inf
        
    if (lc_pars[2] < 0):
        return -np.inf
        
    if (lc_pars[3] < 0) | (lc_pars[3] > lc_pars[2]):
        return -np.inf
        
    if (lc_pars[4] < 0):
        return -np.inf
        
    if (lc_pars[5] < 0) | (lc_pars[5] > (1 + lc_pars[4])):
        return -np.inf

    # Evaluate the model.
    phase, model = transit_circ(lc_pars, ld_pars, time)

    # If LST is given fit the residuals with the long term and LST trend.
    if mat is not None:
        # Evaluate the baseline model.
        pars, fit = baseline(mag, emag, model, mat)
        
        # Compute the log-likelihood.
        lnlike = -.5*np.sum(((mag - model - fit)/emag)**2 + np.log(2*np.pi*emag**2))

    else:
        lnlike = -.5*np.sum(((mag - model)/emag)**2 + np.log(2*np.pi*emag**2))

    return lnlike    
    
def emcee_circ(lc_pars, ld_pars, time, mag, emag, nwalkers, nsteps, nthreads=4):
    
    npars = len(lc_pars)   
    nwalkers = np.maximum(nwalkers, 2*npars)

    # Initial values for the walkers.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = lc_pars[0] + 1e-4*np.random.randn(nwalkers) # F0
    pos[:,1] = lc_pars[1] + 1e-4*np.random.randn(nwalkers) # T0
    pos[:,2] = lc_pars[2] + 1e-4*np.random.randn(nwalkers) # P
    pos[:,3] = lc_pars[3] + 1e-4*np.random.randn(nwalkers) # T14
    pos[:,4] = lc_pars[4] + 1e-4*np.random.randn(nwalkers) # p
    pos[:,5] = lc_pars[5] + 1e-4*np.random.randn(nwalkers) # b
    
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnlike_circ, args=(ld_pars, time, mag, emag), threads=nthreads)
    sampler.run_mcmc(pos, nsteps)
    
    return sampler
 
def transit_ecc(lc_pars, ld_pars, time=None):        
    
    F0, T0, P, a, p, b, esinw, ecosw = lc_pars    
        
    i = impact2inc(b, a)    
    e = ecosw**2 + esinw**2
    
    if e < 1e-3:
        e = 0.
        w = 90.
    else:
        w = np.arctan2(esinw, ecosw)*180/np.pi     
    
    params = batman.TransitParams()
    params.t0 = T0
    params.per = P
    params.rp = p      
    params.a = a                       
    params.inc = i
    params.ecc = e
    params.w = w

    params.limb_dark = ld_pars[0]        
    params.u = ld_pars[1]
    
    if time is None:
        time = params.t0 + np.linspace(0, params.per, 100)
    
    m = batman.TransitModel(params, time)        
    
    phase = (time - params.t0)/params.per    
    model =  -2.5*np.log10(F0*m.light_curve(params))
    
    return phase, model
 
def lnlike_ecc(lc_pars, ld_pars, time, mag, emag):

    e = lc_pars[6]**2 + lc_pars[7]**2

    if (lc_pars[0] < 0):
        return -np.inf
        
    if (lc_pars[2] < 0):
        return -np.inf
        
    if (lc_pars[3]*(1. - e) < (1 + lc_pars[5])):
        return -np.inf        
        
    if (lc_pars[4] < 0):
        return -np.inf
        
    if (lc_pars[5] < 0) | (lc_pars[5] > (1 + lc_pars[4])):
        return -np.inf
        
    if (e > 1):
        return -np.inf

    phase, model = transit_ecc(lc_pars, ld_pars, time)
    lnlike = -.5*np.sum(((mag - model)/emag)**2 + np.log(2*np.pi*emag**2))

    return lnlike    
    
def emcee_ecc(lc_pars, ld_pars, time, mag, emag, nwalkers, nsteps, nthreads=4):
    
    npars = len(lc_pars)   
    nwalkers = np.maximum(nwalkers, 2*npars)

    # Initial values for the walkers.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = lc_pars[0] + 1e-4*np.random.randn(nwalkers) # F0
    pos[:,1] = lc_pars[1] + 1e-4*np.random.randn(nwalkers) # T0
    pos[:,2] = lc_pars[2] + 1e-4*np.random.randn(nwalkers) # P
    pos[:,3] = lc_pars[3] + 1e-4*np.random.randn(nwalkers) # a
    pos[:,4] = lc_pars[4] + 1e-4*np.random.randn(nwalkers) # p
    pos[:,5] = lc_pars[5] + 1e-4*np.random.randn(nwalkers) # b
    pos[:,6] = lc_pars[6] + 1e-4*np.random.randn(nwalkers) # esinw
    pos[:,7] = lc_pars[7] + 1e-4*np.random.randn(nwalkers) # ecosw
    
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnlike_ecc, args=(ld_pars, time, mag, emag), threads=nthreads)
    sampler.run_mcmc(pos, nsteps)
    
#    for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
#        if (i+1) % 100 == 0:
#            print("{0:5.1%}".format(float(i) / nsteps))    
        
    return sampler 
 
###############################################################################
### Transit model.
###############################################################################

def transit_blend(lc_pars, ld_pars, time=None):        
    
    F1, F0, T0, P, T14, p, b = lc_pars    
    
    a = duration2axis(T14, P, p, b)    
    i = impact2inc(b, a)    
    
    params = batman.TransitParams()
    params.t0 = T0
    params.per = P
    params.rp = p      
    params.a = a                       
    params.inc = i
    params.ecc = 0.
    params.w = 90.
                      
    params.limb_dark = ld_pars[0]        
    params.u = ld_pars[1]
    
    if time is None:
        time = params.t0 + np.linspace(0, params.per, 100)
    
    m = batman.TransitModel(params, time)        
    
    phase = (time - params.t0)/params.per    
    model =  -2.5*np.log10(F0*((1. - F1)*m.light_curve(params) + F1))
    
    return phase, model
 
def lnlike_blend(lc_pars, ld_pars, time, mag, emag):

    if (lc_pars[0] < 0) | (lc_pars[0] > 1):
        print 'F0'
        return -np.inf
        
    if (lc_pars[1] < 0):
        print 'F1'
        return -np.inf
        
    if (lc_pars[3] < 0):
        print 'P'
        return -np.inf
        
    if (lc_pars[4] < 0) | (lc_pars[4] > lc_pars[3]):
        print 'T14'
        return -np.inf
        
    if (lc_pars[5] < 0):
        print 'p'
        return -np.inf
        
    if (lc_pars[6] < 0) | (lc_pars[6] > (1 + lc_pars[5])):
        print 'b'
        return -np.inf

    phase, model = transit_blend(lc_pars, ld_pars, time)
    lnlike = -.5*np.sum(((mag - model)/emag)**2 + np.log(2*np.pi*emag**2))

    return lnlike    
    
def emcee_blend(lc_pars, ld_pars, time, mag, emag, nwalkers, nsteps, nthreads=4):
    
    npars = len(lc_pars)   
    nwalkers = np.maximum(nwalkers, 2*npars)

    # Initial values for the walkers.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = np.random.rand(nwalkers) # F1
    pos[:,1] = lc_pars[1] + 1e-4*np.random.randn(nwalkers) # F0
    pos[:,2] = lc_pars[2] + 1e-4*np.random.randn(nwalkers) # T0
    pos[:,3] = lc_pars[3] + 1e-4*np.random.randn(nwalkers) # P
    pos[:,4] = lc_pars[4] + 1e-4*np.random.randn(nwalkers) # T14
    pos[:,5] = lc_pars[5] + 1e-4*np.random.randn(nwalkers) # p
    pos[:,6] = np.random.rand(nwalkers) # b
    
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnlike_blend, args=(ld_pars, time, mag, emag), threads=nthreads)
    sampler.run_mcmc(pos, nsteps)
    
    return sampler 
 
def main(): 
    
    import time    
    
    lc_pars = [1., 0., 2.15, 5, .1, .1, .05, .05]  
    ld_pars = ['quadratic', [0., 0.]]    
    
    start = time.time()    
    
    for i in range(500*7000):    
    
        transit_ecc(lc_pars, ld_pars)    
    
    print time.time() - start    
    
    return
    
if __name__ == '__main__':
    
    main()

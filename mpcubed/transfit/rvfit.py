# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:49:00 2017

@author: talens
"""

import numpy as np

import emcee

###############################################################################
### Helper functions.
###############################################################################

def eccentric_anomaly(M, e, dtol=1e-5):
    
    E0, E1 = M, M    
    while np.any(np.abs(M - E1 + e*np.sin(E1)) > dtol):
        
        E0 = E1
        E1 = E0 - (E0 - e*np.sin(E0) - M)/(1. - e*np.cos(E0))
    
    return E1
    
def true_anomaly(E, e):
    
    cos_nu = (np.cos(E) - e)/(1. - e*np.cos(E))  
    sin_nu = (np.sqrt(1 - e**2)*np.sin(E))/(1 - e*np.cos(E)) 
    
    nu = np.arctan2(sin_nu, cos_nu) 
    
    return nu

###############################################################################
### Radial velocity model.
###############################################################################

def rv_model(rv_pars, time=None):
    
    T0, P, K, gamma, esinw, ecosw = rv_pars   
    
    #
    e = ecosw**2 + esinw**2
    if e < 1e-3:
        w = np.pi/2.
    else:
        w = np.arctan2(esinw, ecosw)    
    
    nu0 = np.pi/2 - w
    E0 = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(nu0/2))
    M0 = E0 - e*np.sin(E0)
    Tp = T0 - P/(2*np.pi)*M0    
    
    if time is None:
        time = T0 + np.linspace(0, P, 500)
    
    phase = (time - T0)/P
    
    # Compute the various anomalies.
    M = 2*np.pi*(time - Tp)/P
    E = eccentric_anomaly(M, e)
    nu = true_anomaly(E, e)
    
    model = K*(np.cos(w + nu) + e*np.cos(w)) + gamma

    return phase, model
    
def lnlike_rv(rv_pars, time, vel, evel):
    
    if (rv_pars[1] < 0):
        return -np.inf, np.inf
    
    if (rv_pars[2] < 0):
        return -np.inf, np.inf
    
    if ((rv_pars[4]**2 + rv_pars[5]**2) > 1):
        return -np.inf, np.inf

    phase, model = rv_model(rv_pars, time)
    chisq = np.sum((vel - model)**2/evel**2)
    lnlike = -.5*np.sum(((vel - model)/evel)**2 + np.log(2*np.pi*evel**2))
    
    return lnlike, chisq
    
def emcee_rv(rv_pars, time, vel, evel, nwalkers, nsteps, nthreads=4):
    
    npars = len(rv_pars)   
    nwalkers = np.maximum(nwalkers, 2*npars)

    # Initial values for the walkers.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = rv_pars[0] + 1e-4*np.random.randn(nwalkers) # T0
    pos[:,1] = rv_pars[1] + 1e-4*np.random.randn(nwalkers) # P
    pos[:,2] = rv_pars[2] + 1e-4*np.random.randn(nwalkers) # K
    pos[:,3] = rv_pars[3] + 1e-4*np.random.randn(nwalkers) # gamma
    pos[:,4] = (2*np.random.rand(nwalkers) - 1)/np.sqrt(2) # ecosw
    pos[:,5] = (2*np.random.rand(nwalkers) - 1)/np.sqrt(2) # esinw
    
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnlike_rv, args=(time, vel, evel), threads=nthreads)
    sampler.run_mcmc(pos, nsteps)
    
    return sampler

###############################################################################
### Fitting a circular orbit model.
###############################################################################

def rv_circ(rv_pars, time=None):
    
    rv_pars = np.asarray(rv_pars)
    rv_pars = np.append(rv_pars, [0., 0.])
    
    return rv_model(rv_pars, time)
    
def lnlike_circ(rv_pars, time, vel, evel):
    
    rv_pars = np.asarray(rv_pars)
    rv_pars = np.append(rv_pars, [0., 0.])
    
    return lnlike_rv(rv_pars, time, vel, evel)
    
def emcee_circ(rv_pars, time, vel, evel, nwalkers, nsteps, nthreads=4):
    
    npars = len(rv_pars)   
    nwalkers = np.maximum(nwalkers, 2*npars)

    # Initial values for the walkers.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = rv_pars[0] + 1e-4*np.random.randn(nwalkers) # T0
    pos[:,1] = rv_pars[1] + 1e-4*np.random.randn(nwalkers) # P
    pos[:,2] = rv_pars[2] + 1e-4*np.random.randn(nwalkers) # K
    pos[:,3] = rv_pars[3] + 1e-4*np.random.randn(nwalkers) # gamma
    
    sampler = emcee.EnsembleSampler(nwalkers, npars, lnlike_circ, args=(time, vel, evel), threads=nthreads)
    sampler.run_mcmc(pos, nsteps)
    
    return sampler
    
def main():

    return    
    
if __name__ == '__main__':
    
    main()

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 15:14:19 2018

@author: talens
"""

import h5py
import numpy as np

import emcee
import batman

import corner
import matplotlib.pyplot as plt

from . import misc, statistics
from .detection import detrend

###############################################################################
### Parameter helper functions.
###############################################################################

def T142axis(T14, P, p, b):
    
    sin_sq = np.sin(np.pi*T14/P)**2
    a_sq = ((1 + p)**2 - b**2*(1 - sin_sq))/sin_sq
    
    return np.sqrt(a_sq)

def axis2T14(a, P, p, b):
    
    sin_sq = ((1 + p)**2 - b**2)/(a**2 - b**2)
    T14 = P/np.pi*np.arcsin(np.sqrt(sin_sq))
    
    return T14

def T232axis(T23, P, p, b):
    
    sin_sq = np.sin(np.pi*T23/P)**2
    a_sq = ((1 - p)**2 - b**2*(1 - sin_sq))/sin_sq
    
    return np.sqrt(a_sq)

def axis2T23(a, P, p, b):
    
    sin_sq = ((1 - p)**2 - b**2)/(a**2 - b**2)
    T23 = P/np.pi*np.arcsin(np.sqrt(sin_sq))
    
    return T23

def impact2inc(b, a):
    
    inc = np.arccos(b/a)
    
    return inc

def inc2impact(inc, a):
    
    b = a*np.cos(inc)
    
    return b

def xy2eccomg(x, y):
    
    ecc = x**2. + y**2.
    omg = np.arctan2(y, x)
    
    return ecc, omg

def eccomg2xy(ecc, omg):
    
    x = np.sqrt(ecc)*np.cos(omg)
    y = np.sqrt(ecc)*np.sin(omg)
    
    return x, y

def axis2dens(a, P):
    
    G = 498.211    
    
    rho = (3.*np.pi*a**3)/(G*P**2)
    
    return rho
    
def dens2axis(rho, P):
    
    G = 498.211      
    
    a = (rho*G*P**2/(3.*np.pi))**(1./3)
    
    return a

def time2phase(time, P, T0):
    
    phase = (time - T0)/P
    phase = np.mod(phase+.5, 1.)-.5
    
    return phase

def mean_anomaly(time, P, Tp):
    
    M = 2.*np.pi*time2phase(time, P, Tp)
    
    return M

def _mean2eccentric(M, ecc, E, dtol=1e-10):
    """ Use Newton's method to solve Kepler's equation. """
    
    dE = 2*dtol
    while np.any(np.abs(dE) > dtol):

        dE = (M + ecc*np.sin(E) - E)/(1. - ecc*np.cos(E))
        E = E + dE

    return E

def mean2eccentric(M, ecc, dtol=1e-10):
    """ Solve Kepler's equation following Mikola (1987) """
    
    # For circular orbits.
    if ecc == 0.:
        return M
    
    # Obtain a first order approximation for E.
    aux = 4.*ecc + .5
    alpha = (1. - ecc)/aux
    beta = M/(2.*aux)
    
    aux = np.sqrt(beta**2 + alpha**3)
    z = beta + aux
    z = np.where(z <= 0.0, beta - aux, z)

    z = np.abs(z)**(1./3)
    s0 = z - alpha/z
    s1 = s0 - (0.078*s0**5)/(1. + ecc)
    E0 = M + ecc*(3.*s1 - 4.*s1**3)

    # Compute a higher order correction.
    se0 = np.sin(E0)
    ce0 = np.cos(E0)

    f = E0 - ecc*se0 - M
    f1 = 1. - ecc*ce0
    f2 = ecc*se0
    f3 = ecc*ce0
    f4 = -f2
    u1 = -f/f1
    u2 = -f/(f1 + 1./2*f2*u1)
    u3 = -f/(f1 + 1./2*f2*u2 + 1./6*f3*u2**2)
    u4 = -f/(f1 + 1./2*f2*u3 + 1./6*f3*u3**2 + 1./24*f4*u3**3)

    E = E0 + u4

    # If necessary, perform further refinement.
    M_ = eccentric2mean(E, ecc)
    if np.any(np.abs(M - M_) > dtol):
        E = _mean2eccentric(M, ecc, E, dtol)
    
    return E

def eccentric2mean(E, ecc):
    
    M = E - ecc*np.sin(E)
    
    return M

def eccentric2true(E, ecc):
    
    x = np.sqrt(1. - ecc)*np.cos(E/2.)
    y = np.sqrt(1. + ecc)*np.sin(E/2.)
    
    nu = 2.*np.arctan2(y, x)

    return nu
    
def true2eccentric(nu, ecc):
    
    x = np.sqrt(1. + ecc)*np.cos(nu/2.)
    y = np.sqrt(1. - ecc)*np.sin(nu/2.)
    
    E = 2.*np.arctan2(y, x) 
    
    return E

def transit2elliptic(lc_pars):
    
    lc_pars = np.asarray(lc_pars)
    
    # Add sqrt(e)cos(w) and sqrt(e)sin(w).
    lc_pars = np.append(lc_pars, np.zeros(2))
    
    # Convert T14 to a/R_*.
    lc_pars[2] = T142axis(lc_pars[2], lc_pars[1], lc_pars[3], lc_pars[4])
    
    return lc_pars

def radvel2elliptic(rv_pars):
    
    rv_pars = np.asarray(rv_pars)
    
    # Add sqrt(e)cos(w) and sqrt(e)sin(w).
    rv_pars = np.append(rv_pars, np.zeros(2))
    
    return rv_pars

###############################################################################
### The transit model and fitting functions.
###############################################################################

def transit_model(time, lc_pars, ld_pars=[0.6], ld_type='linear'):
    """A full transit model."""
    
    if (len(lc_pars) == 5):
        lc_pars = transit2elliptic(lc_pars)
    elif (len(lc_pars) != 7):
        raise ValueError('lc_pars must have length 5 or 7')
    
    T0, P, a, p, b, x, y = lc_pars

    # Convert parameters.
    inc = impact2inc(b, a)
    ecc, omg = xy2eccomg(x, y)

    # Create the batman TransitParams object.
    params = batman.TransitParams()
    params.t0 = 0.
    params.per = P
    params.rp = p      
    params.a = a                       
    params.inc = inc*180/np.pi
    params.ecc = ecc
    params.w = omg*180/np.pi
                      
    params.limb_dark = ld_type     
    params.u = ld_pars
    
    # Evaluate the model.
    m = batman.TransitModel(params, time - T0)        
    
    phase = time2phase(time, P, T0)
    model = m.light_curve(params)
    model = -2.5*np.log10(model)
    
    return phase, model

def transit_evaluate(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}):
    
    # Evaluate the transit model.
    phase, model = transit_model(lc['jd'], lc_pars, ld_pars, ld_type)

    # Evaluate the baseline model.
    lc = detrend.remove_trend(lc, nobs, model=model, method=method, options=options)
    
    return phase, model, lc

def transit_chisq(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}):  
    """Compute the chi-squared value."""

    # Evaluate the transit and baseline models.
    phase, model, lc = transit_evaluate(lc_pars, lc, nobs, ld_pars, ld_type, method, options)
    
    # Compute the chi-square and log-likelihood.
    chisq = np.sum(((lc['mag'] - lc['trend'] - model)/lc['emag'])**2)

    return chisq

def transit_lstsq(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}):
    """Find the best-fit parameters."""
    
    from scipy.optimize import minimize
    
    if (len(lc_pars) == 5):
        bounds = [(None, None), (0, None), (0, None), (0, None), (0, 1)]
    else:
        ValueError('lc_pars must have length 5, transit_lstsq only solves circular orbits')

    res = minimize(transit_chisq, lc_pars, args=(lc, nobs, ld_pars, ld_type, method, options), bounds=bounds)
    
    return res.x

def transit_prior(lc_pars):
    
    if (len(lc_pars) == 5):
        
        # Unpack the transit parameters.
        T0, P, T14, p, b = lc_pars
        
        # Priors on the parameters.
        if (P < 0):
            return True
            
        if (T14 > P/2.): # Equivalent to a < (1 + p).
            return True
            
        if (p < 0.) | (p > 1.):
            return True
            
        if (b < 0.) | (b > (1. + p)):
            return True
        
    elif (len(lc_pars) == 7):
        
        # Unpack the transit parameters.
        T0, P, a, p, b, x, y = lc_pars  
        
        # Convert parameters.
        ecc, omg = xy2eccomg(x, y)
        
        # Priors on the parameters.
        if (P < 0):
            return True
            
        if (a*(1 - ecc) < (1 + p)):
            return True
            
        if (p < 0) | (p > 1):
            return True
            
        if (b < 0) | (b > (1 + p)):
            return True
    
        if (ecc < 0) | (ecc > 0.999): # Avoids e=1. Not clean, but works. 
            return True
        
    else:
        raise ValueError('lc_pars must have length 5 or 7')

    return False

def transit_lnlike(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}):  
    """Compute the log-likelihood value."""

    if transit_prior(lc_pars):
        return -np.inf
        
    # Evaluate the transit and baseline models.
    phase, model, lc = transit_evaluate(lc_pars, lc, nobs, ld_pars, ld_type, method, options)
    
    # Compute the chi-square and log-likelihood.
    lnlike = -0.5*np.sum(np.log(2.*np.pi*lc['emag']**2) + ((lc['mag'] - lc['trend'] - model)/lc['emag'])**2)

    return lnlike

def transit_emcee(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}, **kwargs):
    """Find the best-fit parameters."""

    # Set up the model.
    npars = len(lc_pars) 
    nwalkers = kwargs.pop('nwalkers', 4*npars)
    nsteps = kwargs.pop('nsteps', 300*npars)
    threads = kwargs.pop('threads', 6)
    
    nwalkers = np.maximum(nwalkers, 2*npars)
    
    # Initial walker states.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = lc_pars[0] + 2e-3*np.random.randn(nwalkers) # T0 [JD]
    pos[:,1] = lc_pars[1] + 2e-5*np.random.randn(nwalkers) # P [days]
    pos[:,2] = lc_pars[2] + 2e-3*np.random.randn(nwalkers) # T14 [days] or a/R_*
    pos[:,3] = lc_pars[3] + 2e-2*np.random.randn(nwalkers) # R_p/R_*
    pos[:,4] = np.sqrt(0.75)*np.random.rand(nwalkers) # b = a/R_*cos(i)
    
    if (npars == 7):
    
        ecc = np.sqrt(0.75)*np.random.rand(nwalkers) # e
        omg = 2*np.pi*np.random.rand(nwalkers) # w
        
        x, y = eccomg2xy(ecc, omg)
        pos[:,5] = x # sqrt(e)cos(w)
        pos[:,6] = y # sqrt(e)sin(w)
    
    # Run the model.
    sampler = emcee.EnsembleSampler(nwalkers, npars, transit_lnlike, args=(lc, nobs, ld_pars, ld_type, method, options), threads=threads)
    sampler.run_mcmc(pos, nsteps)
    
    f_acc = sampler.acceptance_fraction
    chain = sampler.chain
    lnprob = sampler.lnprobability
    
    return chain, lnprob, f_acc

###############################################################################
### The radial velocity model and fitting functions.
###############################################################################

def radvel_model(time, rv_pars):
    
    if (len(rv_pars) == 4):
        rv_pars = radvel2elliptic(rv_pars)
    elif (len(rv_pars) != 6):
        raise ValueError('rv_pars must have length 4 or 6')
    
    Tp, P, K, gamma, x, y = rv_pars

    # Convert parameters.
    ecc, omg = xy2eccomg(x, y)

    M = mean_anomaly(time, P, Tp)
    E = mean2eccentric(M, ecc)
    nu = eccentric2true(E, ecc)
    
    phase = time2phase(time, P, Tp)
    model = K*(np.cos(omg + nu) + ecc*np.cos(omg)) + gamma
    
    return phase, model

def radvel_chisq(rv_pars, time, rv, erv):  
    """Compute the chi-squared value."""

    # Evaluate the radial velocity model.
    phase, model = radvel_model(time, rv_pars)
    
    # Compute the chi-square and log-likelihood.
    chisq = np.sum(((rv - model)/erv)**2)

    return chisq

def radvel_lstsq(rv_pars, time, rv, erv):
    """Find the best-fit parameters."""
    
    from scipy.optimize import minimize
    
    if (len(rv_pars) == 4):
        bounds = [(None, None), (0, None), (None, None), (None, None)]
    else:
        ValueError('rv_pars must have length 4, radvel_lstsq only solves circular orbits')

    res = minimize(radvel_chisq, rv_pars, args=(time, rv, erv), bounds=bounds)
    
    return res.x

def radvel_prior(rv_pars):

    if (len(rv_pars) == 4):
        rv_pars = radvel2elliptic(rv_pars)
    elif (len(rv_pars) != 6):
        raise ValueError('rv_pars must have length 4 or 6')

    # Unpack the parameters.
    Tp, P, K, gamma, x, y = rv_pars
    
    # Convert parameters.
    ecc, omg = xy2eccomg(x, y)
    
    # Priors on the parameters.
    if (P < 0):
        return True

    if (ecc < 0) | (ecc > 0.999): # Avoids e=1. Not clean, but works. 
        return True
    
    return False

def radvel_lnlike(rv_pars, time, rv, erv):  
    """Compute the log-likelihood value."""

    if radvel_prior(rv_pars):
        return -np.inf
        
    # Evaluate the radial velocity model.
    phase, model = radvel_model(time, rv_pars)
    
    # Compute the chi-square and log-likelihood.
    lnlike = -0.5*np.sum(np.log(2.*np.pi*erv**2) + ((rv - model)/erv)**2)

    return lnlike

def radvel_emcee(rv_pars, time, rv, erv, **kwargs):
    """Find the best-fit parameters of a circular transit model."""

    # Set up the model.
    npars = len(rv_pars) 
    nwalkers = kwargs.pop('nwalkers', 4*npars)
    nsteps = kwargs.pop('nsteps', 300*npars)
    threads = kwargs.pop('threads', 6)
    
    nwalkers = np.maximum(nwalkers, 2*npars)
    
    # Initial walker states.
    pos = np.zeros((nwalkers, npars))
    pos[:,0] = rv_pars[0] + 2e-3*np.random.randn(nwalkers) # Tp [JD]
    pos[:,1] = rv_pars[1] + 2e-5*np.random.randn(nwalkers) # P [days]
    pos[:,2] = rv_pars[2] + 2e-3*np.random.randn(nwalkers) # K [km/s]
    pos[:,3] = rv_pars[3] + 2e-2*np.random.randn(nwalkers) # gamma [km/s] 
    
    if (npars == 6):
    
        ecc = np.sqrt(0.75)*np.random.rand(nwalkers) # e
        omg = 2*np.pi*np.random.rand(nwalkers) # w
        
        x, y = eccomg2xy(ecc, omg)
        pos[:,4] = x # sqrt(e)cos(w)
        pos[:,5] = y # sqrt(e)sin(w)
    
    # Run the model.
    sampler = emcee.EnsembleSampler(nwalkers, npars, radvel_lnlike, args=(time, rv, erv), threads=threads)
    sampler.run_mcmc(pos, nsteps)
    
    f_acc = sampler.acceptance_fraction
    chain = sampler.chain
    lnprob = sampler.lnprobability
    
    return chain, lnprob, f_acc

###############################################################################
### For reading MASCARA photometry for fitting.
###############################################################################
    
def _apply_mask(lc, nobs):
    
    mask = lc['mask']
    nobs_new = np.zeros_like(nobs)

    strides = np.append(0, np.cumsum(nobs))

    # Loop over blocks of data.
    for i in range(nobs.size):
        
        i1 = strides[i]
        i2 = strides[i+1]
        
        nobs_new[i] = np.sum(mask[i1:i2])
        
    lc_new = lc[mask]
    nobs_new = nobs_new[nobs_new > 0]
    
    return lc_new, nobs_new

def read_mascara(filelist, ascc, ra, dec, aper=0, method='legendre', options={'maxiter':1}):
    
    from .detection.boxlstsq import read_data
    
    # Read the data.
    time, lc2d, nobs = read_data(filelist, [ascc], aper=aper)
    
    lc = lc2d[:,0]
    
    # Perform barycentric correction.
    lc['jd'] = misc.barycentric_dates(lc['jd'], ra, dec)
    
    # Initial trend and remove masked values.
    lc = detrend.remove_trend(lc, nobs, method=method, options=options)
    lc, nobs = _apply_mask(lc, nobs)
    
    return lc, nobs

###############################################################################
### To prepare non-MASCARA photometry for fitting.
###############################################################################
    
def phot2lc(jd, mag, emag):
    
    nobs = np.array([jd.size])
    
    dtype = [('jd', 'float64'), ('mag', 'float64'), ('emag', 'float64'), ('trend', 'float64'), ('mask', 'bool')]

    lc = np.recarray((jd.size,), dtype=dtype)
    lc[:] = 0
    
    lc['jd'] = jd
    lc['mag'] = mag
    lc['emag'] = emag

    return lc, nobs

###############################################################################
### For reading and writing MCMC results.
###############################################################################

def write_emcee(filename, chain, lnprob, f_acc, attrs={}):
    
    with h5py.File(filename) as f:
        
        grp = f.create_group('header')
        
        for key in attrs.keys():
            grp.attrs[key] = attrs[key]
        
        grp = f.create_group('data')
        
        grp.create_dataset('chain', data=chain)
        grp.create_dataset('lnprob', data=lnprob)
        grp.create_dataset('f_acc', data=f_acc)
    
    return

def read_emcee(filename):
    
    with h5py.File(filename, 'r') as f:
        
        grp = f['header']
        
        attrs = {}
        for key in grp.attrs.keys():
            attrs[key] = grp.attrs[key]
        
        grp = f['data']
        
        chain = grp['chain'][()]
        lnprob = grp['lnprob'][()]
        f_acc = grp['f_acc'][()]
        
    return chain, lnprob, f_acc, attrs

###############################################################################
# For determining the best-fit values.
###############################################################################

def round2significance(value, error1, error2=None):
    """Round values to significant digits based on the uncertainties."""
    
    if error2 is None:
        error2 = -error1
    
    ndigits1 = -np.floor(np.log10(error1)).astype('int')  
    ndigits2 = -np.floor(np.log10(-error2)).astype('int')
    
    ndigits = np.maximum(ndigits1, ndigits2) 
    
    value = np.around(value, ndigits)
    error1 = np.around(error1, ndigits)
    error2 = np.around(error2, ndigits)
    
    return ndigits, value, error1, error2

def best_values(flatchain):
    """Get the best fit values from a flattened MCMC chain."""
    
    m50 = np.percentile(flatchain, 50, axis=0)
    m68 = np.percentile(np.abs(flatchain - m50), 68, axis=0) 
    
    m16 = np.percentile(flatchain, 16, axis=0)
    m84 = np.percentile(flatchain, 84, axis=0)   
    
    return m50, m68, m84 - m50, m16 - m50

###############################################################################
### For visualization of the data and model.
###############################################################################

def plot_radvel(rv_pars, time, rv, erv):
    
    phase, model = radvel_model(time, rv_pars)
    
    plt.figure(figsize=(5, 3))

    plt.subplot(111)
    
    # Plot the calibrated data.
    plt.errorbar(phase, rv, erv, marker='o', c='k', ls='none')
    
    # Plot the model.
    phase = np.linspace(0, 1, 1001)
    time = rv_pars[0] + rv_pars[1]*phase
    phase, model = radvel_model(time, rv_pars)
    
    sort = np.argsort(phase)
    plt.plot(phase[sort], model[sort], c=(146./255, 0., 0.), lw=2, zorder=20)
    
    plt.xlim(-.5, .5)
    
    plt.xlabel('Phase')
    plt.ylabel('Radial Velocity [km/s]')
    
    plt.tight_layout()

    plt.show()
    plt.close()
    
    return

def plot_transit(lc_pars, lc, nobs, ld_pars=[0.6], ld_type='linear', method='legendre', options={'maxiter':1}):
    
    if len(lc_pars) == 5:
        T0, P, T14, p, b = lc_pars
    elif len(lc_pars) == 7:
        T0, P, a, p, b, x, y = lc_pars
        T14 = axis2T14(a, P, p, b) # Not exact, but good enough for plotting.
    else:
        raise ValueError('lc_pars must have length 5 or 7')
    
    phase, model, lc = transit_evaluate(lc_pars, lc, nobs, ld_pars, ld_type, method, options)
    
    plt.figure(figsize=(8, 3))

    plt.subplot(111)
    
    # Plot the calibrated data.
    plt.scatter(phase, lc['mag'] - lc['trend'], c='grey', marker='.', edgecolors='none', alpha=0.5)
    
    # Plot the calibrated data, binned in phase.
    nbins = np.ceil(9*P/T14)
    bins = np.linspace(-.5, .5, nbins+1)   
    xbin, ybin, eybin = statistics.bin_data(phase, lc['mag'] - lc['trend'], bins, lc['emag'])
    
    plt.errorbar(xbin, ybin, eybin, fmt='o', c='k')
    
    # Plot the model.
    phase = np.linspace(0, 1, 1001)
    time = lc_pars[1] + lc_pars[0]*phase
    phase, model = transit_model(time, lc_pars, ld_pars, ld_type)
    
    sort = np.argsort(phase)
    plt.plot(phase[sort], model[sort], c=(146./255, 0., 0.), lw=2, zorder=20)
        
    plt.xlim(-.5, .5)
    plt.ylim(0.05, -0.02)
    
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    plt.tight_layout()

    plt.show()
    plt.close()
    
    return

###############################################################################
### For visualization of emcee results.
###############################################################################

def plot_walkers(data, label):
    
    plt.plot(data, c='k', alpha=0.5)
    
    plt.xlabel('Walker Step')
    plt.ylabel(label)
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return

def plot_acceptance(f_acc, f_min):
    
    plt.hist(f_acc, bins=np.linspace(0,1,21))
    plt.axvline(f_min, c='k')
    
    plt.xlabel('$f_{acc}$')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return

def plot_corner(array, labels, figname=None):
        
    corner.corner(array, labels=labels)
    
    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname)
    plt.show()
    plt.close()
    
    return

def results_emcee(labels, chain, lnprob, f_acc, f_min=0.20, fburn=0.4):
    
    nwalkers, nsteps, npars = chain.shape
    nburn = int(fburn*nsteps)
    
    # Plot the walkers for inspection.
    for i in range(npars):
        plot_walkers(chain[:,:,i].T, labels[i])
        
    plot_walkers(lnprob.T, '$\ln L$')
    
    # Plot the acceptance fractions.
    plot_acceptance(f_acc, f_min)
    
    # Remove burn-in and bad walkers.
    good = f_acc > f_min
    flatchain = chain[good,nburn:].reshape((-1,chain.shape[2]))
    flatlnprob = lnprob[good,nburn:].ravel()

    # Make the corner plot.
    array = np.column_stack([flatchain, flatlnprob])
    plot_corner(array, labels + ['$\ln L$'])
    
    pars, p68, p16, p84 = best_values(flatchain)
    
    # Best-fit values and uncertainties.
    for i in range(npars):
        
        ndigits, best, epos, eneg = round2significance(pars[i], p16[i], p84[i])    
    
        print labels[i], ndigits, best, epos, eneg
    
    return pars, p68, p16, p84

def main():
    
    return

if __name__ == '__main__':
    
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np

import matplotlib.pyplot as plt

import emcee
import corner
import scipy.optimize as op

import IO
import models

sol_c = 2.9979e5

def regular_velgrid(wave, spec, window, deltav, pad=200.):
    """ Put a spectrum on a wavelength grid with constant velocity spacing. """   
    
    # Compute the edges of the window in velocity space.
    cwave = (window[0] + window[1])/2.
    v1 = sol_c*np.log(window[0]/cwave)    
    v2 = sol_c*np.log(window[1]/cwave)    
    
    # Round to the nearest multiple of dv and pad.
    vmin = np.floor(v1/deltav)*deltav - pad
    vmax = np.ceil(v2/deltav)*deltav + pad  
    
    # Compute the new wavelength grid.
    v = np.linspace(vmin, vmax, (vmax - vmin)/deltav + 1)  
    new_wave = cwave*np.exp(v/sol_c)
    
    # Interpolate the spectrum to the new grid.
    new_spec = np.interp(new_wave, wave, spec)

    return new_wave, new_spec    
    
def fit_phoenix(specfile, modeldir, window, temps, deltav=.5, order=3):
    
    # Get the model spectra.
    modlist = glob.glob(os.path.join(modeldir, '*'))
    modlist = np.sort(modlist)    
    
    # Read the data.
    wave_data, spec_data, hdr_data = IO.read_hermes(specfile)
    mask = (wave_data > window[0]) & (wave_data < window[1])
    wave_data, spec_data = wave_data[mask], spec_data[mask]
    
    phxteff = np.zeros(len(modlist))
    phxlogg = np.zeros(len(modlist))
    phxZ = np.zeros(len(modlist))
    chisq = np.zeros(len(modlist))
    pars = np.zeros((len(modlist), 2))    
    
    for i in range(len(modlist)):    
    
        modfile = modlist[i]  
    
        # Read the model.
        wave_mod, spec_mod, hdr_mod = IO.read_phoenix(modfile)
        wave_mod, spec_mod = regular_velgrid(wave_mod, spec_mod, window, deltav)
        
        phxteff[i] = hdr_mod['PHXTEFF']
        phxlogg[i] = hdr_mod['PHXLOGG']
        phxZ[i] = hdr_mod['PHXM_H']        
        
        # Check that the temperature of the model.
        if not (temps[0] <= hdr_mod['PHXTEFF'] <= temps[1]):
            chisq[i] = -1
            continue        
        
        # Fit the model to the data.
        pars0 = [0., 10.]
        bounds = [(-150, 150), (1.5, 150)]
        res = op.minimize(models.spec_chisq, pars0, args=(wave_data, spec_data, wave_mod, spec_mod, deltav, order), bounds=bounds)
        
        chisq[i] = res['fun']
        pars[i] = res['x']
        
    mask = (chisq > 0)    
    modlist = modlist[mask]
    chisq = chisq[mask]
    pars = pars[mask]
    phxteff = phxteff[mask]
    phxlogg = phxlogg[mask]
    phxZ = phxZ[mask]        
        
    return modlist, chisq, pars, phxteff, phxlogg, phxZ

def refine_phoenix(specfile, modfile, pars, window, deltav=.5, order=3):
    
    # Read the data.
    wave_data, spec_data, hdr_data = IO.read_hermes(specfile)
    mask = (wave_data > window[0]) & (wave_data < window[1])
    wave_data, spec_data = wave_data[mask], spec_data[mask]
    
    # Read the model.
    wave_mod, spec_mod, hdr_mod = IO.read_phoenix(modfile)
    wave_mod, spec_mod = regular_velgrid(wave_mod, spec_mod, window, deltav)
    
    # Run the MCMC chain.
    ndim, nwalkers = 2, 4
    pos = [pars + 1e-1*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, models.spec_lnlike, args=(wave_data, spec_data, wave_mod, spec_mod, deltav, order))
    sampler.run_mcmc(pos, 1000)
    
    return sampler

def plot_phoenix(specfile, modfile, pars, window, deltav=.5, order=3):
    
    # Read the data.
    wave_data, spec_data, hdr_data = IO.read_hermes(specfile)
    mask = (wave_data > window[0]) & (wave_data < window[1])
    wave_data, spec_data = wave_data[mask], spec_data[mask]
    
    # Read the model.
    wave_mod, spec_mod, hdr_mod = IO.read_phoenix(modfile)
    wave_mod, spec_mod = regular_velgrid(wave_mod, spec_mod, window, deltav)    
    
    # Evaluate the fit.
    fit = models.spec_model(pars, wave_data, spec_data, wave_mod, spec_mod, deltav, order)
    
    plt.figure(figsize=(16,5))
    plt.plot(wave_data, spec_data, 'k')
    plt.fill_between(wave_data, spec_data - np.sqrt(spec_data), spec_data + np.sqrt(spec_data), facecolor='k', alpha=.5)
    plt.plot(wave_data, fit, 'r')
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel('Flux [Counts]')
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return

if __name__ == '__main__':
    
    specfile = '/data2/talens/Spectra/HERMES/HD201585/00717636_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'
    modeldir = '/data2/talens/PHOENIX/Z-0.0/'
    
    modlist, chisq, pars, phxteff, phxlogg, phxZ = fit_phoenix(specfile, modeldir, [4500, 6500], [7000, 8000])
    
    plt.scatter(phxteff + 50*phxZ, phxlogg, c=np.log10(chisq/np.amin(chisq)), cmap=plt.cm.Greys, s=40)
    plt.xlabel(r'$T_{\mathrm{eff}}$ [K]', size='x-large')
    plt.ylabel(r'$\log g$', size='x-large')
    cb = plt.colorbar()
    cb.set_label(r'$\log_{10}(\chi^2/\chi_{\mathrm{min}}^2)$', size='x-large')
    plt.tight_layout()
    plt.show()    
    
    arg = np.argmin(chisq)
    print modlist[arg]
    print pars[arg]
    print phxteff[arg], phxlogg[arg], phxZ[arg]    
    
    exit()    
    
    arg = np.argmin(chisq)
    print pars[arg]
    plot_phoenix(specfile, modlist[arg], pars[arg], window=[5000, 5300])    
    plot_phoenix(specfile, modlist[arg], pars[arg], window=[4700, 5000]) 
    plot_phoenix(specfile, modlist[arg], pars[arg], window=[5300, 5600]) 
    
    sampler = refine_phoenix(specfile, modlist[arg], pars[arg], window=[5000, 5300])    
    
    print sampler.chain.shape    
    
    plt.subplot(211)
    plt.plot(sampler.chain[:,:,0].T, 'k')

    plt.subplot(212)
    plt.plot(sampler.chain[:,:,1].T, 'k')

    plt.tight_layout()
    plt.show()    
    
    samples = sampler.chain[:, 100:, :].reshape((-1, 2))  
    pars = np.percentile(samples, 50, axis=0)
    print pars    
    
    labels = [r'$v_{\mathrm{rad}}$', r'$v\sin i$']
    corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True)
    plt.tight_layout()
    plt.show()  

    plot_phoenix(specfile, modlist[arg], pars, window=[5000, 5300])    
    plot_phoenix(specfile, modlist[arg], pars, window=[4700, 5000]) 
    plot_phoenix(specfile, modlist[arg], pars, window=[5300, 5600])     
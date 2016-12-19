#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

sol_c = 2.9979e5

def blackbody(wave, temp):
    """ Given a wavelength in angstrom and temperature in Kelvin
        returns the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/a.
    """

    w = wave/1e8 # Angstroms to cm.
        
    # constants appropriate to cgs units.
    c1 =  3.7417749e-5 # 2*pi*h*c*c       
    c2 =  1.4387687 # h*c/k
    val =  c2/w/temp  
    
    bbflux =  c1/((w**5.)*(np.exp(val) - 1.))

    return bbflux*1e-8 # Convert to ergs/cm2/s/A

def dshift(wave, vrad):
    """ Given the wavelength and the radial velocity in km/s returns the
        doppler shifted wavelength.
    """
    
    wave = wave*(1. + vrad/sol_c)

    return wave

def lsf_func(vel, vr, vsini, epsilon=.6):
    """ Evaluate the rotational broadening function."""    
    
    # Limb darkening.
    e1 = 2.*(1. - epsilon)
    e2 = np.pi*epsilon/2.
    e3 = np.pi*(1. - epsilon/3.)
    
    # Compute the broadening function.
    x = (vel - vr)/vsini
    x1 = np.abs(1. - x**2)
    
    lsf = (e1*np.sqrt(x1) + e2*x1)/e3 
    lsf = np.where(np.abs(vel - vr) < vsini, lsf, 0.)    
    
    return lsf

def lsf_kernel(deltav, vsini, epsilon=.6):
    """ Returns the kernel for rotational broadening by vsini.""" 
    
    # Limb darkening.
    e1 = 2.*(1. - epsilon)
    e2 = np.pi*epsilon/2.
    e3 = np.pi*(1. - epsilon/3.)

    # Size of kernel, set to be uneven.
    npts = np.ceil(2*vsini/deltav).astype('int')
    if (npts%2 == 0):
        npts = npts + 1
       
    # Create velocity array.
    nwid = npts/2
    x = (np.arange(npts) - nwid)
    x = x*deltav/vsini  
    
    # Return the broadening kernel.
    x1 = np.abs(1. - x**2)
    
    return x*vsini, (e1*np.sqrt(x1) + e2*x1)/e3

def broaden(spec, deltav, vsini, epsilon=.6):
    """ Apply rotational broadening to a spectrum."""     
    
    velgrid, kernel = lsf_kernel(deltav, vsini, epsilon)
    spec = np.convolve(spec, kernel/kernel.sum(), mode='same')
    
    return spec

def vac2air(wave_vac):
    """ Convert wavelength in vacuum to wavelength in air."""
    
    sigma2 = (1e4/wave_vac)**2. # Convert to wavenumber squared

    # Compute conversion factor
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
    
    # Compute air wavelengths.
    wave_air = np.where(wave_vac >= 2000, wave_vac/fact, wave_vac) 

    return wave_air
    
def regular_velgrid(wave, spec, window, deltav, pad=200.):
    """ Put a spectrum on a wavelength grid with constant velocity spacing."""     
    
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
    
def regular_wavegrid(wave, spec, window, dw, pad=200.):
    
    # Compute the new wavelength grid.
    new_wave = np.linspace(window[0] - pad, window[1] + pad, (window[1] - window[0] + 2*pad)/dw)
    
    # Interpolate the spectrum to the new grid.
    new_spec = np.interp(new_wave, wave, spec)
    
    return new_wave, new_spec
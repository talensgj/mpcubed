#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob

import numpy as np
import pyfits as pf

import core

def index_hermes(directory, mode='wavelength_merged'):
    """ List all HERMES spectra in the specified location."""

    mode_values = ['log_merged', 'wavelength_merged']

    # Check that mode is set to a valid value.
    if mode not in mode_values:
        raise ValueError('mode must be one of {}'.format(mode_values))

    # List the spectra.
    filelist = glob.glob(directory + '*{}_c.fits'.format(mode))
    filelist = np.sort(filelist)    
    
    # Create arrays.
    nfiles = len(filelist)
    unseq = np.zeros(nfiles, dtype='int')
    bjd = np.zeros(nfiles)
    bvcor = np.zeros(nfiles)
    exptime = np.zeros(nfiles)
    Object = np.zeros(nfiles, dtype='S20')
    
    for i in range(nfiles):
        filename = filelist[i]        
        
        # Read the header information.
        header = pf.getheader(filename)
        unseq[i] = header['UNSEQ']
        bjd[i] = header['BJD']
        bvcor[i] = header['BVCOR']
        exptime[i] = header['EXPTIME']
        Object[i] = header['OBJECT']
    
    return filelist, unseq, bjd, bvcor, exptime, Object

def wave_hermes(header):
    """ Create the wavelength axis of a HERMES spectrum.""" 

    npoints = int(header['NAXIS1'])
    ref_pixc = int(header['CRPIX1'])
    ref_valc = float(header['CRVAL1'])
    ref_delc = float(header['CDELT1'])
    ref_type = header['CTYPE1']    
    
    wave = np.arange(npoints)*ref_delc + ref_valc - (ref_pixc - 1)*ref_delc
    
    if (ref_type == 'log(wavelength)'):
        wave = np.exp(wave)
    
    return wave

def read_hermes(filename):
    """ Read a HERMES spectrum."""    
    
    spec, header = pf.getdata(filename, header=True)
    wave = wave_hermes(header)

    return wave, spec, header
    
def read_phoenix(filename, air=True, wavefile='/net/boschwijde/data2/talens/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'):
    """ Read a PHOENIX spectrum."""    
    
    spec, header = pf.getdata(filename, header=True)
    wave = pf.getdata(wavefile) 
    
    if air:
        wave = core.vac2air(wave)    
    
    return wave, spec, header
    
def read_kurucz(filename):
    """ Read a Kurucz spectrum."""    
    
    data, header = pf.getdata(filename, header=True)
    wave = data[0]
    spec = data[1]
    
    return wave, spec, header
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from fLCfile import fLCfile

class Data():

    def __init__()

        f = fLCfile('/data2/talens/3mEast/June1.hdf5')
        
        ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        nobs = Mnobs.astype('int')

        staridx = np.arange(len(ascc))
        decidx = pg.find_decidx(dec)
        skyidx = hg.find_gridpoint(ra, dec)
        
        return staridx, decidx, skyidx

    def get_decidx():
        return self.decidx
        
    def get_camdata(idx):
        
        here = (self.decidx == idx)
        ascc = self.ascc[here]
        nobs = self.nobs[here]
        staridx = self.staridx[here]
        skyidx = self.skyidx[here]
        
        jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
        lstidx = lstidx.astype('int')

        # Build indices.    
        staridx = np.repeat(staridx, nobs)
        
        ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
        camtransidx = pg1.find_gridpoint(ha, np.repeat(dec, nobs))
        intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))
        
        skyidx = np.repeat(skyidx, nobs)
        
        # To be replaced with lstseq?
        dayidx = np.floor(jdmid).astype('int')
        #dayidx = dayidx - 2457175 #June1
        dayidx = dayidx - 2457190
        skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
            
        # Flag bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
        flux = flux[here]
        eflux = eflux[here]
        x = x[here]
        y = y[here]

        ## What if there is no data?
        #if len(flux) == 0: continue

        staridx = staridx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        skyidx = skyidx[here]
        skytransidx = skytransidx[here]

        # Convert flux to magnitudes
        mag = 25 - 2.5*np.log10(flux)
        emag = 2.5/np.log(10)*eflux/flux
        
        return mag, emag, x, y, staridx, camtransidx, intrapixidx, skyidx, skytransidx

def CameraTransmission():
    
    ## Function call returns header and related indices.
    
    # Determine the number of declination bins.
    decidx, decuni = np.unique(decidx, return_inverse=True)
    nbins = len(decidx)
    
    ## Function call returns the sky parameters.
    
    # Create arrays for the results.
    niter = np.zeros(nbins)
    chisq = np.zeros(nbins)
    npoints = np.zeros(nbins)
    npars = np.zeros(nbins)
    
    m = np.full(len(ascc), fill_value=np.nan)
    z = np.full((13502*722,), fill_value=np.nan)
    a = np.full((272*722,), fill_value=np.nan)
    b = np.full((272*722,), fill_value=np.nan)
    c = np.full((272*722,), fill_value=np.nan)
    d = np.full((272*722,), fill_value=np.nan)
    
    # Iterate over the declination bins.
    for idx in range(nbins):
        
        ## Function call returns data and related indices for the bin.
        
        if (skypars == True):
            mag = mag - s[skyidx, skytransidx]
            emag = np.sqrt(emag**2 + (sigma1[staridx])**2 + (sigma2[skyidx, skytransidx])**2)
        
        # Make the indices unique.
        staridx, staruni = np.unique(staridx, return_inverse=True)
        camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
        intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)
        
        # Compute the best fit transmission.
        m[staridx], z[camtransidx], a[intrapixidx], b[intrapixidx], c[intrapixidx], d[intrapixidx], niter[idx], chisq[idx], npoints[idx], npars[idx] = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag, emag, x, y, verbose=True, use_weights=False)
    
        # Make a normalization attempt.
        offset = np.nanmedian(m[staridx] - vmag[staridx])
        m[staridx] = m[staridx] - offset
        z[camtransidx] = z[camtransidx] + offset
    
    ## Function call to write results.
    
def SkyTransmission():
    
    ## Function call returns header and related indices.
    
    # Determine the number of sky bins.
    skyidx, skyuni = np.unique(skyidx, return_inverse=True)
    nbins = len(skyidx)
    
    ## Function call returns the camera parameters.
    
    # Create arrays for the results.
    niter = np.zeros(nbins)
    chisq = np.zeros(nbins)
    npoints = np.zeros(nbins)
    npars = np.zeros(nbins)
    
    m = np.full(len(ascc), fill_value=np.nan)
    sigma1 = np.full(len(ascc), fill_value=np.nan)
    s = np.full((hg.npix, 15*13500), fill_value=np.nan)
    sigma2 = np.full((hg.npix, 15*13500), fill_value=np.nan)
    
    # Iterate over the declination bins.
    for idx in range(nbins):
        
        ## Function call returns data and related indices for the bin.
        
        if (campars == True):
            mag = mag - z[camtransidx] - a[intrapixidx]*np.sin(2*np.pi*x) - b[intrapixidx]*np.cos(2*np.pi*x) - c[intrapixidx]*np.sin(2*np.pi*y) - d[intrapixidx]*np.cos(2*np.pi*y)
        
        # Make the indices unique.
        staridx, staruni = np.unique(staridx, return_inverse=True)
        skytransidx, skytransuni = np.unique(skytransidx, return_inverse=True)
        
        # Compute the best fit transmission.
        m[staridx], sigma1[staridx], s[skyidx, skytransidx], sigma2[skyidx, skytransidx], niter[idx], chisq[idx], npoints[idx], npars[idx] = systematics_dev.trans(staruni, skytransuni, mag, emag, verbose=True, use_weights=True)
    
        # Make a normalization attempt.
        offset = np.nanmedian(m[staridx] - vmag[staridx])
        m[staridx] = m[staridx] - offset
        s[skyidx, skytransidx] = s[skyidx, skytransidx] + offset
    
    ## Function call to write results.


def CoarseDecorrelation():

    for i in range(5):
        
        CameraTransmission()
        SkyTransmission()

if __name__ == '__main__':
    

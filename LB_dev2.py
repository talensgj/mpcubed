#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from fLCfile import fLCfile

class CoarseDecorrelation():
    
    def __init__(self, LBfile):
        
        self.LBfile = LBfile
        
        self.camgrid = 'polar'
        self.camnx = 13500
        self.camny = 720
        
        self.ipxgrid = 'polar'
        self.ipxnx = 270
        self.ipxny = 720
        
        self.skygrid = 'healpix'
        self.skynx = 8
    
    def spatial(self):
        
        ascc = self.ascc
        nobs = self.nobs
        decidx = self.decidx
        
        decidx, decuni = np.unique(decidx, return_inverse=True)
        
        nbins = len(decidx)
        for idx in range(nbins):
            
            # Read data.
            here = (decuni == idx)
            flux, eflux, sky, x, y, lst, lstseq = self.f.read_data(['flux0', 'eflux0', 'sky', 'x', 'y', 'lst', 'lstseq'], ascc[here], nobs[here])
            
            ha = np.mod(lst*15. - ra, 360.)
            
            # Create indices.
            staridx = self.staridx[here]
            camtransidx = self.camgrid.find_gridpoint(ha, dec)
            intrapixidx = self.ipxgrid.find_gridpoint(ha, dec)
            skyidx = self.skyidx[here]
            
            # Remove bad data.
            here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
            
            flux = flux[here]
            eflux = eflux[here]
            x = x[here]
            y = y[here]
            lst = lst[here]
            lstseq = lstseq[here]
            
            staridx = staridx[here]
            camtransidx = camtransidx[here]
            intrapixidx = intrapixidx[here]
            skyidx = skyidx[here]
            
            # Convert flux to magnitudes:
            mag, emag = flux2mag(flux, eflux)
            
            # Apply known temporal correction.
            mag = mag - s[skyidx, lstseq]
            emag = np.sqrt(emag**2 + sigma1[staridx]**2 + sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            staridx, ind1 = np.unique(staridx, return_inverse=True)
            camtransidx, ind2 = np.unique(camtransidx, return_inverse=True)
            intrapixidx, ind3 = np.unique(intrapixidx, return_inverse=True)
            
            # Calculate new spatial correction.
            m[staridx], z[camtransidx], A[intrapixidx], niter[idx], chisq[idx], npoints[idx], npars[idx] = trans_ipx(ind1, ind2, ind3, mag, emag, x, y)
            
        return
            
    def temporal(self):
            
        skyidx, skyuni = np.unique(skyidx, return_inverse=True)
        
        nbins = len(skyidx)
        for idx in range(nbins):
            
            # Read data.
            here = (skyuni == idx)
            flux, eflux, sky, x, y, lst, lstseq = self.f.read_data(['flux0', 'eflux0', 'sky', 'x', 'y', 'lst', 'lstseq'], ascc[here], nobs[here])
            
            ha = np.mod(lst*15. - ra, 360.)
            
            # Create indices.
            staridx = self.staridx[here]
            camtransidx = self.camgrid.find_gridpoint(ha, dec)
            intrapixidx = self.ipxgrid.find_gridpoint(ha, dec)
            skyidx = self.skyidx[here]
            
            # Remove bad data.
            here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
            
            flux = flux[here]
            eflux = eflux[here]
            x = x[here]
            y = y[here]
            lst = lst[here]
            lstseq = lstseq[here]
            
            staridx = staridx[here]
            camtransidx = camtransidx[here]
            intrapixidx = intrapixidx[here]
            skyidx = skyidx[here]
            
            # Convert flux to magnitudes:
            mag, emag = flux2mag(flux, eflux)
            
            # Apply known spatial correction.
            mag = mag - z[camtransidx]
            mag = mag - a[intrapixidx]*np.sin(2*np.pi*x) - b[intrapixidx]*np.cos(2*np.pi*x) - c[intrapixidx]*np.sin(2*np.pi*y) - d[intrapixidx]*np.cos(2*np.pi*y)
            
            # Create unique indices.
            staridx, ind1 = np.unique(staridx, return_inverse=True)
            lstseq, ind2 = np.unique(lstseq, return_inverse=True)
            
            # Calculate new temporal correction.
            m[staridx], sigma1[staridx], s[skyidx, lstseq], sigma2[skyidx, lstseq], niter[idx], chisq[idx], npoints[idx], npars[idx] = sky(ind1, ind2, mag, emag, use_weights=True)
            
        return
    
    def calculate():

        
        

if __name__ == '__main__':
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coordinate_grids import PolarGrid, HealpixGrid
from core.LB_decor import spatial_decor, temporal_decor
from usefull_functions_dev import flux2mag

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
    
    def read_data(self, here):
        
        ascc = self.ascc[here]
        ra = self.ra[here]
        dec = self.dec[here]
        nobs = self.nobs[here]
        
        staridx = self.staridx[here]
        skyidx = self.skyidx[here]
        
        # Read data.
        flux, eflux, sky, x, y, lst, lstseq, flags = self.f.read_data(['flux0', 'eflux0', 'sky', 'x', 'y', 'lst', 'lstseq', 'flag'], ascc, nobs)
        lstseq = lstseq.astype('int') - self.lstmin
        
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = np.mod(lst*15. - ra, 360.)
        
        # Create indices.
        staridx = np.repeat(staridx, nobs)
        camtransidx = self.camgrid.find_gridpoint(ha, dec)
        intrapixidx = self.ipxgrid.find_gridpoint(ha, dec)
        skyidx = np.repeat(skyidx, nobs)
        
        # Remove bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
        
        flux = flux[here]
        eflux = eflux[here]
        x = x[here]
        y = y[here]
        lstseq = lstseq[here]
        
        staridx = staridx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        skyidx = skyidx[here]
        
        # Convert flux to magnitudes:
        mag, emag = flux2mag(flux, eflux)
        
        return mag, emag, x, y, staridx, camtransidx, intrapixidx, skyidx, lstseq
    
    def spatial(self):
        
        decidx, decuni = np.unique(self.decidx, return_inverse=True)
        
        nbins = len(decidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        for idx in range(nbins):
            
            # Read data.
            here = (decuni == idx)
            mag, emag, x, y, staridx, camtransidx, intrapixidx, skyidx, lstseq = self.read_data(here)
            
            # Apply known temporal correction.
            if self.got_sky == True:
                mag = mag - self.s[skyidx, lstseq]
                emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + self.sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            camtransidx, ind2, z_nobs = np.unique(camtransidx, return_inverse=True, return_counts=True)
            intrapixidx, ind3, A_nobs = np.unique(intrapixidx, return_inverse=True, return_counts=True)
            
            # Calculate new spatial correction.
            m, z, A, niter[idx], chisq[idx], npoints[idx], npars[idx] = spatial_decor(ind1, ind2, ind3, mag, emag, x, y)
            
            self.m[staridx] = m
            self.z[camtransidx] = z
            self.A[intrapixidx] = A
            
            self.m_nobs[staridx] = m_nobs
            self.z_nobs[camtransidx] = z_nobs
            self.A_nobs[intrapixidx] = A_nobs
            
        return
            
    def temporal(self):
            
        skyidx, skyuni = np.unique(self.skyidx, return_inverse=True)
        
        nbins = len(skyidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        for idx in range(nbins):
            
            # Read data.
            here = (skyuni == idx)
            mag, emag, x, y, staridx, camtransidx, intrapixidx, _, lstseq = self.read_data(here)
            
            # Apply known spatial correction.
            mag = mag - self.z[camtransidx]
            mag = mag - np.sum(self.A[intrapixidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            lstseq, ind2, s_nobs = np.unique(lstseq, return_inverse=True, return_counts=True)
            
            # Calculate new temporal correction.
            m, s, sigma1, sigma2, niter[idx], chisq[idx], npoints[idx], npars[idx] = temporal_decor(ind1, ind2, mag, emag, use_weights=True)
            
            self.m[staridx] = m
            self.s[skyidx[idx], lstseq] = s
            
            self.sigma1[staridx] = sigma1
            self.sigma2[skyidx[idx], lstseq] = sigma2
            
            self.m_nobs[staridx] = m_nobs
            self.s_nobs[skyidx[idx], lstseq] = s_nobs
           
        self.got_sky = True
            
        return
    
    def calculate(self):

        self.lstmin = 12163500
        lstlen = 27000

        self.f = fLCfile(self.LBfile)
        self.camgrid = PolarGrid(self.camnx, self.camny)
        self.ipxgrid = PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = HealpixGrid(self.skynx)
        
        self.ascc, self.ra, self.dec, self.nobs, self.vmag = self.f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.nobs = self.nobs.astype('int')
        
        self.staridx = np.arange(len(self.ascc))
        self.decidx = self.camgrid.find_decidx(self.dec)
        self.skyidx = self.skygrid.find_gridpoint(self.ra, self.dec)
        
        self.m = np.full(len(self.ascc), fill_value=np.nan)
        self.z = np.full(self.camgrid.npix, fill_value=np.nan)
        self.A = np.full((self.ipxgrid.npix, 4), fill_value=np.nan)
        self.s = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.sigma1 = np.full(len(self.ascc), fill_value=np.nan)
        self.sigma2 = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.m_nobs = np.full(len(self.ascc), fill_value=np.nan)
        self.z_nobs = np.full(self.camgrid.npix, fill_value=np.nan)
        self.A_nobs = np.full(self.ipxgrid.npix, fill_value=np.nan)
        self.s_nobs = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.got_sky = False
        
        import matplotlib.pyplot as plt
        
        for niter in range(5):
        
            self.spatial()
            self.temporal()
            
            #plt.plot(self.vmag, self.m, '.')
            #plt.show()
            
            plt.imshow(self.z.reshape(13502, 722).T, interpolation='None', aspect='auto', vmin = np.nanpercentile(self.z, 1), vmax = np.nanpercentile(self.z, 99))
            plt.colorbar()
            plt.show()
            plt.close()
            
            #plt.imshow(self.z_nobs.reshape(13502, 722).T, interpolation='None', aspect='auto')
            #plt.colorbar()
            #plt.show()
            
            plt.imshow(self.A[:,0].reshape(272, 722).T, interpolation='None', aspect='auto', vmin = np.nanpercentile(self.A[:,0], 1), vmax = np.nanpercentile(self.A[:,0], 99))
            plt.colorbar()
            plt.show()
            plt.close()
            
            #plt.imshow(self.A_nobs.reshape(272, 722).T, interpolation='None', aspect='auto')
            #plt.colorbar()
            #plt.show()
            
            plt.imshow(self.s, interpolation='None', aspect='auto', vmin = np.nanpercentile(self.s, 1), vmax = np.nanpercentile(self.s, 99))
            plt.colorbar()
            plt.show()
            plt.close()
        
if __name__ == '__main__':
    
    obj = CoarseDecorrelation('/data2/talens/Orientation/fLC_20150618LPE.hdf5')
    obj.calculate()
    

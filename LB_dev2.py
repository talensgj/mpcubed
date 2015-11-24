#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coordinate_grids import PolarGrid, HealpixGrid
from core.LB_decor import spatial_decor, temporal_decor, spatial_decor_new
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
            
            mag = mag - self.m[staridx]
            # Apply known temporal correction.
            if self.got_sky == True:
                mag = mag - self.s[skyidx, lstseq]
                #emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + self.sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            #staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            camtransidx, ind2, z_nobs = np.unique(camtransidx, return_inverse=True, return_counts=True)
            intrapixidx, ind3, A_nobs = np.unique(intrapixidx, return_inverse=True, return_counts=True)
            
            # Calculate new spatial correction.
            #m, z, A, niter[idx], chisq[idx], npoints[idx], npars[idx] = spatial_decor(ind1, ind2, ind3, mag, emag, x, y)
            z, A, niter[idx], chisq[idx], npoints[idx], npars[idx] = spatial_decor_new(ind2, ind3, mag, emag, x, y)
            #offset = np.nanmedian(m - self.vmag[staridx])
            #m = m - offset
            #z = z + offset
            
            #self.m[staridx] = m
            self.z[camtransidx] = z
            self.A[intrapixidx] = A
            
            #self.m_nobs[staridx] = m_nobs
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
            #m, s, sigma1, sigma2, niter[idx], chisq[idx], npoints[idx], npars[idx] = temporal_decor(ind1, ind2, mag, emag, use_weights=True)
            m, s, niter[idx], chisq[idx], npoints[idx], npars[idx] = temporal_decor(ind1, ind2, mag, emag, use_weights=False)
            
            #offset = np.nanmedian(m - self.vmag[staridx])
            #m = m - offset
            #s = s + offset
            
            self.m[staridx] = m
            self.s[skyidx[idx], lstseq] = s
            
            #self.sigma1[staridx] = sigma1
            #self.sigma2[skyidx[idx], lstseq] = sigma2
            
            self.m_nobs[staridx] = m_nobs
            self.s_nobs[skyidx[idx], lstseq] = s_nobs
           
        self.got_sky = True
            
        return
    
    def calculate(self):

        self.f = fLCfile(self.LBfile)
        self.camgrid = PolarGrid(self.camnx, self.camny)
        self.ipxgrid = PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = HealpixGrid(self.skynx)
        
        self.ascc, self.ra, self.dec, self.nobs, self.vmag = self.f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.nobs = self.nobs.astype('int')
        
        self.staridx = np.arange(len(self.ascc))
        self.decidx = self.camgrid.find_decidx(self.dec)
        self.skyidx = self.skygrid.find_gridpoint(self.ra, self.dec)
        
        #lstseq, = self.f.read_data(['lstseq'], self.ascc, self.nobs) # May be too memory intensive?
        #lstseq = lstseq.astype('int')
        
        #self.lstmin = np.amin(lstseq)
        #self.lstlen = np.amax(lstseq) - np.amin(lstseq) + 1
        
        with h5py.File(self.LBfile, 'r') as f:
            lstmin = f['data'].attrs['lstmin']
            lstmax = f['data'].attrs['lstmax']
            
        self.lstmin = lstmin
        self.lstlen = lstmax - lstmin + 1
        
        #self.m = np.full(len(self.ascc), fill_value=np.nan)
        self.m = np.copy(self.vmag)
        self.z = np.full(self.camgrid.npix, fill_value=np.nan)
        self.A = np.full((self.ipxgrid.npix, 4), fill_value=np.nan)
        self.s = np.full((self.skygrid.npix, self.lstlen), fill_value=np.nan)
        
        self.sigma1 = np.full(len(self.ascc), fill_value=np.nan)
        self.sigma2 = np.full((self.skygrid.npix, self.lstlen), fill_value=np.nan)
        
        self.m_nobs = np.full(len(self.ascc), fill_value=np.nan)
        self.z_nobs = np.full(self.camgrid.npix, fill_value=np.nan)
        self.A_nobs = np.full(self.ipxgrid.npix, fill_value=np.nan)
        self.s_nobs = np.full((self.skygrid.npix, self.lstlen), fill_value=np.nan)
        
        self.got_sky = False
        
        for niter in range(2):
        
            self.spatial()
            self.temporal()
            
        with h5py.File('/data2/talens/3mEast/LBtests/June2_profiled.hdf5') as f:
    
            #hdr = f.create_group('header')
            #hdr.create_dataset('decidx', data=decidx)
            #hdr.create_dataset('niter', data=niter)
            #hdr.create_dataset('chisq', data=chisq)
            #hdr.create_dataset('npoints', data=npoints)
            #hdr.create_dataset('npars', data=npars)
    
            #hdr = f.create_group('header')
            #hdr.create_dataset('skyidx', data=skyidx)
            #hdr.create_dataset('niter', data=niter) 
            #hdr.create_dataset('chisq', data=chisq)
            #hdr.create_dataset('npoints', data=npoints)
            #hdr.create_dataset('npars', data=npars)
            
            grp = f.create_group('data')
            
            grp.create_dataset('magnitudes/ascc', data=self.ascc)
            grp.create_dataset('magnitudes/m', data=self.m)
            grp.create_dataset('magnitudes/sigma', data=self.sigma1)
            grp.create_dataset('magnitudes/nobs', data=self.m_nobs)
            
            idx, = np.where(~np.isnan(self.z))
            grp.create_dataset('camtrans/idx', data=idx)
            grp.create_dataset('camtrans/z', data=self.z[idx])
            grp.create_dataset('camtrans/nobs', data=self.z_nobs[idx])
            
            grp['camtrans'].attrs['grid'] = 'polar'
            grp['camtrans'].attrs['nx'] = self.camnx
            grp['camtrans'].attrs['ny'] = self.camny
            
            idx, = np.where(~np.isnan(self.A[:,0]))
            grp.create_dataset('intrapix/idx', data=idx)
            grp.create_dataset('intrapix/a', data=self.A[idx,0])
            grp.create_dataset('intrapix/b', data=self.A[idx,1])
            grp.create_dataset('intrapix/c', data=self.A[idx,2])
            grp.create_dataset('intrapix/d', data=self.A[idx,3])
            grp.create_dataset('intrapix/nobs', data=self.A_nobs[idx])
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            idx, lstseq = np.where(~np.isnan(self.s))
            grp.create_dataset('skytrans/idx', data=idx)
            grp.create_dataset('skytrans/lstseq', data=lstseq)
            grp.create_dataset('skytrans/s', data=self.s[idx, lstseq])
            grp.create_dataset('skytrans/sigma', data=self.sigma2[idx, lstseq])
            grp.create_dataset('skytrans/nobs', data=self.s_nobs[idx, lstseq])
            
            grp['skytrans'].attrs['grid'] = 'healpix'
            grp['skytrans'].attrs['nx'] = self.skynx
            grp['skytrans'].attrs['lstmin'] = self.lstmin
            grp['skytrans'].attrs['lstlen'] = self.lstlen
            
if __name__ == '__main__':
    
    obj = CoarseDecorrelation('/data2/talens/3mEast/LBtests/June2_fLC_auto.hdf5')
    obj.calculate()
    

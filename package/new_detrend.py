#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from progressbar import ProgressBar

import misc
import IO
from coordinates import grids
from core import cdecor
from core import statistics

def _spatial(self):
    """Solve for the transmission and intrapixel variations."""
    
    decidx, decuni = np.unique(self.decidx, return_inverse=True)
    
    nbins = len(decidx)
    for idx in range(nbins):
        
        # Read data.
        here = (decuni == idx)
        mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self.read_data(here) # Unclear function call + plu how to to split it over multiple lines?
        
        #mag = data[0] if I first capure output as single entity.
        #emag = data.emag if I make output named tuple.
        
        if (len(mag) == 0): continue
        
        # Apply known temporal correction.
        if self.got_sky:
            mag = mag - self.s[skyidx, lstseq]
        
        if self.got_sky & self.sigmas:
            emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + 
                           self.sigma2[skyidx, lstseq]**2)
        
        # Create unique indices.
        staridx, ind1 = np.unique(staridx, return_inverse=True)
        camtransidx, ind2 = np.unique(camtransidx, return_inverse=True)
        intrapixidx, ind3 = np.unique(intrapixidx, return_inverse=True)
        
        # Calculate new spatial correction.
        m, z, A, quality = cdecor.cdecor_intrapix(ind1, ind2, ind3, mag, emag, x, y, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
        
        # Simple magnitude calibration.
        offset = np.nanmedian(m - self.vmag[staridx])
        
        # Store results.
        self.m[staridx] = m - offset
        self.z[camtransidx, decidx[idx]] = z + offset
        self.A[intrapixidx, decidx[idx]] = A
        
        self.spatial_niter[idx] = quality.niter #switched name around _
        self.spatial_chisq[idx] = quality.chisq
        self.spatial_npoints[idx] = quality.npoints
        self.spatial_npars[idx] = quality.npars 
        
        # The count might be moved into the reader, make it easier to do it only once.
        self.nobs_m[staridx] = np.bincount(ind1) #switched name around _
        self.nobs_z[camtransidx, decidx[idx]] = np.bincount(ind2)
        self.nobs_A[intrapixidx, decidx[idx]] = np.bincount(ind3)
        
    return
    

def coarse_decorrelation(LBfile, aperture, sysfile = None, **kwargs):
    
    if sysfile is None:
        head, tail = os.path.split(self.LBfile)
        prefix = 'sys%i_'%self.aper
        tail = prefix + tail.rsplit('_')[-1]
        sysfile = os.path.join(head, tail)
        
    # Initialize with defualt parameters unless arguments were given.
    sigmas = kwargs.pop('sigmas', True)
    outer_maxiter = kwargs.pop('outer_maxiter', 5)
    
    camgrid = 'polar'
    camnx = kwargs.pop('camnx', 13500)
    camny = kwargs.pop('camny', 720)
    
    ipxgrid = 'polar'
    ipxnx = kwargs.pop('ipxnx', 270)
    ipxny = self.camny
    
    skygrid = 'healpix'
    skynx = kwargs.pop('skynx', 8)

    # Options to coarse_decor functions.
    inner_maxiter = kwargs.pop('inner_maxiter', 100)
    dtol = kwargs.pop('dtol', 1e-3)
    verbose = kwargs.pop('verbose', False)
    
    # Perform the coarse decorrelation.
    print 'Performing coarse decorrelation for:', LBfile
    for niter in range(outer_maxiter):
    
        print 'Iteration %i out of %i:'%(niter + 1, outer_maxiter)
        
        print '    Calculating camera systematics...'
        _spatial(LBfile, aperture, sysfile)
        
        print '    Calculating atmospheric systematics...'
        _temporal(LBfile, aperture, sysfile)
    
    return sysfile
    
class CoarseDecorrelation():
    
    def __init__(self, LBfile, aperture, sysfile = None, **kwargs):
        """
            Performs a coarse decorrelation on all data in a given Long Baseline
            fLC file. Removes camera transmission, intrapixel variations and
            sky transmission.
        """
        
        # fLC file and aperture to work on.
        self.LBfile = LBfile
        self.aper = aperture
        
        if not os.path.isfile(self.LBfile):
            print 'File not found:', self.LBfile
            print 'exiting...'
            exit()
        else:
            print 'Calculating corrections for aperture %i of file:'%self.aper, self.LBfile
        
        # The systematics file.
        if sysfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'sys%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            sysfile = os.path.join(head, tail)
        
        self.sysfile = sysfile
        
        if os.path.isfile(self.sysfile):
            print 'Systematics file already exists:', self.sysfile
            print 'exiting...'
            exit()
        else:
            print 'Writing results to:', self.sysfile
        
        # Initialize with defualt parameters unless arguments were given.
        self.sigmas = kwargs.pop('sigmas', True)
        self.outer_maxiter = kwargs.pop('outer_maxiter', 5)
        
        self.camgrid = 'polar'
        self.camnx = kwargs.pop('camnx', 13500)
        self.camny = kwargs.pop('camny', 720)
        
        self.ipxgrid = 'polar'
        self.ipxnx = kwargs.pop('ipxnx', 270)
        self.ipxny = self.camny
        
        self.skygrid = 'healpix'
        self.skynx = kwargs.pop('skynx', 8)
    
        # Options to coarse_decor functions.
        self.inner_maxiter = kwargs.pop('inner_maxiter', 100)
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.verbose = kwargs.pop('verbose', False)
    
        return
    
    def read_data(self, here):
        """
            Reads portions of the data from file, creates appropriate indices,
            removes flagged data and converts fluxes to magnitudes.
        """
        
        ascc = self.ascc[here]
        ra = self.ra[here]
        dec = self.dec[here]
        nobs = self.nobs[here]
        
        staridx = self.staridx[here]
        skyidx = self.skyidx[here]
        
        # Read data.
        fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'lst', 'lstseq', 'flag']
        flux, eflux, sky, x, y, lst, lstseq, flags = self.f.read_data(fields, ascc, nobs)
        lstseq = lstseq.astype('int') - self.lstmin
        
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = np.mod(lst*15. - ra, 360.)
        
        # Create indices.
        staridx = np.repeat(staridx, nobs)
        camtransidx, decidx = self.camgrid.radec2idx(ha, dec)
        intrapixidx, decidx = self.ipxgrid.radec2idx(ha, dec)
        skyidx = np.repeat(skyidx, nobs)
        
        # Remove bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
        
        flux = flux[here]
        eflux = eflux[here]
        x = x[here]
        y = y[here]
        lstseq = lstseq[here]
        
        staridx = staridx[here]
        decidx = decidx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        skyidx = skyidx[here]
        
        # Convert flux to magnitudes:
        mag, emag = misc.flux2mag(flux, eflux)
        
        return mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, skyidx, lstseq
    
    def spatial(self):
        """
            Solves for the camera transmission and intrapixel variations which
            are assumed to be time independent.
        """
        
        decidx, decuni = np.unique(self.decidx, return_inverse=True)
        
        nbins = len(decidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        
        pbar = ProgressBar(maxval = nbins).start()
        for idx in range(nbins):
            
            # Read data.
            here = (decuni == idx)
            mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self.read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known temporal correction.
            if self.got_sky == True:
                mag = mag - self.s[skyidx, lstseq]
                
                if self.sigmas:
                    emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + self.sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            camtransidx, ind2, z_nobs = np.unique(camtransidx, return_inverse=True, return_counts=True)
            intrapixidx, ind3, A_nobs = np.unique(intrapixidx, return_inverse=True, return_counts=True)
            
            # Calculate new spatial correction.
            m, z, A, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor_intrapix(ind1, ind2, ind3, mag, emag, x, y, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Simple magnitude calibration.
            offset = np.nanmedian(m - self.vmag[staridx])
            
            self.m[staridx] = m - offset
            self.z[camtransidx, decidx[idx]] = z + offset
            self.A[intrapixidx, decidx[idx]] = A
            
            self.m_nobs[staridx] = m_nobs
            self.z_nobs[camtransidx, decidx[idx]] = z_nobs
            self.A_nobs[intrapixidx, decidx[idx]] = A_nobs
            
            pbar.update(idx + 1)
            
        pbar.finish()
            
        self.niter_spatial = niter
        self.chisq_spatial = chisq
        self.npoints_spatial = npoints
        self.npars_spatial = npars 
            
        return
        
    def temporal(self):
        """
            Solves for the sky transmission as a function of position and time.
        """
        
        skyidx, skyuni = np.unique(self.skyidx, return_inverse=True)
        
        nbins = len(skyidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        
        pbar = ProgressBar(maxval = nbins).start()
        for idx in range(nbins):
            
            # Read data.
            here = (skyuni == idx)
            mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, _, lstseq = self.read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known spatial correction.
            mag = mag - self.z[camtransidx, decidx]
            mag = mag - np.sum(self.A[intrapixidx, decidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            lstseq, ind2, s_nobs = np.unique(lstseq, return_inverse=True, return_counts=True)
            
            # Calculate new temporal correction.
            if self.sigmas:
                m, s, sigma1, sigma2, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor_sigmas(ind1, ind2, mag, emag, self.sigma1[staridx], self.sigma2[skyidx[idx], lstseq], maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            else:
                m, s, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor(ind1, ind2, mag, emag, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Simple magnitude calibration.
            offset = np.nanmedian(m - self.vmag[staridx])
            
            self.m[staridx] = m - offset
            self.s[skyidx[idx], lstseq] = s + offset
            
            if self.sigmas:
                self.sigma1[staridx] = sigma1
                self.sigma2[skyidx[idx], lstseq] = sigma2
            
            self.m_nobs[staridx] = m_nobs
            self.s_nobs[skyidx[idx], lstseq] = s_nobs
            
            pbar.update(idx + 1)
           
        pbar.finish()
           
        self.niter_temporal = niter
        self.chisq_temporal = chisq
        self.npoints_temporal = npoints
        self.npars_temporal = npars 
           
        self.got_sky = True
            
        return
    
    def run(self):
        """
            Performs the coarse decorrelation on the given Long Baseline fLC file
            and writes the result to the given output location.
        """

        # Set up the IO and coordinate grids.
        self.f = IO.fLCfile(self.LBfile)
        self.camgrid = grids.PolarGrid(self.camnx, self.camny)
        self.ipxgrid = grids.PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = grids.HealpixGrid(self.skynx)
        
        # Read the required header data.
        self.ascc, self.ra, self.dec, self.nobs, self.vmag = self.f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.nobs = self.nobs.astype('int')
        
        # Create indices.
        self.staridx = np.arange(len(self.ascc))
        _, self.decidx = self.camgrid.radec2idx(self.ra, self.dec)
        self.skyidx = self.skygrid.radec2idx(self.ra, self.dec)
        
        # Read global information.
        with h5py.File(self.LBfile, 'r') as f:
            
            station = f['global'].attrs['station']
            camera = f['global'].attrs['camera']
            
            alt0 = f['global/alt0'].value
            az0 = f['global/az0'].value
            th0 = f['global/th0'].value
            x0 = f['global/x0'].value
            y0 = f['global/y0'].value
            
            lstmin = f['global'].attrs['lstmin'].astype('int')
            lstmax = f['global'].attrs['lstmax'].astype('int')
        
        self.lstmin = lstmin
        lstlen = lstmax - lstmin + 1
        
        # Create arrays to hold the results.
        self.m = np.full(len(self.ascc), fill_value=np.nan)
        self.z = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.A = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan)
        self.s = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        if self.sigmas:
            self.sigma1 = np.full(len(self.ascc), fill_value=0)
            self.sigma2 = np.full((self.skygrid.npix, lstlen), fill_value=0)
        
        self.m_nobs = np.full(len(self.ascc), fill_value=np.nan)
        self.z_nobs = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.A_nobs = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2), fill_value=np.nan)
        self.s_nobs = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.got_sky = False
        
        # Perform the coarse decorrelation.
        print 'Performing coarse decorrelation for:', self.LBfile
        for niter in range(self.outer_maxiter):
        
            print 'Iteration %i out of %i:'%(niter + 1, self.outer_maxiter)
            
            print '    Calculating camera systematics...'
            self.spatial()
            
            print '    Calculating atmospheric systematics...'
            self.temporal()
        
        # Write the results to file.
        with h5py.File(self.sysfile) as f:
    
            # Write the header.
            hdr = f.create_group('header')
            
            hdr.attrs['station'] = station
            hdr.attrs['camera'] = camera
            
            hdr.attrs['alt0'] = np.mean(alt0)
            hdr.attrs['az0'] = np.mean(az0)
            hdr.attrs['th0'] = np.mean(th0)
            hdr.attrs['x0'] = np.mean(x0)
            hdr.attrs['y0'] = np.mean(y0)
            
            hdr.attrs['outer_maxiter'] = self.outer_maxiter
            hdr.attrs['inner_maxiter'] = self.inner_maxiter
            hdr.attrs['sigmas'] = self.sigmas
            hdr.attrs['dtol'] = self.dtol
           
            hdr.create_dataset('spatial/niter', data = self.niter_spatial, dtype = 'uint32')
            hdr.create_dataset('spatial/chisq', data = self.chisq_spatial, dtype = 'float64')
            hdr.create_dataset('spatial/npoints', data = self.npoints_spatial, dtype = 'uint32')
            hdr.create_dataset('spatial/npars', data = self.npars_spatial, dtype = 'uint32')
            
            hdr.create_dataset('temporal/niter', data = self.niter_temporal, dtype = 'uint32')
            hdr.create_dataset('temporal/chisq', data = self.chisq_temporal, dtype = 'float64')
            hdr.create_dataset('temporal/npoints', data = self.npoints_temporal, dtype = 'uint32')
            hdr.create_dataset('temporal/npars', data = self.npars_temporal, dtype = 'uint32')
            
            # Write the data.
            grp = f.create_group('data')
            
            # Write the magnitudes.
            grp.create_dataset('magnitudes/ascc', data = self.ascc)
            grp.create_dataset('magnitudes/vmag', data = self.vmag, dtype = 'float32')
            grp.create_dataset('magnitudes/nobs', data = self.m_nobs, dtype = 'uint32')
            grp.create_dataset('magnitudes/mag', data = self.m, dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('magnitudes/sigma', data = self.sigma1, dtype = 'float32')
            
            # Write the camera transmission.
            self.z = np.ravel(self.z)
            idx, = np.where(~np.isnan(self.z))
            grp.create_dataset('trans/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('trans/nobs', data = self.z_nobs[idx], dtype = 'uint32')
            grp.create_dataset('trans/trans', data = self.z[idx], dtype = 'float32')
            
            grp['trans'].attrs['grid'] = 'polar'
            grp['trans'].attrs['nx'] = self.camnx
            grp['trans'].attrs['ny'] = self.camny
            
            # Write the intrapixel variations.
            self.A = np.reshape(self.A, ((self.ipxgrid.nx+2)*(self.ipxgrid.ny+2), 4))
            idx, = np.where(~np.isnan(self.A[:,0]))
            grp.create_dataset('intrapix/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('intrapix/nobs', data = self.A_nobs[idx], dtype = 'uint32')
            grp.create_dataset('intrapix/sinx', data = self.A[idx,0], dtype = 'float32')
            grp.create_dataset('intrapix/cosx', data = self.A[idx,1], dtype = 'float32')
            grp.create_dataset('intrapix/siny', data = self.A[idx,2], dtype = 'float32')
            grp.create_dataset('intrapix/cosy', data = self.A[idx,3], dtype = 'float32')
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            # Write the sky transmission.
            idx, lstseq = np.where(~np.isnan(self.s))
            grp.create_dataset('clouds/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('clouds/lstseq', data = lstseq + self.lstmin, dtype = 'uint32')
            grp.create_dataset('clouds/nobs', data = self.s_nobs[idx, lstseq], dtype = 'uint32')
            grp.create_dataset('clouds/clouds', data = self.s[idx, lstseq], dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('clouds/sigma', data = self.sigma2[idx, lstseq], dtype = 'float32')
            
            grp['clouds'].attrs['grid'] = 'healpix'
            grp['clouds'].attrs['nx'] = self.skynx
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = lstmax
            grp['clouds'].attrs['lstlen'] = lstlen
        
        return

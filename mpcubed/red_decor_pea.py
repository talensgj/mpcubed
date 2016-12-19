#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import misc
import IO
from coordinates import grids
from systematics import cdecor_pea


class CoarseDecorrelation(object):
    
    def __init__(self, LBfile, aperture, sysfile = None, **kwargs):
        """ Perform a coarse decorrelation on all data in a given file."""
        
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
            prefix = 'sys%i_pea_test_'%self.aper
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
        self.sigmas = True # Not an option in polar equal area reduction.
        
        self.camgrid = 'polar'
        self.camnx = kwargs.pop('camnx', 13500)
        self.camny = kwargs.pop('camny', 720)
        
        self.ipxgrid = 'polar'
        self.ipxnx = kwargs.pop('ipxnx', 270)
        self.ipxny = self.camny
        
        self.skygrid = 'polar_eqarea'
        self.skynx = kwargs.pop('skynx', 23)
    
        # Options to coarse_decor functions.
        self.maxiter = kwargs.pop('maxiter', 20)
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.verbose = kwargs.pop('verbose', True)
    
        # Perform the coarse decorrelation.
        self._calculate()
    
        return
    
    def _read_data(self, here): # Clean up the funcion call.
        """ Read a portion of the data, create indices and remove flagged
        datapoints."""
        
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
        _, _, skyidx = self.skygrid.radec2idx(ra, dec)
        
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
        
        ### FIXED V
        mag = mag - self.vmag[staridx]
        
        return mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, skyidx, lstseq
    
    def _calculate(self):
        """ Perform the coarse decorrelation."""

        # Set up the IO and coordinate grids.
        self.f = IO.fLCfile(self.LBfile)
        self.camgrid = grids.PolarGrid(self.camnx, self.camny)
        self.ipxgrid = grids.PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = grids.PolarEAGrid(self.skynx)
        
        # Read the required header data.
        self.ascc, self.ra, self.dec, self.nobs, self.vmag = self.f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.nobs = self.nobs.astype('int')
        
        # Create indices.
        self.staridx = np.arange(len(self.ascc))
        _, self.decidx = self.camgrid.radec2idx(self.ra, self.dec)
        ring, cell, self.skyidx = self.skygrid.radec2idx(self.ra, self.dec)
        
        # Read global information.
        with h5py.File(self.LBfile, 'r') as f:
            
            filelist = f['global/filelist'].value
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
        self.m = np.copy(self.vmag)
        self.z = np.full((self.camgrid.npix,), fill_value=np.nan)
        self.A = np.full((self.ipxgrid.npix, 4), fill_value=np.nan)
        self.s = np.array([])
        
        if self.sigmas:
            self.sigma1 = np.full(len(self.ascc), fill_value=0)
            self.sigma2 = np.array([])
            
        widx = np.array([], dtype='int')
        
        self.nobs_m = np.zeros(len(self.ascc), dtype='uint32')
        self.nobs_z = np.zeros((self.camgrid.npix,), dtype='uint32')
        self.nobs_A = np.zeros((self.ipxgrid.npix,), dtype='uint32')
        self.nobs_s = np.array([], dtype='int')
        
        niter = np.zeros(self.skynx + 2)
        chisq = np.zeros(self.skynx + 2)
        npoints = np.zeros(self.skynx + 2)
        npars = np.zeros(self.skynx + 2)
        
        # Perform the coarse decorrelation.
        print 'Performing coarse decorrelation for:', self.LBfile
        for i in range(self.skynx + 2):

            here = (ring == i)
            if np.sum(here) == 0:
                continue
        
            mag, emag, x, y, staridx, decidx, camidx, ipxidx, skyidx, lstseq = self._read_data(here)
        
            camidx = np.ravel_multi_index((camidx, decidx), (self.camgrid.nx+2, self.camgrid.ny+2))
            skyidx = np.ravel_multi_index((skyidx, lstseq), (self.skygrid.npix, lstlen)) 
            ipxidx = np.ravel_multi_index((ipxidx, decidx), (self.ipxgrid.nx+2, self.ipxgrid.ny+2))

            staridx, idx1 = np.unique(staridx, return_inverse=True)
            camidx, idx2 = np.unique(camidx, return_inverse=True)
            skyidx, idx3 = np.unique(skyidx, return_inverse=True)
            ipxidx, idx4 = np.unique(ipxidx, return_inverse=True)

            self.nobs_m[staridx] = np.bincount(idx1)
            self.nobs_z[camidx] = np.bincount(idx2)
            self.nobs_A[ipxidx] = np.bincount(idx4)
            self.nobs_s = np.append(self.nobs_s, np.bincount(idx3))

            sort = np.argsort(idx4)
            idx1 = idx1[sort]
            idx2 = idx2[sort]
            idx3 = idx3[sort]
            idx4 = idx4[sort]
            mag = mag[sort]
            emag = emag[sort]
            x = x[sort]
            y = y[sort]

            print 'Computing solution for ring', i

            self.z[camidx], s, self.A[ipxidx], self.sigma1[staridx], sigma2, quality = cdecor_pea.cdecor(idx1, idx2, idx3, idx4, mag, emag, x, y, maxiter=self.maxiter, dtol=self.dtol, verbose=self.verbose)
        
            widx = np.append(widx, skyidx)
            self.s = np.append(self.s, s)
            self.sigma2 = np.append(self.sigma2, sigma2)
        
            print 'niter = ', quality.niter
        
            niter[i] = quality.niter
            chisq[i] = quality.chisq
            npoints[i] = quality.npoints
            npars[i] = quality.npars
        
        self.z = self.z.reshape((self.camgrid.nx+2, self.camgrid.ny+2))
        self.A = self.A.reshape((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4))
        
        self.nobs_z = self.nobs_z.reshape((self.camgrid.nx+2, self.camgrid.ny+2))
        self.nobs_A = self.nobs_A.reshape((self.ipxgrid.nx+2, self.ipxgrid.ny+2))
        
        # Write the results to file.
        with h5py.File(self.sysfile) as f:
    
            # Write the header.
            hdr = f.create_group('header')
            
            hdr.create_dataset('filelist', data=filelist)
            hdr.attrs['station'] = station
            hdr.attrs['camera'] = camera
            
            hdr.attrs['alt0'] = np.mean(alt0)
            hdr.attrs['az0'] = np.mean(az0)
            hdr.attrs['th0'] = np.mean(th0)
            hdr.attrs['x0'] = np.mean(x0)
            hdr.attrs['y0'] = np.mean(y0)
            
            #hdr.attrs['outer_maxiter'] = self.outer_maxiter
            hdr.attrs['maxiter'] = self.maxiter
            hdr.attrs['sigmas'] = self.sigmas
            hdr.attrs['dtol'] = self.dtol
            
            hdr.create_dataset('niter', data = niter, dtype = 'uint32')
            hdr.create_dataset('chisq', data = chisq, dtype = 'float64')
            hdr.create_dataset('npoints', data = npoints, dtype = 'uint32')
            hdr.create_dataset('npars', data = npars, dtype = 'uint32')
            
            # Write the data.
            grp = f.create_group('data')
            
            # Write the magnitudes.
            grp.create_dataset('magnitudes/ascc', data = self.ascc)
            grp.create_dataset('magnitudes/vmag', data = self.vmag, dtype = 'float32')
            grp.create_dataset('magnitudes/nobs', data = self.nobs_m, dtype = 'uint32')
            grp.create_dataset('magnitudes/mag', data = self.m, dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('magnitudes/sigma', data = self.sigma1, dtype = 'float32')
            
            # Write the camera transmission.
            idx1, idx2 = np.where(~np.isnan(self.z))
            grp.create_dataset('trans/idx1', data = idx1, dtype = 'uint32')
            grp.create_dataset('trans/idx2', data = idx2, dtype = 'uint32')
            grp.create_dataset('trans/nobs', data = self.nobs_z[idx1, idx2], dtype = 'uint32')
            grp.create_dataset('trans/trans', data = self.z[idx1, idx2], dtype = 'float32')
            
            grp['trans'].attrs['grid'] = 'polar'
            grp['trans'].attrs['nx'] = self.camnx
            grp['trans'].attrs['ny'] = self.camny
            
            # Write the intrapixel variations.
            idx1, idx2 = np.where(~np.isnan(self.A[:,:,0]))
            grp.create_dataset('intrapix/idx1', data = idx1, dtype = 'uint32')
            grp.create_dataset('intrapix/idx2', data = idx2, dtype = 'uint32')
            grp.create_dataset('intrapix/nobs', data = self.nobs_A[idx1, idx2], dtype = 'uint32')
            grp.create_dataset('intrapix/sinx', data = self.A[idx1, idx2, 0], dtype = 'float32')
            grp.create_dataset('intrapix/cosx', data = self.A[idx1, idx2, 1], dtype = 'float32')
            grp.create_dataset('intrapix/siny', data = self.A[idx1, idx2, 2], dtype = 'float32')
            grp.create_dataset('intrapix/cosy', data = self.A[idx1, idx2, 3], dtype = 'float32')
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            # Write the sky transmission.
            idx, lstseq = np.unravel_index(widx, (self.skygrid.npix, lstlen))
            grp.create_dataset('clouds/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('clouds/lstseq', data = lstseq + self.lstmin, dtype = 'uint32')
            grp.create_dataset('clouds/nobs', data = self.nobs_s, dtype = 'uint32')
            grp.create_dataset('clouds/clouds', data = self.s, dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('clouds/sigma', data = self.sigma2, dtype = 'float32')
            
            grp['clouds'].attrs['grid'] = 'polarea'
            grp['clouds'].attrs['nx'] = self.skynx
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = lstmax
            grp['clouds'].attrs['lstlen'] = lstlen
        
        return

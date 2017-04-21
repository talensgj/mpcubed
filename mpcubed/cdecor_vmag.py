#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

import misc
import IO
from coordinates import grids
from systematics import sigmas

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def cdecor(idx2, value, error, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation.

    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        value (float): Values to fit.
        error (float): Measurement errors corresponding to the values.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.

    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.

    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars = npars2
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*value)/np.bincount(idx2, weights)
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit2 < dtol):
                break
        
        par2_old = np.copy(par2)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par2[idx2])**2        
    chisq = np.sum(chisq)
    
    return par2, Quality(niter, chisq, npoints, npars)

def cdecor_intrapix(idx2, idx3, value, error, x, y, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with intrapixel variations.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        idx3 (int): Indices along which to calculate the intrapixel variations.
        value (float): Values to fit.
        error (float): Measurement errors corresponding to the values.
        x (float): The x position corresponding to the values.
        y (float): The y position corresponding to the values. 
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        par3 (float): The amplitudes of the intrapixel variations corresponding
            to idx3.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    sort = np.argsort(idx3)
    idx2 = idx2[sort]
    idx3 = idx3[sort]
    value = value[sort]
    error = error[sort]
    x = x[sort]
    y = y[sort]
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars3 = 4*(np.amax(idx3) + 1)
    npars = npars2 + npars3
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    par3 = np.zeros((npars3/4, 4))
    
    strides = np.cumsum(np.bincount(idx3))
    strides = np.append(0, strides)
    mat = np.vstack([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T
    ipx = np.zeros(len(value))
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*(value - ipx))/np.bincount(idx2, weights)
        
        res = value - par2[idx2]
        wsqrt = np.sqrt(weights)
        for i in range(npars3/4):
            par3[i] = np.linalg.lstsq(mat[strides[i]:strides[i+1],:]*wsqrt[strides[i]:strides[i+1]:,None], res[strides[i]:strides[i+1]]*wsqrt[strides[i]:strides[i+1]])[0]
            ipx[strides[i]:strides[i+1]] = np.sum(mat[strides[i]:strides[i+1],:]*par3[i], axis=1)   
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            dcrit3 = np.nanmax(np.abs(par3 - par3_old))
            
            if (dcrit2 < dtol) & (dcrit3 < dtol):
                break
        
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par2[idx2] - ipx)**2        
    chisq = np.sum(chisq)
    
    return par2, par3, Quality(niter, chisq, npoints, npars)
    
def cdecor_sigmas(idx1, idx2, value, error, sigma1, sigma2, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with extra error terms.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        error (float): Measurement errors corresponding to the values.
        sigma1 (float): Initial value for the extra error corresponding to
            idx1.
        sigma2 (float): Initial value for the extra error corresponding to
            idx2.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        sigma1 (float): The extra error corresponding to idx1.
        sigma2 (float): The extra error corresponding to idx2.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars = npars2
    
    # Create arrays.
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
            
        # Compute the parameters.
        sigma1 = sigmas.find_sigma(idx1, value - par2[idx2], error**2 + (sigma2**2)[idx2])
        par2, sigma2 = sigmas.find_par_sigma(idx2, value, error**2 + (sigma1**2)[idx1])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit2 < dtol):
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_older))
            
            if (dcrit2 < dtol):
                break
        
        if (niter > 0):
            par2_older = np.copy(par2_old)
        
        par2_old = np.copy(par2)
        
    # Compute the chi-square of the fit.
    chisq = (value - par2[idx2])**2/(error**2 + (sigma1**2)[idx1] + (sigma2**2)[idx2])     
    chisq = np.sum(chisq)
    
    return par2, sigma1, sigma2, Quality(niter, chisq, npoints, npars)
    
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
            prefix = 'sys%i_vmag_'%self.aper
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
        skyidx = np.repeat(skyidx, nobs)
        
        # Remove bad data.
        here = (flux > 0) & (eflux > 0) & (flags < 1) & (sky > 0)
        
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
    
    def _spatial(self):
        """ Solve for the time-independent camera transmission and intrapixel 
        variations.
        """
        
        decidx, decuni = np.unique(self.decidx, return_inverse=True)
    
        nbins = len(decidx)
        for idx in range(nbins):
            
            # Read data.
            here = (decuni == idx)
            mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self._read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known temporal correction.
            if self.got_sky:
                mag = mag - self.s[skyidx, lstseq]
            
            if self.got_sky & self.sigmas:
                emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + self.sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            camtransidx, ind2 = np.unique(camtransidx, return_inverse=True)
            intrapixidx, ind3 = np.unique(intrapixidx, return_inverse=True)
            
            # Calculate new spatial correction.
            z, A, quality = cdecor_intrapix(ind2, ind3, mag, emag, x, y, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Store results.
            self.z[camtransidx, decidx[idx]] = z
            self.A[intrapixidx, decidx[idx]] = A
            
            self.spatial_niter[idx] = quality.niter
            self.spatial_chisq[idx] = quality.chisq
            self.spatial_npoints[idx] = quality.npoints
            self.spatial_npars[idx] = quality.npars 
            
            self.nobs_z[camtransidx, decidx[idx]] = np.bincount(ind2)
            self.nobs_A[intrapixidx, decidx[idx]] = np.bincount(ind3)
            
        return
        
    def _temporal(self):
        """ Solve for the time-dependent sky transmission."""
        
        skyidx, skyuni = np.unique(self.skyidx, return_inverse=True)
        
        nbins = len(skyidx)
        for idx in range(nbins):
            
            # Read data.
            here = (skyuni == idx)
            mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, _, lstseq = self._read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known spatial correction.
            mag = mag - self.z[camtransidx, decidx]
            mag = mag - np.sum(self.A[intrapixidx, decidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            # Create unique indices.
            staridx, ind1 = np.unique(staridx, return_inverse=True)
            lstseq, ind2 = np.unique(lstseq, return_inverse=True)
            
            # Calculate new temporal correction.
            if self.sigmas:
                s, sigma1, sigma2, quality = cdecor_sigmas(ind1, ind2, mag, emag, self.sigma1[staridx], self.sigma2[skyidx[idx], lstseq], maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            else:
                s, quality = cdecor(ind2, mag, emag, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Store results.
            self.s[skyidx[idx], lstseq] = s
            
            if self.sigmas:
                self.sigma1[staridx] = sigma1
                self.sigma2[skyidx[idx], lstseq] = sigma2
            
            self.temporal_niter[idx] = quality.niter
            self.temporal_chisq[idx] = quality.chisq
            self.temporal_npoints[idx] = quality.npoints
            self.temporal_npars[idx] = quality.npars 
            
            self.nobs_m[staridx] = np.bincount(ind1)
            self.nobs_s[skyidx[idx], lstseq] = np.bincount(ind2)
            
        self.got_sky = True
            
        return
    
    def _calculate(self):
        """ Perform the coarse decorrelation."""

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
        nbins = len(np.unique(self.decidx))
        self.spatial_niter = np.full(nbins, fill_value = np.nan)
        self.spatial_chisq = np.full(nbins, fill_value = np.nan)
        self.spatial_npoints = np.full(nbins, fill_value = np.nan)
        self.spatial_npars = np.full(nbins, fill_value = np.nan)
        
        nbins = len(np.unique(self.skyidx))
        self.temporal_niter = np.full(nbins, fill_value = np.nan)
        self.temporal_chisq = np.full(nbins, fill_value = np.nan)
        self.temporal_npoints = np.full(nbins, fill_value = np.nan)
        self.temporal_npars = np.full(nbins, fill_value = np.nan)
        
        self.m = np.copy(self.vmag)
        self.z = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.A = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan)
        self.s = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        if self.sigmas:
            self.sigma1 = np.full(len(self.ascc), fill_value=0.)
            self.sigma2 = np.full((self.skygrid.npix, lstlen), fill_value=0.)
        
        self.nobs_m = np.full(len(self.ascc), fill_value=np.nan)
        self.nobs_z = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.nobs_A = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2), fill_value=np.nan)
        self.nobs_s = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.got_sky = False
        
        # Perform the coarse decorrelation.
        print 'Performing coarse decorrelation for:', self.LBfile
        for niter in range(self.outer_maxiter):
        
            print 'Iteration %i out of %i:'%(niter + 1, self.outer_maxiter)
            
            print '    Calculating camera systematics...'
            self._spatial()
            
            print '    Calculating atmospheric systematics...'
            self._temporal()
        
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
            
            hdr.attrs['outer_maxiter'] = self.outer_maxiter
            hdr.attrs['inner_maxiter'] = self.inner_maxiter
            hdr.attrs['sigmas'] = self.sigmas
            hdr.attrs['dtol'] = self.dtol
           
            hdr.create_dataset('spatial/niter', data = self.spatial_niter, dtype = 'uint32')
            hdr.create_dataset('spatial/chisq', data = self.spatial_chisq, dtype = 'float64')
            hdr.create_dataset('spatial/npoints', data = self.spatial_npoints, dtype = 'uint32')
            hdr.create_dataset('spatial/npars', data = self.spatial_npars, dtype = 'uint32')
            
            hdr.create_dataset('temporal/niter', data = self.temporal_niter, dtype = 'uint32')
            hdr.create_dataset('temporal/chisq', data = self.temporal_chisq, dtype = 'float64')
            hdr.create_dataset('temporal/npoints', data = self.temporal_npoints, dtype = 'uint32')
            hdr.create_dataset('temporal/npars', data = self.temporal_npars, dtype = 'uint32')
            
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
            idx, lstseq = np.where(~np.isnan(self.s))
            grp.create_dataset('clouds/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('clouds/lstseq', data = lstseq + self.lstmin, dtype = 'uint32')
            grp.create_dataset('clouds/nobs', data = self.nobs_s[idx, lstseq], dtype = 'uint32')
            grp.create_dataset('clouds/clouds', data = self.s[idx, lstseq], dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('clouds/sigma', data = self.sigma2[idx, lstseq], dtype = 'float32')
            
            grp['clouds'].attrs['grid'] = 'healpix'
            grp['clouds'].attrs['nx'] = self.skynx
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = lstmax
            grp['clouds'].attrs['lstlen'] = lstlen
        
        return

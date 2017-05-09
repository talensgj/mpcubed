#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

import multiprocessing as mp

from . import misc
from . import IO
from .coordinates import grids
from .systematics import sigmas

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def cdecor_spatial(idx2, idx3, value, error, x, y, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with intrapixel variations.
    
    Args:
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
    
def cdecor_temporal(idx1, idx2, value, error, sigma1, sigma2, maxiter=100, dtol=1e-3, verbose=True):
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
    
def spatial_worker(in_queue, out_queue, maxiter, dtol, verbose):
    
    while True:
        
        item = in_queue.get()
    
        if (item == 'DONE'):
            break
        else:
            idx, camtransidx, intrapixidx, idx2, idx4, mag, emag, x, y = item
            trans, amp, quality = cdecor_spatial(idx2, idx4, mag, emag, x, y, maxiter, dtol, verbose)
            out_queue.put((idx, camtransidx, intrapixidx, trans, amp, quality))

    return    
    
def temporal_worker(in_queue, out_queue, maxiter, dtol, verbose):
    
    while True:
        
        item = in_queue.get()
        
        if (item == 'DONE'):
            break
        else:
            idx, staridx, lstseq, idx1, idx3, mag, emag, sig1, sig2 = item
            clouds, sig1, sig2, quality = cdecor_temporal(idx1, idx3, mag, emag, sig1, sig2, maxiter, dtol, verbose)
            out_queue.put((idx, staridx, lstseq, clouds, sig1, sig2, quality))

    return 
    
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
            prefix = 'sys%i_vmag_pea_'%self.aper
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
        self.outer_maxiter = kwargs.pop('outer_maxiter', 5)
        self.nprocs = kwargs.pop('nprocs', 4)
        
        self.camgrid = 'polar'
        self.camnx = kwargs.pop('camnx', 13500)
        self.camny = kwargs.pop('camny', 720)
        
        self.ipxgrid = 'polar'
        self.ipxnx = kwargs.pop('ipxnx', 270)
        self.ipxny = self.camny
        
        self.skygrid = 'polarea'
        self.skynx = kwargs.pop('skynx', 23)
    
        # Options to coarse_decor functions.
        self.inner_maxiter = kwargs.pop('inner_maxiter', 100)
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.verbose = kwargs.pop('verbose', False)
    
        # Perform the coarse decorrelation.
        self._calculate()
    
        return
    
    def _read_data(self, here):
        """ Read a portion of the data, create indices and remove flagged
        datapoints."""
        
        ascc = self.stars['ascc'][here]
        ra = self.stars['ra'][here]
        dec = self.stars['dec'][here]
        nobs = self.stars['nobs'][here]
        
        staridx = self.staridx[here]
        skyidx = self.skyidx[here]
        
        # Read data.
        fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'lst', 'lstseq', 'flag']
        lightcurves = self.f.read_lightcurves(ascc, fields, perstar=False)
        
        flux = lightcurves['flux%i'%self.aper]
        eflux = lightcurves['eflux%i'%self.aper]
        sky = lightcurves['sky']
        x = lightcurves['x']
        y = lightcurves['y']
        lst = lightcurves['lst']
        lstseq = lightcurves['lstseq']
        flags = lightcurves['flag']      
        
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
        mag = mag - self.stars['vmag'][staridx]
        
        return mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, skyidx, lstseq
    
    def _spatial(self):
        """ Solve for the time-independent camera transmission and intrapixel 
        variations.
        """
        
        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, spatial_worker, (in_queue, out_queue, self.inner_maxiter, self.dtol, self.verbose))        
        
        decidx = np.unique(self.decidx)
        
        for idx in decidx:
            
            # Read data.
            here = (self.decidx == idx)
            mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self._read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known temporal correction.
            if self.got_sky:
                mag = mag - self.clouds['clouds'][skyidx, lstseq]
                emag = np.sqrt(emag**2 + self.magnitudes['sigma'][staridx]**2 + self.clouds['sigma'][skyidx, lstseq]**2)
            
            # Create unique indices.
            camtransidx, idx2 = np.unique(camtransidx, return_inverse=True)
            intrapixidx, idx3 = np.unique(intrapixidx, return_inverse=True)
            
            self.trans['nobs'][camtransidx, idx] = np.bincount(idx2)
            self.intrapix['nobs'][intrapixidx, idx] = np.bincount(idx3)
            
            in_queue.put((idx, camtransidx, intrapixidx, idx2, idx3, mag, emag, x, y))            
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')            
            
        for item in iter(out_queue.get, 'DONE'):
            
            idx, camtransidx, intrapixidx, trans, amplitudes, quality = item        
            
            # Store results.
            self.trans['trans'][camtransidx, idx] = trans
            self.intrapix['amplitudes'][intrapixidx, idx] = amplitudes
            
            self.spatial['niter'][idx] = quality.niter
            self.spatial['chisq'][idx] = quality.chisq
            self.spatial['npoints'][idx] = quality.npoints
            self.spatial['npars'][idx] = quality.npars 
            
        return
        
    def _temporal(self):
        """ Solve for the time-dependent sky transmission."""
        
        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, temporal_worker, (in_queue, out_queue, self.inner_maxiter, self.dtol, self.verbose))         
        
        skyidx = np.unique(self.skyidx)
        
        for idx in skyidx:
            
            # Read data.
            here = (self.skyidx == idx)
            mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, _, lstseq = self._read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known spatial correction.
            mag = mag - self.trans['trans'][camtransidx, decidx]
            mag = mag - np.sum(self.intrapix['amplitudes'][intrapixidx, decidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            # Create unique indices.
            staridx, idx1 = np.unique(staridx, return_inverse=True)
            lstseq, idx2 = np.unique(lstseq, return_inverse=True)
            
            self.magnitudes['nobs'][staridx] = np.bincount(idx1)
            self.clouds['nobs'][idx, lstseq] = np.bincount(idx2)

            in_queue.put((idx, staridx, lstseq, idx1, idx2, mag, emag, self.magnitudes['sigma'][staridx], self.clouds['sigma'][idx, lstseq]))
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')              
              
        for item in iter(out_queue.get, 'DONE'):
            
            idx, staridx, lstseq, clouds, sigma1, sigma2, quality = item            
            
            # Store results.
            self.magnitudes['sigma'][staridx] = sigma1
            
            self.clouds['clouds'][idx, lstseq] = clouds
            self.clouds['sigma'][idx, lstseq] = sigma2
            
            self.temporal['niter'][idx] = quality.niter
            self.temporal['chisq'][idx] = quality.chisq
            self.temporal['npoints'][idx] = quality.npoints
            self.temporal['npars'][idx] = quality.npars 
            
        self.got_sky = True
            
        return
    
    def _calculate(self):
        """ Perform the coarse decorrelation."""

        # Set up the IO and coordinate grids.
        self.f = IO.PhotFile(self.LBfile)
        self.camgrid = grids.PolarGrid(self.camnx, self.camny)
        self.ipxgrid = grids.PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = grids.PolarEAGrid(self.skynx)
        
        # Read the required header data.
        self.stars = self.f.read_stars(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.stars['nobs'] = self.stars['nobs'].astype('int')
        
        # Create indices.
        self.staridx = np.arange(len(self.stars['ascc']))
        _, self.decidx = self.camgrid.radec2idx(self.stars['ra'], self.stars['dec'])
        _, _, self.skyidx = self.skygrid.radec2idx(self.stars['ra'], self.stars['dec'])
        
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
        
        # The spatial calculation statistics.
        self.spatial = dict()
        self.spatial['niter'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        self.spatial['chisq'] = np.full(self.camgrid.ny+2, fill_value=np.nan)
        self.spatial['npoints'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        self.spatial['npars'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        
        # The temporal calculation statistics.
        self.temporal = dict()
        self.temporal['niter'] = np.zeros(self.skygrid.npix, dtype='uint32')
        self.temporal['chisq'] = np.full(self.skygrid.npix, fill_value=np.nan)
        self.temporal['npoints'] = np.zeros(self.skygrid.npix, dtype='uint32')
        self.temporal['npars'] = np.zeros(self.skygrid.npix, dtype='uint32')
        
        # The magnitudes.
        self.magnitudes = dict()
        self.magnitudes['ascc'] = self.stars['ascc'] 
        self.magnitudes['vmag'] = self.stars['vmag']
        self.magnitudes['nobs'] = np.zeros(len(self.stars['ascc']), dtype='uint32')
        self.magnitudes['mag'] = np.copy(self.stars['vmag'])
        self.magnitudes['sigma'] = np.zeros(len(self.stars['ascc']))
        
        # The transmission map.
        self.trans = dict()
        self.trans['nobs'] = np.zeros((self.camgrid.nx+2, self.camgrid.ny+2), dtype='uint32')
        self.trans['trans'] = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        
        # The intrapixel amplitudes.
        self.intrapix = dict()
        self.intrapix['nobs'] = np.zeros((self.ipxgrid.nx+2, self.ipxgrid.ny+2), dtype='uint32')
        self.intrapix['amplitudes'] = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan)
        
        # The cloud corrections.
        self.clouds = dict()
        self.clouds['nobs'] = np.zeros((self.skygrid.npix, lstlen), dtype='uint32')
        self.clouds['clouds'] = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        self.clouds['sigma'] = np.zeros((self.skygrid.npix, lstlen))
        
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
            hdr.attrs['sigmas'] = True
            hdr.attrs['dtol'] = self.dtol
           
            idx1 = np.unique(self.decidx)
            for key in self.spatial.keys():
                hdr.create_dataset('spatial/' + key, data=self.spatial[key][idx1])
            
            idx1 = np.unique(self.skyidx)
            for key in self.temporal.keys():
                hdr.create_dataset('temporal/' + key, data=self.temporal[key][idx1])
            
            # Write the data.
            grp = f.create_group('data')
            
            # Write the magnitudes.
            grp.create_dataset('magnitudes/ascc', data=self.magnitudes['ascc'])
            grp.create_dataset('magnitudes/vmag', data=self.magnitudes['vmag'], dtype='float32')
            grp.create_dataset('magnitudes/nobs', data=self.magnitudes['nobs'])
            grp.create_dataset('magnitudes/mag', data=self.magnitudes['mag'], dtype='float32')
            grp.create_dataset('magnitudes/sigma', data=self.magnitudes['sigma'], dtype='float32')

            # Write the camera transmission.
            idx1, idx2 = np.where(self.trans['nobs'] > 0)            
            grp.create_dataset('trans/idx1', data=idx1, dtype='uint32')
            grp.create_dataset('trans/idx2', data=idx2, dtype='uint32')          
            grp.create_dataset('trans/nobs', data=self.trans['nobs'][idx1,idx2])
            grp.create_dataset('trans/trans', data=self.trans['trans'][idx1,idx2], dtype='float32')
            
            grp['trans'].attrs['grid'] = 'polar'
            grp['trans'].attrs['nx'] = self.camnx
            grp['trans'].attrs['ny'] = self.camny
            
            # Write the intrapixel variations.
            idx1, idx2 = np.where(self.intrapix['nobs'] > 0)            
            grp.create_dataset('intrapix/idx1', data=idx1, dtype='uint32')
            grp.create_dataset('intrapix/idx2', data=idx2, dtype='uint32')
            grp.create_dataset('intrapix/nobs', data=self.intrapix['nobs'][idx1,idx2])
            grp.create_dataset('intrapix/sinx', data=self.intrapix['amplitudes'][idx1,idx2,0], dtype='float32')
            grp.create_dataset('intrapix/cosx', data=self.intrapix['amplitudes'][idx1,idx2,1], dtype='float32')
            grp.create_dataset('intrapix/siny', data=self.intrapix['amplitudes'][idx1,idx2,2], dtype='float32')
            grp.create_dataset('intrapix/cosy', data=self.intrapix['amplitudes'][idx1,idx2,3], dtype='float32')
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            # Write the sky transmission.
            idx, lstseq = np.where(self.clouds['nobs'] > 0)            
            grp.create_dataset('clouds/idx', data=idx, dtype='uint32')
            grp.create_dataset('clouds/lstseq', data=lstseq+self.lstmin, dtype='uint32')
            grp.create_dataset('clouds/nobs', data=self.clouds['nobs'][idx, lstseq])
            grp.create_dataset('clouds/clouds', data=self.clouds['clouds'][idx, lstseq], dtype='float32')
            grp.create_dataset('clouds/sigma', data=self.clouds['sigma'][idx, lstseq], dtype='float32')
            
            grp['clouds'].attrs['grid'] = 'polarea'
            grp['clouds'].attrs['nx'] = self.skynx
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = lstmax
            grp['clouds'].attrs['lstlen'] = lstlen
        
        return

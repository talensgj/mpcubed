#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob
import datetime

import numpy as np

import multiprocessing as mp

from .. import io, misc, statistics
from . import grids, sigmas

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

###############################################################################
### Functions for reading the raw data.
###############################################################################

def baseline(date, part, camera, source, dest):
    
    # Check that the date is valid.
    try:
        datetime.datetime.strptime(date, '%Y%m')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMM")
        
    # Get the photometry files matching the date, part and camera.
    if (part == 'A'): 
     
        filelist = glob.glob(os.path.join(source, '*/*/f*_{}0?{}.hdf5'.format(date, camera)))
        filelist = filelist + glob.glob(os.path.join(source, '*/*/f*_{}1[0-5]{}.hdf5'.format(date, camera)))
    
    elif (part == 'B'):
    
        filelist = glob.glob(os.path.join(source, '*/*/f*_{}1[6-9]{}.hdf5'.format(date, camera)))
        filelist = filelist + glob.glob(os.path.join(source, '*/*/f*_{}[23]?{}.hdf5'.format(date, camera)))
    
    elif (part == 'C'):

        filelist = glob.glob(os.path.join(source, '*/*/f*_{}??{}.hdf5'.format(date, camera)))         
    
    else:
        raise ValueError("Part should be either 'A', 'B' or 'C'.")
    
    # Check that there are valid files.
    if len(filelist) == 0:
        print 'No valid data found.'
        return None

    # Sort the filelist.
    filelist = np.sort(filelist)    
    
    print 'Combining files:'
    for filename in filelist:
        print '    {}'.format(filename) 

    # Create the destination directory.    
    outpath = os.path.join(dest, camera)    
    misc.ensure_dir(outpath)
    
    # Determine the prefix.
    f = io.PhotFile(filelist[0])
    siteid = f.get_siteid()
    
    if siteid == 'LP':
        prefix = 'fLC'
    else:
        prefix = 'fast'
    
    # Check that the destination file does not exist.
    photfile = os.path.join(outpath, '{}_{}{}{}.hdf5'.format(prefix, date, part, camera)) 
    if os.path.isfile(photfile):
        print 'Output file already exists: {}'.format(photfile)
        return photfile
    else:
        print 'Writing results to: {}'.format(photfile)
    
    # Combine the files.
    io.make_baseline(photfile, filelist)
    
    return photfile

def _read_header(filename, camgrid, skygrid):
    
    f = io.PhotFile(filename)
    
    data = f.read_header()
    
    settings = dict()
    settings['filelist'] = data['filelist']
    
    if f.get_siteid() == 'LP':
        settings['station'] = data['station']
        settings['camera'] = data['camera']
        settings['alt0'] = np.mean(data['alt0'])
        settings['az0'] = np.mean(data['az0'])
        settings['th0'] = np.mean(data['th0'])
        settings['x0'] = np.mean(data['x0'])
        settings['y0'] = np.mean(data['y0'])
    else:
        settings['station'] = data['site-obs']
        settings['camera'] = data['cam-obs']
    
    # Ensure lstmin and lstmax are integer.
    lstmin = data['lstmin'].astype('int')
    lstmax = data['lstmax'].astype('int')
    
    stars = f.read_stars(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
    
    # Ensure nobs is integer.
    stars['nobs'] = stars['nobs'].astype('int')
    
    # Compute star dependent indices.
    stars['i'] = np.arange(len(stars['ascc']))
    _, stars['n'] = camgrid.radec2idx(stars['ra'], stars['dec'])
    stars['q'] = skygrid.radec2idx(stars['ra'], stars['dec'])
    
    return stars, settings, lstmin, lstmax

def _read_data(filename, aper, stars, index, method, camgrid, ipxgrid):
    
    # Select stars.
    if method == 'spatial':
        mask = stars['n'] == index
    elif method == 'temporal':
        mask = stars['q'] == index
    else:
        raise ValueError('Unknown reading method "{}"'.format(method))
    
    ascc = stars['ascc'][mask]
    ra = stars['ra'][mask]
    dec = stars['dec'][mask]
    nobs = stars['nobs'][mask]
    stars_i = stars['i'][mask]
    stars_q = stars['q'][mask]
    
    # Read data.
    f = io.PhotFile(filename)
    
    if f.get_siteid() == 'LP':
        fields = ['lstseq', 'flux{}'.format(aper), 'eflux{}'.format(aper), 'x', 'y', 'sky', 'lst', 'flag']
    else:
        fields = ['lstseq', 'flux{}'.format(aper), 'eflux{}'.format(aper), 'x', 'y', 'pflag', 'aflag']
    
    lightcurves = f.read_lightcurves(ascc, fields, perstar=False)
    
    lstseq = lightcurves['lstseq'].astype('int')
    x, y = lightcurves['x'], lightcurves['y']
    flux, eflux = lightcurves['flux{}'.format(aper)], lightcurves['eflux{}'.format(aper)]
    
    if f.get_siteid() == 'LP':
        
        sky = lightcurves['sky']
        lst = lightcurves['lst']
        flags = lightcurves['flag']
        
        mask = (flux > 0) & (eflux > 0) & (flags < 1) & (sky > 0)
        
        x, y = x.astype('float64'), y.astype('float64')
        
    else:
        
        aflag = lightcurves['aflag']
        pflag = lightcurves['pflag']
    
        mask = (aflag == 0) & (pflag == 0)
    
        station = f.read_station(['lst', 'exptime'], lstseq)        
        lst = station['lst']
        exptime = station['exptime']
        
        flux, eflux = flux/exptime, eflux/exptime
        
    ra = np.repeat(ra, nobs)
    dec = np.repeat(dec, nobs)
    ha = np.mod(lst*15. - ra, 360.)
    
    # Create indices.
    data_i = np.repeat(stars_i, nobs)
    data_t = lstseq
    data_k, data_n = camgrid.radec2idx(ha, dec)
    data_l, data_n = ipxgrid.radec2idx(ha, dec)
    data_q = np.repeat(stars_q, nobs)
    
    # Remove bad data.
    flux = flux[mask]
    eflux = eflux[mask]
    x = x[mask]
    y = y[mask]
    
    data_i = data_i[mask]
    data_t = data_t[mask]
    data_n = data_n[mask]
    data_k = data_k[mask]
    data_l = data_l[mask]
    data_q = data_q[mask]
    
    # Convert flux to magnitudes:
    mag, emag = misc.flux2mag(flux, eflux)
    
    # Create matrix for intrapixel variations.
    sinx, cosx = np.sin(2*np.pi*x), np.cos(2*np.pi*x)
    siny, cosy = np.sin(2*np.pi*y), np.cos(2*np.pi*y)
    b_mat = np.column_stack([sinx, cosx, siny, cosy])
    
    return mag, emag, b_mat, data_i, data_t, data_n, data_k, data_l, data_q
    
###############################################################################
### Code for calibrating the raw data, inspired by Collier Cameron+ 2006.
###############################################################################

def cdecor_spatial(idx_k, idx_l, mag, emag, b_mat, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with intrapixel variations.
    
    Args:
        idx_k (int): Indices for calculating the transmission corrections.
        idx_l (int): Indices for calculating the intrapixel corrections.
        mag (float): Magnitude points.
        emag (float): The uncertainties corresponding to the magnitude points.
        b_mat (float): Matrix of sine functions for the intrapixel corrections. 
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is False.
        
    Returns:
        T (float): The transmission corrections.
        a (float): The intrapixel corrections.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    sort = np.argsort(idx_l)
    idx_k = idx_k[sort]
    idx_l = idx_l[sort]
    mag = mag[sort]
    emag = emag[sort]
    b_mat = b_mat[sort]
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars_T = np.amax(idx_k) + 1
    npars_a = np.amax(idx_l) + 1
    npars = npars_T + 4*npars_a
    
    # Create arrays.
    weights = 1./emag**2
    T = np.zeros((npars_T,))
    a = np.zeros((npars_a, 4))
    ipx = np.zeros(npoints)
    
    strides = np.cumsum(np.bincount(idx_l))
    strides = np.append(0, strides)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        T = np.bincount(idx_k, weights*(mag - ipx))/np.bincount(idx_k, weights)
        
        res = mag - T[idx_k]
        wsqrt = np.sqrt(weights)
        for i in range(npars_a):
            a[i] = np.linalg.lstsq(b_mat[strides[i]:strides[i+1],:]*wsqrt[strides[i]:strides[i+1]:,None], res[strides[i]:strides[i+1]]*wsqrt[strides[i]:strides[i+1]])[0]
            ipx[strides[i]:strides[i+1]] = np.sum(a[i]*b_mat[strides[i]:strides[i+1],:], axis=1)
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit_T = np.nanmax(np.abs(T - T_old))
            dcrit_a = np.nanmax(np.abs(a - a_old))
            
            if (dcrit_T < dtol) & (dcrit_a < dtol):
                break
        
        T_old = np.copy(T)
        a_old = np.copy(a)
    
    # Compute the chi-square of the fit.
    chisq = weights*(mag - T[idx_k] - ipx)**2        
    chisq = np.sum(chisq)
    
    return T, a, Quality(niter, chisq, npoints, npars)
    
def cdecor_temporal(idx_i, idx_t, mag, emag, sig_m, sig_c, m=None, maxiter=100, dtol=1e-3, verbose=False):
    """ Perform a coarse decorrelation with extra error terms.
    
    Args:
        idx_i (int): Indices for calculating the stellar magnitudes.
        idx_t (int): Indices for calculating the atmospheric corrections.
        mag (float): Magnitude points.
        emag (float): The uncertainties corresponding to the magnitude points.
        sig_m (float): Initial value for the extra error corresponding to
            idx1.
        sig_c (float): Initial value for the extra error corresponding to
            idx2.
        m (float): Values for the stellar magnitudes. If provided they are
            assumed to be fixed.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        m (float): The stellar magnitudes.
        c (float): The atmospheric corrections.
        sig_m (float): Extra error corresponding to the stellar magnitudes.
        sig_c (float): Extra error corresponding to the atmospheric corrections.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars_m = np.amax(idx_i) + 1
    npars_c = np.amax(idx_t) + 1
    
    if m is not None:
        npars = npars_c
        fixed = True
    else:
        npars = npars_m + npars_c
        fixed = False
    
    # Create arrays.
    c = np.zeros(npars_c)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
            
        # Compute the parameters.
        if fixed:
            sig_m = sigmas.find_sigma(idx_i, mag - m[idx_i] - c[idx_t], emag**2 + (sig_c**2)[idx_t])
        else:
            m, sig_m = sigmas.find_par_sigma(idx_i, mag - c[idx_t], emag**2 + (sig_c**2)[idx_t])
            
        c, sig_c = sigmas.find_par_sigma(idx_t, mag - m[idx_i], emag**2 + (sig_m**2)[idx_i])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit_m = np.nanmax(np.abs(m - m_old))
            dcrit_c = np.nanmax(np.abs(c - c_old))
            
            if (dcrit_m < dtol) & (dcrit_c < dtol):
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            dcrit_m = np.nanmax(np.abs(m - m_older))
            dcrit_c = np.nanmax(np.abs(c - c_older))
            
            if (dcrit_m < dtol) & (dcrit_c < dtol):
                break
        
        if (niter > 0):
            m_older = np.copy(m_old)
            c_older = np.copy(c_old)
        
        m_old = np.copy(m)
        c_old = np.copy(c)
        
    # Compute the chi-square of the fit.
    chisq = (mag - m[idx_i] - c[idx_t])**2/(emag**2 + (sig_m**2)[idx_i] + (sig_c**2)[idx_t])     
    chisq = np.sum(chisq)
    
    return m, c, sig_m, sig_c, Quality(niter, chisq, npoints, npars)
    
def worker_spatial(in_queue, out_queue, kwargs):
    
    while True:
        
        item = in_queue.get()
    
        if (item == 'DONE'):
            break
        else:
            n, k, l, idx_k, idx_l, mag, emag, b_mat = item
            T_kn, a_ln, quality = cdecor_spatial(idx_k, idx_l, mag, emag, b_mat, **kwargs)
            out_queue.put((n, k, l, T_kn, a_ln, quality))

    return    
    
def worker_temporal(in_queue, out_queue, kwargs):
    
    while True:
        
        item = in_queue.get()
        
        if (item == 'DONE'):
            break
        else:
            q, i, t, idx_i, idx_t, mag, emag, sig_i, sig_qt, m_i = item
            m_i, c_qt, sig_i, sig_qt, quality = cdecor_temporal(idx_i, idx_t, mag, emag, sig_i, sig_qt, m_i, **kwargs)
            out_queue.put((q, i, t, m_i, c_qt, sig_i, sig_qt, quality))

    return

class CoarseDecorrelation(object):
    
    def __init__(self, num_k=13500, num_l=270, num_n=720, num_q=8, **kwargs):
        """ Perform a coarse decorrelation on all data in a given file."""
        
        # Initialize the coordinate grids.
        self.camgrid = grids.PolarGrid(num_k, num_n)        
        self.ipxgrid = grids.PolarGrid(num_l, num_n)
        self.skygrid = grids.HealpixGrid(num_q)    
    
        # General parameters.
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.fixed = kwargs.pop('fixed', True) 
        self.nprocs = kwargs.pop('nprocs', 4)
        self.outer_maxiter = kwargs.pop('outer_maxiter', 5)
        self.inner_maxiter = kwargs.pop('inner_maxiter', 100)        
        
        # Parameters of the spatial and temporal solvers.
        self.kwargs = dict()
        self.kwargs['dtol'] = self.dtol
        self.kwargs['maxiter'] = self.inner_maxiter        
        self.kwargs['verbose'] = kwargs.pop('verbose', False)
    
        return
    
    def _spatial(self):
        """ Solve for the time-independent camera transmission and intrapixel 
        variations.
        """
        
        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, worker_spatial, (in_queue, out_queue, self.kwargs))        
        
        for n in self.spatial['n']:
            
            # Read data.
            mag, emag, b_mat, data_i, data_t, data_n, data_k, data_l, data_q = _read_data(self.photfile, self.aper, self.stars, n, 'spatial', self.camgrid, self.ipxgrid)
            data_t = data_t - self.clouds['lstmin']
            
            if (len(mag) == 0): continue
        
            mag = mag - self.magnitudes['mag'][data_i]
        
            # Apply known temporal correction.
            if self.got_sky:
                mag = mag - self.clouds['clouds'][data_q,data_t]
                emag = np.sqrt(emag**2 + self.magnitudes['sigma'][data_i]**2 + self.clouds['sigma'][data_q,data_t]**2)
            
            # Create unique indices.
            k, idx_k = np.unique(data_k, return_inverse=True)
            l, idx_l = np.unique(data_l, return_inverse=True)
            
            self.trans['nobs'][k,n] = np.bincount(idx_k)
            self.intrapix['nobs'][l,n] = np.bincount(idx_l)
            
            in_queue.put((n, k, l, idx_k, idx_l, mag, emag, b_mat))            
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')            
            
        for item in iter(out_queue.get, 'DONE'):
            
            n, k, l, T_kn, a_ln, quality = item        
            
            # Store results.
            self.trans['trans'][k,n] = T_kn
            self.intrapix['amplitudes'][l,n] = a_ln
            
            self.spatial['niter'][n] = quality.niter
            self.spatial['chisq'][n] = quality.chisq
            self.spatial['npoints'][n] = quality.npoints
            self.spatial['npars'][n] = quality.npars
            
        return
        
    def _temporal(self):
        """ Solve for the time-dependent sky transmission."""
        
        mngr = mp.Manager()
        in_queue = mp.Queue(2*self.nprocs)
        out_queue = mngr.Queue()
        the_pool = mp.Pool(self.nprocs, worker_temporal, (in_queue, out_queue, self.kwargs))         
        
        for q in self.temporal['q']:
            
            # Read data.
            mag, emag, b_mat, data_i, data_t, data_n, data_k, data_l, data_q = _read_data(self.photfile, self.aper, self.stars, q, 'temporal', self.camgrid, self.ipxgrid)
            data_t = data_t - self.clouds['lstmin']
            
            if (len(mag) == 0): continue
            
            # Apply known spatial correction.
            mag = mag - self.trans['trans'][data_k,data_n]
            mag = mag - np.sum(self.intrapix['amplitudes'][data_l,data_n]*b_mat, axis=1)
            
            # Create unique indices.
            i, idx_i = np.unique(data_i, return_inverse=True)
            t, idx_t = np.unique(data_t, return_inverse=True)
            
            self.magnitudes['nobs'][i] = np.bincount(idx_i)
            self.clouds['nobs'][q,t] = np.bincount(idx_t)

            if self.fixed:
                in_queue.put((q, i, t, idx_i, idx_t, mag, emag, self.magnitudes['sigma'][i], self.clouds['sigma'][q,t], self.magnitudes['mag'][i]))
            else:
                in_queue.put((q, i, t, idx_i, idx_t, mag, emag, self.magnitudes['sigma'][i], self.clouds['sigma'][q,t], None))
            
        for i in range(self.nprocs):
            in_queue.put('DONE')
            
        the_pool.close()
        the_pool.join()
                        
        out_queue.put('DONE')              
              
        for item in iter(out_queue.get, 'DONE'):
            
            q, i, t, m_i, c_qt, sig_i, sig_qt, quality = item            

            # Store results.
            self.magnitudes['mag'][i] = m_i
            self.magnitudes['sigma'][i] = sig_i
            
            self.clouds['clouds'][q,t] = c_qt
            self.clouds['sigma'][q,t] = sig_qt
            
            self.temporal['niter'][q] = quality.niter
            self.temporal['chisq'][q] = quality.chisq
            self.temporal['npoints'][q] = quality.npoints
            self.temporal['npars'][q] = quality.npars 
            
        self.got_sky = True
            
        return
    
    def __call__(self, photfile, aper, sysfile=None):
        """ Perform the coarse decorrelation."""

        self.photfile = photfile
        self.aper = aper
        self.sysfile = sysfile

        # Check if the photometry file exists.
        if not os.path.isfile(self.photfile):
            raise IOError('Photometry file not found: {}'.format(self.photfile))

        # The systematics file.
        if self.fixed:
            prefix = 'sys{}_vmag_'.format(self.aper)
        else:
            prefix = 'sys{}_'.format(self.aper)        
        
        if self.sysfile is None:
            head, tail = os.path.split(self.photfile)
            tail = prefix + tail.rsplit('_')[-1]
            self.sysfile = os.path.join(head, tail)
        
        # Check of the systematics file exists.
        if os.path.isfile(self.sysfile):
            raise IOError('Systematics file already exists: {}'.format(self.sysfile))

        # Read the required header data.
        self.stars, settings, lstmin, lstmax = _read_header(self.photfile, self.camgrid, self.skygrid)
        
        settings['outer_maxiter'] = self.outer_maxiter
        settings['inner_maxiter'] = self.inner_maxiter
        settings['dtol'] = self.dtol
        
        # The spatial calculation statistics.
        self.spatial = dict()
        self.spatial['n'] = np.unique(self.stars['n'])
        self.spatial['niter'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        self.spatial['chisq'] = np.full(self.camgrid.ny+2, fill_value=np.nan)
        self.spatial['npoints'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        self.spatial['npars'] = np.zeros(self.camgrid.ny+2, dtype='uint32')
        
        # The temporal calculation statistics.
        self.temporal = dict()
        self.temporal['q'] = np.unique(self.stars['q'])
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
        self.trans['gridtype'] = 'polar'
        self.trans['num_k'] = self.camgrid.nx
        self.trans['num_n'] = self.camgrid.ny
        self.trans['nobs'] = np.zeros((self.camgrid.nx+2, self.camgrid.ny+2), dtype='uint32')
        self.trans['trans'] = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        
        # The intrapixel amplitudes.
        self.intrapix = dict()
        self.intrapix['gridtype'] = 'polar'
        self.intrapix['num_l'] = self.ipxgrid.nx
        self.intrapix['num_n'] = self.ipxgrid.ny
        self.intrapix['nobs'] = np.zeros((self.ipxgrid.nx+2, self.ipxgrid.ny+2), dtype='uint32')
        self.intrapix['amplitudes'] = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan)
        
        # The cloud corrections.
        self.clouds = dict()
        self.clouds['gridtype'] = 'healpix'
        self.clouds['num_q'] = self.skygrid.nside
        self.clouds['lstmin'] = lstmin
        self.clouds['lstmax'] = lstmax
        self.clouds['lstlen'] = lstmax - lstmin + 1
        self.clouds['nobs'] = np.zeros((self.skygrid.npix, lstmax - lstmin + 1), dtype='uint32')
        self.clouds['clouds'] = np.full((self.skygrid.npix, lstmax - lstmin + 1), fill_value=np.nan)
        self.clouds['sigma'] = np.zeros((self.skygrid.npix, lstmax - lstmin + 1))
        
        self.got_sky = False
        
        # Perform the coarse decorrelation.
        print 'Performing coarse decorrelation for aperture {} of: {}'.format(self.aper, self.photfile)
        print 'Writing results to: {}'.format(self.sysfile) 
        
        for niter in range(self.outer_maxiter):
        
            print 'Iteration {} out of {}:'.format(niter + 1, self.outer_maxiter)
            
            print '    Calculating spatial systematics...'
            self._spatial()
            
            print '    Calculating temporal systematics...'
            self._temporal()
        
        # Write the results to file.
        io.write_calibration(self.sysfile, settings, self.spatial, self.temporal, self.magnitudes, self.trans, self.intrapix, self.clouds)
        
        return self.sysfile

class ApplyDecorrelation(object):

    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        f = io.SysFile(self.sysfile)
    
        self.magnitudes = f.read_magnitudes()    
        self.camgrid, self.trans = f.read_trans()
        self.ipxgrid, self.intrapix = f.read_intrapix()  
        self.skygrid, self.clouds = f.read_clouds()        
        
        return
        
    def get_magnitudes(self, ascc):

        mask = self.magnitudes['ascc'] == ascc
        mag = self.magnitudes['mag'][mask]
        sigma = self.magnitudes['sigma'][mask]
        nobs = self.magnitudes['nobs'][mask]
        
        return mag, sigma, nobs        
        
    def get_transmission(self, ra, dec, lst):
        
        ha = np.mod(lst*15 - ra, 360.)
        dec = np.repeat(dec, len(lst))

        k, n = self.camgrid.radec2idx(ha, dec)

        trans = self.trans['trans'][k,n]
        nobs = self.trans['nobs'][k,n]
        
        return trans, nobs
        
    def get_intrapix(self, ra, dec, lst, x, y):
        
        ha = np.mod(lst*15 - ra, 360.)
        dec = np.repeat(dec, len(lst))

        sinx, cosx = np.sin(2*np.pi*x), np.cos(2*np.pi*x)
        siny, cosy = np.sin(2*np.pi*y), np.cos(2*np.pi*y)
        b_mat = np.column_stack([sinx, cosx, siny, cosy])

        l, n = self.ipxgrid.radec2idx(ha, dec)
        
        ipx = np.sum(b_mat*self.intrapix['amplitudes'][l,n], axis=1)
        nobs = self.intrapix['nobs'][l,n] 
        
        return ipx, nobs
        
    def get_clouds(self, ra, dec, lstseq):
      
        q = self.skygrid.radec2idx(ra, dec)
        t = lstseq - self.clouds['lstmin']
        
        clouds = self.clouds['clouds'][q,t]
        sigma = self.clouds['sigma'][q,t]
        nobs = self.clouds['nobs'][q,t]
        
        return clouds, sigma, nobs
        
    def get_systematics(self, ascc, ra, dec, lstseq, lst, x, y):

        flag = np.zeros(len(lstseq), dtype='uint8')
        
        mag, sigma, nobs = self.get_magnitudes(ascc)        
        
        trans, nobs = self.get_transmission(ra, dec, lst)
        flag = np.where(nobs < 25, flag+2, flag)
        
        ipx, nobs = self.get_intrapix(ra, dec, lst, x, y)
        flag = np.where(nobs < 25, flag+4, flag)
        
        clouds, sigma, nobs = self.get_clouds(ra, dec, lstseq)
        flag = np.where(nobs < 25, flag+8, flag)
        flag = np.where(sigma > .05, flag+16, flag)
        
        systematics = trans + ipx + clouds
        flag = np.where(np.isnan(systematics), flag + 1, flag)

        return mag, trans, ipx, clouds, flag 

    def __call__(self, photfile, aper, redfile=None):
        
        # Check if the photometry file exists.
        if not os.path.isfile(photfile):
            raise IOError('Photometry file not found: {}'.format(photfile))
         
        # The reduced lightcurves file.
        if 'vmag' in self.sysfile:
            prefix = 'red{}_vmag_'.format(aper)
        else:
            prefix = 'red{}_'.format(aper)
        
        if redfile is None:
            head, tail = os.path.split(photfile)
            tail = prefix + tail.rsplit('_')[-1]
            redfile = os.path.join(head, tail)
        
        # Check if the reduced lightcurves file exists.
        if os.path.isfile(redfile):
            raise IOError('Reduced lightcurves file already exists: {}'.format(redfile))
        
        print 'Reading corrections from: {}'.format(self.sysfile) 
        print 'Applying corrections to aperture {} of file: {}'.format(aper, photfile)
        print 'Writing results to: {}'.format(redfile)
        
        # Read the raw photometry.
        f = io.PhotFile(photfile)
    
        settings = f.read_header()
    
        if f.get_siteid() == 'LP':
            
            stars = f.read_stars(['ascc', 'ra', 'dec', 'vmag', 'bmag', 'spectype'])
            
            stars['nobs'] = np.zeros(len(stars['ascc']), dtype='uint32')  
            stars['lstseqmin'] = np.zeros(len(stars['ascc']), dtype='uint32')     
            stars['lstseqmax'] = np.zeros(len(stars['ascc']), dtype='uint32')
            
        else:
            
            stars = f.read_stars()
            station = f.read_station()
        
        # Fields and datatypes for the binned lightcurves.
        lightcurves = dict()
        
        if f.get_siteid() == 'LP':
            names = ['lstseq', 'nobs', 'lst', 'jdmid', 'x', 'y', 'sky', 'esky', 'mag{}'.format(aper), 'emag{}'.format(aper), 'trans{}'.format(aper), 'etrans{}'.format(aper), 'clouds{}'.format(aper), 'eclouds{}'.format(aper)]
            formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']            
        else:
            names = ['lstseq', 'nobs', 'lst', 'jd', 'exptime', 'x', 'y', 'sky', 'esky', 'mag{}'.format(aper), 'emag{}'.format(aper), 'trans{}'.format(aper), 'etrans{}'.format(aper), 'clouds{}'.format(aper), 'eclouds{}'.format(aper)]
            formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']
        
        for i in range(len(stars['ascc'])):                    
                
            # Read the lightcurve.
            lc = f.read_lightcurves(ascc=stars['ascc'][i])
            
            # Remove flagged data.
            if f.get_siteid() == 'LP':
                mask = (lc['flux{}'.format(aper)] > 0) & (lc['eflux{}'.format(aper)] > 0) & (lc['flag'] < 1) & (lc['sky'] > 0)
            else:
                mask = (lc['aflag'] == 0) & (lc['pflag'] == 0)
           
            if np.all(~mask):
                stars['nobs'][i] = 0
                continue   
            
            lc = lc[mask]
    
            # Unpack the data.
            lstseq = lc['lstseq']
            x, y = lc['x'], lc['y']
            flux, eflux = lc['flux{}'.format(aper)], lc['eflux{}'.format(aper)]
            sky = lc['sky']
            
            if f.get_siteid() == 'LP':
                
                jd = lc['jdmid']
                lst = lc['lst']
                
                x, y = x.astype('float64'), y.astype('float64')
                
            else:
                
                idx = np.searchsorted(station['lstseq'], lstseq) 
                jd = station['jd'][idx]
                lst = station['lst'][idx]
                exptime = station['exptime'][idx]
                
                flux, eflux = flux/exptime, eflux/exptime
                
            # Compute the corrected lightcurve.
            mag0, trans, ipx, clouds, cflag = self.get_systematics(stars['ascc'][i], stars['ra'][i], stars['dec'][i], lstseq, lst, x, y)
            
            mag, emag = misc.flux2mag(flux, eflux)
            mag = mag - trans - ipx - clouds  
            
            if f.get_siteid() == 'LP':
                mag = mag - mag0
                trans = trans + mag0
                
            # Remove flagged data.
            mask = (cflag < 1)         
            
            if np.all(~mask):
                stars['nobs'][i] = 0
                continue 
            
            jd = jd[mask]
            lst = lst[mask]
            lstseq = lstseq[mask]
            x = x[mask]
            y = y[mask]
            sky = sky[mask]
            mag = mag[mask]
            emag = emag[mask]
            trans = trans[mask]
            ipx = ipx[mask]
            clouds = clouds[mask] 
            
            if f.get_siteid() != 'LP':
                exptime = exptime[mask]
                
            # Compute the final binned lightcurve.
            lstseq, binidx, nobs = np.unique(lstseq//50, return_inverse=True, return_counts=True)        
            
            lc_bin = np.recarray(len(lstseq), names=names, formats=formats)        
            
            lc_bin['lstseq'] = lstseq
            lc_bin['nobs'] = nobs # Number of raw points used for each binned point.        
            
            if f.get_siteid() == 'LP':
                lc_bin['jdmid'] = statistics.idxstats(binidx, jd, statistic='mean')
                lc_bin['lst'] = statistics.idxstats(binidx, lst, statistic='mean')
            else:
                lc_bin['jd'] = statistics.idxstats(binidx, jd, statistic='mean')
                lc_bin['lst'] = statistics.idxstats(binidx, lst, statistic='mean')
                lc_bin['exptime'] = statistics.idxstats(binidx, exptime, statistic='sum')            

            lc_bin['x'] = statistics.idxstats(binidx, x, statistic='mean')
            lc_bin['y'] = statistics.idxstats(binidx, y, statistic='mean')        
            
            lc_bin['mag{}'.format(aper)] = statistics.idxstats(binidx, mag, statistic='mean')
            lc_bin['emag{}'.format(aper)] = statistics.idxstats(binidx, mag, statistic='std')/np.sqrt(nobs)
            lc_bin['sky'] = statistics.idxstats(binidx, sky, statistic='mean')
            lc_bin['esky'] = statistics.idxstats(binidx, sky, statistic='std')/np.sqrt(nobs)        
            
            lc_bin['trans{}'.format(aper)] = statistics.idxstats(binidx, trans, statistic='mean')
            lc_bin['etrans{}'.format(aper)] = statistics.idxstats(binidx, trans, statistic='std')/np.sqrt(nobs)
            lc_bin['clouds{}'.format(aper)] = statistics.idxstats(binidx, clouds, statistic='mean')
            lc_bin['eclouds{}'.format(aper)] = statistics.idxstats(binidx, clouds, statistic='std')/np.sqrt(nobs)            
            
            lightcurves[stars['ascc'][i]] = lc_bin
        
            stars['nobs'][i] = len(lstseq)
            if f.get_siteid() == 'LP':
                stars['lstseqmin'][i] = lstseq[0]
                stars['lstseqmax'][i] = lstseq[-1] 
            else:
                stars['lstsqmin'][i] = lstseq[0]
                stars['lstsqmax'][i] = lstseq[-1]
               
        io.write_reduced(redfile, settings, stars, lightcurves, f.get_siteid())
            
        return redfile

def run_calibration(date, part, aper, cameras, source, dest):
    
    cal = CoarseDecorrelation()
    
    for cam in cameras:
        
        photfile = baseline(date, part, cam, source, dest)
        
        if photfile is not None:        
        
            sysfile = cal(photfile, aper)    
            
            f = ApplyDecorrelation(sysfile)
            redfile = f(photfile, aper)
    
    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Perform the coarse decorrelation on a baseline.')
    parser.add_argument('date', type=str,
                        help='a date of the form YYYYMM')
    parser.add_argument('part', type=str, choices=['A','B','C'],
                        help='letter indicating which baseline to create') 
    parser.add_argument('dest', type=str,
                        help='Location where the products will be written, e.g. /data3/talens/2017Q1. If the path does not exist it will be created. Subdirectories are generated automatically.')
    parser.add_argument('-c', '--cam', type=str, nargs='+', default=['LPN', 'LPE', 'LPS', 'LPW', 'LPC', 'LSN', 'LSE', 'LSS', 'LSW', 'LSC'],
                        help ='the camera(s) to perform the combination for', dest='cameras')                        
    parser.add_argument('-a', '--aper', type=int, choices=[0,1], default=0,
                        help ='the aperture to perform the coarse decorrelation on', dest='aper')
    parser.add_argument('-d', '--data', type=str, default='/data2/mascara/LaPalma',
                        help='Location of the raw data.', dest='source')
    args = parser.parse_args()
    
    run_calibration(args.date, args.part, args.aper, args.cameras, args.source, args.dest)
    
    return

if __name__ == '__main__':
    
    main() 
    

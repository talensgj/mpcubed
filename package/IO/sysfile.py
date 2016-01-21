#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from ..coordinates import grids

class SysFile(object):
    """ Read data from systematics files.
    
    Attributes:
        sysfile (str): The full path to the file.
        
    """
    
    
    def __init__(self, sysfile):
        """ Initialize a reader of systematics files.
        
        Args:
            sysfile (str): The full path to the file.
            
        """
        
        self.sysfile = sysfile
        
        return
        
    def read_header(self):
        """ Read all header attributes.
        
        Returns:
            data: A list of attribute (key, value) pairs.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            data = f['header'].attrs.items()
        
        return data
        
    def read_pointing(self):
        """ Read the pointing associated with the systematics.
        
        Returns:
            alt0 (float): The pointing altitude in degrees.
            az0 (float): The pointing azimuth in degrees. 
            th0 (float): The pointing orientation in degrees.
            x0 (float): The pointing center in pixels.
            y0 (float): The pointing center in pixels.
            
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            alt0 = f['header'].attrs['alt0']
            az0 = f['header'].attrs['az0']
            th0 = f['header'].attrs['th0']
            x0 = f['header'].attrs['x0']
            y0 = f['header'].attrs['y0']
        
        return alt0, az0, th0, x0, y0
        
    def read_statistics(self, mode):
        """ Read the quality statistics from the header.
        
        Args:
            mode (str): Either 'spatial' or 'temporal' depending
                on which quality parameters should be read.
                
        Returns:
            niter (int): The number of iterations performed.
            chisq (float): The chi-squared value achieved.
            npoints (int): The number of datapoints in the fit.
            npars (int): The number of model parameters used.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['header/' + mode]
            niter = grp['niter'].value
            chisq = grp['chisq'].value
            npoints = grp['npoints'].value
            npars = grp['npars'].value
            
        return niter, chisq, npoints, npars
        
    def read_magnitudes(self):
        """ Read the fitted magnitudes.
        
        Returns:
            ascc (str): The ascc numbers of the stars.
            vmag (float): The catalogue magnitude of the stars.
            mag (float): The fitted magnitude of the stars.
            sigma (float): The extra error term. Set to None if absent.
            nobs (int): The number of datapoints used to compute each magnitude.
            
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/magnitudes']
            ascc = grp['ascc'].value
            vmag = grp['vmag'].value
            mag = grp['mag'].value
            nobs = grp['nobs'].value
            
            try: sigma = grp['sigma'].value
            except: sigma = None
        
        return ascc, vmag, mag, sigma, nobs
        
    def read_trans(self):
        """ Read the fitted transmission.
        
        Returns:
            pg: A polar grid instance corresponding to the transmision.
            trans (float): The transmission values on the grid.
            nobs (int): The number of datapoints used to compute each
                transmission.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/trans']
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            trans = grp['trans'].value
            nobs = grp['nobs'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        trans = pg.values2grid(idx1, idx2, trans, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)

        return pg, trans, nobs
        
    def read_intrapix(self):
        """ Read the fitted intrapixel variations.
        
        Returns:
            pg: A polar grid instance corresponding to the intrapixel
                variations.
            sinx (float): The amplitudes of the sinx term on the grid. 
            cosx (float): The amplitudes of the cosx term on the grid. 
            siny (float): The amplitudes of the siny term on the grid. 
            cosy (float): The amplitudes of the cosy term on the grid. 
            nobs (float): The number of datapoints used to compute the
                amplitudes.
                
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/intrapix']
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            sinx = grp['sinx'].value
            cosx = grp['cosx'].value
            siny = grp['siny'].value
            cosy = grp['cosy'].value
            nobs = grp['nobs'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        sinx = pg.values2grid(idx1, idx2, sinx, np.nan)
        cosx = pg.values2grid(idx1, idx2, cosx, np.nan)
        siny = pg.values2grid(idx1, idx2, siny, np.nan)
        cosy = pg.values2grid(idx1, idx2, cosy, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)
        
        return pg, sinx, cosx, siny, cosy, nobs
    
    def read_clouds(self):
        """ Read the fitted clouds.
        
        Returns:
            hg: A halpix grid instance corresponding to the clouds.
            clouds: The cloud values computed on the grid at each instant.
            sigma: The extra error term. Set to None if absent.
            nobs: The number of datapoints used to compute the clouds.
            lstmin: The lowest lstseq for which clouds were computed.
            lstmax: The highest lstseq for which clouds were computed. 
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/clouds']
            idx = grp['idx'].value
            lstseq = grp['lstseq'].value
            clouds = grp['clouds'].value
            nobs = grp['nobs'].value
            
            try: sigma = grp['sigma'].value
            except: sigma = None
            
            nx = grp.attrs['nx']
            lstmin = grp.attrs['lstmin']
            lstmax = grp.attrs['lstmax']
            lstlen = grp.attrs['lstlen']
    
        lstseq = lstseq - lstmin
    
        hg = grids.HealpixGrid(nx)
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = clouds
        clouds = tmp
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = nobs
        nobs = tmp
        
        if sigma is not None:
            tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
            tmp[idx, lstseq] = sigma
            sigma = tmp

        return hg, clouds, sigma, nobs, lstmin, lstmax

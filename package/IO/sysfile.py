#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from ..coordinates import grids

class SysFile():
    
    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        return
        
    def read_header(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            data = f['header'].attrs.items()
        
        return data
        
    def read_pointing(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            alt0 = f['header'].attrs['alt0']
            az0 = f['header'].attrs['az0']
            th0 = f['header'].attrs['th0']
            x0 = f['header'].attrs['x0']
            y0 = f['header'].attrs['y0']
        
        return alt0, az0, th0, x0, y0
        
    def read_statistics(self, mode):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['header/' + mode]
            niter = grp['niter'].value
            chisq = grp['chisq'].value
            npoints = grp['npoints'].value
            npars = grp['npars'].value
            
        return niter, chisq, npoints, npars
        
    def read_magnitudes(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/magnitudes']
            ascc = grp['ascc'].value
            vmag = grp['vmag'].value
            mag = grp['mag'].value
            nobs = gr['nobs'].value
            
            try: sigma = grp['sigma'].value
            except: sigma = None
        
        return ascc, vmag, mag, sigma, nobs
        
    def read_trans(self):
        
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

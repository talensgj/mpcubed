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
            
            ascc = f['data/magnitudes/ascc'].value
            vmag = f['data/magnitudes/vmag'].value
            mag = f['data/magnitudes/mag'].value
            nobs = f['data/magnitudes/nobs'].value
        
        return ascc, vmag, mag, nobs
        
    def read_trans(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx1 = f['data/trans/idx1'].value
            idx2 = f['data/trans/idx2'].value
            trans = f['data/trans/trans'].value
            nobs = f['data/trans/nobs'].value
            
            nx = f['data/trans'].attrs['nx']
            ny = f['data/trans'].attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        trans = pg.values2grid(idx1, idx2, trans, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)

        return pg, trans, nobs
        
    def read_intrapix(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx1 = f['data/intrapix/idx1'].value
            idx2 = f['data/intrapix/idx2'].value
            sinx = f['data/intrapix/sinx'].value
            cosx = f['data/intrapix/cosx'].value
            siny = f['data/intrapix/siny'].value
            cosy = f['data/intrapix/cosy'].value
            nobs = f['data/intrapix/nobs'].value
            
            nx = f['data/intrapix'].attrs['nx']
            ny = f['data/intrapix'].attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        sinx = pg.values2grid(idx1, idx2, sinx, np.nan)
        cosx = pg.values2grid(idx1, idx2, cosx, np.nan)
        siny = pg.values2grid(idx1, idx2, siny, np.nan)
        cosy = pg.values2grid(idx1, idx2, cosy, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)
        
        return pg, sinx, cosx, siny, cosy, nobs
    
    def read_clouds(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx = f['data/clouds/idx'].value
            lstseq = f['data/clouds/lstseq'].value
            clouds = f['data/clouds/clouds'].value
            nobs = f['data/clouds/nobs'].value
            
            nx = f['data/clouds'].attrs['nx']
            lstmin = f['data/clouds'].attrs['lstmin']
            lstmax = f['data/clouds'].attrs['lstmax']
            lstlen = f['data/clouds'].attrs['lstlen']
    
        lstseq = lstseq - lstmin
    
        hg = grids.HealpixGrid(nx)
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = clouds
        clouds = tmp
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = nobs
        nobs = tmp

        return hg, clouds, nobs, lstmin, lstmax

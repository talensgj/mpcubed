#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from sysrem import sysrem

class SkyTransmission():
    
    def __init__(self):
        self.grid = 'healpix'
        self.nx = 8
        self.ny = 0
        self.margin = 0
        
    def calculate(self, fLC, red):
        
        self.fLC = fLC
        self.red = red
        
        # Read the stellar header information.
        with h5py.File(self.fLC, 'r') as f:
        
            hdr = f['table_header']
            ascc = hdr['ascc'].value
            vmag = hdr['vmag'].value
            ra = hdr['ra'].value
            dec = hdr['dec'].value
            nobs = hdr['nobs'].value.astype('int')
        
        if self.grid == 'healpix':
            
            self.process_healpix()
             
        elif self.grid in ['polar', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid
        
    def process_healpix(self):
        grid = HealpixGrid(self.nx)
        bins, binnum = grid.find_gridpoints(ra, dec, compact=True)
        starcount = index_statistics(binnum, binnum, statistic='count')
        
        for ind in range(len(bins)):
            
            here = (binnum == ind)
            ascc = self.ascc[here]
            vmag = self.vmag[here]
            nobs = self.nobs[here]
            
            nstars = len(ascc)
            ndata = np.sum(nobs)
            select = np.append(0, np.cumsum(nobs))
            
            lstidx = np.zeros(ndata)
            sky = np.zeros(ndata)
            flags1 = np.zeros(ndata)
            
            cflux0 = np.zeros(ndata)
            ecflux0 = np.zeros(ndata)
            flags2 = np.zeros(ndata) 
            
            with h5py.File(self.fLC, 'r') as f, h5py.File(self.red, 'r') as g:
            
                lc = f['data']
                rc = f['data']
            
                for i in range(nstars):
                    
                    lstidx[select[i]:select[i+1]] = lc[ascc_tmp[i]]['lstidx']
                    sky[select[i]:select[i+1]] = lc[ascc_tmp[i]]['sky']
                    flags1[select[i]:select[i+1]] = lc[ascc_tmp[i]]['flag']
                    
                    cflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['cflux0']
                    ecflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['ecflux0']
                    flags2[select[i]:select[i+1]] = rc[ascc_tmp[i]]['flags']
                    
            star_id = np.repeat(np.arange(nstars), nobs)
            
            here = (flags1<1)#&(flags2<1)&(cflux0>0)&(sky>0)&(ecflux0>0)
            lstidx = lstidx[here]
            cflux0 = cflux0[here]
            ecflux0 = ecflux0[here]
            star_id = star_id[here]
            
            times, lstidx = np.unique(lstidx, return_inverse=True)
            # Also unique for star_id 
            
            count = index_statistics(lstidx, lstidx, statistic='count')
            
            a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(lstidx, star_id, cflux0, ecflux0, a2 = (1e7)*10**(vmag/-2.5))

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from sysrem import sysrem

class CameraTransmission():
    
    def __init__(self, grid='polar', nx=13500, ny=720, margin=0):
    
        self.grid = grid
        self.nx = nx
        self.ny = ny
        self.margin = margin

        return 
        
    def _read_data(self, ascc, nobs):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lst = np.zeros(ndata)
        flux0 = np.zeros(ndata)
        eflux0 = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags = np.zeros(ndata)
        
        with h5py.File(self.fLCfile, 'r') as f:
        
            lc = f['data']
        
            for i in range(nstars):
                
                lst[select[i]:select[i+1]] = lc[ascc[i]]['lst']
                flux0[select[i]:select[i+1]] = lc[ascc[i]]['flux0']
                eflux0[select[i]:select[i+1]] = index_statistics(lc[ascc[i]]['lstidx']//50, lc[ascc[i]]['flux0'], statistic='std', keeplength=True)
                sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
                flags[select[i]:select[i+1]] = lc[ascc[i]]['flag']
   
        return lst, flux0, eflux0, sky, flags

    def calculate(self, fLCfile, camfile=None):
        
        self.fLCfile = fLCfile
        
        # Filename of output file.
        if camfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'cam_'+tail.rsplit('_')[-1]
            self.camfile = os.path.join(head, tail)
        else:
            self.camfile = camfile
        
        # Read the stellar header information.
        with h5py.File(self.fLCfile, 'r') as f:
            
            hdr = f['header_table']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
        
        # Calculate the transmission.
        if self.grid == 'polar':
            self._calculate_polar()
             
        elif self.grid in ['healpix', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid

        return
        
    def _calculate_polar(self):
        
        # Polar grid instance.
        pg = PolarGrid(self.nx, self.ny)
        
        # Assign stars to declination bins and count the number of stars in each bin.
        decidx, decuni = pg.find_decidx(self.dec, compact=True)
        starcount = np.bincount(decuni)
        
        # Create arrays.
        niter = np.zeros(len(decidx), dtype='int')
        chisq = np.zeros(len(decidx), dtype='float')
        npoints = np.zeros(len(decidx), dtype='int')
        npars = np.zeros(len(decidx), dtype='int')
        
        for ind in range(len(decidx)):
            
            # Select stars in the current sky bin.
            here = (decuni == ind)
            ascc = self.ascc[here]
            ra = self.ra[here]
            vmag = self.vmag[here]
            nobs = self.nobs[here]
            
            # Read data for these stars.
            lst, flux0, eflux0, sky, flags = self._read_data(ascc, nobs)
            
            # Create the haidx and staridx. 
            ha = np.mod(lst*15.-np.repeat(ra,nobs), 360.)
            haidx = pg.find_raidx(ha)
            staridx = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = (flux0 > 0)&(eflux0 > 0)&(sky > 0)&(flags < 1)
            flux0 = flux0[here]
            eflux0 = eflux0[here]
            haidx = haidx[here]
            staridx = staridx[here]
            
            if len(flux0) == 0: continue
            
            # Make the haidx ascending from 0 and count the number of datapoints at each haidx.
            haidx, hauni = np.unique(haidx, return_inverse=True)
            pointcount = np.bincount(hauni)
            
            # Compute the camera transmission curve.
            flux, camtrans, niter[ind], chisq[ind], npoints[ind], npars[ind] = sysrem(staridx, hauni, flux0, eflux0)
            
            with h5py.File(self.camfile) as f:
                
                grp = f.create_group('data/%i'%decidx[ind])
                grp.create_dataset('haidx', data=haidx)
                grp.create_dataset('pointcount', data=pointcount)
                grp.create_dataset('camtrans', data=camtrans)
        
        with h5py.File(self.camfile) as f:
            
            grp = f.create_group('header')
            
            grp.attrs['grid'] = self.grid
            grp.attrs['nx'] = self.nx
            grp.attrs['ny'] = self.ny
            grp.attrs['margin'] = self.margin
            
            grp.create_dataset('decidx', data = decidx)
            grp.create_dataset('starcount', data = starcount)
            grp.create_dataset('niter', data = niter)
            grp.create_dataset('chisq', data = chisq)
            grp.create_dataset('npoints', data = npoints)
            grp.create_dataset('npars', data = npars)
        
        return
        
class CameraFile():
    
    def __init__(self, camfile):
        
        self.camfile = camfile
        
        return
        
    def _read_camfile(self):
        
        with h5py.File(self.camfile, 'r') as f:
            
            self.grid = f['header'].attrs['grid']
            self.nx = f['header'].attrs['nx']
            self.ny = f['header'].attrs['ny']
            self.margin = f['header'].attrs['margin']
            
            self.starcount = np.full(self.ny+2, fill_value=0)
            self.camtrans = np.full((self.nx+2, self.ny+2), fill_value=np.nan)
            self.pointcount = np.full((self.nx+2, self.ny+2), fill_value=0)
            
            decidx = f['header/decidx'].value
            self.starcount[decidx] = f['header/starcount'].value
            
            for ind in decidx:
                
                try:
                    data = f['data/%i'%ind]
                except:
                    pass
                else:
                    haidx = data['haidx'].value
                    
                    self.camtrans[haidx, ind] = data['camtrans'].value
                    self.pointcount[haidx, ind] = data['pointcount'].value
                    
        return
        
    def visualize(self, wrap=False):
        
        self.wrap = wrap
        self._read_camfile()
            
        if self.grid == 'polar':
            self._visualize_polar()
            
        elif self.grid in ['healpix', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid
            exit()
            
        return
    
    def _visualize_polar(self):
        
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from viridis import viridis 
        
        rcParams['xtick.labelsize'] = 'large'
        rcParams['ytick.labelsize'] = 'large'
        rcParams['axes.labelsize'] = 'x-large'
        rcParams['image.interpolation'] = 'none'
        rcParams['image.origin'] = 'lower'
        
        array = self.camtrans[1:-1,1:-1]
        array = array/np.nanmean(array, axis=0, keepdims=True)
        if self.wrap: array = np.roll(array, self.nx/2, axis=0)
        
        xlim, ylim = np.where(np.isfinite(array))
        
        ax = plt.subplot(111)
        plt.imshow(array.T/np.nanmax(array), aspect='auto', cmap=viridis, vmin=0, vmax=1)
        cb = plt.colorbar()
        plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        plt.xlabel('HA [idx]')
        plt.ylabel('Dec [idx]')
        cb.set_label('Camera')
        plt.show()
 
        return
        
    def correct(self, fLCfile, redfile=None):
        
        self.fLCfile = fLCfile
        
        # Filename of output file.
        if redfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'red_'+tail.rsplit('_')[-1]
            self.redfile = os.path.join(head, tail)
        else:
            self.redfile = redfile
        
        self._read_camfile()
            
        if self.grid == 'polar':

            self._correct_polar()
            
        elif self.grid in ['healpix', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid
            exit()
        
        return
        
    def _correct_polar(self):
            
        # Polar grid instance.
        pg = PolarGrid(self.nx, self.ny)
        
        with h5py.File(self.fLCfile, 'r') as f, h5py.File(self.redfile, 'w-') as g:
            
            ascc = f['header_table/ascc'].value
            ra = f['header_table/ra'].value
            dec = f['header_table/dec'].value
            
            decidx = pg.find_decidx(dec)
            
            # Add header to resulting file.
            grp = g.create_group('header')
            grp.attrs['camfile'] = self.camfile
            
            for i in range(len(ascc)):
                
                lc = f['data/'+ascc[i]]
                
                # Find haidx for the lightcurve.
                ha = np.mod(lc['lst']*15.-ra[i], 360)
                haidx = pg.find_raidx(ha)
            
                # Find the camera transmission curve and correct flux and eflux.
                camtrans0 = self.camtrans[haidx, decidx[i]]
                cflux0 = lc['flux0']/camtrans0
                ecflux0 = lc['eflux0']/camtrans0
                
                # Add also the std of 50 points used in calculating the tranmission.
                sflux0 = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
                scflux0 = sflux0/camtrans0
        
                # Add flags.
                flags = np.where(np.isnan(camtrans0), 1, 0) # Flag data where no transmission coefficient exists.
                flags = flags + np.where(self.pointcount[haidx, decidx[i]]<=5, 2, 0) # Flag data where less than 5 datapoints were used to compute the transmission.
        
                # Combine the reduced data in a record array.
                arlist = [camtrans0, cflux0, ecflux0, sflux0, scflux0, flags]
                names = ['camtrans0', 'cflux0', 'ecflux0', 'sflux0', 'scflux0', 'flags']
                record = np.rec.fromarrays(arlist, names=names)
                
                # Write the result to file.
                g.create_dataset('data/'+ascc[i], data=record)

        return

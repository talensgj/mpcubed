#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import PolarGrid
from sysrem import intrarem

class IntraPixel():
    
    def __init__(self, nx_cam=13500, nx_ipx=270, ny=720):
    
        self.nx_cam = nx_cam
        self.nx_ipx = nx_ipx
        self.ny = ny

        return 
        
    def _read_data(self, ascc, nobs):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lst = np.zeros(ndata)
        y = np.zeros(ndata)
        flux0 = np.zeros(ndata)
        eflux0 = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags = np.zeros(ndata)
        
        with h5py.File(self.fLCfile, 'r') as f:
        
            lc = f['data']
        
            for i in range(nstars):
                
                lst[select[i]:select[i+1]] = lc[ascc[i]]['lst']
                y[select[i]:select[i+1]] = lc[ascc[i]]['y']
                flux0[select[i]:select[i+1]] = lc[ascc[i]]['flux0']
                eflux0[select[i]:select[i+1]] = index_statistics(lc[ascc[i]]['lstidx']//50, lc[ascc[i]]['flux0'], statistic='std', keeplength=True)
                sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
                flags[select[i]:select[i+1]] = lc[ascc[i]]['flag']
   
        return lst, y, flux0, eflux0, sky, flags
        
    def calculate(self, fLCfile, camfile=None):
    
        self.fLCfile = fLCfile
        
        print 'Calculating camera transmission and intrapixel variations for:', os.path.split(self.fLCfile)[1]
        
        # Filename of output file.
        if camfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'camip_'+tail.rsplit('_')[-1]
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
            
        # Polar grid instances.
        pgcam = PolarGrid(self.nx_cam, self.ny)
        pgipx = PolarGrid(self.nx_ipx, self.ny)
        
        # Assign stars to declination bins and count the number of stars in each bin.
        decidx, decuni = pgcam.find_decidx(self.dec, compact=True)
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
            lst, y, flux0, eflux0, sky, flags = self._read_data(ascc, nobs)
            
            # Create the haidx and staridx. 
            ha = np.mod(lst*15.-np.repeat(ra,nobs), 360.)
            haidx_cam = pgcam.find_raidx(ha)
            haidx_ipx = pgipx.find_raidx(ha)
            staridx = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = (flux0 > 0)&(eflux0 > 0)&(sky > 0)&(flags < 1)
            y = y[here]
            flux0 = flux0[here]
            eflux0 = eflux0[here]
            haidx_cam = haidx_cam[here]
            haidx_ipx = haidx_ipx[here]
            staridx = staridx[here]

            if len(flux0) == 0: continue
            
            # Make the haidx ascending from 0 and count the number of datapoints at each haidx.
            haidx_cam, hauni_cam = np.unique(haidx_cam, return_inverse=True)
            pointcount_cam = np.bincount(hauni_cam)
            
            haidx_ipx, hauni_ipx = np.unique(haidx_ipx, return_inverse=True)
            pointcount_ipx = np.bincount(hauni_ipx)
            
            # Compute the camera transmission curve and fit to the intrapixel variations.
            flux, camtrans, a, b, niter[ind], chisq[ind], npoints[ind], npars[ind] = intrarem(staridx, hauni_cam, hauni_ipx, y, flux0, eflux0)
    
            with h5py.File(self.camfile) as f:
                
                grp = f.create_group('data/%i'%decidx[ind])
                grp.create_dataset('haidx_cam', data=haidx_cam)
                grp.create_dataset('pointcount_cam', data=pointcount_cam)
                grp.create_dataset('camtrans', data=camtrans)
                grp.create_dataset('haidx_ipx', data=haidx_ipx)
                grp.create_dataset('pointcount_ipx', data=pointcount_ipx)
                grp.create_dataset('a', data=a)
                grp.create_dataset('b', data=b)

        with h5py.File(self.camfile) as f:

            grp = f.create_group('header')

            grp.attrs['nx_cam'] = self.nx_cam
            grp.attrs['nx_ipx'] = self.nx_ipx
            grp.attrs['ny'] = self.ny

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
            
            self.nx_cam = f['header'].attrs['nx_cam']
            self.nx_ipx = f['header'].attrs['nx_ipx']
            self.ny = f['header'].attrs['ny']
            
            self.starcount = np.full(self.ny+2, fill_value=0)
            
            self.camtrans = np.full((self.nx_cam+2, self.ny+2), fill_value=np.nan)
            self.pointcount_cam = np.full((self.nx_cam+2, self.ny+2), fill_value=0)
            
            self.a = np.full((self.nx_ipx+2, self.ny+2), fill_value=np.nan)
            self.b = np.full((self.nx_ipx+2, self.ny+2), fill_value=np.nan)
            self.pointcount_ipx = np.full((self.nx_ipx+2, self.ny+2), fill_value=0)
            
            decidx = f['header/decidx'].value
            self.starcount[decidx] = f['header/starcount'].value
            
            for ind in decidx:
                
                try:
                    data = f['data/%i'%ind]
                except:
                    pass
                else:
                    haidx_cam = data['haidx_cam'].value
                    haidx_ipx = data['haidx_ipx'].value
                    
                    self.camtrans[haidx_cam, ind] = data['camtrans'].value
                    self.pointcount_cam[haidx_cam, ind] = data['pointcount_cam'].value
                    
                    self.a[haidx_ipx, ind] = data['a'].value
                    self.b[haidx_ipx, ind] = data['b'].value
                    self.pointcount_ipx[haidx_ipx, ind] = data['pointcount_ipx'].value
                    
        return
    
    def visualize(self, wrap=False):
        
        self._read_camfile()
        
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
        if wrap: array = np.roll(array, self.nx_cam/2, axis=0)
        
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
        
        array = np.sqrt(self.a[1:-1,1:-1]**2+self.b[1:-1,1:-1]**2)
        if wrap: array = np.roll(array, self.nx_ipx/2, axis=0)
        
        xlim, ylim = np.where(np.isfinite(array))

        plt.imshow(array.T, aspect='auto', vmin=0, vmax=.1, cmap=viridis)
        cb = plt.colorbar()
        plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        plt.xlabel('HA [idx]')
        plt.ylabel('Dec [idx]')
        cb.set_label('Amplitude')
        plt.show()
        
        array = np.arctan2(self.b[1:-1,1:-1], self.a[1:-1,1:-1])
        if wrap: array = np.roll(array, self.nx_ipx/2, axis=0)

        xlim, ylim = np.where(np.isfinite(array))

        plt.imshow(array.T, aspect='auto', vmin=-np.pi, vmax=np.pi, cmap=viridis)
        cb = plt.colorbar()
        plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        plt.xlabel('HA [idx]')
        plt.ylabel('Dec [idx]')
        cb.set_label('Phase')
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
        
        # Polar grid instance.
        pgcam = PolarGrid(self.nx_cam, self.ny)
        pgipx = PolarGrid(self.nx_ipx, self.ny)
        
        with h5py.File(self.fLCfile, 'r') as f, h5py.File(self.redfile, 'w-') as g:
            
            ascc = f['header_table/ascc'].value
            ra = f['header_table/ra'].value
            dec = f['header_table/dec'].value
            
            decidx = pgcam.find_decidx(dec)
        
            # Add header to resulting file.
            grp = g.create_group('header')
            grp.attrs['camfile'] = self.camfile
    
            for i in range(len(ascc)):
                
                lc = f['data/'+ascc[i]]
        
                # Find haidx for the lightcurve.
                ha = np.mod(lc['lst']*15.-ra[i], 360)
                haidx_cam = pgcam.find_raidx(ha) 
                haidx_ipx = pgipx.find_raidx(ha)
                
                # Find the camera transmission curve and correct flux and eflux.
                camtrans0 = self.camtrans[haidx_cam, decidx[i]]
                cflux0 = lc['flux0']/camtrans0
                ecflux0 = lc['eflux0']/camtrans0
                
                # Find the intrapixel variations and correct cflux and ecflux.
                intrapix0 = self.a[haidx_ipx, decidx[i]]*np.sin(2*np.pi*lc['y']) + self.b[haidx_ipx, decidx[i]]*np.cos(2*np.pi*lc['y']) + 1
                ipcflux0 = cflux0/intrapix0
                eipcflux0 = ecflux0/intrapix0
        
                # Add also the std of 50 points used in calculating the tranmission.
                sflux0 = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
                scflux0 = sflux0/camtrans0
                sipcflux0 = scflux0/intrapix0
        
                # Add flags.
                flags = np.where(np.isnan(camtrans0), 1, 0)
                flags = flags + np.where(self.pointcount_cam[haidx_cam, decidx[i]]<=5, 2, 0)
                
                flags = flags + np.where(np.isnan(intrapix0), 4, 0)
                #flags = flags + np.where(self.pointcount_sky[haidx_ipx, decidx[i]]<=50, 8, 0)
                
                # Combine the reduced data in a record array.
                arlist = [camtrans0, cflux0, ecflux0, intrapix0, ipcflux0, eipcflux0, sflux0, scflux0, sipcflux0, flags]
                names = ['camtrans0', 'cflux0', 'ecflux0', 'intrapix0', 'ipcflux0', 'eipcflux0', 'sflux0', 'scflux0', 'sipcflux0', 'flags']
                record = np.rec.fromarrays(arlist, names=names)
                
                # Write the result to file.
                g.create_dataset('data/'+ascc[i], data=record)
                
        return
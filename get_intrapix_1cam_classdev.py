#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import PolarGrid
from intrarem import trans_intrapixel

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
        
        # Filename of output file.
        if camfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'camip_'+tail.rsplit('_')[-1]
            self.camfile = os.path.join(head, tail)
        else:
            self.camfile = camfile
        
        # Read the stellar header information.
        with h5py.File(self.fLCfile, 'r') as f:
            
            hdr = f['table_header']
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
            flux, camtrans, a, b, niter[ind], chisq[ind], npoints[ind], npars[ind] = trans_intrapixel(staridx, hauni_cam, hauni_ipx, y, flux0, eflux0, verbose=True)
    
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
    
    def visualize(self):
        
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
        array = array.T
        xlim, ylim = np.where(np.isfinite(array))

        plt.imshow(array, aspect='auto', vmin=0, vmax=1, cmap=viridis)
        cb = plt.colorbar()
        cb.set_label('camtrans')
        plt.xlabel('HA')
        plt.ylabel('Dec')
        plt.ylim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.xlim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        array = np.sqrt(self.a[1:-1,1:-1]**2+self.b[1:-1,1:-1]**2)
        array = array.T
        xlim, ylim = np.where(np.isfinite(array))

        plt.imshow(array, aspect='auto', vmin=0, vmax=.1, cmap=viridis)
        cb = plt.colorbar()
        cb.set_label('amplitude')
        plt.xlabel('HA')
        plt.ylabel('Dec')
        plt.ylim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.xlim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        array = np.arctan2(self.b[1:-1,1:-1], self.a[1:-1,1:-1])
        array = array.T
        xlim, ylim = np.where(np.isfinite(array))

        plt.imshow(array, aspect='auto', vmin=-np.pi, vmax=np.pi, cmap=viridis)
        cb = plt.colorbar()
        cb.set_label('phase')
        plt.xlabel('HA')
        plt.ylabel('Dec')
        plt.ylim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.xlim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        
        plt.tight_layout()
        plt.show()
        plt.close()
       
        return
    
    def correct(self, fLCfile, redfile=None):
        
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from viridis import viridis 
        
        rcParams['xtick.labelsize'] = 'large'
        rcParams['ytick.labelsize'] = 'large'
        rcParams['axes.labelsize'] = 'x-large'
        rcParams['image.interpolation'] = 'none'
        rcParams['image.origin'] = 'lower'
        
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
        
        with h5py.File(self.fLCfile, 'r') as f:
            
            ascc = f['table_header/ascc'].value
            ra = f['table_header/ra'].value
            dec = f['table_header/dec'].value
            
            decidx = pgcam.find_decidx(dec)
    
            here = ascc == '605365'
            ascc = ascc[here]
            ra = ra[here]
            dec = dec[here]
            decidx = decidx[here]
    
            for i in range(0, len(ascc), 100):
                
                lc = f['data/'+ascc[i]]
        
                ha = np.mod(lc['lst']*15.-ra[i], 360)
                haidx_cam = pgcam.find_raidx(ha) 
                haidx_ipx = pgipx.find_raidx(ha)
                
                camtrans0 = self.camtrans[haidx_cam, decidx[i]]
                intrapix0 = self.a[haidx_ipx, decidx[i]]*np.sin(2*np.pi*lc['y'])+self.b[haidx_ipx, decidx[i]]*np.cos(2*np.pi*lc['y'])+1
        
                norm = np.median(lc['flux0']/camtrans0)
        
                plt.figure(figsize=(16,8))
                ax = plt.subplot(311)
                plt.title('ASCC %s'%ascc[i])
                plt.plot(lc['jdmid'], lc['flux0'], '.')
                plt.plot(lc['jdmid'], camtrans0*norm, '.')
                plt.ylabel('flux0')
                ax1 = plt.subplot(312, sharex=ax)
                plt.plot(lc['jdmid'], lc['flux0']/(camtrans0*norm), '.')
                plt.plot(lc['jdmid'], intrapix0, '.')
                plt.ylabel('cflux0')
                plt.subplot(313, sharex=ax, sharey=ax1)
                plt.plot(lc['jdmid'], lc['flux0']/(camtrans0*intrapix0*norm), '.')
                plt.ylabel('cipflux0')
                plt.xlabel('Time [JD]')
                plt.tight_layout()
                plt.show()
                plt.close()
        
        return
    
#ip = IntraPixel()
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPN.hdf5')
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPE.hdf5')
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPS.hdf5')
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPW.hdf5')

cf = CameraFile('/data2/talens/Jul2015/camip_20150716LPC.hdf5')
#cf.visualize()
cf.correct('/data2/talens/Jul2015/fLC_20150715LPC.hdf5')

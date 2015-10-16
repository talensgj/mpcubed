#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

class fLCfile():
    
    def __init__(self, fLCfile):
        
        self.fLCfile = fLCfile
    
        return
        
    def read_header(self, fields, ascc=None):
        
        if ascc is None:
            with h5py.File(self.fLCfile, 'r') as f:
                data = [f['header_table/'+field].value for field in fields]
            
        else:
            nstars = len(ascc)
            nfields = len(fields)
            
            data = [np.zeros(nstars) for _ in range(nfields)]
            
            with h5py.File(self.fLCfile, 'r') as f:
                all_ascc = f['header_table/ascc'].value
                here = np.in1d(all_ascc, ascc)
                for j in range(nfields):
                    data[j] = f['header_table/'+fields[j]].value[here]
                    
        return data
        
    def read_data(self, fields, ascc=None, nobs=None):
        
        if ascc is None:
            ascc, nobs = self.read_header(['ascc', 'nobs'])
        elif nobs is None:
            nobs, = self.read_header(['nobs'], ascc)
            
        nstars = len(ascc)
        nfields = len(fields)
        npoints = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        data = [np.zeros(npoints) for _ in range(nfields)]
        with h5py.File(self.fLCfile, 'r') as f:
            for i in range(nstars):
                lc = f['data/'+ascc[i]]
                for j in range(nfields):
                    data[j][select[i]:select[i+1]] = lc[fields[j]]
            
        return data
            
    def visualize(self, ascc, mode='ha'):
        
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from viridis import viridis 

        rcParams['xtick.labelsize'] = 'large'
        rcParams['ytick.labelsize'] = 'large'
        rcParams['axes.labelsize'] = 'x-large'
        rcParams['image.interpolation'] = 'none'
        rcParams['image.origin'] = 'lower'
        
        if len(ascc) == 1:
            self._single(ascc)
        else:
            self._stacked(ascc, mode)
    
    def _single(self, ascc):
        
        ra, dec, vmag, bmag, nobs = self.read_header(['ra', 'dec', 'vmag', 'bmag', 'nobs'], ascc)
        jdmid, flag, flux0, flux1, peak, sky = self.read_data(['jdmid', 'flag', 'flux0', 'flux1', 'peak', 'sky'], ascc)
        
        here = (flag > 0)
        
        jd_day = np.floor(jdmid)[0]
        jdmid = jdmid - jd_day
        
        plt.figure(figsize=(16,8))
        
        ax = plt.subplot(311)
        plt.title(r'RA = %.2f, Dec = %.2f, V = %.2f, B = %.2f, N = %i'%(ra, dec, vmag, bmag, nobs))
        plt.plot(jdmid, flux0, '.', label='flux0')
        plt.plot(jdmid, flux1, '.', label='flux1')
        plt.ylabel('Flux')
        plt.legend()
        
        plt.subplot(312, sharex=ax)
        plt.plot(jdmid, flux0/flux1, '.', label='flux0/flux1')
        plt.plot(jdmid, (peak - sky)/flux0, '.', label='peak/flux0')
        plt.plot(jdmid, (peak - sky)/flux1, '.', label='peak/flux1')
        plt.ylim(0,1)
        plt.ylabel('Fraction')
        plt.legend()
        
        plt.subplot(313, sharex=ax)
        plt.plot(jdmid, sky, '.')
        plt.xlim(np.amin(jdmid)-.02*np.ptp(jdmid), np.amax(jdmid)+.18*np.ptp(jdmid))
        plt.xlabel('Time [JD-%i]'%jd_day)
        plt.ylabel('Sky')
        
        plt.tight_layout()
        plt.show()
        plt.close()
    
    def _stacked(self, ascc, mode):
        
        nstars = len(ascc)
        
        nobs, ra = self.read_header(['nobs', 'ra'], ascc)
        nobs = nobs.astype('int')
        staridx = np.repeat(np.arange(nstars), nobs)
        
        if (mode == 'lst'):
            lstidx, flux0 = self.read_data(['lstidx', 'flux0'], ascc, nobs)
            lstidx = lstidx.astype('int')
            array = np.full((nstars, 13500), fill_value=np.nan)
            array[staridx, lstidx] = flux0
            
        elif (mode == 'ha'):
            from core.coordinate_grids import PolarGrid
            
            lst, flux0 = self.read_data(['lst', 'flux0'], ascc, nobs)
            ha = np.mod(lst*15.-np.repeat(ra, nobs), 360.)
            pg = PolarGrid(13500, 720)
            haidx = pg.find_raidx(ha)
            array = np.full((nstars, 13502), fill_value=np.nan)
            array[staridx, haidx] = flux0
            array = array[:,1:-1]
            
        else:
            exit()
            
        if (np.ptp(ra) > 270.): ra = np.mod(ra+180, 360)
        array = array[np.argsort(ra)]
        array = array/np.nanmedian(array, axis=1, keepdims=True)
        
        ylim, xlim = np.where(np.isfinite(array))  
        
        plt.figure(figsize=(16,8))
                
        plt.imshow(array, aspect='auto', cmap=viridis, vmin=.5, vmax=1.5)
        plt.colorbar().set_label('Normalized Flux')
        plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        
        if (mode == 'lst'): plt.xlabel('Local Sidereal Time [idx]')
        else: plt.xlabel('Hour Angle [idx]')
        plt.ylabel('Stars, RA order')
        
        plt.tight_layout()
        plt.show()
        plt.close()

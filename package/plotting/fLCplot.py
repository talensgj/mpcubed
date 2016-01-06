#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.IO.fLCfile import fLCfile

class fLCplot():
    
    def __init__(self, fLCfile):
        
        self.fLCfile = fLCfile
        
        return
        
    def coverage(self):
        
        f = fLCfile(self.fLCfile)
        ra, dec, nobs = f.read_header(['ra', 'dec', 'nobs'])
        
        plt.figure(figsize = (16, 8))
        
        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0], xticks = np.linspace(0, 360, 13), yticks = np.linspace(-80, 80, 9))

        im = plt.scatter(ra, dec, c = nobs, cmap = viridis)
        plt.xlim(0, 360)
        plt.ylim(-90, 90)
        plt.xlabel('Hour Angle [degrees]')
        plt.ylabel('Declination [degrees]')
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label('nobs')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        return

    def star_overview(self, ascc):
        
        f = fLCfile(self.fLCfile)
        lc = f.read_star(ascc)
        
        here = (lc['flag'] > 0)
        
        offset = np.floor(lc['jdmid'][0])
        lc['jdmid'] = lc['jdmid'] - offset
        
        baseline = np.ptp(lc['jdmid'])
        xmin = np.amin(lc['jdmid']) - .05*baseline
        xmax = np.amax(lc['jdmid']) + .05*baseline
        
        plt.figure(figsize=(16,8))
        
        ax = plt.subplot(311)
        plt.title(r'ASCC %s'%(ascc))
        plt.plot(lc['jdmid'], lc['flux0'], '.', label='flux0')
        plt.plot(lc['jdmid'], lc['flux1'], '.', label='flux1')
        plt.ylabel('Flux')
        plt.legend()
        
        plt.subplot(312, sharex=ax)
        plt.plot(lc['jdmid'], lc['flux0']/lc['flux1'], '.', label='flux0/flux1')
        plt.plot(lc['jdmid'], (lc['peak'] - lc['sky'])/lc['flux0'], '.', label='peak/flux0')
        plt.plot(lc['jdmid'], (lc['peak'] - lc['sky'])/lc['flux1'], '.', label='peak/flux1')
        plt.ylim(0,1)
        plt.ylabel('Fraction')
        plt.legend()
        
        plt.subplot(313, sharex=ax)
        plt.plot(lc['jdmid'], lc['sky'], '.')
        plt.xlim(xmin, xmax)
        plt.xlabel('Time [JD-%i]'%offset)
        plt.ylabel('Sky')
        
        plt.tight_layout()
        plt.show()
        plt.close()

        return

#class fLCPlot():
    
    #def __init__(self, fLCfile):
        
        #self.fLCfile = fLCfile
    
        #return
        
    #def plot_pointing(self):
        
        #with h5py.File(self.fLCfile, 'r') as f:
            
            #alt0 = f['global/alt0'].value
            #az0 = f['global/az0'].value
            #th0 = f['global/th0'].value
            #x0 = f['global/x0'].value
            #y0 = f['global/y0'].value
        
        #plt.figure(figsize = (16, 8))
        
        #plt.subplot(321)
        #plt.plot((alt0 - alt0[0])*60, '.')
        #plt.ylabel(r'$\Delta$alt0 [arcmin]')
        #plt.ylim(-3, 3)
        
        #plt.subplot(323)
        #plt.plot((az0 - az0[0])*60, '.')
        #plt.ylabel(r'$\Delta$az0 [arcmin]')
        #plt.ylim(-3, 3)
        
        #plt.subplot(325)
        #plt.plot((th0 - th0[0])*60, '.')
        #plt.xlabel('Time')
        #plt.ylabel(r'$\Delta$th0 [arcmin]')
        #plt.ylim(-3, 3)
        
        #plt.subplot(322)
        #plt.plot((x0 - x0[0]), '.')
        #plt.ylabel(r'$\Delta$x0 [pixels]')
        #plt.ylim(-3, 3)
        
        #plt.subplot(324)
        #plt.plot((y0 - y0[0]), '.')
        #plt.xlabel('Time')
        #plt.ylabel(r'$\Delta$y0 [pixels]')
        #plt.ylim(-3, 3)
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
        
        #return
        
    #def plot_coverage(self):
        
        #f = fLCfile(self.fLCfile)
        #ra, dec, nobs = f.read_header(['ra', 'dec', 'nobs'])
        
        #plt.figure(figsize = (16, 8))
        
        #gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        #plt.subplot(gs[1,0], xticks = np.linspace(0, 360, 13), yticks = np.linspace(-80, 80, 9))

        #im = plt.scatter(ra, dec, c = nobs, cmap = viridis)
        #plt.xlim(0, 360)
        #plt.ylim(-90, 90)
        #plt.xlabel('Hour Angle [degrees]')
        #plt.ylabel('Declination [degrees]')
        
        #cax = plt.subplot(gs[1,1])
        #cb = plt.colorbar(im, cax = cax)
        #cb.set_label('nobs')
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
        
        #return
        
    #def plot_star(self, ascc, field, mode = 'array'):
        
        #with h5py.File(self.fLCfile, 'r') as f:
            
            #try:
                #lc = f['data/' + ascc].value
            #except:
                #print 'The star was not found in the file, exiting...'
                #exit()
                
        #idx1 = lc['lstseq'] // 13500
        #idx2 = lc['lstseq'] % 13500
                
        #idx1 = idx1 - np.amin(idx1)
                
        #nx = np.amax(idx1) + 1
        #ny = 13500
                
        #image = np.full((nx, ny), fill_value = np.nan)
        #image[idx1, idx2] = lc[field]
        
        #plt.figure(figsize = (16, 8))
        
        #plt.suptitle('ASCC %s'%ascc, size='xx-large')
        
        #gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        #plt.subplot(gs[1,0])
        
        #im = plt.imshow(image, aspect = 'auto', cmap = viridis)
        #plt.xlim(np.amin(idx2) - .5, np.amax(idx2) + .5)
        #plt.xlabel('lstidx')
        #plt.ylabel('lstday')
        
        #cax = plt.subplot(gs[1,1])
        #cb = plt.colorbar(im, cax = cax)
        #cb.set_label(field)
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
        
        #return


#def visualize(self, ascc, mode='ha'):
        
        #if len(ascc) == 1:
            #self._single(ascc)
        #else:
            #self._stacked(ascc, mode)
    
    #def _single(self, ascc):
        
        #ra, dec, vmag, nobs = self.read_header(['ra', 'dec', 'vmag', 'nobs'], ascc)
        #jdmid, flag, flux0, flux1, peak, sky = self.read_data(['jdmid', 'flag', 'flux0', 'flux1', 'peak', 'sky'], ascc)
        
        #here = (flag > 0)
        
        #jd_day = np.floor(jdmid)[0]
        #jdmid = jdmid - jd_day
        
        #plt.figure(figsize=(16,8))
        
        #ax = plt.subplot(311)
        #plt.title(r'ASCC %s, RA = %.2f, Dec = %.2f, V = %.2f, N = %i'%(ascc[0], ra, dec, vmag, nobs))
        #plt.plot(jdmid, flux0, '.', label='flux0')
        #plt.plot(jdmid, flux1, '.', label='flux1')
        #plt.ylabel('Flux')
        #plt.legend()
        
        #plt.subplot(312, sharex=ax)
        #plt.plot(jdmid, flux0/flux1, '.', label='flux0/flux1')
        #plt.plot(jdmid, (peak - sky)/flux0, '.', label='peak/flux0')
        #plt.plot(jdmid, (peak - sky)/flux1, '.', label='peak/flux1')
        #plt.ylim(0,1)
        #plt.ylabel('Fraction')
        #plt.legend()
        
        #plt.subplot(313, sharex=ax)
        #plt.plot(jdmid, sky, '.')
        #plt.xlim(np.amin(jdmid)-.02*np.ptp(jdmid), np.amax(jdmid)+.18*np.ptp(jdmid))
        #plt.xlabel('Time [JD-%i]'%jd_day)
        #plt.ylabel('Sky')
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
    
    #def _stacked(self, ascc, mode):
        
        #nstars = len(ascc)
        
        #nobs, ra = self.read_header(['nobs', 'ra'], ascc)
        #nobs = nobs.astype('int')
        #staridx = np.repeat(np.arange(nstars), nobs)
        
        #if (mode == 'lst'):
            #lstidx, flux0 = self.read_data(['lstidx', 'flux0'], ascc, nobs)
            #lstidx = lstidx.astype('int')
            #array = np.full((nstars, 13500), fill_value=np.nan)
            #array[staridx, lstidx] = flux0
            
        #elif (mode == 'ha'):
            #from core.coordinate_grids import PolarGrid
            
            #lst, flux0 = self.read_data(['lst', 'flux0'], ascc, nobs)
            #ha = np.mod(lst*15.-np.repeat(ra, nobs), 360.)
            #pg = PolarGrid(13500, 720)
            #haidx = pg.find_raidx(ha)
            #array = np.full((nstars, 13502), fill_value=np.nan)
            #array[staridx, haidx] = flux0
            #array = array[:,1:-1]
            
        #else:
            #exit()
            
        #if (np.ptp(ra) > 270.): ra = np.mod(ra+180, 360)
        #array = array[np.argsort(ra)]
        #array = array/np.nanmedian(array, axis=1, keepdims=True)
        
        #ylim, xlim = np.where(np.isfinite(array))  
        
        #plt.figure(figsize=(16,8))
                
        #plt.imshow(array, aspect='auto', cmap=viridis, vmin=.5, vmax=1.5)
        #plt.colorbar().set_label('Normalized Flux')
        #plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        
        #if (mode == 'lst'): plt.xlabel('Local Sidereal Time [idx]')
        #else: plt.xlabel('Hour Angle [idx]')
        #plt.ylabel('Stars, RA order')
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()



if __name__ == '__main__':
    
    obj = fLCPlot('/data2/talens/2015Q2/LPC/fLC_201506ALPC.hdf5')
    obj.star_overview('807144')

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

from fLCfile import fLCfile

class fLCPlot():
    
    def __init__(self, fLCfile):
        
        self.fLCfile = fLCfile
    
        return
        
    def plot_pointing(self):
        
        with h5py.File(self.fLCfile, 'r') as f:
            
            alt0 = f['global/alt0'].value
            az0 = f['global/az0'].value
            th0 = f['global/th0'].value
            x0 = f['global/x0'].value
            y0 = f['global/y0'].value
        
        plt.figure(figsize = (16, 8))
        
        plt.subplot(321)
        plt.plot((alt0 - alt0[0])*60, '.')
        plt.ylabel(r'$\Delta$alt0 [arcmin]')
        plt.ylim(-3, 3)
        
        plt.subplot(323)
        plt.plot((az0 - az0[0])*60, '.')
        plt.ylabel(r'$\Delta$az0 [arcmin]')
        plt.ylim(-3, 3)
        
        plt.subplot(325)
        plt.plot((th0 - th0[0])*60, '.')
        plt.xlabel('Time')
        plt.ylabel(r'$\Delta$th0 [arcmin]')
        plt.ylim(-3, 3)
        
        plt.subplot(322)
        plt.plot((x0 - x0[0]), '.')
        plt.ylabel(r'$\Delta$x0 [pixels]')
        plt.ylim(-3, 3)
        
        plt.subplot(324)
        plt.plot((y0 - y0[0]), '.')
        plt.xlabel('Time')
        plt.ylabel(r'$\Delta$y0 [pixels]')
        plt.ylim(-3, 3)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        return
        
    def plot_coverage(self):
        
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
        
    def plot_star(self, ascc, field, mode = 'array'):
        
        with h5py.File(self.fLCfile, 'r') as f:
            
            try:
                lc = f['data/' + ascc].value
            except:
                print 'The star was not found in the file, exiting...'
                exit()
                
        idx1 = lc['lstseq'] // 13500
        idx2 = lc['lstseq'] % 13500
                
        idx1 = idx1 - np.amin(idx1)
                
        nx = np.amax(idx1) + 1
        ny = 13500
                
        image = np.full((nx, ny), fill_value = np.nan)
        image[idx1, idx2] = lc[field]
        
        plt.figure(figsize = (16, 8))
        
        plt.suptitle('ASCC %s'%ascc, size='xx-large')
        
        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0])
        
        im = plt.imshow(image, aspect = 'auto', cmap = viridis)
        plt.xlim(np.amin(idx2) - .5, np.amax(idx2) + .5)
        plt.xlabel('lstidx')
        plt.ylabel('lstday')
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label(field)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        return


if __name__ == '__main__':
    
    obj = fLCPlot('/data2/talens/2015Q2/LPC/fLC_201506ALPC.hdf5')
    obj.plot_pointing()
    obj.plot_coverage()
    obj.plot_star('807144', 'flux0')

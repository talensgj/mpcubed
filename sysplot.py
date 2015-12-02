#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import mascara
from sysfile import SysFile

class SysPlot():
    
    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        return
    
    def hadec2xy(self, ha, dec):
        
        tmp = ha.shape
        ha, dec = ha.ravel(), dec.ravel()
        
        site = mascara.observer.Site('LaPalma')
        cam = mascara.observer.Camera('east')
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        return x, y
        
    def add_hadecgrid(self):
        
        site = mascara.observer.Site('LaPalma')
        cam = mascara.observer.Camera('east')
        
        ha = np.linspace(0, 360, 360)
        dec = np.linspace(-80, 80, 17)
        ha, dec = np.meshgrid(ha, dec)
        
        tmp = ha.shape
        
        ha, dec = ha.ravel(), dec.ravel()
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        here = (x > -50) & (x < 4008+50) & (y > -50) & (y < 2672+50)
        x[~here] = np.nan
        y[~here] = np.nan
        
        plt.plot(x.T, y.T, c='k')
        
        ha = np.linspace(0, 345, 24)
        dec = np.linspace(-80, 80, 160)
        ha, dec = np.meshgrid(ha, dec)
        
        tmp = ha.shape
        
        ha, dec = ha.ravel(), dec.ravel()
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        here = (x > -50) & (x < 4008+50) & (y > -50) & (y < 2672+50)
        x[~here] = np.nan
        y[~here] = np.nan
        
        plt.plot(x, y, c='k')
        
        return
    
    def plot_statistics(self, mode):
        
        f = SysFile(self.sysfile)
        niter, chisq, npoints, npars = f.read_statistics(mode)
        
        print 'chisq = %.3f'%(np.sum(chisq)/(np.sum(npoints) - np.sum(npars)))
        print 'mean(niter) = %i'%np.mean(niter)
        print 'mean(chisq) = %.3f'%np.mean(chisq/(npoints - npars))
        
        plt.hist(niter, bins=np.linspace(.5, 99.5, 100))
        plt.show()
        
        plt.hist(chisq/(npoints - npars), bins=np.linspace(0, 50, 51))
        plt.show()
        
        return
    
    def plot_camtrans(self):
        
        f = SysFile(self.sysfile)
        z, nobs = f.read_camtrans()
        
        z = z[1:-1,1:-1]
        z = np.ma.masked_invalid(z)
        
        nobs = nobs[1:-1,1:-1]
        nobs = np.ma.masked_invalid(nobs)
        
        ha = np.linspace(0, 360, 13501)
        dec = np.linspace(-90, 90, 721)
        
        ha, dec = np.meshgrid(ha, dec)
        
        x, y = self.hadec2xy(ha, dec)
        
        fig = plt.figure(figsize=(13,9))

        plt.suptitle('Transmission', size='xx-large')

        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0], aspect='equal')
        im = plt.pcolormesh(x, y, z.T, cmap=viridis)
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        plt.subplot(111, aspect='equal')
        plt.pcolormesh(x, y, nobs.T)
        plt.colorbar()
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        plt.show()
        plt.close()
        
    def plot_intrapix(self):
        
        f = SysFile(self.sysfile)
        #niter, station, camera = f.read_header()
        a, b, c, d, nobs = f.read_intrapix()
        
        a = a[1:-1,1:-1]
        a = np.ma.masked_invalid(a)
        
        b = b[1:-1,1:-1]
        b = np.ma.masked_invalid(b)
        
        c = c[1:-1,1:-1]
        c = np.ma.masked_invalid(c)
        
        d = d[1:-1,1:-1]
        d = np.ma.masked_invalid(d)
        
        nobs = nobs[1:-1,1:-1]
        nobs = np.ma.masked_invalid(nobs)
        
        ha = np.linspace(0, 360, 271)
        dec = np.linspace(-90, 90, 721)
        
        ha, dec = np.meshgrid(ha, dec)
        
        x, y = self.hadec2xy(ha, dec)
        
        fig = plt.figure(figsize=(16, 10))

        plt.suptitle('Intrapixel variations.', size='xx-large')

        gs = gridspec.GridSpec(3, 3, width_ratios = [15,15,.5], height_ratios = [1,10,10])
        
        plt.subplot(gs[1,0], aspect='equal')
        plt.title(r'$\sin(2\pi x)$')
        plt.pcolormesh(x, y, a.T, vmin=-.1, vmax=.1, cmap=viridis)
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[1,1], aspect='equal')
        plt.title(r'$\cos(2\pi x)$')
        plt.pcolormesh(x, y, b.T, vmin=-.1, vmax=.1, cmap=viridis)
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[2,0], aspect='equal')
        plt.title(r'$\sin(2\pi y)$')
        plt.pcolormesh(x, y, c.T, vmin=-.1, vmax=.1, cmap=viridis)
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[2,1], aspect='equal')
        plt.title(r'$\cos(2\pi y)$')
        im = plt.pcolormesh(x, y, d.T, vmin=-.1, vmax=.1, cmap=viridis)
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        cax = plt.subplot(gs[1:,2])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label('Amplitude')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        plt.subplot(111, aspect='equal')
        plt.pcolormesh(x, y, nobs.T)
        plt.colorbar()
        self.add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        plt.show()
        plt.close()
        
if __name__ == '__main__':
    
    obj = SysPlot('/data2/talens/2015Q2/LPE/sys_201506BLPE.hdf5')
    obj.plot_camtrans()
        
    
    
    

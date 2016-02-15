#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

import mascara
from .. import IO

class SysPlot(object):
    """ Plot the data from a systematics file.
    
    Attributes:
        sysfile (str): The systematics file to read the data from.
    
    """
    
    def __init__(self, sysfile):
        """ Initialize the plotter of systematics files.
        
        Args:
            sysfile (str): The systematics file to read the data from.
    
        """
        
        self.sysfile = sysfile
        
        return
    
    def _hadec2xy(self, ha, dec):
        
        # Read the pointing from the file.
        f = IO.SysFile(self.sysfile)
        alt0, az0, th0, x0, y0 = f.read_pointing()
        
        # Initialize the coordinate transformations.
        site = mascara.observer.Site('LaPalma')
        #cam = mascara.observer.Camera('east')
        cam = mascara.observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
        
        tmp = ha.shape
        ha, dec = ha.ravel(), dec.ravel()
        
        # Perfrom the coordinate transformations.
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        return x, y
        
    def _xy2hadec(self, x, y):
        
        # Read the pointing from the file.
        f = IO.SysFile(self.sysfile)
        alt0, az0, th0, x0, y0 = f.read_pointing()
        
        # Initialize the coordinate transformations.
        site = mascara.observer.Site('LaPalma')
        #cam = mascara.observer.Camera('east')
        cam = mascara.observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
        
        tmp = x.shape
        x, y = x.ravel(), y.ravel()
        
        # Perfrom the coordinate transformations.
        phi, the = cam.XY2PhiThe(x, y)
        alt, az = cam.PhiThe2Hor(phi, the)
        ha, dec = site.altaz2hadec(alt, az)
        
        ha, dec = ha.reshape(tmp), dec.reshape(tmp)
        
        return ha, dec
        
    def _add_hadecgrid(self):
        
        # Read the pointing from the file.
        f = IO.SysFile(self.sysfile)
        alt0, az0, th0, x0, y0 = f.read_pointing()
        
        # Initialize the coordinate transformations.
        site = mascara.observer.Site('LaPalma')
        #cam = mascara.observer.Camera('east')
        cam = mascara.observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
        
        # Add lines of constant declination.
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
        
        # Add lines of constant hour angle.
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
    
    def plot_magnitudes(self, display=True, savefig=False):
        """ Plot the fitted magnitudes against the catalogue V magnitude. """
        
        # Read the data.
        f = IO.SysFile(self.sysfile)
        ascc, vmag, mag, sigma, nobs = f.read_magnitudes()
        
        # Create the magnitudes plot.
        fig = plt.figure(figsize=(13,9))

        plt.suptitle('Magnitudes', size='xx-large')

        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0])
        
        plt.plot(vmag, mag, '.', alpha = .2)
        plt.xlim(2, 8.4)
        plt.xlabel('V [mag]')
        plt.ylabel('M [mag]')
        
        plt.tight_layout()
        if savefig:
            head, tail = os.path.split(self.sysfile)
            suffix = '_mag.png'
            tail = tail.rsplit('.')[0] + suffix
            plt.savefig(os.path.join(head, tail))
        if display:
            plt.show()
        plt.close()
        
        return
    
    def plot_trans(self, display=True, savefig=False):
        """ Plot the fitted transmission map as a function of x and y. """
        
        # Read the data.
        f = IO.SysFile(self.sysfile)
        pg, trans, nobs = f.read_trans()
        
        # Remove the bounadry and mask NaNs.
        trans = trans[1:-1,1:-1]
        trans = np.ma.masked_invalid(trans)
        
        nobs = nobs[1:-1,1:-1]
        nobs = np.ma.masked_invalid(nobs)
        
        # Create the coordinate grid.
        ha, dec = pg.xedges, pg.yedges
        ha, dec = np.meshgrid(ha, dec)
        x, y = self._hadec2xy(ha, dec)
        
        # Create the transmission plot.
        fig = plt.figure(figsize=(14,9))

        plt.suptitle('Transmission', size='xx-large')

        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0], aspect='equal')
        
        vmin = np.nanpercentile(trans, 1)
        vmax = np.nanpercentile(trans, 99)
        delta = vmax - vmin
        vmin = vmin - .01*delta
        vmax = vmax + .01*delta
        
        im = plt.pcolormesh(x, y, trans.T, cmap=viridis, vmin=vmin, vmax=vmax)
        self._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.ax.invert_yaxis()
        cb.set_label('Magnitude')
        
        plt.tight_layout()
        if savefig:
            head, tail = os.path.split(self.sysfile)
            suffix = '_trans.png'
            tail = tail.rsplit('.')[0] + suffix
            plt.savefig(os.path.join(head, tail))
        if display:
            plt.show()
        plt.close()
        
        ## Create the nobs plot.
        #fig = plt.figure(figsize=(13,9))

        #plt.suptitle('Transmission', size='xx-large')

        #gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        #plt.subplot(gs[1,0], aspect='equal')
        
        #im = plt.pcolormesh(x, y, nobs.T, cmap=viridis)
        #self._add_hadecgrid()
        #plt.xlim(0, 4008)
        #plt.ylim(0, 2672)
        
        #cax = plt.subplot(gs[1,1])
        #cb = plt.colorbar(im, cax = cax)
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
        
        return
        
    def plot_intrapix(self, display=True, savefig=False):
        """ Plot the fitted amplitude maps as a function of x and y. """
        
        # Read the data.
        f = IO.SysFile(self.sysfile)
        pg, a, b, c, d, nobs = f.read_intrapix()
        
        # Remove the boundary and mask NaNs.
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
        
        # Create the coordinate grid.
        ha, dec = pg.xedges, pg.yedges
        ha, dec = np.meshgrid(ha, dec)
        x, y = self._hadec2xy(ha, dec)
        
        # Create the intrapixel plot.
        fig = plt.figure(figsize=(16, 10))

        plt.suptitle('Intrapixel variations.', size='xx-large')

        gs = gridspec.GridSpec(3, 3, width_ratios = [15,15,.5], height_ratios = [1,10,10])
        
        plt.subplot(gs[1,0], aspect='equal')
        plt.title(r'$\sin(2\pi x)$')
        plt.pcolormesh(x, y, a.T, vmin=-.1, vmax=.1, cmap=viridis)
        self._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[1,1], aspect='equal')
        plt.title(r'$\cos(2\pi x)$')
        plt.pcolormesh(x, y, b.T, vmin=-.1, vmax=.1, cmap=viridis)
        self._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[2,0], aspect='equal')
        plt.title(r'$\sin(2\pi y)$')
        plt.pcolormesh(x, y, c.T, vmin=-.1, vmax=.1, cmap=viridis)
        self._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        plt.subplot(gs[2,1], aspect='equal')
        plt.title(r'$\cos(2\pi y)$')
        im = plt.pcolormesh(x, y, d.T, vmin=-.1, vmax=.1, cmap=viridis)
        self._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        cax = plt.subplot(gs[1:,2])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label('Amplitude')
        
        plt.tight_layout()
        if savefig:
            head, tail = os.path.split(self.sysfile)
            suffix = '_ipx.png'
            tail = tail.rsplit('.')[0] + suffix
            plt.savefig(os.path.join(head, tail))
        if display:
            plt.show()
        plt.close()
        
        ## Create the nobs plot.
        #fig = plt.figure(figsize=(13,9))

        #plt.suptitle('Intrapixel variations.', size='xx-large')

        #gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        #plt.subplot(gs[1,0], aspect='equal')
        
        #im = plt.pcolormesh(x, y, nobs.T, cmap=viridis)
        #self._add_hadecgrid()
        #plt.xlim(0, 4008)
        #plt.ylim(0, 2672)
        
        #cax = plt.subplot(gs[1,1])
        #cb = plt.colorbar(im, cax = cax)
        
        #plt.tight_layout()
        #plt.show()
        #plt.close()
        
        return
        
    def plot_clouds(self):
        """ Plot the fitted cloud map as a function of x, y and t."""
        
        import healpy
        
        # Read the data.
        f = IO.SysFile(self.sysfile)
        hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()
        
        #for i in range(0, clouds.shape[1], 50):
            
            #if np.all(np.isnan(clouds[:,i])): continue
            
            #healpy.mollview(clouds[:,i], title = 'Clouds %i'%(i + lstmin), cmap = viridis, min=-.5, max=.5)
            #healpy.graticule(dpar = 10., dmer = 15.)
            
            ##plt.tight_layout()
            #plt.show()
            #plt.close()
            
        xedges = np.linspace(0, 4008, 3*167)
        yedges = np.linspace(0, 2672, 2*167)
            
        xcenter = (xedges[:-1] + xedges[1:])/2.
        ycenter = (yedges[:-1] + yedges[1:])/2.
        
        xcenter, ycenter = np.meshgrid(xcenter, ycenter)
        ha, dec = self._xy2hadec(xcenter, ycenter)
        
        for i in range(12036700-lstmin, 12039236-lstmin+1):
            
            if np.all(np.isnan(clouds[:,i])): continue
            print i + lstmin
            
            lst = ((i + lstmin) % 13500)*24./13500
            ra = np.mod(lst*15. - ha, 360.)
            idx = hg.radec2idx(ra, dec)
            
            array = clouds[idx, i]
            array = np.ma.masked_invalid(array)
            
            fig = plt.figure(figsize=(14,9))

            plt.suptitle('Clouds 2015-06-08 East', size='xx-large')

            gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
            
            plt.subplot(gs[1,0], aspect='equal')
            plt.annotate('seq = {}'.format(i + lstmin), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='large', backgroundcolor='w')
            
            im = plt.pcolormesh(xedges, yedges, array, cmap=cm.Greys, vmin=-.25, vmax=.25)
            self._add_hadecgrid()
            plt.xlim(0, 4008)
            plt.ylim(0, 2672)
            
            cax = plt.subplot(gs[1,1])
            cb = plt.colorbar(im, cax = cax)
            cb.ax.invert_yaxis()
            cb.set_label('Magnitude')
            
            plt.tight_layout()
            plt.savefig('/data2/talens/clouds/img_{}.png'.format(i+lstmin), dpi=160)
            #plt.show()
            plt.close()
            
        return
           
    #def plot_clouds2(self):
    
        #import healpy
        
        ## Read the data.
        #f = IO.SysFile(self.sysfile)
        #hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()
        
        #for i in range(lstmin//50, lstmax//50):
            
            #i1 = np.minimum((i - 1)*50 - lstmin, lstmin)
            #i2 = np.minimum(i*50 - lstmin, lstmin)
            
            #i3 = np.minimum(i*50 - lstmin, lstmin)
            #i4 = np.minimum((i + 1)*50 - lstmin, lstmax)
            
            #i5 = np.minimum((i + 1)*50 - lstmin, lstmax)
            #i6 = np.minimum((i + 2)*50 - lstmin, lstmax)
            
            #fig = plt.figure(figsize=(14, 9))
            
            #plt.suptitle('Clouds')
            
            #gs = gridspec.GridSpec(4, 2, width_ratios=[15,.5], height_ratios=[1,10,10,10])
            
            #plt.subplot(gs[1,0])
            #im = plt.imshow(clouds[:,i1:i2].T, aspect='auto', vmin=-.5, vmax=.5)
            
            #plt.subplot(gs[2,0])
            #im = plt.imshow(clouds[:,i3:i4].T, aspect='auto', vmin=-.5, vmax=.5)
            
            #plt.subplot(gs[3,0])
            #im = plt.imshow(clouds[:,i5:i6].T, aspect='auto', vmin=-.5, vmax=.5)
            
            #cax = plt.subplot(gs[1:,1])
            #cb = plt.colorbar(im, cax = cax)
            #cb.ax.invert_yaxis()
            #cb.set_label('Magnitude')
            
            #plt.tight_layout()
            #plt.show()
            #plt.close()
            
        #return

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os
import glob

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from . import IO

def _hadec2xy(wcspars, ha, dec):
        
    import mascara        
        
    alt0, az0, th0, x0, y0 = wcspars
    
    # Initialize the coordinate transformations.
    site = mascara.observer.Site(28.76025,  -17.8792,  2364.)
    cam = mascara.observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
    
    tmp = ha.shape
    ha, dec = ha.ravel(), dec.ravel()
    
    # Perfrom the coordinate transformations.
    alt, az = site.hadec2altaz(ha, dec, degree=True)
    phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
    x, y = cam.PhiThe2XY(phi, theta)
    
    x, y = x.reshape(tmp), y.reshape(tmp)
    
    return x, y

def wcsgrid(wcspars):

    # Add lines of constant declination.
    ha = np.linspace(0, 360, 360)
    dec = np.linspace(-80, 80, 17)
    ha, dec = np.meshgrid(ha, dec)
    
    tmp = ha.shape
    
    ha, dec = ha.ravel(), dec.ravel()
    
    x, y = _hadec2xy(wcspars, ha, dec)
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
    
    x, y = _hadec2xy(wcspars, ha, dec)
    x, y = x.reshape(tmp), y.reshape(tmp)
    
    here = (x > -50) & (x < 4008+50) & (y > -50) & (y < 2672+50)
    x[~here] = np.nan
    y[~here] = np.nan
    
    plt.plot(x, y, c='k')
    
    return

def plot_polar(grid, data, wcspars, **kwargs):
    
    data = data[1:-1,1:-1]    
    data = np.ma.array(data, mask=np.isnan(data))    
    
    ha, dec = grid.xedges, grid.yedges
    ha, dec = np.meshgrid(ha, dec)
    x, y = _hadec2xy(wcspars, ha, dec)
    
    im = plt.pcolormesh(x, y, data.T, **kwargs)
    wcsgrid(wcspars)

    return im

def fig_transmission(filename, figname):
    
    # Read the transmission map.
    f = IO.SysFile(filename)
    pg, trans, nobs = f.read_trans()
    wcspars = f.read_pointing()    
    
    # Plot the transmission map.
    fig = plt.figure(figsize=(14,9))
    
    plt.suptitle('Transmission', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
    
    plt.subplot(gs[1,0], aspect='equal')
    
    vmin = np.nanpercentile(trans, .1)
    vmax = np.nanpercentile(trans, 99.9)
    im = plot_polar(pg, trans, wcspars, cmap=plt.cm.viridis, vmin=vmin, vmax=vmax)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1,1])
    cb = plt.colorbar(im, cax = cax)
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()
    
    return
    
def fig_intrapix(filename, figname):
    
    # Read the intrapixel amplitudes.
    f = IO.SysFile(filename) 
    pg, sinx, cosx, siny, cosy, nobs = f.read_intrapix() 
    wcspars = f.read_pointing()       
    
    # Plot the amplitudes.
    fig = plt.figure(figsize=(16, 10))
    
    plt.suptitle('Intrapixel Amplitudes', size='xx-large')
    
    gs = gridspec.GridSpec(3, 3, width_ratios = [15,15,.5], height_ratios = [1,10,10])
    
    plt.subplot(gs[1,0], aspect='equal')
    plt.title(r'$\sin(2\pi x)$')
   
    im = plot_polar(pg, sinx, wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)
   
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[1,1], aspect='equal')
    plt.title(r'$\cos(2\pi x)$')
    
    im = plot_polar(pg, cosx, wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)   
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,0], aspect='equal')
    plt.title(r'$\sin(2\pi y)$')
    
    im = plot_polar(pg, siny, wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)  
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,1], aspect='equal')
    plt.title(r'$\cos(2\pi y)$')
    
    im = plot_polar(pg, cosy, wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)   
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1:,2])
    cb = plt.colorbar(im, cax = cax)
    cb.set_label('Amplitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()  
    
    return
    
def fig_clouds(filename, figname):
    
    # Read the data.
    f = IO.SysFile(filename)
    hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Clouds', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,20])
    
    plt.subplot(gs[1,0], xticks=[], yticks=[])
    
    mask = np.isfinite(clouds)
    idx1, = np.where(np.any(mask, axis=1))
    idx2, = np.where(np.any(mask, axis=0))
    clouds = clouds[idx1][:,idx2]    
    im = plt.imshow(clouds.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=-.5, vmax=.5)
    
    idx, = np.where(np.diff(idx2) > 1)
    for i in idx:
        plt.axhline(i+1, xmax=.1, c='k', lw=2)
    
    plt.ylabel('Time')
    plt.xlabel('Sky Patch')
    
    cax = plt.subplot(gs[1,1])
    cb = plt.colorbar(im, cax = cax)
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()    
    
    return
    
def fig_sigma(filename, figname):
    
    # Read the data.
    f = IO.SysFile(filename)
    hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Sigma', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,20])
    
    plt.subplot(gs[1,0], xticks=[], yticks=[])
    
    mask = np.isfinite(sigma)
    idx1, = np.where(np.any(mask, axis=1))
    idx2, = np.where(np.any(mask, axis=0))
    sigma = sigma[idx1][:,idx2]    
    im = plt.imshow(sigma.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=0, vmax=np.nanmax(sigma))
    
    idx, = np.where(np.diff(idx2) > 1)
    for i in idx:
        plt.axhline(i+1, xmax=.1, c='k', lw=2)
    
    plt.ylabel('Time')
    plt.xlabel('Sky Patch')
    
    cax = plt.subplot(gs[1,1])
    cb = plt.colorbar(im, cax = cax)
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()    
    
    return
    
def calibration_summary(filename):
    
    # Split the filepath.
    head, tail = os.path.split(filename)
    
    # Make figure showing the transmission.
    figname = tail.rsplit('.')[0] + '_trans.png'
    figname = os.path.join(head, figname)   
    
    fig_transmission(filename, figname)
    
    # Make a figure showing the intrapixel variations.
    figname = tail.rsplit('.')[0] + '_ipx.png'
    figname = os.path.join(head, figname)  
    
    fig_intrapix(filename, figname)
    
    # Make a figure showing the clouds.
    figname = tail.rsplit('.')[0] + '_clouds.png'
    figname = os.path.join(head, figname)  
    
    fig_clouds(filename, figname)
    
    # Make a figure showing the additional uncertainties.
    figname = tail.rsplit('.')[0] + '_sigma.png'
    figname = os.path.join(head, figname)      
    
    fig_sigma(filename, figname)
    
    return

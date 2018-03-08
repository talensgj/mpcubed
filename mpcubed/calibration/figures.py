# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from .. import io
from lsreduce import io as lsio

def _hadec2xy(wcspars, ha, dec):
        
    import mascara 
    from lsreduce import astrometry       
            
    # Perfrom the coordinate transformations.
    try:
        wcspars['lst']
    except:
        tmp = ha.shape
        ha, dec = ha.ravel(), dec.ravel()
        
        alt0, az0, th0, x0, y0 = wcspars
        site = mascara.observer.Site(28.76025,  -17.8792,  2364.)
        cam = mascara.observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        x, y = x.reshape(tmp), y.reshape(tmp)
    else:
        ra = astrometry.ha2ra(ha, 0.)    
        x, y = astrometry.world2wcs(wcspars, ra, dec, 0.)
        
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

def fig_magnitudes(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_mags.png'
        figname = os.path.join(head, figname)
        
    # Read the magnitudes.
    f = io.SysFile(filename)
    magnitudes = f.read_magnitudes()
    mag, vmag, sigma = magnitudes['mag'], magnitudes['vmag'], magnitudes['sigma']
    
    # Plot the magitudes.
    fig = plt.figure(figsize=(14,9))
    
    plt.suptitle('Magnitudes', size='xx-large')
    
    gs = gridspec.GridSpec(3, 1, height_ratios = [1,5,5])
    
    plt.subplot(gs[1,0])
    
    plt.scatter(vmag, mag - vmag, color='black', marker='.', alpha=.5, edgecolor='none')
    plt.ylabel('V - M')
    
    plt.subplot(gs[2,0])
    
    plt.scatter(vmag, sigma, color='black', marker='.', alpha=.5, edgecolor='none')
    plt.xlabel('V')
    plt.ylabel(r'$\sigma$')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()
    
    return figname

def fig_transmission(filename, astromaster=None, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_trans.png'
        figname = os.path.join(head, figname) 
    
    # Read the transmission map.
    f = io.SysFile(filename)
    
    if astromaster is None:
        wcspars = f.read_pointing()
    else:
        wcspars, polpars, astromask = lsio.read_astromaster(astromaster) 
        
    camgrid, trans = f.read_trans()
        
    trans = trans['trans']
    
    # Plot the transmission map.
    fig = plt.figure(figsize=(14,9))
    
    plt.suptitle('Transmission', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
    
    plt.subplot(gs[1,0], aspect='equal')
    
    vmin = np.nanpercentile(trans, .1)
    vmax = np.nanpercentile(trans, 99.9)
    im = plot_polar(camgrid, trans, wcspars, cmap=plt.cm.viridis, vmin=vmin, vmax=vmax)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1,1])
    cb = plt.colorbar(im, cax = cax)
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()
    
    return figname
    
def fig_intrapix(filename, astromaster=None, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_ipx.png'
        figname = os.path.join(head, figname)   
    
    # Read the intrapixel amplitudes.
    f = io.SysFile(filename) 
    
    if astromaster is None:
        wcspars = f.read_pointing()
    else:
        wcspars, polpars, astromask = lsio.read_astromaster(astromaster)
        
    ipxgrid, intrapix = f.read_intrapix() 
    
    amplitudes = intrapix['amplitudes']       
    
    # Plot the amplitudes.
    fig = plt.figure(figsize=(16, 10))
    
    plt.suptitle('Intrapixel Amplitudes', size='xx-large')
    
    gs = gridspec.GridSpec(3, 3, width_ratios = [15,15,.5], height_ratios = [1,10,10])
    
    plt.subplot(gs[1,0], aspect='equal')
    plt.title(r'$\sin(2\pi x)$')
   
    im = plot_polar(ipxgrid, amplitudes[:,:,0], wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)
   
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[1,1], aspect='equal')
    plt.title(r'$\cos(2\pi x)$')
    
    im = plot_polar(ipxgrid, amplitudes[:,:,1], wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)   
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,0], aspect='equal')
    plt.title(r'$\sin(2\pi y)$')
    
    im = plot_polar(ipxgrid, amplitudes[:,:,2], wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)  
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,1], aspect='equal')
    plt.title(r'$\cos(2\pi y)$')
    
    im = plot_polar(ipxgrid, amplitudes[:,:,3], wcspars, cmap=plt.cm.viridis, vmin=-.05, vmax=.05)   
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1:,2])
    cb = plt.colorbar(im, cax = cax)
    cb.set_label('Amplitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()  
    
    return figname
    
def fig_clouds(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_clouds.png'
        figname = os.path.join(head, figname)
    
    # Read the data.
    f = io.SysFile(filename)
    skygrid, clouds = f.read_clouds()
    
    clouds = clouds['clouds']
    
    # Remove empty rows and columns.
    mask = np.isfinite(clouds)
    q, = np.where(np.any(mask, axis=1))
    t, = np.where(np.any(mask, axis=0))
    clouds = clouds[q][:,t]
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Clouds', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,20])
    
    plt.subplot(gs[1,0], xticks=[], yticks=[])
        
    im = plt.imshow(clouds.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=-.5, vmax=.5)
    
    idx, = np.where(np.diff(t) > 2)
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
    
    return figname
    
def fig_sigma(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_sigma.png'
        figname = os.path.join(head, figname)
    
    # Read the data.
    f = io.SysFile(filename)
    skygrid, clouds = f.read_clouds()
    
    sigma = clouds['sigma']
    
    # Remove empty rows and columns.
    mask = np.isfinite(sigma)
    idx1, = np.where(np.any(mask, axis=1))
    idx2, = np.where(np.any(mask, axis=0))
    sigma = sigma[idx1][:,idx2]  
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Sigma', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,20])
    
    plt.subplot(gs[1,0], xticks=[], yticks=[])
      
    im = plt.imshow(sigma.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=0, vmax=np.nanmax(sigma))
    
    idx, = np.where(np.diff(idx2) > 2)
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
    
    return figname
    
def figs_calibration(filelist, astromaster=None):
    
    for filename in filelist:
        
        fig_magnitudes(filename)
        fig_transmission(filename, astromaster)
        fig_intrapix(filename, astromaster)
        fig_clouds(filename)
        fig_sigma(filename)
    
    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Make figures of calibration terms.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to create the figures for')
    parser.add_argument('-a', '--astro', type=str, default=None,
                        help='the astrometric solution to use when creating the figures')
    args = parser.parse_args()
    
    figs_calibration(args.files, args.astro)
    
    return
    
if __name__ == '__main__':
    
    main()

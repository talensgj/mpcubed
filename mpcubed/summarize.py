# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os
import glob

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

from . import io, misc
from .detection import transit_search as ts

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

def fig_transmission(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_trans.png'
        figname = os.path.join(head, figname) 
    
    # Read the transmission map.
    f = io.SysFile(filename)
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
    
def fig_intrapix(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_ipx.png'
        figname = os.path.join(head, figname)   
    
    # Read the intrapixel amplitudes.
    f = io.SysFile(filename) 
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
    
def fig_clouds(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_clouds.png'
        figname = os.path.join(head, figname)
    
    # Read the data.
    f = io.SysFile(filename)
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
    
def fig_sigma(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        figname = tail.rsplit('.')[0] + '_sigma.png'
        figname = os.path.join(head, figname)
    
    # Read the data.
    f = io.SysFile(filename)
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
    
    fig_transmission(filename)
    fig_intrapix(filename)
    fig_clouds(filename)
    fig_sigma(filename)
    
    return

def plot_periodogram(freq, dchisq, period, zoom=False):
    """ Plot the box least-squares periodogram. """

    plt.annotate(r'$P = {:.5f}$'.format(period), (0, 1), xytext=(10, -10), xycoords='axes fraction', textcoords='offset points', va='top', ha='left', size='x-large', backgroundcolor='w')    
    
    # Plot the box least-squares periodogram.
    plt.plot(freq, dchisq, c='k')
    if zoom:
        plt.xlim(.95/period, 1.05/period)
    else:
        plt.xlim(0, 1.8)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
    
    # Add lines indicating the peak, and integer harmonics.
    freq = 1/period
    plt.axvline(freq, c=(0./255,109./255,219./255))
    for n in range(2, 5):
        plt.axvline(n*freq, c=(0./255,109./255,219./255), ls='--')
        plt.axvline(freq/n, c=(0./255,109./255,219./255), ls='--')
        
    # Add lines indicating the 1 day systematic and harmonics.
    freq = 1/.9972
    plt.axvline(freq, c=(146./255,0,0), lw=2, ymax=.1)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline(freq/n, c=(146./255,0,0), lw=2, ymax=.1)
        
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline((n+1)*freq/n, c=(146./255,0,0), lw=2, ymax=.1)     
    
    return  
    
def plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration, binned=True, zoom=False):
    """ Plot a lightcurve with the box least-squares best fit overplotted. """    
    
    factor = np.ceil(np.abs(depth)/.05)   

    # Plot the phase-folded lightcurve.
    phase = (jdmid - epoch)/period    
    phase = np.mod(phase+.5, 1.)-.5
    plt.scatter(phase, mag, color='black', marker='.', alpha=.5, edgecolor='none')
    
    # Add phase-binned data.
    if binned:
        nbins = np.ceil(9*period/duration)
        bins = np.linspace(-.5, .5, nbins+1)   
        xbin, ybin, eybin = misc.bin_data_err(phase, mag, emag, bins)
        plt.errorbar(xbin, ybin, eybin, fmt='o', c=(0./255,109./255,219./255))
    
    if zoom:
        plt.xlim(-1.5*duration/period, 1.5*duration/period)
        plt.ylim(2*np.abs(depth), -2*np.abs(depth))
    else:
        plt.xlim(-.5, .5)
        plt.ylim(factor*.05, factor*-.03)
    
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$') 
    
    return
    
def plot_boxcurve(period, depth, duration):
    
    # Compute x and y values for plotting a boxfit.
    x = .5*duration/period        
    x = np.array([-.5, -x, -x, x, x, .5])
    y = np.array([0, 0, depth, depth, 0, 0,])
    
    plt.plot(x, y, c=(146./255,0,0), lw=2, zorder=20)    
    
    return

def fig_candidate(ascc, freq, dchisq, jdmid, mag, emag, period, epoch, depth, duration, figdir):
    
    # Create the figure.
    fig = plt.figure(figsize=(16.5, 11.7))
    
    plt.suptitle('ASCC {}'.format(ascc), size='xx-large')    
    
    gs = gridspec.GridSpec(6, 6, height_ratios = [.5, 10, 10, .5, 10, .5])
    
    # Plot the periodogram.
    ax = plt.subplot(gs[1,:4])
    
    plt.title('Periodogram')
    plot_periodogram(freq, dchisq, period)                     
    
    ax = plt.subplot(gs[1,4:])
    
    plt.title('Periodogram zoom')
    plot_periodogram(freq, dchisq, period, zoom=True)    
    
    # Plot the data and the best-fit box-model.
    ax = plt.subplot(gs[2:4,:4])

    plt.title('Photometry')
    plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration)
    plot_boxcurve(period, depth, duration)
    
    ax = plt.subplot(gs[2:4,4:])

    plt.title('Photometry zoom')
    plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration, zoom=True)
    plot_boxcurve(period, depth, duration)

    ax = plt.subplot(gs[4,:3])

    plt.title(r'Half-period, $P = {:.5f}$'.format(.5*period))
    plot_lightcurve(jdmid, mag, emag, .5*period, epoch, depth, duration, binned=False)

    ax = plt.subplot(gs[4,3:])

    plt.title(r'Double-period, $P = {:.5f}$'.format(2.*period))
    plot_lightcurve(jdmid, mag, emag, 2.*period, epoch, depth, duration, binned=False)        
        
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, 'ASCC{}.png'.format(ascc)))
    plt.close()    
    
    return
 
  
def boxlstsq_summary(blsdir, aper=0, method='legendre'): 

    blsfiles = glob.glob(os.path.join(blsdir, 'bls/*'))
    blsfiles = np.sort(blsfiles)
    
    filelist = np.genfromtxt(os.path.join(blsdir, 'data.txt'), dtype='S') 

    figdir = os.path.join(blsdir, 'figures')
    misc.ensure_dir(figdir)

    for blsfile in blsfiles:
        
        print blsfile

        # Read the box least-squares results.
        f = io.blsFile(blsfile)
        hdr = f.read_header(['ascc', 'period', 'epoch', 'depth', 'duration'])  
        data = f.read_data(['freq', 'dchisq'])
    
        # Select stars.
        args, = np.where(hdr['depth'] > 0)
        
        # Read the lightcurves.
        jdmid, lst, mag, emag, trend, mask = ts.read_data(filelist, hdr['ascc'], aper=aper, method=method)
        mask = ~mask
        
        for i in args:
            
            # Make the figure.
            fig_candidate(hdr['ascc'][i], data['freq'], data['dchisq'][:,i], jdmid[mask[i]], mag[i,mask[i]], emag[i,mask[i]], hdr['period'][i], hdr['epoch'][i], hdr['depth'][i], hdr['duration'][i], figdir)     

    return

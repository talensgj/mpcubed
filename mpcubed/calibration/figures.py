# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os
import argparse

import numpy as np

from astropy import wcs

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

fiducial = dict()
fiducial['LSC'] = np.array([-0.1, -32.5, 2004., 1336., 183.2])
fiducial['LSN'] = np.array([1.4, 12.1, 2004., 1336., 1.7])
fiducial['LSE'] = np.array([-45.1, -19.8, 2004., 1336., 291.6])
fiducial['LSS'] = np.array([-3.6, -69.7, 2004., 1336., 184.7])
fiducial['LSW'] = np.array([45.0, -23.9, 2004., 1336., 72.6])
fiducial['LPC'] = np.array([0.24, 28.65, 2054.0, 1356.0, 181.5])
fiducial['LPN'] = np.array([3.9, 69.8, 2056, 1273, 4.3])
fiducial['LPE'] = np.array([314.8, 22.3, 2004, 1329, 252.0])
fiducial['LPS'] = np.array([359.0, -12.1, 2024, 1287, 181.3])
fiducial['LPW'] = np.array([44.5, 20.0, 2030, 1280, 111])


def initial_wcs(pars, scale=9e-3/24., lst=0.):

    # Extract parameters.
    ha0, dec0 = pars[0], pars[1]
    x0, y0 = pars[2], pars[3]
    th0 = pars[4]

    # Convert HA to RA.
    ra0 = np.mod(15. * lst - ha0, 360.)

    # Compute PC matrix.
    sn = np.sin(th0 * np.pi / 180.)
    cs = np.cos(th0 * np.pi / 180.)
    pc0 = np.array([[cs, -sn], [sn, cs]])

    # Convert scale to degrees.
    scale = np.rad2deg(scale)

    # Generate WCS object.
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [x0, y0]
    w.wcs.cdelt = [scale, scale]
    w.wcs.crval = [ra0, dec0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.pc = pc0

    return w


def _hadec2xy(pars, ha, dec):
            
    # Perfrom the coordinate transformations.
    w = initial_wcs(pars)

    ra = np.mod(0. - ha, 360.)
    x, y = w.all_world2pix(ra, dec, 0)
        
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
    
    data = data[1:-1, 1:-1]
    data = np.ma.array(data, mask=np.isnan(data))    
    
    ha, dec = grid.xedges, grid.yedges
    ha, dec = np.meshgrid(ha, dec)
    x, y = _hadec2xy(wcspars, ha, dec)
    
    mask = np.isfinite(x) & np.isfinite(y)
    x[~mask] = 0
    y[~mask] = 0
    
    im = plt.pcolormesh(x, y, data.T, **kwargs)
    wcsgrid(wcspars)

    return im


def fig_magnitudes(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        path = os.path.join(head, 'figures')
        io.ensure_dir(path)
        figname = tail.rsplit('.')[0] + '_mags.png'
        figname = os.path.join(path, figname)
        
    # Read the magnitudes.
    f = io.SysFile(filename)
    magnitudes = f.read_magnitudes()
    mag, vmag, sigma = magnitudes['mag'], magnitudes['vmag'], magnitudes['sigma']
    
    # Plot the magitudes.
    fig = plt.figure(figsize=(14, 9))
    
    plt.suptitle('Magnitudes', size='xx-large')
    
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 5, 5])
    
    plt.subplot(gs[1, 0])
    
    plt.scatter(vmag, mag - vmag, color='black', marker='.', alpha=.5, edgecolor='none')
    
    plt.ylabel('V - M')
    
    plt.subplot(gs[2, 0])
    
    plt.scatter(vmag, sigma, color='black', marker='.', alpha=.5, edgecolor='none')
    
    plt.ylim(0, 0.5)
    
    plt.xlabel('V')
    plt.ylabel(r'$\sigma$')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()
    
    return figname


def fig_transmission(filename, wcspars, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        path = os.path.join(head, 'figures')
        io.ensure_dir(path)
        figname = tail.rsplit('.')[0] + '_trans.png'
        figname = os.path.join(path, figname)
    
    # Read the transmission map.
    f = io.SysFile(filename)
    
    if isinstance(wcspars, str):
        wcspars = fiducial[wcspars]
        
    camgrid, trans = f.read_trans()
    
    trans = trans['trans']
    
    vmin = np.nanpercentile(trans[:, 10:-10], 0.1)
    vmax = np.nanpercentile(trans[:, 10:-10], 99.9)
    
    # Plot the transmission map.
    fig = plt.figure(figsize=(14, 9))
    
    plt.suptitle('Transmission', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios=[15, .5], height_ratios=[1, 10])
    
    plt.subplot(gs[1, 0], aspect='equal')
    
    im = plot_polar(camgrid, trans, wcspars, cmap=plt.cm.viridis, vmin=vmin, vmax=vmax)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1, 1])
    cb = plt.colorbar(im, cax=cax)
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()
    
    return figname


def fig_intrapix(filename, wcspars, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        path = os.path.join(head, 'figures')
        io.ensure_dir(path)
        figname = tail.rsplit('.')[0] + '_ipx.png'
        figname = os.path.join(path, figname)  
    
    # Read the intrapixel amplitudes.
    f = io.SysFile(filename)

    if isinstance(wcspars, str):
        wcspars = fiducial[wcspars]
        
    ipxgrid, intrapix = f.read_intrapix() 
    
    amplitudes = intrapix['amplitudes']       
   
    vlim = 0.1

    # Plot the amplitudes.
    fig = plt.figure(figsize=(16, 10))
    
    plt.suptitle('Intrapixel Amplitudes', size='xx-large')
    
    gs = gridspec.GridSpec(3, 3, width_ratios=[15, 15, .5], height_ratios=[1, 10, 10])
    
    plt.subplot(gs[1, 0], aspect='equal')
    plt.title(r'$\sin(2\pi x)$')
   
    im = plot_polar(ipxgrid, amplitudes[:, :, 0], wcspars, cmap=plt.cm.coolwarm, vmin=-vlim, vmax=vlim)
   
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[1, 1], aspect='equal')
    plt.title(r'$\cos(2\pi x)$')
    
    im = plot_polar(ipxgrid, amplitudes[:, :, 1], wcspars, cmap=plt.cm.coolwarm, vmin=-vlim, vmax=vlim)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2, 0], aspect='equal')
    plt.title(r'$\sin(2\pi y)$')
    
    im = plot_polar(ipxgrid, amplitudes[:, :, 2], wcspars, cmap=plt.cm.coolwarm, vmin=-vlim, vmax=vlim)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2, 1], aspect='equal')
    plt.title(r'$\cos(2\pi y)$')
    
    im = plot_polar(ipxgrid, amplitudes[:, :, 3], wcspars, cmap=plt.cm.coolwarm, vmin=-vlim, vmax=vlim)
    
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1:, 2])
    cb = plt.colorbar(im, cax=cax)
    cb.set_label('Amplitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()  
    
    return figname


def fig_clouds(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        path = os.path.join(head, 'figures')
        io.ensure_dir(path)
        figname = tail.rsplit('.')[0] + '_clouds.png'
        figname = os.path.join(path, figname)
    
    # Read the data.
    f = io.SysFile(filename)
    skygrid, clouds = f.read_clouds()
    
    clouds = clouds['clouds']
    
    # Remove empty rows and columns.
    mask = np.isfinite(clouds)
    q, = np.where(np.any(mask, axis=1))
    t, = np.where(np.any(mask, axis=0))
    clouds = clouds[q][:, t]
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Clouds', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios=[15, .5], height_ratios=[1, 20])
    
    plt.subplot(gs[1, 0], xticks=[], yticks=[])
        
    im = plt.imshow(clouds.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=-.1, vmax=1.5)
    
    idx, = np.where(np.diff(t) > 50)
    for i in idx:
        plt.axhline(i+1, xmax=.1, c='k', lw=2)
    
    plt.ylabel('Time')
    plt.xlabel('Sky Patch')
    
    cax = plt.subplot(gs[1, 1])
    cb = plt.colorbar(im, cax=cax)
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    
    plt.tight_layout()
    
    plt.savefig(figname, dpi=100)
    plt.close()    
    
    return figname


def fig_sigma(filename, figname=None):
    
    if figname is None:
        head, tail = os.path.split(filename)
        path = os.path.join(head, 'figures')
        io.ensure_dir(path)
        figname = tail.rsplit('.')[0] + '_sigma.png'
        figname = os.path.join(path, figname)
    
    # Read the data.
    f = io.SysFile(filename)
    skygrid, clouds = f.read_clouds()
    
    sigma = clouds['sigma']
    
    # Remove empty rows and columns.
    mask = np.isfinite(sigma)
    idx1, = np.where(np.any(mask, axis=1))
    idx2, = np.where(np.any(mask, axis=0))
    sigma = sigma[idx1][:, idx2]
    
    # Plot the clouds.
    fig = plt.figure(figsize=(10, 16))
    
    plt.suptitle('Sigma', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios=[15, .5], height_ratios=[1, 20])
    
    plt.subplot(gs[1, 0], xticks=[], yticks=[])
      
    im = plt.imshow(sigma.T, interpolation='None', aspect='auto', cmap=plt.cm.viridis, vmin=0, vmax=0.5)
    
    idx, = np.where(np.diff(idx2) > 50)
    for i in idx:
        plt.axhline(i+1, xmax=.1, c='k', lw=2)
    
    plt.ylabel('Time')
    plt.xlabel('Sky Patch')
    
    cax = plt.subplot(gs[1, 1])
    cb = plt.colorbar(im, cax=cax)
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

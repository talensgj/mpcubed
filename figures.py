#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

import pandas as pd

# For A4 paper.
rcParams['figure.figsize'] = (11.69,8.27)
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

def transit_duration():
    
    import boxlstsq
    
    Mass = np.array([60., 18., 3.2, 1.7, 1.1, .8, .3])
    Radius = np.array([15., 7., 2.5, 1.3, 1.1, .9, .4])
    
    # O star.
    M = 60.
    R = 15.
    fmax = boxlstsq.freq_max(M, R)
    P = np.linspace(1/fmax, 20., 100)
    q = boxlstsq.phase_duration(1/P, M, R)
    plt.plot(P, q, label='O star', lw=2)
    
    # M star.
    M = .3
    R = .4
    fmax = boxlstsq.freq_max(M, R)
    P = np.linspace(1/fmax, 20., 100)
    q = boxlstsq.phase_duration(1/P, M, R)
    plt.plot(P, q, label='M star', lw=2)
    
    # Sun + code.
    M = 1.
    R = 1.
    fmax = boxlstsq.freq_max(M, R)
    P = np.linspace(1/fmax, 20., 100)
    q = boxlstsq.phase_duration(1/P, M, R)
    plt.plot(P, q, label='Sun', c='k', lw=2)
    plt.plot(P, q/3, ls='--', c='k')
    plt.plot(P, 3*q, ls='--', c='k')
    
    plt.xlim(0, 20)
    plt.ylim(0, .35)
    plt.legend()
    plt.xlabel('P [days]')
    plt.ylabel(r'$\eta/P$')
    plt.show()
    
    return
    
    
def magnitude_parallax():
    
    m = np.linspace(2, 8.5, 500)
    p = np.linspace(0, 30, 500)
    
    m, p = np.meshgrid(m, p)
    d = 1/(.001*p)
    
    tmp = m - 4.83 - 5.*np.log10(d/10.)
    logL = -tmp/2.5
    
    ax = plt.subplot(111)
    ax.invert_xaxis()
    cs = plt.contour(m, p, logL, [-2, -1, 0, 1, 2, 3], colors='k')
    plt.clabel(cs, manual=True)
    plt.xlim(8.4, 2)
    plt.ylim(0, 30)
    plt.xlabel('V')
    plt.ylabel('Parallax [mas]')
        
    plt.show()
        
    return

def depth_radiusstar():
    
    delta = np.linspace(0, .05, 250)
    Rstar = np.linspace(0, 5, 500)
    
    delta, Rstar = np.meshgrid(delta, Rstar)
    
    Rplanet = np.sqrt(delta)*Rstar*10.05
    
    ax = plt.subplot(111)
    cs = plt.contour(Rstar, delta*100, Rplanet, [1, 2, 3, 4, 5], colors='k')
    plt.clabel(cs, manual=True)
    plt.xlim(0, 5)
    plt.ylim(0, 5)
    plt.xlabel(r'$R_*$ [$R_\odot$]')
    plt.ylabel(r'$\delta$ [%]')
    
    plt.show()
   
def coverage_LaPalma():
    
    import healpy
    from mascara import observer
    from package.plotting import viridis
    
    npix = healpy.nside2npix(8)
    
    ccdx = np.hstack([np.zeros(167*2), np.linspace(0, 4008, 167*3), 4008*np.ones(167*2), np.linspace(0, 4008, 167*3)[::-1]])
    ccdy = np.hstack([np.linspace(0, 2672, 167*2), 2672*np.ones(167*3), np.linspace(0, 2672, 167*2)[::-1], np.zeros(167*3)])
    
    alt0 = [90, 49, 49, 49, 49]
    az0 = [0, 0, 90, 180, 270]
    
    site = observer.Site('LaPalma')
    ha0, dec0 = site.altaz2hadec(np.array([90.]), np.array([0.]))
    
    cmap = viridis
    cmap.set_bad('w')
    cmap.set_under('w')
    
    fig = plt.figure(figsize=(8,8))
    plt.subplot(111)
    
    healpy.orthview(np.full(npix, fill_value=np.nan), fig=0, hold=True, rot=(ha0, dec0), half_sky=True, title='', cbar=False, cmap=cmap)
    for alt, az in zip(alt0, az0):
    
        cam = observer.Camera(altitude=alt, azimuth=az, orientation=270., Xo=2004, Yo=1336, nx=4008, ny=2672)
        phi, the = cam.XY2PhiThe(ccdx, ccdy)
        alt, az = cam.PhiThe2Hor(phi, the)
        ha, dec = site.altaz2hadec(alt, az)
        
        healpy.projplot((90-dec)*np.pi/180, ha*np.pi/180, c=(146./255,0,0), lw=2)
    healpy.graticule(10., 15.)
    
    plt.show()
    plt.close()
    
    return
    
def coverage_night(filelist):
    
    import healpy
    from astropy.io import fits
    
    from mascara import observer
    from package.coordinates import grids
    from package.plotting import viridis
    
    cameras = {'N':0, 'E':1, 'S':2, 'W':3, 'C':4}
    mask = np.full((5, 13500), fill_value=np.nan)
    jdmid = np.zeros(13500)
    lst = np.zeros(13500)
    
    for filename in filelist:
        
        print filename
        
        with fits.open(filename) as hdulist:
            cam = hdulist[0].header['CAMERA']
            lstidx = hdulist[0].header['LPSTIDX']
            lst[lstidx] = hdulist[0].header['LPST']
            jdmid[lstidx] = hdulist[0].header['JD']
    
        mask[cameras[cam], lstidx] = 1
    
    hg = grids.HealpixGrid(64)
    coverage = np.zeros(hg.npix)
    
    ccdx = np.linspace(0, 4008, 3*167)
    ccdy = np.linspace(0, 2672, 2*167)
    
    ccdx, ccdy = np.meshgrid(ccdx, ccdy)
    ccdx = ccdx.ravel()
    ccdy = ccdy.ravel()
    
    alt0 = [49, 49, 49, 49, 90]
    az0 = [0, 90, 180, 270, 0]
    
    site = observer.Site('LaPalma')
    ha0, dec0 = site.altaz2hadec(np.array([90.]), np.array([0.]))
    for i in range(5):
        
        cam = observer.Camera(altitude=alt0[i], azimuth=az0[i], orientation=270., Xo=2004, Yo=1336, nx=4008, ny=2672)
        phi, the = cam.XY2PhiThe(ccdx, ccdy)
        alt, az = cam.PhiThe2Hor(phi, the)
        ha, dec = site.altaz2hadec(alt, az)
        
        for j in range(13500):
            
            if (mask[i,j] > 0):
                
                ra = np.mod(lst[j]*15 - ha, 360)
                #idx = hg.radec2idx(ha, dec)
                idx = hg.radec2idx(ra, dec)
                coverage[np.unique(idx)] += 1
    
    coverage[coverage == 0] = np.nan
      
    fig = plt.figure(figsize=(10,6))
    gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,2], width_ratios = [20,.5])
    plt.suptitle('Coverage 30-01-2016', size='xx-large')
      
    cmap = viridis
    cmap.set_under('w')
      
    # Plot the coverage.
    ax = plt.subplot(gs[1,0])
    #healpy.orthview(coverage, rot=(ha0, dec0), half_sky=True, fig=0, sub=211)
    healpy.mollview(coverage, hold=True, title='', cbar=False, cmap=viridis)
    healpy.graticule(10., 15.)
    
    ax = plt.gca()
    image = ax.get_images()[0]
    ax = plt.subplot(gs[1,1])
    cb = plt.colorbar(image, cax=ax)
    cb.set_label('# Images')
    
    plt.subplot(gs[2,0])
    plt.imshow(mask, aspect='auto', interpolation='none', cmap=plt.cm.Greys_r, origin='upper')
    plt.yticks([0, 1, 2, 3, 4], ['N', 'E', 'S', 'W', 'C'])
    plt.xlabel('Time [LST-IDX]')

    #plt.tight_layout()
    plt.show()
    plt.close()
    
    return

    
def main():
    
    #transit_duration()
    #magnitude_parallax()
    #depth_radiusstar()
    coverage_LaPalma()
    
    #filelist = glob.glob('/data2/talens/LaPalma/raw/20160130LP?/*LP?.fits')
    #filelist = filelist[::50]
    #coverage_night(filelist)
    
    return

if __name__ == '__main__':
    main()

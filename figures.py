#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt


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
    
    npix = healpy.nside2npix(8)
    
    ccdx = np.hstack([np.zeros(167*2), np.linspace(0, 4008, 167*3), 4008*np.ones(167*2), np.linspace(0, 4008, 167*3)[::-1]])
    ccdy = np.hstack([np.linspace(0, 2672, 167*2), 2672*np.ones(167*3), np.linspace(0, 2672, 167*2)[::-1], np.zeros(167*3)])
    
    alt0 = [90, 49, 49, 49, 49]
    az0 = [0, 0, 90, 180, 270]
    
    site = observer.Site('LaPalma')
    ha0, dec0 = site.altaz2hadec(np.array([90.]), np.array([0.]))
    
    healpy.orthview(np.zeros(npix), rot=(ha0, dec0), half_sky=True)
    for alt, az in zip(alt0, az0):
    
        cam = observer.Camera(altitude=alt, azimuth=az, orientation=270., Xo=2004, Yo=1336, nx=4008, ny=2672)
        phi, the = cam.XY2PhiThe(ccdx, ccdy)
        alt, az = cam.PhiThe2Hor(phi, the)
        ha, dec = site.altaz2hadec(alt, az)
        
        healpy.projplot((90-dec)*np.pi/180, ha*np.pi/180)
        healpy.graticule(10., 15.)
    
    plt.show()
    
    return
    
def main():
    
    #transit_duration()
    #magnitude_parallax()
    #depth_radiusstar()
    coverage_LaPalma()
    
    return

if __name__ == '__main__':
    main()

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
    
def main():
    
    #transit_duration()
    magnitude_parallax()
    #depth_radiusstar()
    
    return

if __name__ == '__main__':
    main()

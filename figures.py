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
    
    
def main():
    
    transit_duration()
    
    return

if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt

def fourier_fit(t, y, freq):
    
    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)
        
    sin_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        sin_mat[:,i] = np.sin(2*np.pi*freq[i]*t)        
    
    mat = np.hstack([sin_mat, cos_mat])
    pars = np.linalg.lstsq(mat, y)[0]
    
    fit = np.dot(mat, pars)
    
    return pars, fit

def main():
    
    npoints = 500
    BL = 5.
    P = .33
    
    t = np.linspace(0, BL, npoints)
    y = .1*np.sin(2*np.pi*t/P)+.1*np.sin(2*np.pi*4*t/P)
    
    freq = np.fft.rfftfreq(npoints, BL/(npoints-1))
    Y = np.fft.rfft(y)
    Y = np.abs(Y)
    Y = 2*Y/npoints # Factor two because real fft.
    
    pars, fit = fourier_fit(t, y, freq)
    P = np.sqrt(pars[:len(pars)/2]**2 + pars[len(pars)/2:]**2)
    
    print np.amin(1/freq[freq>0])
    
    plt.subplot(311)
    plt.plot(t, y)
    plt.plot(t, fit)
    
    plt.subplot(312)
    plt.plot(freq, Y)
    
    plt.subplot(313)
    plt.plot(freq, P)
    
    plt.show()
    
    return
    
if __name__ == '__main__':
    main()
    

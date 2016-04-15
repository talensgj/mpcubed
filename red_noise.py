#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt

from package.models import transit
import boxlstsq

def rolling_mean(idx, y, window):
    
    nbin = np.bincount(idx)
    ybin = np.bincount(idx, y)
    
    nbin = np.cumsum(nbin)
    ybin = np.cumsum(ybin)
    
    nbin = np.append(0, nbin)
    ybin = np.append(0, ybin)
    
    nbin = (nbin[window:] - nbin[:-window])
    ybin = (ybin[window:] - ybin[:-window])/nbin
    
    return nbin, ybin
    
def transits(condition):
    
    length,count = [],0
    for i in range(len(condition)):

        if condition[i]==True:
            count += 1
        elif condition[i]==False and count>0:
            length.append(count)
            count = 0

        if i==len(condition)-1 and count>0:
            length.append(count)

    return length

def main(args):
    
    np.random.seed(19910909)
    
    idx = np.arange(5000)
    t = np.linspace(0, 40, 5000)
    #y = np.ones(5000) + .01*np.random.randn(5000) + .005*np.sin(2*np.pi*5000.*t)
    y = transit.box_model(t, 2.54, .23, .02, .1) + .01*np.random.randn(5000) + .005*np.sin(2*np.pi*5000.*t)
    yerr = .01*np.ones(5000)
    
    for i in range(20, 5000, 200):
    
        idx = np.delete(idx, np.arange(i, i+50))
        t = np.delete(t, np.arange(i, i+50))
        y = np.delete(y, np.arange(i, i+50))
        yerr = np.delete(yerr, np.arange(i, i+50))
    
    weights = 1/yerr**2
    freq, chisq0, dchisq, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(t, y, weights)
    
    arg = np.argmax(dchisq)
    print dchisq[arg], 1/freq[arg], depth[arg], duration[arg]
    
    phase = freq[arg]*(t - epoch[arg])
    phase = np.mod(phase+.5, 1)-.5
    
    plt.subplot(211)
    plt.plot(freq, dchisq)
    plt.subplot(212)
    plt.errorbar(phase, y, yerr, fmt='.')
    plt.show()
    
    phase = freq[arg]*(t - epoch[arg])
    phase = np.mod(phase + .5, 1)-.5
    mask = np.abs(phase) < .5*freq[arg]*duration[arg]
    
    plt.plot(t[mask], y[mask], '.')
    plt.plot(t[~mask], y[~mask], '.')

    plt.show()
    
    window = np.ceil(5000/20.*duration[arg]).astype('int')
    print window
    npbin, ybin = rolling_mean(idx[~mask], y[~mask], window)
    npbin = np.squeeze(npbin)
    ybin = np.squeeze(ybin)
    
    m0 = np.bincount(npbin)
    m1 = np.bincount(npbin, ybin)
    m2 = np.bincount(npbin, ybin**2)
    
    var = m2/m0 - (m1/m0)**2
    
    plt.plot(npbin, ybin, '.')
    plt.xlim(.5, 25.5)
    plt.show()
    
    plt.plot(np.sqrt(var))
    plt.plot(.01/np.sqrt(np.arange(0,25)))
    plt.show()
    
    nt = transits(mask)
    n = np.sum(nt)
    Nt = len(nt)
    
    nt = np.array(nt)
    print var[nt]
    Sred = (depth[arg])**2*n**2/(np.sum(nt**2.*var[nt]))
    
    print Sred
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

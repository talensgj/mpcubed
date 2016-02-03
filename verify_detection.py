#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m] 

filename = ''
redfile = ''

with h5py.File(filename, 'r') as f:
    ascc = f.keys()
    
for i in range(len(ascc)):
    
    with h5py.File(filename, 'r') as f:
        grp = f[ascc[i]]
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        depth = grp['depth'].value
        nobs = grp.attrs['nobs']
        chisq0 = grp.attrs['chisq0']
    
    with h5py.File(redfile, 'r') as f:
        grp = f['data/'+ascc[i]]
        jdmid = grp['jdmid'].value
        
    
    # Best fit.
    arg = np.argmax(dchisq)
    
    best_chisq = chisq0 - dchisq[arg]
    best_freq = freq[arg]
    best_depth = depth[arg]
    
    if (best_depth < 0): # check signs of depth.
        print 'Best fit is actually a brightening.'
    
    # Variable stars.
    if (best_chisq > 3.5*nobs):
        print 'Large best-fit chisq, possible variable star...'
        
    # Anti-transit ratio:
    args1, = np.where(depth < 0)
    args2, = np.where(depth > 0)
    ratio = np.max(dchisq[args1])/np.max(dchisq[args2])
    
    if (ratio > 1.5): # Check signs of depth.
        print 'Warning poor anti-transit ratio.'    
        
    # Data gaps.
    #q = ((2.*np.pi)**(2./3)*R*freq**(2./3))/(np.pi*(G*M)**(1./3)) # Check units.
    #phase = np.mod(jdmid*best_freq, 1)
    #nbins = 2*np.ceil(1/q)
    #bins = np.linspace(0, 1, nbins)
    #idx = np.searchsorted(bins, phase, 'right') # Has out of range padding.
    #ppbin = np.bincount(idx, minlength=nbins)
    #ranges = misc.zero_runs(ppbin)
    #lengths = ranges[:,1] - ranges[:,0]
    #if np.any(lengths >= 5):
        #print 'Gap of 2.5 times the transit duration folded lightcurve.'
        
    #if (ranges[0,0] == 0) & (ranges[-1,1] == len(jdmid))
        #length = lengths[0] + length[-1]
        #if length >= 5:
            #print 'Gap of 2.5 times the transit duration folded lightcurve.'

    
        
        
        
        
    
    
    
    
    
    
    

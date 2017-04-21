#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

sptype_OC = {'O5':(13.4, -5.1), 'O6':(12.2, -5.1), 'O7':(11., -4.9),
             'O8':(10., -4.6), 'B0':(6.7, -3.4), 'B1':(5.2, -2.6),
             'B2':(4.1, -1.6), 'B3':(3.8, -1.3), 'B5':(3.2, -.5),
             'B6':(2.9, -.1), 'B7':(2.7, .3), 'B8':(2.5, .6),
             'B9':(2.3, .8), 'A0':(2.2, 1.1), 'A1':(2.1, 1.3),
             'A2':(2., 1.5), 'A5':(1.8, 2.2), 'A8':(1.5, 2.7),
             'F0':(1.4, 3.), 'F2':(1.3, 3.4), 'F5':(1.2, 3.9),
             'F8':(1.1, 4.3), 'G0':(1.06, 4.7), 'G2':(1.03, 4.9),
             'G8':(.96, 5.6), 'K0':(.93, 5.7), 'K1':(.91, 6.),
             'K3':(.86, 6.5), 'K4':(.83, 6.7), 'K5':(.8, 7.1),
             'K7':(.74, 7.8), 'M0':(.63, 8.9), 'M1':(.56, 9.6),
             'M2':(.48, 10.4), 'M3':(.41, 11.1), 'M4':(.35, 11.9),
             'M5':(.29, 12.8), 'M6':(.24, 13.8), 'M7':(.20, 14.7)}

sptype_OCinp = {'O5':(13.4, -5.1), 'O6':(12.2, -5.1), 'O7':(11., -4.9),
                'O8':(10., -4.6), 'B0':(6.7, -3.4), 'B1':(5.2, -2.6),
                'B2':(4.1, -1.6), 'B3':(3.8, -1.3), 'B4':(3.5, -.9),
                'B5':(3.2, -.5), 'B6':(2.9, -.1), 'B7':(2.7, .3),
                'B8':(2.5, .6), 'B9':(2.3, .8), 'A0':(2.2, 1.1),
                'A1':(2.1, 1.3), 'A2':(2., 1.5), 'A3':(1.9, 1.7),
                'A4':(1.9, 1.9), 'A5':(1.8, 2.2), 'A6':(1.7, 2.4),
                'A7':(1.6, 2.6), 'A8':(1.5, 2.7), 'F0':(1.4, 3.),
                'F1':(1.35, 3.2), 'F2':(1.3, 3.4), 'F3':(1.3, 3.6),
                'F4':(1.2, 3.8), 'F5':(1.2, 3.9), 'F6':(1.2, 4.0),
                'F7':(1.1, 4.2), 'F8':(1.1, 4.3), 'G0':(1.06, 4.7),
                'G1':(1.04, 4.8), 'G2':(1.03, 4.9), 'G3':(1.02, 5.),
                'G4':(1.01, 5.1), 'G5':(1., 5.2), 'G6':(.99, 5.3),
                'G7':(.98, 5.4), 'G8':(.96, 5.6), 'K0':(.93, 5.7),
                'K1':(.91, 6.), 'K2':(.88, 6.3), 'K3':(.86, 6.5),
                'K4':(.83, 6.7), 'K5':(.8, 7.1), 'K6':(.77, 7.4),
                'K7':(.74, 7.8), 'M0':(.63, 8.9), 'M1':(.56, 9.6),
                'M2':(.48, 10.4), 'M3':(.41, 11.1), 'M4':(.35, 11.9),
                'M5':(.29, 12.8), 'M6':(.24, 13.8), 'M7':(.20, 14.7)}

def flux2mag(flux, eflux=None, m0=25.):
    
    mag = m0 - 2.5*np.log10(flux)
    
    if eflux is not None:
        emag = 2.5/np.log(10.)*np.abs(eflux/flux)
        return mag, emag
    
    return mag

def sigma_clip(array, sigma=5., niter=5):
    
    m0 = np.nanmean(array)
    m1 = np.nanstd(array)
    
    mask = np.abs(array - m0) < sigma*m1
    for i in range(niter):
        m0 = np.nanmean(array[mask])
        m1 = np.nanstd(array[mask])
        
        mask = np.abs(array - m0) < sigma*m1
        
    return m0, m1

def phase(time, period, time_ref=0., fold=True):
    
    phase = (time - time_ref)/period
    
    if fold:
        phase = np.mod(phase, 1.)
        
    return phase

def ensure_dir(path):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    return

def find_ns(lstseq):
    """ Number of sampling points in LST, takes wrapping into account."""
    
    lstidx = (lstseq % 270)
    option1 = np.ptp(lstidx) + 1
    
    lstidx = np.mod(lstidx + 135, 270)
    option2 = np.ptp(lstidx) + 1
    
    if (option2 >= option1):
        return option1, False
    else:
        return option2, True

def transits(jdmid, period, epoch, duration):
    
    phase = (jdmid - epoch)/period
    cycle = np.floor(phase+.5)
    
    phase = np.mod(phase+.5, 1.)-.5
    condition = np.abs(phase) < .5*duration/period
    
    cycles = np.unique(cycle)
    
    length = []
    for i in cycles:
        
        sel = (cycle == i)
        nt = np.sum(condition[sel])

        if (nt > 0):
            length.append(nt)

    return length
    
def get_absmag(mag, d):
    
    Mag = mag - 5.*np.log10(d/10.)
    
    return Mag
    
def get_rstar(Mag, Mag0, R0):
    
    Rstar = R0*10**(-(Mag - Mag0)/5.)
    
    return Rstar
    
def get_rplanet(depth, Rstar):
    
    Rplanet = np.sqrt(depth/.01)*Rstar
    
    return Rplanet
    
def round_to_significance(value, error1, error2=None):
    
    ndigits = -int(np.floor(np.log10(error1)))  
    
    if error2 is None:    
        return ndigits, round(value, ndigits), round(error1, ndigits)
    
    ndigits = np.maximum(ndigits, -int(np.floor(np.log10(-error2)))) 
    
    return ndigits, round(value, ndigits), round(error1, ndigits),  round(error2, ndigits)    
    
#def transits(jdmid, pars):
    
    #phase = (jdmid - pars[1])/pars[0]
    #orbit = np.floor(phase + .5).astype('int')
    
    #phase = np.mod(phase+.5, 1.)-.5
    #sel = (np.abs(phase) < .5*pars[3]/pars[0])
    
    #intransit = np.bincount(orbit[sel])
    #intransit = intransit[intransit > 0]
    
    #return intransit

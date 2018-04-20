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

def barycentric_dates(jdmid, ra, dec, site='Roque de los Muchachos'):
    
    from astropy import time, coordinates, units    
    
    star = coordinates.SkyCoord(ra, dec, frame='icrs', unit=[units.deg, units.deg])
    site = coordinates.EarthLocation.of_site(site)
    
    times = time.Time(jdmid, format='jd', scale='utc', location=site)
    
    ltt_bary = times.light_travel_time(star)
    times_barycentre = times.tdb + ltt_bary
    jdbar = times_barycentre.jd
    
    return jdbar

def flux2mag(flux, eflux=None, m0=25.):
    
    mag = m0 - 2.5*np.log10(flux)
    
    if eflux is not None:
        emag = 2.5/np.log(10.)*np.abs(eflux/flux)
        return mag, emag
    
    return mag

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

def bin_data_noerr(x, y, bins):
    
    npbin, bins = np.histogram(x, bins=bins)
    y_sum, bins = np.histogram(x, weights=y, bins=bins)
    ysq_sum, bins = np.histogram(x, weights=y**2., bins=bins)

    x_bin = (bins[:-1] + bins[1:])/2.
    y_bin = y_sum/npbin
    ey_bin = np.sqrt(ysq_sum/npbin - (y_sum/npbin)**2)/np.sqrt(npbin)

    mask = np.isfinite(y_bin)

    return x_bin[mask], y_bin[mask], ey_bin[mask]

def bin_data_err(x, y, yerr, bins):
    
    weights = 1./yerr**2

    w_sum, bins = np.histogram(x, weights=weights, bins=bins)
    wy_sum, bins = np.histogram(x, weights=weights*y, bins=bins)
    
    x_bin = (bins[:-1] + bins[1:])/2.
    y_bin = wy_sum/w_sum    
    ey_bin = 1/np.sqrt(w_sum)    
    
    mask = np.isfinite(y_bin)    
    
    return x_bin[mask], y_bin[mask], ey_bin[mask]

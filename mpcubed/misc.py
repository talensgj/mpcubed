#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def flux2mag(flux, eflux=None, m0=25.):
    
    mag = m0 - 2.5*np.log10(flux)
    
    if eflux is not None:
        emag = 2.5/np.log(10.)*np.abs(eflux/flux)
        return mag, emag
    
    return mag

def barycentric_dates(jdmid, ra, dec, site='Roque de los Muchachos'):
    
    from astropy import time, coordinates, units    
    
    star = coordinates.SkyCoord(ra, dec, frame='icrs', unit=[units.deg, units.deg])
    site = coordinates.EarthLocation.of_site(site)
    
    times = time.Time(jdmid, format='jd', scale='utc', location=site)
    
    ltt_bary = times.light_travel_time(star)
    times_barycentre = times.tdb + ltt_bary
    jdbar = times_barycentre.jd
    
    return jdbar

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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def altaz2hadec(alt, az, lat):
    
    d2r = np.pi/180.
    alt_r = alt*d2r
    az_r = az*d2r
    lat_r = lat*d2r

    #******************************************************************************
    # find local HOUR ANGLE (in degrees, from 0. to 360.)
    ha = np.arctan2(-np.sin(az_r)*np.cos(alt_r), -np.cos(az_r)*np.sin(lat_r)*np.cos(alt_r) + np.sin(alt_r)*np.cos(lat_r))
    ha = ha/d2r
    ha = np.mod(ha, 360.)

    # Find declination (positive if north of Celestial Equator, negative if south)
    sindec = np.sin(lat_r)*np.sin(alt_r) + np.cos(lat_r)*np.cos(alt_r)*np.cos(az_r)
    dec = np.arcsin(sindec)/d2r  # convert dec to degrees
    
    return ha, dec

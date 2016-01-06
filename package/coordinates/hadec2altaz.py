#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def hadec2altaz(ha, dec, lat, WS=False):
    
    d2r = np.pi/180.

    sh = np.sin(ha*d2r)
    ch = np.cos(ha*d2r)
    sd = np.sin(dec*d2r)
    cd = np.cos(dec*d2r)
    sl = np.sin(lat*d2r)
    cl = np.cos(lat*d2r)

    x = -ch*cd*sl + sd*cl
    y = -sh*cd
    z = ch*cd*cl + sd*sl
    r = np.sqrt(x**2 + y**2)
    # now get Alt, Az

    az = np.arctan2(y, x)/d2r
    alt = np.arctan2(z, r)/d2r

    # correct for negative AZ
    az = np.mod(az, 360.)

    # convert AZ to West from South, if desired
    if WS: az = np.mod(az + 180., 360.)

    return alt, az

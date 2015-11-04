#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def altaz2xy(alt, az, alt0, az0, th0, x0, y0, scale):
    
    dtor = np.pi/180.
    alt = alt*dtor
    az = az*dtor
    alt0 = alt0*dtor
    az0 = az0*dtor
    th0 = th0*dtor
    
    array = np.zeros((3, len(alt)))
    
    array[0] = np.cos(alt)*np.sin(az)
    array[1] = np.cos(alt)*np.cos(az)
    array[2] = np.sin(alt) # 
    
    M = np.zeros((3,3))
    M[0,0] = np.sin(az0)
    M[0,1] = np.cos(az0)
    M[1,0] = -np.cos(az0)
    M[1,1] = np.sin(az0)
    M[2,2] = 1.#
    
    array = np.dot(M, array)
    
    M = np.zeros((3,3))
    M[0,0] = np.sin(alt0)
    M[0,2] = -np.cos(alt0)#
    M[1,1] = 1.
    M[2,0] = np.cos(alt0)
    M[2,2] = np.sin(alt0)#
    
    array = np.dot(M, array)
    
    M = np.zeros((3,3))
    M[0,0] = np.cos(th0)
    M[0,1] = np.sin(th0)
    M[1,0] = -np.sin(th0)
    M[1,1] = np.cos(th0)
    M[2,2] = 1.#
    
    array = np.dot(M, array)
    
    x = scale*array[0]/array[2] + x0
    y = scale*array[1]/array[2] + y0
    z = array[2]
    
    #x = array[0]
    #y = array[1]
    #z = array[2]
    
    return x, y, z

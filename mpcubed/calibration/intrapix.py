#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 16:44:48 2019

@author: talens-irex
"""

import numpy as np       

def ipxmat(x, y):

    sinx, cosx = np.sin(2*np.pi*x), np.cos(2*np.pi*x)
    siny, cosy = np.sin(2*np.pi*y), np.cos(2*np.pi*y)
    b_mat = np.column_stack([sinx, cosx, siny, cosy])

    return b_mat
    
def ipxmod(amps, b_mat):
        
    ipx = np.sum(amps*b_mat, axis=1)
    
    return ipx
    
def ipxsol(mag, weights, b_mat):
    
    wsqrt = np.sqrt(weights)
    amps = np.linalg.lstsq(b_mat*wsqrt[:,np.newaxis], mag*wsqrt, rcond=None)[0]
    ipx = ipxmod(amps, b_mat)

    return amps, ipx

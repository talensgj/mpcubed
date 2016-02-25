#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def frequencies(P, nf, cnst=False):
    
    freq = np.arange(nf)/P
    
    if not cnst:
        freq = freq[1:]
        
    return freq
    
def fftfreq(step, ns, cnst=False):
    
    freq = np.fft.rfftfreq(ns, step)
    
    if not cnst:
        freq = freq[1:]
        
    return freq
    
def sin_mat(time, freq):
    
    smat = np.outer(time, freq)
    smat = np.sin(2*np.pi*smat)
        
    return smat  
      
def fit_sines(time, y, yerr, freq):
    
    mat = sin_mat(time, freq)
    
    pars = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    
    return pars
    
def evaluate_sines(time, pars, freq):
    
    mat = sin_mat(time, freq)
    
    fit = np.dot(mat, pars)
    
    return fit
      
def cos_mat(time, freq):
    
    cmat = np.outer(time, freq)
    cmat = np.cos(2*np.pi*cmat)
        
    return cmat
    
def fit_cosines(time, y, yerr, freq):
    
    mat = cos_mat(time, freq)
        
    pars = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    
    return pars
    
def evaluate_cosines(time, pars, freq):
    
    mat = cos_mat(time, freq)
    
    fit = np.dot(mat, pars)
    
    return fit
    
def fourier_mat(time, freq):
    
    smat = sin_mat(time, freq)
    cmat = cos_mat(time, freq)
    mat = np.hstack([smat, cmat])
    
    return mat
    
def fit_fourier(time, y, yerr, freq):
    
    mat = fourier_mat(time, freq)
    
    pars = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    
    return pars
    
def evaluate_fourier(time, pars, freq):
    
    mat = fourier_mat(time, freq)

    fit = np.dot(mat, pars)

    return fit
    
def fit_mat(y, yerr, mat):
    
    pars = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    
    return pars
    
def main():
    
    return

if __name__ == '__main__':
    main()

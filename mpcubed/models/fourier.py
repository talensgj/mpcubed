#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import linalg

def frequencies(P, nf, cnst=False):
    """ Create an array of harmonics.
    
    Args:
        P (float): A single period.
        nf (int): The number of harmonics.
        cnst (bool): If True freq = 0. is included, default is False.
        
    Returns:
        freq (float): An array of harmonic frequencies.
    
    """
    
    freq = np.arange(nf)/P
    
    if not cnst:
        freq = freq[1:]
        
    return freq
    
def fftfreq(step, ns, cnst=False):
    """ Create an array of fourier frequencies.
    
    Args:
        step (float): The sampling interval.
        ns (int): The number of samples.
        cnst (bool): If True freq = 0. is included, default is False.
    
    Returns:
        freq (float): An array of fourier frequencies.
    
    """
    
    freq = np.fft.rfftfreq(ns, step)
    
    if not cnst:
        freq = freq[1:]
        
    return freq
    
def sin_mat(time, freq):
    """ Create a matrix of sine waves.
    
    Args:
        time (float): Times at which to evaluate the sine waves.
        freq (float): Frequencies for the sine waves.
        
    Returns:
        smat (float): A 2d-array with smat[:,i] = np.sin(2*np.pi*time*freq[i])
        
    """
    
    smat = np.outer(time, freq)
    smat = np.sin(2*np.pi*smat)
        
    return smat  
      
def fit_sines(time, y, weights, freq):
    """ Fit sine waves to data.
    
    Args:
        time (float): Times at which the data were taken.
        y (float): The data.
        weights (float): Weights corresponding to the data.
        freq (float): Frequencies for the sine waves.
        
    Returns:
        pars (float): Best-fit amplitudes of the sine waves.
        
    """
    
    mat = sin_mat(time, freq)
    pars = np.linalg.lstsq(mat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    
    return pars
    
def evaluate_sines(time, pars, freq):
    """ Evaluate the sum of a number of sine waves.
    
    Args:
        time (float): Times at which to evaluate the sine waves.
        pars (float): Amplitudes of the sine waves.
        freq (float): Frequencies of the sine waves.
        
    Returns:
        fit (float): The sum of the sine waves evaluated at the given times.
        
    """
    
    mat = sin_mat(time, freq)
    fit = np.dot(mat, pars)
    
    return fit
      
def cos_mat(time, freq):
    """ Create a matrix of cosine waves.
    
    Args:
        time (float): Times at which to evaluate the cosine waves.
        freq (float): Frequencies for the diffrent cosine waves.
        
    Returns:
        cmat (float): A 2d-array with cmat[:,i] = np.cos(2*np.pi*time*freq[i])
        
    """
    
    cmat = np.outer(time, freq)
    cmat = np.cos(2*np.pi*cmat)
        
    return cmat
    
def fit_cosines(time, y, weights, freq):
    """ Fit cosine waves to data.
    
    Args:
        time (float): Times at which the data were taken.
        y (float): The data.
        weights (float): Weights corresponding to the data.
        freq (float): Frequencies for the cosine waves.
        
    Returns:
        pars (float): Best-fit amplitudes of the cosine waves.
        
    """
    
    mat = cos_mat(time, freq)
    pars = np.linalg.lstsq(mat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    
    return pars
    
def evaluate_cosines(time, pars, freq):
    """ Evaluate the sum of a number of cosine waves.
    
    Args:
        time (float): Times at which to evaluate the cosine waves.
        pars (float): Amplitudes of the cosine waves.
        freq (float): Frequencies of the cosine waves.
        
    Returns:
        fit (float): The sum of the cosine waves evaluated at the given times.
        
    """
    
    mat = cos_mat(time, freq)
    fit = np.dot(mat, pars)
    
    return fit
    
def fourier_mat(time, freq):
    """ Create a matrix of sine and cosine waves.
    
    Args:
        time (float): Times at which to evaluate the waves.
        freq (float): Frequencies for the diffrent waves.
        
    Returns:
        mat (float): A 2d-array with mat[:,i] = np.sin(2*np.pi*time*freq[i])
                     and mat[:,i+len(freq)] = np.cos(2*np.pi*time*freq[i])
        
    """
    
    smat = sin_mat(time, freq)
    cmat = cos_mat(time, freq)
    mat = np.hstack([smat, cmat])
    
    return mat
    
def fit_fourier(time, y, weights, freq):
    """ Fit sine and cosine waves to data.
    
    Args:
        time (float): Times at which the data were taken.
        y (float): The data.
        weights (float): Weights corresponding to the data.
        freq (float): Frequencies for the waves.
        
    Returns:
        pars (float): Best-fit amplitudes of the waves.
        
    """
    
    mat = fourier_mat(time, freq)
    pars = np.linalg.lstsq(mat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    
    return pars
    
def evaluate_fourier(time, pars, freq):
    """ Evaluate the sum of a number of sine and cosine waves.
    
    Args:
        time (float): Times at which to evaluate the waves.
        pars (float): Amplitudes of the waves.
        freq (float): Frequencies of the waves.
        
    Returns:
        fit (float): The sum of the waves evaluated at the given times.
        
    """
    
    mat = fourier_mat(time, freq)
    fit = np.dot(mat, pars)

    return fit

def fit_mat(y, weights, mat):
    """ Fit an arbitrary matrix to data.
    
    Args:
        y (float): The data.
        weights (float): Weights corresponding to the data.
        mat (float): The matrix of basis vectors to fit to the data.
        
    Returns:
        pars (float): The best fit amplitudes of the basis vectors.
    
    """
    
    pars = np.linalg.lstsq(mat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    
    return pars
    

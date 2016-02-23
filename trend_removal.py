#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

import filters
import fourierfuncs

from boxlstsq import boxlstsq

def freqs(P, n):
    """ Compute the harmonics of a period.
    
    P: The base period.
    n: The number of harmonics.
    
    """
    
    freq = (np.arange(n) + 1)/P
    
    return freq
    
def freqs_fft(N, step):
    """ Compute the frequencies in time series.
    
    N: number of points in the series.
    step: Time interval between the points.
    
    """ 
    
    freq = np.fft.rfftfreq(N, step)
    
    return freq
    
def fourier_fit(t, y, yerr, freq):
    
    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)
        
    sin_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        sin_mat[:,i] = np.sin(2*np.pi*freq[i]*t)        
    
    mat = np.hstack([sin_mat, cos_mat])
    pars = np.linalg.lstsq(mat/yerr[:,None], y/yerr)[0]
    
    fit = np.dot(mat, pars)
    
    return pars, fit

def evaluate(t, pars, freq):
    
    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)
        
    sin_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        sin_mat[:,i] = np.sin(2*np.pi*freq[i]*t)
    
    mat = np.hstack([sin_mat, cos_mat])
    fit = np.dot(mat, pars)
    
    return fit

def remove_trend0(jdmid, lst, mag, emag):

    freq_jd = freqs(180., 20)
    freq_lst = freqs(24., 5)
    
    fit1 = 0
    fit2 = 0
    for i in range(5):
        
        pars1, fit1 = fourier_fit(jdmid, mag - fit2, emag, freq_jd)
        pars2, fit2 = fourier_fit(lst, mag - fit1, emag, freq_lst) 

    return pars1, pars2, freq_jd, freq_lst, fit1, fit2

def remove_trend1(lstseq, jdmid, lst, mag, emag):
    
    N = np.ptp(lstseq) + 1 
    step = np.amin(np.diff(jdmid))
    freq_jd = freqs_fft(N, step)
    freq_jd = freq_jd[(freq_jd > 0) & (freq_jd < 1./3)]
    
    N = np.ptp(lstseq%270) + 1
    step = 320./3600.
    freq_lst = freqs_fft(N, step)
    freq_lst = freq_lst[(freq_lst > 0) & (freq_lst < 1./.5)]
    
    fit1 = 0
    fit2 = 0
    for i in range(5):
        
        pars1, fit1 = fourier_fit(jdmid, mag - fit2, emag, freq_jd)
        pars2, fit2 = fourier_fit(lst, mag - fit1, emag, freq_lst) 
        
    return pars1, pars2, freq_jd, freq_lst, fit1, fit2

def cosine_fit(t, y, yerr, P, n):

    freq = freqs(P, n)

    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)

    #cos_mat = np.hstack([np.ones((len(t),1)), cos_mat])

    pars = np.linalg.lstsq(cos_mat/yerr[:,None], y/yerr)[0]
    fit = np.dot(cos_mat, pars)
    
    return pars, fit

def cosine_eval(t, pars, P, n):
    
    freq = freqs(P, n)

    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)
    
    fit = np.dot(cos_mat, pars)
    
    return fit

def main():
    
    with h5py.File('/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5', 'r') as f:
        grp = f['data/807144']
        lstseq = grp['lstseq'].value
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
    
    emag = emag/np.sqrt(nobs)
        
    select = (nobs == 50)
    lstseq = lstseq[select]
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    
    x0 = np.floor(np.amin(jdmid))
    x1 = np.ceil(np.amax(jdmid)) - x0
    jdmid = jdmid - x0
    new_x1 = np.append(jdmid, 2*x1 - jdmid[::-1])
    new_x2 = np.append(lst, lst[::-1])
    new_y = np.append(mag, mag[::-1])
    new_yerr = np.append(emag, emag[::-1])
    
    N = np.ptp(lstseq%270) + 1
    freq = freqs_fft(N, 320./3600)
    
    fit2 = 0
    for i in range(5):
        pars1, fit1 = cosine_fit(new_x1, new_y - fit2, new_yerr, 2*x1, 40)
        pars2, fit2 = fourier_fit(new_x2, new_y - fit1, new_yerr, freq)
        
    print len(new_x1), len(pars1), len(pars2)
        
    chisq = (new_y - fit1 - fit2)**2/new_yerr**2
    print np.sum(chisq)/(len(new_x1) - len(pars1) - len(pars2))
        
    time1 = np.linspace(0, 2*x1, 1000)
    model1 = cosine_eval(time1, pars1, 2*x1, 40)
    
    #time2 = np.linspace(0, 24., 1000)
    #model2 = cosine_eval(time2, pars2, 24., 10)
        
    ax1 = plt.subplot(221)
    plt.plot(new_x1, new_y - fit2, '.')
    plt.plot(time1, model1)
    
    ax2 = plt.subplot(222)
    plt.plot(new_x2, new_y - fit1, '.')
    #plt.plot(time2, model2)
    
    plt.subplot(223, sharex=ax1, sharey=ax1)
    plt.plot(new_x1, new_y - fit1 - fit2, '.')
    
    plt.xlabel('Time [JD]')
    
    plt.subplot(224, sharex=ax2, sharey=ax2)
    plt.plot(new_x2, new_y - fit1 - fit2, '.')
    
    plt.xlabel('Time [LST]')
    
    plt.tight_layout()
    plt.show()
    
    exit()
    #chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 20)
    #pars1, pars2, freq1, freq2, fit1, fit2 = remove_trend0(jdmid, lst, mag, emag)
    pars1, pars2, freq1, freq2, fit1, fit2 = remove_trend1(lstseq, jdmid, lst, mag, emag)
    
    print len(pars1) + len(pars2)
    chisq = (mag - fit1 - fit2)**2/emag**2
    print np.sum(chisq)/(len(jdmid) - len(pars1) - len(pars2))
    
    time1 = np.linspace(np.amin(jdmid), np.amax(jdmid), 1000)
    model1 = evaluate(time1, pars1, freq1)
    
    time2 = np.linspace(np.amin(lst), np.amax(lst), 1000)
    model2 = evaluate(time2, pars2, freq2)
    
    plt.figure(figsize=(16,9))
    
    ax1 = plt.subplot(221)
    plt.errorbar(jdmid, mag - fit2, yerr=emag, fmt='.', c='k')
    plt.plot(time1, model1, c='r', lw=2)
    
    ax2 = plt.subplot(222)
    plt.errorbar(lst, mag - fit1, yerr=emag, fmt='.', c='k')
    plt.plot(time2, model2, c='r', lw=2)
    
    plt.subplot(223, sharex=ax1, sharey=ax1)
    plt.errorbar(jdmid, mag - fit1 - fit2, yerr=emag, fmt='.', c='k')
    
    plt.xlabel('Time [JD]')
    
    plt.subplot(224, sharex=ax2, sharey=ax2)
    plt.errorbar(lst, mag - fit1 - fit2, yerr=emag, fmt='.', c='k')
    
    plt.xlabel('Time [LST]')
    
    plt.tight_layout()
    plt.show()
    
    
    
    
    
    
    
    
    
        
        
        
    
    return

if __name__ == '__main__':
    main()
    

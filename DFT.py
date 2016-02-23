#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

import filters

from boxlstsq import boxlstsq

def weighed_fourier_fit(t, y, yerr, freq):
    
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

def fourier_fit(t, y, freq):
    
    cos_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        cos_mat[:,i] = np.cos(2*np.pi*freq[i]*t)
        
    sin_mat = np.zeros((len(t), len(freq)))
    for i in range(len(freq)):
        sin_mat[:,i] = np.sin(2*np.pi*freq[i]*t)        
    
    mat = np.hstack([sin_mat, cos_mat])
    pars = np.linalg.lstsq(mat, y)[0]
    
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

def fft_vs_lstsq():
    
    npoints = 500
    BL = 5.
    P = .33
    
    t = np.linspace(0, BL, npoints)
    y = .1*np.sin(2*np.pi*t/P)
    
    freq = np.fft.rfftfreq(npoints, BL/(npoints-1))
    Y = np.fft.rfft(y)
    Y = np.abs(Y)
    Y = 2*Y/npoints # Factor two because real fft.
    
    pars, fit = fourier_fit(t, y, freq)
    P = np.sqrt(pars[:len(pars)/2]**2 + pars[len(pars)/2:]**2)
    
    print np.amin(1/freq[freq>0])
    
    
    x = np.linspace(0, 2*BL, 2*npoints)
    fit1 = evaluate(x, pars, freq)
    
    plt.subplot(311)
    plt.plot(t, y)
    plt.plot(t, fit)
    plt.plot(x, fit1)
    
    plt.subplot(312)
    plt.plot(freq, Y)
    
    plt.subplot(313)
    plt.plot(freq, P)
    
    plt.show()
    
    return
    
def lstsq_on_data():
    
    filename = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    with h5py.File(filename, 'r') as f:
        grp = f['data/714995']
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        lstseq = grp['lstseq'].value
        
    select = (nobs == 50)
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    lstseq = lstseq[select]
    
    emag = emag/np.sqrt(50)
    lstidx = lstseq%270
    
    n = np.ptp(lstidx) + 1
    freq = np.fft.rfftfreq(n, 320./3600.)
    
    pars, fit = weighed_fourier_fit(lst, mag, emag, freq)
    
    plt.subplot(211)
    plt.errorbar(lst, mag, yerr=emag, fmt='.')
    plt.plot(lst, fit, '.')
    
    plt.subplot(212)
    plt.errorbar(lst, mag - fit, yerr=emag, fmt='.')
    
    plt.show()
    
    plt.subplot(211)
    plt.errorbar(jdmid, mag, yerr=emag, fmt='.')
    plt.plot(jdmid, fit)
    
    plt.subplot(212)
    plt.errorbar(jdmid, mag - fit, yerr=emag, fmt='.')
    
    plt.show()
    
    return
    
def lstsq_on_data2():
    
    filename = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    with h5py.File(filename, 'r') as f:
        grp = f['data/714995']
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        lstseq = grp['lstseq'].value
        
    select = (nobs == 50)
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    lstseq = lstseq[select]
    
    emag = emag/np.sqrt(50)
    
    n = np.ptp(lstseq) + 1
    step = np.amin(np.diff(jdmid))
    print step, 320./(24*3600.)

    freq = np.fft.rfftfreq(n, step)
    freq = freq[freq<1/3.]
    pars, fit = weighed_fourier_fit(jdmid, mag, emag, freq)
    
    print freq
    time = np.linspace(np.amin(jdmid), np.amax(jdmid), 5000)
    fit1 = evaluate(time, pars, freq)

    plt.subplot(211)
    plt.errorbar(lst, mag, yerr=emag, fmt='.')
    plt.plot(lst, fit, '.')
    
    plt.subplot(212)
    plt.errorbar(lst, mag - fit, yerr=emag, fmt='.')
    
    plt.show()
    
    plt.subplot(211)
    plt.errorbar(jdmid, mag, yerr=emag, fmt='.')
    plt.plot(jdmid, fit, '.')
    plt.plot(time, fit1)
    
    plt.subplot(212)
    plt.errorbar(jdmid, mag - fit, yerr=emag, fmt='.')
    
    plt.show()
    
    return
    
def lstsq_on_data3():
    
    filename = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    with h5py.File(filename, 'r') as f:
        grp = f['data/807144']
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        lstseq = grp['lstseq'].value
        
    select = (nobs == 50)
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    lstseq = lstseq[select]
    
    emag = emag/np.sqrt(50)
    
    fit1 = 0
    fit2 = 0
    for i in range(5):
    
        n = np.ptp(lstseq) + 1
        step = np.amin(np.diff(jdmid))
        freq = np.fft.rfftfreq(n, step)
        freq = freq[freq<1.5]
        pars, fit1 = weighed_fourier_fit(jdmid, mag - fit2, emag, freq)
        
        fit1 = evaluate(jdmid, pars, freq)
        
        plt.subplot(211)
        plt.plot(jdmid, mag - fit2, '.')
        plt.plot(jdmid, fit1)
        
        A = np.sqrt(pars[:len(freq)]**2 + pars[len(freq):]**2)
        
        plt.subplot(212)
        plt.plot(freq, A)
        plt.ylabel('Amplitude')
        plt.xlabel(r'Frequency [day$^{-1}$]')
        
        plt.show()
        
        n = np.ptp(lstseq%270) + 1
        freq = np.fft.rfftfreq(n, 320./3600.)
        freq = freq[1:]
        freq = freq[freq<1./1]
        pars, fit2 = weighed_fourier_fit(lst, mag-fit1, emag, freq)

        plt.subplot(211)
        plt.plot(jdmid, mag - fit1, '.')
        plt.plot(jdmid, fit2)
        
        A = np.sqrt(pars[:len(freq)]**2 + pars[len(freq):]**2)
        
        plt.subplot(212)
        plt.plot(freq, A)
        plt.ylabel('Amplitude')
        plt.xlabel(r'Frequency [hour$^{-1}$]')
        
        plt.show()
        

        fit = fit1 + fit2

        #plt.subplot(211)
        #plt.errorbar(lst, mag, yerr=emag, fmt='.')
        #plt.plot(lst, fit, '.')
        
        #plt.subplot(212)
        #plt.errorbar(lst, mag - fit, yerr=emag, fmt='.')
        
        #plt.show()
        
        #plt.subplot(211)
        #plt.errorbar(jdmid, mag, yerr=emag, fmt='.')
        #plt.plot(jdmid, fit, '.')
        
        #plt.subplot(212)
        #plt.errorbar(jdmid, mag - fit, yerr=emag, fmt='.')
        
        #plt.show()
    
    mag = mag - fit
    
    freq, dchisq, depth, hchisq, chisq0, epoch, duration = boxlstsq(jdmid, mag, 1/emag**2)
    
    plt.plot(freq, dchisq)
    plt.xlim(0, 2)
    plt.ylim(0, 2000)
    plt.show()
    
    return
    
def filter_on_data():
    
    filename = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    with h5py.File(filename, 'r') as f:
        grp = f['data/807144']
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        lstseq = grp['lstseq'].value
        
    select = (nobs == 50)
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    lstseq = lstseq[select]
    
    emag = emag/np.sqrt(50)
    
    chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 20)
    
    mag = mag - fit
    
    freq, dchisq, depth, hchisq, chisq0, epoch, duration = boxlstsq(jdmid, mag, 1/emag**2)
    
    plt.plot(freq, dchisq)
    plt.xlim(0, 2)
    plt.ylim(0, 2000)
    plt.show()
    
    return
    
def mean_on_data():
    
    filename = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    with h5py.File(filename, 'r') as f:
        grp = f['data/807144']
        jdmid = grp['jdmid'].value
        lst = grp['lst'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        lstseq = grp['lstseq'].value
        
    select = (nobs == 50)
    jdmid = jdmid[select]
    lst = lst[select]
    mag = mag[select]
    emag = emag[select]
    lstseq = lstseq[select]
    
    emag = emag/np.sqrt(50)
    lstidx = lstseq%270
    
    fit1 = 0
    fit2 = 0
    
    for i in range(5):
        fit1 = np.bincount(lstidx, (mag - fit2)/emag**2)/np.bincount(lstidx, 1/emag**2)
        fit1 = fit1[lstidx]
        
        plt.subplot(211)
        plt.errorbar(jdmid, mag - fit2, yerr=emag, fmt='.')
        plt.plot(jdmid, fit1, '.')
        
        plt.subplot(212)
        plt.errorbar(jdmid, mag - fit1 - fit2, yerr=emag, fmt='.')
        
        plt.show()
        
        n = np.ptp(lstseq) + 1
        step = np.amin(np.diff(jdmid))
        freq = np.fft.rfftfreq(n, step)
        freq = freq[freq<1/9.]
        pars, fit2 = weighed_fourier_fit(jdmid, mag - fit1, emag, freq)
        
    mag = mag - fit1 - fit2
    
    freq, dchisq, depth, hchisq, chisq0, epoch, duration = boxlstsq(jdmid, mag, 1/emag**2)
    
    plt.plot(freq, dchisq)
    plt.xlim(0, 2)
    plt.ylim(0, 2000)
    plt.show()
    
    
def main():
    
    #filter_on_data()
    mean_on_data()
    #lstsq_on_data3()
    
    return
    
if __name__ == '__main__':
    main()
    

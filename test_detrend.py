#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

# For A4 paper.
rcParams['figure.figsize'] = (11.69,8.27)
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

import boxlstsq
import fourier
import detrend

def read_data(filename, ascc):
    
    with h5py.File(filename, 'r') as f:
        
        try:
            grp = f['data/' + ascc]
        except:
            print 'Star not in file.'
            exit()
            
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
    
    return jdmid, lst, mag, emag, lstseq

def fixed_frequencies():
    
    filename = '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5'
    
    jdmid, lst, mag, emag, lstseq = read_data(filename, '807144')
    
    print np.ptp(jdmid), np.ptp(lst)
    
    # The current model.
    pars1, fit1, chisq1 = detrend.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 21, 24., 6)
    print chisq1, pars1
    
    freq = fourier.frequencies(180., 21)
    lt_model1 = fourier.evaluate_fourier(jdmid, pars1[:40], freq)
    
    freq = fourier.frequencies(24., 6)
    st_model1 = fourier.evaluate_fourier(lst, pars1[40:], freq)
    
    # The new model.
    pars2, fit2, chisq2 = detrend.masc_harmonic(jdmid, lst, mag, 1/emag**2, 1.2*np.ptp(jdmid), 20, 2*np.ptp(lst), 5)
    print chisq2, pars2
    
    freq = fourier.frequencies(1.2*np.ptp(jdmid), 20)
    lt_model2 = fourier.evaluate_fourier(jdmid, pars2[:38], freq)
    
    freq = fourier.frequencies(2*np.ptp(lst), 5)
    st_model2 = fourier.evaluate_fourier(lst, pars2[38:], freq)
    arg = np.argsort(lst)
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag - st_model1, emag, fmt='.', c='k')
    plt.plot(jdmid, lt_model1, c='r')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212)
    ax2.invert_yaxis()
    plt.errorbar(lst, mag - lt_model1, emag, fmt='.', c='k')
    plt.plot(lst[arg], st_model1[arg], c='r')
    plt.xlabel('Time [LST]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag - st_model2, emag, fmt='.', c='k')
    plt.plot(jdmid, lt_model2, c='r')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212)
    ax2.invert_yaxis()
    plt.errorbar(lst, mag - lt_model2, emag, fmt='.', c='k')
    plt.plot(lst[arg], st_model2[arg], c='r')
    plt.xlabel('Time [LST]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag, emag, fmt='.', c='k')
    plt.plot(jdmid, fit1, c='r')
    plt.plot(jdmid, fit2, c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    ax2.invert_yaxis()
    plt.errorbar(jdmid, mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(jdmid, mag - fit2, emag, fmt='.', c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    return

def free_frequencies():
    
    filename = '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5'
    
    jdmid, lst, mag, emag, lstseq = read_data(filename, '714995')
    #jdmid = jdmid - np.amin(jdmid)
    
    step = np.amin(np.diff(jdmid))
    ns = np.ptp(lstseq) + 1
    
    # The current model.
    pars1, fit1, chisq1 = detrend.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 21, 24., 6)
    print chisq1, pars1
    
    pars2, fit2, chisq2 = detrend.iterative_detrending(jdmid, mag, 1/emag**2, step, ns, maxiter=17)
    print chisq2, pars2
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag, emag, fmt='.', c='k')
    plt.plot(jdmid, fit1, c='r')
    plt.plot(jdmid, fit2, c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    ax2.invert_yaxis()
    plt.errorbar(jdmid, mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(jdmid, mag - fit2, emag, fmt='.', c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    return
    
def hybrid():
    
    filename = '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5'
    
    with h5py.File(filename, 'r') as f:
        ascc = f['/header/ascc'].value
    
    #for i in range(0, len(ascc), 50):
        
    jdmid, lst, mag, emag, lstseq = read_data(filename, '344250')
    
    if len(jdmid) == 0: exit()
    
    jdmid = jdmid - np.amin(jdmid)
    
    if np.ptp(jdmid) < 60.: exit()
    
    step = np.amin(np.diff(jdmid))
    ns = np.ptp(lstseq) + 1
    
    # The current model.
    pars1, fit1, chisq1 = detrend.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 21, 24., 6)
    print chisq1, len(pars1)
    
    freq = fourier.frequencies(180., 21)
    lt_model1 = fourier.evaluate_fourier(jdmid, pars1[:40], freq)
    
    freq = fourier.frequencies(24., 6)
    st_model1 = fourier.evaluate_fourier(lst, pars1[40:], freq)
    
    pars2, lt_model2, st_model2, chisq2 = detrend.hybrid(jdmid, lstseq%270, mag, 1/emag**2, step, ns)
    print chisq2, len(pars2)
    print pars2[:10]
    
    fit2 = lt_model2 + st_model2
    
    arg = np.argsort(lst)
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    #plt.title(ascc[i])
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag - st_model1, emag, fmt='.', c='k')
    plt.plot(jdmid, lt_model1, c='r')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212)
    ax2.invert_yaxis()
    plt.errorbar(lst, mag - lt_model1, emag, fmt='.', c='k')
    plt.plot(lst[arg], st_model1[arg], c='r')
    plt.xlabel('Time [LST]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    #plt.title(ascc[i])
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag - st_model2, emag, fmt='.', c='k')
    plt.plot(jdmid, lt_model2, c='r')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212)
    ax2.invert_yaxis()
    plt.errorbar(lst, mag - lt_model2, emag, fmt='.', c='k')
    plt.plot(lst[arg], st_model2[arg], c='r')
    plt.xlabel('Time [LST]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()

    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    #plt.title(ascc[i])
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag, emag, fmt='.', c='k')
    plt.plot(jdmid, fit1, c='r')
    plt.plot(jdmid, fit2, c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    #ax2.invert_yaxis()
    plt.errorbar(jdmid, mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(jdmid, mag - fit2, emag, fmt='.', c='g')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    #freq1, chisq0, dchisq1 = boxlstsq.boxlstsq(jdmid, mag - fit1, 1/emag**2)[:3]
    #freq2, chisq0, dchisq2 = boxlstsq.boxlstsq(jdmid, mag - fit2, 1/emag**2)[:3]
    
    #phase = np.mod(jdmid/2.21857, 1)
    
    #plt.subplot(211)
    #plt.plot(freq1, dchisq1)
    
    #plt.subplot(212)
    #plt.errorbar(phase, mag - fit1, emag, fmt='.')
    #plt.ylim(.1, -.1)
    
    #plt.show()
    
    #plt.subplot(211)
    #plt.plot(freq1, dchisq2)
    
    #plt.subplot(212)
    #plt.errorbar(phase, mag - fit2, emag, fmt='.')
    #plt.ylim(.1, -.1)
    
    #plt.show()
    
    return

def main():
    
    filename = '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5'
        
    jdmid, lst, mag, emag, lstseq = read_data(filename, '807144')
    jdmid = jdmid - np.amin(jdmid)
    
    step = [np.amin(np.diff(jdmid)), 320./3600.]
    ns = [np.ptp(lstseq) + 1, np.ptp(lstseq%270) + 1]
    
    # The current model.
    pars1, fit1, chisq1 = detrend.masc_harmonic(jdmid, lst, mag, 1/emag**2, 180., 21, 24., 6)
    print chisq1, len(pars1)
    
    pars2, fit2, chisq2 = detrend.new_harmonic(jdmid, lst, mag, 1/emag**2, step, ns, True)
    print chisq2, len(pars2)
    
    pars3, fit3, chisq3 = detrend.hybrid(jdmid, lstseq%270, mag, 1/emag**2, step[0], ns[0])
    print chisq3, len(pars3)
    
    arg = np.argsort(lst)

    fig = plt.figure()
    
    ax1 = plt.subplot(211)
    #plt.title(ascc[i])
    ax1.invert_yaxis()
    plt.errorbar(jdmid, mag, emag, fmt='.', c='k')
    plt.plot(jdmid, fit1, c='r')
    plt.plot(jdmid, fit2, c='g')
    plt.plot(jdmid, fit3, c='y')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    #ax2.invert_yaxis()
    plt.errorbar(jdmid, mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(jdmid, mag - fit2, emag, fmt='.', c='g')
    plt.errorbar(jdmid, mag - fit3, emag, fmt='.', c='y')
    plt.xlabel('Time [JD]')
    plt.ylabel(r'$\Delta m$')
    
    plt.show()
    
    freq1, chisq0, dchisq1 = boxlstsq.boxlstsq(jdmid, mag - fit1, 1/emag**2)[:3]
    freq2, chisq0, dchisq2 = boxlstsq.boxlstsq(jdmid, mag - fit2, 1/emag**2)[:3]
    
    phase = np.mod(jdmid/2.21857, 1)
    
    plt.subplot(211)
    plt.plot(freq1, dchisq1)
    
    plt.subplot(212)
    plt.errorbar(phase, mag - fit1, emag, fmt='.')
    plt.ylim(.1, -.1)
    
    plt.show()
    
    plt.subplot(211)
    plt.plot(freq1, dchisq2)
    
    plt.subplot(212)
    plt.errorbar(phase, mag - fit2, emag, fmt='.')
    plt.ylim(.1, -.1)
    
    plt.show()
    
    return

if __name__ == '__main__':
    main()

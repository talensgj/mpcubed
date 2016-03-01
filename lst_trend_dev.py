#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import fourier_dev

def masc_harmonic1(jdmid, lst, value, error, Pjd, njd, Plst, nlst, cnst=False):
    
    freq_jd = fourier_dev.frequencies(Pjd, njd, False)
    freq_lst = fourier_dev.frequencies(Plst, nlst, False)
    
    mat_jd = fourier_dev.fourier_mat(jdmid, freq_jd)
    mat_lst = fourier_dev.fourier_mat(lst, freq_lst)
    
    mat = np.hstack([mat_jd, mat_lst])
    if cnst: 
        mat = np.hstack([np.ones((len(value),1)), mat])
    
    pars = fourier_dev.fit_mat(value, error, mat)
    fit = np.dot(mat, pars)
    
    chisq = (value - fit)**2/error**2
    chisq = np.sum(chisq)
    
    return chisq, pars, fit, mat

def masc_harmonic2(jdmid, lst, value, error, stepjd, njd, steplst, nlst, cnst=False):
    
    freq_jd = fourier_dev.fftfreq(stepjd, njd, False)
    freq_lst = fourier_dev.fftfreq(steplst, nlst, False)

    freq_jd = freq_jd[freq_jd < 1/2.]
    freq_lst = freq_lst#[freq_lst < 1/.5]

    mat_jd = fourier_dev.fourier_mat(jdmid, freq_jd)
    mat_lst = fourier_dev.fourier_mat(lst, freq_lst)
    
    mat = np.hstack([mat_jd, mat_lst])
    if cnst: 
        mat = np.hstack([np.ones((len(value),1)), mat])
    
    pars = fourier_dev.fit_mat(value, error, mat)
    fit = np.dot(mat, pars)
    
    chisq = (value - fit)**2/error**2
    chisq = np.sum(chisq)
    
    return chisq, pars, fit, mat


def iterative_detrending(jdmid, mag, emag, step, ns, maxiter=10):
    
    from scipy.signal import lombscargle
    
    freqs = fourier_dev.fftfreq(step, ns)
    freqs = freqs[freqs <= 24./2.]
    normval = jdmid.shape[0]
    
    use_args = np.array([], dtype='int')
    fit = 0
    
    for niter in range(maxiter):
        
        # Compute the normalized Lomb-Scargle periodogram.
        pgram = lombscargle(jdmid, mag.astype('float64') - fit, 2*np.pi*freqs)
        pgram = np.sqrt(4.*pgram/normval)
        
        # Find the peak in the periogram.
        arg = np.argmax(pgram)
        
        # Add the corresponding frequency to the fit.
        if arg in use_args:
            print 'Frequency already in use.'
            break
        elif (pgram[arg] > .001):
            use_args = np.append(use_args, arg)
        else:
            print 'Small amplitude.'
            break
        
        # Fit the data.
        mat = fourier_dev.fourier_mat(jdmid, freqs[use_args])
        pars = fourier_dev.fit_mat(mag, emag, mat)
        fit = np.dot(mat, pars)
        
    chisq = (mag - fit)**2/emag**2
    chisq = np.sum(chisq)
    
    return chisq, pars, fit


def main():
    
    import matplotlib.pyplot as plt
    
    #jdmid = np.linspace(0, 10, 500)
    #mag = .4*np.sin(2*np.pi*3.*jdmid) + .23*np.cos(2*np.pi*.2*jdmid+.6) + .05*np.random.randn(500)
    #emag = .05*np.ones(500)
    
    #r = np.random.rand(500)
    #jdmid = jdmid[r < .9]
    #mag = mag[r < .9]
    #emag = emag[r < .9]
    
    #freqs, pars, fit = iterative_detrending(jdmid, mag, emag, 10./499, 500)
    
    #print freqs, pars
    
    #plt.subplot(211)
    #plt.errorbar(jdmid, mag, emag, fmt='.')
    #plt.plot(jdmid, fit)
    
    #plt.subplot(212)
    #plt.errorbar(jdmid, mag - fit, emag, fmt='.')
    
    #plt.tight_layout()
    #plt.show()
    
    import h5py
    
    with h5py.File('/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5', 'r') as f:
        
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
        
    step = np.amin(np.diff(jdmid))
    ns = np.ptp(lstseq) + 1
        
    chisq1, pars1, fit1, mat = masc_harmonic1(jdmid, lst, mag, emag, 180., 21, 24., 6)
    chisq2, pars2, fit2 = iterative_detrending(jdmid, mag, emag, step, ns, maxiter=17)
    
    print len(pars1), len(pars2)
    print chisq1, chisq2
    
    print pars1
    print pars2
    
    ax = plt.subplot(311)
    ax.invert_yaxis()
    plt.errorbar(jdmid, mag, emag, fmt='.', c='k')
    plt.plot(jdmid, fit1, c='r')
    plt.plot(jdmid, fit2, c='g')
    
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.errorbar(jdmid, mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(jdmid, mag - fit2, emag, fmt='.', c='g')
    
    plt.subplot(313, sharey=ax)
    plt.errorbar(np.mod(jdmid/2.21857, 1), mag - fit1, emag, fmt='.', c='r')
    plt.errorbar(np.mod(jdmid/2.21857, 1), mag - fit2, emag, fmt='.', c='g')
    
    plt.tight_layout()
    plt.show()
    
    return

if __name__ == '__main__':
    main()

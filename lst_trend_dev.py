#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import fourier_dev

def masc_harmonic(jdmid, lst, value, weights, Pjd, njd, Plst, nlst, cnst=False):
    """ Compute the best fit model for the long-term and LST trends. """
    
    # Compute the frequencies and matrix for the JD times.
    freq_jd = fourier_dev.frequencies(Pjd, njd, False)
    mat_jd = fourier_dev.fourier_mat(jdmid, freq_jd)
    
    # Compute the frequencies and matrix for the LST times.
    freq_lst = fourier_dev.frequencies(Plst, nlst, False)
    mat_lst = fourier_dev.fourier_mat(lst, freq_lst)
    
    # Compute the full matrix of basis functions.
    if cnst: 
        ones = np.ones((len(value),1))
        mat = np.hstack([ones, mat_jd, mat_lst])
    else:
        mat = np.hstcak([mat_jd, mat_lst])
    
    # Calculate the best fit.
    pars = fourier_dev.fit_mat(value, weights, mat)
    fit = np.dot(mat, pars)
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit)**2
    chisq = np.sum(chisq)
    
    return pars, fit, chisq

#def masc_harmonic2(jdmid, lst, value, error, stepjd, njd, steplst, nlst, cnst=False):
    
    #freq_jd = fourier_dev.fftfreq(stepjd, njd, False)
    #freq_lst = fourier_dev.fftfreq(steplst, nlst, False)

    #freq_jd = freq_jd[freq_jd < 1/2.]
    #freq_lst = freq_lst#[freq_lst < 1/.5]

    #mat_jd = fourier_dev.fourier_mat(jdmid, freq_jd)
    #mat_lst = fourier_dev.fourier_mat(lst, freq_lst)
    
    #mat = np.hstack([mat_jd, mat_lst])
    #if cnst: 
        #mat = np.hstack([np.ones((len(value),1)), mat])
    
    #pars = fourier_dev.fit_mat(value, error, mat)
    #fit = np.dot(mat, pars)
    
    #chisq = (value - fit)**2/error**2
    #chisq = np.sum(chisq)
    
    #return chisq, pars, fit, mat

#def iterative_detrending(jdmid, mag, emag, step, ns, maxiter=10):
    
    #from scipy.signal import lombscargle
    
    #freqs = fourier_dev.fftfreq(step, ns)
    #freqs = freqs[freqs <= 24./2.]
    #normval = jdmid.shape[0]
    
    #use_args = np.array([], dtype='int')
    #fit = 0
    
    #for niter in range(maxiter):
        
        ## Compute the normalized Lomb-Scargle periodogram.
        #pgram = lombscargle(jdmid, mag.astype('float64') - fit, 2*np.pi*freqs)
        #pgram = np.sqrt(4.*pgram/normval)
        
        ## Find the peak in the periogram.
        #arg = np.argmax(pgram)
        
        ## Add the corresponding frequency to the fit.
        #if arg in use_args:
            #print 'Frequency already in use.'
            #break
        #elif (pgram[arg] > .001):
            #use_args = np.append(use_args, arg)
        #else:
            #print 'Small amplitude.'
            #break
        
        ## Fit the data.
        #mat = fourier_dev.fourier_mat(jdmid, freqs[use_args])
        #pars = fourier_dev.fit_mat(mag, emag, mat)
        #fit = np.dot(mat, pars)
        
    #chisq = (mag - fit)**2/emag**2
    #chisq = np.sum(chisq)
    
    #return chisq, pars, fit


def main():
    
    return

if __name__ == '__main__':
    main()

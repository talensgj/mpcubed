#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import fourier

def masc_harmonic(jdmid, lst, value, weights, Pjd, njd, Plst, nlst, cnst=False):
    """ Compute the best fit model for the long-term and LST trends. """
    
    # Compute the frequencies and matrix for the JD times.
    freq_jd = fourier.frequencies(Pjd, njd, False)
    mat_jd = fourier.fourier_mat(jdmid, freq_jd)
    
    # Compute the frequencies and matrix for the LST times.
    freq_lst = fourier.frequencies(Plst, nlst, False)
    mat_lst = fourier.fourier_mat(lst, freq_lst)
    
    # Compute the full matrix of basis functions.
    if cnst: 
        ones = np.ones((len(value),1))
        mat = np.hstack([ones, mat_jd, mat_lst])
    else:
        mat = np.hstack([mat_jd, mat_lst])
    
    # Calculate the best fit.
    pars = fourier.fit_mat(value, weights, mat)
    fit = np.dot(mat, pars)
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit)**2
    chisq = np.sum(chisq)
    
    return pars, fit, chisq
    
    
def hybrid(jdmid, lstidx, mag, weights, step, ns, maxiter=10):
    
    freqs = fourier.fftfreq(step, ns)
    freqs = freqs[freqs < 1/6.]
    freqs = np.append(np.amin(freqs)/2, freqs)
    mat = fourier.fourier_mat(jdmid, freqs)
    
    #x = jdmid[:,None]/np.amax(jdmid)
    #mat = np.hstack([x, x**2, mat])
    
    fit2 = 0
    for niter in range(maxiter):
    
        pars1 = fourier.fit_mat(mag - fit2, weights, mat)
        fit1 = np.dot(mat, pars1)
        
        x0 = np.bincount(lstidx, weights)
        x1 = np.bincount(lstidx, weights*(mag - fit1))

        pars2 = x1/x0
        fit2 = pars2[lstidx]
        
    chisq = weights*(mag - fit1 - fit2)**2.
    chisq = np.sum(chisq)
    
    pars = np.append(pars1, pars2[np.unique(lstidx)]) 
    
    return pars, fit1 + fit2, chisq

def new_harmonic(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.), cnst=False):
    
    # Compute the frequencies and matrix for the JD times.
    freq1 = fourier.fftfreq(step[0], ns[0], False)
    freq1 = np.append(np.amin(freq1)/2, freq1)
    freq1 = freq1[freq1 < 1/6.]
    if (len(freq1) != 0):
        mat1 = fourier.fourier_mat(jdmid, freq1)
    else:
        mat1 = np.array([[]]*len(jdmid))
        
    # Compute the frequencies and matrix for the LST times.
    freq2 = fourier.fftfreq(step[1], ns[1], False)
    freq2 = np.append(np.amin(freq2)/2, freq2)
    freq2 = freq2[freq2 < 1/.5]
    if (len(freq2) != 0):
        mat2 = fourier.fourier_mat(lst, freq2)
    else:
        mat2 = np.array([[]]*len(jdmid))
        
    # Compute the full matrix of basis functions.
    if cnst: 
        ones = np.ones((len(value),1))
        mat = np.hstack([ones, mat1, mat2])
    else:
        mat = np.hstack([mat1, mat2])
    
    if (mat.shape[1] == 0):
        return [], np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    # Calculate the best fit.
    pars = fourier.fit_mat(value, weights, mat)
    fit = np.dot(mat, pars)
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit)**2
    chisq = np.sum(chisq)
    
    return pars, fit, chisq

def new_harmonic2(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.), cnst=False):
    
    # Compute the frequencies and matrix for the JD times.
    freq1 = fourier.fftfreq(step[0], ns[0], False)
    freq1 = np.append(np.amin(freq1)/2, freq1)
    freq1 = freq1[freq1 < 1/6.]
    if (len(freq1) != 0):
        mat1 = fourier.fourier_mat(jdmid, freq1)
    else:
        mat1 = np.array([[]]*len(jdmid))
        
    # Compute the frequencies and matrix for the LST times.
    freq2 = fourier.fftfreq(step[1], ns[1], False)
    freq2 = np.append(np.amin(freq2)/2, freq2)
    freq2 = freq2[freq2 < 1/.5]
    if (len(freq2) != 0):
        mat2 = fourier.fourier_mat(lst, freq2)
    else:
        mat2 = np.array([[]]*len(jdmid))
        
    # Compute the full matrix of basis functions.
    if cnst: 
        ones = np.ones((len(value),1))
        mat = np.hstack([ones, mat1, mat2])
    else:
        mat = np.hstack([mat1, mat2])
    
    if (mat.shape[1] == 0):
        return [], np.zeros(len(jdmid)), np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    # Calculate the best fit.
    pars = fourier.fit_mat(value, weights, mat)
    fit = np.dot(mat, pars)
    
    n = 2*len(freq1)
    fit1 = np.dot(mat[:,:n], pars[:n])
    fit2 = np.dot(mat[:,n:], pars[n:])
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit)**2
    chisq = np.sum(chisq)
    
    return pars, fit1, fit2, chisq

def new_harmonic3(jdmid, lst, value, weights, ns, step=(.003693591 ,320./3600.), cnst=False):
    
    # Compute the frequencies and matrix for the JD times.
    freq1 = fourier.fftfreq(step[0], ns[0], False)
    freq1 = np.append(np.amin(freq1)/2, freq1)
    freq1 = freq1[freq1 < 1/6.]
    if (len(freq1) != 0):
        mat1 = fourier.fourier_mat(jdmid, freq1)
    else:
        mat1 = np.array([[]]*len(jdmid))
        
    # Compute the frequencies and matrix for the LST times.
    freq2 = fourier.fftfreq(step[1], ns[1], False)
    freq2 = np.append(np.amin(freq2)/2, freq2)
    freq2 = freq2[freq2 < 1/.5]
    if (len(freq2) != 0):
        mat2 = fourier.fourier_mat(lst, freq2)
    else:
        mat2 = np.array([[]]*len(jdmid))
        
    # Compute the full matrix of basis functions.
    if cnst: 
        ones = np.ones((len(value),1))
        mat = np.hstack([ones, mat1, mat2])
    else:
        mat = np.hstack([mat1, mat2])
    
    if (mat.shape[1] == 0):
        return [], np.zeros(len(jdmid)), np.zeros(len(jdmid)), np.sum(weights*value**2)
    
    # Calculate the best fit.
    pars = fourier.fit_mat(value, weights, mat)
    fit = np.dot(mat, pars)
    
    n = 2*len(freq1)
    fit1 = np.dot(mat[:,:n], pars[:n])
    fit2 = np.dot(mat[:,n:], pars[n:])
    
    # Calculate the chi-square value of the fit.
    chisq = weights*(value - fit)**2
    chisq = np.sum(chisq)
    
    return freq1, freq2, pars[:n], pars[n:], fit1, fit2, chisq

def iterative_detrending(jdmid, mag, weights, step, ns, maxiter=10):
    
    from scipy.signal import lombscargle
    
    freqs1 = fourier.fftfreq(step, ns)
    freqs2 = fourier.frequencies(.9972695664, 20)
    freqs = np.append(freqs1, freqs2)
    #freqs = fourier.fftfreq(step, ns)
    freqs = freqs[freqs < 12.]
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
            print 'Frequency already in use.', 1./freqs[arg]
            break
        elif (pgram[arg] > np.mean(pgram)):
            use_args = np.append(use_args, arg)
        else:
            print 'Small amplitude.'
            break
        
        # Fit the data.
        mat = fourier.fourier_mat(jdmid, freqs[use_args])
        pars = fourier.fit_mat(mag, weights, mat)
        fit = np.dot(mat, pars)
        
    chisq = weights*(mag - fit)**2
    chisq = np.sum(chisq)
    
    print 1/freqs[use_args]
    
    return pars, fit, chisq


def main():
    
    return

if __name__ == '__main__':
    main()

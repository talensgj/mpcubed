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
    
    return chisq, pars, fit

def masc_harminic2(jdmid, lst, value, error, stepjd, njd, steplst, nlst, cnst=False):
    
    freq_jd = fourier_dev.fftfreq(stepjd, njd, False)
    freq_lst = fourier_dev.fftfreq(steplst, nlst, False)

    freq_jd = freq_jd[freq_jd < 1/9.]
    freq_lst = freq_lst[freq_lst < 1/.5]

    mat_jd = fourier_dev.fourier_mat(jdmid, freq_jd)
    mat_lst = fourier_dev.fourier_mat(lst, freq_lst)
    
    mat = np.hstack([mat_jd, mat_lst])
    if cnst: 
        mat = np.hstack([np.ones((len(value),1)), mat])
    
    pars = fourier_dev.fit_mat(value, error, mat)
    fit = np.dot(mat, pars)
    
    chisq = (value - fit)**2/error**2
    chisq = np.sum(chisq)
    
    return chisq, pars, fit

def main():
    
    import h5py
    
    import matplotlib.pyplot as plt
    
    import trend_removal
    
    with h5py.File('/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5', 'r') as f:
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
    
    njd = np.ptp(lstseq) + 1
    stepjd = np.amin(np.diff(jdmid))
    nlst = np.ptp(lstseq%270) + 1
    steplst = 320./3600
    
    chisq1, pars1, fit1 = remove_trend1(jdmid, lst, mag, emag, stepjd, njd, steplst, nlst, True)
    print chisq1, len(pars1)
    print pars1
    
    pars2, fit2 = trend_removal.remove_trend1(lstseq, jdmid, lst, mag, emag)
    print len(pars2)
    print pars2
    
    #print np.allclose(pars1, pars2)
    
    plt.errorbar(jdmid, mag, emag, fmt='.')
    plt.plot(jdmid, fit1, '.')
    plt.plot(jdmid, fit2, '.')
    plt.show()
    
    
    return

if __name__ == '__main__':
    main()

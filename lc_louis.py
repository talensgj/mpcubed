#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

def read_star(data, ascc):
    
    import detrend
    
    jdmid = np.array([])
    lst = np.array([])
    mag = np.array([])
    emag = np.array([])
    nobs = np.array([])
    trend = np.array([])
    
    for filename in data:
        
        with h5py.File(filename, 'r') as f:
        
            try:
                grp = f['data/'+ascc]
            except:
                continue
            
            lstseq = grp['lstseq'].value
            jdmid_ = grp['jdmid'].value
            lst_ = grp['lst'].value
            mag_ = grp['mag0'].value
            emag_ = grp['emag0'].value
            nobs_ = grp['nobs'].value
        
        emag_ = emag_/np.sqrt(nobs_)
        
        sel = (nobs_ > 49)
        lstseq = lstseq[sel]
        jdmid_ = jdmid_[sel]
        lst_ = lst_[sel]
        mag_ = mag_[sel]
        emag_ = emag_[sel]
        nobs_ = nobs_[sel]
        
        if (len(jdmid_) == 0):
            continue
        
        n1 = np.ptp(lstseq) + 1
        n2 = np.ptp(lstseq%270) + 1
        
        pars, trend_, chisq = detrend.new_harmonic(jdmid_, lst_, mag_, 1/emag_**2, [n1, n2])
        
        jdmid = np.append(jdmid, jdmid_)
        lst = np.append(lst, lst_)
        mag = np.append(mag, mag_)
        emag = np.append(emag, emag_)
        nobs = np.append(nobs, nobs_)
        trend = np.append(trend, trend_)
        
    return jdmid, lst, mag, emag, nobs, trend

def main():
    
    import matplotlib.pyplot as plt
    
    data = glob.glob('/data3/talens/2015Q?/LPS/red0_vmag_2015Q?LPS.hdf5')
    data = np.sort(data)
    
    ascc = ['1696259', '1696851', '1697548', '1697570', '1697656', '1698533', '1699174', '1793269', '1793437', '1795064']
    
    for i in range(len(ascc)):
        jdmid, lst, mag, emag, nobs, trend = read_star(data, ascc[i])
    
        plt.plot(jdmid, mag, '.', c='k')
        plt.plot(jdmid, trend, '.', c='r')
        plt.show()
        plt.close()
        
        with h5py.File('/home/talens/Nova_2015LPS.hdf5') as f:
            grp = f.create_group('data/'+ascc[i])
            grp.create_dataset('jdmid', data=jdmid)
            grp.create_dataset('lst', data=lst)
            grp.create_dataset('mag', data=mag)
            grp.create_dataset('emag', data=emag)
            grp.create_dataset('nobs', data=nobs)
            grp.create_dataset('trend', data=trend)
    
    return

if __name__ == '__main__':
    main()

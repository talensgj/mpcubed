#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package.statistics import statistics
import transit_candidates as tc 


def inspect_candidates(data, filelist):

    ascc = np.array([])
    period = np.array([])
    flag = np.array([])
    depth = np.array([])
    
    for filename in filelist:
        
        f = tc.blsFile(filename)
        hdr = f.read_header(['ascc', 'period', 'flag', 'depth'])
    
        ascc = np.append(ascc, hdr['ascc'])
        period = np.append(period, hdr['period'])
        flag = np.append(flag, hdr['flag'])
        depth = np.append(depth, hdr['depth'])

    sel = (flag == 8)
    print np.sum(sel)
    
    tc.candidates(data, filelist, ascc0=ascc[sel])
    
    return
 
def running_sde(dchisq, window=101):
    
    from scipy.ndimage.filters import median_filter
    
    m0 = median_filter(dchisq, [window, 1])
    m1 = 1.4826*median_filter(np.abs(dchisq - m0), [window, 1])
        
    sde = (dchisq - m0)/m1
    
    return sde, m0, m1
    
def simple_sde(dchisq):
    
    m0 = np.mean(dchisq, axis=0)
    m1 = np.std(dchisq, axis=0)
    
    sde = (dchisq - m0)/m1
    
    return sde, m0, m1

def main():
    
    data = glob.glob('/data3/talens/2015Q?/LP?/red0_vmag_2015Q?LP?.hdf5')
    filelist = glob.glob('/data3/talens/boxlstsq/2015Q1234/bls/*')
    
    for filename in filelist:
        
        f = tc.blsFile(filename)
        
        hdr = f.read_header(['ascc'])
        ascc = hdr['ascc']
        
        if '500717' in ascc:
            
            bls = f.read_data(['freq', 'dchisq', 'nt'])
            freq = bls['freq']
            dchisq = bls['dchisq']
            nt = bls['nt']
        
            ssde, sm0, sm1 = simple_sde(dchisq)
            msde, mm0, mm1 = running_sde(dchisq)
            #sde = dchisq
            
            #arg1 = np.argmax(sde, axis=0)
            #arg2 = np.arange(sde.shape[1])
            
            #peaks = sde[arg1, arg2]
            
            for i in range(dchisq.shape[1]):
                
                if ascc[i] == '500717':
                
                    ax = plt.subplot(311)
                    plt.plot(freq, dchisq[:,i])
                    plt.plot(freq, mm0[:,i])
                    
                    plt.subplot(312, sharex=ax)
                    plt.plot(freq, dchisq[:,i]-mm0[:,i])
                    plt.plot(freq, mm1[:,i])
                    
                    plt.subplot(313, sharex=ax)
                    plt.plot(freq, ssde[:,i])
                    plt.plot(freq, msde[:,i])
                    
                    plt.show()
                    plt.close()
    
        #tc.candidates(data, [filename], outdir='/data3/talens/boxlstsq/2015Q1234/sde_candidates/', ascc0=ascc[peaks > 20])
    
    return

if __name__ == '__main__':
    main()

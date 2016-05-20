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
    
def running_median(x, hwindow):

    med = np.zeros(len(x))
    for i in range(hwindow, len(x)-hwindow):

        med[i] = np.median(x[i-hwindow:i+hwindow+1])
        
    return med
    
def running_mad(x, hwindow):
    
    mad = np.zeros(len(x))
    for i in range(hwindow, len(x)-hwindow):

        #mad[i] = statistics.mad(x[i-hwindow:i+hwindow+1])
        mad[i] = np.median(np.abs(np.diff(x[i-hwindow:i+hwindow+1])))
        
    return mad
    
def SDE(dchisq):
    
    m0 = np.mean(dchisq, axis=0)
    m1 = np.std(dchisq, axis=0)
    
    sde = (dchisq - m0)/m1
    
    #m0r = np.zeros(dchisq.shape)
    #m1r = np.zeros(dchisq.shape)
    
    #for i in range(dchisq.shape[1]):
        #m0r[:,i] = running_median(dchisq[:,i], 50)
        #m1r[:,i] = running_mad(dchisq[:,i] - m0r[:,i], 50)
    
    #sde = (dchisq - m0r)/m1r
    
    return sde

def main():
    
    data = glob.glob('/data3/talens/2015Q?/LP?/red0_vmag_2015Q?LP?.hdf5')
    filelist = glob.glob('/data3/talens/boxlstsq/2015Q1234/bls/*')
    
    for filename in filelist:
        
        f = tc.blsFile(filename)
        
        hdr = f.read_header(['ascc'])
        ascc = hdr['ascc']
        
        bls = f.read_data(['freq', 'dchisq', 'nt'])
        freq = bls['freq']
        dchisq = bls['dchisq']
        nt = bls['nt']
    
        sde = SDE(dchisq)
        #sde = dchisq
        
        arg1 = np.argmax(sde, axis=0)
        arg2 = np.arange(sde.shape[1])
        
        peaks = sde[arg1, arg2]
        
        for i in range(sde.shape[1]):
            
            if ascc[i] in ['1339109']:
            
                ax = plt.subplot(311)
                plt.plot(freq, dchisq[:,i])
                
                plt.subplot(312, sharex=ax)
                plt.plot(freq, sde[:,i])
                
                plt.subplot(313, sharex=ax)
                plt.plot(freq, nt[:,i])
                
                plt.show()
                plt.close()
    
        #tc.candidates(data, [filename], outdir='/data3/talens/boxlstsq/2015Q1234/sde_candidates/', ascc0=ascc[peaks > 20])
    
    return

if __name__ == '__main__':
    main()

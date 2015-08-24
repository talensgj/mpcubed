#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/Jul2015/sky_20150716LPC.hdf5', 'r') as f:
    
    plt.subplot(221)
    plt.hist(f['header/nstars'])
    plt.subplot(222)
    plt.hist(f['header/niter'])
    plt.subplot(223)
    plt.hist(f['header/chisq'])
    plt.subplot(224)
    plt.hist(f['header/npoints'].value-f['header/npars'].value)
    plt.show()
    
    bins = f['header/binnum'].value
    
    for ind in bins:
        
        try:
            sky = f['data/%i'%ind]
        except:
            continue
        
        count = sky['count'].value
        
        ax = plt.subplot(311)
        plt.plot(sky['lstidx'][count>5], sky['sky'][count>5], '.', c='k')
        plt.plot(sky['lstidx'][count<=5], sky['sky'][count<=5], '.', c='r')
        plt.ylabel('Sky', size='x-large')
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.subplot(312, sharex=ax)
        plt.plot(sky['lstidx'][count>5], sky['count'][count>5], '.', c='k')
        plt.plot(sky['lstidx'][count<=5], sky['count'][count<=5], '.', c='r')
        plt.ylabel('# points', size='x-large')
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.subplot(313, sharex=ax)
        plt.plot(sky['lstidx'][count>5], (sky['chisq_sky'].value/sky['count'].value)[count>5], '.', c='k')
        plt.plot(sky['lstidx'][count<=5], (sky['chisq_sky'].value/sky['count'].value)[count<=5], '.', c='r')
        plt.xlabel('Time [lstidx]', size='x-large')
        plt.ylabel(r'$\chi^2/$(# points)', size='x-large')
        plt.xticks(size='large')
        plt.yticks(size='large')
        
        plt.tight_layout()
        plt.show()
        plt.close()

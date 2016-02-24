#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit

from run_boxlstsq import data_mc

def main():
    
    blsfile = '/home/talens/MASCARA/bls_266.hdf5'
    with h5py.File(blsfile, 'r') as f:
        
        grp = f['header']
        ascc = grp['ascc'].value
        chisq0 = grp['chisq0'].value
        period = grp['period'].value
        best_depth = grp['depth'].value
        best_epoch = grp['epoch'].value
        best_duration = grp['duration'].value
        flag = grp['flag'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        depth = grp['depth'].value
        duration = grp['duration'].value
        epoch = grp['epoch'].value
        nt = grp['nt'].value
    
    data = ['/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPS/red0_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5']
    
    jdmid, lst, mag, emag, mask = data_mc(data, ascc)
    
    jdmid = jdmid - np.amin(jdmid)
    
    for i in range(len(ascc)):
        
        if (flag[i] > 0): continue 
        
        x = jdmid[~mask[i]]
        y = mag[i, ~mask[i]]
        yerr = emag[i, ~mask[i]]
        
        phase = (x - best_epoch[i])/period[i]
        phase = np.mod(phase + .5, 1.) - .5
        
        time = np.linspace(0, period[i], 9*period[i]/best_duration[i] + 1)
        mphase = (time - best_epoch[i])/period[i]
        mphase = np.mod(mphase + .5, 1.) - .5
        model = transit.box_model(time, period[i], best_epoch[i], -best_depth[i], best_duration[i])
        
        sort = np.argsort(mphase)
        mphase = mphase[sort]
        model = model[sort]
        
        fig  = plt.figure(figsize=(16,9))
        
        gs = gridspec.GridSpec(3, 1, height_ratios = [1,10,10])
        
        plt.suptitle('ASCC {}'.format(ascc[i]), size='xx-large')
        
        ax = plt.subplot(gs[1])
        ax.invert_yaxis()
        plt.errorbar(phase, y, yerr, fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.ylabel('Magnitude')
        
        pbls = [period[i], best_epoch[i], best_depth[i], best_duration[i], 10.]
        try:
            popt, pcov = curve_fit(transit.softened_box_model, x, y, pbls, yerr, absolute_sigma=True)
        except:
            popt = pbls
        
        phase = (x - popt[1])/popt[0]
        phase = np.mod(phase + .5, 1.) - .5
        
        time = np.linspace(0, popt[0], 9*popt[0]/popt[3] + 1)
        mphase = (time - popt[1])/popt[0]
        mphase = np.mod(mphase + .5, 1.) - .5
        model = transit.softened_box_model(time, *popt)
        sort = np.argsort(mphase)
        mphase = mphase[sort]
        model = model[sort]
        
        ax = plt.subplot(gs[2])
        ax.invert_yaxis()
        plt.errorbar(phase, y, yerr, fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        fig = plt.figure(figsize=(16, 12))
        
        gs = gridspec.GridSpec(4, 2, height_ratios = [1,10,10,10])
        
        plt.suptitle('ASCC {}'.format(ascc[i]), size='xx-large')
        
        ax = plt.subplot(gs[1,:])
        plt.annotate('$P = {:.2f}$ days'.format(period[i]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, dchisq[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        plt.subplot(gs[2,0], sharex=ax)
        plt.annotate(r'$\delta = {:.2f}$'.format(best_depth[i]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, depth[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.ylabel(r'Depth')
        
        plt.subplot(gs[2,1], sharex=ax)
        plt.annotate(r'$\eta = {:.2f}$ days'.format(best_duration[i]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, duration[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.ylabel('Duration [days]')
        
        plt.subplot(gs[3,0], sharex=ax)
        plt.annotate(r'$T_p = {:.2f}$ days'.format(best_epoch[i]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, epoch[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel('Epoch [days]')
        
        plt.subplot(gs[3,1], sharex=ax)
        plt.plot(freq, nt[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel('nt')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
    return

if __name__ == '__main__':
    main()

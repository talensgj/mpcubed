#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np 

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from boxlstsq_ms_dev import freq_max, phase_duration

def bls_headers(data):
    
    ascc = np.array([])
    flag = np.array([], dtype='int')
    period = np.array([])
    depth = np.array([])
    duration = np.array([])
    nt = np.array([])
    
    for filename in data:
        with h5py.File(filename, 'r') as f:
            grp = f['header']
            ascc_ = grp['ascc'].value
            flag_ = grp['flag'].value
            period_ = grp['period'].value
            depth_ = grp['depth'].value
            duration_ = grp['duration'].value
            nt_ = grp['nt'].value
            
            if '906991' in ascc_:
                print filename
            
            ascc = np.append(ascc, ascc_)
            flag = np.append(flag, flag_)
            period = np.append(period, period_)
            depth = np.append(depth, depth_)
            duration = np.append(duration, duration_)
            nt = np.append(nt, nt_)
            
    return ascc, flag, period, depth, duration, nt

def lookup():
    
    R = np.linspace(.1, 50, 51)
    M = np.linspace(.1, 100, 101)
    
    Rg = (R[:-1] + R[1:])/2
    Mg = (M[:-1] + M[1:])/2
    
    Rg, Mg = np.meshgrid(Rg, Mg)
    
    fmax = freq_max(Mg, Rg)
    
    plt.pcolormesh(R, M, fmax, norm=LogNorm())
    plt.xlabel(r'$R/R_\odot$')
    plt.ylabel(r'$M/M_\odot$')
    cb = plt.colorbar()
    cb.set_label(r'$f_{max}$ [day$^{-1}$]')
    plt.show()
    
    return

def diagnostics():
    
    data = glob.glob('/data2/talens/2015Q2_vmag/boxlstsq/bls0_*.hdf5')
    data = np.sort(data)
    
    ascc, flag, period, depth, duration, nt = bls_headers(data)

    select = (flag < 1)
    
    plt.figure(figsize=(16, 6))

    plt.plot(1/period, depth, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(1/period[select], depth[select], c='r', zorder=1)
    plt.xlim(0, 1.8)
    plt.ylim(-.1, .1)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel('Depth')
    
    plt.tight_layout()
    plt.show()
    
    q = phase_duration(1/period, 1., 1.)
    
    plt.figure(figsize=(16, 6))

    plt.plot(1/period, duration/(period*q), '.', c='k', alpha=.2, zorder=0)
    plt.scatter(1/period[select], duration[select]/(period[select]*q[select]), c='r', zorder=1)
    plt.xlim(0, 1.8)
    plt.ylim(0, 3.5)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel('Duration [q]')
    
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(16, 6))

    plt.plot(depth, duration/(period*q), '.', c='k', alpha=.2, zorder=0)
    plt.scatter(depth[select], duration[select]/(period[select]*q[select]), c='r', zorder=1)
    plt.xlim(-.1, .1)
    plt.ylim(0, 3.5)
    plt.xlabel(r'Depth')
    plt.ylabel('Duration [q]')
    
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(16, 6))
    
    plt.subplot(111, xscale='log')
    plt.plot(nt, depth, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(nt[select], depth[select], c='r', zorder=1)
    plt.ylim(-.1, .1)
    plt.xlabel(r'$n_t$')
    plt.ylabel('Depth')
    
    plt.tight_layout()
    plt.show()

    return 

def common_noise(filename):
    
    with h5py.File(filename, 'r') as f:
        grp = f['header']
        flag = grp['flag'].value
        period = grp['period'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        
    plt.pcolormesh(freq, np.arange(163), dchisq.T, vmin=0, vmax=1000)
    cb = plt.colorbar()
    cb.set_label(r'$\Delta\chi^2$')
    
    plt.scatter(1/period, np.arange(len(period)), c='w', s=50)
    plt.xlim(freq[0], freq[-1])
    plt.ylim(0, 162)

    plt.xlabel('Frequency [day$^{-1}$]')
    
    plt.show()
        
    return
    
def periodogram(filename, ascc):
    
    with h5py.File(filename, 'r') as f:
        grp = f['header']
        sID = grp['ascc'].value
        flag = grp['flag'].value
        period = grp['period'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        
    i, = np.where(sID == ascc)
    
    print flag[i]
    print period[i].dtype
    plt.figure(figsize=(16,9))
    
    plt.title('ASCC {}'.format(ascc))
    plt.annotate('$P = {:.2f}$ days'.format(float(period[i])), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
    plt.plot(freq, dchisq[:,i], c='k')
    plt.axvline(1/period[i], c='r', ls='--')
    plt.xlabel('Frequency [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
        
    plt.tight_layout()
    plt.show()
        
    return

def plot_candidates(filename, data, outdir):
    
    patch = filename[-8:-5]
    print patch
    
    from package.models import transit
    from run_boxlstsq import data_mc
    
    with h5py.File(filename, 'r') as f:
        
        grp = f['header']
        ascc = grp['ascc'].value
        chisq0 = grp['chisq0'].value
        period = grp['period'].value
        best_depth = grp['depth'].value
        best_epoch = grp['epoch'].value
        best_duration = grp['duration'].value
        best_nt = grp['nt'].value
        flag = grp['flag'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        depth = grp['depth'].value
        duration = grp['duration'].value
        epoch = grp['epoch'].value
        nt = grp['nt'].value
    
    select = (flag < 1)
    if np.sum(select) < 1: 
        return
    
    jdmid, lst, mag, emag, mask, trend = data_mc(data, ascc)
    jdmid = jdmid - np.amin(jdmid)
    
    for i in range(len(ascc)):
        
        if (flag[i] > 0): continue 
        
        x = jdmid[~mask[i]]
        y = mag[i, ~mask[i]]
        yerr = emag[i, ~mask[i]]
        fit = trend[i, ~mask[i]]
        
        fig  = plt.figure(figsize=(16,9))
        
        gs = gridspec.GridSpec(3, 1, height_ratios = [1,10,10])
        
        plt.suptitle('ASCC {}'.format(ascc[i]), size='xx-large')
        
        ax = plt.subplot(gs[1])
        ax.invert_yaxis()
        plt.errorbar(x, y + fit, yerr, fmt='.', c='k')
        plt.plot(x, fit, '.', c='r')
        plt.ylabel('Magnitude')
        
        ax = plt.subplot(gs[2], sharex=ax)
        ax.invert_yaxis()
        plt.errorbar(x, y, yerr, fmt='.', c='k')
        plt.ylim(.1, -.1)
        plt.xlabel('Time [MJD]')
        plt.ylabel('Magnitude')
        
        plt.tight_layout()
        if np.amax(dchisq[:,i] > 1000):
            plt.savefig(outdir + 'prime_candidate_ASCC{}_lst_patch{}.png'.format(ascc[i], patch))
        else:
            plt.savefig(outdir + 'candidate_ASCC{}_lst_patch{}.png'.format(ascc[i], patch))
        #plt.show()
        plt.close()
        
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
        
        pbls = [period[i], best_epoch[i], -best_depth[i], best_duration[i], 10.]
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
        if np.amax(dchisq[:,i] > 1000):
            plt.savefig(outdir + 'prime_candidate_ASCC{}_lc_patch{}.png'.format(ascc[i], patch))
        else:
            plt.savefig(outdir + 'candidate_ASCC{}_lc_patch{}.png'.format(ascc[i], patch))
        #plt.show()
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
        plt.annotate(r'$\delta = {:.1f}$ %'.format(best_depth[i]*100), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
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
        plt.annotate(r'$n_t = {:d}$'.format(int(best_nt[i])), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, nt[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel('nt')
        
        plt.tight_layout()
        if np.amax(dchisq[:,i] > 1000):
            plt.savefig(outdir + 'prime_candidate_ASCC{}_bls_patch{}.png'.format(ascc[i], patch))
        else:
            plt.savefig(outdir + 'candidate_ASCC{}_bls_patch{}.png'.format(ascc[i], patch))
        #plt.show()
        plt.close()
        
    return

def main():
    
    #diagnostics()
    #lookup()
    periodogram('/data2/talens/2015Q2_vmag/boxlstsq/bls0_vmag_2015Q2_patch269.hdf5', '906991')
    
    return

if __name__ == '__main__':
    main()
    

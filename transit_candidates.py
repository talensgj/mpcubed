#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

# For A4 paper.
rcParams['figure.figsize'] = (11.69,8.27)
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit

from boxlstsq import phase_duration
from transit_search import read_header, read_data

def boxlstsq_header(filelist):
    
    ascc = np.array([])
    flag = np.array([], dtype='int')
    period = np.array([])
    depth = np.array([])
    duration = np.array([])
    nt = np.array([])
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            grp = f['header']
            ascc_ = grp['ascc'].value
            flag_ = grp['flag'].value
            period_ = grp['period'].value
            depth_ = grp['depth'].value
            duration_ = grp['duration'].value
            nt_ = grp['nt'].value
            
            ascc = np.append(ascc, ascc_)
            flag = np.append(flag, flag_)
            period = np.append(period, period_)
            depth = np.append(depth, depth_)
            duration = np.append(duration, duration_)
            nt = np.append(nt, nt_)
            
    return ascc, flag, period, depth, duration, nt

def plot_diagnostics(filelist):
    
    ascc, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    q = phase_duration(1/period, 1., 1.)
    select = (flag < 1)
    
    ax = plt.subplot(311)
    plt.plot(1/period, depth, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(1/period[select], depth[select], c='r', zorder=1)
    plt.xlim(0, 1.8)
    plt.ylim(-.1, .1)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel('Depth')
    
    plt.subplot(312, sharex=ax)
    plt.plot(1/period, duration/(period*q), '.', c='k', alpha=.2, zorder=0)
    plt.scatter(1/period[select], duration[select]/(period[select]*q[select]), c='r', zorder=1)
    plt.xlim(0, 1.8)
    plt.ylim(0, 3.5)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel('Duration [q]')
    
    plt.subplot(313, sharex=ax, yscale='log')
    plt.plot(1/period, nt, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(1/period[select], nt[select], c='r', zorder=1)
    plt.xlim(0, 1.8)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel(r'$n_t$')
    
    plt.tight_layout()
    plt.show()
    
    ax = plt.subplot(221, yscale='log')
    plt.plot(duration/(period*q), nt, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(duration[select]/(period[select]*q[select]), nt[select], c='r', zorder=1)
    plt.xlim(0, 3.5)
    plt.xlabel('Duration [q]')
    plt.ylabel(r'$n_t$')
    
    plt.subplot(223)
    plt.plot(duration/(period*q), depth, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(duration[select]/(period[select]*q[select]), depth[select], c='r', zorder=1)
    plt.xlim(0, 3.5)
    plt.ylim(-.1, .1)
    plt.xlabel('Duration [q]')
    plt.ylabel('Depth')
    
    plt.subplot(224, xscale='log')
    plt.plot(nt, depth, '.', c='k', alpha=.2, zorder=0)
    plt.scatter(nt[select], depth[select], c='r', zorder=1)
    plt.ylim(-.1, .1)
    plt.ylabel('Depth')
    plt.xlabel(r'$n_t$')
    
    plt.tight_layout()
    plt.show()
    
    return 

def plot_candidates(filename, data, outdir):
    
    ASCC, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)
    
    # Read the box least-squares results.
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
    
    ntransit = best_nt*(319.1262613/(24*3600))/best_duration
    
    # Determine if there are any transit candidates.
    select = (flag < 1)
    if np.sum(select) < 1: 
        return
    
    # Read the data.
    jdmid, lst, mag, emag, mask, trend = read_data(data, ascc)
    jdmid = jdmid - np.amin(jdmid)
    
    for i in range(len(ascc)):
        
        if (flag[i] > 0): continue 
        
        phase = np.mod((jdmid - best_epoch[i])/period[i] + .5, 1.) - .5
        
        # The best fit box-model.
        model = transit.box_model(jdmid, period[i], best_epoch[i], -best_depth[i], best_duration[i])
        weights = np.where(mask, 0., 1./emag**2)
        mean = np.sum(weights[i]*(mag[i] - model))/np.sum(weights[i])
        model = model + mean
        
        sort = np.argsort(phase)
        mphase = phase[sort]
        model = model[sort]
        
        # Bin the phase-folded lightcurve.
        nbins = np.ceil(9*period[i]/best_duration[i])
        bins = np.linspace(-.5, .5, nbins + 1)
        
        x0, _ = np.histogram(phase, bins, weights=weights[i])
        x1, _ = np.histogram(phase, bins, weights=weights[i]*mag[i])
        
        bphase = (bins[:-1] + bins[1:])/2.
        bmag = x1/x0
        ebmag = 1/np.sqrt(x0)
        
        # Make the figure.
        fig = plt.figure()
        
        gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
    
        plt.suptitle('ASCC {}, {}, $V={:.1f}$\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $N_t = {:.1f}, flag={:d}$'.format(ascc[i], sptype[ASCC==ascc[i]][0], vmag[ASCC==ascc[i]][0], best_depth[i]*100, period[i], best_duration[i], ntransit[i], flag[i]), size='xx-large')
        
        ax = plt.subplot(gs[1,:])
        plt.plot(freq, dchisq[:,i], c='k')
        
        plt.axvline(1/period[i], c='g', ls='--')
        for n in range(2, 5):
            plt.axvline(n/period[i], c='g', ls='--')
            plt.axvline(1/(n*period[i]), c='g', ls='--')
        #for n in range(1, 5):
            #plt.axvline(n/(period[i]*(n+1.)), c='g', ls='--', ymax=.4)
            #plt.axvline((n+1.)/(n*period[i]), c='g', ls='--', ymax=.4)
        
        plt.axvline(1/.9972, c='r', ymax=.1, lw=2)
        for n in range(2, 5):
            plt.axvline(n/.9972, c='r', ymax=.1, lw=2)
            plt.axvline(1/(n*.9972), c='r', ymax=.1, lw=2)
        for n in range(1, 5):
            plt.axvline(n/(.9972*(n+1.)), c='r', ymax=.1, lw=2)
            plt.axvline((n+1.)/(n*.9972), c='r', ymax=.1, lw=2)
            
        plt.xlim(0, 1.8)
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        plt.subplot(gs[2,:])
        plt.errorbar(phase[~mask[i]], mag[i,~mask[i]], emag[i,~mask[i]], fmt='.', c='k')
        #plt.errorbar(bphase, bmag, ebmag, fmt='.', c='g')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        plt.tight_layout()
        if (np.amax(dchisq[:,i]) > 1000):
            plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}.png'.format(ascc[i])))
        else:
            plt.savefig(os.path.join(outdir, 'candidate_ASCC{}.png'.format(ascc[i])))
        #plt.show()
        plt.close()
        
    return

def plot_lightcurve(data, ascc, period):
    
    jdmid, lst, mag, emag, mask, trend = read_data(data, ascc)
    jdmid = jdmid - np.amin(jdmid)
    
    for i in range(len(ascc)):
        phase = np.mod(jdmid/period[i], 1.)
        
        ax = plt.subplot(111)
        ax.invert_yaxis()
        plt.title('ASCC {}'.format(ascc[i]))
        plt.errorbar(phase[~mask[i]], mag[i, ~mask[i]], emag[i, ~mask[i]], fmt='.', c='k')
        plt.xlim(0,1)
        #plt.show()
        plt.savefig('binary_ASCC{}.png'.format(ascc[i]))
        plt.close()
    
    return

def plot_periodogram(filelist, ascc):
    
    for filename in filelist:
        with h5py.File(filename, 'r') as f:
            grp = f['header']
            sID = grp['ascc'].value
            
            if ascc in sID:
                grp = f['data']
                freq = grp['freq'].value
                dchisq = grp['dchisq'].value
                
                arg, = np.where(sID == ascc)
                
                plt.plot(freq, dchisq[:,arg])
                plt.show()
                
                break
    
    return
                

def periodogram_animation(filelist, data, sID):
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            ascc = f['header/ascc'].value
        
        if sID in ascc:
            break
    
    with h5py.File(filename, 'r') as f:
        
        grp = f['header']
        ascc = grp['ascc'].value
        chisq0 = grp['chisq0'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        depth = grp['depth'].value
        duration = grp['duration'].value
        epoch = grp['epoch'].value
    
    # Read the data.
    jdmid, lst, mag, emag, mask, trend = read_data(data, ascc)
    jdmid = jdmid - np.amin(jdmid)
    
    arg, = np.where(ascc == sID)
    arg = arg[0]
    
    for i in range(len(freq)):
        
        print i
    
        phase = np.mod((jdmid - epoch[i,arg])*freq[i] + .5, 1.) - .5
    
        # The best fit box-model.
        model = transit.box_model(jdmid, 1/freq[i], epoch[i,arg], -depth[i,arg], duration[i,arg])
        weights = np.where(mask, 0., 1./emag**2)
        mean = np.sum(weights[arg]*(mag[arg] - model))/np.sum(weights[arg])
        model = model + mean
        
        sort = np.argsort(phase)
        mphase = phase[sort]
        model = model[sort]
        
        # Make the figure.
        fig = plt.figure()
        
        gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
        
        plt.suptitle('ASCC {}'.format(sID))
        
        ax1 = plt.subplot(gs[1,:])
        plt.annotate('$P = {:.2f}$ days'.format(1/freq[i]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
        plt.plot(freq, dchisq[:,arg], c='k')
        plt.axvline(freq[i], c='g', ls='--')
        plt.xlim(freq[i]-.1, freq[i]+.1)
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        ax2 = plt.subplot(gs[2,:])
        ax2.invert_yaxis()
        plt.errorbar(phase[~mask[arg]], mag[arg,~mask[arg]], emag[arg,~mask[arg]], fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        plt.tight_layout()
        plt.savefig('/data2/talens/boxlstsq_movie/mov_{:04d}.png'.format(i))
        #plt.show()
        plt.close()
    
    return


def main():
    
    data = ['/data3/talens/2015Q2/LPN/red0_vmag_2015Q2LPN.hdf5',
            '/data3/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data3/talens/2015Q2/LPS/red0_vmag_2015Q2LPS.hdf5',
            '/data3/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5',
            '/data3/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5']
    
    filelist = glob.glob('/data3/talens/2015Q2/boxlstsq/bls0_*.hdf5')
    
    for filename in filelist:
        plot_candidates(filename, data, '/data3/talens/2015Q2/boxlstsq/figures')
    
    #filelist = glob.glob('/data2/talens/2015Q2_vmag/boxlstsq/bls0_*.hdf5')
    #ascc1, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    #ascc1 = ascc1[flag == 0]
    
    #filelist = glob.glob('/data2/talens/2015Q2_vmag/boxlstsq_new/bls0_*.hdf5')
    #ascc2, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    #ascc2 = ascc2[flag == 0]
    
    #print 'In 1 but not in 2:', np.setdiff1d(ascc1, ascc2)
    #print 'In 2 but not in 1:', np.setdiff1d(ascc2, ascc1)
    
    #for ascc in np.setdiff1d(ascc1, ascc2):
        #os.system('shotwell /data2/talens/2015Q2_vmag/boxlstsq/figures/*{}*'.format(ascc))
    
    #for ascc in np.setdiff1d(ascc2, ascc1):
        #os.system('shotwell /data2/talens/2015Q2_vmag/boxlstsq_new/figures/*{}*'.format(ascc))
    
    return

if __name__ == '__main__':
    main()

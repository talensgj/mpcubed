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

from scipy import optimize

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

def plot_candidates(data, filelist, outdir):
    
    ascc, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    
    select = (flag == 0)
    ascc = ascc[select]
    
    #plot_periodogram(data, filelist, ascc, outdir)
    new_periodogram(data, filelist, ascc, outdir)
    #plot_lightcurve(data, filelist, ascc, outdir)
    
    return

def plot_periodogram(data, filelist, ascc, outdir=None):
    
    ASCC, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            sID = grp['ascc'].value
            period = grp['period'].value
            best_depth = grp['depth'].value
            best_epoch = grp['epoch'].value
            best_duration = grp['duration'].value
            best_nt = grp['nt'].value
            flag = grp['flag'].value
            
            grp = f['data']
            freq = grp['freq'].value
            dchisq = grp['dchisq'].value

        select = np.in1d(sID, ascc)
        args, = np.where(select)
        if (len(args) == 0):
            continue 
            
        ntransit = best_nt*(319.1262613/(24*3600))/best_duration
            
        # Read the data.
        jdmid, lst, mag, emag, mask, trend = read_data(data, sID)
        jdmid = jdmid - np.amin(jdmid)
            
        for i in args:
            
            if (best_duration[i]/period[i] > .15):
                continue
                
            if (np.amax(dchisq[:,i]) < 400):
                continue
            
            phase = np.mod((jdmid - best_epoch[i])/period[i] + .5, 1.) - .5
    
            # The best fit box-model.
            model = transit.box_model(jdmid, period[i], best_epoch[i], -best_depth[i], best_duration[i])
            weights = np.where(mask, 0., 1./emag**2)
            mean = np.sum(weights[i]*(mag[i] - model))/np.sum(weights[i])
            model = model + mean
            
            sort = np.argsort(phase)
            mphase = phase[sort]
            model = model[sort]
            
            # Plot the periodogram.
            fig = plt.figure()
    
            gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
        
            plt.suptitle('ASCC {}, {}, $V={:.1f}$\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $N_t = {:.1f}$, flag={:d}'.format(sID[i], sptype[ASCC==sID[i]][0], vmag[ASCC==sID[i]][0], best_depth[i]*100, period[i], best_duration[i], ntransit[i], flag[i]), size='xx-large')
            
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
            plt.plot(mphase, model, c='r', lw=2)
            plt.xlim(-.5, .5)
            plt.ylim(.1, -.1)
            plt.xlabel('Phase')
            plt.ylabel(r'$\Delta m$')
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            elif (np.amax(dchisq[:,i]) >= 900):
                plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}.png'.format(sID[i])))
            else:
                plt.savefig(os.path.join(outdir, 'candidate_ASCC{}.png'.format(sID[i])))
            
            plt.close()

    return


def minfunc(pars, jdmid, mag, weights):
    
    model = transit.softened_box_model(jdmid, *pars)
    chisq = weights*(mag - model)**2

    return np.sum(chisq)

def plot_periodogram_refined(data, filelist, ascc, outdir=None):
    
    ASCC, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            sID = grp['ascc'].value
            period = grp['period'].value
            best_depth = grp['depth'].value
            best_epoch = grp['epoch'].value
            best_duration = grp['duration'].value
            best_nt = grp['nt'].value
            flag = grp['flag'].value
            
            grp = f['data']
            freq = grp['freq'].value
            dchisq = grp['dchisq'].value

        select = np.in1d(sID, ascc)
        args, = np.where(select)
        if (len(args) == 0):
            continue 
            
        ntransit = best_nt*(319.1262613/(24*3600))/best_duration
            
        # Read the data.
        jdmid, lst, mag, emag, mask, trend = read_data(data, sID)
        jdmid = jdmid - np.amin(jdmid)
            
        weights = np.where(mask, 0, 1/emag**2)
            
        for i in args:
            
            if (best_duration[i]/period[i] > .15):
                continue
                
            if (np.amax(dchisq[:,i]) < 400):
                continue
            
            p_bls = [period[i], best_epoch[i], -best_depth[i], best_duration[i], 10.] 
            #res = optimize.minimize(minfunc, p_bls, args=(jdmid, mag[i], weights[i]), jac=False)
            #res = optimize.minimize(minfunc, p_bls, method='SLSQP', args=(jdmid, mag[i], weights[i]), bounds=[(period[i]-.1, period[i]+.1), (best_epoch[i]-.1, best_epoch[i]+.1), (.005, .03), (0, .2), (1., 50.)])
            #p_opt = res.x
            
            try:
                p_opt, pcov = optimize.curve_fit(transit.softened_box_model, jdmid[~mask[i]], mag[i, ~mask[i]], p_bls, sigma=emag[i,~mask[i]], absolute_sigma=True)
            except:
                p_opt = p_bls
            
            
            phase = np.mod((jdmid - p_opt[1])/p_opt[0] + .5, 1.) - .5
    
            # Bin the data.
            nbins = 9*np.ceil(p_opt[0]/p_opt[3])
            bins = np.linspace(-.5, .5, nbins)
            m0, bins = np.histogram(phase, bins=bins, weights=weights[i])
            m1, bins = np.histogram(phase, bins=bins, weights=weights[i]*mag[i])
            m2, bins = np.histogram(phase, bins=bins, weights=weights[i]*mag[i]**2)
            
            phase_bin = (bins[:-1] + bins[1:])/2.
            mag_bin = m1/m0
            emag_bin = np.sqrt(1./m0)
    
            # The best fit box-model.
            model = transit.softened_box_model(jdmid, *p_opt)
            
            sort = np.argsort(phase)
            mphase = phase[sort]
            mmodel = model[sort]
            
            # Plot the periodogram.
            fig = plt.figure()
    
            gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
        
            plt.suptitle('ASCC {}, {}, $V={:.1f}$\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $N_t = {:.1f}$, flag={:d}'.format(sID[i], sptype[ASCC==sID[i]][0], vmag[ASCC==sID[i]][0], best_depth[i]*100, period[i], best_duration[i], ntransit[i], flag[i]), size='xx-large')
            
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
            plt.plot(phase[~mask[i]], mag[i,~mask[i]], '.', alpha=.5, c='k')
            plt.errorbar(phase_bin, mag_bin, emag_bin, fmt='o', c='g')
            plt.plot(mphase, mmodel, c='r', lw=2)
            plt.xlim(-.5, .5)
            plt.ylim(.05, -.02)
            #plt.xlabel('Phase')
            plt.ylabel(r'$\Delta m$')
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            elif (np.amax(dchisq[:,i]) >= 900):
                plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}.png'.format(sID[i])))
            else:
                plt.savefig(os.path.join(outdir, 'candidate_ASCC{}.png'.format(sID[i])))
            
            plt.close()

    return
    
def add_periodogram(ax, freq, pgram):
    
    ax.plot(freq, pgram, c='k')
    
    ax.set_xlim(0, 1.8)
    ax.set_xlabel('Frequency [day$^{-1}$]')
    ax.set_ylabel(r'$\Delta\chi^2$')
    
    return
    
def add_harmonics1(ax, freq, **kwargs):
    
    plt.axvline(freq, **kwargs)
    for n in range(2, 5):
        plt.axvline(n*freq, **kwargs)
        plt.axvline(freq/n, **kwargs)
    
    return
    
def add_harmonics2(ax, freq, **kwargs):
    
    plt.axvline(freq, **kwargs)
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), **kwargs)
        plt.axvline((n+1)*freq/n, **kwargs)
    
    return
    
def add_data(ax, jdmid, mag, emag, mask, pars):
    
    phase = (jdmid - pars[1])/pars[0]
    phase = np.mod(phase + .5, 1) - .5
    
    plt.errorbar(phase[~mask], mag[~mask], emag[~mask], fmt='.', c='k')
    
    plt.xlim(-.5, .5)
    plt.ylim(.1, -.1)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')

    return
    
def add_bindata(ax1, ax2, jdmid, mag, emag, mask, pars):
    
    phase = (jdmid - pars[1])/pars[0]
    phase = np.mod(phase + .5, 1) - .5
    
    weights = np.where(mask, 0., 1/emag**2)
    
    npoints = 9*np.ceil(pars[0]/pars[3])
    bins = np.linspace(-.5, .5, npoints+1)
    m0, bins = np.histogram(phase, bins=bins, weights=weights)
    m1, bins = np.histogram(phase, bins=bins, weights=weights*mag)
    
    phase_bin = (bins[:-1] + bins[1:])/2
    mag_bin = m1/m0
    emag_bin = np.sqrt(1/m0)
    
    ax1.set_title(r'$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days'.format(-100*pars[2], pars[0], pars[3]))
    
    ax1.plot(phase[~mask], mag[~mask], '.', c='k', alpha=.5)
    ax1.errorbar(phase_bin, mag_bin, emag_bin, fmt='o', c='g')
    
    model = transit.softened_box_model(jdmid, *pars)
    model_bin = transit.softbox_p(phase_bin, *pars)
    
    ax2.plot(phase[~mask], (mag - model)[~mask], '.', c='k', alpha=.5)
    ax2.errorbar(phase_bin, mag_bin - model_bin, fmt='o', c='g')
    
    ax1.set_xlim(-.5, .5)
    ax1.set_ylim(.05, -.02)
    ax1.set_xlabel('Phase')
    ax1.set_ylabel(r'$\Delta m$')
    
    ax2.set_xlim(-.5, .5)
    ax2.set_ylim(.02, -.02)
    ax2.set_yticks([-.02, -.01, 0., .01, .02])
    ax2.set_xlabel('Phase')
    
    return
    
def add_boxmodel(ax, pars):
    
    npoints = 9*np.ceil(pars[0]/pars[3])
    x = np.linspace(0, pars[0], npoints)
    y = transit.box_model(x, *pars)
    
    x = (x - pars[1])/pars[0]
    x = np.mod(x + .5, 1) - .5
    
    sort = np.argsort(x)
    x = x[sort]
    y = y[sort]
    
    plt.plot(x, y, c='r', lw=2)
    
    return
    
def add_sboxmodel(ax, jdmid, mag, emag, mask, pars):
    
    pars.append(10.)
    try:
        pars, pcov = optimize.curve_fit(transit.softened_box_model, jdmid[~mask], mag[~mask], pars, sigma=emag[~mask], absolute_sigma=True)
    except:
        pass
    
    npoints = 9*np.ceil(pars[0]/pars[3])
    x = np.linspace(0, pars[0], npoints)
    y = transit.softened_box_model(x, *pars)
    
    x = (x - pars[1])/pars[0]
    x = np.mod(x + .5, 1) - .5
    
    sort = np.argsort(x)
    x = x[sort]
    y = y[sort]
    
    plt.plot(x, y, c='r', lw=2, zorder=10)
    
    return pars
    
def new_periodogram(data, filelist, ascc, outdir=None):
    
    ASCC, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            sID = grp['ascc'].value
            period = grp['period'].value
            best_depth = grp['depth'].value
            best_epoch = grp['epoch'].value
            best_duration = grp['duration'].value
            best_nt = grp['nt'].value
            flag = grp['flag'].value
            
            grp = f['data']
            freq = grp['freq'].value
            dchisq = grp['dchisq'].value

        select = np.in1d(sID, ascc)
        args, = np.where(select)
        if (len(args) == 0):
            continue 
    
        ntransit = best_nt*(319.1262613/(24*3600))/best_duration
                
        # Read the data.
        jdmid, lst, mag, emag, mask, trend = read_data(data, sID)
        jdmid = jdmid - np.amin(jdmid)
            
        weights = np.where(mask, 0, 1/emag**2)
        
        for i in args:
            
            if (best_duration[i]/period[i] > .15):
                print 'eta/P > .15'
                continue
            
            pars = [period[i], best_epoch[i], -best_depth[i], best_duration[i]]
            
            # Plot the periodogram.
            fig = plt.figure(figsize=(8.27, 11.69))

            gs = gridspec.GridSpec(5, 2, height_ratios = [.5,10,10,10,5])
        
            plt.suptitle('ASCC {}, {}, $V={:.1f}$\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $N_t = {:.1f}$, flag={:d}'.format(sID[i], sptype[ASCC==sID[i]][0], vmag[ASCC==sID[i]][0], best_depth[i]*100, period[i], best_duration[i], ntransit[i], flag[i]), size='xx-large')
            
            ax = plt.subplot(gs[1,:])
            add_periodogram(ax, freq, dchisq[:,i])
            add_harmonics1(ax, 1./period[i], c='g', ls='--')
            add_harmonics1(ax, 1./.9972, c='r', lw=2, ymax=.1)
            add_harmonics2(ax, 1./.9972, c='r', lw=2, ymax=.1)
            
            ax = plt.subplot(gs[2,:])
            add_data(ax, jdmid, mag[i], emag[i], mask[i], pars)
            add_boxmodel(ax, pars)
            
            ax1 = plt.subplot(gs[3,:])
            new_pars = add_sboxmodel(ax1, jdmid, mag[i], emag[i], mask[i], pars)
            ax2 = plt.subplot(gs[4,:])
            add_bindata(ax1, ax2, jdmid, mag[i], emag[i], mask[i], new_pars)
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            elif (np.amax(dchisq[:,i]) < 400):
                plt.savefig(os.path.join(outdir, 'maybe_candidate_ASCC{}.png'.format(sID[i])))
            elif (np.amax(dchisq[:,i]) >= 900):
                plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}.png'.format(sID[i])))
            else:
                plt.savefig(os.path.join(outdir, 'candidate_ASCC{}.png'.format(sID[i])))
            
            plt.close()

    return

def plot_lightcurve(data, filelist, ascc, outdir=None):
    
    ASCC, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header']
            sID = grp['ascc'].value
            period = grp['period'].value
            best_depth = grp['depth'].value
            best_epoch = grp['epoch'].value
            best_duration = grp['duration'].value
            best_nt = grp['nt'].value
            flag = grp['flag'].value
            
            grp = f['data']
            freq = grp['freq'].value
            dchisq = grp['dchisq'].value

        select = np.in1d(sID, ascc)
        args, = np.where(select)
        if (len(args) == 0):
            continue 
            
        ntransit = best_nt*(319.1262613/(24*3600))/best_duration
            
        # Read the data.
        jdmid, lst, mag, emag, mask, trend = read_data(data, sID)
        jdmid = jdmid - np.amin(jdmid)
            
        for i in args:
            
            if (best_duration[i]/period[i] > .15):
                continue
                
            if (np.amax(dchisq[:,i]) < 400):
                continue
            
            # Plot the lightcurve.
            fig = plt.figure()
    
            gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
        
            plt.suptitle('ASCC {}, {}, $V={:.1f}$'.format(sID[i], sptype[ASCC==sID[i]][0], vmag[ASCC==sID[i]][0]), size='xx-large')
            
            ax1 = plt.subplot(gs[1,:])
            ax1.invert_yaxis()
            plt.errorbar(jdmid[~mask[i]], mag[i,~mask[i]] + trend[i,~mask[i]], emag[i,~mask[i]], fmt='.', c='k')
            plt.plot(jdmid[~mask[i]], trend[i,~mask[i]], c='r')
            plt.xlabel('Time [days]')
            plt.ylabel(r'$\Delta m$')
            plt.ylim(.5, -.5)
            
            ax2 = plt.subplot(gs[2,:], sharex=ax1)
            ax2.invert_yaxis()
            plt.errorbar(jdmid[~mask[i]], mag[i,~mask[i]], emag[i,~mask[i]], fmt='.', c='k')
            plt.xlabel('Time [days]')
            plt.ylabel(r'$\Delta m$')
            plt.ylim(.1, -.1)
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            else:
                plt.savefig(os.path.join(outdir, 'lightcurve_ASCC{}.png'.format(sID[i])))
        
            plt.close()

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
    
    data = ['/data3/talens/2015Q1/LPN/red0_vmag_2015Q1LPN.hdf5',
            '/data3/talens/2015Q1/LPE/red0_vmag_2015Q1LPE.hdf5',
            '/data3/talens/2015Q1/LPS/red0_vmag_2015Q1LPS.hdf5',
            '/data3/talens/2015Q1/LPW/red0_vmag_2015Q1LPW.hdf5',
            '/data3/talens/2015Q1/LPC/red0_vmag_2015Q1LPC.hdf5',
            '/data3/talens/2015Q2/LPN/red0_vmag_2015Q2LPN.hdf5',
            '/data3/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data3/talens/2015Q2/LPS/red0_vmag_2015Q2LPS.hdf5',
            '/data3/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5',
            '/data3/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5',
            '/data3/talens/2015Q3/LPN/red0_vmag_2015Q3LPN.hdf5',
            '/data3/talens/2015Q3/LPE/red0_vmag_2015Q3LPE.hdf5',
            '/data3/talens/2015Q3/LPS/red0_vmag_2015Q3LPS.hdf5',
            '/data3/talens/2015Q3/LPW/red0_vmag_2015Q3LPW.hdf5',
            '/data3/talens/2015Q3/LPC/red0_vmag_2015Q3LPC.hdf5',
            '/data3/talens/2015Q4/LPN/red0_vmag_2015Q4LPN.hdf5',
            '/data3/talens/2015Q4/LPE/red0_vmag_2015Q4LPE.hdf5',
            '/data3/talens/2015Q4/LPS/red0_vmag_2015Q4LPS.hdf5',
            '/data3/talens/2015Q4/LPW/red0_vmag_2015Q4LPW.hdf5',
            '/data3/talens/2015Q4/LPC/red0_vmag_2015Q4LPC.hdf5']
    
    
    filelist = glob.glob('/data3/talens/boxlstsq/2015Q1234/bls0*.hdf5')
    
    plot_candidates(data, filelist, '/data3/talens/boxlstsq/2015Q1234/new_figures')

    
    
    #data = ['/data3/talens/2015Q1/LPN/red0_vmag_2015Q1LPN.hdf5',
            #'/data3/talens/2015Q1/LPE/red0_vmag_2015Q1LPE.hdf5',
            #'/data3/talens/2015Q1/LPS/red0_vmag_2015Q1LPS.hdf5',
            #'/data3/talens/2015Q1/LPW/red0_vmag_2015Q1LPW.hdf5',
            #'/data3/talens/2015Q1/LPC/red0_vmag_2015Q1LPC.hdf5',
            #'/data3/talens/2015Q2/LPN/red0_vmag_2015Q2LPN.hdf5',
            #'/data3/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5',
            #'/data3/talens/2015Q2/LPS/red0_vmag_2015Q2LPS.hdf5',
            #'/data3/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5',
            #'/data3/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5',
            #'/data3/talens/2015Q3/LPN/red0_vmag_2015Q3LPN.hdf5',
            #'/data3/talens/2015Q3/LPE/red0_vmag_2015Q3LPE.hdf5',
            #'/data3/talens/2015Q3/LPS/red0_vmag_2015Q3LPS.hdf5',
            #'/data3/talens/2015Q3/LPW/red0_vmag_2015Q3LPW.hdf5',
            #'/data3/talens/2015Q3/LPC/red0_vmag_2015Q3LPC.hdf5',
            #'/data3/talens/2015Q4/LPN/red0_vmag_2015Q4LPN.hdf5',
            #'/data3/talens/2015Q4/LPE/red0_vmag_2015Q4LPE.hdf5',
            #'/data3/talens/2015Q4/LPS/red0_vmag_2015Q4LPS.hdf5',
            #'/data3/talens/2015Q4/LPW/red0_vmag_2015Q4LPW.hdf5',
            #'/data3/talens/2015Q4/LPC/red0_vmag_2015Q4LPC.hdf5']
    
    #filelist = glob.glob('/data3/talens/boxlstsq/2015Q1234/bls0*.hdf5')
    
    #new_periodogram(data, filelist, '444741')
    
    exit()
    
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

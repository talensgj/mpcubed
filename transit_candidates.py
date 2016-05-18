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

import pandas as pd

# For A4 paper.
rcParams['figure.figsize'] = (11.69,8.27)
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit
from package.plotting import viridis

import boxlstsq
from transit_search import read_header, read_data

from scipy import optimize

sptype_OC = {'O5':(13.4, -5.1), 'O6':(12.2, -5.1), 'O7':(11., -4.9), 'O8':(10., -4.6),
             'B0':(6.7, -3.4), 'B1':(5.2, -2.6), 'B2':(4.1, -1.6), 'B3':(3.8, -1.3), 'B5':(3.2, -.5), 'B6':(2.9, -.1), 'B7':(2.7, .3), 'B8':(2.5, .6), 'B9':(2.3, .8),
             'A0':(2.2, 1.1), 'A1':(2.1, 1.3), 'A2':(2., 1.5), 'A5':(1.8, 2.2), 'A8':(1.5, 2.7),
             'F0':(1.4, 3.), 'F2':(1.3, 3.4), 'F5':(1.2, 3.9), 'F8':(1.1, 4.3),
             'G0':(1.06, 4.7), 'G2':(1.03, 4.9), 'G8':(.96, 5.6),
             'K0':(.93, 5.7), 'K1':(.91, 6.), 'K3':(.86, 6.5), 'K4':(.83, 6.7), 'K5':(.8, 7.1), 'K7':(.74, 7.8),
             'M0':(.63, 8.9), 'M1':(.56, 9.6), 'M2':(.48, 10.4), 'M3':(.41, 11.1), 'M4':(.35, 11.9), 'M5':(.29, 12.8), 'M6':(.24, 13.8), 'M7':(.20, 14.7)}

sptype_OCinp = {'O5':(13.4, -5.1), 'O6':(12.2, -5.1), 'O7':(11., -4.9), 'O8':(10., -4.6),
                'B0':(6.7, -3.4), 'B1':(5.2, -2.6), 'B2':(4.1, -1.6), 'B3':(3.8, -1.3), 'B4':(3.5, -.9), 'B5':(3.2, -.5), 'B6':(2.9, -.1), 'B7':(2.7, .3), 'B8':(2.5, .6), 'B9':(2.3, .8),
                'A0':(2.2, 1.1), 'A1':(2.1, 1.3), 'A2':(2., 1.5), 'A3':(1.9, 1.7), 'A4':(1.9, 1.9), 'A5':(1.8, 2.2), 'A6':(1.7, 2.4), 'A7':(1.6, 2.6), 'A8':(1.5, 2.7),
                'F0':(1.4, 3.), 'F1':(1.35, 3.2), 'F2':(1.3, 3.4),  'F3':(1.3, 3.6), 'F4':(1.2, 3.8), 'F5':(1.2, 3.9), 'F6':(1.2, 4.0), 'F7':(1.1, 4.2), 'F8':(1.1, 4.3),
                'G0':(1.06, 4.7), 'G1':(1.04, 4.8), 'G2':(1.03, 4.9), 'G3':(1.02, 5.), 'G4':(1.01, 5.1), 'G5':(1., 5.2), 'G6':(.99, 5.3), 'G7':(.98, 5.4), 'G8':(.96, 5.6),
                'K0':(.93, 5.7), 'K1':(.91, 6.), 'K2':(.88, 6.3), 'K3':(.86, 6.5), 'K4':(.83, 6.7), 'K5':(.8, 7.1), 'K6':(.77, 7.4), 'K7':(.74, 7.8),
                'M0':(.63, 8.9), 'M1':(.56, 9.6), 'M2':(.48, 10.4), 'M3':(.41, 11.1), 'M4':(.35, 11.9), 'M5':(.29, 12.8), 'M6':(.24, 13.8), 'M7':(.20, 14.7)}

#def boxlstsq_header(filelist):
    
    #ascc = np.array([])
    #flag = np.array([], dtype='int')
    #period = np.array([])
    #depth = np.array([])
    #duration = np.array([])
    #nt = np.array([])
    
    #for filename in filelist:
        
        #with h5py.File(filename, 'r') as f:
            #grp = f['header']
            #ascc_ = grp['ascc'].value
            #flag_ = grp['flag'].value
            #period_ = grp['period'].value
            #depth_ = grp['depth'].value
            #duration_ = grp['duration'].value
            #nt_ = grp['nt'].value
            
            #ascc = np.append(ascc, ascc_)
            #flag = np.append(flag, flag_)
            #period = np.append(period, period_)
            #depth = np.append(depth, depth_)
            #duration = np.append(duration, duration_)
            #nt = np.append(nt, nt_)
            
    #return ascc, flag, period, depth, duration, nt

class StarCatalogue(object):
    
    def __init__(self, catalogue='/data3/talens/Catalogue/I_280B.hdf5'):
        
        cat = dict()
        
        with h5py.File(catalogue, 'r') as f:
        
            grp = f['data']
            
            ascc = grp['ASCC'].value 
            cat['ra'] = grp['RAhour'].value
            cat['dec'] = grp['DEdeg'].value
            cat['plx'] = grp['Plx'].value
            cat['sptype'] = grp['SpType'].value
            cat['vmag'] = grp['Vmag'].value
            cat['bmag'] = grp['Bmag'].value
            cat['tyc1'] = grp['TYC1'].value
            cat['tyc2'] = grp['TYC2'].value
            cat['tyc3'] = grp['TYC3'].value
            cat['hd'] = grp['HD'].value
    
        cat['plx'] = cat['plx']/1e2 # [mas]
        cat['vmag'] = cat['vmag']/1e3 # [mag]
        cat['bmag'] = cat['bmag']/1e3 # [mag]
    
        ascc = ascc.astype('|S32')
        self.cat = pd.DataFrame(cat, index=ascc)
    
        return
        
    def get_star(self, ascc0):
        
        return self.cat.loc[ascc0]

class RedFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
        
        return
        
    def read_header(self, fields):
        
        hdr = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['header']
            
            for field in fields:
                hdr[field] = grp[field].value
        
        return hdr
        
    def read_lightcurve(self, ascc, fields):
        
        lc = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data/'+ascc]
            
            for field in fields:
                lc[field] = grp[field].value
                
        return lc

class blsFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
    
        return
        
    def read_header(self, fields):
        
        hdr = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['header']
            
            for field in fields:
                hdr[field] = grp[field].value
                
        return hdr
        
    def read_data(self, fields):
        
        data = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data']
            
            for field in fields:
                data[field] = grp[field].value
        
        return data

#def plot_diagnostics(filelist):
    
    #ascc, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    #q = boxlstsq.phase_duration(1/period, 1., 1.)
    #select = (flag < 1)
    
    #ax = plt.subplot(311)
    #plt.plot(1/period, depth, '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(1/period[select], depth[select], c='r', zorder=1)
    #plt.xlim(0, 1.8)
    #plt.ylim(-.1, .1)
    #plt.xlabel(r'Frequency [day$^{-1}$]')
    #plt.ylabel('Depth')
    
    #plt.subplot(312, sharex=ax)
    #plt.plot(1/period, duration/(period*q), '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(1/period[select], duration[select]/(period[select]*q[select]), c='r', zorder=1)
    #plt.xlim(0, 1.8)
    #plt.ylim(0, 3.5)
    #plt.xlabel(r'Frequency [day$^{-1}$]')
    #plt.ylabel('Duration [q]')
    
    #plt.subplot(313, sharex=ax, yscale='log')
    #plt.plot(1/period, nt, '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(1/period[select], nt[select], c='r', zorder=1)
    #plt.xlim(0, 1.8)
    #plt.xlabel(r'Frequency [day$^{-1}$]')
    #plt.ylabel(r'$n_t$')
    
    #plt.tight_layout()
    #plt.show()
    
    #ax = plt.subplot(221, yscale='log')
    #plt.plot(duration/(period*q), nt, '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(duration[select]/(period[select]*q[select]), nt[select], c='r', zorder=1)
    #plt.xlim(0, 3.5)
    #plt.xlabel('Duration [q]')
    #plt.ylabel(r'$n_t$')
    
    #plt.subplot(223)
    #plt.plot(duration/(period*q), depth, '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(duration[select]/(period[select]*q[select]), depth[select], c='r', zorder=1)
    #plt.xlim(0, 3.5)
    #plt.ylim(-.1, .1)
    #plt.xlabel('Duration [q]')
    #plt.ylabel('Depth')
    
    #plt.subplot(224, xscale='log')
    #plt.plot(nt, depth, '.', c='k', alpha=.2, zorder=0)
    #plt.scatter(nt[select], depth[select], c='r', zorder=1)
    #plt.ylim(-.1, .1)
    #plt.ylabel('Depth')
    #plt.xlabel(r'$n_t$')
    
    #plt.tight_layout()
    #plt.show()
    
    #return 

def get_absmag(mag, d):
    
    Mag = mag - 5.*np.log10(d/10.)
    
    return Mag
    
def get_rstar(Mag, Mag0, R0):
    
    Rstar = R0*10**(-(Mag - Mag0)/5.)
    
    return Rstar
    
def get_rplanet(delta, Rstar):
    
    Rplanet = np.sqrt(delta/.01)*Rstar
    
    return Rplanet

def find_ns(lstseq):
    """ Number of sampling points in LST, takes wrapping into account."""
    
    lstidx = lstseq%270
    option1 = np.ptp(lstidx) + 1
    
    lstidx = np.mod(lstidx + 135, 270)
    option2 = np.ptp(lstidx) + 1
    
    if (option2 >= option1):
        return option1, False
    else:
        return option2, True
        
def transits(jdmid, period, epoch, duration):
    
    phase = (jdmid - epoch)/period
    cycle = np.floor(phase+.5)
    
    phase = np.mod(phase+.5, 1.)-.5
    condition = np.abs(phase) < .5*duration/period
    
    cycles = np.unique(cycle)
    
    length = []
    for i in cycles:
        
        sel = (cycle == i)
        nt = np.sum(condition[sel])

        if (nt > 0):
            length.append(nt)

    return length
    
#def transits(jdmid, pars):
    
    #phase = (jdmid - pars[1])/pars[0]
    #orbit = np.floor(phase + .5).astype('int')
    
    #phase = np.mod(phase+.5, 1.)-.5
    #sel = (np.abs(phase) < .5*pars[3]/pars[0])
    
    #intransit = np.bincount(orbit[sel])
    #intransit = intransit[intransit > 0]
    
    #return intransit
    
    
def rolling_mean(idx, y, window):
    
    nbin = np.bincount(idx)
    ybin = np.bincount(idx, y)
    
    nbin = np.cumsum(nbin)
    ybin = np.cumsum(ybin)
    
    nbin = np.append(0, nbin)
    ybin = np.append(0, ybin)
    
    nbin = (nbin[window:] - nbin[:-window])
    ybin = (ybin[window:] - ybin[:-window])/nbin
    
    return nbin, ybin


def plot_lightcurve(filename, ascc, outdir=None):
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    import detrend
    
    fields = ['lstseq', 'lst', 'jdmid', 'mag0', 'emag0', 'nobs']
    
    f = RedFile(filename)
    try:
        lc = f.read_lightcurve(ascc, fields)
    except:
        print 'Failed to read lightcurve.'
        return
    
    lstseq = lc['lstseq']
    lst = lc['lst']
    jdmid = lc['jdmid']
    mag0 = lc['mag0']
    emag0 = lc['emag0']
    nobs = lc['nobs']
    
    jdmid = jdmid - 2400000.5
    emag0 = emag0/np.sqrt(nobs)
    
    f_g = sum(nobs==50)/float(sum(nobs>25))

    sel = (nobs == 50)
    lstseq = lstseq[sel]
    lst = lst[sel]
    jdmid = jdmid[sel]
    mag0 = mag0[sel]
    emag0 = emag0[sel]
    
    if (len(jdmid) == 0):
        print 'No good datapoints.'
        return
        
    n1 = np.ptp(lstseq)
    n2, wrap = find_ns(lstseq)
    
    if wrap:
        lst = np.mod(lst+12, 24)-12
    
    n1 = np.maximum(n1, 2)
    n2 = np.maximum(n2, 2)
    
    weights = 1/(emag0)**2
    pars, trend1, trend2, chisq = detrend.new_harmonic2(jdmid, lst, mag0, weights, [n1, n2])
    
    head, tail = os.path.split(filename)
    tail = tail.rsplit('.')[0]
    tail = tail.rsplit('_')[-1]
    
    majorLocator1 = MultipleLocator(1)
    majorLocator2 = MultipleLocator(20)
    
    fig = plt.figure(figsize=(11.69, 8.27))

    gs = gridspec.GridSpec(3, 2, height_ratios = [1, 5, 5], width_ratios = [2, 3])

    plt.suptitle('{}, ASCC {}, $f_g = {:.2f}$ \n $\\chi^2_{{\\nu}} = {:.3f}$, npoints = {}, npars = {}'.format(tail, ascc, f_g, chisq/(len(jdmid) - len(pars)), len(jdmid), len(pars)), size='xx-large')
    
    ax1 = plt.subplot(gs[1,0])
    ax1.invert_yaxis()
    plt.errorbar(lst, mag0 - trend1, emag0, fmt='.', c='k')
    plt.plot(lst, trend2, '.', c='r')
    ax1.xaxis.set_major_locator(majorLocator1)
    plt.ylabel('$\Delta m$')
    
    ax2 = plt.subplot(gs[1,1], sharey=ax1)
    plt.errorbar(jdmid, mag0 - trend2, emag0, fmt='.', c='k')
    plt.plot(jdmid, trend1, '.', c='r')
    ax2.xaxis.set_major_locator(majorLocator2)
    
    ax3 = plt.subplot(gs[2,0], sharex=ax1)
    ax3.invert_yaxis()
    plt.errorbar(lst, mag0 - trend1 - trend2, emag0, fmt='.', c='k')
    plt.xlabel('Time [LST]')
    plt.ylabel('$\Delta m$')
    
    ax4 = plt.subplot(gs[2,1], sharex=ax2, sharey=ax3)
    plt.errorbar(jdmid, mag0 - trend1 - trend2, emag0, fmt='.', c='k')
    plt.xlabel('Time [MJD]')
    
    plt.tight_layout()
    
    if outdir is None:
        plt.show()
    else:
        outfile = 'lightcurve_ASCC{}_{}.png'.format(ascc, tail)
        outfile = os.path.join(outdir, outfile)
        plt.savefig(outfile)
    plt.close()

    return
    
    
def plot_periodogram(ascc, star, Nt, flag, jdmid, mag, emag, mask, freq, dchisq, pars, outdir=None):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    mask = mask[~mask]
    
    # Create the figure.
    fig = plt.figure(figsize=(8.27, 11.69))
    gs = gridspec.GridSpec(5, 2, height_ratios = [.5,10,10,10,5])

    # Create the title.
    line1 = r'ASCC {}, {}, $V={:.1f}$'.format(ascc, star['sptype'], star['vmag'])
    if (star['hd'] != 0):
        line1 = line1 + ', HD {}'.format(star['hd'])
    else:
        line1 = line1 + ', TYC {}-{}-{}'.format(star['tyc1'], star['tyc2'], star['tyc3'])
    line2 = '\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $N_t = {:.1f}$, flag={}'.format(pars[2]*100, pars[0], pars[3], Nt, flag)
    
    plt.suptitle(line1 + line2, size='xx-large')
    
    # Plot the periodogram.
    ax = plt.subplot(gs[1,:])
    
    plt.plot(freq, dchisq, c='k')
    
    freq = 1/pars[0]
    plt.axvline(freq, c='g', ls='--')
    for n in range(2, 5):
        plt.axvline(n*freq, c='g', ls='--')
        plt.axvline(freq/n, c='g', ls='--')
        
    freq = 1/.9972
    plt.axvline(freq, c='r', lw=2, ymax=.1)
    for n in range(2, 5):
        plt.axvline(n*freq, c='r', lw=2, ymax=.1)
        plt.axvline(freq/n, c='r', lw=2, ymax=.1)
        
    plt.axvline(freq, c='r', lw=2, ymax=.1)
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), c='r', lw=2, ymax=.1)
        plt.axvline((n+1)*freq/n, c='r', lw=2, ymax=.1)
        
    plt.xlim(0, 1.8)
    plt.xlabel('Frequency [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
    
    # Plot the data and the best-fit box-model.
    ax = plt.subplot(gs[2,:])
    
    phase = (jdmid - pars[1])/pars[0]
    phase = np.mod(phase + .5, 1) - .5
    mod = transit.box(jdmid, *pars)
    m = np.sum((mag - mod)/emag**2)/np.sum(1/emag**2)
    mod = mod + m
    
    phase_plot = np.linspace(-.5, .5, 9*pars[0]/pars[3])
    mod_plot = transit.box(phase_plot*pars[0]+pars[1], *pars)
    
    plt.errorbar(phase, mag, emag, fmt='.', c='k')
    plt.plot(phase_plot, mod_plot, c='r', lw=2)
    plt.xlim(-.5, .5)
    plt.ylim(.1, -.1)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    # Plot the binned data and the best-fit refined model.
    ax = plt.subplot(gs[3,:])
    
    sbox_pars, chisq, flag = refine_fit(jdmid, mag, emag, mask, pars)
    phase = (jdmid - sbox_pars[1])/sbox_pars[0]
    phase = np.mod(phase+.5, 1.)-.5
        
    npoints = 9*np.ceil(sbox_pars[0]/sbox_pars[3])
    if (npoints < 0):
        npoints = 9
    
    bins = np.linspace(-.5, .5, npoints+1)
    weights = 1/emag**2
    
    m0, bins = np.histogram(phase, bins=bins, weights=weights)
    m1, bins = np.histogram(phase, bins=bins, weights=weights*mag)
    
    phase_bin = (bins[:-1] + bins[1:])/2
    mag_bin = m1/m0
    emag_bin = np.sqrt(1/m0)
    
    phase_plot = np.linspace(-.5, .5, 9*sbox_pars[0]/sbox_pars[3])
    mod_plot = transit.softbox(phase_plot*sbox_pars[0]+sbox_pars[1], *sbox_pars)
        
    mod = transit.softbox(jdmid, *sbox_pars)
    mod_bin = transit.softbox(phase_bin*sbox_pars[0]+sbox_pars[1], *sbox_pars) 
        
    plt.title(r'$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days'.format(sbox_pars[2]*100, sbox_pars[0], sbox_pars[3]))
    plt.plot(phase, mag, '.', c='k', alpha=.5)
    plt.errorbar(phase_bin, mag_bin, emag_bin, fmt='o', c='g')
    plt.plot(phase_plot, mod_plot, c='r', lw=2)
    plt.xlim(-.5, .5)
    plt.ylim(.05, -.02)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    # Plot the residuals.
    ax = plt.subplot(gs[4,:])
    
    plt.plot(phase, mag-mod, '.', c='k', alpha=.5)
    plt.errorbar(phase_bin, mag_bin-mod_bin, emag_bin, fmt='o', c='g')
    plt.xlim(-.5, .5)
    plt.ylim(.02, -.02)
    plt.yticks([-.02, -.01, 0., .01, .02])
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
        
    plt.tight_layout()
    if outdir is None:
        plt.show()
    else:
        outfile = 'candidate_ASCC{}.png'.format(ascc)
        outfile = os.path.join(outdir, outfile)
        plt.savefig(outfile)
    plt.close()
    
    return

def plot_transit(ascc, star, Nt, flag, jdmid, mag, emag, mask, freq, dchisq, pars, outdir=None):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    mask = mask[~mask]
    
    # Create the figure.
    fig = plt.figure(figsize=(16, 5))
    gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

    # Plot the data and the best-fit box-model.
    ax = plt.subplot(gs[1,:])
    
    #phase = (jdmid - pars[1])/pars[0]
    #phase = np.mod(phase + .5, 1) - .5
    #mod = transit.box(jdmid, *pars)
    #m = np.sum((mag - mod)/emag**2)/np.sum(1/emag**2)
    #mod = mod + m
    
    #x0 = .5*pars[3]/pars[0]
    #phase_plot = np.array([-.5, -x0, -x0, x0, x0, .5])
    #mod_plot = np.array([0., 0., -pars[2], -pars[2], 0., 0.])
    
    #plt.errorbar(phase, mag, emag, fmt='.', c='k')
    #plt.plot(phase_plot, mod_plot, c='r', lw=2)
    #plt.xlim(-.5, .5)
    #plt.ylim(.1, -.1)
    #plt.xlabel('Phase')
    #plt.ylabel(r'$\Delta m$')
        
    sbox_pars, chisq, flag = refine_fit(jdmid, mag, emag, mask, pars)
    phase = (jdmid - sbox_pars[1])/sbox_pars[0]
    phase = np.mod(phase+.5, 1.)-.5
        
    npoints = 9*np.ceil(sbox_pars[0]/sbox_pars[3])
    bins = np.linspace(-.5, .5, npoints+1)
    weights = 1/emag**2
    
    m0, bins = np.histogram(phase, bins=bins, weights=weights)
    m1, bins = np.histogram(phase, bins=bins, weights=weights*mag)
    
    phase_bin = (bins[:-1] + bins[1:])/2
    mag_bin = m1/m0
    emag_bin = np.sqrt(1/m0)
    
    phase_plot = np.linspace(-.5, .5, 9*sbox_pars[0]/sbox_pars[3])
    mod_plot = transit.softbox(phase_plot*sbox_pars[0]+sbox_pars[1], *sbox_pars)
        
    mod = transit.softbox(jdmid, *sbox_pars)
    mod_bin = transit.softbox(phase_bin*sbox_pars[0]+sbox_pars[1], *sbox_pars) 
        
    print sbox_pars, flag
        
    plt.suptitle('ASCC {}, $P={:.3f}$ days, $\delta={:.1f}$ %'.format(ascc, sbox_pars[0], sbox_pars[2]*100), size='xx-large')
    plt.plot(phase, mag, '.', c='k', alpha=.8)
    plt.errorbar(phase_bin, mag_bin, emag_bin, fmt='o', c='r', ms=8)
    plt.plot(phase_plot, mod_plot, c='g', lw=4)
    plt.xlim(-.5, .5)
    plt.ylim(.05, -.02)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
        
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return

def phase_orbit(lstseq, phase, orbit, data, dx, dy, **kwargs):
    """ Make a phase-orbit plot of the given data."""
    
    # Segement the data fro plotting.
    args1, = np.where(np.diff(lstseq) > 1)
    args2, = np.where(np.diff(np.floor(orbit)) > 0)
    
    args = np.append(args1, args2)
    args = np.unique(args)
        
    args = args + 1
    args = np.append(0, args)
    args = np.append(args, len(data))
   
    # Plot each data segment.
    for i in range(len(args)-1):
    
        # Get the segment.
        x = phase[args[i]:args[i+1]]
        y = orbit[args[i]:args[i+1]]
        z = data[args[i]:args[i+1]]
       
        # Convert the center-cooridnates to edge-coordinates. 
        x = np.append(x - dx/2., x[-1] + dx/2.)
        y = np.append(y - dy/2., y[-1] - dy/2.)
        
        x = np.vstack([x, x])
        y = np.vstack([y, y + dy])
        
        # Give the data the appropriate shape.
        z = z[:, None].T
        
        # Plot the data.
        plt.pcolormesh(x, y, z, **kwargs)
    
    return

def plot_transit_coverage(lstseq, jdmid, lst, mag, emag, trend, mask, pars, ascc, star):
    
    import ephem
    
    lstseq = lstseq[~mask]
    jdmid = jdmid[~mask]
    lst = lst[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    trend = trend[~mask]
    
    # Make sure the arrays are sorted in time.
    sort = np.argsort(lstseq)
    lstseq = lstseq[sort]
    jdmid = jdmid[sort]
    lst = lst[sort]
    mag = mag[sort]
    emag = emag[sort]
    trend = trend[sort]
    
    # Compute the phase.
    phase = (jdmid - pars[1])/pars[0]
    orbit = phase + .5
    phase = np.mod(phase+.5, 1.)-.5
    
    # Compute the moon phase.
    m = ephem.Moon()
    mphase = np.zeros(len(jdmid))
    for i in range(len(jdmid)):
        m.compute(jdmid[i])
        mphase[i] = m.moon_phase
    
    # Set the cell sizes for the plotting the data.    
    dx = 319.1262613/(pars[0]*24*60*60)
    dy = 1.
    
    # Create the figure title.
    title = r'ASCC {}, $P = {:.4f}$ days'.format(ascc, pars[0])
    
    if (star['hd'] != 0):
        title = title + ', HD {}'.format(star['hd'])
    else:
        title = title + ', TYC {}-{}-{}'.format(star['tyc1'], star['tyc2'], star['tyc3'])
    
    # Create the figure.
    fig = plt.figure(figsize=(11.69, 8.27))
    gs = gridspec.GridSpec(3, 2, height_ratios=[.5,10,10])
    
    plt.suptitle(title, size='xx-large')
    
    plt.subplot(gs[1,0])
    phase_orbit(lstseq, phase, orbit, mag, dx, dy, vmin=-.02, vmax=.05, cmap=viridis)
    cbar = plt.colorbar()
    cbar.set_label(r'$\Delta m$')
    cbar.set_ticks(np.linspace(-.02, .05, 8))
    plt.axvline(.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.axvline(-.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.xlim(-.5, .5)
    plt.ylim(-.5, np.amax(orbit)+.5)
    plt.xlabel('Phase')
    plt.ylabel('Orbit')
    
    plt.subplot(gs[1,1])
    phase_orbit(lstseq, phase, orbit, trend, dx, dy, vmin=-.2, vmax=.2, cmap=viridis)
    cbar = plt.colorbar()
    cbar.set_label(r'Trend')
    cbar.set_ticks(np.linspace(-.2, .2, 9))
    plt.axvline(.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.axvline(-.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.xlim(-.5, .5)
    plt.ylim(-.5, np.amax(orbit)+.5)
    plt.xlabel('Phase')
    plt.ylabel('Orbit')
    
    plt.subplot(gs[2,0])
    phase_orbit(lstseq, phase, orbit, mphase, dx, dy, vmin=0, vmax=1, cmap=viridis)
    cbar = plt.colorbar()
    cbar.set_label('Moon Phase')
    cbar.set_ticks(np.linspace(0, 1, 6))
    plt.axvline(.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.axvline(-.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.xlim(-.5, .5)
    plt.ylim(-.5, np.amax(orbit)+.5)
    plt.xlabel('Phase')
    plt.ylabel('Orbit')
    
    plt.subplot(gs[2,1])
    phase_orbit(lstseq, phase, orbit, lst, dx, dy, vmin=0, vmax=24, cmap=viridis)
    cbar = plt.colorbar()
    cbar.set_label('Time [LST]')
    cbar.set_ticks(np.linspace(0, 24, 7))
    plt.axvline(.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.axvline(-.5*pars[3]/pars[0], c='k', lw=2, ls='--')
    plt.xlim(-.5, .5)
    plt.ylim(-.5, np.amax(orbit)+.5)
    plt.xlabel('Phase')
    plt.ylabel('Orbit')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return

def refine_fit(jdmid, mag, emag, mask, box_pars, display=False):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    
    # Evaluate the box-model.
    box_mod = transit.box(jdmid, *box_pars)
    
    # Find the out-of-transit level.
    m = np.sum((mag - box_mod)/emag**2)/np.sum(1/emag**2)
    box_mod = box_mod + m
    box_pars = np.append(box_pars, m)
    
    box_chisq = transit.box_chisq(box_pars, jdmid, mag, emag)
    
    # Initial guess for the refined fit.
    sbox_pars = np.insert(box_pars, -1, 10.)
    
    # Perform the refined fit.
    try:
        time0 = np.amin(jdmid)
        x = jdmid - time0
        sbox_pars[1] = sbox_pars[1] - time0
        sbox_pars, pcov = optimize.curve_fit(transit.softbox, x, mag, sbox_pars, sigma=emag, absolute_sigma=True)
        sbox_pars[1] = sbox_pars[1] + time0
    except:
        flag = 1
        sbox_mod = box_mod
        sbox_chisq = box_chisq
    else:
        flag = 0
        sbox_mod = transit.softbox(jdmid, *sbox_pars)
        sbox_chisq = transit.softbox_chisq(sbox_pars, jdmid, mag, emag)
    
    # Show the result.
    if display:
        ax = plt.subplot(211)
        ax.invert_yaxis()
        
        phase = (jdmid - box_pars[1])/box_pars[0]
        phase = np.mod(phase+.5, 1.)-.5 
        
        plt.plot(phase, mag, '.', c='k')
        plt.plot(phase, box_mod, '.', c='r')
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        ax = plt.subplot(212)
        ax.invert_yaxis()
        
        phase = (jdmid - sbox_pars[1])/sbox_pars[0]
        phase = np.mod(phase+.5, 1.)-.5 
        
        plt.plot(phase, mag, '.', c='k')
        plt.plot(phase, sbox_mod, '.', c='r')
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        plt.tight_layout()
        plt.show()
        plt.close()
    
    return sbox_pars, sbox_chisq, flag

def ellipsoidal_variations(jdmid, mag, emag, mask, pars, display=False):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    
    # Evaluate the box-model. 
    box_mod = transit.box_model(jdmid, *pars)
    
    # Find the out-of-transit level.
    m = np.sum((mag - box_mod)/emag**2)/np.sum(1/emag**2)
    box_mod = box_mod + m
    
    # Find in-transit points.
    phase = (jdmid - pars[1])/pars[0]
    phase = np.mod(phase+.5, 1.)-.5
    sel = (np.abs(phase) < .5*pars[3]/pars[0])
    
    # Create the model.
    theta = np.arccos(np.sin(2*np.pi*phase))
    p = -np.cos(2.*theta)
    
    # Compute the weights.
    weights = 1/emag**2
    weights[sel] = 0.
    
    # Find the amplitude and S/N of the model.
    u = np.sum(weights*(mag - box_mod)*p)
    v = np.sum(weights*p**2.)
    eps = u/v
    SNe = u/np.sqrt(v)
    
    if display:
        mod = np.where(sel, box_mod, eps*p)
        
        ax = plt.subplot(211)
        ax.invert_yaxis()
        
        plt.plot(phase, mag, '.', c='k')
        plt.plot(phase, mod, '.', c='r')
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        ax = plt.subplot(212)
        ax.invert_yaxis()
        
        plt.plot(phase, mag - mod, '.', c='k')
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
    return SNe, eps

def new_flags(star, pars):
    
    new_flag = 0
    
    # Check the transit duration.
    if (pars[3]/pars[0] > .15):
        new_flag += 1 
    
    # Check the parallax and planet radius.
    sptype = star['sptype']
    plx = star['plx']
    vmag = star['vmag']
    
    try:
        R, Mv = sptype_OCinp[sptype[:2]]
    except:
        Rstar = -1
        Rplanet = -1
    else:
        
        if (plx > 0):
            # Compute derived quantities.
            d = 1/(plx*1e-3) # [pc]
            Vmag = get_absmag(vmag, d)
            Rstar = get_rstar(Vmag, Mv, R)
            Rplanet = get_rplanet(-pars[2], Rstar)
            
            if (Rplanet > 10.):
                new_flag += 2
            
        else:
            Rstar = -1
            Rplanet = -1 
            new_flag += 4
    
    return new_flag, Rstar, Rplanet

def scale_sigma(jdmid, mag, emag, mask, pars):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    
    # Find the out-of-transit level.
    mod = transit.box(jdmid, *pars)
    m = np.sum((mag - mod)/emag**2)/np.sum(1/emag**2)
    pars = np.append(pars, m)
    
    # Compute the chi-square.
    chisq = transit.box_chisq(pars, jdmid, mag, emag)
    scale = chisq/(len(jdmid) - len(pars))
    
    # Find the in-transit points.
    phase = (jdmid - pars[1])/pars[0]
    phase = np.mod(phase+.5, 1.)-.5
    intransit = (np.abs(phase) < .5*pars[3]/pars[0])
    
    # Compute the S/N with the scaled errors.
    weights = 1/(scale*emag**2)
    
    t = np.sum(weights[~intransit])
    m = np.sum((weights*mag)[~intransit])
    
    r = np.sum(weights[intransit])
    s = np.sum((weights*mag)[intransit])
    
    SNr = (m/t - s/r)/np.sqrt(1/t + 1/r)
            
    return np.abs(SNr), scale

def red_noise(lstseq, jdmid, mag, emag, mask, pars):
    
    lstseq = lstseq[~mask]
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    
    # Subtract the box-model.
    mod = transit.box_model(jdmid, *pars)
    m = np.sum((mag - mod)/emag**2)/np.sum(1/emag**2)
    mod = mod + m
    mag = mag - mod
    
    # Number of points per transit.
    nt = transits(jdmid, pars[0], pars[1], pars[3])
    nt = np.array(nt)
    
    if (len(nt) == 0):
        return 0.
    
    # Window of the rolling mean.
    window = pars[3]/(319.1262613/(24*3600))
    window = np.ceil(window)

    # The rolling mean.
    #npbin, xbin = rolling_mean(lstseq_, jdmid_, window)
    npbin, ybin = rolling_mean(lstseq, mag, window)
    
    # Compute the variance.
    m0 = np.bincount(npbin)
    m1 = np.bincount(npbin, ybin)
    m2 = np.bincount(npbin, ybin**2)
    var = m2/m0 - (m1/m0)**2
    
    # Compute the red noise.
    SNr = (pars[2]*np.sum(nt))**2/np.sum(nt**2.*var[nt])
    
    return SNr

def candidates(data, filelist, outdir=None, ascc0=None):
    
    if outdir is not None:
        
        import csv
        from astropy.coordinates import ICRS
        from astropy import units as u
        
        head, tail = os.path.split(filelist[0])
        tail = tail.rsplit('_')[-2]
        outfile = tail + '.csv'
        outfile = os.path.join(outdir, outfile)
    
    # Create an instance of the ASCC catalogue.
    cat = StarCatalogue()
    
    for filename in filelist:
        print filename
        # Read the box least-squares results.
        f = blsFile(filename)
        
        fields = ['ascc', 'flag', 'period', 'depth', 'epoch', 'duration', 'nt']
        hdr = f.read_header(fields)
        ascc, flag, period, depth, epoch, duration, nt = tuple(hdr[field] for field in fields)
        
        fields = ['freq', 'dchisq']
        bls = f.read_data(fields)
        freq, dchisq = bls['freq'], bls['dchisq']

        # Check if there are good candidates in the file.
        if ascc0 is not None:
            if np.any(np.in1d(ascc, ascc0)):
                pass
            else:
                continue
        else:
            if np.any(flag == 0):
                pass
            else:
                continue

        # Compute some derived quantities.
        q = duration/period
        Nt = nt*(319.1262613/(24*3600))/duration
        SN = np.sqrt(np.amax(dchisq, axis=0))

        # Read the lightcurves.
        lstseq, jdmid, lst, mag, emag, mask, trend, cam = read_data(data, ascc)
        lstseq = lstseq.astype('int')
        #jdmid = jdmid - np.amin(jdmid)
        
        for i in range(len(ascc)):
            
            if ascc0 is not None:
                if (ascc[i] in ascc0):
                    pass
                else:
                    continue
            else:
                if (flag[i] == 0):
                    pass
                else:
                    continue
                
            star = cat.get_star(ascc[i])
            pars = np.array([period[i], epoch[i], -depth[i], duration[i]])
            
            # Plot the periodogram.
            plot_periodogram(ascc[i], star, Nt[i], flag[i], jdmid, mag[i], emag[i], mask[i], freq, dchisq[:,i], pars, outdir)
            #plot_transit(ascc[i], star, Nt[i], flag[i], jdmid, mag[i], emag[i], mask[i], freq, dchisq[:,i], pars, outdir)
            
            # Plot the transit coverage.
            for j in range(5):
                
                if np.any(cam == j):
                    print j
                    sel = (cam == j)
                    plot_transit_coverage(lstseq[sel], jdmid[sel], lst[sel], mag[i][sel], emag[i][sel], trend[i][sel], mask[i][sel], pars, ascc[i], star)
            
            # Get the star and planet radii.
            new_flag, Rstar, Rplanet = new_flags(star, pars)
            
            # Compute the ellipsoidal variations.
            SNe, eps = ellipsoidal_variations(jdmid, mag[i], emag[i], mask[i], pars)
            if (SNe > 10.):
                new_flag += 8
            
            # Compute signal-to-red-noise.
            SNr, scale = scale_sigma(jdmid, mag[i], emag[i], mask[i], pars)
            if (SNr < 5.): 
                new_flag += 16
        
            if outdir is not None:  
                
                # Format the coordinates.
                c = ICRS(star['ra']*u.hour, star['dec']*u.degree)
                dec = c.dec.to_string(u.degree, alwayssign=True, precision=0)
                ra = c.ra.to_string(u.hour, precision=0)
                
                # Build up the row.
                row = [ascc[i], ra, dec]
                fields = ['plx', 'sptype', 'vmag', 'bmag', 'hd', 'tyc1', 'tyc2', 'tyc3']
                row = row + [star[field] for field in fields]   
                row = row + [flag[i], period[i], epoch[i], depth[i], duration[i]]
                row = row + [new_flag, Rstar, Rplanet, SN[i], SNr, SNe]

                # Write the row to file.
                with open(outfile, 'a') as csvfile:
                    w = csv.writer(csvfile)
                    w.writerow(row)

    return

#def read_star(data, ascc):
    
    #import detrend
    
    #jdmid = np.array([])
    #lst = np.array([])
    #mag = np.array([])
    #emag = np.array([])
    #nobs = np.array([])
    #trend = np.array([])
    
    #for filename in data:
        
        #with h5py.File(filename, 'r') as f:
        
            #try:
                #grp = f['data/'+ascc]
            #except:
                #continue
            
            #lstseq = grp['lstseq'].value
            #jdmid_ = grp['jdmid'].value
            #lst_ = grp['lst'].value
            #mag_ = grp['mag0'].value
            #emag_ = grp['emag0'].value
            #nobs_ = grp['nobs'].value
        
        #emag_ = emag_/np.sqrt(nobs_)
        
        #sel = (nobs_ > 49)
        #lstseq = lstseq[sel]
        #jdmid_ = jdmid_[sel]
        #lst_ = lst_[sel]
        #mag_ = mag_[sel]
        #emag_ = emag_[sel]
        #nobs_ = nobs_[sel]
        
        #if (len(jdmid_) == 0):
            #continue
        
        #n1 = np.ptp(lstseq) + 1
        #n2 = np.ptp(lstseq%270) + 1
        
        #pars, trend_, chisq = detrend.new_harmonic(jdmid_, lst_, mag_, 1/emag_**2, [n1, n2])
        
        #jdmid = np.append(jdmid, jdmid_)
        #lst = np.append(lst, lst_)
        #mag = np.append(mag, mag_)
        #emag = np.append(emag, emag_)
        #nobs = np.append(nobs, nobs_)
        #trend = np.append(trend, trend_)
        
    #return jdmid, lst, mag, emag, nobs, trend

#def boxlstsq_refine(data, filelist):
    
    #for filename in filelist:
        
        #with h5py.File(filename, 'r') as f:
            
            #grp = f['header']
            #ascc = grp['ascc'].value
            #flag = grp['flag'].value
        
            #grp = f['data']
            #freq0 = grp['freq'].value
            #dchisq0 = grp['dchisq'].value
        
        #for i in range(len(ascc)):
            
            #if (flag[i] > 0): continue
            
            #jdmid, lst, mag, emag, nobs, trend = read_star(data, ascc[i])
            
            #if len(jdmid) == 0: continue
            
            ##freq1, chisq0, dchisq1, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag - trend, 1/emag**2, OS=1)
            ##freq2, chisq0, dchisq2, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag - trend, 1/emag**2, OS=3)
            ##freq3, chisq0, dchisq3, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag - trend, 1/emag**2, OS=5)
    
            ##plt.title('ASCC {}'.format(ascc[i]))
            ###plt.plot(freq0, dchisq0[:,i])
            ##plt.plot(freq1, dchisq1, label='OS=1')
            ##plt.plot(freq2, dchisq2, label='OS=3')
            ##plt.plot(freq3, dchisq3, label='OS=5')
            ##plt.legend()
            ##plt.xlabel(r'Frequency [day$^{-1}$]')
            ##plt.ylabel(r'$\Delta\chi^2$')
            ##plt.show()
            ##plt.close()
    
            #freq, chisq0, dchisq, hchisq, depth, epoch, duration, nt = boxlstsq.boxlstsq(jdmid, mag - trend, 1/emag**2, M=.5, R=.01, fmax=24.)
    
            #plt.title('ASCC {}'.format(ascc[i]))
            #plt.plot(freq, dchisq)
            #plt.xlabel(r'Frequency [day$^{-1}$]')
            #plt.ylabel(r'$\Delta\chi^2$')
            #plt.show()
            #plt.close()
    
    #return


def main():
    
    #filelist = glob.glob('/data3/talens/boxlstsq/2015Q2/bls/*')
    #filelist = np.sort(filelist)
    
    #ascc0 = np.array([])
    #for filename in filelist:
        #f = blsFile(filename)
        #hdr = f.read_header(['ascc', 'flag'])
        #ascc0 = np.append(ascc0, hdr['ascc'][hdr['flag']==0])
    
    data = glob.glob('/data3/talens/2015Q2/LP?/red0_vmag_2015Q?LP?.hdf5')
    data = np.sort(data)
    
    filelist = glob.glob('/data3/talens/boxlstsq/2015Q2/bls/*')
    filelist = np.sort(filelist)
    
    candidates(data, filelist, ascc0=['31753'])
    
    #data = glob.glob('/data3/talens/2015Q2_pea_ra/LP?/red0_pea_ra_2015Q?LP?.hdf5')
    #data = np.sort(data)
    
    #print data
    
    #for filename in data:
        #plot_lightcurve(filename, '31753')
    
    return

if __name__ == '__main__':
    main()

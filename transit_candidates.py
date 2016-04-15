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
    
    def __init__(self, redfile):
        
        self.redfile = redfile
        
        return
        
    def read_lightcurve(self, ascc, fields):
        
        lc = dict()
        
        with h5py.File(self.redfile, 'r') as f:
            grp = f['data/'+ascc]
            for field in fields:
                lc[field] = grp[field].value
                
        return lc

class BlsFiles(object):
    
    def __init__(self, filelist):
        
        filenames = np.sort(filelist)
        
        self._read_header(filelist)
        
        return
        
    def _read_header(self, filelist):
    
        ascc = np.array([])
        flag = np.array([], dtype='int')
        period = np.array([])
        epoch = np.array([])
        depth = np.array([])
        duration = np.array([])
        nt = np.array([])
        
        for filename in filelist:
            
            with h5py.File(filename, 'r') as f:
                
                grp = f['header']
                ascc_ = grp['ascc'].value
                flag_ = grp['flag'].value
                period_ = grp['period'].value
                epoch_ = grp['epoch'].value
                depth_ = grp['depth'].value
                duration_ = grp['duration'].value
                nt_ = grp['nt'].value
                
                ascc = np.append(ascc, ascc_)
                flag = np.append(flag, flag_)
                period = np.append(period, period_)
                epoch = np.append(epoch, epoch_)
                depth = np.append(depth, depth_)
                duration = np.append(duration, duration_)
                nt = np.append(nt, nt_)
            
        dic = dict()
        dic['flag'] = flag
        dic['period'] = period
        dic['epoch'] = epoch
        dic['depth'] = depth
        dic['duration'] = duration
        dic['nt'] = nt
        
        self.header = pd.DataFrame(dic, index=ascc)
        
        return
    
    def get_header(self):
        
        return self.header
    
    def get_candidates(self):
        
        ascc = self.header[self.header['flag'] == 0].index.tolist()
        
        return ascc
        
    def get_star(self, ascc0):
        
        return self.header.loc[ascc0]
        
    #def get_data(self, ascc0):
        
        #si = self.get_star(ascc0)
        #filename = si.filename

        #with h5py.File(filename, 'r') as f:
            
            #grp = f['data']
            #freq = grp['freq'].value
            #dchisq = grp['dchisq'].value
            
        #return freq, dchisq
   
def get_absmag(mag, d):
    
    Mag = mag - 5.*np.log10(d/10.)
    
    return Mag
    
def get_rstar(Mag, Mag0, R0):
    
    Rstar = R0*10**(-(Mag - Mag0)/5.)
    
    return Rstar
    
def get_rplanet(delta, Rstar):
    
    Rplanet = np.sqrt(delta/.01)*Rstar
    
    return Rplanet
    
def verify_candidates(filelist, outfile=None, verbose=True):
    
    import csv
    from astropy.coordinates import ICRS
    from astropy import units as u
    
    cat = StarCatalogue()
    bls = BlsFiles(filelist)
    ascc = bls.get_candidates()
    
    for ascc0 in ascc:
    
        # Get catalogue entry.
        cat_res = cat.get_star(ascc0)
        ra, dec, plx, sptype, vmag, bmag, tyc, hd = cat_res['ra'], cat_res['dec'], cat_res['plx'], cat_res['sptype'], cat_res['vmag'], cat_res['bmag'], (cat_res['tyc1'], cat_res['tyc2'], cat_res['tyc3']), cat_res['hd']
        
        # Get header bls result.
        bls_res = bls.get_star(ascc0)
        flag, period, epoch, depth, duration = bls_res['flag'], bls_res['period'], bls_res['epoch'], bls_res['depth'], bls_res['duration']
        
        new_flag = 0
        
        q = duration/period
        if (q > .15):
            new_flag += 1 
        
        # Try to get radius and absolute magnitude.
        try:
            R, Mv = sptype_OCinp[sptype[:2]]
        except:
            print 'Unknown type:', sptype
            Rstar = -1
            Rplanet = -1
        else:
            
            if (plx > 0):
                # Compute derived quantities.
                d = 1/(plx*1e-3) # [pc]
                Vmag = get_absmag(vmag, d)
                Rstar = get_rstar(Vmag, Mv, R)
                Rplanet = get_rplanet(depth, Rstar)
                
                if (Rplanet > 10.):
                    new_flag += 2
                
            else:
                Rstar = -1
                Rplanet = -1 
                new_flag += 4
        
        # Format the coordinates.
        c = ICRS(ra*u.hour, dec*u.degree)
        dec = c.dec.to_string(u.degree, alwayssign=True, precision=0)
        ra = c.ra.to_string(u.hour, precision=0)
        
        # Display the results.
        if verbose:
            print 'ASCC {}:'.format(ascc0) 
            print '  RA = {}, Dec = {}, Plx={} mas'.format(ra, dec, plx)
            print '  {}, V = {}, B = {}'.format(sptype, vmag, bmag)
            print '  HD {}'.format(hd)
            print '  TYC {}'.format(tyc)
            print
            print '  flag     = {}'.format(flag)
            print '  period   = {:.4f} days'.format(period)
            print '  epoch    = {:.4f}'.format(epoch)
            print '  depth    = {:.1f} %'.format(depth*100)
            print '  duration = {:.4f} days'.format(duration)    
            print 
            print '  flag = {}'.format(new_flag)
            print '  Rs = {:.2f} Rsun, Rp = {:.2f} Rjup'.format(Rstar, Rplanet)
            print '  eta/P = {:.2f}'.format(q)
            print
            
        if outfile is not None:   
            with open(outfile, 'a') as csvfile:
                w = csv.writer(csvfile)
                w.writerow([ascc0, ra, dec, plx, sptype, vmag, bmag, hd, tyc, flag, period, epoch, depth, duration, new_flag, Rstar, Rplanet])
    
    return

def plot_diagnostics(filelist):
    
    ascc, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    q = boxlstsq.phase_duration(1/period, 1., 1.)
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
    
    #new_periodogram(data, filelist, ascc, outdir)
    #plot_lightcurve(data, filelist, ascc, outdir)
    #phase_lst(data, filelist, ascc, outdir)
    
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

def find_ns(lstseq):
    """ Number of sampling points in LST, takes wrapping into account."""
    
    lstidx = lstseq%270
    option1 = np.ptp(lstidx) + 1
    
    lstidx = np.mod(lstidx + 135, 270)
    option2 = np.ptp(lstidx) + 1
    
    if (option2 > option1):
        return option1, False
    else:
        return option2, True

def plot_lightcurve(filename, ascc, outdir=None):
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    import detrend
    
    fields = ['lstseq', 'lst', 'jdmid', 'mag0', 'emag0', 'nobs']
    
    f = RedFile(filename)
    try:
        lc = f.read_lightcurve(ascc, fields)
    except:
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
        return
        
    n1 = np.ptp(lstseq)
    n2, wrap = find_ns(lstseq)
    
    n1 = np.maximum(n1, 2)
    n2 = np.maximum(n2, 2)
    
    weights = 1/(emag0)**2
    pars, trend1, trend2, chisq = detrend.new_harmonic2(jdmid, lst, mag0, weights, [n1, n2])
    
    head, tail = os.path.split(filename)
    tail = tail.rsplit('.')[0]
    tail = tail.rsplit('_')[-1]
    
    majorLocator1 = MultipleLocator(1)
    majorLocator2 = MultipleLocator(20)
    
    if wrap:
        lst = np.mod(lst+12, 24)-12
    
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

def phase_lst(data, filelist, ascc, outdir=None):
    
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
                
            phase = (jdmid - best_epoch[i])/period[i]
            phase = np.mod(phase + .5, 1) -.5
            
            fig = plt.figure(figsize=(8.27, 4))
            plt.title(r'ASCC {}, $P={:.4f}$ days'.format(sID[i], period[i]))
            plt.plot(phase[~mask[i]], lst[~mask[i]], '.', c='k', alpha=.5)
            plt.xlim(-.5, .5)
            plt.ylim(0, 24)
            plt.xlabel('Phase')
            plt.ylabel('Time [LST]')
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            elif (np.amax(dchisq[:,i]) < 400):
                plt.savefig(os.path.join(outdir, 'maybe_candidate_ASCC{}_p2.png'.format(sID[i])))
            elif (np.amax(dchisq[:,i]) >= 900):
                plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}_p2.png'.format(sID[i])))
            else:
                plt.savefig(os.path.join(outdir, 'candidate_ASCC{}_p2.png'.format(sID[i])))
            plt.close()
                
    return

def phase_sigma(data, filelist, ascc, outdir=None):
    
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
                
            phase = (jdmid - best_epoch[i])/period[i]
            phase = np.mod(phase + .5, 1) -.5
            
            fig = plt.figure(figsize=(8.27, 4))
            plt.title(r'ASCC {}, $P={:.4f}$ days'.format(sID[i], period[i]))
            plt.plot(phase[~mask[i]], emag[i,~mask[i]], '.', c='k', alpha=.5)
            plt.xlim(-.5, .5)
            plt.ylim(0, .05)
            plt.xlabel('Phase')
            plt.ylabel(r'$\sigma$')
            
            plt.tight_layout()
            
            if outdir is None:
                plt.show()
            elif (np.amax(dchisq[:,i]) < 400):
                plt.savefig(os.path.join(outdir, 'maybe_candidate_ASCC{}_p2.png'.format(sID[i])))
            elif (np.amax(dchisq[:,i]) >= 900):
                plt.savefig(os.path.join(outdir, 'prime_candidate_ASCC{}_p2.png'.format(sID[i])))
            else:
                plt.savefig(os.path.join(outdir, 'candidate_ASCC{}_p2.png'.format(sID[i])))
            plt.close()

    return


def main():
    
    data = glob.glob('/data3/talens/2015Q?/LP?/red0_vmag_2015Q?LP?.hdf5')
    data = np.sort(data)
    
    filelist = glob.glob('/data3/talens/boxlstsq/2015Q1234/bls0*.hdf5')
    filelist = np.sort(filelist)
    
    verify_candidates(filelist, outfile = '/data3/talens/2015Q1234.csv')
    
    #ascc, flag, period, depth, duration, nt = boxlstsq_header(filelist)
    
    #ascc = ascc[flag == 0] 
    
    #for filename in data:
        #for ascc0 in ascc:
            #plot_lightcurve(filename, ascc0, '/data3/talens/boxlstsq/2015Q1234/test')
    
    return

if __name__ == '__main__':
    main()

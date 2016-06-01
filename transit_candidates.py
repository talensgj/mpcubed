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

from mpcubed import IO
from mpcubed import misc
from mpcubed.models import transit

from transit_search import read_header, read_data
    
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
        R, Mv = misc.sptype_OCinp[sptype[:2]]
    except:
        Rstar = -1
        Rplanet = -1
    else:
        
        if (plx > 0):
            # Compute derived quantities.
            d = 1/(plx*1e-3) # [pc]
            Vmag = misc.get_absmag(vmag, d)
            Rstar = misc.get_rstar(Vmag, Mv, R)
            Rplanet = misc.get_rplanet(-pars[2], Rstar)
            
            if (Rplanet > 10.):
                new_flag += 2
            
        else:
            Rstar = -1
            Rplanet = -1 
            new_flag += 4
    
    return new_flag, Rstar, Rplanet

def ellipsoidal_variations(jdmid, mag, emag, mask, pars, display=False):
    
    jdmid = jdmid[~mask]
    mag = mag[~mask]
    emag = emag[~mask]
    
    # Evaluate the box-model. 
    box_mod = transit.box(jdmid, *pars)
    
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

def refine_fit(jdmid, mag, emag, mask, box_pars, display=False):
    
    from scipy import optimize
    
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

def plot_candidate(ascc, star, freq, dchisq, jdmid, mag, emag, mask, pars, Nt, flag, outdir=None):
    
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
    
    plt.plot(freq, dchisq, c=(0,0,0))
    
    freq = 1/pars[0]
    plt.axvline(freq, c=(0./255,109./255,219./255), ls='--')
    for n in range(2, 5):
        plt.axvline(n*freq, c=(0./255,109./255,219./255), ls='--')
        plt.axvline(freq/n, c=(0./255,109./255,219./255), ls='--')
        
    freq = 1/.9972
    plt.axvline(freq, c=(146./255,0,0), lw=2, ymax=.1)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline(freq/n, c=(146./255,0,0), lw=2, ymax=.1)
        
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline((n+1)*freq/n, c=(146./255,0,0), lw=2, ymax=.1)
        
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
    
    plt.errorbar(phase, mag, emag, fmt='.', c=(0,0,0))
    plt.plot(phase_plot, mod_plot, c=(146./255,0,0), lw=2)
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
    plt.plot(phase, mag, '.', c=(0,0,0), alpha=.5)
    plt.errorbar(phase_bin, mag_bin, emag_bin, fmt='o', c=(0./255,109./255,219./255))
    plt.plot(phase_plot, mod_plot, c=(146./255,0,0), lw=2)
    plt.xlim(-.5, .5)
    plt.ylim(.05, -.02)
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$')
    
    # Plot the residuals.
    ax = plt.subplot(gs[4,:])
    
    plt.plot(phase, mag-mod, '.', c=(0,0,0), alpha=.5)
    plt.errorbar(phase_bin, mag_bin-mod_bin, emag_bin, fmt='o', c=(0./255,109./255,219./255))
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
    
def candidates(directory, outdir='candidates'):
    
    import csv
    from astropy.coordinates import ICRS
    from astropy import units as u
    
    # Check that the directory exists.
    if not os.path.exists(directory):
        print 'Error: path {} does not exist.'
        exit()
    
    # Find out on which data the bls was run.
    data = np.loadtxt(os.path.join(directory, 'data.txt'), dtype='S')
    
    # Check that the bls directory exists.
    blsdir = os.path.join(directory, 'bls')
    if not os.path.exists(blsdir):
        print 'Error: box least-squares directory does not exist.'
        exit()
    
    # Check that there are bls files.
    filelist = os.listdir(blsdir)
    if (len(filelist) == 0):
        print 'Error: No bls files found.'    
        exit()
    
    # Create the candidates directory.
    outdir = os.path.join(directory, outdir) 
    misc.ensure_dir(outdir)
    if (len(os.listdir(outdir)) > 0):
        print 'Error: the output directory {} is not empty.'.format(outdir)
        exit()
    else:
        print 'Writing results to:', outdir
    
    # Name of the output file.
    head, tail = os.path.split(directory)
    outfile = os.path.join(outdir, tail + '.csv')
    
    # Create an instance of the ASCC catalogue.
    cat = IO.StarCat()
    
    for filename in filelist:
        
        filename = os.path.join(blsdir, filename)
        
        # Read the box least-squares results.
        f = IO.blsFile(filename)
        
        fields = ['ascc', 'flag', 'period', 'depth', 'epoch', 'duration', 'nt']
        hdr = f.read_header(fields)
        ascc, flag, period, depth, epoch, duration, nt = tuple(hdr[field] for field in fields)
        
        if np.all(flag > 0):
            continue
        
        fields = ['freq', 'dchisq']
        bls = f.read_data(fields)
        freq, dchisq = bls['freq'], bls['dchisq']

        # Compute some derived quantities.
        q = duration/period
        Nt = nt*(319.1262613/(24*3600))/duration
        SN = np.sqrt(np.amax(dchisq, axis=0))

        # Read the lightcurves.
        lstseq, jdmid, lst, mag, emag, mask, trend, cam = read_data(data, ascc)
        lstseq = lstseq.astype('int')
        
        for i in range(len(ascc)):
            
            if (flag[i] != 0):
                continue
            
            star = cat.get_star(ascc[i])
            pars = np.array([period[i], epoch[i], -depth[i], duration[i]])
            
            # Plot the periodogram.
            plot_candidate(ascc[i], star, freq, dchisq[:,i], jdmid, mag[i], emag[i], mask[i], pars, Nt[i], flag[i], outdir)
            
            # Get the star and planet radii.
            new_flag, Rstar, Rplanet = new_flags(star, pars)
            
            # Compute the ellipsoidal variations.
            SNe, eps = ellipsoidal_variations(jdmid, mag[i], emag[i], mask[i], pars)
            if (SNe > 10.):
                new_flag += 8
            
            # Compute signal-to-red-noise.
            SNr = 0
            scale = 0
            #SNr, scale = scale_sigma(jdmid, mag[i], emag[i], mask[i], pars)
            #if (SNr < 5.): 
                #new_flag += 16
                
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

def main():
    
    candidates('/data3/talens/boxlstsq/test')
    
    return

if __name__ == '__main__':
    main()

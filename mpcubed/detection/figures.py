# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os
import glob

import numpy as np
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from .. import io, misc
from . import boxlstsq

def plot_periodogram(freq, dchisq, period, zoom=False):
    """ Plot the box least-squares periodogram. """

    plt.annotate(r'$P = {:.5f}$'.format(period), (0, 1), xytext=(10, -10), xycoords='axes fraction', textcoords='offset points', va='top', ha='left', size='x-large', backgroundcolor='w')    
    
    # Plot the box least-squares periodogram.
    plt.plot(freq, dchisq, c='k')
    if zoom:
        plt.xlim(.95/period, 1.05/period)
    else:
        plt.xlim(0, 1.8)
    plt.xlabel(r'Frequency [day$^{-1}$]')
    plt.ylabel(r'$\Delta\chi^2$')
    
    # Add lines indicating the peak, and integer harmonics.
    freq = 1/period
    plt.axvline(freq, c=(0./255,109./255,219./255), zorder=-10)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(0./255,109./255,219./255), ls='--', zorder=-10)
        plt.axvline(freq/n, c=(0./255,109./255,219./255), ls='--', zorder=-10)
        
    # Add lines indicating the 1 day systematic and harmonics.
    freq = 1/.9972
    plt.axvline(freq, c=(146./255,0,0), lw=2, ymax=.1)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline(freq/n, c=(146./255,0,0), lw=2, ymax=.1)
        
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), c=(146./255,0,0), lw=2, ymax=.1)
        plt.axvline((n+1)*freq/n, c=(146./255,0,0), lw=2, ymax=.1)     
    
    return  
    
def plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration, binned=True, zoom=False):
    """ Plot a lightcurve with the box least-squares best fit overplotted. """    
    
    factor = np.ceil(np.abs(depth)/.05)   

    # Plot the phase-folded lightcurve.
    phase = (jdmid - epoch)/period    
    phase = np.mod(phase+.5, 1.)-.5
    plt.scatter(phase, mag, color='black', marker='.', alpha=.5, edgecolor='none')
    
    # Add phase-binned data.
    if binned:
        nbins = np.ceil(9*period/duration)
        bins = np.linspace(-.5, .5, nbins+1)   
        xbin, ybin, eybin = misc.bin_data_err(phase, mag, emag, bins)
        plt.errorbar(xbin, ybin, eybin, fmt='o', c=(0./255,109./255,219./255))
    
    if zoom:
        plt.xlim(-1.5*duration/period, 1.5*duration/period)
        plt.ylim(2*np.abs(depth), -2*np.abs(depth))
    else:
        plt.xlim(-.5, .5)
        plt.ylim(factor*.05, factor*-.03)
    
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$') 
    
    return
    
def plot_boxcurve(period, depth, duration):
    
    # Compute x and y values for plotting a boxfit.
    x = .5*duration/period        
    x = np.array([-.5, -x, -x, x, x, .5])
    y = np.array([0, 0, depth, depth, 0, 0,])
    
    plt.plot(x, y, c=(146./255,0,0), lw=2, zorder=20)    
    
    return

def fig_candidate(ascc, freq, dchisq, jdmid, mag, emag, period, epoch, depth, duration, figdir):
    
    # Create the figure.
    fig = plt.figure(figsize=(16.5, 11.7))
    
    plt.suptitle('ASCC {}'.format(ascc), size='xx-large')    
    
    gs = gridspec.GridSpec(6, 6, height_ratios = [.5, 10, 10, .5, 10, .5])
    
    # Plot the periodogram.
    ax = plt.subplot(gs[1,:4])
    
    plt.title('Periodogram')
    plot_periodogram(freq, dchisq, period)                     
    
    ax = plt.subplot(gs[1,4:])
    
    plt.title('Periodogram zoom')
    plot_periodogram(freq, dchisq, period, zoom=True)    
    
    # Plot the data and the best-fit box-model.
    ax = plt.subplot(gs[2:4,:4])

    plt.title('Photometry')
    plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration)
    plot_boxcurve(period, depth, duration)
    
    ax = plt.subplot(gs[2:4,4:])

    plt.title('Photometry zoom')
    plot_lightcurve(jdmid, mag, emag, period, epoch, depth, duration, zoom=True)
    plot_boxcurve(period, depth, duration)

    ax = plt.subplot(gs[4,:3])

    plt.title(r'Half-period, $P = {:.5f}$'.format(.5*period))
    plot_lightcurve(jdmid, mag, emag, .5*period, epoch, depth, duration, binned=False)

    ax = plt.subplot(gs[4,3:])

    plt.title(r'Double-period, $P = {:.5f}$'.format(2.*period))
    plot_lightcurve(jdmid, mag, emag, 2.*period, epoch, depth, duration, binned=False)        
        
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, 'ASCC{}.png'.format(ascc)))
    plt.close()    
    
    return
 
def worker(queue, figdir):

    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:
            
            hdr, data, jd, mag, emag, mask = item
            
            args, = np.where(hdr['depth'] > 0)
            
            for i in args:
                fig_candidate(hdr['ascc'][i], data['freq'], data['dchisq'][:,i], jd[mask[i]], mag[i,mask[i]], emag[i,mask[i]], hdr['period'][i], hdr['epoch'][i], hdr['depth'][i], hdr['duration'][i], figdir)
    
    return    

def figs_boxlstsq(blsdir, aper=0, method='legendre', nprocs=6): 

    # Get the box least-squares files.
    blsfiles = glob.glob(os.path.join(blsdir, 'bls/*'))
    blsfiles = np.sort(blsfiles)
    
    # Get the photometry files.
    filelist = np.genfromtxt(os.path.join(blsdir, 'data.txt'), dtype='S') 

    # Create the output directory.
    figdir = os.path.join(blsdir, 'figures')
    misc.ensure_dir(figdir)

    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, worker, (the_queue, figdir))

    for blsfile in blsfiles:
        
        print blsfile

        # Read the box least-squares results.
        f = io.blsFile(blsfile)
        hdr = f.read_header(['ascc', 'period', 'epoch', 'depth', 'duration'])  
        data = f.read_data(['freq', 'dchisq'])
    
        # Read the lightcurves.
        jd, lst, mag, emag, trend, mask = boxlstsq.read_data(filelist, hdr['ascc'], np.zeros(len(hdr['ascc'])), aper=aper, method=method)
        
        the_queue.put((hdr, data, jd, mag, emag, mask))
        
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()

    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Make figures for the database.')
    parser.add_argument('blsdir', type=str,
                        help='the BLS run to create the figures for')
    parser.add_argument('-a', '--aper', type=int, default=0,
                        help ='the aperture to use', dest='aper')
    parser.add_argument('-m', '--method', type=str, default='legendre',
                        help ='detrending method', dest='method')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()

    figs_boxlstsq(args.blsdir, args.aper, args.method, args.nprocs)
    
    return

if __name__ == '__main__':
    
    main()

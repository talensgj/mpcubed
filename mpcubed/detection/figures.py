# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 16:17:22 2017

@author: talens
"""

import os
import glob
import argparse
import multiprocessing as mp

import numpy as np

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

from .. import io, statistics
from . import detrend


def plot_periodogram(freq, dchisq, period, zoom=False):
    """ Plot the box least-squares periodogram. """

    plt.annotate(r'$P = {:.5f}$'.format(period), (0, 1), xytext=(10, -10), xycoords='axes fraction',
                 textcoords='offset points', va='top', ha='left', size='x-large', backgroundcolor='w')
    
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
    plt.axvline(freq, c=(0./255, 109./255, 219./255), zorder=-10)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(0./255, 109./255, 219./255), ls='--', zorder=-10)
        plt.axvline(freq/n, c=(0./255, 109./255, 219./255), ls='--', zorder=-10)
        
    # Add lines indicating the 1 day systematic and harmonics.
    freq = 1/.9972
    plt.axvline(freq, c=(146./255, 0, 0), lw=2, ymax=.1)
    for n in range(2, 5):
        plt.axvline(n*freq, c=(146./255, 0, 0), lw=2, ymax=.1)
        plt.axvline(freq/n, c=(146./255, 0, 0), lw=2, ymax=.1)
        
    for n in range(1, 5):
        plt.axvline(n*freq/(n+1), c=(146./255, 0, 0), lw=2, ymax=.1)
        plt.axvline((n+1)*freq/n, c=(146./255, 0, 0), lw=2, ymax=.1)
    
    return  


def plot_lightcurve(jd, mag, emag, box_pars, offset=0.5, binned=True, zoom=False):
    """ Plot a lightcurve with the box least-squares best fit overplotted. """    
    
    T0, P, T14, delta = box_pars
    
    factor = np.ceil(np.abs(delta)/.05)   

    # Plot the phase-folded lightcurve.
    phase = (jd - T0)/P    
    phase = np.mod(phase+offset, 1.)-offset
    plt.scatter(phase, mag, color='black', marker='.', alpha=.5, edgecolor='none')
    
    # Add phase-binned data.
    if binned:
        nbins = np.ceil(9*P/T14)
        bins = np.linspace(-offset, 1-offset, nbins+1)   
        xbin, ybin, eybin = statistics.bin_data(phase, mag, bins, emag)
        plt.errorbar(xbin, ybin, eybin, fmt='o', c=(0./255, 109./255, 219./255))
    
    if zoom:
        plt.xlim(-1.5*T14/P, 1.5*T14/P)
        plt.ylim(2*np.abs(delta), -2*np.abs(delta))
    else:
        plt.xlim(-offset, 1-offset)
        plt.ylim(factor*.05, factor*-.03)
    
    plt.xlabel('Phase')
    plt.ylabel(r'$\Delta m$') 
    
    return


def plot_boxcurve(box_pars):
    
    T0, P, T14, delta = box_pars
    
    # Compute x and y values for plotting a boxfit.
    x = .5*T14/P        
    x = np.array([-.5, -x, -x, x, x, .5])
    y = np.array([0, 0, -delta, -delta, 0, 0])
    
    plt.plot(x, y, c=(146./255, 0, 0), lw=2, zorder=20)
    
    return


def figs_candidates(hdr, data, lc2d, nobs, method, figdir):
    
    args, = np.where(hdr['sde'] >= 5)
    
    box_pars = np.column_stack([hdr['epoch'], hdr['period'], hdr['duration'], hdr['depth']])
    
    for i in args:
                
        # Perform the secondary calibration.
        lc = detrend.remove_trend(lc2d[:, i], nobs, method=method)
        lc = lc[lc['mask']]
    
        # Create the figure.
        plt.figure(figsize=(16.5, 11.7))
        
        plt.suptitle('ASCC {}'.format(hdr['ascc'][i]), size='xx-large')    
        
        gs = gridspec.GridSpec(6, 6, height_ratios=[.5, 10, 10, .5, 10, .5])
        
        # Plot the periodogram.
        plt.subplot(gs[1, :4])
        
        plt.title('Periodogram')
        plot_periodogram(data['freq'], data['dchisq_dec'][:, i], box_pars[i, 1])
        
        plt.subplot(gs[1, 4:])
        
        plt.title('Periodogram zoom')
        plot_periodogram(data['freq'], data['dchisq_dec'][:, i], box_pars[i, 1], zoom=True)
        
        # Plot the data and the best-fit box-model.
        plt.subplot(gs[2:4, :4])
    
        plt.title('Photometry')
        plot_lightcurve(lc['jd'], lc['mag'] - lc['trend'], lc['emag'], box_pars[i])
        plot_boxcurve(box_pars[i])
        
        plt.subplot(gs[2:4, 4:])
    
        plt.title('Photometry zoom')
        plot_lightcurve(lc['jd'], lc['mag'] - lc['trend'], lc['emag'], box_pars[i], zoom=True)
        plot_boxcurve(box_pars[i])
    
        plt.subplot(gs[4, :3])
    
        box_pars[i, 1] = .5*box_pars[i, 1]  # Half original period.
        plt.title(r'Half-period, $P = {:.5f}$'.format(box_pars[i, 1]))
        plot_lightcurve(lc['jd'], lc['mag'] - lc['trend'], lc['emag'], box_pars[i], binned=False)
    
        plt.subplot(gs[4, 3:])
    
        box_pars[i, 1] = 4.*box_pars[i, 1]  # Twice original period.
        plt.title(r'Double-period, $P = {:.5f}$'.format(box_pars[i, 1]))
        plot_lightcurve(lc['jd'], lc['mag']-lc['trend'], lc['emag'], box_pars[i], offset=0.25, binned=False)        
            
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, 'ASCC{}.png'.format(hdr['ascc'][i])))
        plt.close()    
    
    return


def worker(queue, method, figdir):

    while True:
        
        item = queue.get()
        
        if item == 'DONE':
            break
        else:

            # Unpack the arguments.
            hdr, data, tmpfile = item

            # Read the temporary file.
            npzfile = np.load(tmpfile)
            ascc, time, lc2d, nobs = npzfile['ascc'], npzfile['time'], npzfile['lc2d'], npzfile['nobs']

            # Remove the temporary file.
            os.remove(tmpfile)

            # Make the diagnostic figures.
            figs_candidates(hdr, data, lc2d, nobs, method, figdir)
    
    return    


def figs_boxlstsq(blsdir, aper=0, method='legendre', nprocs=6): 

    from .boxlstsq import read_data
    
    # Get the box least-squares files.
    blsfiles = glob.glob(os.path.join(blsdir, 'bls/*'))
    blsfiles = np.sort(blsfiles)
    
    # Get the photometry files.
    filelist = np.genfromtxt(os.path.join(blsdir, 'data.txt'), dtype='S', delimiter=',')

    # Create the output directory.
    figdir = os.path.join(blsdir, 'figures')
    io.ensure_dir(figdir)

    # Set up the multiprocessing.
    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, worker, (the_queue, method, figdir))

    print('Making boxlstsq diagnostic figures for:')

    for blsfile in blsfiles:
        
        print(' ', blsfile)

        # Read the box least-squares results.
        f = io.BoxlstsqFile(blsfile)
        hdr = f.read_header(['ascc', 'period', 'epoch', 'depth', 'duration', 'sde'])
        data = f.read_data(['freq', 'dchisq_dec'])
    
        # Read the lightcurves.
        time, lc2d, nobs = read_data(filelist, hdr['ascc'], aper=aper)

        # Save the lightcurves to a temporary file.
        # Data volumes are too large to put directly in the queue.
        tmpfile = os.path.splitext(blsfile)[0] + '.npz'
        np.savez(tmpfile, ascc=hdr['ascc'], time=time, lc2d=lc2d, nobs=nobs)

        the_queue.put((hdr, data, tmpfile))
        
    # End the multiprocessing.
    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()

    return


def main():

    parser = argparse.ArgumentParser(description='Make figures for the database.')
    parser.add_argument('blsdir', type=str,
                        help='the BLS run to create the figures for')
    parser.add_argument('-a', '--aper', type=int, default=0,
                        help='the aperture to use', dest='aper')
    parser.add_argument('-m', '--method', type=str, default='legendre',
                        help='detrending method', dest='method')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()

    figs_boxlstsq(args.blsdir, args.aper, args.method, args.nprocs)
    
    return


if __name__ == '__main__':
    
    main()

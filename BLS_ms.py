#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from time import time

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m]  

def cumsum_to_grid(time, bin_edges, weights):

    nt = len(time)
    nedges = len(bin_edges)
    nw = len(weights)
    
    if (np.diff(bin_edges) < 0).any():
        raise AttributeError('bins must increase monotonically.')
    
    if nt != nw:
        raise ValueError('weights should have the same first dimension as a.')
    
    if weights.ndim>1:
        nstars = np.shape(weights)[1]
        n = np.zeros((nedges,nstars), weights.dtype)
        zero = np.array([0]*nstars, weights.dtype)
    else:
        n = np.zeros((nedges,), weights.dtype)
        zero = np.array(0, weights.dtype)
    
    block = 65536
    for i in np.arange(0, len(time), block):
        sa = time[i:i+block]
        sw = weights[i:i+block]
        cw = np.concatenate(([zero,], sw.cumsum(axis=0)))
        bin_index = np.r_[sa.searchsorted(bin_edges[:-1], 'left'), sa.searchsorted(bin_edges[-1], 'right')]
        n += cw[bin_index]
    
    return n

def BLS_ms(time, flux, flux_err, **options):
    """
    Perform Box Least-Squares on a given lightcurve.

    Returns the maximum Signal Residue as a function of orbital frequency.

    Parameters
    ----------
    time : array_like
        Array contaning the measurement times in units of [days].
    flux : array_like
        Fluxes measured at the times in time in arbitrary units.
    flux_err : array_like
        Measurement errors on the fluxes.
    M : float, optional
        Mass of the target star in [solar units]. Default is solar value.
    R : float, optional
        Radius of the target star in [solar units]. Defaults is solar value.
    freq : array_like
        Array of freqencies to search in units of [1/days]. If this is set
        fmin, fmax and OS are ignored.
    fmin : float, optional
        The lower boundary of the frequency range to search in units of [1/days].
        If not set it is determined from the temporal span of the data.
    fmax : float, optional
        The lower boundary of the frequency range to search in units of [1/days].
        If not set it is determined from the Hill radius of the star.
    OS : int
        OverSampling parameter. Sets the number of frequencies searched when 
        fmin and fmax determine the frequency boundaries. Default is 3.
    ES : int 
        EpochSampling parameter. Sets the step in transit epoch as 1/ES times
        the minimum transit duration searched. Default is 3.
    Qmin : int
        Sets the minimum transit duration sampled as 1/Qmin times the expected
        duration. Default is 3.
    Qmax : int
        Sets the maximum transit duration sampled as Qmax times the expected
        duration. Default is 3.
        
    Returns
    -------
    freq : array
        Contains the frequencies searched in units of [1/days]. Is the same
        as the freq parameter when set. Otherwise it depends on M, R, fmin, fmax
        and OS.
    dchisq : array
        The maximum Signal Residue at each frequency.
        
    Notes 
    -----
    For details on the principle of BLS see Kovacs+02
    For details on the chisq formulation of BLS see Collier Cameron+06
    For details on the optimzed parameter choices see Ofir+14 
    """
    time = np.asarray(time)
    flux = np.asarray(flux)
    flux_err = np.asarray(flux_err)
    
    if len(time)!=len(flux):
        raise ValueError('time and flux must have the same first dimension.')
    
    if np.shape(flux)!=np.shape(flux_err):
        raise ValueError('flux and flux_err must be the same shape.')
        
    time = time*SecInDay # [s]
    S = np.max(time)-np.min(time) # [s]  
    weights = 1./flux_err**2.
    
    t = np.sum(weights, axis=0)
    flux = flux - np.sum(flux*weights, axis=0)/t # Subtract average
    chisq0 = np.sum(flux**2.*weights, axis=0) # Best fit constant model.
    
    # First check stellar parameters.
    if options.get("M") is not None:
        M = options.get("M")
    else: M = 1.

    if options.get("R") is not None:
        R = options.get("R")
    else: R = 1.
    
    print 'Using M=%.2f Msun and R=%.2f Rsun for the stellar parameters.'%(M, R)
    
    M *= Msun # [kg]
    R *= Rsun # [m]
    
    # Defaults for remaining keyword arguments.
    fmin = 2./S # Two transits in timespan.
    fmax = 1./(2.*np.pi)*np.sqrt(G*M/(3.*R)**3.) # Smallest orbit determined by the Hill radius.
    OS = 3 # Frequency oversampling.
    ES = 3 # Transit epoch sampling.
    Qmin = 3 # Minimum transit duration sampling.
    Qmax = 3 # Maximum transit duration sampling.
    
    # Frequency sampling.
    if options.get('freq') is not None:
        freq = options.get('freq')/SecInDay #[s]
        nfreq = len(freq)
        print 'Using user specified frequencies, ignoring fmin, fmax, OS.'
    
    else:
        if options.get("fmin") is not None:
            if  options.get("fmin")/SecInDay < fmin:
                print 'Warning: minimum frequency contains <2 events.'
            fmin = options.get("fmin")/SecInDay
                
        if options.get('fmax') is not None:
            if  options.get('fmax')/SecInDay < fmax:
                print 'Warning: maximum frequency is inside the Hill radius.'
            fmax = options.get("fmax")/SecInDay
                
        if fmax<fmin:
            raise AttributeError('fmax must be larger than fmin.')
            
        if options.get('OS') is not None:
            OS = options.get('OS')
        
        print 'Sampling frequencies between fmin=%.2e, fmax=%.2e per day with OS=%i'%(fmin*SecInDay, fmax*SecInDay, OS)
    
        # Compute sample frequencies in range fmin, fmax
        A = ((2.*np.pi)**(2./3.)*R)/(np.pi*(G*M)**(1./3.)*S*OS)
        xmax = int(np.rint((fmax**(1./3.)-fmin**(1./3.))*3./A+1.))
        xvec = np.array(range(1, xmax))
        freq = (A/3.*(xvec-1.)+fmin**(1./3.))**3.
        nfreq = len(freq)

    if options.get('ES') is not None:
        ES = options.get('ES')

    print 'Sampling transit epochs at %.2f time(s) the smallest sampled transit.'%(1./ES)

    if options.get('Qmin') is not None:
        Qmin = options.get('Qmin')
        
    if options.get('Qmax') is not None:
        Qmax = options.get('Qmax')

    NumSteps = Qmin*Qmax
    
    print 'Sampling transits between %.2f and %.2f time(s) the expected duration in %i steps.'%(1./Qmin, Qmax, NumSteps)
    
    # Calculate the expected transit durations.
    q = ((2.*np.pi)**(2./3)*R*freq**(2./3))/(np.pi*(G*M)**(1./3))
    
    # Calculate the number of bins in the sampling grid.
    nbins = np.ceil((Qmin*ES)/q)
    nbins = np.minimum(nbins, len(time))
    
    time = time - time[0]
    
    #dchisq = np.zeros((nfreq,np.shape(flux)[1]))
    depth = np.zeros((nfreq,np.shape(flux)[1]))
    dchisq = np.zeros((nfreq,np.shape(flux)[1]))
    hchisq = np.zeros((nfreq,np.shape(flux)[1]))
    # Loop over frequency domain.
    for i in xrange(nfreq):
        
        # Phase fold the lightcurve. 
        phase = (time*freq[i])%1.
        phase_sort = np.argsort(phase)
        
        phase_fold = np.take(phase, phase_sort, axis=0)
        flux_fold = np.take(flux, phase_sort, axis=0)
        weights_fold = np.take(weights, phase_sort, axis=0)

        # Sum the data on a regular grid.
        bins = np.linspace(0., 1., nbins[i]+1)
        r_cum = cumsum_to_grid(phase_fold, bins, weights_fold)
        s_cum = cumsum_to_grid(phase_fold, bins, flux_fold*weights_fold)
        q_cum = cumsum_to_grid(phase_fold, bins, flux_fold**2*weights_fold)
        
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES+1] + r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES+1] + s_cum[-1], axis=0)
        q_cum = np.append(q_cum, q_cum[1:NumSteps*ES+1] + q_cum[-1], axis=0)
        
        # Calculate the indices of the start and end of transit for all epochs
        # and durations.
        i1, i2 = np.indices((nbins[i], NumSteps))
        i1 += 1
        i2 = i1+(i2+1)*ES
       
        # Compute arrays of r and s values.
        r = r_cum[i2] - r_cum[i1 - 1]
        s = s_cum[i2] - s_cum[i1 - 1]
        q = q_cum[i2] - q_cum[i1 - 1]
        
        # Find the quality of the best fit.
        depth_tmp = s*t/(r*(t - r))
        dchisq_tmp = s*depth_tmp
        hchisq_tmp = chisq0 - s**2/(t - r) - q
        
        for j in range(np.shape(flux)[1]):
            if np.all(np.isnan(dchisq_tmp[:,:,j])): continue
            args = np.nanargmax(dchisq_tmp[:,:,j])

            depth[i,j] = depth_tmp[:,:,j].ravel()[args]
            dchisq[i,j] = dchisq_tmp[:,:,j].ravel()[args]
            hchisq[i,j] = hchisq_tmp[:,:,j].ravel()[args]
    
    return freq*SecInDay, dchisq, depth, hchisq

def BLS(time, flux, flux_err, **options):
    """
    Perform Box Least-Squares on a given lightcurve.

    Returns the maximum Signal Residue as a function of orbital frequency.

    Parameters
    ----------
    time : array_like
        Array contaning the measurement times in units of [days].
    flux : array_like
        Fluxes measured at the times in time in arbitrary units.
    flux_err : array_like
        Measurement errors on the fluxes.
    M : float, optional
        Mass of the target star in [solar units]. Default is solar value.
    R : float, optional
        Radius of the target star in [solar units]. Defaults is solar value.
    freq : array_like
        Array of freqencies to search in units of [1/days]. If this is set
        fmin, fmax and OS are ignored.
    fmin : float, optional
        The lower boundary of the frequency range to search in units of [1/days].
        If not set it is determined from the temporal span of the data.
    fmax : float, optional
        The lower boundary of the frequency range to search in units of [1/days].
        If not set it is determined from the Hill radius of the star.
    OS : int
        OverSampling parameter. Sets the number of frequencies searched when 
        fmin and fmax determine the frequency boundaries. Default is 3.
    ES : int 
        EpochSampling parameter. Sets the step in transit epoch as 1/ES times
        the minimum transit duration searched. Default is 3.
    Qmin : int
        Sets the minimum transit duration sampled as 1/Qmin times the expected
        duration. Default is 3.
    Qmax : int
        Sets the maximum transit duration sampled as Qmax times the expected
        duration. Default is 3.
        
    Returns
    -------
    freq : array
        Contains the frequencies searched in units of [1/days]. Is the same
        as the freq parameter when set. Otherwise it depends on M, R, fmin, fmax
        and OS.
    dchisq : array
        The maximum Signal Residue at each frequency.
        
    Notes 
    -----
    For details on the principle of BLS see Kovacs+02
    For details on the chisq formulation of BLS see Collier Cameron+06
    For details on the optimzed parameter choices see Ofir+14 
    """
    time = np.asarray(time)
    flux = np.asarray(flux)
    flux_err = np.asarray(flux_err)
    
    if len(time)!=len(flux):
        raise ValueError('time and flux must have the same first dimension.')
    
    if np.shape(flux)!=np.shape(flux_err):
        raise ValueError('flux and flux_err must be the same shape.')
        
    time = time*SecInDay # [s]
    S = np.max(time)-np.min(time) # [s]  
    weights = 1./flux_err**2.
    
    t = np.sum(weights, axis=0)
    flux = flux - np.sum(flux*weights, axis=0)/t # Subtract average
    chisq0 = np.sum(flux**2.*weights, axis=0) # Best fit constant model.
    
    # First check stellar parameters.
    if options.get("M") is not None:
        M = options.get("M")
    else: M = 1.

    if options.get("R") is not None:
        R = options.get("R")
    else: R = 1.
    
    print 'Using M=%.2f Msun and R=%.2f Rsun for the stellar parameters.'%(M, R)
    
    M *= Msun # [kg]
    R *= Rsun # [m]
    
    # Defaults for remaining keyword arguments.
    fmin = 2./S # Two transits in timespan.
    fmax = 1./(2.*np.pi)*np.sqrt(G*M/(3.*R)**3.) # Smallest orbit determined by the Hill radius.
    OS = 3 # Frequency oversampling.
    ES = 3 # Transit epoch sampling.
    Qmin = 3 # Minimum transit duration sampling.
    Qmax = 3 # Maximum transit duration sampling.
    
    # Frequency sampling.
    if options.get('freq') is not None:
        freq = options.get('freq')/SecInDay #[s]
        nfreq = len(freq)
        print 'Using user specified frequencies, ignoring fmin, fmax, OS.'
    
    else:
        if options.get("fmin") is not None:
            if  options.get("fmin")/SecInDay < fmin:
                print 'Warning: minimum frequency contains <2 events.'
            fmin = options.get("fmin")/SecInDay
                
        if options.get('fmax') is not None:
            if  options.get('fmax')/SecInDay < fmax:
                print 'Warning: maximum frequency is inside the Hill radius.'
            fmax = options.get("fmax")/SecInDay
                
        if fmax<fmin:
            raise AttributeError('fmax must be larger than fmin.')
            
        if options.get('OS') is not None:
            OS = options.get('OS')
        
        print 'Sampling frequencies between fmin=%.2e, fmax=%.2e per day with OS=%i'%(fmin*SecInDay, fmax*SecInDay, OS)
    
        # Compute sample frequencies in range fmin, fmax
        A = ((2.*np.pi)**(2./3.)*R)/(np.pi*(G*M)**(1./3.)*S*OS)
        xmax = int(np.rint((fmax**(1./3.)-fmin**(1./3.))*3./A+1.))
        xvec = np.array(range(1, xmax))
        freq = (A/3.*(xvec-1.)+fmin**(1./3.))**3.
        nfreq = len(freq)

    if options.get('ES') is not None:
        ES = options.get('ES')

    print 'Sampling transit epochs at %.2f time(s) the smallest sampled transit.'%(1./ES)

    if options.get('Qmin') is not None:
        Qmin = options.get('Qmin')
        
    if options.get('Qmax') is not None:
        Qmax = options.get('Qmax')

    NumSteps = Qmin*Qmax
    
    print 'Sampling transits between %.2f and %.2f time(s) the expected duration in %i steps.'%(1./Qmin, Qmax, NumSteps)
    
    # Calculate the expected transit durations.
    q = ((2.*np.pi)**(2./3)*R*freq**(2./3))/(np.pi*(G*M)**(1./3))
    
    # Calculate the number of bins in the sampling grid.
    nbins = np.ceil((Qmin*ES)/q)
    #nbins = np.minimum(nbins, len(time))

    time = time - time[0]
    
    #dchisq = np.zeros((nfreq,np.shape(flux)[1]))
    depth = np.zeros((nfreq))
    dchisq = np.zeros((nfreq))
    hchisq = np.zeros((nfreq))
    # Loop over frequency domain.
    for i in xrange(nfreq):
        
        # Phase fold the lightcurve. 
        phase = (time*freq[i])%1.
        phase_sort = np.argsort(phase)
        
        phase_fold = np.take(phase, phase_sort, axis=0)
        flux_fold = np.take(flux, phase_sort, axis=0)
        weights_fold = np.take(weights, phase_sort, axis=0)

        # Sum the data on a regular grid.
        bins = np.linspace(0., 1., nbins[i]+1)
        r_cum = cumsum_to_grid(phase_fold, bins, weights_fold)
        s_cum = cumsum_to_grid(phase_fold, bins, flux_fold*weights_fold)
        q_cum = cumsum_to_grid(phase_fold, bins, flux_fold**2*weights_fold)
        
        r_cum = np.append(r_cum, r_cum[1:NumSteps*ES+1] + r_cum[-1], axis=0)
        s_cum = np.append(s_cum, s_cum[1:NumSteps*ES+1] + s_cum[-1], axis=0)
        q_cum = np.append(q_cum, q_cum[1:NumSteps*ES+1] + q_cum[-1], axis=0)
        
        # Calculate the indices of the start and end of transit for all epochs
        # and durations.
        i1, i2 = np.indices((nbins[i], NumSteps))
        i1 += 1
        i2 = i1+(i2+1)*ES
       
        # Compute arrays of r and s values.
        r = r_cum[i2] - r_cum[i1 - 1]
        s = s_cum[i2] - s_cum[i1 - 1]
        q = q_cum[i2] - q_cum[i1 - 1]
        
        # Find the quality of the best fit.
        depth_tmp = s*t/(r*(t - r))
        dchisq_tmp = s*depth_tmp
        hchisq_tmp = chisq0 - s**2/(t - r) - q
        
        args = np.nanargmax(dchisq_tmp.ravel())

        depth[i] = depth_tmp.ravel()[args]
        dchisq[i] = dchisq_tmp.ravel()[args]
        hchisq[i] = hchisq_tmp.ravel()[args]
    
    return freq*SecInDay, dchisq, depth, hchisq
        

def main():
    from astropy.io import fits
    
    f = fits.open('kplr010666592-2009166043257_llc.fits')
    times = f[1].data['TIME']
    flux = f[1].data['PDCSAP_FLUX']
    flux_err = f[1].data['PDCSAP_FLUX_ERR']
    f.close()
    
    here = np.isfinite(flux)
    times = times[here]
    flux = flux[here]
    flux_err = flux_err[here]
    
    med = np.median(flux)
    std = np.std(flux)
    
    here = (flux>med-2*std)&(flux<med+2*std)
    
    coef = np.polyfit(times[here], flux[here], 6)
    fit = np.polyval(coef, times)
    
    plt.title('HaTP-7 b')
    plt.errorbar(times, flux, yerr=flux_err, fmt='.', label='Kepler data')
    plt.plot(times, fit, c='r', label='Polynomial')
    plt.xlabel('Time [BJD-2454833]')
    plt.ylabel('Flux')
    plt.legend()
    plt.tight_layout()
    plt.savefig('HATP7_flux.png', dpi=300)
    plt.show()

    flux = flux/fit
    flux_err = flux_err/fit

    start = time()
    freq, dchisq = BLS(times, flux, flux_err)
    
    plt.title('HaTP-7 b')
    plt.plot(1/freq, dchisq, label='Periodogram')
    plt.xlabel(r'Period [days]')
    plt.ylabel(r'$\Delta\chi^2$')
    plt.axvline(2.2047298, c='k', ls='--', label='Measured period')
    plt.tight_layout()
    plt.legend()
    plt.savefig('HATP7_BLS.png', dpi=300)
    plt.show()
    
    plt.title('HaTP-7 b')
    plt.plot(1/freq, dchisq, label='Periodogram')
    plt.xlabel(r'Period [days]')
    plt.ylabel(r'$\Delta\chi^2$')
    plt.axvline(2.2047298, c='k', ls='--', label='Measured period')
    plt.xlim(1,3.5)
    plt.tight_layout()
    plt.legend()
    plt.savefig('HATP7_BLS_zoom.png', dpi=300)
    plt.show()
    
    arg = np.argmax(dchisq)
    period = 1/freq[arg]
    
    phase = ((times-times[0])/period)%1

    plt.title('HaTP-7 b')
    plt.errorbar(phase, flux, yerr=flux_err, fmt='.')
    plt.xlabel(r'Reduced Phase')
    plt.ylabel('Relative Flux')
    plt.tight_layout()
    plt.savefig('HATP7_folded.png', dpi=300)
    plt.show()
    
    return 0

if __name__ == '__main__':
	main()


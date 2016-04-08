#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

def reduced_lc(filelist, ascc):
    
    colors = {'LPN':'b', 'LPE':'r', 'LPS':'g', 'LPW':'y', 'LPC':'c'}
    
    plt.figure(figsize=(16,6))
    
    ax1 = plt.subplot(211)
    ax1.set_title(ascc)
    ax2 = plt.subplot(212, sharey=ax1)
    
    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            try:
                grp = f['data/'+ascc]
            except:
                continue
            
            lst = grp['lst'].value
            jdmid = grp['jdmid'].value
            mag = grp['mag0'].value
            
        cam = filename[-8:-5]
        
        ax1.plot(jdmid, mag, '.', c=colors[cam])
        ax2.plot(lst, mag, '.', c=colors[cam])
        
    ax1.set_ylabel('Magnitude')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Magnitude')
    ax2.set_ylim(.5, -.5)
    
    plt.tight_layout()
    plt.show()
    
    return

def systematics(filename):
    
    from package import plotting

    systematics = plotting.SysPlot(filename)
    systematics.plot_trans()
    systematics.plot_intrapix()

    return
    
def coverage_nobs(filelist):

    from package.plotting import viridis

    filelist = np.sort(filelist)

    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    nobs = np.array([])

    for filename in filelist:
        
        with h5py.File(filename, 'r') as f:
            
            try:
                grp = f['header_table']
            except:
                pass
            else:
                ascc_ = grp['ascc'].value
                ra_ = grp['ra'].value
                dec_ = grp['dec'].value
                nobs_ = grp['nobs'].value
                
        ascc = np.append(ascc, ascc_)
        ra = np.append(ra, ra_)
        dec = np.append(dec, dec_)
        nobs = np.append(nobs, nobs_)
        
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    ra = ra[args]
    dec = dec[args]
    nobs = np.bincount(idx, nobs)

    plt.figure(figsize=(16, 9))
    plt.subplot(111, projection='mollweide')
    plt.scatter((ra - 180)*np.pi/180, dec*np.pi/180, c=nobs, cmap=viridis, edgecolor='None')
    plt.colorbar()
    plt.tight_layout()
    plt.show()

    return
    
def coverage_baseline(filename):
    
    from package.plotting import viridis
    
    with h5py.File(filename, 'r') as f:
        
        try:
            grp = f['header']
        except:
            pass
        else:
            ascc = grp['ascc'].value
            ra = grp['ra'].value
            dec = grp['dec'].value
            nobs = grp['nobs'].value
            jdmin = grp['jdmin'].value
            jdmax = grp['jdmax'].value

    plt.figure(figsize=(16,9))
    plt.subplot(111, projection='mollweide')
    plt.scatter((ra - 180)*np.pi/180, dec*np.pi/180, c=(jdmax - jdmin), cmap=viridis, edgecolor='None')
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    
    return

def sigma_fraction(filename):
    
    from package import IO
    
    f = IO.SysFile(filename)
    hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()

    sigma = sigma[np.isfinite(sigma)]
    print sigma.size, sum(sigma == 0), sum(sigma == 2)

    return

def clouds(filename):
    
    from package import IO
    from package.plotting import viridis
    from matplotlib import colors
    
    f = IO.SysFile(filename)
    hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()

    select = np.where(~np.isnan(clouds))

    plt.subplot(121)
    plt.hist(clouds[select], bins=np.linspace(-5,5,1001), log=True, histtype='step')
    plt.xlim(-5,5)
    plt.ylim(1,)

    plt.subplot(122)
    plt.hist(sigma[select], bins=np.linspace(0,2,1001), log=True, histtype='step')
    plt.xlim(0,2)
    plt.ylim(1,)

    plt.tight_layout()
    plt.show()

    plt.hist2d(clouds[select], sigma[select], bins=[np.linspace(-5,5,1001), np.linspace(0,2,1001)], norm=colors.LogNorm(), cmap=viridis)
    cb = plt.colorbar()

    plt.xlim(-5,5)
    plt.ylim(0,2)

    cb.set_label('Count')
    plt.xlabel('clouds')
    plt.ylabel('sigma')

    plt.tight_layout()
    plt.show()

    return
    
def astrometry(tablename='LPEastrometry'):
    
    from mascara.reduction import dbquery

    myq = dbquery.MascDB()
    myq.connect()

    cur = myq.dbconnect.cursor()
    cur.execute("SELECT lstseq, alt, az, tho, xo, yo FROM {}".format(tablename))
    data = cur.fetchall()

    myq.disconnect()

    data = np.array(data)
    data = data.T

    plt.subplot(511)
    plt.plot(data[1], '.')

    plt.subplot(512)
    plt.plot(data[2], '.')

    plt.subplot(513)
    plt.plot(data[3], '.')

    plt.subplot(514)
    plt.plot(data[4], '.')

    plt.subplot(515)
    plt.plot(data[5], '.')

    plt.tight_layout()
    plt.show()

    return
    
def std_dec(filename):
    
    import detrend
    
    with h5py.File(filename, 'r') as f:
        
        grp = f['header']
        ascc = grp['ascc'].value
        dec = grp['dec'].value
        Nobs = grp['nobs'].value
        
        grp = f['data']
        
        std0 = np.zeros(len(ascc))
        std = np.zeros(len(ascc))
        
        X = np.array([])
        Y = np.array([])
        M = np.array([])
        
        for i in range(len(ascc)):
            
            lstseq = grp[ascc[i] + '/lstseq'].value
            jdmid = grp[ascc[i] + '/jdmid'].value
            lst = grp[ascc[i] + '/lst'].value
            mag = grp[ascc[i] + '/mag0'].value
            emag = grp[ascc[i] + '/emag0'].value
            nobs = grp[ascc[i] + '/nobs'].value
            x = grp[ascc[i] + '/x'].value
            y = grp[ascc[i] + '/y'].value

            sel = (nobs == 50)
            x = x[sel]
            y = y[sel]
            lstseq = lstseq[sel]
            jdmid = jdmid[sel]
            lst = lst[sel]
            mag = mag[sel]
            emag = emag[sel]/np.sqrt(50.)
            
            if len(jdmid) == 0:
                continue
                
            weights = 1/emag**2
            
            n1 = np.ptp(lstseq) + 1
            n2 = np.ptp(lstseq%270) + 1
            
            n1 = np.maximum(n1, 2)
            n2 = np.maximum(n2, 2)
            
            pars, fit, chisq = detrend.new_harmonic(jdmid, lst, mag, weights, [n1, n2])
            
            std0[i] = np.std(mag)
            std[i] = np.std(mag - fit)
            
            X = np.append(X, x)
            Y = np.append(Y, y)
            M = np.append(M, mag - fit)
    
    xedges = np.linspace(0, 4008, 3*167+1)
    yedges = np.linspace(0, 2672, 2*167+1)
    bins = [xedges, yedges]
    m0, xedges, yedges = np.histogram2d(X, Y, bins=bins)
    m1, xedges, yedges = np.histogram2d(X, Y, bins=bins, weights=M)
    m2, xedges, yedges = np.histogram2d(X, Y, bins=bins, weights=M**2)
    
    std = np.sqrt(m2/m0 - (m1/m0)**2)
    
    plt.subplot(111, aspect='equal')
    plt.pcolormesh(xedges, yedges, std.T, vmin=0, vmax=.2)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.savefig('/home/talens/stdmap_2015Q2LPS.png')
    #plt.show()
    
    #sel = (Nobs > 270)
            
    #plt.plot(dec[sel], std0[sel], '.', c='r', alpha=.5)
    #plt.plot(dec[sel], std[sel], '.', c='k', alpha=.5)
    #plt.ylim(0, .2)
    #plt.xlabel('Declination [deg]')
    #plt.ylabel('STD')
    #plt.show()
    
    return
    
def main(args):
    
    import glob
    
    #filelist = glob.glob('/data3/talens/2015Q2/*/red0_vmag_2015Q2LP?.hdf5')
    #reduced_lc(filelist, '807144')
    
    #filename = '/data3/talens/2015Q2/LPE/sys0_vmag_201506ALPE.hdf5'
    #systematics(filename)
    
    #filelist = glob.glob('/data2/mascara/LaPalma/201506??LPE/fLC/fLC_*.hdf5')
    #coverage_nobs(filelist)
    
    #filename = '/data3/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    #coverage_baseline(filename)
    
    #filename = '/data3/talens/2015Q2/LPE/sys0_vmag_201506ALPE.hdf5'
    #sigma_fraction(filename)
    
    #filename = '/data3/talens/2015Q2/LPE/sys0_vmag_201506ALPE.hdf5'
    #clouds(filename)
    
    filename = '/data3/talens/2015Q2/LPS/red0_vmag_2015Q2LPS.hdf5'
    std_dec(filename)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

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
    
    from pckage.plotting import viridis
    
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

    plt.hist2d(clouds[select], sigma[select], bins=[np.linspace(-5,5,1001), np.linspace(0,2,1001)], norm=colors.LogNorm())
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
    
def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

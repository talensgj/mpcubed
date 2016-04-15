#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table, vstack

def read_catalog(path='/data2/talens/I_280B'):
    
    filelist = glob.glob(os.path.join(path, 'cc*.dat.gz'))
    filelist = np.sort(filelist)
    
    colnames = ['RAhour', 'DEdeg', 'e_RAhour', 'e_DEdeg', 'Plx', 'e_Plx', 'pmRA', 'pmDE', 'e_pmRA', 'e_pmDE', 'Bmag', 'Vmag', 'e_Bmag', 'e_Vmag', 'Scat', 'v1', 'v2', 'v3', 'v4', 'd12', 'd3', 'd4', 'd5', 'd6', 'SpType', 'TYC1', 'TYC2', 'TYC3', 'HIP', 'HD', 'DM', 'ASCC', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag']
    colstarts = [0, 13, 26, 31, 36, 43, 49, 57, 65, 71, 77, 83, 89, 94, 99, 104, 105, 106, 107, 108, 110, 111, 112, 113, 115, 136, 140, 145, 147, 154, 161, 170, 178, 184, 189, 195, 200, 206]
    
    first = True
    for filename in filelist:
        
        print filename
        
        tmp = Table.read(filename, format='ascii.fixed_width_no_header',
                         names=colnames, col_starts=colstarts)
        
        if (tmp['v4'].dtype != '|S1'):
            print 'Fixing datatype.'
            v4 = tmp['v4'].astype('|S1')
            tmp.replace_column('v4', v4)
        
        if first:
            data = tmp
            first = False
        else:
            data = vstack([data, tmp])
    
    with h5py.File(os.path.join(path, 'ASCC.hdf5')) as f:
        
        grp = f.create_group('data')
        for name in colnames:
            grp.create_dataset(name, data=data[name])
    
    return data
   
def stellar_radius(VMAG, R, Mv):
    
    Rstar = R*10**((Mv - VMAG)/5.)
    
    return Rstar
    
def parallax(catalogue='/data2/talens/I_280B/ASCC.hdf5'):
    
    type_dict = {'O5':(13.4, -5.1), 'O6':(12.2, -5.1), 'O7':(11., -4.9), 'O8':(10., -4.6),
                 'B0':(6.7, -3.4), 'B1':(5.2, -2.6), 'B2':(4.1, -1.6), 'B3':(3.8, -1.3), 'B5':(3.2, -.5), 'B6':(2.9, -.1), 'B7':(2.7, .3), 'B8':(2.5, .6), 'B9':(2.3, .8),
                 'A0':(2.2, 1.1), 'A1':(2.1, 1.3), 'A2':(2., 1.5), 'A5':(1.8, 2.2), 'A8':(1.5, 2.7),
                 'F0':(1.4, 3.), 'F2':(1.3, 3.4), 'F5':(1.2, 3.9), 'F8':(1.1, 4.3),
                 'G0':(1.06, 4.7), 'G2':(1.03, 4.9), 'G8':(.96, 5.6),
                 'K0':(.93, 5.7), 'K1':(.91, 6.), 'K3':(.86, 6.5), 'K4':(.83, 6.7), 'K5':(.8, 7.1), 'K7':(.74, 7.8),
                 'M0':(.63, 8.9), 'M1':(.56, 9.6), 'M2':(.48, 10.4), 'M3':(.41, 11.1), 'M4':(.35, 11.9), 'M5':(.29, 12.8), 'M6':(.24, 13.8), 'M7':(.20, 14.7)}
    
    with h5py.File(catalogue, 'r') as f:
        
        grp = f['data']
        ascc = grp['ASCC'].value
        plx = grp['Plx'].value
        vmag = grp['Vmag'].value
        sptype = grp['SpType'].value
        
    sel = (plx != 999999) & (plx > 0)
    
    plx = plx*1e-2
    vmag = vmag*1e-3
        
    d = np.where(sel, 1/(plx*1e-3), -1) # [pc]
    VMAG = np.where(sel, vmag - 5*np.log10(d/10.), 99.999)
    
    Rstar = np.zeros(len(ascc))-1
    for i in range(len(ascc)):
        
        stype = sptype[i][:2]
        
        try:
            pars = type_dict[stype]
        except:
            print stype, d[i]
        else:
            Rstar[i] = stellar_radius(VMAG[i], *pars)
            print stype, d[i], Rstar[i]
    
    with h5py.File('ASCC_radii.hdf5') as f:
        
        grp = f.create_group('data')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('d', data=d)
        grp.create_dataset('rstar', data=Rstar)
    
    return
    
def main(args):
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

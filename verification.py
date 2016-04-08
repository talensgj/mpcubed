#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

def fLC_lstseq(filename):
    
    with h5py.File(filename, 'r') as f:
        ascc = f['header_table/ascc'].value
        
        for i in range(0, len(ascc), 50):
            lc = f['data/' + ascc[i]].value
            
            plt.plot(lc['jdmid'], lc['lstseq'])
            plt.show()
            plt.close()
    
    return 

def fLC_nobs(filelist):
    
    filelist = np.sort(filelist)
    
    for filename in filelist:
        
        try: 
            with h5py.File(filename, 'r') as f:
                grp = f['header_table']
                
        except:
            print 'Broken file', filename
            continue
        
        with h5py.File(filename, 'r') as f:
        
            ascc = f['header_table/ascc'].value
            nobs = f['header_table/nobs'].value
            
            for i in range(len(ascc)):
                
                try:
                    lc = f['data/'+ascc[i]].value
                except:
                    print 'Star not in data/', filename
                    continue
                
                if len(lc) != nobs[i]:
                    print filename
                    break
    
    return

def quarter(filename):
    
    with h5py.File(filename, 'r') as f:

        grp = f['header']
        ascc = grp['ascc'].value
        nobs = grp['nobs'].value
        jdmin = grp['jdmin'].value
        jdmax = grp['jdmax'].value
        lstseqmin = grp['lstseqmin'].value
        lstseqmax = grp['lstseqmax'].value

        print nobs.dtype
        print jdmin.dtype
        print lstseqmin.dtype

        grp = f['data']
        
        ascc1 = grp.keys()
        
        print len(ascc), len(ascc1)
        print np.setdiff1d(ascc, ascc1)
        print np.setdiff1d(ascc1, ascc)
        
        for i in range(len(ascc)):
            
            jdmid = grp[ascc[i] + '/jdmid'].value
            mag = grp[ascc[i] + '/mag0'].value
            lstseq = grp[ascc[i] + '/lstseq'].value
            
            if (len(jdmid) != nobs[i]):
                print 'Bad nobs.'
                
            if (np.amin(jdmid) != jdmin[i]):
                print 'Bad jdmin.'
                
            if (np.amax(jdmid) != jdmax[i]):
                print 'Bad jdmax.'
                
            if (np.amin(lstseq) != lstseqmin[i]):
                print 'Bad lstseqmin.'
                
            if (np.amax(lstseq) != lstseqmax[i]):
                print 'Bad lstseqmax.'

    return

def main(args):
    
    import glob
    
    filename = '/data3/talens/2015Q3/LPW/fLC_201507BLPW.hdf5'
    fLC_lstseq(filename)
    
    filelist = glob.glob('/data2/mascara/LaPalma/2015060[1-5]LPE/fLC/fLC_*.hdf5')
    fLC_nobs(filelist)
    
    filename = '/data3/talens/2015Q3/LPW/red0_vmag_2015Q3LPW.hdf5'
    quarter(filename)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

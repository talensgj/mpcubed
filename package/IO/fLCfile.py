#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

class fLCfile():
    
    def __init__(self, fLCfile):
        """
        Read data from a standard fLC file.
        """
        
        self.fLCfile = fLCfile
    
        return
        
    def read_global(self):
        """
        Read the data from the global group.
        """
        
        with h5py.File(self.fLCfile, 'r') as f:
            data = f['global'].attrs.items()
        
        return data
        
    def read_header(self, fields, ascc=None):
        """
        Read data from the header.
        """
        
        if ascc is None:
            
            # If no ascc are given read data for all stars.
            with h5py.File(self.fLCfile, 'r') as f:
                data = [f['header_table/'+field].value for field in fields]
            
        else:
            nstars = len(ascc)
            nfields = len(fields)
            
            # Create arrays.
            data = [np.zeros(nstars) for _ in range(nfields)]
            
            # Read the data for the specified stars.
            with h5py.File(self.fLCfile, 'r') as f:
                
                all_ascc = f['header_table/ascc'].value
                here = np.in1d(all_ascc, ascc)
                
                for j in range(nfields):
                    data[j] = f['header_table/'+fields[j]].value[here]
                    
        return data
    
    def read_data(self, fields, ascc=None, nobs=None):
        """
        Read data from the lightcurves of multiple stars.
        """
        
        if ascc is None:
            # If no stars are specified read stars and nobs from the header.
            ascc, nobs = self.read_header(['ascc', 'nobs'])
        elif nobs is None:
            # If stars are specified but no nobs read nobs from the header.
            nobs, = self.read_header(['nobs'], ascc)
            
        nstars = len(ascc)
        nfields = len(fields)
        npoints = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        # Create arrays.
        data = [np.zeros(npoints) for _ in range(nfields)]
        
        # Read the data.
        with h5py.File(self.fLCfile, 'r') as f:
            for i in range(nstars):
                lc = f['data/'+ascc[i]].value
                for j in range(nfields):
                    data[j][select[i]:select[i+1]] = lc[fields[j]]
            
        return data
           
    def read_star(self, ascc):
        """
        Read the full lightcurve of a particular star.
        """
        
        with h5py.File(self.fLCfile, 'r') as f:
            data = f['data/' + ascc].value
        
        return data
        
if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Outputs the global group for a given fLC file.')
    parser.add_argument('fLCfile', help = 'The input fLC file.')
    args = parser.parse_args()
    
    f = fLCfile(args.fLCfile)
    data = f.read_global()

    for key, value in data:
        print key, value
        
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import red_merge_dev

filelist = ['/data2/talens/2015Q2/LPE/test_tmp_06A.hdf5']*2
outfile = '/data2/talens/2015Q2/LPE/test_qfile.hdf5'

red_merge_dev.make_quarterfile(filelist, outfile)




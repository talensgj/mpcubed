#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "mpcubed",
    version = "3.0.4",
    author = "Geert Jan Talens",
    author_email = "talens@strw.leidenuniv.nl",
    description = "The MASCARA Post-Processing Pipeline.",
    license = "CCC",
    packages = find_packages(),
    long_description = read('README'),
    entry_points={'console_scripts': 
        ['merge_files = mpcubed.io:main',
         'run_calibration = mpcubed.calibration.cdecor:main',
         'figs_calibration = mpcubed.calibration.figures:main',
         'run_boxlstsq = mpcubed.detection.boxlstsq:main',
         'figs_boxlstsq = mpcubed.detection.figures:main']}
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "mpcubed",
    version = "1.0.16",
    author = "Geert Jan Talens",
    author_email = "talens@strw.leidenuniv.nl",
    description = "The MASCARA Post-Processing Pipeline.",
    license = "CCC",
    packages = find_packages(),
    long_description = read('README'),
)

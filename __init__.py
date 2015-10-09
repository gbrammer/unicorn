#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Created by Gabriel Brammer on 2011-05-18.

"""
__version__ = "1.0"

import os
from socket import gethostname as hostname

if hostname().startswith('uni') or hostname().startswith('hyp'):
    GRISM_HOME = '/3DHST/Spectra/Work/'
else:
    GRISM_HOME = '/research/HST/GRISM/3DHST/'

if hostname().startswith('850dhcp8'):
    GRISM_HOME = '/3DHST/Spectra/Work/'
    #threedhst.sex.RUN_MODE='direct'

if hostname().lower().startswith('gabriel') | hostname().startswith('vpn-brammer'):
    GRISM_HOME = '/Users/brammer/3DHST/Spectra/Work/'

GRISM_HOME = os.getenv('THREEDHST')

if GRISM_HOME[-1] != '/':
    GRISM_HOME += '/'
        
import threedhst

try:
    import utils_c #as utils_c
except:
    print """Couldn't import "utils_c" """

import plotting
import prepare
import reduce
import candels
import analysis
import go_3dhst
import galfit
import catalogs

import survey_paper
import go_acs
import fast

import interlace_fit
import intersim

import prepare_acs

noNewLine = '\x1b[1A\x1b[1M'

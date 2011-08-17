#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Created by Gabriel Brammer on 2011-05-18.

$URL$
$Author$
$Date$

"""
__version__ = "$Rev$"

import threedhst

import prepare
import reduce
import candels
import analysis
import go_3dhst
import galfit
import catalogs
import survey_paper
import go_acs

from socket import gethostname as hostname

if hostname().startswith('uni'):
    GRISM_HOME = '/3DHST/Spectra/Work/'
else:
    GRISM_HOME = '/research/HST/GRISM/3DHST/'

if hostname().startswith('850dhcp8'):
    GRISM_HOME = '/3DHST/Spectra/Work/'
    #threedhst.sex.RUN_MODE='direct'
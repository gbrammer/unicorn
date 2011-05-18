#!/usr/bin/env python
# encoding: utf-8
"""
prepare.py

Scripts to "prepare" the FLT files for the 3D-HST reduction:

Direct images:
    1) Align to WCS reference
    2) Subtract sky backgrounds
    3) Run multidrizzle to flag cosmic rays
    4) Make mosaics combining pointings for a given field
    
Grism images:
    1) Copy shifts from direct images
    2) Divide by the F140W flat
    3) Subtract sky background
    4) Run multidrizzle to flag cosmic rays
    5) ...
    
$ URL: https://subversion.assembla.com/svn/threedhst_internal/trunk/prepare.py $
$ Author: gbrammer $
$ Date: 2011-05-18 00:37:27 -0500 (Wed, 18 May 2011) $

"""
__version__ = " $ Rev: 2 $"

def COSMOS():
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files

    os.chdir('/3DHST/Spectra/Work/COSMOS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    gris = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True)
    
    
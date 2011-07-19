"""

Run the full reduction scripts on the ACS parallels

"""
import threedhst
import threedhst.prep_flt_files

import unicorn
import glob
import os
import shutil

import pyraf
from pyraf import iraf

from threedhst.prep_flt_files import process_3dhst_pair as pair

def test():
    """
    COSMOS-23 ACS overlaps with COSMOS-25 WFC3
    """
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/PREP_FLT')
    
    root='jbhm51'
    files=glob.glob('../RAW/%s*fits*' %(root))
    for file in files:
        os.system('cp %s .' %(file))
        os.system('gunzip %s' %(os.path.basename(file)))
        print file
    
    #### Destripe + CTE correction
    threedhst.prep_flt_files.prep_acs(force=True)
    
    os.system('rm *flt.fits')
    
    asn_direct_file = 'jbhm51010_asn.fits'
    asn_grism_file = 'jbhm51020_asn.fits'
    
    #### Copy corrected FLT files to . 
    asn = threedhst.utils.ASNFile(asn_direct_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    #
    asn = threedhst.utils.ASNFile(asn_grism_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    
    #### Flag CRs, subtract background and align WCS    
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    ALIGN = '/3DHST/Spectra/Work/COSMOS/MOSAIC/COSMOS-F140w_11-05-23_sci.fits'
    
    threedhst.shifts.run_tweakshifts(asn_direct_file, verbose=True)

    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.03, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)

    threedhst.shifts.refine_shifts(ROOT_DIRECT=asn_direct_file.split('_as')[0].upper(),
              ALIGN_IMAGE=ALIGN, 
              ALIGN_EXTENSION = 0,
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.03, pixfrac=1, driz_cr=False,
        updatewcs=False, clean=True, median=False)
        
    ### Grism
    threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)

    threedhst.prep_flt_files.startMultidrizzle(asn_grism_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.03, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)
    
    ## Check
    threedhst.gmap.makeImageMap(['JBHM51010_drz.fits[1]*4', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-25-F140W_drz.fits', 'JBHM51020_drz.fits[1]', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-25-G141_drz.fits'], aper_list=[15])
    
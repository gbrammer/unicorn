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
    
    threedhst.prep_flt_files.prep_acs(force=True)
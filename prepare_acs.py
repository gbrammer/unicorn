#!/usr/bin/env python
# encoding: utf-8
"""
prepare_acs.py

Created by Ivelina Momcheva on 2012-08-06.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import numpy as np
import glob
import threedhst
import unicorn

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

def GOODSS_acs (FORCE=False):
    import unicorn
    import threedhst
    from unicorn.go_acs import process_acs_pair as pair
    import glob
    import os

    direct=glob.glob('jbhj*10_asn.fits')
    grism = glob.glob('jbhj*20_asn.fits')
    
    direct.remove('jbhj09010_asn.fits')
    direct.remove('jbhj12010_asn.fits')
    direct.remove('jbhj17010_asn.fits')
    direct.remove('jbhj25010_asn.fits')
    direct.remove('jbhj30010_asn.fits')    
    grism.remove('jbhj09020_asn.fits')
    grism.remove('jbhj12020_asn.fits')
    grism.remove('jbhj17020_asn.fits')
    grism.remove('jbhj25020_asn.fits')
    grism.remove('jbhj30020_asn.fits')    

    
    FORCE=True
    
    loop_list = range(len(direct))
    for i in loop_list:
        root=direct[i].split('_asn.fits')[0][0:6]
    
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/GOODS-S/PREP_FLT')
    
        files=glob.glob('../RAW/%s*flt.fits*' %(root))
        for file in files:
            os.system('cp %s .' %(file))
            os.system('gunzip %s' %(os.path.basename(file)))
            print file

        threedhst.prep_flt_files.prep_acs(root=root,force=True)
    
        #os.system('rm %s*flt.fits' %(root))

        ALIGN_FILES = ['/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits', '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-2_F160W_v1_sci.fits', '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-1_F160W_v1_sci.fits']
    
    
        #### Main preparation loop
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False, field='GOODS-S')
        ALIGN = ALIGN_FILES[0]
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], field = 'GOODS-S', ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, align_geometry='shift')

    #return False        
    
    #make mosaic
    
    #import threedhst.prep_flt_files
    #SCALE=0.05
    #PIXFRAC=1
    #NX, NY, ra, dec = 22000,22000, 53.098403, -27.807325

    #direct_files = glob.glob('GOODS-S-*-F814W_asn.fits')
    #threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-F814W',
    #                  path_to_FLT='./', run_multidrizzle=False)

    #threedhst.prep_flt_files.startMultidrizzle('GOODS-S-F814W_asn.fits',
    #         use_shiftfile=True, skysub=False,
    #         final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
    #         updatewcs=False, clean=True, median=False,
    #         ra=ra, dec=dec,
    #         final_outnx = NX, final_outny=NY)

    #extract spectra

    direct = glob.glob('GOODS-S-0*F814W_asn.fits')
    
    for i in range(len(direct)):
        root = direct[i].split('-F814W_asn.fits')[0]
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/GOODS-S/DATA/')
        #master_catalog = '../GOODS-S_fake.cat'
        #master_segmentation = '../GOODS-S_fake.seg.fits'
        master_catalog = '/3DHST/Photometry/Work/GOODS-S/v2/sextr/catalogs/GOODS-S_F814W.cat'
        master_segmentation = '/3DHST/Photometry/Work/GOODS-S/v2/sextr/checkimages/GOODS-S_F814W.seg.fits'    
        unicorn.go_acs.make_external_catalog(root =root,master_segmentation=master_segmentation, \
            master_catalog=master_catalog, reference_image='../PREP_FLT/'+root+'-F814W_drz.fits')
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/GOODS-S/')
        #hyperion: WFC3_DIR = '/Users/ivastar/3DHST/Data.2.0/GOODS-S/'
        WFC3_DIR = '/3DHST/Spectra/Release/v2.0/GOODS-S/'
        unicorn.go_acs.reduce_acs(root=root, LIMITING_MAGNITUDE=20.,match_wfc3 = True, WFC3_DIR=WFC3_DIR)


def COSMOS_acs (FORCE=False):
    import unicorn
    import threedhst
    from unicorn.go_acs import process_acs_pair as pair
    import glob
    import os

    direct=glob.glob('jbhm5*10_asn.fits')
    grism = glob.glob('jbhm5*20_asn.fits')
    
    #direct.remove('jbhm29010_asn.fits')
    direct.remove('jbhm35010_asn.fits')
    direct.remove('jbhm36010_asn.fits')
    direct.remove('jbhm39010_asn.fits')
    #grism.remove('jbhm29020_asn.fits')
    grism.remove('jbhm35020_asn.fits')
    grism.remove('jbhm36020_asn.fits')
    grism.remove('jbhm39020_asn.fits')
    
    FORCE=True
    
    loop_list = range(len(direct))
    for i in loop_list:
        root=direct[i].split('_asn.fits')[0][0:6]
    
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/PREP_FLT')
    
        files=glob.glob('../RAW/%s*flt.fits*' %(root))
        for file in files:
            os.system('cp %s .' %(file))
            os.system('gunzip %s' %(os.path.basename(file)))
            print file

        threedhst.prep_flt_files.prep_acs(root=root,force=True)
    
        #os.system('rm %s*flt.fits' %(root))

        ALIGN = '/3DHST/Ancillary/COSMOS/CANDELS/UCSC/cosmos_sect*_wfc3ir_F160W_wfc3ir_drz_sci.fits'

    
    
        #### Main preparation loop
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False, field='COSMOS')
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], field='COSMOS',ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, align_geometry='shift')

    #return False        
    
    #make mosaic
    
    #import threedhst.prep_flt_files
    #SCALE=0.05
    #PIXFRAC=1
    #NX, NY, ra, dec = 10700,22000, 150.12634, 2.3336697

    #direct_files = glob.glob('GOODS-S-*-F814W_asn.fits')
    #threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-F814W',
    #                  path_to_FLT='./', run_multidrizzle=False)

    #threedhst.prep_flt_files.startMultidrizzle('COSMOS-F814W_asn.fits',
    #         use_shiftfile=True, skysub=False,
    #         final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
    #         updatewcs=False, clean=True, median=False,
    #         ra=ra, dec=dec,
    #         final_outnx = NX, final_outny=NY)


    
    
    



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


    #make reference catalog
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/GOODS-S/ALIGN')
    plt.clf()
    f606_lim = 24.5
    f160_lim = 22.
    f606 = threedhst.sex.mySexCat('GOODS-S_acsv_conv.cat')
    f160 = threedhst.sex.mySexCat('GOODS-S_F160W.cat')
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'b.')
    plt.xlim([10,40])
    plt.ylim([10,40])
    plt.xlabel('F160W')
    plt.ylabel('F606W')    
    plt.plot([f160_lim,f160_lim],[0,50],color='red')
    plt.plot([0,50],[f606_lim,f606_lim],color='red')
    MAG_F606W = f606['MAG_AUTO']
    MAG_F606W[(f160['MAG_AUTO']<f160_lim) & (f606['MAG_AUTO'] >f606_lim)] = f606_lim-0.01
    f606.renameColumn(original='MAG_AUTO', new='MAG_F606W', verbose=True)
    f606.addColumn(data=MAG_F606W, format='%f', name='MAG_AUTO',verbose=True)
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'r.')
    f606.write(outfile='GOODS-S_F606W_ref.cat')
 
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/GOODS-S/PREP_FLT')    
    direct = glob.glob('GOODS-S-*F814W_asn.fits')
    direct.remove('GOODS-S-09-F814W_asn.fits')
    direct.remove('GOODS-S-12-F814W_asn.fits')
    direct.remove('GOODS-S-25-F814W_asn.fits')
    
    for i in range(len(direct)):
        root = direct[i].split('-F814W_asn.fits')[0]
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/GOODS-S/DATA/')
        master_catalog = '/3DHST/Spectra/Work/ACS_PARALLEL/GOODS-S/ALIGN/GOODS-S_F606W_ref.cat'
        master_segmentation = '/3DHST/Spectra/Work/ACS_PARALLEL/GOODS-S/ALIGN/GOODS-S_F125W_F140W_F160W.seg.fits'    
        unicorn.go_acs.make_external_catalog(root =root,master_segmentation=master_segmentation, \
            master_catalog=master_catalog, reference_image='../PREP_FLT/'+root+'-F814W_drz.fits')
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/GOODS-S/')
        #WFC3_DIR = '/3DHST/Spectra/Release/v2.0/GOODS-S/'
        unicorn.go_acs.reduce_acs(root=root, LIMITING_MAGNITUDE=24.5,match_wfc3 = False, WFC3_DIR='')


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

    # make reference catalog
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/COSMOS/ALIGN')
    plt.clf()
    f606_lim = 24.5
    f160_lim = 22.
    f606 = threedhst.sex.mySexCat('F606W.cat')
    f160 = threedhst.sex.mySexCat('F160W.cat')
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'b.')
    plt.xlim([10,40])
    plt.ylim([10,40])
    plt.xlabel('F160W')
    plt.ylabel(['F606W'])    
    plt.plot([f160_lim,f160_lim],[0,50],color='red')
    plt.plot([0,50],[f606_lim,f606_lim],color='red')
    MAG_F606W = f606['MAG_AUTO']
    MAG_F606W[(f160['MAG_AUTO']<f160_lim) & (f606['MAG_AUTO'] >f606_lim)] = f606_lim-0.01
    f606.renameColumn(original='MAG_AUTO', new='MAG_F606W', verbose=True)
    f606.addColumn(data=MAG_F606W, format='%f', name='MAG_AUTO',verbose=True)
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'r.')
    f606.write(outfile='COSMOS_F606W_ref.cat')
        
        #extract spectra
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/COSMOS/PREP_FLT')
    direct = glob.glob('COSMOS-*F814W_asn.fits')
    #direct.remove('COSMOS-29-F814W_asn.fits')
    #direct.remove('COSMOS-35-F814W_asn.fits')
    #direct.remove('COSMOS-36-F814W_asn.fits')
    #direct.remove('COSMOS-39-F814W_asn.fits')
    
    direct = ['COSMOS-54-F814W_asn.fits', 'COSMOS-55-F814W_asn.fits', 'COSMOS-56-F814W_asn.fits']
    for i in range(len(direct)):
        root = direct[i].split('-F814W_asn.fits')[0]
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/DATA/')
        master_catalog = '/3DHST/Spectra/Work/ACS_PARALLEL/COSMOS/ALIGN/COSMOS_F606W_ref.cat'
        master_segmentation = '/3DHST/Spectra/Work/ACS_PARALLEL/COSMOS/ALIGN/F160W_seg.fits'    
        unicorn.go_acs.make_external_catalog(root =root,master_segmentation=master_segmentation, \
            master_catalog=master_catalog, reference_image='../PREP_FLT/'+root+'-F814W_drz.fits')
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/')
        #WFC3_DIR = '/3DHST/Spectra/Release/v2.0/GOODS-S/'
        unicorn.go_acs.reduce_acs(root=root, LIMITING_MAGNITUDE=24.5,match_wfc3 = False, WFC3_DIR='')
        
    
        

def UDS_acs():
    
    #make reference catalog
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/UDS/ALIGN')
    plt.clf()
    f606_lim = 24.5
    f160_lim = 22.
    f606 = threedhst.sex.mySexCat('UDS_F606W_conv.cat')
    f160 = threedhst.sex.mySexCat('UDS_F160W.cat')
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'b.')
    plt.xlim([10,40])
    plt.ylim([10,40])
    plt.xlabel('F160W')
    plt.ylabel(['F606W'])    
    plt.plot([f160_lim,f160_lim],[0,50],color='red')
    plt.plot([0,50],[f606_lim,f606_lim],color='red')
    MAG_F606W = f606['MAG_AUTO']
    MAG_F606W[(f160['MAG_AUTO']<f160_lim) & (f606['MAG_AUTO'] >f606_lim)] = f606_lim-0.01
    f606.renameColumn(original='MAG_AUTO', new='MAG_F606W', verbose=True)
    f606.addColumn(data=MAG_F606W, format='%f', name='MAG_AUTO',verbose=True)
    plt.plot(f160['MAG_AUTO'],f606['MAG_AUTO'],'r.')
    f606.write(outfile='UDS_F606W_ref.cat')
    
    
    #extract spectra
    #on hyperion
    
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/UDS/PREP_FLT')
    direct = glob.glob('UDS-*F814W_asn.fits')
    direct.remove('UDS-15-F814W_asn.fits')
    direct.remove('UDS-18-F814W_asn.fits')
    direct.remove('UDS-20-F814W_asn.fits')
    direct.remove('UDS-21-F814W_asn.fits')
    direct.remove('UDS-22-F814W_asn.fits')
    
    for i in range(len(direct)):
        root = direct[i].split('-F814W_asn.fits')[0]
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/UDS/DATA/')
        master_catalog = '/3DHST/Spectra/Work/ACS_PARALLEL/UDS/ALIGN/UDS_F606W_ref.cat'
        master_segmentation = '/3DHST/Spectra/Work/ACS_PARALLEL/UDS/ALIGN/UDS_F125W_F140W_F160W.seg.fits'    
        unicorn.go_acs.make_external_catalog(root =root,master_segmentation=master_segmentation, \
            master_catalog=master_catalog, reference_image='../PREP_FLT/'+root+'-F814W_drz.fits')
        os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/UDS/')
        #WFC3_DIR = '/3DHST/Spectra/Release/v2.0/GOODS-S/'
        unicorn.go_acs.reduce_acs(root=root, LIMITING_MAGNITUDE=24.5,match_wfc3 = False, WFC3_DIR='')
        
    
    
    
    



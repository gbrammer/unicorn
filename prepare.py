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
    
$URL$
$Author$
$Date$

"""
__version__ = " $Rev$"

import threedhst
import unicorn

import numpy as np
import os

def fix_GOODSN_asn():
    """
    Replace the bad flt files with the new fixed ones.
    """
    import threedhst.catIO as catIO
    
    list = catIO.Readfile('files.info')
    
    bad = ['ib3701s4q_flt.fits','ib3701skq_flt.fits','ib3702u8q_flt.fits','ib3702uoq_flt.fits','ib3703uzq_flt.fits','ib3703vfq_flt.fits','ib3703vmq_flt.fits','ib3704wpq_flt.fits','ib3704x8q_flt.fits','ib3705xzq_flt.fits','ib3705y1q_flt.fits','ib3705ylq_flt.fits','ib3706b2q_flt.fits','ib3706biq_flt.fits','ib3706bpq_flt.fits','ib3707c8q_flt.fits','ib3707cqq_flt.fits','ib3708i5q_flt.fits','ib3708ipq_flt.fits','ib3709j3q_flt.fits','ib3709joq_flt.fits','ib3728d9q_flt.fits']
    
    asn = threedhst.utils.ASNFile(threedhst.utils.find_fits_gz('../RAW/ib3701050_asn.fits'))
    
    targets = np.unique(list.targname)
    for target in targets:
        root='GOODS-N-'+target.split('GRISM')[1]
                
        direct_x = list.file[(list.targname == target) & (list.filter == 'F140W')]
        direct = []
        for file in direct_x:
            direct.append(file.split('_flt')[0])
                
        grism_x = list.file[(list.targname == target) & (list.filter == 'G141')]
        grism = []
        for file in grism_x:
            if not (file.replace('.gz','') in bad):
                grism.append(file.split('_flt')[0])
        #
        if root == 'GOODS-N-23':
            ### Use only the new grism images because at different angle
            while grism[0].startswith('ib3707'):
                grism.pop(0)
        
        print root, len(direct), len(grism)
                    
        asn.exposures = direct
        asn.product = root+'-F140W'
        asn.write(root+'-F140W_asn.fits')
        
        asn.exposures = grism
        asn.product = root+'-G141'
        asn.write(root+'-G141_asn.fits')
                
        
        
def GOODSN(FORCE=False, GET_SHIFT=True):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
    
    os.chdir(unicorn.GRISM_HOME+'PREP_FLT')
    
    #### Use ACS alignment images
    ALIGN = '/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_nz*drz*fits'
    
    #### Process direct + grism pairs
    direct=glob.glob('GOODS-N-*-F140W_asn.fits')
    grism = glob.glob('GOODS-N-*-G141_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=False)
    
    # threedhst.gmap.makeImageMap(['GOODS-N-42-F140W_drz.fits', 'GOODS-N-42-G141_drz.fits*3', 'GOODS-N-42-F140W_align.fits[0]*4'], aper_list=[15], polyregions=["GOODS-N-42-F140W_asn.pointing.reg"])
    
    #### Direct+grism mosaics
    direct_files = glob.glob('GOODS-N-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-N-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    direct_files = glob.glob('GOODS-N-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-N-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.128254
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False) #,
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False)

    #### Make full alignment mosaic
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN), matchImage='GOODS-N-F140W_drz.fits', match_extension=1, output='GOODS-N-F850LP.fits')
    
    threedhst.gmap.makeImageMap(['GOODS-N-F140W_drz.fits', 'GOODS-N-F850LP.fits[0]*6','GOODS-N-G141_drz.fits*2'], aper_list=[13,14,15], polyregions=glob.glob('*F140W_asn.pointing.reg'), zmin=-0.1, zmax=3)
    
    
def COSMOS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
    
    os.chdir(unicorn.GRISM_HOME+'COSMOS/NEW_PREP')
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    #ALIGN = '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits'
    #ALIGN = '/3DHST/Ancillary/COSMOS/NMBS/COSMOS-1_K_sci.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True, GET_SHIFT=False, SKIP_DIRECT=False)
        threedhst.gmap.makeImageMap(['COSMOS-6-F140W_drz.fits','/tmp/junk.fits','COSMOS-6-G141_drz.fits','/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'][0:2], aper_list=[14], polyregions=["COSMOS-6-F140W_asn.pointing.reg"])
    
    #### Direct mosaic
    direct_files = glob.glob('COSMOS-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Direct mosaic
    direct_files = glob.glob('COSMOS-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.128254
    #SCALE = 0.5
    PIXFRAC=0.6
    NX, NY = int(9670*0.06/SCALE), int(18890*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12634, dec=2.3336697,
             final_outnx = NX, final_outny=NY)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12634, dec=2.3336697,
             final_outnx = NX, final_outny=NY)
    
    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits', 'COSMOS-G141_drz.fits'], aper_list=[13,14], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))
    
    ### Temporary release of the direct mosaic
    # iraf.imcopy('COSMOS-F140W_drz.fits[1]','COSMOS-F140w_11-05-23_sci.fits')
    # iraf.imcopy('COSMOS-F140W_drz.fits[2]','COSMOS-F140w_11-05-23_wht.fits')
    # !tar czvf COSMOS-F140w_11-05-23.tar.gz COSMOS-*-F140W_shifts.txt COSMOS-*-F140W_tweak.fits COSMOS-*-F140W_asn.fits COSMOS-F140W_shifts.txt COSMOS-F140W_asn.fits
    # !mv COSMOS-F140w_11-05-23* ../MOSAIC
        
    threedhst.gmap.makeImageMap(['../MOSAIC/COSMOS-F140w_11-05-23_sci.fits[0]', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/Chandra/CC0570_img.fits[0]*1.5', '/3DHST/Ancillary/COSMOS/Spitzer/mips_24_GO3_sci_10.fits[0]*400', '/3DHST/Ancillary/COSMOS/VLA/vla_20cm_dp_sin_10.fits[0]*40000'], aper_list=[13,14,15,16,17], tileroot=['F140W','F814W','WIRDS-K','0.5-7keV', 'MIPS-24', 'VLA-20cm'])
    
    ### don't need high-res tiles of lo-res images sitting around
    os.system('remove ~/Sites/FITS/tiles/MIPS*1[67].png')
    os.system('remove ~/Sites/FITS/tiles/VLA*1[7].png')
    os.system('remove ~/Sites/FITS/tiles/0.5*17.png')
    
    
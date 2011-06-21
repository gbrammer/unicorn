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
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('GOODS-N-[0-9]*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='GOODS-N-*[0-9]-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)
    
  
def COSMOS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
    
    os.chdir(unicorn.GRISM_HOME+'COSMOS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    #ALIGN = '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits'
    #ALIGN = '../NMBS/COSMOS-1.v4.K_nosky.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    
    threedhst.gmap.makeImageMap(['COSMOS-11-F140W_drz.fits', 'COSMOS-11-G141_drz.fits','/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'][0:], aper_list=[16], polyregions=glob.glob("COSMOS-*-F140W_asn.pointing.reg"))
            
    threedhst.gmap.makeImageMap(['COSMOS-28-F140W_drz.fits','/tmp/junk.fits', 'COSMOS-28-G141_drz.fits','/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'][0:2], aper_list=[14], polyregions=glob.glob("COSMOS-*-F140W_asn.pointing.reg"))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('COSMOS-[0-9]*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='COSMOS-*[0-9]-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)
    
    # #### Test, combine three pointings with the same orientation, but use the single pointing as the detection image
    # threedhst.utils.combine_asn_shifts(['COSMOS-14-G141_asn.fits', 'COSMOS-2-G141_asn.fits', 'COSMOS-6-G141_asn.fits'], out_root='COSMOS-14o',
    #                    path_to_FLT='./', run_multidrizzle=True)
    #                    
    # os.system('cp COSMOS-14o* ../DATA')
    # os.chdir('../')
    # 
    # unicorn.go_3dhst.set_parameters(direct='F140W', LIMITING_MAGNITUDE=24)
    # threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/COSMOS-14-F140W_drz.fits'
    # threedhst.options['OTHER_BANDS'] = []
    # threedhst.process_grism.reduction_script(asn_grism_file='COSMOS-14o_asn.fits')
    # unicorn.analysis.make_SED_plots(grism_root='COSMOS-14o')
    # unicorn.go_3dhst.clean_up()
    
    #### Direct mosaic
    direct_files = glob.glob('COSMOS-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Direct mosaic
    direct_files = glob.glob('COSMOS-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.06
    #SCALE = 0.5
    PIXFRAC=0.6
    NX, NY = int(9670*0.06/SCALE), int(18890*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12634, dec=2.3336697,
             final_outnx = NX, final_outny=NY, ivar_weights=False)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12634, dec=2.3336697,
             final_outnx = NX, final_outny=NY, ivar_weights=False)
    
    ### Check the mosaic
    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits', 'COSMOS-G141_drz.fits', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04'], aper_list=[14], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))
    
    ### Temporary release of the direct mosaic
    iraf.imcopy('COSMOS-F140W_drz.fits[1]','COSMOS-F140w_11-06-17_sci.fits')
    iraf.imcopy('COSMOS-F140W_drz.fits[2]','COSMOS-F140w_11-06-17_wht.fits')
    #!tar czvf COSMOS-F140w_11-06-17.tar.gz COSMOS-*-F140W_shifts.txt COSMOS-*-F140W_tweak.fits COSMOS-*-F140W_asn.fits COSMOS-F140W_shifts.txt COSMOS-F140W_asn.fits
    #!mv COSMOS-F140w_11-06-17* ../MOSAIC
        
    threedhst.gmap.makeImageMap(['../MOSAIC/COSMOS-F140w_11-06-17_sci.fits[0]', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/Chandra/CC0570_img.fits[0]*1.5', '/3DHST/Ancillary/COSMOS/Spitzer/mips_24_GO3_sci_10.fits[0]*400', '/3DHST/Ancillary/COSMOS/VLA/vla_20cm_dp_sin_10.fits[0]*40000'], aper_list=[13,14,15,16,17], tileroot=['F140W','F814W','WIRDS-K','0.5-7keV', 'MIPS-24', 'VLA-20cm'])
    
    ### don't need high-res tiles of lo-res images sitting around
    os.system('rm ~/Sites/FITS/tiles/MIPS*1[67].png')
    os.system('rm ~/Sites/FITS/tiles/VLA*1[7].png')
    os.system('rm ~/Sites/FITS/tiles/0.5*17.png')

def GOODSS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/PREP_FLT')
    
    ALIGN_FILES = ['/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits', '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-2_F160W_v1_sci.fits', '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-1_F160W_v1_sci.fits']
    
    #### Main preparation loop
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        ALIGN = ALIGN_FILES[0]
        if pointing.startswith('GOODS-S-28'):
            ALIGN = ALIGN_FILES[1]
        #
        if pointing.startswith('GOODS-S-1'):
            ALIGN = ALIGN_FILES[2]
        #
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    # test
    threedhst.gmap.makeImageMap(['GOODS-S-23-F140W_drz.fits', 'GOODS-S-23-F140W_align.fits[0]*4', 'GOODS-S-23-G141_drz.fits'], aper_list=[13,14,15,16], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'))
    
    #### Direct mosaic
    direct_files = glob.glob('GOODS-S-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Direct mosaic
    direct_files = glob.glob('GOODS-S-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.128254
    #SCALE = 0.5
    PIXFRAC=0.6
    NX, NY = int(16210*0.06/SCALE), int(18100*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-S-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=53.154223, dec=-27.807325,
             final_outnx = NX, final_outny=NY)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-S-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=53.154223, dec=-27.807325,
             final_outnx = NX, final_outny=NY)
    
    threedhst.gmap.makeImageMap(['GOODS-S-F140W_drz.fits', 'GOODS-S-G141_drz.fits'], aper_list=[13,14,15], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'))
    
def AEGIS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'AEGIS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/AEGIS/ACS/mos_i_scale1*drz.fits'
    ALIGN = '/3DHST/Ancillary/AEGIS/WIRDS/WIRDS_Ks_141927+524056_T0002.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    
    threedhst.gmap.makeImageMap(['AEGIS-9-F140W_drz.fits', 'AEGIS-9-F140W_align.fits[0]*0.04', 'AEGIS-9-G141_drz.fits'][0:], aper_list=[14,15,16], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'))
    
    #### Direct mosaic
    direct_files = glob.glob('AEGIS-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='AEGIS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Direct mosaic
    direct_files = glob.glob('AEGIS-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='AEGIS-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    ## some changes
    
    SCALE = 0.128254
    #SCALE = 0.5
    PIXFRAC=0.6
    NX, NY = int(18000*0.06/SCALE), int(19000*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('AEGIS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=214.88705, dec=52.85442,
             final_outnx = NX, final_outny=NY)
    
    threedhst.prep_flt_files.startMultidrizzle('AEGIS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=214.88705, dec=52.85442,
             final_outnx = NX, final_outny=NY)
    #
    threedhst.gmap.makeImageMap(['AEGIS-F140W_drz.fits', 'AEGIS-G141_drz.fits'], aper_list=[13,14,15], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'))

def SN_MARSHALL():
    ####********************************************####
    ####              SN-MARSHALL (UDS)
    ####********************************************####
    
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair

    os.chdir(unicorn.GRISM_HOME+'SN-MARSHALL/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_drz.fits'
    
    direct = glob.glob('../RAW/ib*80_asn.fits')
    grism = glob.glob('../RAW/ib*70_asn.fits')
    
    for i in range(1,6):        
        dir_i=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        gris_i=threedhst.prep_flt_files.make_targname_asn(grism[i], newfile=True)
        pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    
    ## Other direct images
    direct = glob.glob('../RAW/*[43]0_asn.fits')
    olds = glob.glob('MARSHALL-[BC]-F*shifts.txt')
    for old in olds:
        os.remove(old)
    
    for i in range(1,len(direct)):
        dir_i=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if os.path.exists(dir_i.replace('_asn.fits','_shifts.txt')):
            shutil.move(dir_i.replace('_asn.fits','_shifts.txt'),'old_shifts.txt')
            shutil.move(dir_i, 'old_asn.fits')
            combine=True
        else:
            combine=False
        #    
        dir_i=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        pair(dir_i, None, ALIGN_IMAGE = ALIGN, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False)
        #
        if combine:
            threedhst.utils.combine_asn_shifts([dir_i, 'old_asn.fits'],
                        out_root=dir_i.split('_asn')[0], path_to_FLT='./',
                        run_multidrizzle=False)
        
    #### Direct images
    asn_list = glob.glob('MARSHALL-?-F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-F160W',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    asn_list = glob.glob('MARSHALL-?-F125W_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-F125W',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    #### First orientation
    asn_list = glob.glob('MARSHALL-[135]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-225-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-225-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    #### Second orientation
    asn_list = glob.glob('MARSHALL-[246]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-245-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-245-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    #### Check
    # threedhst.gmap.makeImageMap(['MARSHALL-1-F160W_drz.fits', 'MARSHALL-1-F160W_align.fits[0]', 'MARSHALL-1-G141_drz.fits'], aper_list=[13,14,15])
    threedhst.gmap.makeImageMap(['MARSHALL-F125W_drz.fits', 'MARSHALL-F160W_drz.fits', 'MARSHALL-1-F160W_align.fits[0]', 'MARSHALL-225-G141_drz.fits', 'MARSHALL-245-G141_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[13,14,15,16])
    
#
def DADDI():
    ####******************************************************####
    ####              HIGHZ-CLUSTER (z>2?) from Daddi
    ####******************************************************####
    
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair

    os.chdir(unicorn.GRISM_HOME+'DADDI/PREP_FLT')
    ALIGN = None
    
    direct = glob.glob('../RAW/ib*10_asn.fits')
    grism = glob.glob('../RAW/ib*20_asn.fits')
    
    for i in range(len(direct)):        
        dir_i=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        gris_i=threedhst.prep_flt_files.make_targname_asn(grism[i], newfile=True)
        pair(direct[i], grism[i], ALIGN_IMAGE = None, TWEAKSHIFTS_ONLY=True, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    #
    threedhst.gmap.makeImageMap(['HIGHZ-CLUSTER-4-F140W_drz.fits', 'HIGHZ-CLUSTER-4-G141_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[14, 15, 16])
    
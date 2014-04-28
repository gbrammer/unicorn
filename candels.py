import os
import glob
import shutil

import numpy as np
import pyfits
import scipy.linalg
import matplotlib.pyplot as plt

import threedhst
import threedhst.prep_flt_files
import threedhst.catIO as catIO
import threedhst.eazyPy as eazy

import unicorn

def goods_ers():
    """
    GOODS-ERS field (not candels)
    """
    import unicorn.candels
    import threedhst
    import threedhst.prep_flt_files
    
    unicorn.candels.make_asn_files()
    
    #ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'
        
    files=glob.glob('WFC3-*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    #
    #### Make detection images for the grism
    os.system('cp '+unicorn.GRISM_HOME+'/ERS/PREP_FLT/WFC3-ERSII-G01-F140W_drz.fits .')
    
    #### F098M
    os.system('cp WFC3-ERSII-G01-F140W_drz.fits WFC3-ERSII-G01-F098M_drz.fits')
    threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='WFC3-ERSII-IR*-F098M',
                                    pointing='WFC3-ERSII-G01-F098M',
                                    run_multidrizzle=True, grow=200)
     
    #### F125W
    os.system('cp WFC3-ERSII-G01-F140W_drz.fits WFC3-ERSII-G01-F125W_drz.fits')
    threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='WFC3-ERSII-IR*-F125W',
                                    pointing='WFC3-ERSII-G01-F125W',
                                    run_multidrizzle=True, grow=200)
    
    #### F160W
    os.system('cp WFC3-ERSII-G01-F140W_drz.fits WFC3-ERSII-G01-F160W_drz.fits')
    threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='WFC3-ERSII-IR*-F160W',
                                    pointing='WFC3-ERSII-G01-F160W',
                                    run_multidrizzle=True, grow=200)
    
    ### Test
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN_IMAGE), matchImage='WFC3-ERSII-G01-F160W_drz.fits', match_extension=1, output='ERS-F850LP.fits')
    threedhst.gmap.makeImageMap(['ERS-F850LP.fits[0]*4','WFC3-ERSII-G01-F098M_drz.fits[1]*2', 'WFC3-ERSII-G01-F140W_drz.fits', 'WFC3-ERSII-G01-F125W_drz.fits', 'WFC3-ERSII-G01-F160W_drz.fits',], aper_list=[14,15], tileroot=['z850','f098m','f140w','f125w','f160w'])
    
def egs():
    import unicorn.candels
    
    os.chdir('/Users/gbrammer/CANDELS/EGS/PREP_FLT/')
    unicorn.candels.make_asn_files()
    
    #ALIGN_IMAGE = 'AEGIS-N2_K_sci.fits'
    ALIGN_IMAGE = '/3DHST/Ancillary/AEGIS/WIRDS/WIRDS_Ks_141927+524056_T0002.fits'
    ALIGN_IMAGE = '/3DHST/Ancillary/AEGIS/ACS/mos_i_scale2_drz.fits'
    ALIGN_IMAGE = '/3DHST/Ancillary/AEGIS/CANDELS/hlsp_candels_hst_wfc3_egsa01_f160w_v0.5_drz.fits'

    files=glob.glob('EGS-V*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=False, DIRECT_HIGHER_ORDER=2,
                SCALE=0.06, geometry='rxyscale,shift')
    
    # V0T-F125 71/1I-F160Wsatellite trail, V54-125 - bad shift
    root = 'EGS-V72-F125W'

    threedhst.shifts.refine_shifts(ROOT_DIRECT=root,
      ALIGN_IMAGE=ALIGN_IMAGE,
      fitgeometry='shift', clean=True,
      ALIGN_EXTENSION=0)
    #
    threedhst.prep_flt_files.startMultidrizzle(root+'_asn.fits',
         use_shiftfile=True, skysub=False,
         final_scale=0.06, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False)
    #
    
    for filt in ['F125W', 'F160W']:
        files=glob.glob('EGS-V*-'+filt+'_asn.fits')
        threedhst.utils.combine_asn_shifts(files, out_root='EGS-'+filt, path_to_FLT='./', run_multidrizzle=False)
        #
        SCALE = 0.06
        threedhst.prep_flt_files.startMultidrizzle( 'EGS-'+filt+'_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=214.92061, dec=52.878457, final_outnx=7547, 
             final_outny=31900, final_rot=42)
    #
    threedhst.shifts.plot_shifts('EGS-F125W', ALIGN_IMAGE, skip_swarp=False)
    
    threedhst.gmap.makeImageMap(['/3DHST/Ancillary/AEGIS/ACS/mos_i_scale2_drz.fits[0]*2','EGS-epoch2-F125W_drz.fits', 'EGS-epoch2-F160W_drz.fits'], aper_list=[14], tileroot=['F814W','F125W','F160W'], polyregions=glob.glob('*F125W_asn.pointing.reg'))
    
    threedhst.gmap.makeImageMap(['/3DHST/Ancillary/AEGIS/WIRDS/WIRDS_Ks_141927+524056_T0002.fits[0]*0.04', '/3DHST/Ancillary/AEGIS/ACS/mos_i_scale2_drz.fits[0]*5', 'PREP_FLT/EGS-epoch2-F125W_drz.fits', 'PREP_FLT/EGS-epoch2-F160W_drz.fits'], aper_list=[12], tileroot=['WIRDS','F814W', 'F125W', 'F160W'])

def tile41_prep():
    
    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()    
    
    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/COSMOS/CANDELS/ASTRODRIZZLE/cosmos-f160w-astrodrizzle-v0.1_drz_sci.fits'
    
    files = glob.glob('TILE*asn.fits')

    for file in files:
        #if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')    
    
def cosmos_prep():
    
    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    #ALIGN_IMAGE = '/3DHST/Ancillary/COSMOS/CANDELS/candels_public/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits'
    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/COSMOS/CANDELS/ASTRODRIZZLE/cosmos-f160w-astrodrizzle-v0.1_drz_sci.fits'
    
    for filter in ['F125W','F160W']:
        files=glob.glob('COSMOS-V*'+filter+'_asn.fits')
        threedhst.utils.combine_asn_shifts(files, out_root='COSMOS-'+filter,
                           path_to_FLT='./', run_multidrizzle=False)
        for file in files:
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                    SCALE=0.06, geometry='rotate,shift')    

def cosmos():
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 

    INDEF = iraf.INDEF
    no = iraf.no
    yes = iraf.yes
    
    import unicorn.candels
    
    os.chdir('/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/')
    unicorn.candels.make_asn_files()
    
    #ALIGN_IMAGE = 'AEGIS-N2_K_sci.fits'
    ALIGN_IMAGE = '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits'

    files=glob.glob('COSMOS-V*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                SCALE=0.06, geometry='rotate,shift')
    #
    threedhst.gmap.makeImageMap(['COSMOS-V12-F125W_drz.fits','COSMOS-V12-F160W_drz.fits', 'COSMOS-V12-F160W_align.fits[0]*4',], aper_list=[15,16])
    
    # V0T-F125 71/1I-F160Wsatellite trail, V54-125 - bad shift
    root = 'EGS-V72-F125W'

    threedhst.shifts.refine_shifts(ROOT_DIRECT=root,
      ALIGN_IMAGE=ALIGN_IMAGE,
      fitgeometry='shift', clean=True,
      ALIGN_EXTENSION=0)
    #
    threedhst.prep_flt_files.startMultidrizzle(root+'_asn.fits',
         use_shiftfile=True, skysub=False,
         final_scale=0.06, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False)
    #
    
    ##### DQ flag: F125W - 19,23,27
    #####          F160W - 7,11,13,15,23,24,25,49
    ## done: 19,23,24,27 ,7
    for visit in [11,13,15,25,49]:
        threedhst.dq.checkDQ(asn_direct_file='COSMOS-V%02d-F125W_asn.fits' %(visit), asn_grism_file='COSMOS-V%02d-F160W_asn.fits' %(visit), path_to_flt='../RAW/', SIMPLE_DS9=False)
    
    bad_visits = [19,23,27] # epoch1
    bad_visits = [67,51] #epoch2
    for visit in bad_visits:
        threedhst.dq.checkDQ(asn_direct_file='COSMOS-V%02d-F125W_asn.fits' %(visit), asn_grism_file='COSMOS-V%02d-F125W_asn.fits' %(visit), path_to_flt='../RAW/', SIMPLE_DS9=False)
        os.remove('COSMOS-V%02d-F125W_drz.fits' %(visit))
    #
    bad_visits = [7,11,13,15,23,24,25,49]  #epoch1
    bad_visits = [94,93,89,79,78,52]  #epoch2
    for visit in bad_visits:
        threedhst.dq.checkDQ(asn_direct_file='COSMOS-V%02d-F160W_asn.fits' %(visit), asn_grism_file='COSMOS-V%02d-F160W_asn.fits' %(visit), path_to_flt='../RAW/', SIMPLE_DS9=False)
        os.remove('COSMOS-V%02d-F160W_drz.fits' %(visit))
    
    # Redo BG subtraction for visits with trails masked
    files=glob.glob('COSMOS-V*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=False, DIRECT_HIGHER_ORDER=2,
                SCALE=0.06, geometry='rotate,shift')
    
    #### Align everything to the F160W image from epoch1
    #TEST
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN_IMAGE), matchImage='COSMOS-V45-F160W_drz.fits', match_extension=1, output='tmp.fits')
    
    files=glob.glob('COSMOS-V[0-4]*-F160W_asn.fits')
    files.append('COSMOS-V50-F160W_asn.fits')
    for file in files:
        visit = int(file.split('-')[1][1:])
        print visit
        ### Redo because messed up with a typo earlier
        threedhst.shifts.refine_shifts(ROOT_DIRECT='COSMOS-V%02d-F160W' %(visit), ALIGN_IMAGE='COSMOS-V%02d-F125W_drz.fits' %(visit), fitgeometry='shift', clean=True, ALIGN_EXTENSION=1)
        threedhst.prep_flt_files.startMultidrizzle('COSMOS-V%02d-F160W_asn.fits' %(visit), use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
        # threedhst.shifts.refine_shifts(ROOT_DIRECT='COSMOS-V%02d-F125W' %(visit), ALIGN_IMAGE=file.replace('asn','drz'), fitgeometry='shift', clean=True, ALIGN_EXTENSION=1)
        # threedhst.prep_flt_files.startMultidrizzle('COSMOS-V%02d-F125W_asn.fits' %(visit), use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
        # #
        # threedhst.shifts.refine_shifts(ROOT_DIRECT='COSMOS-V%02d-F125W' %(visit+44), ALIGN_IMAGE=file.replace('asn','drz'), fitgeometry='shift', clean=True, ALIGN_EXTENSION=1)
        # threedhst.prep_flt_files.startMultidrizzle('COSMOS-V%02d-F125W_asn.fits' %(visit+44), use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
        #
        # threedhst.shifts.refine_shifts(ROOT_DIRECT='COSMOS-V%02d-F160W' %(visit+44), ALIGN_IMAGE=file.replace('asn','drz'), fitgeometry='shift', clean=True, ALIGN_EXTENSION=1)
        # threedhst.prep_flt_files.startMultidrizzle('COSMOS-V%02d-F160W_asn.fits' %(visit+44), use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
    
    for filt in ['F125W', 'F160W']:
        # files=glob.glob('COSMOS-V[0-4]*-'+filt+'_asn.fits')
        # files.append('COSMOS-V50-'+filt+'_asn.fits')
        # out='COSMOS-epoch1-'+filt
        #
        # files=glob.glob('COSMOS-V[5-9]*-'+filt+'_asn.fits')
        # xx = files.pop(0)
        # out='COSMOS-epoch2-'+filt
        #
        files=glob.glob('COSMOS-V*-'+filt+'_asn.fits')
        out='COSMOS-full-'+filt
        #
        threedhst.utils.combine_asn_shifts(files, out_root=out, path_to_FLT='./', run_multidrizzle=False)
        SCALE = 0.06
        threedhst.prep_flt_files.startMultidrizzle(out+'_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.1291915, dec=2.36638225, final_outnx=4667*0.1283/SCALE, 
             final_outny=11002*0.1283/SCALE, build_drz=False)
    
    ###  Combined detection image
    files=glob.glob('COSMOS-V*-F???W_asn.fits')
    out='COSMOS-J+H'
    threedhst.utils.combine_asn_shifts(files, out_root=out, path_to_FLT='./', run_multidrizzle=False)
    SCALE = 0.06
    threedhst.prep_flt_files.startMultidrizzle(out+'_asn.fits',
         use_shiftfile=True, skysub=False,
         final_scale=SCALE, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False,
         ra=150.1291915, dec=2.36638225, final_outnx=4667*0.1283/SCALE, 
         final_outny=11002*0.1283/SCALE, build_drz=False)
    #
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN_IMAGE), matchImage='COSMOS-full-F160W_drz_sci.fits', match_extension=0, output='COSMOS-full-ACSi.fits')
    
    #
    iraf.imcalc('COSMOS-epoch2-F125W_drz_sci.fits,COSMOS-epoch1-F125W_drz_sci.fits', "diff_F125W.fits", "im1-im2")
    iraf.imcalc('COSMOS-epoch2-F160W_drz_sci.fits,COSMOS-epoch1-F160W_drz_sci.fits', "diff_F160W.fits", "im1-im2")
    
    threedhst.shifts.plot_shifts('COSMOS-F125W', ALIGN_IMAGE, skip_swarp=False)
    
    
    # p!sh ~/Sites/FITS/keys.sh

def uds_prep():
    
    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    ALIGN_IMAGE = '/3DHST/Ancillary/UDS/CANDELS/candels_public/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'

    for filter in ['F125W','F160W']:
        files = glob.glob('UDS-V*'+filter+'_asn.fits')
        files.remove('UDS-V5K-'+filter+'_asn.fits')
        threedhst.utils.combine_asn_shifts(files, out_root='UDS-'+filter,
                           path_to_FLT='./', run_multidrizzle=False)
        for file in files:
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.06, geometry='rotate,shift')       
 
    filter = 'F140W'
    files = glob.glob('UDS-*-'+filter+'_asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                SCALE=0.06, geometry='rotate,shift')            
    threedhst.utils.combine_asn_shifts(files, out_root='UDS-'+filter,
                           path_to_FLT='./', run_multidrizzle=False)
    
def uds_marshall_prep():
    
    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    ALIGN_IMAGE = '/3DHST/Ancillary/UDS/CANDELS/candels_public/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'
    
    for filter in ['F125W','F160W']:
        files=glob.glob('MARSHALL-*'+filter+'_asn.fits')
        for file in files:
            #if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                    SCALE=0.06, geometry='rotate,shift')       

def uds():
    import unicorn.candels
    
    #unicorn.candels.make_asn_files()
    
    os.chdir('/Users/gbrammer/CANDELS/UDS/PREP_FLT')
    
    list = catIO.Readfile('files.info')
    epoch1 = list.date_obs < '2010-12-15'
    epoch2 = list.date_obs > '2010-12-15'
    
    ALIGN_IMAGE='../UKIDSS/UDS_K.fits'
    ALIGN_IMAGE = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f125w_v0.5_drz.fits'
        
    ########################
    ### Epoch1
    ########################
    visits = np.unique(list.targname[epoch1])
    for visit in visits:
        filters = np.unique(list.filter[list.targname == visit])
        for filter in filters:
            file = visit+'-'+filter+'_asn.fits'
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.06, geometry='rxyscale,shift')
    #
    root='UDS-V4Z-F160W'
    threedhst.shifts.refine_shifts(ROOT_DIRECT=root,
      ALIGN_IMAGE=ALIGN_IMAGE,
      fitgeometry='shift', clean=False,
      ALIGN_EXTENSION=0)
    
    threedhst.prep_flt_files.startMultidrizzle(root+'_asn.fits',
         use_shiftfile=True, skysub=False,
         final_scale=0.06, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False)
    
    ## Mosaic
    f125w_list = []
    f160w_list = []
    for visit in visits:
        f125w_list.append(visit+'-F125W_asn.fits')
        f160w_list.append(visit+'-F160W_asn.fits')
    
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='UDS-epoch1-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='UDS-epoch1-F160W', path_to_FLT='./', run_multidrizzle=False)
    
    for filt in ['F125W', 'F160W']:
        SCALE = 0.06
        ra, dec, nx, ny = 34.4081405556, -5.20018201235, 10422, 4349
        threedhst.prep_flt_files.startMultidrizzle( 'UDS-epoch1-'+filt+'_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec, final_outnx=nx*0.128254/SCALE, final_outny=ny*0.128254/SCALE)
             
    ########################
    ### Epoch2
    ########################
    visits = np.unique(list.targname[epoch2])
    for visit in visits:
        filters = np.unique(list.filter[list.targname == visit])
        for filter in filters:
            file = visit+'-'+filter+'_asn.fits'
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.06, geometry='rxyscale,shift')
    ## Mosaic
    f125w_list = []
    f160w_list = []
    for visit in visits:
        f125w_list.append(visit+'-F125W_asn.fits')
        f160w_list.append(visit+'-F160W_asn.fits')
    
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='UDS-epoch2-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='UDS-epoch2-F160W', path_to_FLT='./', run_multidrizzle=False)
    
    for filt in ['F125W', 'F160W']:
        SCALE = 0.06
        ra, dec, nx, ny = 34.4081405556, -5.20018201235, 10422, 4349
        threedhst.prep_flt_files.startMultidrizzle( 'UDS-epoch2-'+filt+'_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec, final_outnx=nx*0.128254/SCALE, final_outny=ny*0.128254/SCALE)

    ########################
    ### Full Mosaic
    ########################
    f125w_list = glob.glob('UDS-V*-F125W_asn.fits')
    f160w_list = glob.glob('UDS-V*-F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='UDS-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='UDS-F160W', path_to_FLT='./', run_multidrizzle=False)
    
    for filt in ['F125W', 'F160W']:
        SCALE = 0.06
        ra, dec, nx, ny = 34.4081405556, -5.20018201235, 10422, 4349
        threedhst.prep_flt_files.startMultidrizzle( 'UDS-'+filt+'_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec, final_outnx=nx*0.128254/SCALE, final_outny=ny*0.128254/SCALE)
    
    #### Check pointings
    threedhst.shifts.plot_shifts('UDS-F125W', ALIGN_IMAGE, skip_swarp=False)
    
    
    #### Matched ACS images
    os.chdir('/3DHST/Ancillary/UDS/CANDELS')
    for band in ['f814w', 'f606w']:
        threedhst.shifts.matchImagePixels(input= glob.glob('hlsp_candels_hst_acs_uds-tot_%s_v1.0_drz.fits' %(band)), matchImage='hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits', match_extension=0, output='UDS-%s.fits' %(band))
    #
    os.chdir('/3DHST/Ancillary/UDS/CANDELS')

    threedhst.shifts.matchImagePixels(input= glob.glob('../UKIDSS/UDS_K.fits'), matchImage='hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits', match_extension=0, output='UDS-K.fits')

def cdfs_prep():

    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'

    if not only_f140w and not only_flat:
        for filter in ['F125W','F160W']:
            files=glob.glob('GOODS-*-*-'+filter+'_asn.fits')
            files.remove('GOODS-SD2-V7G-'+filter+'_asn.fits')
            files.remove('GOODS-SD5-VGQ-'+filter+'_asn.fits')
            files.remove('GOODS-SD5-VGX-'+filter+'_asn.fits')
            files.append('GOODS-TESTORB-'+fiter+'_asn.fits')
            threedhst.utils.combine_asn_shifts(files, out_root='GOODS-S-'+filter,
                path_to_FLT='./', run_multidrizzle=False)
            for file in files:
                if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file, 
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                        SCALE=0.06, geometry='rotate,shift')  

    if only_f140w:
        for filter in ['F140W']:
            files = glob.glob('GOODS-SOUTH-*'+filter+'_asn.fits')
            threedhst.utils.combine_asn_shifts(files, out_root='GOODS-S-'+filter,
                path_to_FLT='./', run_multidrizzle=False)
            for file in files:
                if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file, 
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                        SCALE=0.06, geometry='rotate,shift')    

    if only_flat:
        for filter in ['F125W','F160W']:
            files=glob.glob('GOODS-*-*'+filter+'_asn.fits')
            files.remove('GOODS-SD2-V7G-'+filter+'_asn.fits')
            files.remove('GOODS-SD5-VGQ-'+filter+'_asn.fits')
            files.remove('GOODS-SD5-VGX-'+filter+'_asn.fits')
            for file in files:
                unicorn.candels.prep_candels(asn_file=file, 
                    GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
                    redo_segmentation=False)    
                    
                    
    ### Making GOODS-S-ALL asn tables
    ### On hyperion in /Volumes/3D-HST/MOSAICS/GOODS-S/ASN
    files = glob.glob('*F125W*asn.fits')
    ### Contains the following files:
    ### ['ERS-F125W_asn.fits', 'GEORGE1-F125W_asn.fits', 'GEORGE2-F125W_asn.fits', 'GOODS-S-F125W_asn.fits', 'HUDF09-01-DEEP-F125W_asn.fits', 'HUDF09-02-DEEP-F125W_asn.fits', 'HUDF09-DEEP-F125W_asn.fits']
    threedhst.utils.combine_asn_shifts(files, out_root='GOODS-S-ALL-F125W',path_to_FLT='../ALL_FLT/', run_multidrizzle=False)
    
    files = glob.glob('*F160W*asn.fits')
    ### ['ERS-F160W_asn.fits', 'GEORGE1-F160W_asn.fits', 'GEORGE2-F160W_asn.fits', 'GOODS-S-F160W_asn.fits', 'HUDF09-01-DEEP-F160W_asn.fits', 'HUDF09-02-DEEP-F160W_asn.fits', 'HUDF09-DEEP-F160W_asn.fits', 'HUDF12-DEEP-F160W_asn.fits']
    threedhst.utils.combine_asn_shifts(files, out_root='GOODS-S-ALL-F160W',path_to_FLT='../ALL_FLT/', run_multidrizzle=False)
    
    files=glob.glob('*F140W*asn.fits')
    ### ['CDFS-AGN1-F140W_asn.fits', 'GOODS-S-F140W_asn.fits', 'HUDF12-DEEP-F140W_asn.fits', 'WFC3-ERSII-G01-F140W_asn.fits']
    threedhst.utils.combine_asn_shifts(files, out_root='GOODS-S-ALL-F140W',path_to_FLT='../ALL_FLT/', run_multidrizzle=False)

def george_prep():

    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'

    files = glob.glob('GEORGE*asn.fits')
    
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')  
    
def cdfs_agn1_prep():

    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'
     
    files = glob.glob('CDFS*asn.fits')
    
    for file in files:
        #if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')  
        
def hudf12_prep():

    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()
    
    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'
    
    for filter,order in zip(['F140W','F160W'],[2,1]):
        files=glob.glob('HUDF-DEEP-WFC3-*-*-'+filter+'_asn.fits')
        for file in files:
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=order,
                    SCALE=0.06, geometry='rotate,shift')  

    ### make ASN tables for each band:
    
    for filter in ['F140W','F160W']:
        files = glob.glob('HUDF-DEEP-WFC3-*-%s_asn.fits'%(filter))
        threedhst.utils.combine_asn_shifts(files, out_root='HUDF12-DEEP-%s'%(filter), 
            path_to_FLT='./', run_multidrizzle=False)
        threedhst.prep_flt_files.startMultidrizzle('HUDF12-DEEP-%s_asn.fits'%(filter),
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)
        
    

def hudf_prep():
    
    import unicorn.candels
    import glob
    from threedhst import catIO
    
    ### Reduction for the pointing which matches GOODS-S-01 F140W

    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-1_F160W_v1_sci.fits'  

    files = glob.glob('HUDF05-01-DEEP-WFC3*asn.fits')

    for file in files:
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')  

    for filter in ['F125W','F160W']:
        files = glob.glob('HUDF05-01-DEEP-WFC3-*-%s_asn.fits'%(filter))
        threedhst.utils.combine_asn_shifts(files, out_root='HUDF09-01-DEEP-%s'%(filter), 
            path_to_FLT='./', run_multidrizzle=False)
        threedhst.prep_flt_files.startMultidrizzle('HUDF09-01-DEEP-%s_asn.fits'%(filter),
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)


    ### Reduction for the pointing which matches GOODS-S-28 in F140W
    ### removed ib5x5b*, ib5x5c* and ib5x5f form the files.info file --> lost pointing
    ### Having problems processing it so break it down by date

    info = catIO.Readfile('files.info')
    asn = threedhst.utils.ASNFile('../RAW/ib5x59020_asn.fits')
    
    ### Make ASN files for each visit (by date_obs)
    hudf02 = info.targname == 'HUDF05-02-DEEP-WFC3'
    for filter in ['F125W', 'F160W']:
        sub = hudf02 & (info.filter == filter)
        dates = np.unique(info.date_obs[sub])
        for date in dates:
            visit = sub & (info.date_obs == date)
            files = info.file[visit]
            asn.exposures = [file.split('_flt')[0] for file in files]
            name = 'HUDF-02-%s-%s' %(''.join(date.split('-')[1:]), filter)
            asn.product = name
            asn.write(name+'_asn.fits')
            print name, asn.exposures
    #
    
    ##### Process visits separately
    files=glob.glob('HUDF-02-0*F1*W_asn.fits')
    
    FORCE=False
    
    #### Just run tweakshifts
    for file in files:
        threedhst.process_grism.fresh_flt_files(file)
        threedhst.prep_flt_files.threedhst.shifts.run_tweakshifts(file)
    
    #### For some reason, the second exposure in some visits has a large
    #### offset.  MANUALLY change it to zero or something close and then run the prep
    #### scripts.
    for file in files:
        if ((not os.path.exists(file.replace('asn','drz'))) & (not os.path.exists(file.replace('asn','drz')+'.gz'))) | FORCE:
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = None, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
        
    #### These visits are just aligned internally to check the 
    #### tweakshifts.  You could then do "refine_shifts" to align them
    #### to the reference, something like:
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-2_F160W_v1_sci.fits'  
    
    for file in files:
        for geom in ['rotate','shift']:
            threedhst.shifts.refine_shifts(ROOT_DIRECT=file.split('_asn')[0], 
              ALIGN_IMAGE=ALIGN_IMAGE, ALIGN_EXTENSION=0,  
              fitgeometry=geom, clean=True)
            #
            threedhst.prep_flt_files.startMultidrizzle(file,
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
                 
    for filter in ['F125W','F160W']:
        files = glob.glob('HUDF05-02-DEEP-WFC3-*-%s_asn.fits'%(filter))
        threedhst.utils.combine_asn_shifts(files, out_root='HUDF09-02-DEEP-%s'%(filter), 
            path_to_FLT='./', run_multidrizzle=False)
        threedhst.prep_flt_files.startMultidrizzle('HUDF09-02-DEEP-%s_asn.fits'%(filter),
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)


    ### Reduce the cental HUDF pointing
    
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09_F160W_v1_sci.fits'
    
    ## removed ib5x0ei6q_flt.fits, ib5x0ei8q_flt.fits and ib5x0eicq_flt.fits 
    ## from files.info --> weird persistence pattern
    unicorn.candels.make_asn_files()
    files = glob.glob('HUDF-DEEP-WFC3-*_asn.fits')

    for file in files:
        unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')  

    for filter in ['F125W','F160W']:
        files = glob.glob('HUDF-DEEP-WFC3-*-%s_asn.fits'%(filter))
        threedhst.utils.combine_asn_shifts(files, out_root='HUDF09-DEEP-%s'%(filter), 
            path_to_FLT='./', run_multidrizzle=False)
        threedhst.prep_flt_files.startMultidrizzle('HUDF09-DEEP-%s_asn.fits'%(filter),
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)

def hudf_v3_prep():
    
    ### HUDF
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09_F160W_v1_sci.fits'
    
    files_f125w = glob.glob('ib5x0[cde]020_asn.fits')
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f125w[0])
    for file in files_f125w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f125w[0])
    
    files_f160w = glob.glob('ib5x1[efg]020_asn.fits')
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f160w[0])
    for file in files_f160w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f160w[0])
    
    
    for file in [files_f125w[0],files_f160w[0]]:
        unicorn.candels.prep_candels(asn_file=file,
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
            GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
            
    ### Flanking pointing 1
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-1_F160W_v1_sci.fits'  
    files_f125w, files_f160w = [],[]
    
    files_f125w = glob.glob('ib5x3[9ad]020_asn.fits')
    
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f125w[0])
    for file in files_f125w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        print 'Appending '+files_f125w[0]+' with '+file+'.'
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f125w[0])
    
    #files_f160w = glob.glob('ib5x4[578]020_asn.fits')
    files_f160w = glob.glob('ib5x4[578cf]020_asn.fits')
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f160w[0])
    for file in files_f160w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        print 'Appending '+files_f160w[0]+' with '+file+'.'
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f160w[0])
    
    
    for file in [files_f125w[0],files_f160w[0]]:
        unicorn.candels.prep_candels(asn_file=file,
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
            GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
            
    ### Flanking pointing 2
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/HUDF09/hlsp_hudf09_hst_wfc3ir_hudf09-2_F160W_v1_sci.fits'  
    files_f125w, files_f160w = [],[]
    
    files_f125w = glob.glob('ib5x5[ade]020_asn.fits')
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f125w[0])
    for file in files_f125w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        print 'Appending '+files_f125w[0]+' with '+file+'.'
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f125w[0])
    
    files_f160w = glob.glob('ib5x6[5678]020_asn.fits')
    asn1 = threedhst.utils.ASNFile('../RAW/'+files_f160w[0])
    for file in files_f160w[1:]:
        asn_tmp = threedhst.utils.ASNFile('../RAW/'+file)
        print 'Appending '+files_f160w[0]+' with '+file+'.'
        asn1.append(asn_tmp)
    asn1.write(out_file = files_f160w[0])
    
    
    for file in [files_f125w[0],files_f160w[0]]:
        unicorn.candels.prep_candels(asn_file=file,
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
            GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
            
    files = ['ib5x39020_asn.fits','ib5x5a020_asn.fits']
    root = 'HUDF-DEEP'
    filter = 'F125W'
    threedhst.utils.combine_asn_shifts(files, out_root='%s-%s' %(root, filter), path_to_FLT='./', run_multidrizzle=False)

    files = ['ib5x45020_asn.fits','ib5x65020_asn.fits']
    root = 'HUDF-DEEP'
    filter = 'F160W'
    threedhst.utils.combine_asn_shifts(files, out_root='%s-%s' %(root, filter), path_to_FLT='./', run_multidrizzle=False)

def cdfs():
    import unicorn.candels

    os.chdir('/Users/gbrammer/CANDELS/GOODS-S/PREP_FLT')
    
    unicorn.candels.make_asn_files()
    
    list = catIO.Readfile('files.info')
    
    list.group = []
    for targ in list.targname:
        list.group.append('-'.join(targ.split('-')[0:2]))
    
    list.group = np.array(list.group)
    
    epoch1 = list.group == 'GOODS-SD5'
    epoch2 = list.group == 'GOODS-SD2'
    epoch3 = list.group == 'GOODS-S075'
    epoch4 = list.group == 'GOODS-SDW'

    wide1 = list.group == 'GOODS-W1'
    wide2 = list.group == 'GOODS-W2'

    ALIGN_IMAGE='/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz*fits'
    #ALIGN_IMAGE='../PREP_FLT/goods-s_f160w_drz.fits'
    
    ########################
    ### Deep epochs
    ########################
    epochs = [epoch1, epoch2, epoch3, epoch4]
    
    visits = np.unique(list.targname)
    for visit in visits:
        filters = np.unique(list.filter[list.targname == visit])
        for filter in filters:
            file = visit+'-'+filter+'_asn.fits'
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.06, geometry='rxyscale,shift')

        ## Mosaic
        f125w_list = []
        f160w_list = []
        for visit in visits:
            f125w_list.append(visit+'-F125W_asn.fits')
            f160w_list.append(visit+'-F160W_asn.fits')
        #
        threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-epoch%d-F125W' %(e), path_to_FLT='./', run_multidrizzle=False)
        threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-epoch%d-F160W' %(e), path_to_FLT='./', run_multidrizzle=False)
        #
        SCALE = 0.128254
        for filt in ['F125W', 'F160W']:
            ra, dec, nx, ny = 53.127353, -27.798616, int(5670*0.128254/SCALE), int(5220*0.128254/SCALE)
            threedhst.prep_flt_files.startMultidrizzle( 'GOODS-epoch%d-' %(e)+filt+'_asn.fits',
                use_shiftfile=True, skysub=False,
                final_scale=SCALE, pixfrac=1.0, driz_cr=False,
                updatewcs=False, clean=True, median=False,
                ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    #
    files=glob.glob('GOODS-S*drz.fits')
    files.extend(glob.glob('*align.fits'))
    for file in files: 
        os.remove(file)
        
    #### Full deep mosaic
    f105w_list = glob.glob('GOODS-S*F105W_asn.fits')
    f125w_list = glob.glob('GOODS-S*F125W_asn.fits')
    f160w_list = glob.glob('GOODS-S*F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(f105w_list,out_root='GOODS-deep-F105W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-deep-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-deep-F160W', path_to_FLT='./', run_multidrizzle=False)

    SCALE = 0.06
    ra, dec, nx, ny = 53.127353, -27.798616, int(5670*0.128254/SCALE), int(5220*0.128254/SCALE)
    for filt in ['F105W','F125W', 'F160W']:
        threedhst.prep_flt_files.startMultidrizzle('GOODS-deep-'+filt+ '_asn.fits',
        use_shiftfile=True, skysub=False,
        final_scale=SCALE, pixfrac=0.8, driz_cr=False,
        updatewcs=False, clean=True, median=False,
        ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    
    #### Full wide mosaic
    f105w_list = glob.glob('GOODS-W*F105W_asn.fits')
    f125w_list = glob.glob('GOODS-W*F125W_asn.fits')
    f160w_list = glob.glob('GOODS-W*F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(f105w_list,out_root='GOODS-wide-F105W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-wide-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-wide-F160W', path_to_FLT='./', run_multidrizzle=False)

    SCALE = 0.06
    ra, dec, nx, ny = 53.166591, -27.898528, int(11320*0.06/SCALE), int(7900*0.06/SCALE)
    for filt in ['F105W','F125W', 'F160W']:
        threedhst.prep_flt_files.startMultidrizzle('GOODS-wide-'+filt+ '_asn.fits',
        use_shiftfile=True, skysub=False,
        final_scale=SCALE, pixfrac=0.8, driz_cr=False,
        updatewcs=False, clean=True, median=False, ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    
    #### Check pointings
    threedhst.shifts.plot_shifts('GOODS-deep-F125W', ALIGN_IMAGE, skip_swarp=False)
    threedhst.shifts.plot_shifts('/3DHST/Ancillary/GOODS-S/CANDELS/hlsp_candels_hst_wfc3_gsd01_f125w_v0.5_drz.fits', 'GOODS-deep-F125W_align.fits', WEIGHT_IMAGE='/3DHST/Ancillary/GOODS-S/CANDELS/hlsp_candels_hst_wfc3_gsd01_f125w_v0.5_wht.fits', drz=False, skip_swarp=False)

    threedhst.shifts.plot_shifts('GOODS-wide-F125W', ALIGN_IMAGE, skip_swarp=False)
    
    #### Check PSF
    files=glob.glob('GOODS-[wd]???*drz.fits')
    for file in files[3:]:
        unicorn.candels.easy_sextractor(drz_file=file)
    
    d105 = threedhst.sex.mySexCat('GOODS-deep-F105W_drz.cat')
    d125 = threedhst.sex.mySexCat('GOODS-deep-F125W_drz.cat')
    d160 = threedhst.sex.mySexCat('GOODS-deep-F160W_drz.cat')
    
    w105 = threedhst.sex.mySexCat('GOODS-wide-F105W_drz.cat')
    w125 = threedhst.sex.mySexCat('GOODS-wide-F125W_drz.cat')
    w160 = threedhst.sex.mySexCat('GOODS-wide-F160W_drz.cat')
    
    #plt.plot(d105['MAG_AUTO'], d105['FLUX_RADIUS'], marker='o', linestyle='None', color='blue', alpha=0.1)
    plt.plot(d125['MAG_AUTO'], d125['FLUX_RADIUS'], marker='o', linestyle='None', color='blue', alpha=0.1)
    plt.plot(d160['MAG_AUTO'], d160['FLUX_RADIUS'], marker='o', linestyle='None', color='red', alpha=0.1)
    plt.plot(w125['MAG_AUTO'], w125['FLUX_RADIUS'], marker='o', linestyle='None', color='purple', alpha=0.2)
    plt.plot(w160['MAG_AUTO'], w160['FLUX_RADIUS'], marker='o', linestyle='None', color='orange', alpha=0.2)
    
    plt.xlim(16,27)
    plt.ylim(0,20)
    
    ########################
    ### Wide epochs
    ########################
    ALIGN_IMAGE = '../ECDFS/MUSYC_ECDFS_BVR.fits'
    epochs = [wide1, wide2]
    
    for i,epoch in enumerate(epochs):
        e=i+1
        visits = np.unique(list.targname[epoch])
        for visit in visits:
            filters = np.unique(list.filter[list.targname == visit])
            for filter in filters:
                file = visit+'-'+filter+'_asn.fits'
                if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file, 
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                        SCALE=0.128254)
        ## Mosaic
        f125w_list = []
        f160w_list = []
        for visit in visits:
            f125w_list.append(visit+'-F125W_asn.fits')
            f160w_list.append(visit+'-F160W_asn.fits')
        #
        threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-wide%d-F125W' %(e), path_to_FLT='./', run_multidrizzle=False)
        threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-wide%d-F160W' %(e), path_to_FLT='./', run_multidrizzle=False)
        #
        SCALE = 0.128254
        for filt in ['F125W', 'F160W']:
            ra, dec, nx, ny = 53.165233, -27.900847, int(5350*0.128254/SCALE), int(3680*0.128254/SCALE)
            threedhst.prep_flt_files.startMultidrizzle( 'GOODS-wide%d-' %(e)+filt+'_asn.fits',
                use_shiftfile=True, skysub=False,
                final_scale=SCALE, pixfrac=1.0, driz_cr=False,
                updatewcs=False, clean=True, median=False,
                ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    #
    files=glob.glob('GOODS-W*drz.fits')
    files.extend(glob.glob('*align.fits'))
    for file in files: 
        os.remove(file)
    
    #### Full wide mosaic
    f125w_list = glob.glob('GOODS-W?-*F125W_asn.fits')
    f160w_list = glob.glob('GOODS-W?-*F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-wide-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-wide-F160W', path_to_FLT='./', run_multidrizzle=False)

    SCALE = 0.128254
    ra, dec, nx, ny = 53.165233, -27.900847, int(5670*0.128254/SCALE), int(5220*0.128254/SCALE)
    for filt in ['F125W', 'F160W']:
        threedhst.prep_flt_files.startMultidrizzle('GOODS-wide-'+filt+ '_asn.fits',
        use_shiftfile=True, skysub=False,
        final_scale=SCALE, pixfrac=1.0, driz_cr=False,
        updatewcs=False, clean=True, median=False,
        ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    
    ########################
    ###   ERS imaging
    ########################
    ALIGN_IMAGE='../ACS/h_sz*drz*fits'
    ers = list.group == 'WFC3-ERSII'
    
    visits = np.unique(list.targname[ers])
    for visit in visits:
        filters = np.unique(list.filter[list.targname == visit])
        for filter in filters:
            if filter.startswith('G'):
                continue
            file = visit+'-'+filter+'_asn.fits'
            if not os.path.exists(file.replace('asn','drz')):
                unicorn.candels.prep_candels(asn_file=file, 
                    ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.128254)
    ## Mosaic
    f098m_list = []
    f125w_list = []
    f160w_list = []
    for visit in visits:
        if not 'G01' in visit:
            f098m_list.append(visit+'-F098M_asn.fits')
            f125w_list.append(visit+'-F125W_asn.fits')
            f160w_list.append(visit+'-F160W_asn.fits')
    #
    threedhst.utils.combine_asn_shifts(f098m_list,out_root='GOODS-ERS-F098M', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f125w_list,out_root='GOODS-ERS-F125W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(f160w_list,out_root='GOODS-ERS-F160W', path_to_FLT='./', run_multidrizzle=False)
    #
    SCALE = 0.128254
    for filt in ['F098M','F125W', 'F160W']:
        ra, dec, nx, ny = 53.098662, -27.715892, int(5410*0.128254/SCALE), int(3460*0.128254/SCALE)
        threedhst.prep_flt_files.startMultidrizzle( 'GOODS-ERS-'+filt+ '_asn.fits',
            use_shiftfile=True, skysub=False,
            final_scale=SCALE, pixfrac=1.0, driz_cr=False,
            updatewcs=False, clean=True, median=False,
            ra=ra, dec=dec, final_outnx=nx, final_outny=ny)
    #
    files=glob.glob('WFC3-ERSII-IR*drz.fits')
    files.extend(glob.glob('*align.fits'))
    for file in files: 
        os.remove(file)
#
def ers():
    import unicorn.candels
    
    os.chdir('/Users/gbrammer/CANDELS/ERS/PREP_FLT/')
    #unicorn.candels.make_asn_files()
    
    ALIGN_IMAGE = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    
    filter = 'F125W'
    
    os.chdir('/Users/gbrammer/CANDELS/ERS/%s/' %(filter))
    files=glob.glob('WFC3*%s_asn.fits' %(filter))
    
    #files=glob.glob('GOODS*F125W_asn.fits')
    #files=glob.glob('GOODS*F105W_asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    #
    redo = ['WFC3-ERSII-IR02-F125W']
    
    for root in redo:
        threedhst.dq.checkDQ(root+'_asn.fits',root+'_asn.fits', size=900)
    
    for root in redo:
        unicorn.candels.prep_candels(asn_file=root+'_asn.fits', 
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=1,
            GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
    
def goodss():
    import unicorn.candels
    
    os.chdir('/Users/gbrammer/CANDELS/GOODS-S/PREP_FLT/')
    #unicorn.candels.make_asn_files()
    
    #ALIGN_IMAGE = 'AEGIS-N2_K_sci.fits'
    ALIGN_IMAGE = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    
    files=glob.glob('GOODS*F160W_asn.fits')
    #files=glob.glob('GOODS*F125W_asn.fits')
    #files=glob.glob('GOODS*F105W_asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    
    ### F125W
    redo_shifts = ['GOODS-S205-VI3-F125W']
    bad = ['GOODS-SD2-V7G-F125W', 'GOODS-SD5-VGQ-F125W', 'GOODS-SD5-VGX-F125W']
    redo = ['GOODS-S075-VJ3-F125W', 'GOODS-S080-VDV-F125W', 'GOODS-S080-VDX-F125W', 'GOODS-S080-VDY-F125W', 'GOODS-SD2-V7F-F125W', 'GOODS-SD2-V8J-F125W', 'GOODS-SD5-VG0-F125W', 'GOODS-SD5-VGU-F125W', 'GOODS-SDW-V02-F125W', 'GOODS-SDW-VH5-F125W', 'GOODS-SDW-VHB-F125W', 'GOODS-SDW-VHE-F125W', 'GOODS-W1-VIY-VIY-F125W', 'GOODS-W1-VJ0-VJ0-F125W', 'GOODS-W2-VIO-VIO-F125W']
    
    ### F160W
    redo_shifts = ['GOODS-S205-VI3-F160W']
    bad = ['GOODS-SD2-V7G-F160W', 'GOODS-SD5-VGQ-F160W', 'GOODS-SD5-VGX-F160W']
    redo = ['GOODS-S075-VJ3-F160W', 'GOODS-S075-VJ5-F160W', 'GOODS-S080-VDV-F160W', 'GOODS-S080-VDW-F160W', 'GOODS-S080-VDX-F160W', 'GOODS-S080-VDY-F160W', 'GOODS-S205-VHZ-F160W', 'GOODS-S205-VI7-F160W', 'GOODS-SD2-V7A-F160W', 'GOODS-SD2-V7D-F160W', 'GOODS-SD3-VEN-F160W', 'GOODS-SD5-V8Q-F160W', 'GOODS-SD5-VGC-F160W', 'GOODS-SD5-VGU-F160W', 'GOODS-SDW-V02-F160W', 'GOODS-SDW-VH4-F160W', 'GOODS-SDW-VHB-F160W', 'GOODS-SDW-VHD-F160W', 'GOODS-SDW-VHE-F160W', 'GOODS-W1-VIV-VIV-F160W', 'GOODS-WIDE115-V4D-F160W']
    
    for root in redo:
        threedhst.dq.checkDQ(root+'_asn.fits',root+'_asn.fits', size=900)
    
    for root in redo:
        unicorn.candels.prep_candels(asn_file=root+'_asn.fits', 
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=1,
            GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
    
    #########   XXXXXXXXXXXXXXXX todo: ERS fields
    
    os.chdir('/Users/gbrammer/CANDELS/GOODS-S/F105W')
    files=glob.glob('GOODS*asn.fits')
    threedhst.utils.combine_asn_shifts(files, out_root='GS-F105W', path_to_FLT='./', run_multidrizzle=False)
    #
    SCALE = 0.06
    threedhst.prep_flt_files.startMultidrizzle('GS-F105W_asn.fits',
         use_shiftfile=True, skysub=False,
         final_scale=SCALE, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False) 
    #
    os.chdir('/3DHST/Ancillary/GOODS-S/GOODS_ACS')
    for band in ['i','b','v','z']:
        threedhst.shifts.matchImagePixels(input= glob.glob('/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_s%s*drz_img.fits' %(band)), matchImage='/Volumes/Crucial/3DHST/Ancillary/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits', match_extension=0, output='GS-ACS%s.fits' %(band))
    
    ### Unicorn, match ACS to 0.06" WFC3 image
    os.chdir('/3DHST/Ancillary/GOODS-S/GOODS_ACS')
    for band in ['i','b','v','z'][:1]:
        threedhst.shifts.matchImagePixels(input= glob.glob('/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_s%s*drz_img.fits' %(band)), matchImage='../CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits', match_extension=0, output='GS-ACS%s.fits' %(band))
    
    ######## Big RGB map of the whole field
    ### scaling is related to the filter AB zeropoint
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    rgb = '/Volumes/Crucial/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits[0]*%.3f, /Volumes/Crucial/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F125W_wfc3ir_drz_sci.fits[0]*%.3f, /Volumes/Crucial/3DHST/Ancillary/GOODS-S/GS-ACSi.fits[0]*%.3f' %(scales[0], scales[1], scales[2])
    
    #### SWarp G141 images to put in same map
    ref = '/Volumes/Crucial/3DHST/Ancillary/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    ext = 0
    
    ref = '/3DHST/Spectra/Work/GOODS-S/PREP_FLT/GOODS-S-34-F140W_drz.fits'
    ext = 1
    
    #### Subimage around interesting sources
    #10820.215 6150.307
    #5633.7106 11078.678
    #6125.7472 8651.4943 
    imcopy('/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits[0][5000:11400,5500:11600]', 'GS_ref.fits')
    ref = 'GS_ref.fits'
    ext = 0
    
    os.chdir('/Users/gbrammer/Sites/GS_MAP/')
    threedhst.shifts.matchImagePixels( input=glob.glob('/Users/gbrammer/Sites/GS_MAP/GOODS-S-G141_drz_v0.5.fits'), matchImage=ref, match_extension=ext, input_extension=1, output='GS-G141.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/GOODS-S/GOODS_ACS/GS-ACSi.fits'), matchImage=ref, match_extension=ext, input_extension=0, output='GS-i.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F125W_wfc3ir_drz_sci.fits'), matchImage=ref, match_extension=ext, input_extension=0, output='GS-J.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits'), matchImage=ref, match_extension=ext, input_extension=0, output='GS-H.fits')
    
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    
    scales = np.array(scales)*2
    
    rgb = 'GS-H.fits[0]*%.3f, GS-J.fits[0]*%.3f, GS-i.fits[0]*%.3f' %(scales[0], scales[1], scales[2])
    
    threedhst.gmap.makeImageMap(['GS-i.fits[0]*%.3f' %(scales[2]), 'GS-H.fits[0]*%.3f' %(scales[0]), rgb, 'GS-G141.fits[0]*3'], aper_list=[14,16,17], tileroot=['ACS-i','WFC3-H', 'iJH', 'G141'], extension=1, path='/Users/gbrammer/Sites/GS_MAP/')
    
    #### Equirect image
    im_r = pyfits.open('/Volumes/Crucial/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits')
    im_g = pyfits.open('/Volumes/Crucial/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F125W_wfc3ir_drz_sci.fits')
    im_b = pyfits.open('/Volumes/Crucial/3DHST/Ancillary/GOODS-S/GS-ACSi.fits')
    xc, yc = 8193, 10055
    NX, NY = 6000, 6000
    #NX, NY = 500, 500
    #NX, NY = 1500, 800
    
    xc, yc = 16384/2, 20108/2
    
    ### F160W
    sub_r = im_r[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.96-25.96))
    ### F125W
    sub_g = im_g[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(26.25-25.96))
    ### F814W
    sub_b = im_b[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.94-25.96))*1.5
    
    Q, alpha, m0 = 10.,8.,-0.02
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='goodss-rgb.jpg', shape=(sub_r.shape[1]/2, sub_r.shape[0]/2))
    

def goodsn_prep(only_f140w=False, only_flat=False):
    
    import unicorn.candels
    import glob
    
    unicorn.candels.make_asn_files()

    #ALIGN_IMAGE = '/3DHST/Spectra/Work/GOODS-N/PREP_FLT/GOODS-N-F140W_drz.fits'
    ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f160w-astrodrizzle-v0.1_drz_sci.fits'

    bad_f125w = ['GOODSN-ORI090-V7P-F125W', 'GOODSN-ORI180-V1C-F125W', 'GOODSN-ORI180-V1F-F125W', 'GOODSN-ORI180-V1K-F125W', 'GOODSN-ORI180-V1N-F125W', 'GOODSN-ORI180-V1R-F125W']
    bad_f160w = ['GOODSN-ORI326-V5L-F160W','GOODSN-ORI090-V7P-F160W','GOODSN-ORI180-V1C-F160W', 'GOODSN-ORI180-V1F-F160W', 'GOODSN-ORI180-V1K-F160W', 'GOODSN-ORI180-V1N-F160W', 'GOODSN-ORI180-V1R-F160W']

    if not only_f140w and not only_flat:
        for filter, bad_files in zip(['F125W','F160W'],[bad_f125w,bad_f160w]):
            files=glob.glob('GOODSN-ORI*-V*'+filter+'_asn.fits')
            for bad_file in bad_files:
                print 'Will skip ', bad_file+'_asn.fits'
                files.remove(bad_file+'_asn.fits')
            threedhst.utils.combine_asn_shifts(files, out_root='GOODS-N-'+filter,
                path_to_FLT='./', run_multidrizzle=False)
            for file in files:
                if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file, 
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                        SCALE=0.06, geometry='rotate,shift')  

    if only_f140w:
        for filter in ['F140W']:
            files = glob.glob('GOODS-N-*-'+filter+'_asn.fits')
            threedhst.utils.combine_asn_shifts(files, out_root='GOODS-N-'+filter,
                path_to_FLT='./', run_multidrizzle=False)
            for file in files:
                if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file, 
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                        SCALE=0.06, geometry='rotate,shift')    
            
    if only_flat:
        for filter, bad_files in zip(['F125W','F160W'],[bad_f125w,bad_f160w]):
            files=glob.glob('GOODSN-ORI*-V*'+filter+'_asn.fits')
            for bad_file in bad_files:
                print 'Will skip ', bad_file+'_asn.fits'
                files.remove(bad_file+'_asn.fits')
            for file in files:
                unicorn.candels.prep_candels(asn_file=file, 
                    GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
                    redo_segmentation=False)    

def colfax_prep():
       
       import unicorn.candels
       import glob
    
       unicorn.candels.make_asn_files()
       
       ALIGN_IMAGE = '/Volumes/robot/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f160w-astrodrizzle-v0.1_drz_sci.fits'
       
       files = glob.glob('*asn.fits')
       
       for file in files:                
           if not os.path.exists(file.replace('asn','drz')):
                    unicorn.candels.prep_candels(asn_file=file,
                        ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                        SCALE=0.06, geometry='rotate,shift')       
    
def goodsn():
    import unicorn.candels
    
    os.chdir('/Users/gbrammer/CANDELS/GOODS-N/PREP_FLT/')
    unicorn.candels.make_asn_files()

    ### Some images fall outside of the 3D-HST area
    #ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_nz*drz*fits'   
    
    ALIGN_IMAGE = '/3DHST/Spectra/Work/GOODS-N/PREP_FLT/GOODS-N-F140W_drz.fits'
    
    files=glob.glob('GOOD*F160W*asn.fits')
    files=glob.glob('GOOD*F125W*asn.fits')
    #files=glob.glob('GOOD*F105W*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=1,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    
    ### Tiles outside of the F140W coverage
    filter='F160W'
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_nz*drz*fits'   
    acs_align = ['GOODSN-ORI326-V5Q', 'GOODSN-ORI326-V5R', 'GOODSN-ORI326-V5L', 'GOODSN-ORI326-V5M', 'GOODSN-ORI326-V5S', 'GOODSN-ORI326-V5V', 'GOODSN-ORI326-V5Y']
    for acs in acs_align:
            unicorn.candels.prep_candels(asn_file=acs+'-%s_asn.fits' %(filter), 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    
    ### F125W flagged
    bad = ['GOODSN-ORI180-V1C-F125W', 'GOODSN-ORI180-V1F-F125W', 'GOODSN-ORI180-V1K-F125W', 'GOODSN-ORI180-V1N-F125W', 'GOODSN-ORI180-V1R-F125W']
    redo = ['GOODSN-ORI180-V1X-F125W', 'GOODN-ORI090-V3M-F125W', 'GOODN-ORI090-V3M-F125W',  'GOODSN-ORI135-V2J-F125W', 'GOODSN-ORI135-V2P-F125W', 'GOODSN-ORI135-V2X-F125W', 'GOODSN-ORI135-V2Y-F125W', 'GOODSN-ORI180-V1D-F125W', 'GOODSN-ORI180-V1P-F125W', 'GOODSN-ORI180-V1U-F125W', 'GOODN-ORI090-V3O-F125W', 'GOODSN-ORI135-V2A-F125W', 'GOODSN-ORI135-V2C-F125W', 'GOODSN-ORI020-V4A-F125W', 'GOODSN-ORI020-V4C-F125W', 'GOODSN-ORI020-V4F-F125W', 'GOODSN-ORI020-V4H-F125W', 'GOODSN-ORI020-V4I-F125W', 'GOODSN-ORI020-V4J-F125W', 'GOODSN-ORI020-V4L-F125W', 'GOODSN-ORI020-V4W-F125W', 'GOODSN-ORI326-V5Y-F125W', 'GOODSN-ORI326-V5L-F125W', 'GOODSN-ORI326-V5O-F125W']
    
    ### F160W flagged pointings
    bad = ['GOODSN-ORI180-V1C-F160W', 'GOODSN-ORI180-V1F-F160W', 'GOODSN-ORI180-V1K-F160W', 'GOODSN-ORI180-V1N-F160W', 'GOODSN-ORI180-V1R-F160W']
    redo = ['GOODSN-ORI135-V2O-F160W', 'GOODSN-ORI135-V2Q-F160W', 'GOODSN-ORI135-V2U-F160W', 'GOODSN-ORI180-V1B-F160W', 'GOODSN-ORI180-V1I-F160W', 'GOODSN-ORI180-V1O-F160W', 'GOODSN-ORI180-V1Q-F160W', 'GOODSN-ORI180-V1S-F160W', 'GOODSN-ORI180-V1T-F160W', 'GOODSN-ORI180-V1V-F160W', 'GOODSN-ORI180-V1W-F160W', 'GOODSN-ORI180-V1X-F160W', 'GOODSN-ORI020-V4A-F160W', 'GOODSN-ORI020-V4W-F160W', 'GOODSN-ORI020-V4X-F160W', 'GOODSN-ORI326-V5F-F160W', 'GOODSN-ORI326-V5O-F160W', 'GOODSN-ORI326-V5T-F160W']
    
    for root in redo:
        threedhst.dq.checkDQ(root+'_asn.fits',root+'_asn.fits', size=900)
    
    for root in redo:
        unicorn.candels.prep_candels(asn_file=root+'_asn.fits', 
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=1,
            GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
    
    #### Mosaics
    filter = 'F160W'
    
    files, root = glob.glob('GOODS*020*asn.fits'), 'GN-020'
    
    files, root = glob.glob('GOOD*asn.fits'), 'GN'
    
    files, root = glob.glob('GOOD*asn.fits'), 'GN-v2'

    files, root = glob.glob('GOOD*asn.fits'), 'GN-v3'
    
    keep_list = []
    for file in files:
        if file.split('_asn.fits')[0] not in bad:
            keep_list.append(file)
            
    threedhst.utils.combine_asn_shifts(keep_list, out_root='%s-%s' %(root, filter), path_to_FLT='./', run_multidrizzle=False)
    SCALE = 0.06
    #### Match 3D-HST F140W
    NX, NY = int(6840*0.128254/SCALE), int(8042*0.128254/SCALE)
    ra0, dec0 = 189.17736, 62.23892
    #### Full CANDELS Mosaic
    NX, NY = 18500, 18500
    ra0, dec0 = 189.236, 62.246087
    threedhst.prep_flt_files.startMultidrizzle('%s-%s_asn.fits' %(root, filter),
         use_shiftfile=True, skysub=False,
         final_scale=SCALE, pixfrac=0.8, driz_cr=False,
         updatewcs=False, clean=True, median=False,
         ra=ra0, dec=dec0,
         final_outnx = NX, final_outny=NY, ivar_weights=True, build_drz=False) 
    
    #### Matched ACS images
    os.chdir('/Users/gbrammer/CANDELS/GOODS-N/ACS_MATCH')
    for band in ['i','b','v','z'][1:]:
        threedhst.shifts.matchImagePixels(input= glob.glob('/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_n%s*drz_img.fits' %(band)), matchImage='/3DHST/Spectra/Work/GOODS-N/PREP_FLT/GOODS-N-F140W_drz.fits', match_extension=1, output='GN-ACS%s.fits' %(band))

    #
    for band in ['i','b','v','z'][1:]:
        threedhst.shifts.matchImagePixels(input= glob.glob('/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_n%s*drz_img.fits' %(band)), matchImage='/Users/gbrammer/CANDELS/GOODS-N/F160W/GN-v3-F160W_drz_sci.fits', match_extension=0, output='GN-ACS%s_v3.fits' %(band))
    
    ### Unicorn
    os.chdir('/Volumes/robot/3DHST/Ancillary/GOODS-N/GOODS_ACS')
    for band in ['i','b','v','z'][:1]:
        threedhst.shifts.matchImagePixels(input= glob.glob('h_n%s*drz_img.fits' %(band)), matchImage='../3DHST_F140W/GOODS-N-F140W_drz_v0.5.fits', match_extension=1, output='GN-ACS%s.fits' %(band))
    
    #### Check for rejected stars
    filter = 'F160W'
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()

    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['FILTER']    = 'Y'
    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '3'
    se.options['ANALYSIS_THRESH']  = '3' 

    #### Run SExtractor
    threedhst.options['MAG_ZEROPOINT'] = '25.96'
    se.options['CATALOG_NAME']    = 'GN-%s.cat' %(filter)
    se.options['CHECKIMAGE_NAME'] = 'GN-%s_seg.fits' %(filter)
    se.options['WEIGHT_IMAGE']    = 'GN-%s_drz_weight.fits' %(filter)
    status = se.sextractImage('GN-v2-%s_drz_sci.fits' %(filter))
    
    ZPs = {}
    ZPs['F160W'] = 25.96
    ZPs['F125W'] = 26.25
    
    cat = threedhst.sex.mySexCat('GN-%s.cat' %(filter))
    r50 = np.cast[float](cat.FLUX_RADIUS)
    mag = np.cast[float](cat.MAG_AUTO)+ZPs[filter]
    stars = (mag < 20) & (r50 < 3.5) #& (mag > 18)
    
    #wht = pyfits.open('GN-F160W_drz_weight.fits')
    wht = pyfits.open('GN-v2-%s_drz_weight.fits' %(filter))
    NX = 25
    im = np.zeros((2*NX, 2*NX))
    x = np.cast[int](np.round(cat['X_IMAGE'][stars]))
    y = np.cast[int](np.round(cat['Y_IMAGE'][stars]))
    for i in range(stars.sum()):
        print unicorn.noNewLine + '%d of %d' %(i+1, stars.sum())
        sub = wht[0].data[y[i]-NX:y[i]+NX, x[i]-NX:x[i]+NX]
        im += sub / np.median(sub)
    
    pyfits.writeto('star.fits', data=im, clobber=True)
        
def make_asn_files(force=False, make_region=False, uniquename=True):
    """
    Read a files.info file and make ASN files for each visit/filter[/date].
    """
    list = catIO.Readfile('files.info')
    list.pa_v3 = np.cast[int](np.round(list.pa_v3))
    asn = threedhst.utils.ASNFile(glob.glob('../RAW/i*asn.fits')[0])
    
    dates = np.array([''.join(date.split('-')[1:]) for date in list.date_obs])
    targets = np.unique(list.targname)
    
    visits = np.array([file[4:6] for file in list.file])
    for visit in np.unique(visits):
        angle = list.pa_v3[visits == visit][0]
        for target in np.unique(list.targname[visits == visit]):
            #### 3D-HST targname translations
            translate = {'AEGIS-':'aegis-', 'COSMOS-':'cosmos-', 'GNGRISM':'goodsn-', 'GOODS-SOUTH-':'goodss-', 'UDS-':'uds-'}
            target_use = target
            for key in translate.keys():
                target_use = target_use.replace(key, translate[key])
            ## pad i < 10 with zero
            for key in translate.keys():
                if translate[key] in target_use:
                    spl = target_use.split('-')
                    if int(spl[-1]) < 10:
                        spl[-1] = '%02d' %(int(spl[-1]))
                        target_use = '-'.join(spl)
            #
            for filter in np.unique(list.filter[(list.targname == target) & (visits == visit)]):
                if uniquename:
                    product='%s-%s-%03d-%s' %(target_use, visit.upper(), int(angle), filter)
                else:
                    product='%s-%s' %(target_use, filter)             
                #       
                use = (list.targname == target) & (list.filter == filter) & (visits == visit)
                exposure_list = []
                for file in list.file[use]:
                    f = file.split('.gz')[0]
                    if f not in exposure_list:
                        exposure_list.append(f)
                
                print product, len(exposure_list)
                        
                if (os.path.exists('%s_asn.fits' %(product))) & (not force):
                    continue
                #
                asn.exposures = []
                asn.product=product
                for fits in exposure_list:
                    root = os.path.basename(fits).split('_flt')[0]
                    asn.exposures.append(root)
                #
                asn.write(asn.product+'_asn.fits')
                if make_region:
                    threedhst.regions.asn_region(asn.product+'_asn.fits', path_to_flt='../RAW')
                    
    # for target in targets:        
    #     filters = np.unique(list.filter[list.targname == target])
    #     print target, filters
    #     for filter in filters:
    #         angles = np.unique(list.pa_v3[(list.targname == target) & (list.filter == filter)])
    #         for angle in angles:
    #             if (os.path.exists('%s-%03d-%s_asn.fits' %(target, int(angle), filter))) & (not force):
    #                 continue
    # 
    #             #print asn.in_fits[1].columns
    #             use = (list.targname == target) & (list.filter == filter) & (list.pa_v3 == angle)
    #             asn.exposures = []
    #             asn.product='%s-%03d-%s' %(target, int(angle), filter)
    #             for fits in list.file[use]:
    #                 asn.exposures.append(os.path.basename(fits).split('_flt')[0])
    # 
    #             #try:
    #             asn.write(asn.product+'_asn.fits')
    #             #except:
    #             #    continue
    # 
    #             threedhst.regions.asn_region(asn.product+'_asn.fits', path_to_flt='../RAW')
    #             print asn.product
            
    
def prep_candels(asn_file='ib3706050_asn.fits',
                       ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
                       ALIGN_EXTENSION=0,
                       PATH_TO_RAW='../RAW',
                       GET_SHIFT=True,
                       DIRECT_HIGHER_ORDER=2,
                       SCALE=0.06,
                       bg_skip=False,
                       geometry='rxyscale,shift',
                       clean=True,
                       redo_segmentation=True,
                       clean_zeros=True):
    
    import threedhst
    import threedhst.prep_flt_files
    
    threedhst.process_grism.fresh_flt_files(asn_file, from_path=PATH_TO_RAW)

    ##### Make region files for the pointing
    if not os.path.exists(asn_file.replace('fits','pointing.reg')):
        threedhst.regions.asn_region(asn_file)
    
    threedhst.prep_flt_files.prep_flt(asn_file=asn_file,
                    get_shift=GET_SHIFT, 
                    bg_only=False, bg_skip=bg_skip, redo_background=True,
                    ALIGN_IMAGE=ALIGN_IMAGE, 
                    ALIGN_EXT=ALIGN_EXTENSION,
                    skip_drz=False, final_scale=SCALE, pixfrac=0.8,
                    IMAGES=[],
                    align_geometry=geometry, clean=clean,
                    initial_order=0, save_fit=False, TWEAKSHIFTS_ONLY=(ALIGN_IMAGE is None), redo_segmentation=redo_segmentation)
    
    if DIRECT_HIGHER_ORDER > 0:
        threedhst.prep_flt_files.prep_flt(asn_file=asn_file,
                    get_shift=False, first_run=False, 
                    bg_only=False, bg_skip=bg_skip, redo_background=False,
                    skip_drz=False, final_scale=SCALE, pixfrac=0.8,
                    IMAGES=[], clean=clean,
                    initial_order=DIRECT_HIGHER_ORDER, save_fit=False,
                    redo_segmentation=redo_segmentation)
    
    #### Set pixels with zero weight to zero value
    if clean_zeros:
        print 'Clean pixels with zero weight.'
        drz = pyfits.open(asn_file.replace('asn','drz'), mode='update')
        drz[1].data[drz[2].data == 0] = 0
        drz.flush()
        
def make_test_catalog():
    """
    Run SExtractor on some CANDELS images to see what the mag vs. size relation looks like
    """
    import threedhst
    os.chdir('/Users/gbrammer/CANDELS/JUNK')
    
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()

    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['FILTER']    = 'Y'
    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '3'
    se.options['ANALYSIS_THRESH']  = '3' 

    #### Run SExtractor
    threedhst.options['MAG_ZEROPOINT'] = '26.25'
    se.options['CATALOG_NAME']    = 'UDS01_F125W_drz.cat'
    se.options['CHECKIMAGE_NAME'] = 'UDS01_F125W_seg.fits'
    se.options['WEIGHT_IMAGE']    = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f125w_v0.5_wht.fits'
    status = se.sextractImage('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f125w_v0.5_drz.fits')

    threedhst.options['MAG_ZEROPOINT'] = '25.96'
    se.options['CATALOG_NAME']    = 'UDS01_F160W_drz.cat'
    se.options['CHECKIMAGE_NAME'] = 'UDS01_F160W_seg.fits'
    se.options['WEIGHT_IMAGE']    = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_wht.fits'
    status = se.sextractImage('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_drz.fits')
    
    ######### Locally
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/FIRST_PAPER/GRISM_v1.6')
    uds125 = threedhst.sex.mySexCat('UDS01_F125W_drz.cat')
    uds160 = threedhst.sex.mySexCat('UDS01_F160W_drz.cat')
    threed = catIO.Readfile('full_sextractor.cat')
    
    fig = unicorn.catalogs.plot_init(xs=4,aspect=2)
    
    #fig = plt.figure(figsize=(5.5,5), dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.15,
                        bottom=0.10,right=0.99,top=0.97)

    ## zphot
    ax = fig.add_subplot(311)   
    plt.plot(np.cast[float](uds125.MAG_AUTO)+26.25, np.cast[float](uds125.FLUX_RADIUS), marker='o', alpha=0.2, linestyle='None', markersize=3, color='blue')
    plt.plot([0,100],[2.5,2.5], alpha=0.2, color='black', linewidth=10)
    plt.xlim(14.1,25)
    plt.ylim(0,11)
    plt.text(15,9,'UDS F125W')
    
    ax = fig.add_subplot(312)
    plt.plot(np.cast[float](uds160.MAG_AUTO)+26.25, np.cast[float](uds160.FLUX_RADIUS), marker='o', alpha=0.2, linestyle='None', markersize=3, color='blue')
    plt.plot([0,100],[2.5,2.5], alpha=0.2, color='black', linewidth=10)
    plt.xlim(14.1,25)
    plt.ylim(0,11)
    plt.text(15,9,'UDS F160W')
    #plt.close()
    
    ax = fig.add_subplot(313)
    plt.plot(th.mag_f1392w, th.flux_radius, marker='o', color='red', alpha=0.1, linestyle='None', markersize=3)
    plt.plot([0,100],[2.5,2.5], alpha=0.2, color='black', linewidth=10)
    plt.xlim(14.1,25)
    plt.ylim(0,11)
    plt.xlabel('MAG_AUTO')
    plt.ylabel('FLUX_RADIUS')
    plt.text(15,9,'3D-HST F140W')
    
    plt.savefig('size_candels_3dhst.png')
    
    plt.close()

def easy_sextractor(drz_file='UDS-F125W_drz.fits', threshold=3, zeropoint=26.25, auto_zp=True, mode='direct'):
    
    header = pyfits.getheader(drz_file, 0)
    filter = header['FILTER']

    if auto_zp:
        if filter == 'F105W':
            zeropoint = 26.27
        if filter == 'F140W':
            zeropoint = 26.46
        if filter == 'F125W':
            zeropoint = 26.25
        if filter == 'F160W':
            zeropoint = 25.96
    
    ROOT_GRISM = drz_file.split('_drz.fits')[0]
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = drz_file+'[1]'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '%.2f' %(threshold) 
    se.options['ANALYSIS_THRESH']  = '%.2f' %(threshold) 
    se.options['MAG_ZEROPOINT'] = '%.2f' %(zeropoint)
    
    status = se.sextractImage(drz_file+'[0]', mode=mode)
    
def flt_rejected_stars():
    pass
    
def star_stacks():
    """ 
    Make stacks of the images and weights around stars, 3DHST and CANDELS
    """
    import unicorn
    import threedhst.catIO as catIO
    import threedhst
    import pyfits
    import os
    import numpy as np
    import glob
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    
    NX = 50
    phot = catIO.Readfile('full_sextractor.cat')
    uds125 = threedhst.sex.mySexCat('UDS01_F125W_drz.cat')
    uds160 = threedhst.sex.mySexCat('UDS01_F160W_drz.cat')
    
    stars_125 = (np.cast[float](uds125.MAG_AUTO)+26.25 > 18) & (np.cast[float](uds125.MAG_AUTO)+26.25 < 19) & (np.cast[float](uds125.FLUX_RADIUS) < 3)
    stars_160 = (np.cast[float](uds160.MAG_AUTO)+26.25 > 18) & (np.cast[float](uds160.MAG_AUTO)+26.25 < 19) & (np.cast[float](uds160.FLUX_RADIUS) < 3)
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    
    # d_125 = pyfits.open('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f125w_v0.5_drz.fits')
    #d_125 = pyfits.open('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f125w_v0.5_wht.fits')

    #d_125 = pyfits.open('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_drz.fits')
    d_125 = pyfits.open('/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_wht.fits')
    
    #### Have to do wht, drz separately for memory issues
    
    data_125 = d_125[0].data
    
    drz_125 = np.zeros((2*NX,2*NX))
    wht_125 = np.zeros((2*NX,2*NX))
    
    idx = np.arange(len(stars_125))
    count=0
    for i in idx[stars_125]:
        print str(i)
        xc = np.round(float(uds125.X_IMAGE[i]))
        yc = np.round(float(uds125.Y_IMAGE[i]))
        drz_125 += data_125[yc-NX:yc+NX,xc-NX:xc+NX]
        count+=1
        #wht_125 += w_125_data[yc-NX:yc+NX,xc-NX:xc+NX]
    
    pyfits.writeto('star_drz_160.fits', data=drz_125/count, clobber=True)
    
    pyfits.writeto('star_wht_160.fits', data=drz_125/count, clobber=True)
    
    ####### Threedhst
    stars_3d = (phot.mag_f1392w > 18) & (phot.mag_f1392w < 19) & (phot.flux_radius < 3) & ((phot.field == 'AEGIS') | (phot.field == 'COSMOS') | (phot.field == 'GOODS-S') | (phot.field == 'GOODS-N'))
    idx = np.arange(len(stars_3d))
    old_pointing = 'xxx'
    
    drz_125 = np.zeros((2*NX,2*NX))
    wht_125 = np.zeros((2*NX,2*NX))
    
    count=0.
    for obj in phot.id[stars_3d]:
        path = unicorn.analysis.get_grism_path(obj)
        pointing = obj.split('G141')[0]+'F140W'
        print unicorn.noNewLine+pointing
        if pointing != old_pointing:
            drz = pyfits.open(path+'PREP_FLT/'+pointing+'_drz.fits')
            old_pointing = pointing
        #
        i = phot.id == obj
        xc = np.round(float(phot.x_image[i]))
        yc = np.round(float(phot.y_image[i]))
        #
        try:
            drz_125 += drz[1].data[yc-NX:yc+NX,xc-NX:xc+NX]
            wht_125 += drz[2].data[yc-NX:yc+NX,xc-NX:xc+NX]
            count+=1
        except:
            pass
        
    #
    pyfits.writeto('star_drz_3d.fits', data=drz_125/count, clobber=True)
    pyfits.writeto('star_wht_3d.fits', data=wht_125/count, clobber=True)
    
def check_drz_dq():
    import time
    import threedhst.dq
    fp = open('dq_check.txt','a')
    fp.write('# %s\n' %(time.ctime()))
    
    files=glob.glob('*drz.fits')
    ds9 = threedhst.dq.myDS9()
    
    for file in files:
        im= pyfits.open(file)
        ds9.view(im[1])
        ds9.scale(-0.02,0.05)
        ds9.set('zoom to fit')
        im.close()
        #
        inp = raw_input('%s 1/[0]: ' %(file))
        if inp != '':
            fp.write(file+'\n')
    
    fp.close()
    
def flag_dq_problems():
    import threedhst.dq
    bad_images = np.unique(np.loadtxt('dq_check.txt', dtype=np.str))
    for bad in bad_images:
        asn = bad.replace('_drz','_asn')
        threedhst.dq.checkDQ(asn_direct_file=asn, asn_grism_file=asn, path_to_flt='../RAW/')
    
def redo_processing_for_bad_dq():
    import threedhst.dq
    bad_images = np.unique(np.loadtxt('dq_check.txt', dtype=np.str))
    for bad in bad_images:
        asn = bad.replace('_drz','_asn')
        unicorn.candels.prep_candels(asn_file=asn, 
            ALIGN_IMAGE = None, ALIGN_EXTENSION=0,
            GET_SHIFT=False, DIRECT_HIGHER_ORDER=2,
            SCALE=0.06, geometry='rxyscale,shift')

#########################  
########## Things for the CLASH survey
#########################

def clash_make_rms_map(image='hlsp_clash_hst_wfc3ir_macs1206_total_v1_wht.fits', include_poisson=False):
    """
    Make an RMS image by simply taking the inverse of an inverse variance image.
    """
    wht = pyfits.open(image)
    out = image.replace('wht','rms')
    rms = 1./np.sqrt(wht[0].data)
    if include_poisson:
        drz = pyfits.open(image.replace('wht','drz'))
        texp = drz[0].header['EXPTIME']
        var = rms**2*texp**2 + drz[0].data*texp
        rms = np.sqrt(var)/texp
    
    wht[0].data = rms
    wht.writeto(out, clobber=True)

#
def clash_all_clusters():
    import unicorn
    
    unicorn.candels.clash_run_SExtractor(cluster='a2261')
    unicorn.candels.clash_run_SExtractor(cluster='a383')
    unicorn.candels.clash_run_SExtractor(cluster='macs1149')
    unicorn.candels.clash_run_SExtractor(cluster='macs1206')
    unicorn.candels.clash_run_SExtractor(cluster='macs2129')
    
    ### Redo catalogs
    for cluster in ['a2261','a383','macs1149','macs1206','macs2129']:
        os.chdir('/Users/gbrammer/CLASH/%s/CAT/' %(cluster))
        unicorn.candels.clash_make_catalog()
        
def clash_run_SExtractor(cluster='macs2129'):
    
    import threedhst
    
    try:
        os.mkdir('/Users/gbrammer/CLASH/%s/CAT' %(cluster))
    except:
        pass
        
    os.chdir('/Users/gbrammer/CLASH/%s/CAT' %(cluster))
    
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    threedhst.sex.USE_CONVFILE = 'gauss_2.0_5x5.conv'
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()

    #cluster = 'macs2129'
    
    detect_image = glob.glob('../hlsp*hst_wfc3ir*%s*total*drz.fits' %(cluster))[-1]
        
    se.overwrite = True
    #se.options['CATALOG_TYPE'] = 'FITS_1.0'
    
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '2.5'
    se.options['ANALYSIS_THRESH']  = '2.5' 
    
    se.options['DETECT_MINAREA']    = '8'
    se.options['DEBLEND_NTHRESH']    = '64'
    se.options['DEBLEND_MINCONT']    = '0.00001'
    se.options['DEBLEND_MINCONT']    = '0.000002'
    se.options['GAIN']    = '1.0'
    
    se.options['MEMORY_OBJSTACK'] = '8000'
    se.options['MEMORY_PIXSTACK'] = '800000'
    
    scale = 0.065
    apertures = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,2.5,3])/scale
    pix = '%.2f' %(apertures[0])
    for ap in apertures[1:]:
        pix += ',%.2f' %(ap)
        
    se.options['PHOT_APERTURES'] = pix
    #### Run SExtractor
    se.options['MAG_ZEROPOINT'] = '25'
    se.options['CATALOG_NAME']    = cluster+'_total_drz.cat'
    se.options['CHECKIMAGE_NAME'] = cluster+'_total_seg.fits'
    # se.options['WEIGHT_IMAGE']    = detect_image.replace('drz','wht')
    # se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    
    if not os.path.exists(detect_image.replace('drz','rms')):
        unicorn.candels.clash_make_rms_map(image=detect_image.replace('drz','wht'), include_poisson=False)
    se.options['WEIGHT_IMAGE']    = detect_image.replace('drz','rms')
    se.options['WEIGHT_TYPE']     = 'MAP_RMS'
    se.options['WEIGHT_GAIN']     = 'N'
    
    status = se.sextractImage(detect_image, mode='direct')
    #threedhst.sex.sexcatRegions(cluster+'_total_drz_xx.cat', cluster+'_total_drz_xx.reg', format=2)
    
    
    all_bands = glob.glob('../*drz.fits')
    for band in all_bands:
        se.options['CATALOG_NAME'] = os.path.basename(band).replace('fits','cat')
        se.options['WEIGHT_IMAGE']= '%s,%s' %(detect_image.replace('drz','rms'), band.replace('drz','wht'))
        se.options['WEIGHT_TYPE']     = 'MAP_RMS,MAP_WEIGHT'
        #
        fp = open('auto.sex','w')
        fp.write(se._makeOptionStr())
        fp.close()
        os.system('sex %s %s -c auto.sex' %(detect_image, band))
    #
    base = os.path.basename(detect_image).split('.fits')[0]
    threedhst.sex.sexcatRegions(base+'.cat', base+'.reg', format=2)
    
    unicorn.candels.clash_make_catalog()
    
def clash_get_coords():
    for cluster in ['a2261','a383','macs1149','macs1206','macs2129']:
        os.chdir('/Users/gbrammer/CLASH/%s/CAT/' %(cluster))
        catfile = glob.glob('*eazy.cat')[-1]
        cat = catIO.Readfile(catfile)
        print '%s %.6f %.6f' %(cluster, np.mean(cat.ra), np.mean(cat.dec))
        
def clash_make_catalog():
    
    files=glob.glob('*_f*_*cat')
    NFILT = len(files)
    filters = []
    fluxes = []
    errors = []
    
    flux_column = 'APER3' # 0.4" diameter 
    
    total_file = glob.glob('*_wfc3ir*total*cat')[0]
    tcat = threedhst.sex.mySexCat(total_file)
    ftot = tcat['FLUX_AUTO']/tcat['FLUX_'+flux_column]
    
    cluster = os.getcwd().split('/')[-2]
    ebv = {}
    ebv['a2261'] = 0.0430
    ebv['a383'] = 0.0309
    ebv['macs1149'] = 0.0227
    ebv['macs1206'] = 0.0616
    ebv['macs2129'] = 0.0793
    
    RES = eazy.FilterFile('/usr/local/share/EAZY/eazy-photoz/filters/FILTER.RES.latest')
    
    dered_str = '# E(B-V) = %.4f\n' %(ebv[cluster])
    
    for file in files:
        cat = threedhst.sex.mySexCat(file)
        filter=file.split('_')[5]
        filters.append(filter)
        #
        idx = RES.search(filter, verbose=False)[-1]
        try:
            dered = RES.filters[idx].extinction_correction(ebv[cluster], mag=False)
            dered_mag = RES.filters[idx].extinction_correction(ebv[cluster], mag=True)
        except:
            print '\n Problem with extinction correction, %s' %(filter)
            dered = 1.0
            dered_mag = -1
        #
        dered_str += '# %s %.3f\n' %(filter, dered_mag)
        head = pyfits.getheader('../%s' %(file.replace('cat','fits')))
        zp=-2.5*np.log10(head['PHOTFLAM']) - 21.10 - 5 *np.log10(head['PHOTPLAM']) + 18.6921
        #
        flux = np.cast[float](cat['FLUX_'+flux_column])*10**(-0.4*(zp-25))*ftot
        if os.path.exists('../%s_empty.fits' %(file.replace('.cat',''))):
            empty = pyfits.open('../%s_empty.fits' %(file.replace('.cat','')))
            texp = head['EXPTIME']
            error_color2 = threedhst.utils.biweight(empty[2].data[:,0])**2*texp**2+np.cast[float](cat['FLUX_'+flux_column])*texp
            error = np.sqrt(error_color2)/texp*10**(-0.4*(zp-25))*ftot
            depth = ' %.2f' %(zp-2.5*np.log10(5*threedhst.utils.biweight(empty[2].data[:,0])))
        else:
            error = np.cast[float](cat['FLUXERR_'+flux_column])*10**(-0.4*(zp-25))*ftot
            depth = ''
        #
        fluxes.append(flux*dered)
        errors.append(error*dered)
        print unicorn.noNewLine+file+depth
        
    fluxes = np.array(fluxes)
    errors = np.array(errors)
    
    bad = fluxes == 0
    fluxes[bad] = -99
    errors[bad] = -99
    
    bad = (ftot < 0) #| (ftot > 20)
    fluxes[:,bad] = -99
    errors[:,bad] = -99
    
    NFILT, NOBJ = fluxes.shape
    for i in range(NOBJ):
        sn = fluxes[:,i]/errors[:,i]
        ok = sn > 1
        if len(ok[ok]) < 3:
            fluxes[:,i] = -99
    
    #### Write the ascii catalog
    id, ra, dec = tcat['NUMBER'], tcat['X_WORLD'], tcat['Y_WORLD']
    
    header = '# id ra dec'
    for filt in filters:
        header = '%s f_%s e_%s' %(header, filt, filt)
    
    
    fp = open(total_file.split('_total')[0]+'_eazy.cat','w')
    fp.write(header+'\n')
    fp.write(dered_str)
    for i in range(NOBJ):
        print unicorn.noNewLine+'%d' %(i)
        flux = fluxes[:,i]
        error = errors[:,i]
        line = '%-6d %.6f %.6f' %(id[i], ra[i], dec[i])
        for j in range(NFILT):
            line = '%s %.3e %.3e' %(line, flux[j], error[j])
        #    
        fp.write(line+'\n')
    
    fp.close()
    
    if (1 == 0):
        for filt in filters:
            print 'f_%-6s Fxx' %(filt)
            print 'e_%-6s Exx' %(filt)
        #
        cat = catIO.Readfile(glob.glob('*eazy.cat')[0])
        zout = catIO.Readfile('OUTPUT/photz.zout')
        hmag = 25-2.5*np.log10(cat.f_f160w)
        
        hlim = 25
        plt.rcParams['text.usetex'] = False
        plt.hist(zout.z_peak[(hmag < hlim) & (zout.q_z < 5)], range=(0,5), bins=400)
        
        test = (hmag < hlim) & (zout.z_peak > 1.5) & (zout.q_z < 5)
        idx = np.arange(len(hmag))
        for i in idx[test]:
            print '%-6d %.3f %.2f %.6f %.6f' %(cat.id[i], zout.z_peak[i], hmag[i], cat.ra[i], cat.dec[i])
            
            
def clash_red_galaxies():
    import threedhst
    jcat = threedhst.sex.mySexCat('hlsp_clash_hst_wfc3ir_macs2129_f125w_v1_drz.cat')
    hcat = threedhst.sex.mySexCat('hlsp_clash_hst_wfc3ir_macs2129_f160w_v1_drz.cat')
    
    jflux = jcat['FLUX_APER4']*10**(-0.4*(25-26.25))
    hflux = hcat['FLUX_APER4']*10**(-0.4*(25-25.96))
    
    jhmag = -2.5*np.log10(jflux/hflux)
    
    plt.plot(hcat['MAG_AUTO'], jhmag, marker='o', linestyle='none', color='blue',alpha=0.5)
    plt.xlim(20,30)
    plt.ylim(-0.4,3)
    
    #
    drz = pyfits.open('../hlsp_clash_hst_wfc3ir_macs2129_f125w_v1_drz.fits')
    zp=-2.5*np.log10(drz[0].header['PHOTFLAM']) - 21.10 - 5 *np.log10(drz[0].header['PHOTPLAM']) + 18.6921

def get_full_candels_region(field='COSMOS'):
    import unicorn
    
    if field == 'AEGIS':
        field = 'EGS'
        
    reg_files = glob.glob(os.getenv('CANDELS')+'/%s/PREP_FLT/*F160W*pointing.reg' %(field))
    #out_file = os.getenv('CANDELS')+'/COSMOS_full.reg'
    
    regions = []
    for file in reg_files:
        print unicorn.noNewLine+file
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        for line in lines:
            if 'polygon' in line:
                regions.append(unicorn.survey_paper.polysplit(region=line[:-1], get_shapely=True))
    #
    un = regions[0]
    for pp in regions[1:]:
        un = un.union(pp)
    
    #
    if un.geometryType() is 'MultiPolygon':
        x, y = [], []
        for sub_poly in un.geoms:
            xi, yi = sub_poly.exterior.xy
            x.append(np.array(xi)), y.append(np.array(yi))
            #ax.plot(x,y, alpha=candels_alpha, color=candels_color, linewidth=1)
            #ax.fill(x,y, alpha=0.1, color='0.7')
    else:        
        x, y = un.exterior.xy
        x, y = np.array(x), np.array(y)
        #plt.plot(x,y, alpha=0.5, color='blue', linewidth=1)
        #plt.fill(x,y, alpha=0.1, color='0.7')
    
    return un, x, y
    
    
        
def massive_galaxies_morphologies(field='COSMOS', masslim=11.2, zlim=(1.7,6), save_fits=False, id=None, thumb_size = 8):
    import threedhst
    import threedhst.eazyPy as eazy
    import unicorn
    import shapely
    from shapely.geometry import Point,Polygon
    try:
        import stsci.tools.wcsutil as wcsutil
    except:
        import pytools.wcsutil as wcsutil
        
    cat, zout, fout = unicorn.analysis.read_catalogs(root='%s-1' %(field))
    try:
        test = cat.use
    except:
        cat.use = np.ones(len(cat.id), dtype=np.int)
    
    full_geom, xreg, yreg = unicorn.candels.get_full_candels_region(field=field)
    
    in_candels = cat.id < 0
    print 'Find objects within CANDELS'

    idx = np.arange(len(cat.id))
    for i in idx:
        #print unicorn.noNewLine+'%d' %(i)
        in_candels[i] = full_geom.contains(Point(cat.ra[i], cat.dec[i]))

    print unicorn.noNewLine+'Find objects within CANDELS - Done.'
        
    #
    eazy_path = {}
    eazy_path['AEGIS'] = '/3DHST/Ancillary/AEGIS/NMBS/Photometry//aegis-n2.deblend.redshifts'
    eazy_path['COSMOS'] = '/3DHST/Ancillary/COSMOS/NMBS/Photometry//cosmos-1.deblend.redshifts'
    eazy_path['GOODS-S'] = '/3DHST/Ancillary/GOODS-S/FIREWORKS/FAST'
    eazy_path['GOODS-N'] = '/3DHST/Ancillary/GOODS-N/MODS/Photometry/FAST/../EAZY/OUTPUT'
    eazy_path['UDS'] = '/3DHST/Ancillary/UDS/UKIDSS/FAST'
    
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    OUTPUT_DIRECTORY = eazy_path[field]
    CACHE_FILE='Same'
    
    F160W = {}
    F160W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F160W_drz_sci.fits'
    F160W['AEGIS'] = '/Users/gbrammer/CANDELS/EGS/PREP_FLT/EGS-F160W_drz_sci.fits'
    F160W['UDS'] = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'
    F160W['GOODS-S'] = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    
    F125W = {}
    F125W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F125W_drz_sci.fits'
    F125W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F125W_drz_sci.fits'
    F125W['AEGIS'] = '/Users/gbrammer/CANDELS/EGS/PREP_FLT/EGS-F125W_drz_sci.fits'
    F125W['UDS'] = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f125w_v1.0_drz.fits'
    F125W['GOODS-S'] = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F125W_wfc3ir_drz_sci.fits'
    
    im_f125w = pyfits.open(F125W[field])
    im_f160w = pyfits.open(F160W[field])
    
    wcs_f125w = wcsutil.WCSObject('f125w',header=im_f125w[0].header)
    wcs_f160w = wcsutil.WCSObject('f160w',header=im_f160w[0].header)
        
    # thumb_size: arcsec
    ps_f125w = np.sqrt(wcs_f125w.cd11**2+wcs_f125w.cd12**2)*3600.
    ps_f160w = np.sqrt(wcs_f160w.cd11**2+wcs_f160w.cd12**2)*3600.
    
    N_125 = int(np.round(thumb_size / ps_f125w / 2))
    N_160 = int(np.round(thumb_size / ps_f160w / 2))
    ticks_125 = np.arange(0,thumb_size+1)/ps_f125w
    ticks_160 = np.arange(0,thumb_size+1)/ps_f160w
    
    plt.gray()
    
    massive = in_candels & (cat.use == 1) & (fout.lmass > masslim) & (zout.z_peak > zlim[0]) & (zout.z_peak < zlim[1])
    
    if id is not None:
        massive = cat.id == id
        
    for ii in idx[massive]:
        print ii
        
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(ii, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
        zgrid, pz = eazy.getEazyPz(ii, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
        
        fig = unicorn.catalogs.plot_init(xs=8, aspect=1./4*0.95, left=0.08, right=0.02, bottom=0.05, top=0.04, fontsize=8, NO_GUI=True)
        
        #### SED
        ax = fig.add_subplot(141)
        ymax = (fobs)[(efobs > 0) & (fobs > -99)].max()
        ax.errorbar(lci, fobs, efobs, marker='o', linestyle='None', color='black', ecolor='black', ms=4, alpha=0.7)
        ax.plot(lambdaz, temp_sed, color='blue', alpha=0.3)
        ax.set_ylim(-0.05*ymax, 1.1*ymax)
        ax.set_xlim(3000,8.2e4)
        ax.semilogx()
        
        #### p(z)
        ax = fig.add_subplot(142)
        ax.plot(zgrid, pz, color='orange')
        ax.fill(zgrid, pz, color='orange', alpha=0.5)
        if zout.z_spec[ii] > 0:
            ax.plot(np.array([1,1])*zout.z_spec[ii], [0, 1.05*pz.max()], color='red', linewidth=1)
        
        #### Label
        ax.text(1.1, 1.05, '%s-%d, K=%.2f, z=%.2f, logM=%.2f, Av=%.2f' %(field, cat.id[ii], cat.kmag[ii], zout.z_peak[ii], fout.lmass[ii], fout.Av[ii]), transform=ax.transAxes, horizontalalignment='center')
        
        #### F125W thumbnail
        x125, y125 = wcs_f125w.rd2xy((cat.ra[ii], cat.dec[ii]), hour=False)
        x125, y125 = np.int(np.round(x125))-1, np.int(np.round(y125))-1
        
        x160, y160 = wcs_f160w.rd2xy((cat.ra[ii], cat.dec[ii]), hour=False)
        x160, y160 = np.int(np.round(x160))-1, np.int(np.round(y160))-1
        
        thumb125 = im_f125w[0].data[y125-N_125:y125+N_125, x125-N_125:x125+N_125]
        thumb160 = im_f160w[0].data[y160-N_160:y160+N_160, x160-N_160:x160+N_160]
        
        if save_fits:
            hdu = [pyfits.PrimaryHDU()]
            header = pyfits.Header()
            header.update('RA', cat.ra[ii], comment='Right Ascension')
            header.update('DEC', cat.dec[ii], comment='Declination')
            
            header.update('FILTER','F125W')
            header.update('SCALE', np.abs(wcs_f125w.cd11)*3600., comment='Pixel scale')
            hdu.append(pyfits.ImageHDU(data=thumb125, header=header))

            header.update('FILTER','F160W')
            header.update('SCALE', np.abs(wcs_f160w.cd11)*3600., comment='Pixel scale')
            hdu.append(pyfits.ImageHDU(data=thumb160, header=header))
            
            hdul = pyfits.HDUList(hdu)
            hdul.writeto('%s-%05d_%3.1f.fits' %(field, cat.id[ii], zout.z_peak[ii]), clobber=True)
            
        l125 = np.log10(thumb125+0.1); l125[~np.isfinite(l125)]
        l160 = np.log10(thumb160+0.1); l160[~np.isfinite(l160)]
        
        scl = (-0.3,1)
        scl = (0.1,1)
        ax = fig.add_subplot(143)
        #ax.imshow(l125, vmin=-1.2, vmax=0.3, interpolation='nearest')
        ax.imshow(0-l125, vmin=scl[0], vmax=scl[1], interpolation='nearest')
        xtick = ax.set_xticks(ticks_125);ytick = ax.set_yticks(ticks_125) 
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_ylabel('F125W')

        
        ax = fig.add_subplot(144)
        #ax.imshow(l160, vmin=-1.2, vmax=0.3, interpolation='nearest')
        ax.imshow(0-l160, vmin=scl[0], vmax=scl[1], interpolation='nearest')
        xtick = ax.set_xticks(ticks_160);ytick = ax.set_yticks(ticks_160) 
        ax.set_ylabel('F160W')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        unicorn.catalogs.savefig(fig, '%s-%05d_%3.1f.pdf' %(field, cat.id[ii], zout.z_peak[ii]))
        
        plt.close()
        
def massive_galaxies_spectra(field='COSMOS', masslim=11.2, zlim=[1.7,6]):
    
    import threedhst
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import unicorn
    import shapely
    from shapely.geometry import Point,Polygon
    try:
        import stsci.tools.wcsutil as wcsutil
    except:
        import pytools.wcsutil as wcsutil
    
    cat, zout, fout = unicorn.analysis.read_catalogs(root='%s-1' %(field))
    try:
        test = cat.use
    except:
        cat.use = np.ones(len(cat.id), dtype=np.int)
        
    full_geom, xreg, yreg = unicorn.candels.get_full_candels_region(field=field)
    
    in_candels = cat.id < 0
    print 'Find objects within CANDELS'

    idx = np.arange(len(cat.id))
    for i in idx:
        #print unicorn.noNewLine+'%d' %(i)
        in_candels[i] = full_geom.contains(Point(cat.ra[i], cat.dec[i]))

    print unicorn.noNewLine+'Find objects within CANDELS - Done.'
        
    #
    eazy_path = {}
    eazy_path['AEGIS'] = '/3DHST/Ancillary/AEGIS/NMBS/Photometry//aegis-n2.deblend.redshifts'
    eazy_path['COSMOS'] = '/3DHST/Ancillary/COSMOS/NMBS/Photometry//cosmos-1.deblend.redshifts'
    eazy_path['GOODS-S'] = '/3DHST/Ancillary/GOODS-S/FIREWORKS/FAST'
    eazy_path['GOODS-N'] = '/3DHST/Ancillary/GOODS-N/MODS/Photometry/FAST/../EAZY/OUTPUT'
    eazy_path['UDS'] = '/3DHST/Ancillary/UDS/UKIDSS/FAST'
    
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    OUTPUT_DIRECTORY = eazy_path[field]
    CACHE_FILE='Same'
    
    F160W = {}
    F160W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F160W_drz_sci.fits'
    F160W['AEGIS'] = '/Users/gbrammer/CANDELS/EGS/PREP_FLT/EGS-F160W_drz_sci.fits'
    F160W['UDS'] = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'
    F160W['GOODS-S'] = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    
    im_f160w = pyfits.open(F160W[field])
    
    wcs_f160w = wcsutil.WCSObject('f160w',header=im_f160w[0].header)
        
    thumb_size = 8 # arcsec
    ps_f160w = np.sqrt(wcs_f160w.cd11**2+wcs_f160w.cd12**2)*3600.
    
    N_160 = int(np.round(thumb_size / ps_f160w / 2))
    ticks_160 = np.arange(0,thumb_size+1)/ps_f160w
    
    #massive = in_candels & (cat.use == 1) & (fout.lmass > 10.98) & (zout.z_peak > 1.5)
    
    massive = in_candels & (cat.use == 1) & (fout.lmass > masslim) & (zout.z_peak > zlim[0]) & (zout.z_peak < zlim[1])
    
    #### 3D-HST interlaced objects
    specs = glob.glob(unicorn.GRISM_HOME+'/%s/PREP_FLT/*wcsfix' %(field))
    specs_id, specs_ra, specs_dec, specs_pointing = [], [], [], []
    for spec in specs:
        print unicorn.noNewLine+spec
        wcs = catIO.Readfile(spec)
        pointing = int(spec.split('FLT/%s-' %(field))[1].split('_inter')[0])
        specs_id.extend(list(wcs.id))
        specs_ra.extend(list(wcs.ra))
        specs_dec.extend(list(wcs.dec))
        specs_pointing.extend([pointing]*len(wcs.id))
    
    specs_id, specs_ra, specs_dec, specs_pointing = np.array(specs_id), np.array(specs_ra), np.array(specs_dec), np.array(specs_pointing)
    
    for ii in idx[massive]:
        print ii
        
        dr = np.sqrt((cat.ra[ii]-specs_ra)**2/np.cos(cat.ra[ii]/360.*2*np.pi)**2+(cat.dec[ii]-specs_dec)**2)*3600.
        match = dr < 1
        if np.sum(match) == 0:
            continue
        
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(ii, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
        zgrid, pz = eazy.getEazyPz(ii, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
        
        fig = unicorn.catalogs.plot_init(xs=8, aspect=2./3, left=0.08, right=0.02, bottom=0.05, top=0.04, fontsize=8, NO_GUI=True)
                        
        #### SED
        ax = fig.add_subplot(231)
        ymax = (fobs)[(efobs > 0) & (fobs > -99)].max()
        ax.errorbar(lci, fobs, efobs, marker='o', linestyle='None', color='black', ecolor='black', ms=4, alpha=0.7)
        ax.plot(lambdaz, temp_sed, color='blue', alpha=0.3)
        ax.set_ylim(-0.05*ymax, 1.1*ymax)
        ax.set_xlim(3000,8.2e4)
        ax.semilogx()
        
        #### p(z)
        ax = fig.add_subplot(232)
        ax.plot(zgrid, pz, color='orange')
        ax.fill(zgrid, pz, color='orange', alpha=0.5)
        if zout.z_spec[ii] > 0:
            ax.plot(np.array([1,1])*zout.z_spec[ii], [0, 1.05*pz.max()], color='red', linewidth=1)
        
        #### Label
        ax.text(0.5, 1.03, '%s-%d, K=%.2f, z=%.2f, logM=%.2f, Av=%.2f' %(field, cat.id[ii], cat.kmag[ii], zout.z_peak[ii], fout.lmass[ii], fout.Av[ii]), transform=ax.transAxes, horizontalalignment='center')
        
        #### F160W thumbnail
        x160, y160 = wcs_f160w.rd2xy((cat.ra[ii], cat.dec[ii]), hour=False)
        x160, y160 = np.int(np.round(x160))-1, np.int(np.round(y160))-1
        
        thumb160 = im_f160w[0].data[y160-N_160:y160+N_160, x160-N_160:x160+N_160]
        
        l160 = np.log10(thumb160+0.1); l160[~np.isfinite(l160)]
        
        scl = (-0.3,1)
        scl = (0.1,1)
        ax = fig.add_subplot(233)
        #ax.imshow(l160, vmin=-1.2, vmax=0.3, interpolation='nearest')
        ax.imshow(0-l160, vmin=scl[0], vmax=scl[1], interpolation='nearest')
        xtick = ax.set_xticks(ticks_160);ytick = ax.set_yticks(ticks_160) 
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_ylabel('F160W')
                    
        icount = 0
        NMAT = np.maximum(np.sum(match),2)
        for id, pointing in zip(specs_id[match], specs_pointing[match]):
            try:
                spc2D = pyfits.open(unicorn.GRISM_HOME+'/%s/PREP_FLT/%s-%d_%05d.2D.fits' %(field, field, pointing, id))            
            except:
                NMAT -= 1
        #
        NMAT = np.maximum(NMAT,2)
        
        for id, pointing in zip(specs_id[match], specs_pointing[match]):
            try:
                spc2D = pyfits.open(unicorn.GRISM_HOME+'/%s/PREP_FLT/%s-%d_%05d.2D.fits' %(field, field, pointing, id))            
            except:
                continue
            
            spc1D = unicorn.reduce.Interlace1D(unicorn.GRISM_HOME+'/%s/PREP_FLT/%s-%d_%05d.1D.fits' %(field, field, pointing, id), PNG=False)
            #spc1D.show()
            #
            icount +=1 
            ax = fig.add_subplot(200+10*NMAT+NMAT+icount)
            ax.plot(lambdaz/1.e4, temp_sed, color='blue', linewidth=2, alpha=0.4)
            ax.errorbar(lci/1.e4, fobs, efobs, marker='o', linestyle='None', color='red', ecolor='red', ms=6, alpha=0.7)
            plt.xlim(0.8, 2.1)
            ax.set_title('%s-%d_%05d.1D.fits' %(field, pointing, id))
            #            
            keep = (spc1D.data.wave > 1.12e4) & (spc1D.data.wave < 1.6e4)
            if np.sum(keep) > 10:
                wave = np.cast[np.double](spc1D.data.wave)
                flux = (spc1D.data.flux-spc1D.data.contam)/spc1D.data.sensitivity
                err = spc1D.data.error/spc1D.data.sensitivity
                #
                NW = len(wave)
                NBIN = 3
                wbin = np.zeros(NW/NBIN)
                fbin, ebin = wbin*1., wbin*1.
                for i in range(NW/NBIN):
                    wbin[i] = wave[i*NBIN+1]
                    fbin[i] = np.sum(flux[i*NBIN:i*NBIN+NBIN])/3.
                    ebin[i] = np.sqrt(np.sum(err[i*NBIN:i*NBIN+NBIN]**2))/3.
                #    
                #
                yint = unicorn.utils_c.interp_conserve_c(wave, lambdaz, temp_sed) 
                anorm = np.sum(flux[keep]*yint[keep]/err[keep]**2)/np.sum(yint[keep]**2/err[keep]**2)
                #ax.errorbar(spc1D.data.wave/1.e4, flux/anorm, err/anorm, marker=',', color='black', alpha=0.5, linestyle='None')
                ax.plot(spc1D.data.wave/1.e4, flux/anorm, color='black', alpha=0.2)
                ax.plot(spc1D.data.wave/1.e4, spc1D.data.contam/spc1D.data.sensitivity/anorm, color='orange', alpha=0.5, linewidth=2)
                ax.errorbar(wbin/1.e4, fbin/anorm, ebin/anorm, marker='o', linestyle='None', ms=3, alpha=0.4, capsize=1.1, color='black', ecolor='black')
                ymax = flux[keep].max()/anorm
                #print ymax
                ax.set_ylim(-0.05*ymax, 1.1*ymax)
                ax.set_xlim(0.9,1.8)
                
        unicorn.catalogs.savefig(fig, '%s-%05d_spec.pdf' %(field, cat.id[ii]))

        plt.close()
        
def proposal_figure():
    
    ### Compact: COSMOS-01966, mass=11.28, z=2.28
    ### Merger: COSMOS-5024, mass=11.32, z=2.61
    ### Extended: COSMOS-16561
    import threedhst
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import unicorn
    import shapely
    from shapely.geometry import Point,Polygon
    try:
        import stsci.tools.wcsutil as wcsutil
    except:
        import pytools.wcsutil as wcsutil
    
    field = 'COSMOS'
    compact = 1966
    merger = 5024
    extended = 16561 # 12867
    
    cat, zout, fout = unicorn.analysis.read_catalogs(root='COSMOS')
    cat.Jmag = 25-2.5*np.log10(cat.J)
    cat.Jmag = 25-2.5*np.log10(cat.J1*0.5+cat.J2*0.5)
    cat.Hmag = 25-2.5*np.log10(cat.H)
    
    F160W = {}
    F160W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F160W_drz_sci.fits'
    F160W['AEGIS'] = '/Users/gbrammer/CANDELS/EGS/PREP_FLT/EGS-F160W_drz_sci.fits'
    F160W['UDS'] = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'
    F160W['GOODS-S'] = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    
    F125W = {}
    F125W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F125W_drz_sci.fits'
    F125W['COSMOS'] = '/Users/gbrammer/CANDELS/COSMOS/PREP_FLT/COSMOS-full-F125W_drz_sci.fits'
    F125W['AEGIS'] = '/Users/gbrammer/CANDELS/EGS/PREP_FLT/EGS-F125W_drz_sci.fits'
    F125W['UDS'] = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f125w_v1.0_drz.fits'
    F125W['GOODS-S'] = '/Users/gbrammer/CANDELS/GOODS-S/UCSC/GOODS-S_F125W_wfc3ir_drz_sci.fits'
    
    im_f125w = pyfits.open(F125W[field])
    im_f160w = pyfits.open(F160W[field])
    
    wcs_f125w = wcsutil.WCSObject('f125w',header=im_f125w[0].header)
    wcs_f160w = wcsutil.WCSObject('f160w',header=im_f160w[0].header)
        
    thumb_size = 8 # arcsec
    ps_f125w = np.sqrt(wcs_f125w.cd11**2+wcs_f125w.cd12**2)*3600.
    ps_f160w = np.sqrt(wcs_f160w.cd11**2+wcs_f160w.cd12**2)*3600.
    
    N_125 = int(np.round(thumb_size / ps_f125w / 2))
    N_160 = int(np.round(thumb_size / ps_f160w / 2))
    ticks_125 = np.arange(0,thumb_size+1)/ps_f125w
    ticks_160 = np.arange(0,thumb_size+1)/ps_f160w
        
    labels = ['Compact','Extended','Merger/clump']
    ids = [compact, extended, merger]
    scl = (-0.3,1)
    scl = (0.1,1)
    scales = [(0.1,1),(0.3,1),(0.3,1)]
    
    
    idx= np.arange(len(cat.id))
    
    vert = True
    if vert:
        aspect = 3./2
        xs=4
    else:
        aspect = 2./3
        xs=4*3/2
        
    fig = unicorn.catalogs.plot_init(xs=xs, aspect=aspect, left=0.01, right=0.01, bottom=0.01, top=0.01, fontsize=10)
    fig.subplots_adjust(wspace=0.00,hspace=0.00)
    
    for i in range(3):
               
        ii = idx[cat.id == ids[i]][0]
        
        x125, y125 = wcs_f125w.rd2xy((cat.ra[ii], cat.dec[ii]), hour=False)
        x125, y125 = np.int(np.round(x125))-1, np.int(np.round(y125))-1
        
        x160, y160 = wcs_f160w.rd2xy((cat.ra[ii], cat.dec[ii]), hour=False)
        x160, y160 = np.int(np.round(x160))-1, np.int(np.round(y160))-1
        
        thumb125 = im_f125w[0].data[y125-N_125:y125+N_125, x125-N_125:x125+N_125]*1.5
        thumb160 = im_f160w[0].data[y160-N_160:y160+N_160, x160-N_160:x160+N_160]
        
        l125 = np.log10(thumb125+0.1); l125[~np.isfinite(l125)]
        l160 = np.log10(thumb160+0.1); l160[~np.isfinite(l160)]
        
        scl = scales[i]
        if vert:
            ax = fig.add_subplot(321+i*3)
        else:
            ax = fig.add_subplot(231+i)
            
        #ax.imshow(l125, vmin=-1.2, vmax=0.3, interpolation='nearest')
        ax.imshow(0-l125, vmin=scl[0], vmax=scl[1], interpolation='nearest')
        xtick = ax.set_xticks(ticks_125);ytick = ax.set_yticks(ticks_125) 
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(0.1, 0.1, r'F125W / $J$=%.1f' %(cat.Jmag[ii]), transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom')
        ax.text(0.1, 0.9, r'C1-%05d' %(ids[i]), transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')
        ax.text(0.1, 0.8, r'$z=$%.2f, log M/M$_\odot$=%.1f' %(zout.z_peak[ii], fout.lmass[ii]), transform=ax.transAxes, horizontalalignment='left', verticalalignment='top')

        if vert:
            ax = fig.add_subplot(321+i*3)
        else:
            ax = fig.add_subplot(231+i+3)
            
        #ax.imshow(l160, vmin=-1.2, vmax=0.3, interpolation='nearest')
        ax.imshow(0-l160, vmin=scl[0], vmax=scl[1], interpolation='nearest')
        xtick = ax.set_xticks(ticks_160);ytick = ax.set_yticks(ticks_160) 
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(0.1, 0.1, r'F160W / $H$=%.1f' %(cat.Hmag[ii]), transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom')
        ax.text(0.5, 0.9, r'$%s$' %(labels[i]), transform=ax.transAxes, horizontalalignment='center', verticalalignment='top', fontsize=12)

def go_make_thumbnails():
    """
    Make thunbnails with a range of sizes for the 3D-HST fields
    """
    import unicorn.candels
    import glob
    
    box_sizes = [1.5, 3.0, 6.0] # 12
    for field in ['COSMOS','GOODS-N','GOODS-S']:
        for box in box_sizes:
            unicorn.candels.threedhst_RGB_thumbnails(field=field, box_size=box, skip=True)
    
    ### Make separate directories for each size to make easier to copy/download
    for field in ['COSMOS','GOODS-N','GOODS-S']:
        for box in box_sizes:
            print '%s/%04.1f' %(field, box)
            try:
                os.mkdir('%s/%04.1f' %(field, box))
            except:
                pass
            #
            os.system('mv %s/*%04.1f.png %s/%04.1f/' %(field, box, field, box))

            
def threedhst_RGB_thumbnails(field='COSMOS', box_size=3, skip=True, rgb_channel=0):
    """
    Make thumbnails for all grism objects in a given 3D-HST field
    
    Thumbnail `box_size` is the radius, given in arcsec
    
    rgb_channel = 0: All RGB
    rgb_chanel = (1,2,3): only (R,G,B)
    
    GBB, 10/10/2012
    """
    import pywcs
    import threedhst.catIO as catIO
    
    #### Working directory
    os.chdir('/3DHST/Spectra/Work/RGB/v2.1')
    try:
        os.mkdir(field)
    except:
        pass
        
    catalogs = {}
    catalogs['COSMOS'] = 'cosmos_3dhst.v2.1.cat'
    catalogs['GOODS-N'] = 'goodsn_3dhst.v2.1.cat'
    catalogs['GOODS-S'] = 'goodss_3dhst.v2.1.cat'
    catalogs['AEGIS'] = 'aegis_3dhst.v2.1.cat'
    catalogs['UDS'] = 'uds_3dhst.v2.1.cat'
    
    zfit = catIO.Readfile('/3DHST/Spectra/Release/v2.1/%s/WFC3/%s.zfit.linematched.dat' %(field, field), save_fits=False)
    cat = catIO.Readfile('/3DHST/Photometry/Release/v2.1/%s/Catalog/%s' %(field, catalogs[field]), save_fits=False)
    
    #### Images
    if field == 'GOODS-N':
        mag = 25-2.5*np.log10(cat.f_f140w)
        PATH = '/3DHST/Ancillary/GOODS-N/CANDELS/Gabe'
        im_r = pyfits.open(os.path.join(PATH, 'GN-v3-F160W_drz_sci.fits'))
        im_g = pyfits.open(os.path.join(PATH, 'GN-v3-F125W_drz_sci.fits'))
        im_b = pyfits.open(os.path.join(PATH, '../../GOODS_ACS/GN-ACSi_v3.fits'))
    
    if field == 'COSMOS':
        mag = 25-2.5*np.log10(cat.f_f140w)
        PATH = '/3DHST/Ancillary/COSMOS/CANDELS/Gabe/'
        im_r = pyfits.open(os.path.join(PATH, 'COSMOS-full-F160W_drz_sci.fits'))
        im_g = pyfits.open(os.path.join(PATH, 'COSMOS-full-F125W_drz_sci.fits'))
        im_b = pyfits.open(os.path.join(PATH, 'COSMOS-full-ACSi.fits'))
    
    if field == 'GOODS-S':
        mag = 25-2.5*np.log10(cat.f_f140w)
        PATH = '/3DHST/Ancillary/GOODS-S/'
        im_r = pyfits.open(os.path.join(PATH, 'CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits'))
        im_g = pyfits.open(os.path.join(PATH, 'CANDELS/ucsc_mosaics/GOODS-S_F125W_wfc3ir_drz_sci.fits'))
        im_b = pyfits.open(os.path.join(PATH, 'GOODS_ACS/GS-ACSi.fits'))
    
    if field == 'UDS':
        mag = 25-2.5*np.log10(cat.f_f140w)
        PATH = '/Volumes/robot/3DHST/Photometry/Work/UDS/v2/images'
        im_r = pyfits.open(os.path.join(PATH, 'UDS_F160W_sci.fits'))
        im_g = pyfits.open(os.path.join(PATH, 'UDS_F125W_sci.fits'))
        im_b = pyfits.open(os.path.join(PATH, 'UDS_F814W_sci.fits'))

    if field == 'AEGIS':
        mag = 25-2.5*np.log10(cat.f_f140w)
        PATH = '/3DHST/Photometry/Work/AEGIS/IMAGES/SMALL_PIXSCL/'
        im_r = pyfits.open(os.path.join(PATH, 'AEGIS_F160W_sci.fits'))
        im_g = pyfits.open(os.path.join(PATH, 'AEGIS_F125W_sci.fits'))
        im_b = pyfits.open(os.path.join(PATH, 'AEGIS_F814W_sci.fits'))
        im_14 = pyfits.open(os.path.join(PATH, 'AEGIS_F140W_sci.fits'))
        
        ### Make two-color F814W/F140W where CANDELS not available
        fill = (im_r[0].data == 0) & (im_14[0].data != 0)
        r_scale = 10**(-0.4*(26.46-25.96))
        im_r[0].data[fill] = im_14[0].data[fill]*r_scale
        b_scale = 10**(-0.4*(25.94-25.96))*1.5
        im_g[0].data[fill] = (im_r[0].data[fill]*0.5 + im_b[0].data[fill]*0.5*b_scale) / 10**(-0.4*(26.25-25.96))
    
    print 'Reading large images....'
    
    ### Image WCS
    shape = im_r[0].data.shape
    wcs = pywcs.WCS(im_r[0].header)
    
    ### Objects with grism spectrum IDs
    keep = (zfit.spec_id != '00000')
    idx = np.arange(len(keep))[keep]
    idx = idx[np.argsort(mag[idx])]
    
    ### Box size
    pix_scale = im_r[0].header['CD1_1']**2
    if 'CD1_2' in im_r[0].header.keys():
        pix_scale += im_r[0].header['CD1_2']**2
    
    pix_scale = np.sqrt(pix_scale)*3600.
        
    NX = int(np.round(box_size/pix_scale))
    NY = NX
    
    ### View in DS9
    use_ds9 = False
    
    #### Scale parameters
    Q, alpha, m0 = 5.,3.,-0.02
  
    #### for fainter galaxies, lower SB features
    Q, alpha, m0 = 3.5, 5, -0.01
    
    #### Interpolate scale parameters with mag
    Q_i = [5,3.5]
    alpha_i = [3,5]
    m0_i = [-0.02,-0.01]
    mag_i = [18,21]
    
    for i in range(len(idx)):
        obj = zfit.spec_id[idx][i]
        out_image = '%s/%s_rgb_%04.1f.png' %(field, obj, box_size)
        if os.path.exists(out_image) & skip:
            continue
        #
        ra, dec = cat.ra[idx][i], cat.dec[idx][i]
        #ra, dec = np.cast[float](ds9.get('pan fk5').split())
        xy = np.round(wcs.wcs_sky2pix(ra, dec,0))
        xc, yc = int(xy[0]), int(xy[1])
        if (xc < 0) | (yc < 0) | (xc > shape[1]) | (yc > shape[0]):
            continue
        #
        # Browse with DS9 one by one
        if use_ds9:
            xy = np.round(np.cast[float](ds9.get('pan').split()))
            xc, yc = int(xy[0]), int(xy[1])
            obj = 'tmp'
        
        ### F160W
        sub_r = im_r[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.96-25.96))
        ### F125W
        sub_g = im_g[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(26.25-25.96))
        ### F814W
        sub_b = im_b[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.94-25.96))*1.5
        
        if rgb_channel == 1:
            sub_g*=0
            sub_b*=0
            out_image = out_image.replace('.png','_R.png')
        
        if rgb_channel == 2:
            sub_r*=0
            sub_b*=0
            out_image = out_image.replace('.png','_G.png')
        
        if rgb_channel == 3:
            sub_r*=0
            sub_g*=0
            out_image = out_image.replace('.png','_B.png')

        #
        #### Interpolate scale parameters
        Q = np.interp(mag[idx][i], mag_i, Q_i, left=Q_i[0], right=Q_i[1])
        alpha = np.interp(mag[idx][i], mag_i, alpha_i, left=alpha_i[0], right=alpha_i[1])
        m0 = np.interp(mag[idx][i], mag_i, m0_i, left=m0_i[0], right=m0_i[1])
        #
        unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename=out_image, shape=sub_r.shape)
        #labels = ['%s %d' %(root, cat2.cat.id[idx][i]), 'z= %.2f' %(cat2.zout.z_peak[idx][i]), 'log M= %.1f' %(cat2.fout.lmass[idx][i])]
        #old_im = 'RGB/%s_0.png' %(obj)
        #unicorn.candels.thumbnail_annotate(old_im, label=labels, italic=True, fit=True).save(old_im.replace('.png','_t.png'))
        print unicorn.noNewLine + obj + ' (%d of %d)' %(i+1, len(idx))
    
def test_lupton():
    import pywcs
    import unicorn.catalogs2 as cat2
    
    cat2.read_catalogs('GOODS-N')
    mag = 25-2.5*np.log10(cat2.cat.f_f140)
    PATH = '/Volumes/Crucial/3DHST/Ancillary/GOODS-N/HST_MOSAICS'
    im_r = pyfits.open(os.path.join(PATH, 'GN-v2-F160W_drz_sci.fits'))
    im_g = pyfits.open(os.path.join(PATH, 'GN-v2-F125W_drz_sci.fits'))
    im_b = pyfits.open(os.path.join(PATH, 'GN-ACSi.fits'))
    root='GOODS-N'

    cat2.read_catalogs('COSMOS')
    mag = 25-2.5*np.log10(cat2.cat.F160W)
    PATH = '/research/HST/CANDELS/COSMOS/PREP_FLT'
    im_r = pyfits.open(os.path.join(PATH, 'COSMOS-full-F160W_drz_sci.fits'))
    im_g = pyfits.open(os.path.join(PATH, 'COSMOS-full-F125W_drz_sci.fits'))
    im_b = pyfits.open(os.path.join(PATH, 'COSMOS-full-ACSi.fits'))
    root='COSMOS'

    cat2.read_catalogs('GOODS-S')
    mag = 25-2.5*np.log10(cat2.cat.f_f160)
    PATH = '/Volumes/Crucial/3DHST/Ancillary/GOODS-S'
    im_r = pyfits.open(os.path.join(PATH, 'UCSC/GOODS-S_F160W_wfc3ir_drz_sci.fits'))
    im_g = pyfits.open(os.path.join(PATH, 'UCSC/GOODS-S_F125W_wfc3ir_drz_sci.fits'))
    im_b = pyfits.open(os.path.join(PATH, 'GS-ACSi.fits'))
    root='GOODS-S'
    
    cat, zout, fout = unicorn.analysis.read_catalogs('UDS-17')
    cat2.cat, cat2.zout, cat2.fout = cat, zout, fout
    cat2.zfit = zout
    cat2.zfit.z_max_spec = cat2.zfit.z_peak
    
    mag = 25-2.5*np.log10(cat2.cat.k_totf)
    PATH = '/3DHST/Ancillary/UDS/CANDELS'
    im_r = pyfits.open(os.path.join(PATH, 'hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'))
    im_g = pyfits.open(os.path.join(PATH, 'hlsp_candels_hst_wfc3_uds-tot_f125w_v1.0_drz.fits'))
    im_b = pyfits.open(os.path.join(PATH, 'UDS-f814w.fits'))
    root='UDS'
    
    shape = im_r[0].data.shape
    wcs = pywcs.WCS(im_r[0].header)
    
    #### selection
    test = ((cat2.zout.z_peak > 0.1) & (cat2.fout.lmass > 10.) & (cat2.fout.Av > -2) & (mag > 16) & (mag < 24)) | (cat2.zfit.z_max_spec >= 0) #& (cat.f160W_flux_radius > 3)
    
    #test = (cat2.zfit.z_max_spec >= 0)
    
    idx = np.arange(len(test))[test]
    idx = idx[np.argsort(mag[idx])]
    
    i = 0
    NX, NY = 83, 83 # 5" radius

    #### Run it
    i = i+1
    skip=True
    
    use_ds9 = False

    Q, alpha, m0 = 5.,3.,-0.05
    
    if root == 'UDS':
        Q = 3.5
        m0 = -0.04
        
    for i in range(len(idx)):
        obj = '%s-%05d' %(root, cat2.cat.id[idx][i])
        if (os.path.exists(os.path.join(PATH,'RGB/')+obj+'_0.png') & skip):
            continue
        #
        xy = np.round(wcs.wcs_sky2pix(cat2.cat.ra[idx][i], cat2.cat.dec[idx][i],0))
        xc, yc = int(xy[0]), int(xy[1])
        if (xc < 0) | (yc < 0) | (xc > shape[1]) | (yc > shape[0]):
            continue
        #
        # Browse with DS9 one by one
        if use_ds9:
            xy = np.round(np.cast[float](ds9.get('pan').split()))
            xc, yc = int(xy[0]), int(xy[1])
            obj = 'tmp'
        #
        sub_r = im_r[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.96-25.96))
        sub_g = im_g[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(26.25-25.96))
        sub_b = im_b[0].data[yc-NY:yc+NY, xc-NX:xc+NX]*10**(-0.4*(25.94-25.96))*1.5
        #
        unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename=os.path.join(PATH,'RGB/')+obj+'_0.png')
        #unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q/2., alpha=alpha*1.05, m0=m0, filename=os.path.join(PATH,'RGB/')+obj+'_2.png')
        #unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q/4., alpha=alpha*1.1, m0=m0, filename=os.path.join(PATH,'RGB/')+obj+'_4.png')
        labels = ['%s %d' %(root, cat2.cat.id[idx][i]), 'z= %.2f' %(cat2.zout.z_peak[idx][i]), 'log M= %.1f' %(cat2.fout.lmass[idx][i])]
        old_im = 'RGB/%s_0.png' %(obj)
        unicorn.candels.thumbnail_annotate(old_im, label=labels, italic=True, fit=True).save(old_im.replace('.png','_t.png'))
        print unicorn.noNewLine + obj + ' (%d of %d)' %(i+1, len(idx))
        
    ###  Selection for the bar fraction
    test = ((cat2.zout.z_peak > 0.1) & (cat2.fout.lmass > 10.4) & (cat2.fout.Av > -2) & (mag > 18) & (mag < 24)) & (cat2.cat.f160W_flux_radius > 3)
    
    idx = np.arange(len(test))[test]
    
    idx = idx[np.argsort((cat2.cat.aimage/cat2.cat.bimage)[idx])]
    idx = idx[np.argsort((cat2.zout.z_peak)[idx])]
    
    #### sequence of theta in cosmos
    test = ((cat2.zout.z_peak > 0.1) & (cat2.fout.lmass > 10.) & (cat2.fout.Av > -2) & (mag > 18) & (mag < 22)) & (cat2.cat.f160W_flux_radius > 3) & (cat2.cat.f160W_ellip > 0.4)
    
    idx = np.arange(len(test))[test]
    idx = idx[np.argsort((cat2.cat.f160W_theta)[idx])]
    
    os.system('rm RGB/seq*png')
    
    for i in range(len(idx)):
        obj = '%s-%05d' %(root, cat.id[idx][i])
        print unicorn.noNewLine + obj + ' (%d of %d)' %(i+1, len(idx))
        #### make sequence files
        shutil.copy('RGB/%s_0.png' %(obj), 'RGB/seq_%d.png' %(i+1))
        #### Add labels
        labels = ['%s %d' %(root, cat2.cat.id[idx][i]), 'z= %.2f' %(cat2.zout.z_peak[idx][i]), 'log M= %.1f' %(cat2.fout.lmass[idx][i])]
        old_im = 'RGB/%s_0.png' %(obj)
        unicorn.candels.thumbnail_annotate(old_im, label=labels, italic=True, fit=True).save(old_im.replace('.png','_t.png'))
        #### Copy to BARS directory
        #shutil.copy('RGB/%s_0_t.png' %(obj), '/research/HST/BARS/RGB/')
        #

    #im.rotate(45).show()
    
def thumbnail_annotate(filename, label=None, cross=0.05, font='cmunss.otf', italic=False, fit=True, fontsize=16, fontdir='STD'):
    """
    Add text and maybe a crosshair to an image thumbnail.
    
    The label text can either be a string or a list of strings that will
    each be printed on subsequent lines.
    
    The files for the LaTeX-like Computer modern font can be downloaded here:
    http://cm-unicode.sourceforge.net/download.html
    
    """
    import Image
    import ImageFont, ImageDraw, ImageOps
    
    #### Computer modern, italic font.  Check the "Fontbook" application
    #### to find your available fonts and their paths.  Appears to work 
    #### with at least 'ttf' and 'otf' formatted font files.
    if italic:
        font = 'cmunit.otf'
    
    #### Input is either a filename or a PIL image object
    if isinstance(filename, str):
        im = Image.open(filename)
    else:
        im = filename
    
    #### Select the font, use default if font file not found
    if fontdir == 'STD':
        fontpath = '%s/Library/Fonts/%s' %(os.getenv('HOME'), font)
    else:
        fontpath = '/%s' %(fontdir, font)
        
    if os.path.exists(fontpath):
        f = ImageFont.truetype(fontpath, fontsize)
    else: 
        print 'Font: %s not found, using default' %(fontpath)
        f = ImageFont.load_default()
        fit = False
        
    sh = im.size
    draw = ImageDraw.Draw(im)
    
    #### Add text labels
    if label is not None:
        #### Make sure the label will fit with the specified fontsize
        if not isinstance(label, list):
            label = [label]
        for line in label:
            fs = draw.textsize(line, font=f)
            if (fs[0] > 0.9*sh[0]) & fit:
                fontsize = int(np.round(fontsize*0.9*sh[0]/fs[0]))
                f = ImageFont.truetype('/Users/gbrammer/Library/Fonts/%s' %(font),  fontsize)
        
        #### Add the label(s)
        for il, line in enumerate(label):
            fs = draw.textsize(line, font=f)    
            draw.text((int(0.02*sh[0]), int(0.02*sh[1])+fs[1]*il), line, font=f, fill=(255,255,255))
    
    #### Add the crosshair        
    if cross > 0:
        dx, dy = 0.05*sh[0], 0.05*sh[1]
        draw.line((sh[0]/2., sh[1]/2.+3*dy, sh[0]/2., sh[1]/2.+4*dy), fill=(255,255,255))
        draw.line((sh[0]/2.-4*dx, sh[1]/2., sh[0]/2.-3*dx, sh[1]/2.), fill=(255,255,255))
    
    #im.save('junk.png')
    return im
    
def luptonRGB(imr, img, imb, Q=5, alpha=3, m0=-0.05, m1=1, shape=(300,300), filename='junk.png', ds9=None, verbose=False, rgb_clip=True):
    """
    Make a 3 color image scaled with the color clipping and 
    asinh scaling from Lupton et al. (2004)
    """   
    import Image
    
    I = (imr+img+imb-3*m0)/3.
    fI = np.arcsinh(alpha*Q*I)/Q
    M = m0 + np.sinh(Q*1.)/(alpha*Q)
    #ds9.v(fI, vmin=0, vmax=1)
    if verbose:
        print 'min, max = %f, %f' %(m0, M)
        
    fI[I < m0] = 0
    R = np.maximum(imr-m0, 0)*fI/I
    G = np.maximum(img-m0, 0)*fI/I
    B = np.maximum(imb-m0, 0)*fI/I
        
    min_RGB = np.minimum(np.minimum(R,G),B)
    zero = min_RGB < 0
    zero = fI < 0
    R[zero] = 0.
    G[zero] = 0.
    B[zero] = 0.
    
    R[R < 0] = 0
    G[G < 0] = 0
    B[B < 0] = 0
    
    max_RGB = np.maximum(np.maximum(R,G),B)
    if rgb_clip:
        clip = max_RGB > 1
        R[clip] = R[clip]/max_RGB[clip]
        G[clip] = G[clip]/max_RGB[clip]
        B[clip] = B[clip]/max_RGB[clip]
    else:
        R[R > 1] = 1.
        G[G > 1] = 1.
        B[B > 1] = 1.

    if ds9 is not None:
        #ds9.set('rgb True')
        v1=1
        ds9.set('rgb lock colorbar')
        ds9.set('rgb red'); ds9.v(R, vmin=0, vmax=v1); ds9.set('scale linear')
        ds9.set('rgb green'); ds9.v(G, vmin=0, vmax=v1); ds9.set('scale linear')
        ds9.set('rgb blue'); ds9.v(B, vmin=0, vmax=v1); ds9.set('scale linear')
        return True
        
    #rgb = np.array([R,G,B]).T
    #im = Image.merge('RGB', (Image.fromarray((R[::-1,:]*255).astype('uint8')), Image.fromarray((G[::-1,:]*255).astype('uint8')), Image.fromarray((B[::-1,:]*255).astype('uint8'))))
    im = Image.merge('RGB', (Image.fromarray((R[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((G[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((B[::-1,:]*255).astype('int8'), mode='L')))
    
    im = im.resize(shape)
    im.save(filename)
    
def logRGB(imr, img, imb, lexp=0.8e4, cmap=[1.39165, 0.628952], zs=[-0.01, 2], scale=[1.4,1.,0.3], shape=None, filename=None):
    """
    Default colors for pixel-matched F140W / F814W / F435W
    """
    import Image
    
    # imgr = pyfits.open('../hlsp_clash_hst_wfc3ir_macs0717_f140w_v1_drz.fits')
    # imgg = pyfits.open('../hlsp_clash_hst_acs_macs0717_f814w_v1_drz.fits')
    # imgb = pyfits.open('../hlsp_clash_hst_acs_macs0717_f435w_v1_drz.fits')
    # 
    # imr = imgr[0].data[2000:2900, 2000:2900]
    # img = imgg[0].data[2000:2900, 2000:2900]
    # imb = imgb[0].data[2000:2900, 2000:2900]

    # zs = [-0.01,2]
    # lexp = 0.8e4
    # cmap=[1.39165, 0.628952]
    # scale = [1,1,1]
    # 
    # scale = [1.4,1.,0.3]
    
    contrast, bias = cmap
    
    clipped = []
    ims = [imr, img, imb]
    for i in range(3):
        scl = np.array(zs)*scale[i]
        clip = (np.clip(ims[i], scl[0], scl[1])-scl[0])/(scl[1]-scl[0])
        clip_log = np.clip((np.log10(lexp*clip+1)/np.log10(lexp)-bias)*contrast+0.5, 0, 1)
        clipped.append(clip_log)
    #
    im = Image.merge('RGB', (Image.fromarray((clipped[0][::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((clipped[1][::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((clipped[2][::-1,:]*255).astype('int8'), mode='L')))
    
    if shape is None:
        shape = imr.shape
    #
    im = im.resize(shape)
    
    if filename is None:
        filename='test.png'
    
    im.save(filename)
    
    return im
    
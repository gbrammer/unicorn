import os
import glob
import shutil

import numpy as np
import pyfits
import scipy.linalg
import matplotlib.pyplot as plt

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe 

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

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
    
    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
        
    files=glob.glob('WFC3-*asn.fits')
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=True, DIRECT_HIGHER_ORDER=2)
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
    
def cosmos():
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
            
def make_asn_files(force=False):
    """
    Read a files.info file and make ASN files for each visit/filter.
    """
    list = catIO.Readfile('files.info')
    asn = threedhst.utils.ASNFile(glob.glob('../RAW/i*asn.fits')[1])
    
    visits = np.unique(list.targname)
    for visit in visits:        
        filters = np.unique(list.filter[list.targname == visit])
        print visit, filters
        for filter in filters:
            if (os.path.exists(visit+'-'+filter+'_asn.fits')) & (not force):
                continue
            
            #print asn.in_fits[1].columns
            
            use = (list.targname == visit) & (list.filter == filter)
            asn.exposures = []
            asn.product=visit+'-'+filter
            for fits in list.file[use]:
                asn.exposures.append(os.path.basename(fits).split('_flt')[0])
            
            #try:
            asn.write(asn.product+'_asn.fits')
            #except:
            #    continue
                
            threedhst.regions.asn_region(asn.product+'_asn.fits', path_to_flt='../RAW')
            print asn.product
            
    
def prep_candels(asn_file='ib3706050_asn.fits',
                       ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
                       ALIGN_EXTENSION=0,
                       PATH_TO_RAW='../RAW',
                       GET_SHIFT=True,
                       DIRECT_HIGHER_ORDER=2,
                       SCALE=0.06,
                       bg_skip=False,
                       geometry='rxyscale,shift'):
    
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
                    align_geometry=geometry, clean=True,
                    initial_order=0, save_fit=False, TWEAKSHIFTS_ONLY=(ALIGN_IMAGE is None))
    
    if DIRECT_HIGHER_ORDER > 0:
        threedhst.prep_flt_files.prep_flt(asn_file=asn_file,
                    get_shift=False, first_run=False, 
                    bg_only=False, bg_skip=bg_skip, redo_background=False,
                    skip_drz=False, final_scale=SCALE, pixfrac=0.8,
                    IMAGES=[], clean=True,
                    initial_order=DIRECT_HIGHER_ORDER, save_fit=False)

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
    
    
        
def massive_galaxies_morphologies(field='COSMOS'):
    import threedhst
    import threedhst.eazyPy as eazy
    import unicorn
    import shapely
    from shapely.geometry import Point,Polygon
    import stsci.tools.wcsutil as wcsutil
    
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
        
    thumb_size = 8 # arcsec
    ps_f125w = np.sqrt(wcs_f125w.cd11**2+wcs_f125w.cd12**2)*3600.
    ps_f160w = np.sqrt(wcs_f160w.cd11**2+wcs_f160w.cd12**2)*3600.
    
    N_125 = int(np.round(thumb_size / ps_f125w / 2))
    N_160 = int(np.round(thumb_size / ps_f160w / 2))
    ticks_125 = np.arange(0,thumb_size+1)/ps_f125w
    ticks_160 = np.arange(0,thumb_size+1)/ps_f160w
    
    plt.gray()
    
    massive = in_candels & (cat.use == 1) & (fout.lmass > 11.5) & (zout.z_peak > 2.3) & (zout.z_peak < 2.7)
    
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
        
        unicorn.catalogs.savefig(fig, '%s-%05d_%3.1f.png' %(field, cat.id[ii], zout.z_peak[ii]))
        
        plt.close()
        
def massive_galaxies_spectra(field='COSMOS'):
    
    import threedhst
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import unicorn
    import shapely
    from shapely.geometry import Point,Polygon
    import stsci.tools.wcsutil as wcsutil
    
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
    
    massive = in_candels & (cat.use == 1) & (fout.lmass > 11.5) & (zout.z_peak > 1.7)
    
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
            spc1D.show()
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
                ax.plot(spc1D.data.wave/1.e4, spc1D.data.contam/anorm, color='orange', alpha=0.5, linewidth=2)
                ax.errorbar(wbin/1.e4, fbin/anorm, ebin/anorm, marker='o', linestyle='None', ms=3, alpha=0.4, capsize=1.1, color='black', ecolor='black')
                ymax = flux[keep].max()/anorm
                #print ymax
                ax.set_ylim(-0.05*ymax, 1.1*ymax)
                ax.set_xlim(0.9,1.8)
                
        unicorn.catalogs.savefig(fig, '%s-%05d_spec.png' %(field, cat.id[ii]))

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
    import stsci.tools.wcsutil as wcsutil
    
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
        
        
    

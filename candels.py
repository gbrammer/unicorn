import os
import glob
import shutil

import numpy as np
import pyfits
import scipy.linalg
# from scipy import polyfit, polyval
import matplotlib.pyplot as plt

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe 

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

import threedhst
import threedhst.prep_flt_files
import threedhst.catIO as catIO
import unicorn

noNewLine = '\x1b[1A\x1b[1M'

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
    #
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
    noNewLine = '\x1b[1A\x1b[1M'
    
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
        print noNewLine+pointing
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
        
    
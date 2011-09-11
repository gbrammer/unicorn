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
import glob
import shutil

def check_COSMOS_stars():
    """
    Run SExtractor and try to get the FWHM of the new reduction compared to the old,
    which had many stars with central pixels rejected by the CR rejection.
    """    
    import matplotlib.pyplot as plt
    import unicorn.go_3dhst as go
    
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=26)
    os.chdir('/3DHST/Spectra/Work/COSMOS/PREP_FLT_UNICORN')
    files=glob.glob('COSMOS-*-F140W_drz.fits')
        
    for file in files:
        ROOT_GRISM = file.split('_drz.fits')[0]
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
        se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = file+'[1]'
        se.options['FILTER']    = 'Y'
        se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
        se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
        se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
        status = se.sextractImage(file+'[0]')
            
    os.system('grep "#" COSMOS-1-F140W_drz.cat > cosmos_old.cat')
    os.system('cat COSMOS-*-F140W_drz.cat |grep -v "#" >> cosmos_old.cat')
    
    #### Full mosaic
    file='COSMOS-F140W_drz.fits'
    ROOT_GRISM = file.split('_drz.fits')[0]
    se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
    se.options['WEIGHT_IMAGE']    = file+'[1]'
    status = se.sextractImage(file+'[0]')
    
    
    os.chdir('../NEW')
    files=glob.glob('../PREP_FLT/COSMOS-2[678]-F140W_drz.fits')
    for file in files:
        ROOT_GRISM = os.path.basename(file).split('_drz.fits')[0]
        iraf.imcopy(file+'[1]',ROOT_GRISM+'_SCI.fits')
        iraf.imcopy(file+'[2]',ROOT_GRISM+'_WHT.fits')
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
        se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = ROOT_GRISM+'_WHT.fits'
        se.options['FILTER']    = 'Y'
        se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
        se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
        se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
        status = se.sextractImage(ROOT_GRISM+'_SCI.fits')
        os.system('rm '+ROOT_GRISM+'_[SW]?[IT].fits')
    
    os.system('grep "#" COSMOS-1-F140W_drz.cat > cosmos_new.cat')
    os.system('cat COSMOS-*-F140W_drz.cat |grep -v "#" >> cosmos_new.cat')
    
    old = threedhst.sex.mySexCat('../PREP_FLT_UNICORN/cosmos_old.cat')
    new = threedhst.sex.mySexCat('cosmos_new.cat')
    
    ######### GOODS-N
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=26)
    os.chdir('/3DHST/Spectra/Work/GOODS-N/PREP_FLT/')
    files=glob.glob('GOODS-N-*-F140W_drz.fits')
        
    for file in files:
        ROOT_GRISM = file.split('_drz.fits')[0]
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
        se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = file+'[1]'
        se.options['FILTER']    = 'Y'
        se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
        se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
        se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
        status = se.sextractImage(file+'[0]')
            
    os.system('grep "#" GOODS-N-11-F140W_drz.cat > goods-n.cat')
    os.system('cat GOODS-N-*-F140W_drz.cat |grep -v "#" >> goods-n.cat')
    
    #### Full mosaic
    file='GOODS-N-F140W_drz.fits'
    ROOT_GRISM = file.split('_drz.fits')[0]
    se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
    se.options['WEIGHT_IMAGE']    = file+'[1]'
    status = se.sextractImage(file+'[0]')
    
    ##### GOODS-S
    os.chdir('/3DHST/Spectra/Work/GOODS-S/PREP_FLT')
    
    files=glob.glob('GOODS-S-*-F140W_drz.fits')
    for file in files:
        ROOT_GRISM = os.path.basename(file).split('_drz.fits')[0]
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
        se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = file+'[1]'
        se.options['FILTER']    = 'Y'
        se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
        se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
        se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
        status = se.sextractImage(file+'[0]')
    
    os.system('grep "#" GOODS-S-6-F140W_drz.cat > goods-s.cat')
    os.system('cat GOODS-S-*-F140W_drz.cat |grep -v "#" >> goods-s.cat')
    
    goodss = threedhst.sex.mySexCat('GOODS-N-F140W_drz.cat')
    
    
    
    ###### Make plots ################
    
    os.chdir('/3DHST/Spectra/Work/COSMOS/NEW/')
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=1.5)
    
    ax = fig.add_subplot(211)
    plt.plot(old.MAG_AUTO, old.FLUX_RADIUS, marker='o', linestyle='None', alpha=0.1, color='blue', markersize=3)
    plt.text(15,10,'GRISM_v1.6')
    plt.semilogy()
    plt.xlim(14,26)
    plt.ylim(1,20)

    ax = fig.add_subplot(212)
    plt.plot(new.MAG_AUTO, new.FLUX_RADIUS, marker='o', linestyle='None', alpha=0.1, color='blue', markersize=3)
    plt.plot(goodsn.MAG_AUTO, goodsn.FLUX_RADIUS, marker='o', linestyle='None', alpha=0.1, color='purple', markersize=3)
    plt.text(15,10,'August 16')
    plt.semilogy()
    plt.xlim(14,26)
    plt.ylim(1,20)
    plt.xlabel('MAG_AUTO')
    plt.ylabel('FLUX_RADIUS')
    
    fig.savefig('f140w_sizes_aug16.png')
    
    #### Signal to noise
    SN = np.cast[float](goodsn.FLUX_AUTO)/np.cast[float](goodsn.FLUXERR_AUTO)
    plt.plot(goodsn.MAG_AUTO, SN, marker='o', linestyle='None', alpha=0.1, color='purple', markersize=3)
    plt.semilogy()
    plt.xlim(14,30)
    plt.ylim(0.1,500)
    
    ##### Stack images and weights

    os.chdir('../PREP_FLT_UNICORN')

    NX = 100
    img_old = np.zeros((2*NX, 2*NX))
    wht_old = np.zeros((2*NX, 2*NX))
    
    files=glob.glob('COSMOS-*-F140W_drz.cat')
    for file in files:
        print file
        cat = threedhst.sex.mySexCat(file)
        drz = pyfits.open(file.replace('cat','fits'))
        radius = np.cast[float](cat.FLUX_RADIUS)
        mag = np.cast[float](cat.MAG_AUTO)
        xpix, ypix = np.round(np.cast[float](cat.X_IMAGE)), np.round(np.cast[float](cat.Y_IMAGE))
        keep = (radius < 3.5) & (mag >= 18) & (mag < 18.5)
        NOBJ = len(cat.NUMBER)
        idx = np.arange(NOBJ)[keep]
        if len(idx) > 0:
            for i in idx:
                try:
                    img_old += drz[1].data[ypix[i]-NX:ypix[i]+NX, xpix[i]-NX:xpix[i]+NX]
                    wht_old += drz[2].data[ypix[i]-NX:ypix[i]+NX, xpix[i]-NX:xpix[i]+NX]
                except:
                    pass
                    
    img_old /= img_old.max()
    wht_old /= np.median(wht_old)
    
    #
    os.chdir('../NEW')
    NX = 100
    img_new = np.zeros((2*NX, 2*NX))
    wht_new = np.zeros((2*NX, 2*NX))
    
    files=glob.glob('COSMOS-*-F140W_drz.cat')
    for file in files:
        print file
        cat = threedhst.sex.mySexCat(file)
        drz = pyfits.open('../PREP_FLT/'+file.replace('cat','fits'))
        radius = np.cast[float](cat.FLUX_RADIUS)
        mag = np.cast[float](cat.MAG_AUTO)
        xpix, ypix = np.round(np.cast[float](cat.X_IMAGE)), np.round(np.cast[float](cat.Y_IMAGE))
        keep = (radius < 3.5) & (mag >= 18) & (mag < 18.5)
        NOBJ = len(cat.NUMBER)
        idx = np.arange(NOBJ)[keep]
        if len(idx) > 0:
            for i in idx:
                try:
                    img_new += drz[1].data[ypix[i]-NX:ypix[i]+NX, xpix[i]-NX:xpix[i]+NX]
                    wht_new += drz[2].data[ypix[i]-NX:ypix[i]+NX, xpix[i]-NX:xpix[i]+NX]
                except:
                    pass
                    
    img_new /= img_new.max()
    wht_new /= np.median(wht_new)
    
    import threedhst.dq
    ds9 = threedhst.dq.myDS9()
    
    ds9.frame(1)
    ds9.v(img_old, vmin=0, vmax=1)
    ds9.set('scale log')
    ds9.set('cmap value 1.08889 0.270588')

    ds9.frame(2)
    ds9.v(img_new, vmin=0, vmax=1)
    ds9.set('scale log')
    ds9.set('cmap value 1.08889 0.270588')

    ds9.frame(3)
    ds9.v(wht_old, vmin=0.5, vmax=1.1)
    ds9.set('scale log')
    ds9.set('cmap value 6.7325 0.960993')

    ds9.frame(4)
    ds9.v(wht_new, vmin=0.5, vmax=1.1)
    ds9.set('scale log')
    ds9.set('cmap value 6.7325 0.960993')
    
    
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
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/PREP_FLT')
    
    #### Use ACS alignment images
    ALIGN = '/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_nz*drz*fits'
    
    #### Process direct + grism pairs
    direct=glob.glob('GOODS-N-*-F140W_asn.fits')
    grism = glob.glob('GOODS-N-*-G141_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rxyscale')
    
    #### Refine one set of shifts
    threedhst.shifts.refine_shifts(ROOT_DIRECT='GOODS-N-31-F140W', 
              ALIGN_IMAGE='GOODS-N-21-F140W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-31-F140W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    threedhst.gmap.makeImageMap(['GOODS-N-21-F140W_drz.fits', 'GOODS-N-21-F140W_align.fits[0]*4','GOODS-N-31-F140W_drz.fits'], aper_list=[16], polyregions=glob.glob("GOODS-N-*-F140W_asn.pointing.reg"))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('GOODS-N-[0-9]*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.startMultidrizzle(file, 
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        #
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='GOODS-N-*[0-9]-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)
    
def GOODSN_mosaic():
    import threedhst.prep_flt_files
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/PREP_FLT')
    
    #### Direct+grism mosaics
    direct_files = glob.glob('GOODS-N-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-N-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    direct_files = glob.glob('GOODS-N-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-N-G141',
                       path_to_FLT='./', run_multidrizzle=False)
        
    SCALE = 0.06
    PIXFRAC=0.8
    NX, NY = int(6840*0.128254/SCALE), int(8042*0.128254/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=189.17736, dec=62.23892,
             final_outnx = NX, final_outny=NY, ivar_weights=False)
             
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False)

    #### Make full alignment mosaic
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN), matchImage='GOODS-N-F140W_drz.fits', match_extension=1, output='GOODS-N-F850LP.fits')
    
    threedhst.gmap.makeImageMap(['GOODS-N-F140W_drz.fits','GOODS-N-F850LP.fits[0]*4'], aper_list=[15], polyregions=glob.glob('GOODS-N-*-F140W_asn.pointing.reg')) #, zmin=-0.1, zmax=3)
    
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['GOODS-N-F140W_drz.fits', 'GOODS-N-G141_drz.fits*2', 'GOODS-N-F850LP.fits[0]*6','/3DHST/Ancillary/GOODS-N/CDFN/paper13-cdfn-figure3-FB-binned1pix-smooth.fits[0]'], aper_list=zooms, tileroot=['F140W', 'G141', 'GOODS-z850','0.5-8keV'], polyregions=glob.glob('GOODS-N*F140W_asn.pointing.reg')) #, zmin=-0.1, zmax=3)
    
    iraf.imcopy('GOODS-N-F140W_drz.fits[1]', '../MOSAIC/GOODS-N-F140w_11-09-08_sci.fits')
    iraf.imcopy('GOODS-N-F140W_drz.fits[2]', '../MOSAIC/GOODS-N-F140w_11-09-08_wht.fits')
    # !tar czvf GOODS-N-F140w_11-09-08.tar.gz GOODS-N-*-F140W_shifts.txt GOODS-N-*-F140W_tweak.fits GOODS-N-*-F140W_asn.fits GOODS-N-F140W_shifts.txt GOODS-N-F140W_asn.fits
    # !mv GOODS-N-F140w_11-09-08* ../MOSAIC
    
  
def COSMOS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
        
    os.chdir(unicorn.GRISM_HOME+'COSMOS/PREP_FLT')
    #ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    ALIGN = '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits'
    #ALIGN = '../NMBS/COSMOS-1.v4.K_nosky.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rxyscale')
    
    #### Fix shifts for COSMOS-18
    threedhst.shifts.refine_shifts(ROOT_DIRECT='COSMOS-18-F140W', 
              ALIGN_IMAGE='COSMOS-17-F140W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='rxyscale', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-18-F140W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    threedhst.gmap.makeImageMap(['COSMOS-12-F140W_drz.fits', 'COSMOS-12-G141_drz.fits','/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'][0:2], aper_list=[14,15,16], polyregions=glob.glob("COSMOS-*-F140W_asn.pointing.reg"))

    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits'], aper_list=[16], polyregions=glob.glob("COSMOS-*-F140W_asn.pointing.reg"))
            
    threedhst.gmap.makeImageMap(['COSMOS-28-F140W_drz.fits','/tmp/junk.fits', 'COSMOS-28-G141_drz.fits','/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'][0:2], aper_list=[14], polyregions=glob.glob("COSMOS-*-F140W_asn.pointing.reg"))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('COSMOS-*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.startMultidrizzle(file, 
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        #
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='COSMOS-*-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)
    
def COSMOS_mosaic():
    import threedhst.prep_flt_files
    os.chdir(unicorn.GRISM_HOME+'COSMOS/PREP_FLT')

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
    PIXFRAC=0.8
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
    iraf.imcopy('COSMOS-F140W_drz.fits[1]','../MOSAIC/COSMOS-F140w_11-09-08_sci.fits')
    iraf.imcopy('COSMOS-F140W_drz.fits[2]','../MOSAIC/COSMOS-F140w_11-09-08_wht.fits')
    # !tar czvf COSMOS-F140w_11-09-08.tar.gz COSMOS-*-F140W_shifts.txt COSMOS-*-F140W_tweak.fits COSMOS-*-F140W_asn.fits COSMOS-F140W_shifts.txt COSMOS-F140W_asn.fits
    # !mv COSMOS-F140w_11-09-08* ../MOSAIC
        
    zooms = [13,14,15]
    threedhst.gmap.makeImageMap(['../MOSAIC/COSMOS-F140w_11-07-31_sci.fits[0]', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/Chandra/CC0570_img.fits[0]*1.5', '/3DHST/Ancillary/COSMOS/Spitzer/mips_24_GO3_sci_10.fits[0]*400', '/3DHST/Ancillary/COSMOS/VLA/vla_20cm_dp_sin_10.fits[0]*40000'], aper_list=zooms, tileroot=['F140W','ACS','WIRDS-K','0.5-7keV', 'MIPS-24', 'VLA-20cm'], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))
    
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
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rxyscale')
            
    # test
    threedhst.gmap.makeImageMap(['GOODS-S-24-F140W_drz.fits','GOODS-S-23-F140W_drz.fits', 'GOODS-S-24-F140W_align.fits[0]*4'], aper_list=[16], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('GOODS-S-*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.startMultidrizzle(file, 
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        #
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='GOODS-S-*-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)

def GOODSS_mosaic():
    import threedhst.prep_flt_files
    os.chdir(unicorn.GRISM_HOME+'GOODS-S/PREP_FLT')
    
    #### Direct mosaic
    direct_files = glob.glob('GOODS-S-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Direct mosaic
    direct_files = glob.glob('GOODS-S-*-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.1
    #SCALE = 0.5
    PIXFRAC=0.8
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
    
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['GOODS-S-F140W_drz.fits', 'GOODS-S-G141_drz.fits','/3DHST/Ancillary//GOODS-S/CANDELS/hlsp_candels_hst_wfc3_gsd01_f125w_v0.5_drz.fits[0]','/3DHST/Ancillary/GOODS-S/CDFS/CDFS-4Ms-0p5to8-asca-im-bin1.fits[0]'], aper_list=zooms, tileroot=['F140W','G141','F125W','0.5-8keV'], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'))
    
def UDF():
    import threedhst
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/PREP_FLT')

    #### Take GOODS-S-34 as the reference and align the other pointings
    files = glob.glob('GOODS-S-3[6-8]-F140W*asn.fits')
    for file in files:
        threedhst.shifts.refine_shifts(ROOT_DIRECT=file.split('_asn')[0], 
              ALIGN_IMAGE='GOODS-S-34-F140W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='rxyscale', clean=True)
        #
        threedhst.prep_flt_files.startMultidrizzle(file,
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
        #
        threedhst.shifts.make_grism_shiftfile(file, file.replace('F140W','G141'))
        #
        threedhst.prep_flt_files.startMultidrizzle(file.replace('F140W','G141'),
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
        
    

    os.chdir(unicorn.GRISM_HOME+'UDF/PREP_FLT')
    
    ##### Make a dumy combined image of 3D-HST + PRIMO to get the total area covered 
    ##### to make a common detection image
    # files = glob.glob('*G141_asn.fits')
    # threedhst.utils.combine_asn_shifts(files, out_root='dummy-G141',
    #                     path_to_FLT='./', run_multidrizzle=False)
    # #
    # threedhst.prep_flt_files.startMultidrizzle('dummy-G141_asn.fits',
    #          use_shiftfile=True, skysub=False,
    #          final_scale=0.2, pixfrac=0.8, driz_cr=False,
    #          updatewcs=False, clean=True, median=False)
    
    ra0, dec0 = 53.158733, -27.785316
    NX, NY = 1083*0.2/0.06, 1106*0.2/0.06
    
    ######## Combine the 3D-HST pointings
    #### Direct mosaic
    direct_files = glob.glob('GOODS-S-3[4678]-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='UDF-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### GRISM mosaic
    direct_files = glob.glob('GOODS-S-3[4678]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='UDF-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.06
    PIXFRAC=0.8
    
    threedhst.prep_flt_files.startMultidrizzle('UDF-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra0, dec=dec0, final_outnx=NX*0.06/SCALE, final_outny=NY*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('UDF-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra0, dec=dec0, final_outnx=NX*0.06/SCALE, final_outny=NY*0.06/SCALE)
             
    threedhst.gmap.makeImageMap(['GOODS-S-34-F140W_drz.fits', 'UDF-F140W_drz.fits',  'GOODS-S-34-G141_drz.fits', 'UDF-G141_drz.fits'], aper_list=[15,16], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'), zmin=-0.01, zmax=0.1)
     
    #### Make a second UDF G141 file for testing with the flux cube
    direct_files = glob.glob('GOODS-S-3[4678]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='UDF-FC-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #
    ###################
    ####  Make detection image from CANDELS
    ###################
    os.chdir('/Users/gbrammer/CANDELS/GOODS-S/PREP_FLT')
    unicorn.candels.make_asn_files()
    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    ### Find visits that overlap with the grism pointing
    fp = open('/3DHST/Spectra/Work/UDF/PREP_FLT/grism.reg')
    lines = fp.readlines()
    fp.close()
    
    region = lines[4]
    spl = np.float_(np.array(region[region.find('(')+1:region.find(')')].split(',')))
    grismx, grismy = spl[0::2],spl[1::2]
    
    region_files = glob.glob('GOODS*F160W*pointing.reg')
    visits = []
    for region_file in region_files:
        fp = open(region_file)
        lines = fp.readlines()
        fp.close()
        #
        region = lines[1]
        spl = np.float_(np.array(region[region.find('(')+1:region.find(')')].split(',')))
        regx, regy = spl[0::2],spl[1::2]
        #
        if threedhst.regions.polygons_intersect(grismx, grismy, regx, regy):
            visits.append(region_file.split('-F160W')[0])
          
    import tarfile
    fptar = tarfile.open('udf_regions.tar.gz','w|gz')
    for visit in visits:
        fptar.add(visit+'-F160W_asn.pointing.reg')
    
    fptar.close()
      
    #visits = ['GOODS-SD5-VMX','GOODS-SD5-VGX','GOODS-S205-VHS','GOODS-SDW-VH7','GOODS-S075-VCM','GOODS-SD2-V71','GOODS-SD5-VGU','GOODS-S205-VHT','GOODS-S075-VJ3','GOODS-SD5-VGV','GOODS-SD2-V7J','GOODS-S1-V3J']
    
    FORCE=True
    for visit in visits:
        if os.path.exists(visit+'-F125W_asn.fits') & (not os.path.exists(visit+'-F125W_drz.fits')) | FORCE:
            print visit
            for filter in ['F125W','F160W']:
                try:
                    unicorn.candels.prep_candels(asn_file=visit+'-'+filter+'_asn.fits', 
                        ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                        SCALE=0.06, geometry='rxyscale')
                except:
                    pass
    #
    ##### bad shifts, redo with CANDELS as alignment
    bad_shifts = ['GOODS-S075-VJ9','GOODS-S205-V60','GOODS-S205-VHY','GOODS-SD2-V7G','GOODS-SD3-VEL','GOODS-SD5-VGQ','GOODS-SDW-V02','GOODS-SDW-VHB']
    bad_shifts = ['GOODS-SD2-V7G','GOODS-SD5-VGQ']
    ALIGN = '/3DHST/Ancillary//GOODS-S/CANDELS/hlsp_candels_hst_wfc3_gsd01_f125w_v0.5_drz.fits'
    for bad in bad_shifts[0:1]:
        for filter in ['F125W','F160W']:
            try:
                unicorn.candels.prep_candels(asn_file=bad+'-'+filter+'_asn.fits', 
                    ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=0,
                    GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                    SCALE=0.06, geometry='rxyscale')
            except:
                pass
    
    ##### Something very wrong with VGQ, V7G, delete the first visit from the 
    ##### files.info file, the exposures were repeated with a HOPR.
    os.system('rm GOODS-SD5-VGQ*')
    os.system('rm GOODS-SD2-V7G*')
    
    threedhst.gmap.makeImageMap( ['GOODS-SD2-V7G-F125W_drz.fits','GOODS-SD2-V7G-F125W_align.fits[0]*4','GOODS-SD2-V7G-F160W_drz.fits','GOODS-SD2-V7G-F160W_align.fits[0]*4'], zmin=-0.06, zmax=0.6, aper_list=[15])
    
    ### Combine CANDELS visits into a single image                
    for filter in ['F125W','F160W']:
        drz_files = glob.glob('GOODS-*-'+filter+'_drz.fits')
        direct_files = []
        for drz in drz_files:
            direct_files.append(drz.replace('drz','asn'))
        #
        threedhst.utils.combine_asn_shifts(direct_files, out_root='UDF-'+filter,
            path_to_FLT='./', run_multidrizzle=False)
        #
        threedhst.prep_flt_files.startMultidrizzle('UDF-'+ filter+'_asn.fits',
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.6, driz_cr=False,
                updatewcs=False, clean=True, median=False,
                ra=ra0, dec=dec0,
                final_outnx=NX, final_outny=NY)
    #
    os.system('cp UDF-F???W_drz.fits '+unicorn.GRISM_HOME+'UDF/PREP_FLT')
    
    os.chdir(unicorn.GRISM_HOME+'UDF/PREP_FLT')
    threedhst.gmap.makeImageMap( ['UDF-F140W_drz.fits','UDF-F125W_drz.fits','UDF-F160W_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[15,16])
    
    #### Fill in empty areas of F140W image with the average of F125W and F160W
    f140 = pyfits.open('UDF-F140W_drz.fits')
    f125 = pyfits.open('UDF-F125W_drz.fits')
    f160 = pyfits.open('UDF-F160W_drz.fits')
    
    avg = 0.5*f125[1].data*10**(-0.4*(26.25-26.46))+0.5*f160[1].data*10**(-0.4*(25.96-26.46))
    wht = 0.5*f125[2].data+0.5*f160[2].data

    mask = f140[2].data < 900
    f140[1].data[mask] = avg[mask]
    f140[2].data[mask] = wht[mask]
    f140.writeto('UDF-fill-F140W_drz.fits', clobber=True)

    pyfits.writeto('udf-candels-f125w.fits',f125[1].data, header=f160[1].header)
    pyfits.writeto('udf-candels-f160w.fits',f160[1].data, header=f160[1].header)
    
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
    
    threedhst.gmap.makeImageMap(['AEGIS-28-F140W_drz.fits', 'AEGIS-28-F140W_align.fits[0]*0.04', 'AEGIS-28-G141_drz.fits', 'AEGIS-9-G141_drz.fits'][0:], aper_list=[13, 14,15], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings, not coded yet but need to regenerate the original 
    #### single drz pointings to avoid growing them each time.
    files = glob.glob('AEGIS-*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.startMultidrizzle(file, 
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        #
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='AEGIS-*-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)
    
def AEGIS_mosaic():
    import threedhst.prep_flt_files
    
    os.chdir(unicorn.GRISM_HOME+'AEGIS/PREP_FLT')
    
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
    PIXFRAC=0.8
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
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['AEGIS-F140W_drz.fits', 'AEGIS-G141_drz.fits','/3DHST/Ancillary/AEGIS/ACS/mos_i_scale2_drz.fits[0]','/3DHST/Ancillary/AEGIS/NMBS/AEGIS-N2_K_sci.fits[0]'], aper_list=zooms, tileroot=['F140W','G141','ACS-i','NMBS-K'], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'))

#
def UDS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'UDS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_drz.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rxyscale')
    
    #### Check
    threedhst.gmap.makeImageMap(['UDS-23-F140W_drz.fits', 'UDS-23-F140W_align.fits[0]','UDS-23-G141_drz.fits'], aper_list=[14,15,16],  polyregions=glob.glob('UDS-*-F140W_asn.pointing.reg'))
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings, not coded yet but need to regenerate the original 
    #### single drz pointings to avoid growing them each time.
    files = glob.glob('UDS-*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.startMultidrizzle(file, 
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        #
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='UDS-*-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=200)

def SN_GEORGE():
    ####********************************************####
    ####              SN-GEORGE (GOODS-S)
    ####********************************************####
    
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair

    os.chdir(unicorn.GRISM_HOME+'SN-GEORGE/PREP_FLT')
    
    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    ##### Direct + grism
    pair('ibfug1030_asn.fits','ibfug1040_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False)
    status = os.system('mv GEORGE-1-F160W_drz.fits GEORGE-2-F160W_drz.fits')
    status = os.system('mv GEORGE-1-F160W_tweak.fits GEORGE-2-F160W_tweak.fits')
    status = os.system('mv GEORGE-1-F160W_asn.fits GEORGE-2-F160W_asn.fits')
    status = os.system('mv GEORGE-1-F160W_shifts.txt GEORGE-2-F160W_shifts.txt')

    pair('ibfug2030_asn.fits','ibfug2040_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False, GET_SHIFT=True)
    status = os.system('mv GEORGE-1-F125W_drz.fits GEORGE-2-F125W_drz.fits')
    status = os.system('mv GEORGE-1-F125W_tweak.fits GEORGE-2-F125W_tweak.fits')
    status = os.system('mv GEORGE-1-F125W_asn.fits GEORGE-2-F125W_asn.fits')
    status = os.system('mv GEORGE-1-F125W_shifts.txt GEORGE-2-F125W_shifts.txt')
    
    
    ##### Direct only
    pair('ibfug3020_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F125W_drz.fits GEORGE-3-F125W_drz.fits')
    status = os.system('mv GEORGE-1-F125W_tweak.fits GEORGE-3-F125W_tweak.fits')
    status = os.system('mv GEORGE-1-F125W_asn.fits GEORGE-3-F125W_asn.fits')
    status = os.system('mv GEORGE-1-F125W_shifts.txt GEORGE-3-F125W_shifts.txt')

    pair('ibfug4020_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F125W_drz.fits GEORGE-4-F125W_drz.fits')
    status = os.system('mv GEORGE-1-F125W_tweak.fits GEORGE-4-F125W_tweak.fits')
    status = os.system('mv GEORGE-1-F125W_asn.fits GEORGE-4-F125W_asn.fits')
    status = os.system('mv GEORGE-1-F125W_shifts.txt GEORGE-4-F125W_shifts.txt')
    
    pair('ibfug5020_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F125W_drz.fits GEORGE-5-F125W_drz.fits')
    status = os.system('mv GEORGE-1-F125W_tweak.fits GEORGE-5-F125W_tweak.fits')
    status = os.system('mv GEORGE-1-F125W_asn.fits GEORGE-5-F125W_asn.fits')
    status = os.system('mv GEORGE-1-F125W_shifts.txt GEORGE-5-F125W_shifts.txt')
    
    #### later F160W don't have ASN files
    asn = threedhst.utils.ASNFile('ibfug5020_asn.fits')

    asn.exposures = ['ibfug3hes','ibfug3hns']
    asn.product = 'ibfug3010'.upper()
    asn.write('ibfug3010_asn.fits')
    
    pair('ibfug3010_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F160W_drz.fits GEORGE-3-F160W_drz.fits')
    status = os.system('mv GEORGE-1-F160W_tweak.fits GEORGE-3-F160W_tweak.fits')
    status = os.system('mv GEORGE-1-F160W_asn.fits GEORGE-3-F160W_asn.fits')
    status = os.system('mv GEORGE-1-F160W_shifts.txt GEORGE-3-F160W_shifts.txt')
    
    
    asn.exposures = ['ibfug4g0q','ibfug4g9q']
    asn.product = 'ibfug4010'.upper()
    asn.write('ibfug4010_asn.fits')

    pair('ibfug4010_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F160W_drz.fits GEORGE-4-F160W_drz.fits')
    status = os.system('mv GEORGE-1-F160W_tweak.fits GEORGE-4-F160W_tweak.fits')
    status = os.system('mv GEORGE-1-F160W_asn.fits GEORGE-4-F160W_asn.fits')
    status = os.system('mv GEORGE-1-F160W_shifts.txt GEORGE-4-F160W_shifts.txt')

    
    asn.exposures = ['ibfug5jrq','ibfug5k0q']
    asn.product = 'ibfug5010'.upper()
    asn.write('ibfug5010_asn.fits')

    pair('ibfug5010_asn.fits', None, ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    status = os.system('mv GEORGE-1-F160W_drz.fits GEORGE-5-F160W_drz.fits')
    status = os.system('mv GEORGE-1-F160W_tweak.fits GEORGE-5-F160W_tweak.fits')
    status = os.system('mv GEORGE-1-F160W_asn.fits GEORGE-5-F160W_asn.fits')
    status = os.system('mv GEORGE-1-F160W_shifts.txt GEORGE-5-F160W_shifts.txt')
       
    ###### Combine everything
    asn_files = glob.glob('GEORGE-?-F125W_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_files, out_root='GEORGE-F125W',
                    path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('GEORGE-F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    asn_files = glob.glob('GEORGE-?-F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_files, out_root='GEORGE-F160W',
                    path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('GEORGE-F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    asn_files = glob.glob('GEORGE-?-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_files, out_root='GEORGE-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('GEORGE-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    #### Original F125W not quite aligned to F160W, which looks closer to ACS
    threedhst.shifts.refine_shifts(ROOT_DIRECT='GEORGE-F125W', 
              ALIGN_IMAGE='GEORGE-F160W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle('GEORGE-F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    threedhst.gmap.makeImageMap(['GEORGE-F125W_drz.fits', 'GEORGE-F160W_drz.fits', 'GEORGE-1-F160W_align.fits[0]*4', 'GEORGE-G141_drz.fits*2'], aper_list=[15], polyregions=glob.glob('GEORGE-1-F125W_asn.pointing.reg'))
    
    
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
    threedhst.gmap.makeImageMap(['MARSHALL-F125W_drz.fits', 'MARSHALL-F160W_drz.fits', 'MARSHALL-1-F160W_align.fits[0]', 'MARSHALL-225-G141_drz.fits', 'MARSHALL-245-G141_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[15], polyregions=glob.glob('MARSHALL-1-F160W_asn.pointing.reg'))
    
def GOODS_ERS():
    ####******************************************************####
    ####              GOODS-S ERS Field
    ####******************************************************####
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import pyfits
    # reload threedhst.grism_sky   #### Reset the flat-field for F140W/G141
    
    os.chdir(unicorn.GRISM_HOME+'ERS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    direct = glob.glob('../RAW/ib*20_asn.fits')
    grism = glob.glob('../RAW/ib*10_asn.fits')
    
    for i in range(len(direct))[:-1]:        
        dir_i=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        gris_i=threedhst.prep_flt_files.make_targname_asn(grism[i], newfile=True)
    
    #### F140W / G141
    threedhst.grism_sky.set_G141()
    pair(direct[1], grism[1], ALIGN_IMAGE = ALIGN, TWEAKSHIFTS_ONLY=False, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False)
    
    #### F098M / G102
    threedhst.grism_sky.set_G102()
    pair(direct[0], grism[0], ALIGN_IMAGE = ALIGN, TWEAKSHIFTS_ONLY=False, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=True, sky_images=['../CONF/G102_master_flatcorr.fits'], second_pass=False, overall=True)
    
    ##### Check results
    threedhst.gmap.makeImageMap(['WFC3-ERSII-G01-F140W_drz.fits', 'WFC3-ERSII-G01-F140W_align.fits[0]*4', 'WFC3-ERSII-G01-F098M_drz.fits', 'WFC3-ERSII-G01-G141_drz.fits','WFC3-ERSII-G01-G102_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[14, 15])
    
#
def SN_PRIMO():
    ####********************************************####
    ####              SN-PRIMO (GOODS-S / UDF)
    ####********************************************####
    
    import os
    import threedhst
    import unicorn
    import unicorn.candels
    
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair

    os.chdir(unicorn.GRISM_HOME+'SN-PRIMO/PREP_FLT')
    
    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    ##### Direct + grism
    pair('ibfup1020_asn.fits','ibfup1030_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False)

    pair('ibfup2020_asn.fits','ibfup2030_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False)
    
    #### Make empty ASNs for 11-01, only has single direct image for each set of 
    #### grism exposures
    asn = threedhst.utils.ASNFile('ibfup1020_asn.fits')
    asn.exposures = ['ibfupacsq','ibfupbdkq']
    asn.product = 'ibfupa010'.upper()
    asn.write('ibfupa010_asn.fits')
    
    pair('ibfupa010_asn.fits','ibfupa020_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False)
    
    files=glob.glob('PRIMO-1101-G141*')
    for file in files:
        out = file.replace('1101','1101-a')
        status = os.system('mv %s %s' %(file, out))
    
    #### Second set
    asn = threedhst.utils.ASNFile('ibfup1020_asn.fits')
    asn.exposures = ['ibfupbdkq', 'ibfupacsq']
    asn.product = 'ibfupb010'.upper()
    asn.write('ibfupb010_asn.fits')
    
    pair('ibfupb010_asn.fits','ibfupb020_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=False)
    
    files=glob.glob('PRIMO-1101-G141*')
    for file in files:
        out = file.replace('1101','1101-b')
        status = os.system('mv %s %s' %(file, out))
       
    ##### Other direct images
    files=glob.glob('ibfup[5-7]*asn.fits')
    files.extend(glob.glob('ibfup[cd]*asn.fits'))
    
    for file in files:
        pair(file,'ibfup2030_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_DIRECT=False, SKIP_GRISM=True)
    
    ##### Combine different sets of direct / grism images
    
    ## something wrong with these
    status = os.system('rm PRIMO-112*')
    status = os.system('rm PRIMO-1027*')
    
    direct_files = glob.glob('PRIMO-*-F125W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='PRIMO-F125W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('PRIMO-F125W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False)
    
    direct_files = glob.glob('PRIMO-*-F160W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='PRIMO-F160W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('PRIMO-F160W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=53.159487, dec=-27.776371,
             final_outnx=2980, final_outny=2380)
    
    direct_files = glob.glob('PRIMO-1101-?-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='PRIMO-1101-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('PRIMO-1101-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False)
    
    
    threedhst.gmap.makeImageMap(['PRIMO-F160W_drz.fits','PRIMO-F125W_drz.fits', 'PRIMO-1026-F160W_align.fits[0]','PRIMO-1026-G141_drz.fits', 'PRIMO-1101-G141_drz.fits'], zmin=-0.06, zmax=0.6, aper_list=[15], polyregions=glob.glob('PRIMO-*-G141_asn.pointing.reg'))
    
    ###################
    ####  Make detection image from CANDELS
    ###################
    os.chdir('/Users/gbrammer/CANDELS/GOODS-S/PREP_FLT')
    unicorn.candels.make_asn_files()
    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    visits = ['GOODS-SD5-VMX','GOODS-SD5-VGX','GOODS-S205-VHS','GOODS-SDW-VH7','GOODS-S075-VCM','GOODS-SD2-V71','GOODS-SD5-VGU','GOODS-S205-VHT','GOODS-S075-VJ3','GOODS-SD5-VGV','GOODS-SD2-V7J','GOODS-S1-V3J']
    for visit in visits:
        if os.path.exists(visit+'-F125W_asn.fits') & (not os.path.exists(visit+'-F125W_drz.fits')):
            print visit
            for filter in ['F125W','F160W']:
                try:
                    unicorn.candels.prep_candels(asn_file=visit+'-'+filter+'_asn.fits', 
                        ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=0,
                        GET_SHIFT=True, DIRECT_HIGHER_ORDER=2,
                        SCALE=0.06)
                except:
                    pass
    
    #### Shifts not right for this visit
    threedhst.shifts.refine_shifts(ROOT_DIRECT='GOODS-SD5-VMX-F160W', 
              ALIGN_IMAGE='GOODS-SD5-VMX-F125W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='shift', clean=True)
    #
    threedhst.prep_flt_files.startMultidrizzle('GOODS-SD5-VMX-F160W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False)
    #
    threedhst.gmap.makeImageMap( ['GOODS-SD5-VMX-F160W_drz.fits','GOODS-SD5-VMX-F125W_drz.fits','GOODS-SD5-VMX-F160W_align.fits[0]*4'], zmin=-0.06, zmax=0.6, aper_list=[15], polyregions=glob.glob('PRIMO-*-G141_asn.pointing.reg'))
    
    ### Combine CANDELS visits into a single image                
    for filter in ['F125W','F160W']:
        drz_files = glob.glob('GOODS-*-'+filter+'_drz.fits')
        direct_files = []
        for drz in drz_files:
            direct_files.append(drz.replace('drz','asn'))
        #
        threedhst.utils.combine_asn_shifts(direct_files, out_root='PRIMO-'+filter,
            path_to_FLT='./', run_multidrizzle=False)
    
        threedhst.prep_flt_files.startMultidrizzle('PRIMO-'+ filter+'_asn.fits',
                use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.6, driz_cr=False,
                updatewcs=False, clean=True, median=False,
                ra=53.159487, dec=-27.776371,
                final_outnx=2980, final_outny=2380)
    
    threedhst.gmap.makeImageMap(['PRIMO-F125W_drz.fits', 'PRIMO-F160W_drz.fits','/3DHST/Spectra/Work/SN-PRIMO/PREP_FLT/PRIMO-F160W_drz.fits','/3DHST/Spectra/Work/SN-PRIMO/PREP_FLT/PRIMO-1101-G141_drz.fits'], zmin=-0.08, zmax=2, aper_list=[16], polyregions=glob.glob('PRIMO-*-G141_asn.pointing.reg'))
    
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

def count_pointings():
    import unicorn
    import glob
    os.chdir(unicorn.GRISM_HOME)
    directories = glob.glob('*')
    for dir in directories:
        cats = glob.glob(dir+'/HTML/*drz.cat')
        print '%s: %0d' %(dir, len(cats))

def copy_drz_images():
    """
    cd /Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/DRZ
    
    for field in AEGIS COSMOS GOODS-N GOODS-S ERS; do
        rsync -avz /3DHST/Spectra/Work/${field}/PREP_FLT/*F140W_drz.fits ./
        rsync -avz /3DHST/Spectra/Work/${field}/DATA/*G141_drz.fits ./
    done
    
    gzip *drz.fits
    
    """
    return True
    
def add_links_to_overview_websites():
    """
    cd /Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/MOSAIC
    
    for field in AEGIS COSMOS GOODS-N GOODS-S; do
        echo $field
        cd $field
        cp scripts/threedhst.js scripts/threedhst.js.bkup
        perl -pi -e "s/\/\/ GEvent/GEvent/" scripts/threedhst.js
        perl -pi -e "s/xxx\//unicorn.astro.yale.edu\/P\/GRISM_v1.6\//" scripts/threedhst.js
        perl -pi -e "s/F140W.yyy/G141_index.html/" scripts/threedhst.js
        cd ../
    done
    
    """
    return True
    
def check_weights():
    """
    Check if MultiDrizzle is still flagging stars as cosmic rays
    """
    
    os.chdir(unicorn.GRISM_HOME+'COSMOS/PREP_FLT/')
    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits[1]', 'COSMOS-F140W_drz.fits[2]*0.0008'], aper_list=[15], zmin=-0.1, zmax=0.5, polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))

    threedhst.gmap.makeImageMap(['COSMOS-23-F140W_drz.fits[1]', 'COSMOS-23-F140W_drz.fits[2]*0.0008'], aper_list=[16], zmin=-0.1, zmax=0.5, polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))
    
    
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

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import os
import glob
import shutil

def check_COSMOS_stars():
    """
    Run SExtractor and try to get the FWHM of the new reduction compared to the old,
    which had many stars with central pixels rejected by the CR rejection.
    """    
    from pyraf import iraf

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
    import unicorn
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair
    import threedhst.prep_flt_files
    import drizzlepac
    import glob
    import os
        
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/INTERLACE_v4.1.1')

    ALIGN = '/3DHST/Photometry/Release/v4.0/GOODS-N/HST_Images/goodsn_3dhst.v4.0.F160W_orig_sci.fits'
    CATALOG = '/3DHST/Photometry/Work/GOODS-N/v4/sextr/catalogs/GOODS-N_IR_radec.cat'
    
    unicorn.candels.make_asn_files(uniquename=True)
        
    
    ### Quality checks on the additional pointings:
    direct_files =  glob.glob('goodsn-*-4*-*-F140W_asn.fits')
    grism_files =  glob.glob('goodsn-*-4*-*-G141_asn.fits')
    for direct, grism in zip(direct_files, grism_files):
        print direct, grism
        threedhst.dq.checkDQ(asn_grism_file=grism, asn_direct_file=direct, path_to_flt='./', size=1000)
        
    
    ### The following direct images require masks:
    ### 'ib3702uuq_flt.fits', 'ib3703uxq_flt.fits', 'ib3704wpq_flt.fits', 'ib3705xzq_flt.fits', 'ib3707c8q_flt.fits', 'ib3708ivq_flt.fits', 'ib3710n4q_flt.fits', 'ib3711c8q_flt.fits', 'ib3712m5q_flt.fits', 'ib3715siq_flt.fits', 'ib3715sxq_flt.fits', 'ib3716q3q_flt.fits', 'ib3716qiq_flt.fits', 'ib3717wyq_flt.fits', 'ib3720fbq_flt.fits', 'ib3721xtq_flt.fits', 'ib3721y8q_flt.fits', 'ib3723r0q_flt.fits', 'ib3727u2q_flt.fits', 'ib3728d9q_flt.fits', 'ib3728doq_flt.fits', 'ib3749obq_flt.fits'

    ## Rename assosciation tables
    grism = 'G141'
    filter = 'F140W'

    ### Rename single exposures
    rename = ['goodsn-15','goodsn-16','goodsn-17','goodsn-18','goodsn-25','goodsn-26','goodsn-27','goodsn-28','goodsn-32', 'goodsn-33', 'goodsn-34', 'goodsn-35', 'goodsn-36', 'goodsn-41', 'goodsn-42', 'goodsn-43', 'goodsn-44', 'goodsn-45', 'goodsn-46']
    for target in rename:
        if len(glob.glob('{}*-??-???-{}_asn.fits'.format(target, grism))) > 1:
            print 'MORE THAN ONE MATCH: {}'.target
            continue
        else:
            file = glob.glob('{}*-??-???-{}_asn.fits'.format(target, grism))[0]
        root = file.split('-{}'.format(grism))[0]
        print root
        asn_grism = threedhst.utils.ASNFile(file)
        asn_grism.product = '{}-{}'.format(target,grism)
        asn_grism.write('{}_asn.fits'.format(asn_grism.product))
        asn_direct = threedhst.utils.ASNFile(file.replace(grism, filter))
        asn_direct.product = '{}-{}'.format(target, filter)
        asn_direct.write('{}_asn.fits'.format(asn_direct.product))
        print 'mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(root, grism)
        os.system('mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(root, grism))
        print 'mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(root, filter)
        os.system('mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(root, filter))
        
    rename_repeats = {'goodsn-11-01-345':'goodsn-11', 'goodsn-11-41-345':'goodsn-111', 'goodsn-14-04-342':'goodsn-14', 'goodsn-14-44-341':'goodsn-114', 'goodsn-23-07-339':'goodsn-23', 'goodsn-23-47-345':'goodsn-123'}
    for target in rename_repeats.keys():
        asn_grism = threedhst.utils.ASNFile('{}-{}_asn.fits'.format(target, grism))
        asn_grism.product = '{}-{}'.format(rename_repeats[target],grism)
        asn_grism.write('{}_asn.fits'.format(asn_grism.product))
        asn_direct = threedhst.utils.ASNFile('{}-{}_asn.fits'.format(target, filter))
        asn_direct.product = '{}-{}'.format(rename_repeats[target],filter)
        asn_direct.write('{}_asn.fits'.format(asn_direct.product))
        print 'mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(target, grism)
        os.system('mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(target, grism))
        print 'mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(target, filter)
        os.system('mv {0}-{1}_asn.fits {0}-{1}_asn.fits.X'.format(target, filter))
            

    direct=glob.glob('goodsn-*-F140W_asn.fits')
    grism = glob.glob('goodsn-*-G141_asn.fits')
 
    for i in range(len(direct)):
        print direct[i], grism[i]
        pair(direct_asn=direct[i], grism_asn=grism[i], radec=CATALOG, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, skip_direct=False, ACS=False, align_threshold=6)

    ### Coadd the ASNs for repeat observations:
    combine = ['goodsn-12', 'goodsn-13', 'goodsn-21', 'goodsn-22', 'goodsn-24', 'goodsn-31']
    grism = 'G141'
    filter = 'F140W'
    
    for target in combine:
        files = glob.glob('{}*-??-???-{}_asn.fits'.format(target, grism))
        print files
        rots = []
        for file in files:
            rots.append(file.split('-{}'.format(grism))[0][-3:])
        for rot in np.unique(rots):
            files = glob.glob('{}*-??-{}-{}_asn.fits'.format(target, rot, grism))
            exposures = []
            direct_exposures = []
            for file in files:
                asn_grism = threedhst.utils.ASNFile(file)
                exposures.extend(asn_grism.exposures)
                print 'mv {0} {0}.X'.format(file)
                os.system('mv {0} {0}.X'.format(file))
                direct = glob.glob(file.split('-{}'.format(grism))[0]+'-{}_asn.fits'.format(filter))[0]
                asn_direct = threedhst.utils.ASNFile(direct)
                direct_exposures.extend(asn_direct.exposures)  
                print 'mv {0} {0}.X'.format(direct)
                os.system('mv {0} {0}.X'.format(direct))
            asn_grism.exposures = exposures
            asn_grism.product = '{}-{}'.format(target, grism)
            asn_grism.write('{}_asn.fits'.format(asn_grism.product), clobber=True)
            asn_direct.exposures = direct_exposures
            asn_direct.product = '{}-{}'.format(target, filter)
            asn_direct.write('{}_asn.fits'.format(asn_direct.product), clobber=True)
            drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(asn_grism.product), 
                clean=True, skysub=False, final_pixfrac=0.8, context=False, resetbits=0, 
                driz_sep_bits=576, final_bits=576, preserve=False, driz_separate=False,     
                driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False)
            
    redo_cr = ['goodsn-12-G141_asn.fits', 'goodsn-13-G141_asn.fits', 'goodsn-21-G141_asn.fits', 'goodsn-22-G141_asn.fits', 'goodsn-24-G141_asn.fits', 'goodsn-31-G141_asn.fits']

    for asn in redo_cr:
        drizzlepac.astrodrizzle.AstroDrizzle(asn, clean=True, skysub=False, final_scale=None, 
        final_pixfrac=0.8, context=False, driz_sep_bits=576, final_bits=576, preserve=False, 
        driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7') 
    
            

    #### Diagnostic HTML:
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
             
    #
    threedhst.shifts.plot_shifts('GOODS-N-F140W', '/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_nz*drz*fits', skip_swarp=False)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False,ra=189.17736, dec=62.23892,
             final_outnx = NX, final_outny=NY, ivar_weights=False)

    #### Make full alignment mosaic
    threedhst.shifts.matchImagePixels(input=glob.glob(ALIGN), matchImage='GOODS-N-F140W_drz.fits', match_extension=1, output='GOODS-N-F850LP.fits')
    
    threedhst.gmap.makeImageMap(['GOODS-N-F140W_drz.fits','GOODS-N-F850LP.fits[0]*4'], aper_list=[15], polyregions=glob.glob('GOODS-N-*-F140W_asn.pointing.reg')) #, zmin=-0.1, zmax=3)
    
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['GOODS-N-F140W_drz.fits', 'GOODS-N-G141_drz.fits*2', 'GOODS-N-F850LP.fits[0]*6','/3DHST/Ancillary/GOODS-N/CDFN/paper13-cdfn-figure3-FB-binned1pix-smooth.fits[0]'], aper_list=zooms, tileroot=['F140W', 'G141', 'GOODS-z850','0.5-8keV'], polyregions=glob.glob('GOODS-N*F140W_asn.pointing.reg')) #, zmin=-0.1, zmax=3)

    zooms=[12,13,14,15,16]
    threedhst.gmap.makeImageMap(['/Volumes/robot/3DHST/Spectra/Work/GOODS-N/Interlace_GBB/GOODS-N-F140W_drz.fits[1]', '/Volumes/robot/3DHST/Spectra/Work/GOODS-N/Interlace_GBB/GOODS-N-G141_drz.fits[1]*2', '/3DHST/Photometry/Work/GOODS-N/DATA/mosaics/goodsn_acsb_trim.fits[0]*15','/3DHST/Photometry/Work/GOODS-N/DATA/mosaics/goodsn_acsz_trim.fits[0]*6','/3DHST/Ancillary/GOODS-N/MIPS/n_mips_1_s1_v0.36_sci.fits[0]*100','/3DHST/Ancillary/GOODS-N/HERSCHEL/L2_GOODS-N_All250_DR1/L2_GOODS-N_image_SMAP250_dr1.fits[1]*20000.','/3DHST/Ancillary/GOODS-N/CDFN/paper13-cdfn-figure3-FB-binned1pix-smooth.fits[0]*8','/3DHST/Ancillary/GOODS-N/VLA/GOODSN_1_4GHz.fits[0]*60000.'], aper_list=zooms, tileroot=['F140W', 'G141', 'ACSb-F435W','ACSz-F850LP','MIPS-24','SPIRE-250','0.5-8keV','VLA'], polyregions=glob.glob('GOODS-N*F140W_asn.pointing.reg'),path='/3DHST/Spectra/Work/GOODS-N/MOSAIC_HTML/')
    
    from pyraf import iraf

    iraf.imcopy('GOODS-N-F140W_drz.fits[1]', '../MOSAIC/GOODS-N-F140w_11-10-06_sci.fits')
    iraf.imcopy('GOODS-N-F140W_drz.fits[2]', '../MOSAIC/GOODS-N-F140w_11-10-06_wht.fits')
    # !tar czvf GOODS-N-F140w_11-09-08.tar.gz GOODS-N-*-F140W_shifts.txt GOODS-N-*-F140W_tweak.fits GOODS-N-*-F140W_asn.fits GOODS-N-F140W_shifts.txt GOODS-N-F140W_asn.fits
    # !mv GOODS-N-F140w_11-09-08* ../MOSAIC
    
    ##### Fix associations of GOODS-N revisits to just include FLTs not affected by the background blowout
    # Targets where the revisit angle was the same can be interlaced
    #
    # Target   Angle_offset
    # GNGRISM11: 0.459
    # GNGRISM12: 0.001
    # GNGRISM13: 0.000
    # GNGRISM14: 0.557
    # GNGRISM21: 0.000
    # GNGRISM22: 0.000
    # GNGRISM23: 6.000
    # GNGRISM24: 0.000
    # GNGRISM31: 0.000
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/PREP_FLT/')
    os.chdir('/3DHST/Spectra/Work/GOODS-N/INTERLACE_v4.0_REVISIT')
    #### 
    info = catIO.Readfile('files.info')
    asn = threedhst.utils.ASNFile('../RAW/ib3701050_asn.fits')
    new_pointings = np.unique(info.targname[(info.date_obs > '2011-01-01') & (info.dec_targ > 62.)])
    for pointing in new_pointings:
        #### Only do pointings that were re-observed at the same ORIENT
        test = (info.targname == pointing)
        angles = np.diff(info.pa_v3[test])
        print '%s: %.3f' %(pointing, np.abs(angles).max())
        if np.abs(angles).max() > 0.1:
            continue
        #
        for filter in ['F140W', 'G141']:
            #### Original images and re-visit images
            print filter
            old = (info.targname == pointing) & (info.date_obs < '2011-01-01') & (info.filter == filter)
            new = (info.targname == pointing) & (info.date_obs > '2011-01-01') & (info.filter == filter)
            asn.exposures = ['' for i in range(4)]
            #### Loop through the four exposures and replace with newer if postarg is the same
            for i in range(4):
                match = info.postarg1[old][i] == info.postarg1[new]
                if match.sum() == 0:
                    asn.exposures[i] = info.file[old][i].replace('_flt.fits.gz','')
                    print info.file[old][i].replace('_flt.fits.gz','')
                    continue
                #
                asn.exposures[i] = info.file[new][match][0].replace('_flt.fits.gz','')
                #### Make a differenc image to see if shifts are OK.  Some show offsets of a pixel maybe
                old_img = pyfits.open('../RAW/'+info.file[old][i])
                new_img = pyfits.open('../RAW/'+info.file[new][match][0])
                diff = old_img[1].data-new_img[1].data
                pyfits.writeto(info.file[old][i].replace('.fits.gz','_diff.fits'), diff, header=old_img[1].header, clobber=True)
                print info.file[old][i].replace('_flt.fits.gz','') + '  ->  ' + info.file[new][match][0].replace('_flt.fits.gz','')
                #### Keep direct images that were OK before
                # if (FILTER == 'F140W')  & (np.median(diff[100:400,100:400]) < 0.1):
                #     asn.exposures.append(info.file[old][i].replace('_flt.fits.gz',''))
                #     print info.file[old][i].replace('_flt.fits.gz','') + ' was OK.  Keep it.'
            #
            #### Write the ASN file
            asn.product=pointing.replace('GNGRISM','GOODS-N-')+'-'+filter
            asn.write(pointing.replace('GNGRISM','GOODS-N-')+'-'+filter+'_asn.fits')
    #
    files = glob.glob('*F140W_asn.fits')
    
    ALIGN = '/Users/brammer/3DHST/Ancillary/Mosaics/goods-n-f160w-astrodrizzle-v4.0_drz_sci.fits'
    ALIGN = '/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f160w-astrodrizzle-v4.0_drz_sci.fits'
    
    for file in files:
        threedhst.prep_flt_files.process_3dhst_pair(file, file.replace('F140W', 'G141'), adjust_targname=False, ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
            
    ### Have cleaned up the bad-sky G141 reads of the revisits that can't be interlaced
    threedhst.prep_flt_files.process_3dhst_pair('../RAW/ib3701050_asn.fits', '../RAW/ib3701060_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
    threedhst.prep_flt_files.process_3dhst_pair('../RAW/ib3704050_asn.fits', '../RAW/ib3704060_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
    threedhst.prep_flt_files.process_3dhst_pair('../RAW/ib3707050_asn.fits', '../RAW/ib3707060_asn.fits', ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
    
    ### This one had a big satellite trail that was masked in the MultiAccum file
    file='GOODS-N-31-F140W_asn.fits'
    threedhst.prep_flt_files.process_3dhst_pair(file, file.replace('F140W', 'G141'), adjust_targname=False, ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
    
    files = glob.glob('*F140W_asn.fits')
    for file in files:
        unicorn.reduce.interlace_combine(file.split('_asn')[0], view=False, NGROW=125)
        unicorn.reduce.interlace_combine(file.split('_asn')[0].replace('F140W', 'G141'), view=False, NGROW=125)
    
    #### Shifts of new G141 images from misalignment of revisit
    # files=glob.glob('/tmp/OLD/*G141_inter.fits.gz')
    # for i in range(len(files)):
    #     file = files[i]
    #     old = pyfits.open(file)
    #     new = pyfits.open(os.path.basename(file).split('.gz')[0])
    #     shift = new[1].data*0.
    #     shift[:-1,:] += new[1].data[1:,:]
    #     #shift[1:,:] += new[1].data[:-1,:]
    #     #shift[:,:] += new[1].data[:,:]
    #     ds9.view(old[1].data-shift)
        
    # Just these two need y shifts
    im = pyfits.open('GOODS-N-22-G141_inter.fits', mode='update')
    im[1].data[:-1,:] = im[1].data[1:,:]*1
    im[2].data[:-1,:] = im[2].data[1:,:]*1
    im.flush()
    
    im = pyfits.open('GOODS-N-31-G141_inter.fits', mode='update')
    im[1].data[:-1,:] = im[1].data[1:,:]*1
    im[2].data[:-1,:] = im[2].data[1:,:]*1
    im.flush()
    
def COSMOS(FORCE=False):
    import unicorn
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
        
    #os.chdir(unicorn.GRISM_HOME+'COSMOS/INTERLACE_v4.1.1')
    ALIGN = '/3DHST/Photometry/Release/v4.0/COSMOS/HST_Images/cosmos_3dhst.v4.0.F160W_orig_sci.fits'
    CATALOG = '/3DHST/Photometry/Work/COSMOS/Sex/cosmos_3dhst.v4.0.IR_radec.cat'
    
    unicorn.candels.make_asn_files(uniquename=False)
        
    direct=glob.glob('cosmos-*-F140W_asn.fits')
    grism = glob.glob('cosmos-*-G141_asn.fits')
    
    ### The following direct images have masks:
    ### ibhm40brq_flt.fits, ibhm45xzq_flt.fits, ibhm45ziq_flt.fits, ibhm45zqq_flt.fits, ibhm47h8q_flt.fits, ibhm53onq_flt.fits
    loop_list = range(len(direct))
    FORCE=True
    for i in loop_list:
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct_asn=direct[i], grism_asn=grism[i], radec=CATALOG, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, skip_direct=False, ACS=False, align_threshold=6)

    """
    # Diagnostic HTML maps:
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
    """

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
    
    threedhst.shifts.plot_shifts('COSMOS-F140W', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits', clean=True, verbose=True, ALIGN_EXTENSION=0, skip_swarp=True, threshold=10)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12634, dec=2.3336697,
             final_outnx = NX, final_outny=NY, ivar_weights=False)
    
    ### Check the mosaic
    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits', 'COSMOS-G141_drz.fits', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04'], aper_list=[14], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))
    
    threedhst.shifts.plot_shifts('COSMOS-F140W', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_*_sci.fits', skip_swarp=False)
    
    ### Temporary release of the direct mosaic
    from pyraf import iraf

    iraf.imcopy('COSMOS-F140W_drz.fits[1]',
        '../MOSAIC/COSMOS-F140w_11-10-06_sci.fits')
    iraf.imcopy('COSMOS-F140W_drz.fits[2]',
        '../MOSAIC/COSMOS-F140w_11-10-06_wht.fits')
    # !tar czvf COSMOS-F140w_11-09-08.tar.gz COSMOS-*-F140W_shifts.txt COSMOS-*-F140W_tweak.fits COSMOS-*-F140W_asn.fits COSMOS-F140W_shifts.txt COSMOS-F140W_asn.fits
    # !mv COSMOS-F140w_11-09-08* ../MOSAIC
        
    zooms = [13,14,15]
    threedhst.gmap.makeImageMap(['../MOSAIC/COSMOS-F140w_11-07-31_sci.fits[0]', '/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/Chandra/CC0570_img.fits[0]*1.5', '/3DHST/Ancillary/COSMOS/Spitzer/mips_24_GO3_sci_10.fits[0]*400', '/3DHST/Ancillary/COSMOS/VLA/vla_20cm_dp_sin_10.fits[0]*40000'], aper_list=zooms, tileroot=['F140W','ACS','WIRDS-K','0.5-7keV', 'MIPS-24', 'VLA-20cm'], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'))

    zooms = [12,13,14,15,16]
    threedhst.gmap.makeImageMap(['/3DHST/Spectra/Work/COSMOS/INTERLACE/COSMOS-F140W_drz.fits[1]','/3DHST/Spectra/Work/COSMOS/INTERLACE/COSMOS-G141_drz.fits[1]', '/3DHST/Ancillary/COSMOS/CANDELS/Gabe/COSMOS-full-F160W_drz_sci.fits[0]','/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.04', '/3DHST/Ancillary/COSMOS/Chandra/CC0570_img.fits[0]*1.5', '/3DHST/Ancillary/COSMOS/Spitzer/mips_24_GO3_sci_10.fits[0]*400', '/3DHST/Ancillary/COSMOS/VLA/vla_20cm_dp_sin_10.fits[0]*60000'], aper_list=zooms, tileroot=['F140W','G141','F160W','ACS 814W','WIRDS-Ks','0.5-7keV', 'MIPS-24', 'VLA-20cm'], polyregions=glob.glob('COSMOS-*-F140W_asn.pointing.reg'),path='/3DHST/Spectra/Work/COSMOS/MOSAIC_HTML/')
    
    ### don't need high-res tiles of lo-res images sitting around
    os.system('rm ~/Sites/FITS/tiles/MIPS*1[67].png')
    os.system('rm ~/Sites/FITS/tiles/VLA*1[7].png')
    os.system('rm ~/Sites/FITS/tiles/0.5*17.png')
    
    #### Make full ACS mosaics matched to the UCSC candels image
    os.chdir('/Volumes/robot/3DHST/Ancillary/COSMOS/ACS/MOSAIC')
    MATCH = '../../CANDELS/Gabe/COSMOS-full-F125W_drz_sci.fits'

    input = glob.glob('../*sci.fits')
    threedhst.shifts.matchImagePixels(input=input, matchImage=MATCH, output='COSMOS_ACSi_sci.fits', input_extension=0, match_extension=0)

    input = glob.glob('../*wht.fits')
    threedhst.shifts.matchImagePixels(input=input, matchImage=MATCH, output='COSMOS_ACSi_wht.fits', input_extension=0, match_extension=0)
    
def GOODSS(FORCE=False):
    import unicorn
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair
    import glob
    import os

    #os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE_v4.1.5')
        
    ALIGN = '/3DHST/Photometry/Release/v4.0/GOODS-S/HST_Images/goodss_3dhst.v4.0.F160W_orig_sci.fits'
    CATALOG = '/3DHST/Photometry/Work/GOODS-S/v4/sextr/catalogs/GOODS-S_IR_radec.cat'
    
    #### Main preparation loop
    ### The following direct images have satellite tracks, make sure you have masks:
    ### ibhj04gbq_flt.fits, ibhj07yzq_flt.fits, ibhj28i0q_flt.fits, ibhj29naq_flt.fits
    
    unicorn.candels.make_asn_files(uniquename=False)
    
    direct=glob.glob('goodss-*-F140W_asn.fits')
    grism = glob.glob('goodss-*-G141_asn.fits')
    
    
    loop_list = range(len(direct))
    FORCE=True
    for i in loop_list:
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct_asn=direct[i], grism_asn=grism[i], radec=CATALOG, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=None, skip_direct=False, ACS=False, align_threshold=6)
            

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
    
    SCALE = 0.06
    #SCALE = 0.5
    PIXFRAC=0.8
    NX, NY, ra, dec = int(16210*0.06/SCALE), int(18100*0.06/SCALE), 53.154223, -27.807325
    #matching the Magee CANDELS reduction
    #SCALE = 0.06
    #PIXFRAC=0.8
    #NX, NY, ra, dec = 14867, 18341, 53.131588, -27.808022
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-S-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec,
             final_outnx = NX, final_outny=NY)
    
    threedhst.prep_flt_files.startMultidrizzle('GOODS-S-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=53.154223, dec=-27.807325,
             final_outnx = NX, final_outny=NY)
    
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['GOODS-S-F140W_drz.fits', 'GOODS-S-G141_drz.fits','/3DHST/Ancillary//GOODS-S/CANDELS/hlsp_candels_hst_wfc3_gsd01_f125w_v0.5_drz.fits[0]','/3DHST/Ancillary/GOODS-S/CDFS/CDFS-4Ms-0p5to8-asca-im-bin1.fits[0]'], aper_list=zooms, tileroot=['F140W','G141','F125W','0.5-8keV'], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'))

    zooms=[12,13,14,15,16]
    threedhst.gmap.makeImageMap(['/Volumes/robot/3DHST/Spectra/Work/GOODS-S/INTERLACE/GOODS-S-F140W_drz.fits[1]', '/Volumes/robot/3DHST/Spectra/Work/GOODS-S/INTERLACE/GOODS-S-G141_drz.fits[1]','/3DHST/Ancillary/GOODS-S/CANDELS/candels_release/hlsp_candels_hst_wfc3_gsd04_f160w_v0.5_drz.fits[0]','/3DHST/Ancillary/GOODS-S/ISAAC/GOODS_ISAAC_mosaic_KS_V2.0.fits[0]','/3DHST/Photometry/Work/GOODS-S/DATA/mips/s_mips_1_s1_v0.30_sci.fits[0]*100','/3DHST/Ancillary/GOODS-S/HERSCHEL/L2_ECDFS_All250_DR1/L2_ECDFS_image_SMAP250_dr1.fits[1]*50000.','/3DHST/Ancillary/GOODS-S/CDFS/CDFS-4Ms-0p5to8-asca-im-bin1.fits[0]','/3DHST/Ancillary/GOODS-S/VLA/ECDFS_DR1.FITS[0]*60000.'], aper_list=zooms, tileroot=['F140W','G141','F160W','ISAAC Ks','MIPS-24','SPIRE-250','0.5-8keV','VLA'], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'),path='/3DHST/Spectra/Work/GOODS-S/MOSAIC_HTML/')

    
def UDF():
    import threedhst
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'UDF/PREP_FLT')

    ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    
    direct_files = glob.glob('GOODS-S-3?-F140W*asn.fits')
    grism_files = glob.glob('GOODS-S-3?-G141*asn.fits')
    for direct, grism in zip(direct_files, grism_files):
        pair(direct, grism, ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
    #### Combinations for noise tests
    #direct_files = glob.glob('GOODS-S-3[4678]-F140W_asn.fits')
    asn_files = glob.glob('GOODS-S*asn.fits')
    for file in asn_files:
        threedhst.utils.combine_asn_shifts([file], out_root=file.replace('GOODS-S','UDF').replace('_asn.fits',''), path_to_FLT='./', run_multidrizzle=False)
    #
    threedhst.utils.combine_asn_shifts(asn_files[0:4:2], out_root='UDF-A-F140W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(asn_files[1:4:2], out_root='UDF-A-G141', path_to_FLT='./', run_multidrizzle=False)

    threedhst.utils.combine_asn_shifts(asn_files[4::2], out_root='UDF-B-F140W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(asn_files[5::2], out_root='UDF-B-G141', path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.utils.combine_asn_shifts(asn_files[0::2], out_root='UDF-Full-F140W', path_to_FLT='./', run_multidrizzle=False)
    threedhst.utils.combine_asn_shifts(asn_files[1::2], out_root='UDF-Full-G141', path_to_FLT='./', run_multidrizzle=False)
    
    direct_images = glob.glob('UDF*F140W_asn.fits')
    for direct in direct_images:
        threedhst.prep_flt_files.startMultidrizzle(direct,
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    ######## From here down was for the original quick analysis I did with the UDF 
    threedhst.gmap.makeImageMap(['GOODS-S-34-F140W_drz.fits', 'GOODS-S-34-F140W_align.fits[0]*4'], aper_list=[15,16], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'), zmin=-0.1, zmax=1)
    
    for file in files:
        threedhst.shifts.refine_shifts(ROOT_DIRECT=file.split('_asn')[0], 
              ALIGN_IMAGE='GOODS-S-34-F140W_drz.fits', ALIGN_EXTENSION=1,  
              fitgeometry='shift', clean=True)
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
    
    #ra0, dec0 = 53.158733, -27.785316
    ra0, dec0 = 53.159816, -27.784101
    #NX, NY = 1083*0.2/0.06, 1106*0.2/0.06
    NX, NY = 3380, 3286
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
                        SCALE=0.06, geometry='rotate,shift')
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
                    SCALE=0.06, geometry='rotate,shift')
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
    import unicorn
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir(unicorn.GRISM_HOME+'AEGIS/INTERLACE_v4.1.5')
        
    ALIGN = '/3DHST/Photometry/Release/v4.0/AEGIS/HST_Images/aegis_3dhst.v4.0.F160W_orig_sci.fits'
    CATALOG = '/3DHST/Photometry/Work/AEGIS/Sex/aegis_3dhst.v4.0.IR_radec.cat'
    
    ### The following direct images have masks:
    ### ibhj40hvq_flt.fits, ibhj49agq_flt.fits, ibhj49boq_flt.fits, ibhj50enq_flt.fits ibhj50euq_flt.fits ibhj50f1q_flt.fits ibhj52brq_flt.fits ibhj53jbq_flt.fits ibhj66ddq_flt.fits
    
    unicorn.candels.make_asn_files(uniquename=False)
    
    direct=glob.glob('aegis-*-F140W_asn.fits')
    grism = glob.glob('aegis-*-G141_asn.fits')
    
    loop_list = range(len(direct))
    FORCE=True
    for i in loop_list:
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct_asn=direct[i], grism_asn=grism[i], radec=CATALOG, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, skip_direct=False, ACS=False, align_threshold=6)
                        
          
    ### Coadd the ASNs for repeat observations:
    target = 'aegis-01'
    grism = 'G141'
    filter = 'F140W'
    
    files = glob.glob('{}*-??-???-{}_asn.fits'.format(target, grism))
    print files
    rots = []
    for file in files:
        rots.append(file.split('-{}'.format(grism))[0][-3:])
    print rots
    #   
    for rot in np.unique(rots):
        files = glob.glob('{}*-??-{}-{}_asn.fits'.format(target, rot, grism))
        exposures = []
        direct_exposures = []
        for file in files:
            asn_grism = threedhst.utils.ASNFile(file)
            exposures.extend(asn_grism.exposures)
            print 'mv {0} {0}.X'.format(file)
            os.system('mv {0} {0}.X'.format(file))
            direct = glob.glob(file.split('-{}'.format(grism))[0]+'-{}_asn.fits'.format(filter))[0]
            asn_direct = threedhst.utils.ASNFile(direct)
            direct_exposures.extend(asn_direct.exposures)  
            print 'mv {0} {0}.X'.format(direct)
            os.system('mv {0} {0}.X'.format(direct))
        asn_grism.exposures = exposures
        asn_grism.product = '{}-{}'.format(target, grism)
        asn_grism.write('{}_asn.fits'.format(asn_grism.product), clobber=True)
        asn_direct.exposures = direct_exposures
        asn_direct.product = '{}-{}'.format(target, filter)
        asn_direct.write('{}_asn.fits'.format(asn_direct.product), clobber=True)
        drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(asn_grism.product), 
            clean=True, skysub=False, final_pixfrac=0.8, context=False, resetbits=0, 
            driz_sep_bits=576, final_bits=576, preserve=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False)
            
    redo_cr = ['aegis-01-G141_asn.fits','aegis-01-F140W_asn.fits']

    for asn in redo_cr:
        drizzlepac.astrodrizzle.AstroDrizzle(asn, clean=True, skysub=False, final_scale=None, 
        final_pixfrac=0.8, context=False, driz_sep_bits=576, final_bits=576, preserve=False, 
        driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7') 
          
    ### Fix the backgrounds for pointing which have bad background
    pointings = ['aegis-02-F140W','aegis-07-F140W','aegis-08-F140W','aegis-09-F140W']
    for pointing in pointings:
        threedhst.prep_flt_astrodrizzle.subtract_flt_background(pointing,scattered_light=False, 
            sex_background=True)
        drizzlepac.astrodrizzle.AstroDrizzle('{}_asn.fits'.format(pointing), clean=False, context=False, 
            preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, 
            driz_cr=True, driz_cr_corr=True, driz_combine=True)

    """
    # Make diagnostic HTML maps:
    threedhst.gmap.makeImageMap(['AEGIS-1-F140W_drz.fits', 'AEGIS-1-F140W_align.fits[0]*4', 'AEGIS-1-G141_drz.fits'][0:], aper_list=[15, 16], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'), tileroot=['wfc3','acs','g141'] )
    threedhst.gmap.makeImageMap(['AEGIS-1-F140W_drz.fits', 
'AEGIS-1-F140W_align.fits[0]*4', 
'AEGIS-1-G141_drz.fits'][0:], aper_list=[15, 16], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'), tileroot=['wfc3','acs','g141'] )

    
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
    ##
    os.chdir('/3DHST/Spectra/Work/AEGIS/REVISIT/')
    #### Shifts of REVISIT
    new = pyfits.open('../RAW/ibhj69heq_flt.fits.gz')[1].data
    old = pyfits.open('../RAW/ibhj39usq_flt.fits.gz')[1].data
    ### dx, dy = -2, 0
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/AEGIS/REVISIT/MultiAccum')
    bad=1024
    unicorn.prepare.make_mask_crr_file(bad_value=bad)
    unicorn.prepare.flag_bad_reads('ibhj39uuq_raw.fits', ds9=ds9, bad_value=bad)
    unicorn.prepare.flag_bad_reads('ibhj39viq_raw.fits', ds9=ds9, bad_value=bad)
    unicorn.prepare.flag_bad_reads('ibhj69hgq_raw.fits', ds9=ds9, bad_value=bad)
    unicorn.prepare.flag_bad_reads('ibhj69hnq_raw.fits', ds9=ds9, bad_value=bad)
    
    #### Shift the later visit
    os.system("cp ../../RAW/ibhj69heq_flt.fits.gz ../../RAW/ibhj69hlq_flt.fits.gz .; gunzip *.gz")
    files=glob.glob('ibhj69*flt.fits')
    for file in files:
        flt = pyfits.open(file, mode='update')
        for ext in range(1,6):
            flt[ext].data[:,0:-2] = flt[ext].data[:,2:]
        #
        flt.flush()
    #
    os.system('cp *flt.fits ../../RAW/')
    
    #### Make ASN files
    os.chdir('/3DHST/Spectra/Work/AEGIS/REVISIT/')
    info = catIO.Readfile('files.info')
    asn = threedhst.utils.ASNFile('../RAW/ibhj70030_asn.fits')
    for filter in ['F140W', 'G141']:
        sel = (info.targname == 'AEGIS-1') & (info.filter == filter) & (info.date_obs < '2013-01-01')
        asn.exposures = []
        for file in info.file[sel]:
            asn.exposures.append(file.split('_flt')[0])
        #
        asn.product = 'AEGIS-1-%s' %(filter)
        asn.write('AEGIS-1-%s_asn.fits' %(filter))
        print asn.product
    
    ALIGN = '/Users/brammer/3DHST/Ancillary/Mosaics/aegis-f160w-astrodrizzle-v4.0_drz_sci.fits'
    ALIGN = '/3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f160w-astrodrizzle-v4.0_drz_sci.fits'
    
    threedhst.prep_flt_files.process_3dhst_pair('AEGIS-1-F140W_asn.fits', 'AEGIS-1-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate, shift')
    unicorn.reduce.interlace_combine(root='AEGIS-1-F140W', view=False, use_error=True, make_undistorted=False, pad=60, NGROW=125, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=False)
    unicorn.reduce.interlace_combine(root='AEGIS-1-G141', view=False, use_error=True, make_undistorted=False, pad=60, NGROW=125, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=False)
    """

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
        
    SCALE = 0.06
    #SCALE = 0.5
    PIXFRAC=0.8
    NX, NY = int(20047*0.06/SCALE), int(20900*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('AEGIS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=214.86, dec=52.84,
             final_outnx = NX, final_outny=NY)
    

    threedhst.prep_flt_files.startMultidrizzle('AEGIS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=214.86, dec=52.84,
             final_outnx = NX, final_outny=NY)

    #Mosaic which matches the v0.5 CANDELS image
    #SCALE=0.06
    #PIXFRAC=0.8
    #NX,NY,ra,dec = 38500,12600, 214.98054, 52.91692
    #threedhst.prep_flt_files.startMultidrizzle('AEGIS-F140W_asn.fits',
    #         use_shiftfile=True, skysub=False,
    #         final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
    #        updatewcs=False, clean=True, median=False,
    #         ra=ra, dec=dec, final_rot=-49.5760,
    #         final_outnx = NX, final_outny=NY,
    #         refimage='/3DHST/Ancillary/AEGIS/CANDELS/hlsp_candels_hst_wfc3_egsa01_f160w_v0.5_drz.fits')


    #
    zooms=[13,14,15]
    threedhst.gmap.makeImageMap(['AEGIS-F140W_drz.fits', 'AEGIS-G141_drz.fits','/3DHST/Ancillary/AEGIS/ACS/mos_i_scale2_drz.fits[0]','/3DHST/Ancillary/AEGIS/NMBS/AEGIS-N2_K_sci.fits[0]'], aper_list=zooms, tileroot=['F140W','G141','ACS-i','NMBS-K'], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'))

    zooms=[12,13,14,15,16]
    threedhst.gmap.makeImageMap(['AEGIS-F140W_drz.fits[1]*5.','AEGIS-G141_drz.fits[1]*10.'], aper_list=zooms, tileroot=['F140W','G141'], polyregions=glob.glob('AEGIS-*-F140W_asn.pointing.reg'),path='/3DHST/Spectra/Work/AEGIS/MOSAIC_HTML/')

#
def UDS(FORCE=False):

    import unicorn
    from threedhst.prep_flt_astrodrizzle import prep_direct_grism_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
        
    #os.chdir(unicorn.GRISM_HOME+'UDS/INTERLACE_v4.1.5')

    ALIGN = '/3DHST/Photometry/Release/v4.0/UDS/HST_Images/uds_3dhst.v4.0.F160W_orig_sci.fits'
    CATALOG = '/3DHST/Photometry/Work/UDS/v4/sextr/catalogs/UDS_IR_radec.cat'
    
    unicorn.candels.make_asn_files(uniquename=False)
        
    direct=glob.glob('uds-*-F140W_asn.fits')
    grism = glob.glob('uds-*-G141_asn.fits')
    
    ### The following direct images have masks:
    ### ibhm12ubq_flt.fits, ibhm12udq_flt.fits, ibhm12uiq_flt.fits, ibhm12ukq_flt.fits, ibhm23glq_flt.fits, ibhm24fkq_flt.fits, ibhm26n5q_flt.fits, ibhm27qnq_flt.fits, ibhm28slq_flt.fits
    ### ibhm26n0q_flt.fits, ibhm26n7q_flt.fits, ibhm26neq_flt.fits, ibhm26noq_flt.fits, ibhm25m8q_flt.fits, ibhm26n7q_flt.fits, ibhm25mmq_flt.fits, ibhm25mtq_flt.fits 
    loop_list = range(len(direct))
    FORCE=True
    for i in loop_list:
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=True)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct_asn=direct[i], grism_asn=grism[i], radec=CATALOG, raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=0.06, skip_direct=False, ACS=False)

    """
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
    """
    
def UDS_mosaic():
    import threedhst.prep_flt_files
    
    os.chdir(unicorn.GRISM_HOME+'UDS/PREP_FLT')
    
    #### Direct mosaic
    direct_files = glob.glob('UDS-*-F140W_asn.fits')
    direct_files.remove('UDS-18-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='UDS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Gris mosaic
    direct_files = glob.glob('UDS-*-G141_asn.fits')
    direct_files.remove('UDS-18-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='UDS-G141',
                       path_to_FLT='./', run_multidrizzle=False)
    
    #### Make the image
    PIXFRAC=0.8
    SCALE = 0.06
    ra, dec, nx, ny = 34.4081405556, -5.20018201235, int(10422*0.128254/SCALE), int(4349*0.128254/SCALE)
        
    threedhst.prep_flt_files.startMultidrizzle('UDS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec,
             final_outnx = nx, final_outny=ny, build_drz=False)
    
    threedhst.prep_flt_files.startMultidrizzle('UDS-G141_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=ra, dec=dec,
             final_outnx = nx, final_outny=ny, build_drz=False)

    ### To match CANDELS v1.0 mosaic
    #PIXFRAC = 0.8
    #SCALE = 0.06
    #ra, dec, nx, ny = 34.350027, -5.2000058, 30720, 12800
    #threedhst.prep_flt_files.startMultidrizzle('UDS-F140W_asn.fits',
    #         use_shiftfile=True, skysub=False,
    #         final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
    #         updatewcs=False, clean=True, median=False,
    #         ra=ra, dec=dec,
    #         final_outnx = nx, final_outny=ny)

    zooms=[12,13,14,15,16]
    threedhst.gmap.makeImageMap(['UDS-F140W_drz_sci.fits[0]','UDS-G141_drz_sci.fits[0]','/3DHST/Ancillary/UDS/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits[1]'], aper_list=zooms, tileroot=['F140W','G141','F160W'], polyregions=glob.glob('GOODS-S-*-F140W_asn.pointing.reg'),path='/3DHST/Spectra/Work/UDS/MOSAIC_HTML/')

    #### Problems with UDS alignment pointings
    threedhst.shifts.refine_shifts(ROOT_DIRECT='UDS-16-F140W', 
              ALIGN_IMAGE='UDS_F140W_ref.fits', ALIGN_EXTENSION=0,  
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle('UDS-16-F140W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
def Koekemoer():
    """ High-z AGN from Koekemoer program """
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.catIO as catIO
    import glob
    
    os.chdir(unicorn.GRISM_HOME+'Koekemoer/PREP_FLT')
    
    ALIGN = '/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits'
    ALIGN = '/Users/brammer/3DHST/Ancillary/Mosaics/goods-s-f160w-astrodrizzle-v4.0_drz_sci.fits'
    
    ALIGN_EXT=0
    
    info = catIO.Readfile('files.info')
    
    unicorn.candels.make_asn_files(make_region=False)
    
    #### Make ASN files for different grisms / orientations
    asn = threedhst.utils.ASNFile(glob.glob('../RAW/i*asn.fits')[0])
    files = glob.glob('*F140W_asn.fits')
    
    for file in files[1:]:
        threedhst.prep_flt_files.process_3dhst_pair(file, file.replace('F140W','G141'), adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
        threedhst.prep_flt_files.process_3dhst_pair(file, file.replace('F140W','G102'), adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate,shift', sky_images=['sky_g102_f105w.fits'])
    
    ## Use v4.0 as reference
    REF_ROOT='GS-F125W'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F125W', REFERENCE='goodss_3dhst.v4.0.F125W_orig_sci.fits', SEGM = 'goodss_3dhst.v4.0.F160W_seg.fits', sci_extension=0)
    
    #### Make interlaced images
    for file in files[1:]:
        NGROW=125
        CATALOG='goodss_3dhst.v4.0.F125W_conv.cat'
        CATALOG='goodss_3dhst.v4.0.F125W_F140W_F160W_det.cat'
        pointing = file.split('_asn.fits')[0]
        #unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing, NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing, view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, auto_offsets=True)
        #
        unicorn.reduce.interlace_combine(pointing, pad=60, NGROW=NGROW, auto_offsets=True)
        unicorn.reduce.interlace_combine(pointing.replace('F140W','G141'), pad=60, NGROW=NGROW, auto_offsets=True)
        unicorn.reduce.interlace_combine(pointing.replace('F140W','G102'), pad=60, NGROW=NGROW, auto_offsets=True)
        
    #### Need to change the MAG_AUTO fluxes from combined to F125W
    import research.v4 as cat2
    cat2.read_catalogs('GOODS-S')
    for file in files:
        pointing = file.split('-F140W_asn.fits')[0]
        c = threedhst.sex.mySexCat(pointing+'_inter.cat')
        c.renameColumn('MAG_AUTO', 'MAG_COMBINE')
        jmag = np.clip(25-2.5*np.log10(cat2.cat.f_F125W[c.id-1]), 0,30)
        c.addColumn(jmag, format='%.4f', name='MAG_AUTO')
        os.system('cp %s %s.OLD' %(c.filename, c.filename))
        c.write()
        
    #### G141 and G102
    for grism in ['G102','G141']:
        #os.mkdir(grism+'_Nor')
        os.chdir(grism+'_Nor'*0)
        os.system('ln -s ../CDFS-AGN2*%s* .' %(grism))
        os.system('ln -s ../CDFS-AGN2*F140W* .')
        os.system('ln -s ../CDFS-AGN2*F140W* .')
        os.system('ln -s ../CDFS-AGN2-??-???_*inter* .')
        os.chdir('../')
    
    grism = 'G141'
    
    os.chdir(grism)
    files = glob.glob('*F140W_asn.fits')
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=False, force=True)
    for file in files[0:1]:
        pointing = file.split('-F140W_asn.fits')[0]
        model = unicorn.reduce.process_GrismModel(root=pointing, grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, REFINE_MAG_LIMIT=22, make_zeroth_model=True, use_segm=False, model_slope=0, direct='F140W', grism=grism)
    
    ### Use Nor's test config files
    os.chdir(grism+'_Nor')
    import unicorn
    import glob
    unicorn.reduce.set_grism_config(grism=grism, chip=1, use_new_config=True, force=True)
    files = glob.glob('*F140W_asn.fits')
    for file in files[0:1]:
        pointing = file.split('-F140W_asn.fits')[0]
        model = unicorn.reduce.process_GrismModel(root=pointing, grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, REFINE_MAG_LIMIT=22, make_zeroth_model=True, use_segm=False, model_slope=0, direct='F140W', grism=grism)
    
    
    for targ in [1,2]:
        for angle in [83, 91]:
            pointing='CDFS-AGN%d-%d' %(targ, angle)
            model = unicorn.reduce.process_GrismModel(root=pointing, grow_factor=2, growx=2, growy=2, MAG_LIMIT=27, REFINE_MAG_LIMIT=23, make_zeroth_model=True, use_segm=False, model_slope=0, grism='G141')
    
    ### Run the fits!
    model = unicorn.reduce.process_GrismModel('CDFS-AGN1-91', grism='G141')
    test = model.cat.mag < 22
    for id in model.objects[test]:
        model.twod_spectrum(id)
        gris = unicorn.interlace_fit.GrismSpectrumFit('CDFS-AGN1-91_%05d' %(id), use_mag_prior=True)
        gris.fit_in_steps()
    
    ### Check automated interlace offsets
    # unicorn.reduce.get_interlace_offsets('CDFS-AGN1-83-F140W_asn.fits', growx=2, growy=2, path_to_flt='./', verbose=True, first_zero=False, raw=False)
    
    unicorn.reduce.interlace_combine('CDFS-AGN1-83-F140W', NGROW=125, auto_offsets=True, view=False)
    unicorn.reduce.interlace_combine('CDFS-AGN1-83-G141', NGROW=125, auto_offsets=True, view=False)
    unicorn.reduce.interlace_combine('CDFS-AGN1-83-G102', NGROW=125, auto_offsets=True, view=False)
    
    unicorn.reduce.reduce_pointing(file='CDFS-AGN1-83-G141_asn.fits', clean_all=False, clean_spectra=True, make_images=False, make_model=True, fix_wcs=True, extract_limit=23.5, skip_completed_spectra=True, MAG_LIMIT=27.5, out_path='./')
    
def SN_TILE41():
    """ Supernova in COSMOS """
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.catIO as catIO
    
    os.chdir(unicorn.GRISM_HOME+'SN-TILE41/PREP_FLT')
    
    ALIGN = os.getenv('CANDELS')+'/COSMOS/PREP_FLT/COSMOS-full-F160W_drz_sci.fits'
    ALIGN_EXT=0
    
    ALIGN = '/Users/brammer/3DHST/Ancillary/Mosaics/cosmos-f160w-astrodrizzle-v4.0_drz_sci.fits'
    ALIGN_EXT=0
    
    info = catIO.Readfile('files.info')
    
    #### Make ASN files for different grisms / orientations, combine mixed filters for direct images to get the shifts right
    asn = threedhst.utils.ASNFile(glob.glob('../RAW/i*asn.fits')[0])
    
    match = (info.targname == 'TILE41') & (info.filter == 'F160W') & (info.pa_v3 < 140)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-132-F160W'
    asn.write(asn.product+'_asn.fits', clobber=True)

    match = (info.targname == 'TILE41') & (info.filter == 'G141') & (info.pa_v3 < 140)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-132-G141'
    asn.write(asn.product+'_asn.fits', clobber=True)

    match = (info.targname == 'TILE41') & (info.filter != 'G141') & (info.pa_v3 > 140)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-152-direct'
    asn.write(asn.product+'_asn.fits', clobber=True)

    match = (info.targname == 'TILE41') & (info.filter == 'G141') & (info.pa_v3 > 140)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-152-G141'
    asn.write(asn.product+'_asn.fits', clobber=True)

    match = (info.targname == 'TILE41') & (info.filter != 'G102') & (info.pa_v3 < 130)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-129-direct'
    asn.write(asn.product+'_asn.fits', clobber=True)

    match = (info.targname == 'TILE41') & (info.filter == 'G102') & (info.pa_v3 < 130)
    asn.exposures = []
    for exp in info.file[match]:
        asn.exposures.append(exp.split('_flt')[0])
    asn.product = 'TILE41-129-G102'
    asn.write(asn.product+'_asn.fits', clobber=True)
    
    #unicorn.candels.make_asn_files(force=False)
    
    #### First G141 angle
    ## direct only
    pair('TILE41-132-F160W_asn.fits', 'TILE41-132-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
    threedhst.shifts.make_grism_shiftfile('TILE41-132-F160W_asn.fits', 'TILE41-132-G141_asn.fits')
    sf_direct = threedhst.shifts.ShiftFile('TILE41-132-F160W_shifts.txt')
    sf_grism = threedhst.shifts.ShiftFile('TILE41-132-G141_shifts.txt')
    
    ### Two grism exposures at each position
    for i in range(sf_direct.nrows):
        sf_grism.xshift[i*2] = sf_direct.xshift[i]
        sf_grism.xshift[i*2+1] = sf_direct.xshift[i]
        sf_grism.yshift[i*2] = sf_direct.yshift[i]
        sf_grism.yshift[i*2+1] = sf_direct.yshift[i]
    
    sf_grism.write(sf_grism.filename)
    
    ## Now do grism
    pair('TILE41-132-F160W_asn.fits', 'TILE41-132-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate,shift')
    
    #### Second G141 angle
    ## direct only
    pair('TILE41-152-direct_asn.fits', 'TILE41-152-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
    threedhst.shifts.make_grism_shiftfile('TILE41-152-direct_asn.fits', 'TILE41-152-G141_asn.fits')
    sf_direct = threedhst.shifts.ShiftFile('TILE41-152-direct_shifts.txt')
    sf_grism = threedhst.shifts.ShiftFile('TILE41-152-G141_shifts.txt')
    
    for i in range(sf_direct.nrows):
        sf_grism.xshift[i*2] = sf_direct.xshift[i]
        sf_grism.xshift[i*2+1] = sf_direct.xshift[i]
        sf_grism.yshift[i*2] = sf_direct.yshift[i]
        sf_grism.yshift[i*2+1] = sf_direct.yshift[i]
    
    sf_grism.write(sf_grism.filename)
    
    ## Now do grism
    pair('TILE41-152-direct_asn.fits', 'TILE41-152-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate,shift')

    #### G102 Visit
    ## direct only
    pair('TILE41-129-direct_asn.fits', 'TILE41-129-G102_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=True, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
    threedhst.shifts.make_grism_shiftfile('TILE41-129-direct_asn.fits', 'TILE41-129-G102_asn.fits')
    sf_direct = threedhst.shifts.ShiftFile('TILE41-129-direct_shifts.txt')
    sf_grism = threedhst.shifts.ShiftFile('TILE41-129-G102_shifts.txt')
    
    for i in range(sf_direct.nrows):
        sf_grism.xshift[i*2] = sf_direct.xshift[i]
        sf_grism.xshift[i*2+1] = sf_direct.xshift[i]
        sf_grism.yshift[i*2] = sf_direct.yshift[i]
        sf_grism.yshift[i*2+1] = sf_direct.yshift[i]
    
    sf_grism.write(sf_grism.filename)
    
    ###### Now do grism
    ### G102 Flat
    sky = pyfits.open('../CONF/WFC3.IR.G102.sky.V1.0.fits')
    IREF = os.getenv('iref')
    flat = pyfits.open(IREF+'/uc72113oi_pfl.fits')
    flat_im = flat[1].data[5:-5,5:-5]
    pyfits.writeto('../CONF/sky_g102_f105w.fits', data=sky[0].data/flat_im, header=sky[0].header, clobber=True)
    
    threedhst.grism_sky.flat_f140 = pyfits.open(IREF+'/uc72113oi_pfl.fits')
    threedhst.grism_sky.flat_g141 = pyfits.open(IREF+'/u4m1335li_pfl.fits')
    threedhst.grism_sky.flat = threedhst.grism_sky.flat_g141[1].data[5:1019,5:1019] / threedhst.grism_sky.flat_f140[1].data[5:1019, 5:1019]
    threedhst.grism_sky.flat[flat <= 0] = 5
    threedhst.grism_sky.flat[flat > 5] = 5
    
    pair('TILE41-129-direct_asn.fits', 'TILE41-129-G102_asn.fits', adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate,shift', sky_images=['sky_g102_f105w.fits'])
    
    #### Make correct direct images separated by band
    threedhst.utils.combine_asn_shifts(['TILE41-129-direct_asn.fits', 'TILE41-132-F160W_asn.fits','TILE41-152-direct_asn.fits'], out_root='TILE41-direct', path_to_FLT='./', run_multidrizzle=False)
    
    sf = threedhst.shifts.ShiftFile('TILE41-direct_shifts.txt')
    import copy
    
    for filter in ['F105W','F125W','F160W']:
        print unicorn.noNewLine+filter
        asn = threedhst.utils.ASNFile('TILE41-direct_asn.fits')
        sfout = copy.deepcopy(sf)
        for i in np.arange(sfout.nrows-1,-1,-1):
            head = pyfits.getheader(sfout.images[i])
            if head['FILTER'] != filter:
                p = asn.exposures.pop(i)
                p = sfout.pop(i)
        #
        asn.product = 'TILE41-%s' %(filter)
        asn.write('TILE41-%s_asn.fits' %(filter))
        sfout.write('TILE41-%s_shifts.txt' %(filter))
    
    #
    for filter in ['F105W','F125W','F160W']:
        threedhst.prep_flt_files.startMultidrizzle('TILE41-%s_asn.fits' %(filter), use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.6, driz_cr=False, updatewcs=False, clean=True, median=False, ra=150.0726, dec=2.1947412, final_outnx=2890, final_outny=2680)
    #
    from pyraf import iraf

    for filter in ['F105W','F125W','F160W']:
        try:
            os.remove('../DATA/TILE41-%s_sci.fits' %(filter))
        except:
            pass
        #
        iraf.imcopy('TILE41-%s_drz.fits[1]' %(filter), '../DATA/TILE41-%s_sci.fits' %(filter))
    
    ### 3-Color Display
    ds9 = threedhst.dq.myDS9()
    scale = 'linear'
    ds9.set('rgb yes')
    scl = (-0.02,1.)
    
    ds9.set('rgb blue')
    ds9.set('file TILE41-F105W_drz.fits')
    s = 10**(0.4*(26.27-25.96))/(1.6/1.0)**2
    ds9.set('scale limits %.3f %.3f' %(scl[0]*s, scl[1]*s))
    ds9.set('scale '+scale)
    
    ds9.set('rgb green')
    ds9.set('file TILE41-F125W_drz.fits')
    s = 10**(0.4*(26.25-25.96))/(1.6/1.25)**2
    ds9.set('scale limits %.3f %.3f' %(scl[0]*s, scl[1]*s))
    ds9.set('scale '+scale)

    ds9.set('rgb red')
    ds9.set('file TILE41-F160W_drz.fits')
    ds9.set('scale limits %.3f %.3f' %(scl[0], scl[1]))
    ds9.set('scale '+scale)

    ds9.set('rgb lock colorbar')
    
    ##########
    ###
#
def SN_CLASH():
    """
    Reductions for the Perlmutter MACS1720 grism observations
    """
    import glob
    import unicorn.candels
    import threedhst
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    
    os.chdir("/3DHST/Spectra/Work/SN-PERLMUTTER/PREP_FLT")
    unicorn.candels.make_asn_files()
    
    ALIGN_IMAGE = '/Volumes/WD3_Data/Users/gbrammer/CLASH/macs1720/hlsp_clash_hst_wfc3ir_macs1720_f160w_v1_drz.fits'
    
    ALIGN_IMAGE = '../CLASH/hlsp_clash_hst_wfc3ir_macs1720_f160w_v1_drz.fits'
    files=glob.glob('MACS*F1*asn.fits')
    files.extend(glob.glob('SN-*F1*asn.fits'))
        
    ## Subtract imaging backgrounds
    for file in files:
        if not os.path.exists(file.replace('asn','drz')):
            threedhst.shifts.make_blank_shiftfile(asn_file=file, xshift=0, yshift=0, rot=0., scale=1.0)
            #
            unicorn.candels.prep_candels(asn_file=file, ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0, GET_SHIFT=False, DIRECT_HIGHER_ORDER=2)
    
    ##  Align to clash
    for file in files:
        for geom in ['rotate', 'shift']:
            threedhst.shifts.refine_shifts(ROOT_DIRECT=file.split('_asn')[0], ALIGN_IMAGE=ALIGN_IMAGE, ALIGN_EXTENSION=0, fitgeometry=geom, clean=True)
            #
            threedhst.prep_flt_files.startMultidrizzle(file, use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
    
    ## Make combined reference images
    threedhst.utils.combine_asn_shifts(glob.glob('[SM]*-3*F105W_asn.fits'), 'MACSJ1720+3536-F105W')
    threedhst.prep_flt_files.startMultidrizzle('MACSJ1720+3536-F105W_asn.fits', use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=1, driz_cr=False, updatewcs=False, clean=True, median=False)

    threedhst.utils.combine_asn_shifts(glob.glob('[SM]*-3*F1[46]0W_asn.fits'), 'MACSJ1720+3536-F140W')
    threedhst.prep_flt_files.startMultidrizzle('MACSJ1720+3536-F140W_asn.fits', use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=1, driz_cr=False, updatewcs=False, clean=True, median=False)
    
    ### Subtract grism background (SKIP_DIRECT=True)
    threedhst.prep_flt_files.process_3dhst_pair('SN-L1-PANCHA-307-F160W_asn.fits', 'SN-L1-PANCHA-307-G141_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, adjust_targname=False)
    threedhst.prep_flt_files.process_3dhst_pair('SN-L1-PANCHA-312-F160W_asn.fits', 'SN-L1-PANCHA-312-G141_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, adjust_targname=False)
    
    threedhst.prep_flt_files.process_3dhst_pair('SN-L1-PANCHA-307-F105W_asn.fits', 'SN-L1-PANCHA-307-G102_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, adjust_targname=False, sky_images=['../CONF/G102_master_flatcorr.fits'])
    threedhst.prep_flt_files.process_3dhst_pair('SN-L1-PANCHA-312-F105W_asn.fits', 'SN-L1-PANCHA-312-G102_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, adjust_targname=False, sky_images=['../CONF/G102_master_flatcorr.fits'])
    
    #### Make catalog and segmentation image from F140W
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = 'MACSJ1720+3536-F140W.cat'
    se.options['CHECKIMAGE_NAME'] = 'MACSJ1720+3536-F140W_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'MACSJ1720+3536-F140W_drz.fits[1]'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_MINAREA']    = '10'
    se.options['DETECT_THRESH']    = '1.5'
    se.options['ANALYSIS_THRESH']  = '1.5'
    se.options['MAG_ZEROPOINT'] = '%.5f' %(unicorn.reduce.ZPs['F140W'])
    se.options['DEBLEND_NTHRESH'] = '64'
    se.options['DEBLEND_MINCONT'] = '0.00002'
    se.options['PHOT_APERTURES'] = '6.6667,8.333333,16.667,8.335,20'
    se.options['GAIN']  = '1.0'
    
    status = se.sextractImage('MACSJ1720+3536-F140W_drz.fits[0]')
    
    threedhst.sex.sexcatRegions('MACSJ1720+3536-F140W.cat', 'MACSJ1720+3536-F140W.reg', format=2)
    
    ### xxx regenerate drz images
    files=glob.glob('SN*F1*asn.fits')
    for file in files:
        threedhst.prep_flt_files.startMultidrizzle(file, use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False)
        
    #### Need to make fake grism reference images with four exposures so that blot_combine will work
    #### Then need to run multidrizzle.  The results will be nonsense because direct/grism mixed
    switch = {'F160W':'G141', 'F105W':'G102'}
    for filter in switch.keys():
        files=glob.glob('SN-L1*%s*asn.fits' %(filter))
        for file in files:
            #### Put last three grism exposures in the direct ASN
            asn_grism = threedhst.utils.ASNFile(file.replace(filter, switch[filter]))
            asn_direct = threedhst.utils.ASNFile(file)
            for i in range(1,4):
                asn_direct.exposures.append(asn_grism.exposures[i])
            #
            asn_direct.product = 'X'+asn_direct.product
            asn_direct.write('X'+asn_direct.file)
            ### Make shiftfile for the new list, overwriting the old
            threedhst.shifts.make_grism_shiftfile(asn_grism.file, 'X'+asn_direct.file)
            threedhst.prep_flt_files.startMultidrizzle('X'+file, use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False, refimage='MACSJ1720+3536-F140W_drz.fits[1]')
    

    #### Now set up blotting back
    unicorn.reduce.prepare_blot_reference(REF_ROOT='M1720_F160W', filter='F140W', REFERENCE = 'MACSJ1720+3536-F140W_drz.fits', SEGM = 'MACSJ1720+3536-F140W_seg.fits', sci_extension=1)
    unicorn.reduce.prepare_blot_reference(REF_ROOT='M1720_F105W', filter='F105W', REFERENCE = 'MACSJ1720+3536-F105W_drz.fits', SEGM = 'MACSJ1720+3536-F140W_seg.fits', sci_extension=1)
    
    NGROW=0 #125
    pad=60
    CATALOG = 'MACSJ1720+3536-F140W.cat'
    
    for filter in ['F105W', 'F160W']:
        files=glob.glob('XSN-L1-PANCHA*%s_asn.fits' %(filter))
        REF_ROOT = 'M1720_%s' %(filter)
        for file in files:
            unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT=file.split('_asn')[0], NGROW=NGROW, verbose=True)
    
    #### above not working
    #### 
    for filter in ['F105W']:
        files=glob.glob('XSN-L1-PANCHA*%s_asn.fits' %(filter))
        REF_ROOT = 'M1720_%s' %(filter)
        for file in files:
            unicorn.reduce.interlace_combine_blot(root=file.split('_asn')[0], view=False, pad=pad, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, growx=2, growy=2, auto_offsets=True)
            #
            unicorn.reduce.interlace_combine(root=file.split('_asn')[0][1:], view=False, pad=pad,  NGROW=NGROW, growx=2, growy=2, auto_offsets=True)
            unicorn.reduce.interlace_combine(root=file.split('_asn')[0][1:].replace(filter,'G141'), view=False, pad=pad,  NGROW=NGROW, growx=2, growy=2, auto_offsets=True)
    
    #### have to rename XSN*_ref_inter.fits to SN*_ref_inter.fits 
    #
    model = unicorn.reduce.GrismModel(root='SN-L1-PANCHA-307', grow_factor=2, growx=2, growy=2, MAG_LIMIT=22, use_segm=False, grism='G102', direct='F105W')
    
    model = unicorn.reduce.GrismModel(root='SN-L1-PANCHA-307', grow_factor=2, growx=2, growy=2, MAG_LIMIT=22, use_segm=False, grism='G141', direct='F160W')
    model.compute_full_model(BEAMS=['A', 'B', 'C', 'D'], view=None, MAG_LIMIT=22.0, save_pickle=True, refine=True, model_slope=0)
    

    model2 = unicorn.reduce.GrismModel(root='SN-L1-PANCHA-312', grow_factor=2, growx=2, growy=2, MAG_LIMIT=22, use_segm=False, grism='G141', direct='F160W')
    model2.compute_full_model(BEAMS=['A', 'B', 'C', 'D'], view=None, MAG_LIMIT=22.0, save_pickle=True, refine=True, model_slope=0)
    
    id=622
    for m in [model, model2]:
        m.twod_spectrum(id, miny=-60, USE_REFERENCE_THUMB=True)
        m.show_2d(savePNG=True)
        
    #
    unicorn.hudf.stack(id=id, dy=20, save=True, inverse=False, scale=[1,99], fcontam=100., ref_wave = 1.4e4, root='MACSJ1720', search='SN-L1-PANCHA')
    
    #### Set up redshift fits.  Will read photometric/eazy catalogs defined in unicorn.analysis.read_catalogs
    gris = unicorn.interlace_fit.GrismSpectrumFit('MACSJ1720_%05d' %(id))
    gris = unicorn.interlace_fit.GrismSpectrumFit('SN-L1-PANCHA-307_%05d' %(id))
    gris = unicorn.interlace_fit.GrismSpectrumFit('SN-L1-PANCHA-312_%05d' %(id))
    #gris = unicorn.interlace_fit.GrismSpectrumFit('MACSJ1720_%05d' %(id))
    
    ### Fit it!
    gris.fit_in_steps(zrfirst=[0.4,2.4])
    
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

    # reload threedhst.grism_sky   #### Reset the flat-field for F140W/G141
    
    os.chdir(unicorn.GRISM_HOME+'ERS/PREP_FLT')
    #ALIGN = '/3DHST/Ancillary/GOODS-S/GOODS_ACS/h_sz*drz_img.fits'
    ALIGN_FILES=['/3DHST/Ancillary/GOODS-S/CANDELS/ucsc_mosaics/GOODS-S_F160W_wfc3ir_drz_sci.fits']
    
    ALIGN_FILES=['goodss_3dhst.v2.1.F140W_orig_sci.fits']
    direct = [1,'WFC3-ERSII-G01-F140W_asn.fits']
    grism = [1,'WFC3-ERSII-G01-G141_asn.fits']
    
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
    
    ##### Try interlacing it
    unicorn.reduce.prepare_blot_reference(REF_ROOT='ERS_F140W', filter='F140W', REFERENCE = 'goodss_3dhst.v2.1.F140W_orig_sci.fits', SEGM = 'goodss_3dhst.v2.1.F125W_F140W_F160W_seg.fits')

    NGROW=200
    pad=350
    CATALOG='goodss_3dhst.v2.1.F125W_F140W_F160W_det.cat'
    
    direct='WFC3-ERSII-G01-F140W_asn.fits'

    extract_limit = 24
    skip_completed=False
    REF_ROOT='ERS_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    pointing=threedhst.prep_flt_files.make_targname_asn(direct, newfile=False).split('-F140')[0]
    #
    unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)

    unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=pad, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, growx=1, growy=1, auto_offsets=True)
    unicorn.reduce.interlace_combine(pointing+'-F140W', pad=pad, NGROW=NGROW, growx=1, growy=1, auto_offsets=True)
    unicorn.reduce.interlace_combine(pointing+'-G141', pad=pad, NGROW=NGROW, growx=1, growy=1, auto_offsets=True)
    
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

#
def COOPER(FORCE=False):
    """
    Reduce pointings from a grism program in AEGIS by M. Cooper
    """
    import unicorn
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os

    os.chdir('/3DHST/Spectra/Work/COOPER/PREP_FLT')
    
    ### reduced by hand before???
    
    ## RGB Image
    imr = pyfits.open('EGS12007881.12012083-F140W_drz.fits')
    img = pyfits.open('EGS12007881.12012083-ACSi.fits')
    imb = pyfits.open('EGS12007881.12012083-ACSv.fits')
    
    f = 0.5
    scales = np.array([10**(-0.4*(26.46-26.46)), 10**(-0.4*(25.95928-26.46))*1.5, 10**(-0.4*(26.50512-26.46))*1.5])*f
    #
    xc, yc, NX = 1754, 1924, 50
    sub_r = imr[1].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[0]
    sub_g = img[0].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[1]*1.0
    sub_b = imb[0].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[2]*1.0
    #
    # sub_r = imr[1].data*scales[0]
    # sub_g = img[0].data*scales[1]*1.0
    # sub_b = imb[0].data*scales[2]*1.0
    #
    Q, alpha, m0 = 3., 8, -0.02
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='COOPER_viH.png', shape=sub_r.shape)
    
    #### ACS + F140W photometry
    import pyregion
    filters = ['F140W', 'F814W', 'F606W']
    ZPs = [26.46, 25.95928, 26.50512]
    ## ACS http://www-int.stsci.edu/~koekemoe/egs/
    ZPs = [26.46, 25.937, 26.486] 
    images = [imr[1], img[0], imb[0]]
    
    fluxes = {}
    mags = {}
    for id in [0, 616, 625]:
        fluxes[id] = {}
        mags[id] = {}
        
    for filter, ZP, image in zip(filters, ZPs, images):
        print filter
        for id in [0, 616, 625]:
            r = pyregion.open('blob_%d.reg' %(id)).as_imagecoord(header=image.header)
            aper = r.get_mask(hdu=image)
            fluxes[id][filter] = np.sum(image.data*aper)
            mags[id][filter] = -2.5*np.log10(fluxes[id][filter])+ZP
    
    #### CFHT-LS photometry, color image
    ext = '14561_14592_7791_7822.fits' # 6"
    #ext = '14405_14748_7635_7978.fits' # larger
    
    imr = pyfits.open('CFHTLS/D3.IQ.I.%s' %(ext))
    img = pyfits.open('CFHTLS/D3.IQ.R.%s' %(ext))
    imb = pyfits.open('CFHTLS/D3.IQ.G.%s' %(ext))
    
    #f = 10
    f = 0.4
    scales = np.array([10**(-0.4*(30-26.46))*1, 10**(-0.4*(30-26.46))*1.2, 10**(-0.4*(30-26.46))*1.5])*f
    sub_r = imr[0].data*scales[0]
    sub_g = img[0].data*scales[1]*1.0
    sub_b = imb[0].data*scales[2]*1.0
    
    Q, alpha, m0 = 3., 8, -0.02
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='COOPER_CFHTLS_gri.png', shape=(100,100))
    
    for filter in ['U','G','R','I','Z']:
       image = pyfits.open('CFHTLS/D3.IQ.%s.%s' %(filter, ext))[0]
       for id in [0, 616, 625]:
           r = pyregion.open('blob_%d.reg' %(id)).as_imagecoord(header=image.header)
           aper = r.get_mask(hdu=image)
           fluxes[id][filter] = np.sum(image.data*aper)
           mags[id][filter] = -2.5*np.log10(fluxes[id][filter])+30
       
    import cPickle as pickle
    fp = open('egs_blob.pkl','wb')
    pickle.dump(filters, fp)
    pickle.dump(ZPs, fp)
    pickle.dump(fluxes, fp)
    pickle.dump(mags, fp)
    fp.close()
    
    fp = open('egs_blob.readme','w')
    fp.write('filters ZPs fluxes mags, 616 is the blob, 625 is the companion\n')
    fp.close()
    
def make_new_flats():
    """
    Go and find all 3D-HST flt/seg/mask.reg images
    and copy them to a new directory.  Then generate a master flat-field.
    """
    import shutil
    
    PATH = '/3DHST/Spectra/Work/FLATS/F140W'
    for field in ['AEGIS','COSMOS','GOODS-S','GOODS-N','UDS']:
        os.chdir('/3DHST/Spectra/Work/%s/PREP_FLT/' %(field))
        info = catIO.Readfile('files.info')
        f140 = info.filter == 'F140W'
        print '%s: %d (%d)' %(field, f140.sum(), f140.sum()/4.)
        for file in info.file[f140]:
            root=file.split('.fits')[0]
            seg = glob.glob(root+'.seg.fits')
            if len(seg) > 0:
                print '%s, %s' %(field, root)
                shutil.copy('../RAW/%s' %(file), PATH)
                shutil.copy(seg[0], PATH)
                if os.path.exists('%s.fits.mask.reg' %(root)):
                    shutil.copy('%s.fits.mask.reg' %(root), PATH)
    
    #### Make flats for different fields
    os.chdir(PATH)
    info = catIO.Readfile('files.info')
    info.field_name = []
    for targ in info.targname:
        info.field_name.append(targ[:3])
    
    info.field_name = np.array(info.field_name)
    
    cosmos = info.field_name == 'COS'
    flt_files = info.file[cosmos]
    unicorn.prepare.make_flat(flt_files, output='cosmos_v0.fits', GZ='')

    subset = info.field_name == 'AEG'
    flt_files = info.file[subset]
    unicorn.prepare.make_flat(flt_files, output='aegis_v0.fits', GZ='')

    subset = info.field_name == 'UDS'
    flt_files = info.file[subset]
    unicorn.prepare.make_flat(flt_files, output='uds_v0.fits', GZ='')

    subset = info.field_name == 'GOO'
    flt_files = info.file[subset]
    unicorn.prepare.make_flat(flt_files, output='goodss_v0.fits', GZ='')

    subset = info.field_name == 'GNG'
    flt_files = info.file[subset]
    unicorn.prepare.make_flat(flt_files, output='goodsn_v0.fits', GZ='')
    
    subset = (info.field_name == 'CAM') | (info.field_name == 'COL') | (info.field_name == 'GRB') | (info.field_name == 'HIG') | (info.field_name == 'HUD') 
    flt_files = info.file[subset]
    unicorn.prepare.make_flat(flt_files, output='flat_UDF_F140W_v0.fits', GZ='')
    
    ### By time
    import astropy.time
    t = astropy.time.Time(info.date_obs, format='iso', scale='utc')
    so = np.argsort(t.mjd)
    
    so = so[info.field_name != 'GNG']
    N = 3
    NT = len(so)/N
    for i in range(N):
        subset = so[(i*NT):(i+1)*NT]
        flt_files = info.file[subset]
        unicorn.prepare.make_flat(flt_files, output='flat_3DHST_F140W_t%d_v0.1.fits' %(i+1), GZ='')
    #
    for i in range(N):
        subset = so[(i*NT):(i+1)*NT]
        flt_files = info.file[subset]
        print 'flat_3DHST_F140W_t%d_v0.1.fits  %.1f' %(i+1, t.mjd[subset][0])
    
    # flat_3DHST_F140W_t1_v0.1.fits  55465.0
    # flat_3DHST_F140W_t2_v0.1.fits  55793.0
    # flat_3DHST_F140W_t3_v0.1.fits  55911.0
        
    # unicorn.prepare.make_flat(flt_files, output='xxx.fits', GZ='')
    
    ##### CANDELS on UNICORN
    shell = """
    cd /Users/gbrammer/FLATS/CANDELS
    for field in AEGIS COSMOS GOODS-S GOODS-N UDS; do
        ln -sf /3DHST/Ancillary/${field}/CANDELS/PREP_FLT/*seg.fits .
        ln -sf /3DHST/Ancillary/${field}/CANDELS/PREP_FLT/*mask.reg .
        ln -sf /3DHST/Ancillary/${field}/CANDELS/RAW/*flt.fits.gz .
    done
    
    cd /Users/gbrammer/FLATS/3DHST
    for field in AEGIS COSMOS GOODS-S GOODS-N UDS; do
        ln -sf /3DHST/Spectra/Work/${field}/PREP_FLT/*seg.fits .
        ln -sf /3DHST/Spectra/Work/${field}/PREP_FLT/*mask.reg .
        ln -sf /3DHST/Spectra/Work/${field}/RAW/*flt.fits.gz .
    done
    
    """
    
    os.chdir('/Users/gbrammer/FLATS/CANDELS')
    info = catIO.Readfile('files.info')
    info.field_name = []
    for targ, ra in zip(info.targname, info.ra_targ):
        info.field_name.append(targ[:3]+'%03d' %(int(ra/10)*10))
    
    info.field_name = np.array(info.field_name)
    
    t = astropy.time.Time(info.date_obs, format='iso', scale='utc')
    so = np.argsort(t.mjd)
    
    for filter in ['F125W', 'F160W']:
        ok = info.filter[so] == filter
        N = 2
        NT = ok.sum()/N
        for i in range(N):
            subset = so[ok][(i*NT):(i+1)*NT]
            flt_files = info.file[subset]
            unicorn.prepare.make_flat(flt_files, output='flat_%s_t%d_v0.3.fits' %(filter, i+1), GZ='')
    #
    for filter in ['F125W', 'F160W']:
        ok = info.filter[so] == filter
        N = 2
        NT = ok.sum()/N
        for i in range(N):
            subset = so[ok][(i*NT):(i+1)*NT]
            flt_files = info.file[subset]
            print 'flat_%s_t%d_v0.3.fits  %.1f' %(filter, i+1, t.mjd[subset][0])
     
    # flat_F125W_t1_v0.1.fits  55071.0
    # flat_F125W_t2_v0.1.fits  55770.0
    # flat_F160W_t1_v0.1.fits  55079.0
    # flat_F160W_t2_v0.1.fits  55717.0
     
    
def clean_flat():
    """
    """

def make_flat(flt_files, output='master_flat.fits', GZ=''):
    """
    Given a list of FLT files, make a median average master flat
    
    GZ = '.gz' if segmentation images are gzipped
    
    """    
    import iraf

    import scipy.ndimage as nd
    from iraf import noao
    from iraf import imred
    from iraf import ccdred
    
    pipeline_flats = {}
    
    NF = len(flt_files)
    idx = np.arange(NF)
    X = np.zeros((NF,1014.**2))
    
    root = output.split('.fits')[0]
    
    fp_ic = open('%s.imcombine' %(root), 'w')
    fp_info = open('%s.image_info' %(root), 'w')
    
    for i, file in enumerate(flt_files):
        fi = file.replace('.seg','')
        segfile = fi.split('flt')[0]+'flt.seg.fits%s' %(GZ)
        if not os.path.exists(segfile):
            continue
        #    
        ### Skip masked images that might have problems
        if os.path.exists(fi.replace('.gz','')+'.mask.reg'):
            print 'Mask found, skip...'
            continue
        #
        print '%d %s' %(i, file)
        flt = pyfits.open(fi)
        #### Divide by the pipeline flat
        pfl = flt[0].header['PFLTFILE'].split('iref$')[1]
        if pfl not in pipeline_flats.keys():
            print 'Pipeline flat: %s' %(pfl)
            flat = pyfits.open(os.getenv('iref') + '/%s' %(pfl))[1].data[5:-5,5:-5]
            flat[flat <= 0] = 5
            flat[flat > 5] = 5
            pipeline_flats[pfl] = flat*1.
        else:
            flat = pipeline_flats[pfl]
        #
        ### Segmentation mask
        segm = pyfits.open(segfile)[0].data
        grow = nd.maximum_filter(segm, size=(5,5))
        masked = grow == 0
        ### DQ mask, hot pixels and the "death star"
        dq_ok = (flt[3].data & (4+32+16)) == 0
        #dq_ok = (flt[3].data & (4+32)) == 0
        #
        ok = masked & np.isfinite(flt[1].data) & (dq_ok)
        background =  np.median(flt[1].data[ok])
        flt[1].data *= flat/background
        #flt[1].data /= background
        flt[1].data[(ok == False)] = 0
        X[i,:] = flt[1].data.flatten()
        ### For imcombine
        if not os.path.exists('imcombine'):
            os.mkdir('imcombine')
        #
        flt[1].data[(ok == False)] = -9999
        icfile = 'imcombine/%s.fits' %(fi.split('_flt')[0])
        flt[1].header.update('SKYMED', background)
        flt[1].header.update('SKYCTS', background*flt[0].header['EXPTIME'])
        pyfits.writeto(icfile, data=flt[1].data, header=flt[1].header, clobber=True)
        fp_ic.write(icfile+'\n')
        fp_info.write('%s %f %f\n' %(fi, flt[0].header['EXPTIME'], background))
        
    fp_ic.close()
    fp_info.close()
    
    #### Average
    nsum = np.sum(X != 0, axis=0).reshape(1014,1014)
    avg = np.sum(X, axis=0).reshape(1014,1014)/nsum
     
    ### Fill empty pixels with no input images
    sky = avg
    
    #### Imcombine, much faster than numpy for median 
    files=glob.glob('ic_%s*' %(root))
    for file in files: 
        os.remove(file)
    
    NREJECT = np.maximum(np.minimum(5, (NF/2)-5), 1)
    iraf.imcombine ( input = '@%s.imcombine' %(root), output = 'ic_%s.fits' %(root), 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = 'ic_%s_n.pl' %(root), 
       expmasks = '', sigmas = 'ic_%s_s.fits' %(root), 
       logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'double', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = '!SKYCTS', statsec = '', expname = '', lthreshold = 0.5, 
       hthreshold = 2, nlow = NREJECT, nhigh = NREJECT, nkeep = 5, 
       mclip = iraf.yes, lsigma = 4.0, hsigma = 4.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
    
    #iraf.badpiximage('ic_%s_n.pl' %(root), 'ic_%s.fits' %(root), 'ic_'
    iraf.imcopy('ic_%s_n.pl' %(root), 'ic_%s_n.fits' %(root))
    
    ic = pyfits.open('ic_%s.fits' %(root))
    
    ### Checking stats
    # X2 = X.reshape((len(flt_files), 1014, 1014))
    # so = np.argsort(X2, axis=0)
    # i, j = 500, 500
    # for i in range(495,505):
    #     for j in range(495, 505):
    #         column = X2[:,i,j][so[:,i,j]]
    #         ok = (column > 0.5) & (column < 2)
    #         data = np.mean(column[ok][3:-3])
    #         print '%d %d %f %f %f' %(i,j, data, ic[0].data[i,j], data/ic[0].data[i,j])
            
    sky = ic[0].data
    sky[sky == 0] = np.inf
    #flatim = pyfits.open(os.getenv('iref')+'/uc721143i_pfl.fits')
    #flatim[1].data[5:-5,5:-5] = ic[0].data
    #flatim.writeto('ic2.fits', clobber=True)
    
    ### Median with mask
    # m = np.ma.masked_array(X, X == 0)
    # sky = np.ma.median(m, axis=0)
    # rms = np.ma.std(m)
    # sky = sky.reshape(1014,1014)
    
    sigma = pyfits.open('ic_%s_s.fits' %(root))[0].data
    med = np.median(sigma[sigma != 0])
    std = threedhst.utils.biweight(sigma[sigma != 0])
    noisy = (sigma-med) > (5*std)
    
    pam = pyfits.open('%s/ir_wfc3_map.fits' %(os.getenv('iref')))[1].data
    
    bad = (np.isfinite(sky) == False) | (sky/flat > 1.15) | noisy | (sky/pam > 1.08)
    x,y = np.where(bad)
    NX = len(x)
    pad = 1
    for i in range(NX):
        xi = x[i]
        yi = y[i]
        sub = sky[xi-pad:xi+pad+2,yi-pad:yi+pad+2]
        if (np.sum(sub) != 0.0):
            sky[xi,yi] = np.median(sub[np.isfinite(sub)])
    
    still_bad = (np.isfinite(sky) == False) | (sky <= 0.01)
    sky[still_bad] = flat[still_bad]
    
    #### Some residual scaling W.R.T pipeline flats??
    ratio = sky / pipeline_flats[pipeline_flats.keys()[0]]
    fig = unicorn.plotting.plot_init(xs=4, square=True, NO_GUI=True)
    ax = fig.add_subplot(111)
    med, std = np.median(ratio[~bad]), threedhst.utils.biweight(ratio[~bad])
    ax.hist(ratio[~bad].flatten(), range=(0.9,1.1), bins=100, alpha=0.5, label='%s: %.4f, %.4f' %(pipeline_flats.keys()[0], med, std))
    ax.legend()
    
    sky /= med
    
    unicorn.plotting.savefig(fig, '%s_scale.png' %(root))
    
    
    # bad_flat = (flat < 0.5)
    # sky[bad_flat] = flat[bad_flat]
        
    # im_sky = pyfits.PrimaryHDU(data=sky)
    # im_n = pyfits.ImageHDU(data=nsum)
    # im = pyfits.HDUList([im_sky, im_n])
    # im.writeto('sky.fits', clobber=True)
    
    #### for DIRECT flat
    flatim = pyfits.open(os.getenv('iref')+'/uc721143i_pfl.fits')
    flatim[1].data[5:-5,5:-5] = sky
    flatim[3].data[5:-5,5:-5] = nsum
    flatim.writeto(output, clobber=True)
    
    #flatim.writeto('/research/HST/GRISM/IREF/cosmos_f140w_flat.fits', clobber=True)
               
def redo_all_sky_subtraction():
    import unicorn.prepare
    
    for field in ['AEGIS','COSMOS','UDS','GOODS-S','GOODS-N']:
        os.chdir(unicorn.GRISM_HOME+'/%s/PREP_FLT' %(field))
        asn_files = glob.glob(field+'-[0-9]*G141_asn.fits')
        for asn_file in asn_files:
            print asn_file
            unicorn.prepare.redo_sky_subtraction(asn_file=asn_file, PATH_TO_RAW='../RAW', sky_images=['sky.G141.set001.fits', 'sky.G141.set002.fits', 'sky.G141.set003.fits', 'sky.G141.set004.fits', 'sky.G141.set005.fits', 'sky.G141.set025.fits', 'sky.G141.set120.fits'])
            
def redo_sky_subtraction(asn_file='GOODS-S-30-G141_asn.fits', PATH_TO_RAW='../RAW', sky_images=['sky.G141.set001.fits', 'sky.G141.set002.fits', 'sky.G141.set003.fits', 'sky.G141.set004.fits', 'sky.G141.set005.fits', 'sky.G141.set025.fits', 'sky.G141.set120.fits'], sky_components=False):
    """
    
    sky_images = {'G141':['G141_zodi_0.00_excess_1.00.fits', 'G141_zodi_0.25_excess_0.75.fits', 'G141_zodi_0.50_excess_0.50.fits', 'G141_zodi_0.75_excess_0.25.fits', 'G141_zodi_1.00_excess_0.00.fits'], 
                  'G102':['G102_zodi_0.00_excess_1.00.fits', 'G102_zodi_0.25_excess_0.75.fits', 'G102_zodi_0.50_excess_0.50.fits', 'G102_zodi_0.75_excess_0.25.fits', 'G102_zodi_1.00_excess_0.00.fits']}
    
    
    """
    
    import threedhst
    import threedhst.grism_sky
    
    asn = threedhst.utils.ASNFile(asn_file)
    
    threedhst.process_grism.fresh_flt_files(asn_file, 
                 from_path=PATH_TO_RAW, preserve_dq=False)
    
    ## Run the sky background division             
    asn_grism = threedhst.utils.ASNFile(asn_file)
    for exp in asn_grism.exposures:
        threedhst.grism_sky.remove_grism_sky(flt=exp+'_flt.fits', list=sky_images, path_to_sky='../CONF/', verbose=True, second_pass=True, overall=True, sky_components=sky_components)
    
    ## Run Multidrizzle twice, the first time to flag CRs + hot pixels
    threedhst.prep_flt_files.startMultidrizzle(asn_file, use_shiftfile=True, skysub=False,
            final_scale=0.128, pixfrac=1.0, driz_cr=True,
            updatewcs=True, median=True, clean=True)
            
    threedhst.prep_flt_files.startMultidrizzle(asn_file, use_shiftfile=True, skysub=False,
            final_scale=0.06, pixfrac=0.8, driz_cr=False,
            updatewcs=False, median=False, clean=True)
    
    
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

def make_mask_crr_file(bad_value=1024):
    """
    Need to edit CRR file in $iref to not make flagged reads UNSTABLE
    """
    crr = pyfits.open(os.getenv('iref')+'/u6a1748ri_crr.fits')
    try:
        crr[1].data['BADINPDQ'] *= 0
    except:
        pass
    #
    try:
        crr[1].data['BADINPDQ'] += bad_value
    except:
        pass

    crr.writeto('test_crr.fits', clobber=True)
    
def flag_bad_reads(image='ib3701s4q_raw.fits', flag_reads = None, ds9=None, vmi=0, vma=2, bad_value=1024, satellite_reads=[]):
    """
    Flag individual reads of an RAW image and re-run calwf3
    """
    import mywfc3.bg
    
    raw = pyfits.open(image, mode='update')
    time, ramp, reads = mywfc3.bg.get_bg_ramp(raw)
    raw = pyfits.open(image, mode='update')
    dt = np.diff(time)
    NSAMP = raw[0].header['NSAMP']*1
    ### Display
    if flag_reads is None:
        if ds9 is None:
            ds9 = threedhst.dq.myDS9()
        
        ds9.set('width 867')
        ds9.set('height 756')
        #
        ds9.set('tile yes')
        ds9.set('frame delete all')
        for i in range(NSAMP-1):
            ds9.frame(i+1)
            ds9.v(reads[i,:,:]/dt[i], vmin=vmi, vmax=vma)
        #
        ds9.set('zoom to fit')
        ds9.set('match frames image')
        ds9.set('lock colorbar')
        ds9.set('lock scale')
        #
        read_frames = raw_input('Frames to ignore: ')
        if (read_frames == ' '):
            if len(satellite_reads) == 0:
                return False
            else:
                flag_reads = []
        else:
            flag_reads = np.cast[int](read_frames.split(','))
            fp = open(image.replace('_raw.fits', '_raw.flag_reads'), 'w')
            fp.write(' '.join(np.cast[str](flag_reads))+'\n')
            fp.close()
        
    for read in range(NSAMP):
        raw['DQ', NSAMP-read].header['PIXVALUE'] = 0
    
    # dr = np.diff(flag_reads)
    # i=0
    # if flag_reads[0] <=2 :
    #     raw['DQ', NSAMP-flag_reads[0]].header['PIXVALUE'] = bad_value
    #     for i in range(len(dr)):
    #         if dr[i] > 1:
    #             break
    #         else:
    #             raw['DQ', NSAMP-flag_reads[i+1]].header['PIXVALUE'] = bad_value
    
    for read in flag_reads:
        raw['DQ', NSAMP-read].header['PIXVALUE'] = bad_value
    
    for read in satellite_reads:
        raw['DQ', NSAMP-read].header['PIXVALUE'] = 0
        raw['SCI', NSAMP-read].data *= 10
        
    raw[0].header['BADINPDQ'] = bad_value # 0    
    raw[0].header['CRREJTAB'] = 'test_crr.fits' #'iref$u6a1748ri_crr.fits'
    raw.flush()
    
    from wfc3tools import calwf3
    for ext in ['_flt', '_ima']:
        if os.path.exists(image.replace('_raw', ext)):
            os.remove(image.replace('_raw', ext))
    
    calwf3.calwf3(image, verbose=True)
    
    print '\n*** Take out %d from FLT/DQ array ***' %(bad_value)
    flt = pyfits.open(image.replace('_raw', '_flt'), mode='update')
    flt['DQ'].data -= bad_value
    flt.flush()
    
def flag_all():
    """
    Flag high-background READS in the G141 images
    """
    ds9 = threedhst.dq.myDS9()
    unicorn.prepare.flag_bad_reads('ib3701s4q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3701skq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3702u8q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3702uoq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3703uzq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3703vfq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3703vmq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3704wrq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3704x8q_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3705y1q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3705ylq_raw.fits', ds9=ds9)
    
    unicorn.prepare.flag_bad_reads('ib3706b2q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3706biq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3706bpq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3707caq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3707cqq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3707cuq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3708i5q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3708ipq_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3708i9q_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3709j3q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3709joq_raw.fits', ds9=ds9)
    
    unicorn.prepare.flag_bad_reads('ib3710n6q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3710nmq_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3711bkq_raw.fits', ds9=ds9)
    
    unicorn.prepare.flag_bad_reads('ib3747a5q_raw.fits', ds9=ds9)

    unicorn.prepare.flag_bad_reads('ib3749o5q_raw.fits', ds9=ds9)
    unicorn.prepare.flag_bad_reads('ib3749oqq_raw.fits', ds9=ds9)
    
    ##### GOODS-N F140W
    unicorn.prepare.flag_bad_reads('ib3707c8q_raw.fits', ds9=ds9)
    
    
    #### Testing
    os.chdir('/Users/brammer/3DHST/Spectra/Work/GOODS-N/FixMultiAccum/TEST')
    
    image = 'ib3703vmq_raw.fits'
    bad = 1024
    
    raw = pyfits.open('../%s' %(image))
    raw.writeto(image, clobber=True)
    
    unicorn.prepare.make_mask_crr_file(bad_value=bad)
    unicorn.prepare.flag_bad_reads(image, ds9=ds9, bad_value=bad)
    
def test_masked():
    """
    Something like wf3ir with masked arrays
    """
    import numpy as np
    import numpy.ma as ma
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/GOODS-N/FixMultiAccum/TEST')
    image = 'ib3749o5q_raw.fits'
    
    raw = pyfits.open('../%s' %(image))
    raw.writeto(image, clobber=True)
    
    raw = pyfits.open(image, mode='update')
    time, ramp, reads = mywfc3.bg.get_bg_ramp(raw)
    dt = np.diff(time)
    NSAMP = raw[0].header['NSAMP']
    dt_2D = reads*0+1.
    for i in range(NSAMP-1):
        dt_2D[i,:,:] *= dt[i]
        
    cps = ma.MaskedArray(reads/dt_2D, mask=~np.isfinite(reads), copy=True)
    counts = ma.MaskedArray(reads, mask=~np.isfinite(reads), copy=True)
    GAIN = 2.4
    cps_err = np.sqrt(12**2+np.maximum(reads, 0)*GAIN)/dt_2D
    
    avg = counts.sum(axis=0)/dt_2D.sum(axis=0)
    
    rma.mask = np.abs(rma/std) > 4
    
#
def split_multiaccum(ima, scale_flat=True):
    """
    Pull out the MultiAccum reads of a RAW or IMA file into a single 3D 
    matrix.
    
    Returns cube[NSAMP,1024,1014], time, NSAMP
    """    
    skip_ima = ('ima' in ima.filename()) & (ima[0].header['FLATCORR'] == 'COMPLETE')
    if scale_flat & ~skip_ima:
        FLAT_F140W = pyfits.open(os.path.join(os.getenv('iref'), 'uc721143i_pfl.fits'))[1].data
    else:
        FLAT_F140W = 1
            
    NSAMP = ima[0].header['NSAMP']
    sh = ima['SCI',1].shape
    
    cube = np.zeros((NSAMP, sh[0], sh[1]))
    if 'ima' in ima.filename():
        dq = np.zeros((NSAMP, sh[0], sh[1]), dtype=np.int)
    else:
        dq = 0
        
    time = np.zeros(NSAMP)
    for i in range(NSAMP):
        if ima[0].header['UNITCORR'] == 'COMPLETE':
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data*ima['TIME',i+1].header['PIXVALUE']/FLAT_F140W
        else:
            cube[NSAMP-1-i, :, :] = ima['SCI',i+1].data/FLAT_F140W
        
        if 'ima' in ima.filename():
            dq[NSAMP-1-i, :, :] = ima['DQ',i+1].data
        
        time[NSAMP-1-i] = ima['TIME',i+1].header['PIXVALUE']
    
    return cube, dq, time, NSAMP
    
# def flag_RAW_crs(raw='icbw16e7q_raw.fits'):
#     """
#     Manually subtract variable average background from RAW reads so that
#     pipeline can flag CRs
#     """
#     import shutil
#     
#     if not os.path.exists(raw.replace('raw','raw_BKUP')):
#         shutil.copy(raw, raw.replace('raw','raw_BKUP'))
#     
#     img = pyfits.open(raw, mode='update')
#     
#     cube, dq, time, NSAMP = unicorn.prepare.split_multiaccum(img, scale_flat=False)
#     diff = np.diff(cube, axis=0)
#     ramp_cps = np.median(diff, axis=1)
#     avg_ramp = np.median(ramp_cps, axis=1)
#     
#     cumsum = np.cumsum(avg_ramp)
#     
#     for i in range(NSAMP-1):
#         img['sci', i+1].data -= cumsum[NSAMP-2-i]*10 #- 2.1*img['time',i+1].header['PIXVALUE']
#         
#     for i in range(NSAMP):
#         img['sci',i+1].data = np.cast[np.int16](img['sci',i+1].data)
#         
#     img.flush(output_verify='fix')
#     
#     wfc3tools.calwf3.calwf3(raw)
    
def make_IMA_FLT(raw='ibhj31grq_raw.fits', pop_reads=[], remove_ima=True, fix_saturated=True):
    """
    Run calwf3, if necessary, to generate ima & flt files.  Then put the last
    read of the ima in the FLT SCI extension and let Multidrizzle flag 
    the CRs.
    
    Optionally pop out reads affected by satellite trails or earthshine.  The 
    parameter `pop_reads` is a `list` containing the reads to remove, where
    a value of 1 corresponds to the first real read after the 2.9s flush.
    
    Requires IRAFX for wfc3tools
    """
    import wfc3tools
        
    ### Remove existing products or calwf3 will die
    for ext in ['flt','ima']:
        if os.path.exists(raw.replace('raw', ext)):
            os.remove(raw.replace('raw', ext))
    
    ### Run calwf3
    wfc3tools.calwf3.calwf3(raw)
    
    flt = pyfits.open(raw.replace('raw', 'flt'), mode='update')
    ima = pyfits.open(raw.replace('raw', 'ima'))
    
    cube, dq, time, NSAMP = unicorn.prepare.split_multiaccum(ima,
                                       scale_flat=False)
    
    readnoise_2D = np.zeros((1024,1024))
    readnoise_2D[512: ,0:512] += ima[0].header['READNSEA']
    readnoise_2D[0:512,0:512] += ima[0].header['READNSEB']
    readnoise_2D[0:512, 512:] += ima[0].header['READNSEC']
    readnoise_2D[512: , 512:] += ima[0].header['READNSED']
    readnoise_2D = readnoise_2D**2
        
    if len(pop_reads) > 0:
        ### Pop out reads affected by satellite trails or earthshine
        threedhst.showMessage('Pop reads %s from %s' %(pop_reads, ima.filename()))
        
        diff = np.diff(cube, axis=0)
        dt = np.diff(time)
        final_exptime = time[-1]
        final_sci = cube[-1,:,:]*1
        for read in pop_reads:
            final_sci -= diff[read,:,:]
            final_exptime -= dt[read]
        
        #final_var = ima[0].header['READNSEA']**2 + final_sci        
        final_var = readnoise_2D + final_sci        
        final_err = np.sqrt(final_var)/final_exptime
        final_sci /= final_exptime
        
        flt[0].header['EXPTIME'] = final_exptime
        
    else:
        final_sci = ima['SCI', 1].data*1
        final_err = ima['ERR', 1].data*1
    
    final_dq = ima['DQ', 1].data*1
    
    #### For saturated pixels, look for last read that was unsaturated
    #### Background will be different under saturated pixels but maybe won't
    #### matter so much for such bright objects.
    if fix_saturated:
        print 'Fix Saturated pixels:'
        #### Saturated pixels
        zi, yi, xi = np.indices(dq.shape)
        saturated = (dq & 256) > 0
        # 1024x1024 index array of reads where pixels not saturated
        zi_flag = zi*1
        zi_flag[saturated] = 0
        last_ok_read = np.max(zi_flag, axis=0)

        zi_idx = zi < 0
        for i in range(2, NSAMP-1):
            zi_idx[i,:,:] = zi[i,:,:] == last_ok_read

        time_array = time[zi]
        time_array[0,:,:] = 1.e-3 # avoid divide-by-zero
        # pixels that saturated before the last read
        fix = (last_ok_read < (ima[0].header['NSAMP'] - 1)) & (last_ok_read > 1)
        #err = np.sqrt(ima[0].header['READNSEA']**2 + cube)/time_array
        err = np.sqrt(readnoise_2D + cube)/time_array

        final_sci[fix] = np.sum((cube/time_array)*zi_idx, axis=0)[fix]
        final_err[fix] = np.sum(err*zi_idx, axis=0)[fix]

        fixed_sat = (zi_idx.sum(axis=0) > 0) & ((final_dq & 256) > 0)
        final_dq[fixed_sat] -= 256
        print '  Nsat = %d' %(fixed_sat.sum())
        flt['DQ'].data |= final_dq[5:-5,5:-5] & 256
        
    else:
        #### Saturated pixels
        flt['DQ'].data |= ima['DQ',1].data[5:-5,5:-5] & 256
        
    flt['SCI'].data = final_sci[5:-5,5:-5]
    flt['ERR'].data = final_err[5:-5,5:-5]
    
    
    #### Some earthshine flares DQ masked as 32: "unstable pixels"
    mask = (flt['DQ'].data & 32) > 0
    if mask.sum() > 1.e4:
        threedhst.showMessage('Take out excessive DQ=32 flags (N=%e)' %(mask.sum()))
        flt['DQ'].data[mask] -= 32
        
    ### Update the FLT header
    flt[0].header['IMA2FLT'] = (1, 'FLT extracted from IMA file')
    flt[0].header['IMASAT'] = (fix_saturated*1, 'Manually fixed saturation')

    # if pyfits.__version__ > '3.2':
    #     flt[0].header['IMA2FLT'] = (1, 'FLT extracted from IMA file')
    #     flt[0].header['IMASAT'] = (fix_saturated*1, 'Manually fixed saturation')
    # else:
    #     flt[0].header.update('IMA2FLT', 1, comment='FLT extracted from IMA file')
    #     flt[0].header.update('IMASAT', fix_saturated*1, comment='Manually fixed saturation')
        
    flt.flush()
    
    ### Remove the IMA file
    if remove_ima:
        os.remove(raw.replace('raw', 'ima'))

def show_MultiAccum_reads(raw='ib3701s4q_ima.fits'):
    """
    Make a figure (.ramp.png) showing the individual reads of an 
    IMA or RAW file.
    """    
    import matplotlib.pyplot as plt
    
    img = pyfits.open(raw)
    
    if 'raw' in raw:
        gain=2.5
    else:
        gain=1
        
    cube, dq, time, NSAMP = split_multiaccum(img)
    diff = np.diff(cube, axis=0)
    dt = np.diff(time)
    fig = plt.figure(figsize=[10,10])
    
    #(xs=10, aspect=0.8, wspace=0., hspace=0., left=0.05, NO_GUI=True)
    for j in range(1,NSAMP-1):
        ax = fig.add_subplot(4,4,j)
        ax.imshow(diff[j,:,:]/dt[j]*gain, vmin=0, vmax=4, origin='lower')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(20,5,'%d' %(j), ha='left', va='bottom')
    #
    ax = fig.add_subplot(4,4,16)
    ramp_cps = np.median(diff, axis=1)
    avg_ramp = np.median(ramp_cps, axis=1)
    ax.plot(time[2:], ramp_cps[1:,16:-16:4]/100*gain, alpha=0.1, color='black')
    ax.plot(time[2:], avg_ramp[1:]/100*gain, alpha=0.8, color='red', linewidth=2)
    
    fig.tight_layout(h_pad=0.0, w_pad=0.0, pad=0.0)
    plt.savefig(raw.split('_')[0]+'_ramp.png')
    
    return fig
    
def redo_pointing(pointing='GOODS-S-31', grism_exposures = [2,4]):
    
    os.chdir(os.getenv('THREEDHST')+'/3DHST_VariableBackgrounds/RAW')
    
    asn = threedhst.utils.ASNFile('../PREP_FLT/%s-G141_asn.fits' %(pointing)) 
    for exp in grism_exposures:
        unicorn.prepare.make_IMA_FLT(raw=asn.exposures[exp-1]+'_raw.fits')
    
    #
    os.chdir('../PREP_FLT/')
    
    #threedhst.process_grism.fresh_flt_files(pointing+'-F140W_asn.fits', from_path='../RAW', preserve_dq=False)
    
    #threedhst.prep_flt_files.startMultidrizzle(pointing+'-F140W_asn.fits', use_shiftfile=True, skysub=True, final_scale=0.06, pixfrac=0.8, driz_cr=True, updatewcs=False, clean=True, median=True)
    
    threedhst.prep_flt_files.process_3dhst_pair(pointing+'-F140W_asn.fits', pointing+'-G141_asn.fits', adjust_targname=False, ALIGN_IMAGE = None, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True, align_geometry='rotate, shift')
    
    if False:
        unicorn.reduce.set_grism_config(grism='G141', use_new_config=True, force=True)
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25., make_zeroth_model=False, BEAMS=['A','B','C','D','E'])

        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=21., REFINE_MAG_LIMIT=21, make_zeroth_model=False, BEAMS=['A','B','C','D','E'])
        
        id=27230
        model.twod_spectrum(id, refine=False, CONTAMINATING_MAGLIMIT=18, USE_REFERENCE_THUMB=True)
        model.show_2d(savePNG=True)

        gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(pointing, id), lowz_thresh=0.)
        gris.fit_in_steps(zrfirst=[np.clip(gris.z_peak-0.3, 0,3),gris.z_peak+0.3])

        os.system('open %s_%05d*zfit*png' %(pointing, id))
    
    
    
    
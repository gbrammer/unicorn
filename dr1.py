import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import pyfits

import threedhst
import threedhst.catIO as catIO
import threedhst.eazyPy as eazy
import threedhst.dq

import unicorn

def udf_prepare():
    """
    Make images and catalogs for extracting UDF spectra
    """
    import os
    import threedhst
    import unicorn
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.catIO as catIO
    
    os.chdir(unicorn.GRISM_HOME+'UDF/PREP_FLT')
    
    ALIGN = '../XDF/xdfh_sci.fits'
    ALIGN_EXT=0
    
    info = catIO.Readfile('files.info')
    
    #### Make ASN files for different grisms / orientations
    asn = threedhst.utils.ASNFile(glob.glob('../RAW/i*asn.fits')[0])
    
    ## 10-26
    for d in ['10-26','11-01']:
        match = (info.targname == 'PRIMO') & (info.filter == 'G141') & (info.date_obs == '2010-%s' %(d))
        asn.exposures = []
        for exp in info.file[match]:
            asn.exposures.append(exp.split('_flt')[0])
        #
        asn.product = 'PRIMO-%s-G141' %(d.replace('-',''))
        asn.write(asn.product+'_asn.fits', clobber=True)
        #
        match = (info.targname == 'PRIMO') & (info.filter != 'G141') & (info.date_obs == '2010-%s' %(d))
        filt = info.filter[match][0]
        asn.exposures = []
        for exp in info.file[match]:
            asn.exposures.append(exp.split('_flt')[0])
        #
        asn.product = 'PRIMO-%s-%s' %(d.replace('-',''), filt)
        asn.write(asn.product+'_asn.fits', clobber=True)
        
    for pointing in [34,36,37,38]:
        for filt in ['F140W','G141']:
            match = (info.targname == 'GOODS-SOUTH-%d' %(pointing)) & (info.filter == filt)
            asn.exposures = []
            for exp in info.file[match]:
                asn.exposures.append(exp.split('_flt')[0])
            #
            asn.product = 'GOODS-SOUTH-%d-%s' %(pointing, filt)
            asn.write(asn.product+'_asn.fits', clobber=True)
    
    
    ##### Run background subtraction on all images        
    direct = glob.glob('*[0-9]-F*asn.fits')
    grism = glob.glob('*[0-9]-G141_asn.fits')
    
    for i in range(len(direct)):
        if not os.path.exists(grism[i].replace('asn','drz')):
            pair(direct[i], grism[i], adjust_targname=False, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=ALIGN_EXT, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
    ### Fix offsets for 1026 since aperture combination was different
    files=['ibfup1myq_flt.fits','ibfup1n1q_flt.fits']
    for file in files:
        im = pyfits.open(file, mode='update')
        im[0].header.update('POSTARG1', im[0].header['POSTARG1']+8.814150)
        im[0].header.update('POSTARG2', im[0].header['POSTARG2']+0.025025)
        im.flush()
        
    #### Interlaced combinations
    ## Need to fake a combination for the interlaced direct image for PRIMO
    ## The first image is the direct image and the rest are G141 exposures
    ## to fill a 2x2 interlaced array
    ##
    ## Give them "F140W" filenames to work with interlacing code
    for d, f in zip(['1026', '1101'], ['F160W', 'F125W']):
        asn_im = threedhst.utils.ASNFile('PRIMO-'+d+'-%s_asn.fits' %(f))
        asn = threedhst.utils.ASNFile('PRIMO-'+d+'-G141_asn.fits')
        sf = threedhst.shifts.ShiftFile('PRIMO-'+d+'-G141_shifts.txt')
        #
        asn.exposures[0] = asn_im.exposures[0]
        sf.images[0] = asn.exposures[0]+'_flt.fits'
        #
        ###  Enough images to fill 2x2 grid: (Even-Even, Odd-Odd, OE, EO)
        # xo, yo = unicorn.reduce.get_interlace_offsets('PRIMO-'+d+'-G141_asn.fits', verbose=1, path_to_flt='./')
        keep = [0,1,3,5]
        for i in range(len(asn.exposures))[::-1]:
            if i not in keep:
                p = asn.exposures.pop(i)
                p = sf.images.pop(i)
                p = sf.xshift.pop(i)
                p = sf.yshift.pop(i)
                p = sf.scale.pop(i)
                p = sf.rotate.pop(i)
        #
        ### Image is combination of G141, F140W but need Multidrizzle outputs
        asn.product = 'PRIMO-'+d+'-F140W'
        sf.nrows = 4
        asn.write('%s_asn.fits' %(asn.product), clobber=True)
        sf.write('%s_shifts.txt' %(asn.product))
        threedhst.prep_flt_files.startMultidrizzle(asn.product + '_asn.fits',
                     use_shiftfile=True, skysub=False,
                     final_scale=0.06, pixfrac=0.8, driz_cr=False,
                     updatewcs=False, clean=True, median=False)
    
    #####################
    #### Deep F160W reference for interlaced reductions
    #####################
    
    wht = pyfits.open('xdfh_wht.fits')
    sci = pyfits.open('xdfh_sci.fits')
    texp = 3.e3
    f = 10**(-0.4*(33.4549980163574-25.96))
    #wht[0].data = wht[0].data*1000.-
    wht[0].data = (1.e5*wht[0].data**2+1./((sci[0].data*f+0.5)/texp))*f**2
    wht.writeto('xdfh_VAR.fits', clobber=True)
    wht[0].data = 1./np.sqrt(wht[0].data)
    wht.writeto('xdfh_SIG.fits', clobber=True)
    sci[0].data *= f
    sci.writeto('xdfh_sci_scaled.fits')
    
    ## Make catalog
    os.chdir("/research/HST/GRISM/3DHST/UDF/XDF")
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = 'xdf.cat'
    se.options['CHECKIMAGE_NAME'] = 'xdf_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'xdfh_VAR.fits'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '1.5'
    se.options['ANALYSIS_THRESH']  = '1.5'
    se.options['MAG_ZEROPOINT'] = '33.45499801'
    se.options['DEBLEND_NTHRESH'] = '64'
    se.options['DEBLEND_MINCONT'] = '0.00005'
    
    status = se.sextractImage('xdfh_sci.fits')
    
    threedhst.sex.sexcatRegions('xdf.cat', 'xdf.reg', format=2)
    
    ## Prep blot drizzle images
    REF_ROOT = 'XDF-F160W'
    CATALOG = '../XDF/xdf.cat'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F160W', REFERENCE = '../XDF/xdfh_sci_scaled.fits', SEGM = '../XDF/xdf_seg.fits', sci_extension=0)
    
    REF_ROOT = 'HUDF12-F160W'
    CATALOG = '../HUDF12/hudf12.cat'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F160W', REFERENCE = '../HUDF12/hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits', SEGM = '../HUDF12/hudf12_seg.fits', sci_extension=0)
    
    ## Use new F140W image
    REF_ROOT = 'HUDF12-F140W'
    CATALOG = '../F140W/HUDF12-F140W.cat'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F140W', REFERENCE = '../F140W/HUDF12-F140W_drz_sci.fits', SEGM = '../F140W/HUDF12-F140W_seg.fits', sci_extension=0)
    
    ### Generate DRZ images
    files=glob.glob('*F*asn.fits')
    for file in files: 
        threedhst.prep_flt_files.startMultidrizzle(file,
                     use_shiftfile=True, skysub=False,
                     final_scale=0.06, pixfrac=0.8, driz_cr=False,
                     updatewcs=False, clean=True, median=False)
    
        
    NGROW=125
    ROOT = 'PRIMO-1101'
    ROOT = 'GOODS-SOUTH-34'
    for p in [34, 36, 37, 38]:
        ROOT = 'GOODS-SOUTH-%d' %(p)
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = ROOT+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=ROOT+'-F140W', view=False, pad=60+200*(ROOT=='PRIMO-1026'), REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, auto_offsets=True, NSEGPIX=3)
        # seg = pyfits.open(ROOT+'_inter_seg.fits', mode='update')
        # seg[0].data[seg[0].data < 0] = 0
        # seg.flush()
        #
        unicorn.reduce.interlace_combine(root=ROOT+'-G141', view=False, pad=60+200*(ROOT=='PRIMO-1026'),  NGROW=NGROW, auto_offsets=True)
        unicorn.reduce.interlace_combine(root=ROOT+'-F140W', view=False, pad=60+200*(ROOT=='PRIMO-1026'),  NGROW=NGROW, auto_offsets=True)
        #
        ref = pyfits.open(ROOT+'_ref_inter.fits', mode='update')
        ref[0].header['FILTER'] = 'F140W'
        ref.flush()
    
    ### Shift 1026-G141 image right by 130 pixels 
    #xo, yo = unicorn.reduce.get_interlace_offsets('PRIMO-1026-F140W_asn.fits', verbose=1, path_to_flt='./')
    im = pyfits.open('PRIMO-1026-G141_inter.fits', mode='update')
    fill = im[1].data*0
    fill[:,480:2580] = im[1].data[:,350:2450]
    im[1].data = fill*1
    fill = im[2].data*0
    fill[:,480:2580] = im[2].data[:,350:2450]
    im[2].data = fill*1
    im.flush()
    
    #### Shift 1101-G141 image down by 1 pixel
    im = pyfits.open('PRIMO-1101-G141_inter.fits', mode='update')
    fill = im[1].data*0
    fill[100:-100,:] = im[1].data[101:-99,:]
    im[1].data = fill*1
    fill = im[2].data*0
    fill[100:-100,:] = im[2].data[101:-99,:]
    im[2].data = fill*1
    im.flush()
    
    if ROOT.startswith('PRIMO'):
        ref = pyfits.open(ROOT+'_ref_inter.fits')
        im140 = pyfits.open(ROOT+'-F140W_inter.fits', mode='update')
        im140[0].header['FILTER'] = 'F140W'
        im140[1].data = ref[1].data #/ 10**(-0.4*(26.46-25.96))
        #im140[2].data = im140[2].data # / 10**(-0.4*(26.46-25.96))
        im140.flush()
    
    #
    files = glob.glob('*ref_inter.fits')
    for file in files:
        ROOT=file.split('_ref_inter')[0]
        model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=28, REFINE_MAG_LIMIT=23, make_zeroth_model=False, grism='G141')
        model.make_wcs_region_file()
        
    #model = unicorn.reduce.GrismModel(root=ROOT, MAG_LIMIT=30, grism='G141')
    ### Force use F160W as detection image
    #if ROOT.startswith('PRIMO'):
    #model.direct[1].data = model.im[1].data*1./10**(-0.4*(26.46-25.96))
    
    model.get_corrected_wcs()
    model.make_wcs_region_file()
    
    #for p in [34,36,37,38]:
    #    ROOT = 'GOODS-SOUTH-%d' %(p)
    
    
    ##### Extract all objects
    c = threedhst.sex.mySexCat('../F140W/HUDF12-F140W.cat')
    c.write(c.filename.replace('.cat','.reform.cat'), reformat_header=True)
    c = catIO.Readfile('../F140W/HUDF12-F140W.reform.cat')
    ok = c.mag_auto < 26
    hudf.extract_all(c.number[ok], miny=-80)
    
    ### FOr some reason, some 2D files weren't extracted with 80 pix.  Redo those
    bad = []
    files=glob.glob('*2D.fits')
    for 
    
    for i in range(ok.sum()):
        id = c.number[ok][i]
        twod = glob.glob('[GP]*_%05d.2D.fits' %(id))
        ### temporary
        # if (len(twod) < 2) | os.path.exists('UDF_%05d_stack.png' %(id)):
        #     continue
        # #
        #### Redo background fits
        if len(twod) > 0:
            try:
                hudf.stack(id, dy=40, inverse=True)
                hudf.fix_2d_background('UDF_%05d' %(id), force=True)
                hudf.stack(id, dy=16, inverse=True)
                if c.mag_auto[ok][i] > 24:
                    hudf.fix_2d_background('UDF_%05d' %(id), force=False)
            except:
                pass
        #
        # twod = glob.glob('[GP]*_%05d.2D.fits' %(id))
        # ### Get just center of the trace
        # if len(twod) > 0:
        #     try:
        #         hudf.stack(id, dy=12, inverse=True)
        #     except:
        #         pass
    
    ######### Subtract background on smaller cutout
    for id in c.number[ok]:
        twod = glob.glob('[GP]*_%05d.2D.fits' %(id))
        if len(twod) > 0:
            try:
                hudf.stack(id, dy=18, inverse=True)
                hudf.fix_2d_background('UDF_%05d' %(id), force=False)
            except:
                pass
    
    for id in c.number[ok][::-1]:
        if (os.path.exists('UDF_%05d.new_zfit.png' %(id))) | (not os.path.exists('UDF_%05d.bg.png' %(id))):
            continue
        #
        try:
            unicorn.analysis.FORCE_GOODSS = False
            self = test.SimultaneousFit('UDF_%05d' %(id), RELEASE=False, p_flat=1.e-8, lowz_thresh=0.)
            if self.dr > 0.5:
                unicorn.analysis.FORCE_GOODSS = True
                self = test.SimultaneousFit('UDF_%05d' %(id), RELEASE=False, p_flat=1.e-8, lowz_thresh=0.)
            #
            self.read_master_templates()
            self.new_fit_constrained(faint_limit=25)
            #os.system('open UDF_%05d.new_zfit.png' %(id))
            #self.new_fit_free_emlines()
        except:
            pass
                 
def extract_new_redshifts(PATH='./'):
    """
    UDF 
    
    The fitting code saves a FITS file with p(z) but doesn't print the 
    results.  Pull it out and make a catalog
    """
    import unicorn.interlace_test as test
    import unicorn.catalogs2 as cat2
    
    c = catIO.Readfile('../F140W/HUDF12-F140W.reform.cat')
    zsp = cat2.SpeczCatalog()
    dr, idx = zsp.match_list(c.x_world, c.y_world)
    c.z_spec = zsp.zspec[idx] #[dr < 0.3]
    c.z_spec[dr > 0.3] = -1
    
    ok = c.mag_auto < 26
    
    original_path = os.getcwd()
    
    fp = open('udf_redshifts.dat','w')
    fp.write('# id mag flag z_max z_peak l68 u68 l95 u95 z_spec z_source\n')
    fp.write('# %s' %(PATH))
    
    os.chdir(PATH)
    
    for i in range(ok.sum()):
        id = c.number[ok][i]
        z_spec = c.z_spec[ok][i]
        if z_spec > 0:
            z_source = zsp.source[idx][ok][i]
        else:
            z_source = '-'
        #
        logstr = 'UDS_%05d  %.2f ' %(id, c.mag_auto[ok][i])
        if os.path.exists('UDF_%05d.new_zfit.pz.fits' %(id)):
            self = test.SimultaneousFit('UDF_%05d' %(id))
            self.read_master_templates()
            status = self.new_load_fits()
            if status:
                logstr += ' 1 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(self.z_max_spec, self.z_peak_spec, self.c68[0], self.c68[1], self.c95[0], self.c95[1], z_spec, z_source)
            else:
                logstr += '-1 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(-1, -1, -1, -1, -1, -1, -1, '-')
        else:
            logstr += ' 0 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(-1, -1, -1, -1, -1, -1, -1, '-')

        print logstr
        fp.write(logstr+'\n')
    
    fp.close()
    os.chdir(original_path)
    
def fix_UDF_nans():
    """
    UDF
    
    I had fixed a bug in the UDF stacking script that caused the thumbnails
    to be filled with NaNs if one of the four exposures had no valid pixels 
    in it.  Go through and find them and extract them again
    """
    
    files = glob.glob('UDF*2D.fits')
    for file in files:
        im = pyfits.open(file)
        if np.isfinite(im['DSCI'].data).sum() == 0:
            print file
            id=int(file.split('_')[1].split('.')[0])
            hudf.stack(id, dy=16, inverse=True)
            
def check_backgrounds():
    """
    UDF
    
    Plot the automatically-determined backgrounds as a function of position
    in the UDF frame
    """
    c = catIO.Readfile('../F140W/HUDF12-F140W.reform.cat')
    ok = c.mag_auto < 27.5
    
    fp = open('udf_backgrounds.dat','w')
    fp.write('# id c0 cx cy x y mag\n')
    for i in np.arange(c.N)[ok]:
        id = c.number[i]
        bgfile = 'UDF_FIT/UDF_%05d.bg.dat' %(id)
        if os.path.exists(bgfile):
            line = open(bgfile).readlines()[1][:-1]
            fp.write('%s  %7.1f %7.1f  %.2f\n' %(line, c.x_image[i], c.y_image[i], c.mag_auto[i]))
    
    fp.close()
    
    bg = catIO.Readfile('udf_backgrounds.dat')
    
    ok = (bg.mag > 24) & (bg.c0+0.003 > 0) & (bg.c0+0.003 < 0.004)
    plt.scatter(bg.x[ok], bg.y[ok], c=bg.c0[ok], s=30, vmin=-0.004, vmax=0.002)
    
    #### Try 2D as in sciypy example (need newer version)
    # from scipy.interpolate import griddata
    # 
    # points = [bg.x[ok], bg.y[ok]]
    # values = bg.c0[ok]
    # 
    # grid_z2 = griddata(points, values, (bg.x, bg.y), method='cubic')
    
    from scipy import interpolate
    bg_spline = interpolate.SmoothBivariateSpline(bg.x[ok], bg.y[ok], bg.c0[ok], kx=4, ky=4)
    
    test = bg.c0*0.
    for i in range(bg.N):
        test[i] = bg_spline(bg.x[i], bg.y[i])
    #
    plt.scatter(bg.x+3600, bg.y, c=test, s=30, vmin=-0.004, vmax=0.002)
        
    plt.text(2000,3800,'Data (mag > 24)', ha='center', va='center')
    plt.text(5600,3800,'Interpolated', ha='center', va='center')
    plt.savefig('background_spline.pdf')
    
    import cPickle as pickle
    fp = open('background_spline.pkl','wb')
    pickle.dump(bg_spline, fp)
    fp.close()

def new_background():
    """
    UDF
    
    Take out the old best-fit backgrounds that were applied automatically
    and subtract the interpolated version determined above
    """
    import unicorn.hudf as hudf
    import cPickle as pickle
    
    fp = open('background_spline.pkl','rb')
    bg_spline = pickle.load(fp)
    fp.close()
    
    files = glob.glob('UDF*2D.fits')
    for file in files:
        im = pyfits.open(file)
        header = im[0].header
        if 'BG_C0' in header.keys():
            hudf.fix_2d_background(object=file.split('.2D')[0], force=False, clip=100, remove=False)
        #
        new_bg = bg_spline(header['X_PIX'], header['Y_PIX'])
        im = pyfits.open(file, mode='update')
        im[0].header.update('BG_C0', new_bg)
        im[0].header.update('BG_CX', 0)
        im[0].header.update('BG_CY', 0)
        im['SCI'].data -= new_bg
        im.flush()
        
def rgb_mosaic_browser():
    
    #### UDF browser, note "f" parameter for deeper UDF images
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]

    f140 = pyfits.open('HUDF12-F140W_drz_sci.fits')
    match = ['HUDF12_F160W.fits', 'HUDF12_F125W.fits', 'UDF_ACS_i.fits']
    
    f125 = pyfits.open('HUDF12_F125W.fits')
    mask = f125[0].data == 0
    
    #### Fill empty HUDF12 images with wide-F140W
    for i in range(3):
        print match[i]
        im = pyfits.open(match[i], mode='update')
        scale = 10**(-0.4*(26.46-25.96))/scales[i]
        im[0].data[mask] = f140[0].data[mask]*scale
        im.flush()
    #
    f = 6
    rgb1 = 'HUDF12_F160W.fits[0]*%.3f, HUDF12_F125W.fits[0]*%.3f, UDF_ACS_i.fits[0]*%.3f' %(scales[0]*f, scales[1]*f*1., scales[2]*f*1.)
    f = 20
    rgb2 = 'HUDF12_F160W.fits[0]*%.3f, HUDF12_F125W.fits[0]*%.3f, UDF_ACS_i.fits[0]*%.3f' %(scales[0]*f, scales[1]*f*1., scales[2]*f*1.)
    
    #threedhst.gmap.makeImageMap([rgb1,rgb2], aper_list=[15, 16], tileroot=['HJY','deep'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    threedhst.gmap.makeImageMap([rgb1,rgb2], aper_list=[15,16,17], tileroot=['iJH', 'deep'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    #### iJH
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    
    #### UDS
    #acs_f814w = '/3DHST/Ancillary/UDS/CANDELS/public/hlsp_candels_hst_acs_uds-tot_f814w_v1.0_drz.fits'
    wfc3_f125w = '/3DHST/Ancillary/UDS/CANDELS/ASTRODRIZZLE/uds-f125w-astrodrizzle-v0.1_drz_sci.fits'
    wfc3_f160w = '/3DHST/Ancillary/UDS/CANDELS/ASTRODRIZZLE/uds-f160w-astrodrizzle-v0.1_drz_sci.fits'
    acs_f814w = '/3DHST/Photometry/Work/UDS/v3/images/UDS_F814W_sci.fits'
    
    x0, y0, N = 30720/2, 12800/2, 500  # UDS
    x0, y0, N = 20480/2, 20480/2, 500  # GOODS-N
    x0, y0, N = 23000/2, 20000/2, 500  # GOODS-S
    iraf.imcopy('%s[0][%d:%d,%d:%d]' %(wfc3_f160w, x0-N, x0+N-1, y0-N, y0+N-1), 'sub_f160w.fits')
    iraf.imcopy('%s[0][%d:%d,%d:%d]' %(wfc3_f125w, x0-N, x0+N-1, y0-N, y0+N-1), 'sub_f125w.fits')
    iraf.imcopy('%s[0][%d:%d,%d:%d]' %( acs_f814w, x0-N, x0+N-1, y0-N, y0+N-1), 'sub_f814w.fits')
    #     
    # wfc3_f125w, wfc3_f160w, acs_f814w = 'sub_f125w.fits', 'sub_f160w.fits', 'sub_f814w.fits'
    
    f = 2
    
    #scales[2] = 10**(-0.4*(25.94-25.96))*1.5/2
    rgb = '%s[0]*%.3f, %s[0]*%.3f, %s[0]*%.3f' %(wfc3_f160w, scales[0]*f, wfc3_f125w, scales[1]*f, acs_f814w, scales[2]*f)
    
    threedhst.gmap.makeImageMap([rgb], aper_list=[14,15,16], tileroot=['iJH'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    
    #### GOODS-S
    wfc3_f125w = '/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f125w-astrodrizzle-v0.1_drz_sci.fits'
    wfc3_f160w = '/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v0.1_drz_sci.fits'
    acs_f814w = '/3DHST/Photometry/Work/GOODS-S/v3/images/GOODS-S_candels_f814w_sci.fits'
    
    threedhst.shifts.matchImagePixels(input=glob.glob(acs_f814w), matchImage=wfc3_f160w, output='GOODS-S_F814W_match3.0.fits', input_extension=0, match_extension=0)
    acs_f814w = 'GOODS-S_F814W_match3.0.fits'
    
    f = 2
    
    #scales[2] = 10**(-0.4*(25.94-25.96))*1.5/2
    rgb = '%s[0]*%.3f, %s[0]*%.3f, %s[0]*%.3f' %(wfc3_f160w, scales[0]*f, wfc3_f125w, scales[1]*f, acs_f814w, scales[2]*f)
    
    threedhst.gmap.makeImageMap([rgb], aper_list=[14,15,16], tileroot=['iJH'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    #### GOODS-N
    wfc3_f125w = '/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f125w-astrodrizzle-v0.1_drz_sci.fits'
    wfc3_f160w = '/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f160w-astrodrizzle-v0.1_drz_sci.fits'
    threedhst.shifts.matchImagePixels(input=glob.glob('/3DHST/Ancillary/GOODS-N/GOODS_ACS/h_ni_sect*_v2.0_drz_img.fits'), matchImage=wfc3_f160w, output='GOODS-N_F775W_match3.0.fits', input_extension=0, match_extension=0)
    acs_f814w = 'GOODS-N_F775W_match3.0.fits'
    
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.666-25.96))*1.5]
    
    rgb = '%s[0]*%.3f, %s[0]*%.3f, %s[0]*%.3f' %(wfc3_f160w, scales[0]*f, wfc3_f125w, scales[1]*f, acs_f814w, scales[2]*f)
    
    threedhst.gmap.makeImageMap([rgb], aper_list=[14,15,16], tileroot=['iJH'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    #### cosmos
    wfc3_f125w = '/3DHST/Ancillary/COSMOS/CANDELS/ASTRODRIZZLE/cosmos-f125w-astrodrizzle-v0.1_drz_sci.fits'
    wfc3_f160w = '/3DHST/Ancillary/COSMOS/CANDELS/ASTRODRIZZLE/cosmos-f160w-astrodrizzle-v0.1_drz_sci.fits'
    threedhst.shifts.matchImagePixels(input=glob.glob('/3DHST/Ancillary//COSMOS/ACS/acs_I_030mas_*_sci.fits'), matchImage=wfc3_f160w, output='COSMOS_F814W_match3.0.fits', input_extension=0, match_extension=0)
    acs_f814w = 'COSMOS_F814W_match3.0.fits'
    
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    
    rgb = '%s[0]*%.3f, %s[0]*%.3f, %s[0]*%.3f' %(wfc3_f160w, scales[0]*f, wfc3_f125w, scales[1]*f, acs_f814w, scales[2]*f)
    
    threedhst.gmap.makeImageMap([rgb], aper_list=[14,15,16], tileroot=['iJH'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    #### AEGIS
    wfc3_f125w = '/3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f125w-astrodrizzle-v0.1_drz_sci.fits'
    wfc3_f160w = '/3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f160w-astrodrizzle-v0.1_drz_sci.fits'
        
    # swarp /3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f160w-astrodrizzle-v0.1_drz_sci.fits[0] -c make_square.swarp
    # mv coadd.fits AEGIS_F160W_v0.1_sq.fits; rm coadd.weight.fits
    # swarp /3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f125w-astrodrizzle-v0.1_drz_sci.fits[0] -c make_square.swarp
    # mv coadd.fits AEGIS_F125W_v0.1_sq.fits; rm coadd.weight.fits
    wfc3_f125w, wfc3_f160w = 'AEGIS_F125W_v0.1_sq.fits', 'AEGIS_F160W_v0.1_sq.fits' 
    threedhst.shifts.matchImagePixels(input = glob.glob('/3DHST/Ancillary/AEGIS/ACS/mos_i_scale1_section*[0-9]_drz.fits'), matchImage=wfc3_f160w, output='AEGIS_F814W_match3.0_sq.fits', input_extension=0, match_extension=0)
    acs_f814w = 'AEGIS_F814W_match3.0_sq.fits'
    
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    
    rgb = '%s[0]*%.3f, %s[0]*%.3f, %s[0]*%.3f' %(wfc3_f160w, scales[0]*f, wfc3_f125w, scales[1]*f, acs_f814w, scales[2]*f)
    
    threedhst.gmap.makeImageMap([rgb], aper_list=[14,15,16], tileroot=['iJH'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
            
            
import os
import glob
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

import threedhst
import unicorn
import threedhst.prep_flt_files
from threedhst.prep_flt_files import process_3dhst_pair as pair

def SN_MARSHALL():
    ####********************************************####
    ####              SN-MARSHALL (UDS)
    ####********************************************####
    
    
    os.chdir(unicorn.GRISM_HOME+'SN-MARSHALL/PREP_FLT')
    ALIGN_IMAGE = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds01_f160w_v0.5_drz.fits'
    
    #### NEW process all direct imaging
    unicorn.candels.make_asn_files()
    files=glob.glob('[UM]*-F1[26]*asn.fits')
    FORCE=False
    
    #### First processing step
    for file in files:
        if ((not os.path.exists(file.replace('asn','drz'))) & (not os.path.exists(file.replace('asn','drz')+'.gz'))) | FORCE:
            #threedhst.process_grism.fresh_flt_files(file)
            #threedhst.prep_flt_files.threedhst.shifts.run_tweakshifts(file)
            unicorn.candels.prep_candels(asn_file=file, 
                ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
                GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
                SCALE=0.06, geometry='rotate,shift')
    
    ### DQ mask
    bad = ['MARSHALL-B6-F125W', 'MARSHALL-1-G141', 'MARSHALL-3-G141', 'MARSHALL-4-G141', 'MARSHALL-6-G141', 'MARSHALL-6-F125W']
    for b in bad[1:]:
        threedhst.dq.checkDQ(asn_direct_file='%s_asn.fits' %(b), asn_grism_file='%s_asn.fits' %(b), path_to_flt='../RAW/', SIMPLE_DS9=False)
    
    for b in bad:
        unicorn.candels.prep_candels(asn_file='%s_asn.fits' %(b), 
            ALIGN_IMAGE = ALIGN_IMAGE, ALIGN_EXTENSION=0,
            GET_SHIFT=False, DIRECT_HIGHER_ORDER=1,
            SCALE=0.06, geometry='rotate,shift')
        
    #### Align to CANDELS reference
    files=glob.glob('[UM]*-F1*drz.fits')
    ALIGN_IMAGE = '/3DHST/Ancillary/UDS/CANDELS/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'
    for file in files:
        root=file.split('_drz')[0]
        if (not os.path.exists(root+'_align.map')) | FORCE:
            threedhst.prep_flt_files.startMultidrizzle('%s_asn.fits' %(root),
                     use_shiftfile=True, skysub=False,
                     final_scale=0.06, pixfrac=0.6, driz_cr=False,
                     updatewcs=False, clean=True, median=False)
            #
            for geom in ['rotate', 'shift']:
                threedhst.shifts.refine_shifts(ROOT_DIRECT=root, 
                    ALIGN_IMAGE=ALIGN_IMAGE, ALIGN_EXTENSION=0,  
                    fitgeometry=geom, clean=True)
                #
                threedhst.prep_flt_files.startMultidrizzle('%s_asn.fits' %(root),
                         use_shiftfile=True, skysub=False,
                         final_scale=0.06, pixfrac=0.6, driz_cr=False,
                         updatewcs=False, clean=True, median=False)
    
    ###
    
    files=glob.glob('*-F1[26]*asn.fits')
    for file in files:
        asn = threedhst.utils.ASNFile(file)
        for exp in asn.exposures:
            flt = pyfits.open(exp+'_flt.fits', mode='update')
            bad = (flt['SCI'].data > 1.4e4) & ((flt['DQ'].data & 4096) == 0)
            print exp, bad.sum()
            if bad.sum() == 0:
                continue
            #
            flt['DQ'].data[bad] |= 4096
            flt.flush()
    #
    import shutil
    for filter in ['F125W','F160W'][::-1]:
        threedhst.utils.combine_asn_shifts(glob.glob('[MU]*-%s_asn.fits' %(filter)),
            out_root='SN-MARSHALL-%s' %(filter), path_to_FLT='./', run_multidrizzle=False)
        #
        threedhst.prep_flt_files.startMultidrizzle('SN-MARSHALL-%s_asn.fits' %(filter),
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False, ra=34.439064, dec=-5.2595809,
                 final_outnx=3300, final_outny=3300, ivar_weights=True, build_drz=False)
        #
        shutil.move('SN-MARSHALL-%s_drz_weight.fits' %(filter), 'SN-MARSHALL-%s_drz_ivar.fits' %(filter))
        #
        threedhst.prep_flt_files.startMultidrizzle('SN-MARSHALL-%s_asn.fits' %(filter),
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False, ra=34.439064, dec=-5.2595809,
                 final_outnx=3300, final_outny=3300, ivar_weights=False, build_drz=False)
        #
        shutil.move('SN-MARSHALL-%s_drz_weight.fits' %(filter), 'SN-MARSHALL-%s_drz_exp.fits' %(filter))
        
    #### Process the grism exposures
    direct = glob.glob('../RAW/ib*80_asn.fits')
    grism = glob.glob('../RAW/ib*70_asn.fits')
    
    direct = glob.glob('MARSHALL-?-F*_asn.fits')
    grism = glob.glob('MARSHALL-?-G141_asn.fits')
    
    #for i in range(0,6):  
    for i in [2,3,5]:      
        threedhst.shifts.make_grism_shiftfile(direct[i], grism[i])
        pair(direct[i], grism[i], ALIGN_IMAGE = None, SKIP_GRISM=False, GET_SHIFT=False, SKIP_DIRECT=True)
            
    # #### Direct images
    # asn_list = glob.glob('MARSHALL-?-F160W_asn.fits')
    # threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-F160W',
    #                 path_to_FLT='./', run_multidrizzle=False)
    # threedhst.prep_flt_files.startMultidrizzle('MARSHALL-F160W_asn.fits',
    #              use_shiftfile=True, skysub=False,
    #              final_scale=0.06, pixfrac=0.8, driz_cr=False,
    #              updatewcs=False, clean=True, median=False)
    # 
    # asn_list = glob.glob('MARSHALL-?-F125W_asn.fits')
    # threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-F125W',
    #                 path_to_FLT='./', run_multidrizzle=False)
    # threedhst.prep_flt_files.startMultidrizzle('MARSHALL-F125W_asn.fits',
    #              use_shiftfile=True, skysub=False,
    #              final_scale=0.06, pixfrac=0.8, driz_cr=False,
    #              updatewcs=False, clean=True, median=False)
    
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
    #
    #### Placeholder direct image for each orientation
    asn_list = glob.glob('MARSHALL-[1]-F*_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-225-F140W',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-225-F140W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    #
    asn_list = glob.glob('MARSHALL-[2]-F*_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-245-F140W',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-245-F140W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
def make_catalog():
    
    os.chdir('/3DHST/Spectra/Work/SN-MARSHALL/Catalog')
    
    #### Clean up remaining bad pixels with large values
    drz = pyfits.open('SN-MARSHALL-F160W_drz_sci.fits', mode='update')
    exp = pyfits.open('SN-MARSHALL-F160W_drz_exp.fits', mode='update')
    ivar = pyfits.open('SN-MARSHALL-F160W_drz_ivar.fits', mode='update')
    
    bad = drz[0].data > 1000
    drz[0].data[bad] = 0
    drz.flush()
    ivar[0].data[bad] = 1.e-10
    ivar.flush()
    exp[0].data[bad] = 0.
    exp.flush()
    
    #### Some bad pixels
    # import scipy.ndimage as nd
    # med = nd.median_filter(drz[0].data+0.1, size=(3,3))
    # yi, xi = np.indices(bad.shape)
    # R = np.sqrt((xi-2005)**2+(yi-3331)**2)
    # bad = ((R < 60) & ((drz[0].data+0.1)/med > 1.2))*1.
    # grow = nd.convolve(bad, weights=np.ones((3,3))/9.)
    # drz[0].data[grow > 0] = 0
    # ivar[0].data[grow > 0] /= 10
    # drz.flush()
    # ivar.flush()
    
    #### Make sigma image
    counts = drz[0].data*exp[0].data
    GAIN = 2.5
    # sigma_poisson = np.sqrt(counts/s*t_exp/GAIN)/t_exp
    poisson =  (drz[0].data)/GAIN/exp[0].data # squared for variance
    poisson[counts < 10] = 0
    zero = ivar[0].data == 0
    var = 1./ivar[0].data + poisson
    var[zero] = 0
    var[var == 0] = 1.e6
    
    #### Region of spurious sources
    # import pyregion
    # reg_file = 'rms_mask.reg'
    # r = pyregion.open(reg_file).as_imagecoord(header=drz[0].header)
    # aper = r.get_mask(hdu=drz[0])
    # var[aper > 0] *= 3**2
    
    pyfits.writeto('SN-MARSHALL-F160W_drz_rms.fits', header=drz[0].header, data=np.sqrt(var), clobber=True)
    
    #### Make catalog
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['CHECKIMAGE_NAME'] = 'SN-MARSHALL-F160W_seg.fits'
    se.options['WEIGHT_TYPE']     = 'MAP_RMS'
    se.options['WEIGHT_IMAGE']    = 'SN-MARSHALL-F160W_drz_rms.fits[0]'
    se.options['WEIGHT_GAIN'] = 'N'
    se.options['GAIN'] = '0'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '1.1' 
    se.options['ANALYSIS_THRESH']  = '1.1' 
    
    se.options['DEBLEND_NTHRESH']    = '32' 
    se.options['DEBLEND_MINCONT']  = '0.0001' 
    #
    #apertures = np.arange(0.2,1.5,0.1)
    #apertures = np.append(apertures, np.arange(1.5,3,0.3) )
    #apertures = np.append(apertures, np.arange(3,8,0.6) )
    apertures = [0.2]
    str_ap = ['%.3f' %(ap/0.06) for ap in apertures]
    se.options['PHOT_APERTURES'] = ','.join(str_ap)
    
    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['MAG_ZEROPOINT'] = '%.2f' %(unicorn.reduce.ZPs['F160W'])
    se.options['CATALOG_NAME']    = 'SN-MARSHALL-F160W.cat'
    status = se.sextractImage('SN-MARSHALL-F160W_drz_sci.fits')
    threedhst.sex.sexcatRegions('SN-MARSHALL-F160W.cat', 'SN-MARSHALL-F160W.reg', format=2)
    
    c = threedhst.sex.mySexCat('SN-MARSHALL-F160W.cat')
    c.write(c.filename.replace('.cat','.reform.cat'), reformat_header=True)
    #c = catIO.Readfile('SN-MARSHALL-F160W.reform.cat')
    
def interlaced_products():
    
    os.chdir("/3DHST/Spectra/Work/SN-MARSHALL/PREP_FLT")
    
    # unicorn.reduce.interlace_combine(root='MARSHALL-225-F140W', view=True, use_error=True, make_undistorted=False, pad = 60, NGROW=125, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True)
    # unicorn.reduce.interlace_combine(root='MARSHALL-225-G141', view=True, use_error=True, make_undistorted=False, pad = 60, NGROW=125, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=True)
    
    REF_ROOT = 'MARSHALL-F160W'
    unicorn.reduce.prepare_blot_reference(REF_ROOT=REF_ROOT, filter='F160W', REFERENCE = '../Catalog/SN-MARSHALL-F160W_drz_sci.fits', SEGM = '../Catalog/SN-MARSHALL-F160W_seg.fits', sci_extension=0)
    
    #### Make Fake list of F140W files using the G141 list
    for p in ['225','245']:
        asn = threedhst.utils.ASNFile('../PREP_FLT/MARSHALL-%s-F140W_asn.fits' %(p))
        sf = threedhst.shifts.ShiftFile('MARSHALL-%s-G141_shifts.txt' %(p))
        while(len(sf.images) > 5):
            out = sf.pop(5)
        #
        sf.images[0] = asn.exposures[0]+'_flt.fits'
        sf.headerlines[1] = '# refimage: %s_flt.fits[1]\n' %(asn.exposures[0])
        asn.exposures = []
        for im in sf.images:
            asn.exposures.append(im.split('_flt')[0])
        #
        sf.write('MARSHALL-%s-F140W_shifts.txt' %(p))
        asn.write('MARSHALL-%s-F140W_asn.fits' %(p))
        #
        threedhst.prep_flt_files.startMultidrizzle('MARSHALL-%s-F140W_asn.fits' %(p),
                     use_shiftfile=True, skysub=False,
                     final_scale=0.06, pixfrac=0.8, driz_cr=False,
                     updatewcs=False, clean=True, median=False)
        
    #### Make reference images with the matched catalog
    CATALOG = '../Catalog/SN-MARSHALL-F160W.cat'
    NGROW=125
    for p in ['225','245']:
        ROOT = 'MARSHALL-%s' %(p)
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = ROOT+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=ROOT+'-F140W', view=False, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, auto_offsets=True, NSEGPIX=3)
        # seg = pyfits.open(ROOT+'_inter_seg.fits', mode='update')
        # seg[0].data[seg[0].data < 0] = 0
        # seg.flush()
        #
        unicorn.reduce.interlace_combine(root=ROOT+'-G141', view=False, pad=60,  NGROW=NGROW, auto_offsets=True)
        unicorn.reduce.interlace_combine(root=ROOT+'-F140W', view=False, pad=60,  NGROW=NGROW, auto_offsets=True)
        #
        ref = pyfits.open(ROOT+'_ref_inter.fits', mode='update')
        ref[0].header['FILTER'] = 'F160W'
        ref.flush()
        #
        im = pyfits.open(ROOT+'-F140W_inter.fits', mode='update')
        im[1].data = ref[1].data
        im.flush()
        #
    
    ### 245 off by a pixel for some reason.  Shift the grism image LL by (1,1)
    import scipy.ndimage as nd
    im = pyfits.open('MARSHALL-225-G141_inter.fits', mode='update')
    im[1].data[:,:-1] = im[1].data[:,1:]
    im[2].data[:,:-1] = im[2].data[:,1:]
    im.flush()

    im = pyfits.open('MARSHALL-245-G141_inter.fits', mode='update')
    im[1].data[:-1,:-1] = im[1].data[1:,1:]
    im[2].data[:-1,:-1] = im[2].data[1:,1:]
    im.flush()
    
    ### more shifts
    im = pyfits.open('MARSHALL-225-G141_inter.fits', mode='update')
    im[1].data[:,:-1] = im[1].data[:,1:]*1
    im[2].data[:,:-1] = im[2].data[:,1:]*1
    im.flush()

    im = pyfits.open('MARSHALL-245-G141_inter.fits', mode='update')
    im[1].data[:,2:] = im[1].data[:,:-2]*1
    im[2].data[:,2:] = im[2].data[:,:-2]*1
    im.flush()
    
    ###
    bl = pyfits.open('ibfuw2aaq_blot.fits')[0].data[:,125:-125]
    fl = pyfits.open('ibfuw2aaq_flt.fits')[1].data

    bl = pyfits.open('ibfuw1q4q_blot.fits')[0].data[:,125:-125]
    fl = pyfits.open('ibfuw1q4q_flt.fits')[1].data
    
    ##### cd ../Interlace; cp../PREP_FLT/MARSHALL-2?5* .
    os.chdir("/3DHST/Spectra/Work/SN-MARSHALL/Interlace")
    
    for p in ['225','245']:
        ROOT = 'MARSHALL-%s' %(p)
        model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=27.5, REFINE_MAG_LIMIT=23, make_zeroth_model=False, grism='G141')
        model.make_wcs_region_file()
    
    ### Check background
    for p in ['225','245']:
        gris = pyfits.open('MARSHALL-%s-G141_inter.fits' %(p), mode='update')
        model = pyfits.open('MARSHALL-%s_inter_model.fits' %(p))[0].data
        ### fit background
        ok = (model < 1.e-3) & (gris[1].data < 1.e-2) & (gris[1].data != 0)
        yi, xi = np.indices(ok.shape)
        from scipy import interpolate
        #bg_spline = interpolate.SmoothBivariateSpline(xi[ok], yi[ok], gris[ok], kx=3, ky=3)
        skern = 10
        xkern = np.arange(-5*skern, 5*skern+0.1)
        kernel = 1./np.sqrt(2*np.pi*skern**2)*np.exp(-xkern**2/2/skern**2)
        #
        i = 0 # x
        n = ok.sum(axis=i)
        flux = (gris[1].data*ok).sum(axis=i)
        avg = flux/np.maximum(n,1)
        sm_avg = nd.convolve1d(avg, kernel)
        sm_avg[n == 0] = 0
        plt.plot(avg, alpha=0.5)
        plt.plot(sm_avg, alpha=0.8)
        sm_2d = np.dot(np.ones((ok.shape[0],1)), sm_avg.reshape((1,-1)))
        sm_2d[gris[1].data == 0] = 0
        gris[1].data -= sm_2d
        gris[2].data = np.sqrt(gris[2].data**2+sm_2d**2)
        gris.flush()
        
def extract_spectra():
    ### check model
    # bl = pyfits.open('ibfuw1q4q_blot.fits')[0].data[:,125:-125]
    # fl = pyfits.open('ibfuw1q4q_flt.fits')[1].data
     
    ### test
    import unicorn.catalogs2 as cat2
    zsp = cat2.SpeczCatalog()
    dr, idx = zsp.match_list(model.ra_wcs, model.dec_wcs)
    ok = dr < 0.5
    
    #
    p = '225'
    ROOT = 'MARSHALL-%s' %(p)
    model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=27.5, REFINE_MAG_LIMIT=23, make_zeroth_model=False, grism='G141')
    
    import unicorn.interlace_test as test
    id=723
    # model.twod_spectrum(id, USE_REFERENCE_THUMB=True, miny=-80)
    # hudf.fix_2d_background('%s_%05d' %(model.root, id))
    model.twod_spectrum(id, USE_REFERENCE_THUMB=True, maxy=26)
    #hudf.fix_2d_background('%s_%05d' %(model.root, id))
    #
    self = test.SimultaneousFit('%s_%05d' %(model.root, id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.0, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True)
    self.read_master_templates()
    self.new_load_fits()
    self.new_fit_constrained(faint_limit=25)
    
    #### Extract all 1D/2D spectra
    for p in ['225','245']:
        ROOT = 'MARSHALL-%s' %(p)
        model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=27.5, REFINE_MAG_LIMIT=23, make_zeroth_model=False, grism='G141')
        ok = model.cat.mag < 26.8
        for id in model.objects[ok]:
            if os.path.exists('%s_%05d.2D.fits' %(ROOT, id)):
                continue
            #
            status = model.twod_spectrum(id, USE_REFERENCE_THUMB=True, maxy=26)
            print id
    
    ### Combine the two pointings
    import unicorn.interlace_test as test
    import unicorn.hudf as hudf
    
    #### Stack for all objects with 2D spectra
    files=glob.glob('MARSHALL-2*2D.fits')
    ids = []
    for file in files:
        ids.append(int(file[13:18]))
    #
    ids = np.unique(ids)
    
    for id in ids[2::3]:
        if os.path.exists('MARSHALL_%05d_stack.png' %(id)):
            print id
            continue
        #
        hudf.stack(id=id, dy=12, save=True, inverse=True, scale=[1,90], fcontam=2., ref_wave = 1.4e4, root='MARSHALL', search='MARSHALL-2[24]5')
        #
        self = test.SimultaneousFit('MARSHALL_%05d' %(id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.0, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True)
        if self.status is False:
            continue
        #
        if self.dr > 0.5:
            continue
        #
        self.read_master_templates()
        self.new_load_fits()
        try:
            self.new_fit_constrained(faint_limit=25.5, get_tilt=(self.twod.im[0].header['MAG'] < 25.5))
        except:
            self.read_master_templates()
            self.new_load_fits()
            try:
                self.new_fit_constrained(get_tilt=False)
            except:
                pass
        #    
    
    ### test single extraction    
    # id=150
    # 
    # hudf.stack(id=id, dy=12, save=True, inverse=True, scale=[1,90], fcontam=2., ref_wave = 1.4e4, root='MARSHALL', search='MARSHALL-2[24]5')
    # #
    # self = test.SimultaneousFit('MARSHALL_%05d' %(id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.0, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True)
    # self.read_master_templates()
    # self.new_load_fits()
    # try:
    #     self.new_fit_constrained(faint_limit=25.5, get_tilt=(self.twod.im[0].header['MAG'] < 25.5))
    # except:
    #     self.read_master_templates()
    #     self.new_load_fits()
    #     self.new_fit_constrained(get_tilt=False)
    # #    
    # os.system('open MARSHALL_%05d*[tk].*png' %(id))
    # 
    # 
    # if status is False:
    #         #
    #         self = test.SimultaneousFit('%s_%05d' %(model.root, id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.0, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True)
    #         # if self.status is False:
    #         #     continue
    #         # #
    #         # if self.dr > 0.5:
    #         #     continue
    #         # #
    #         self.read_master_templates()
    #         self.new_load_fits()
    #         self.new_fit_constrained(faint_limit=24.5, get_tilt=(self.twod.im[0].header['MAG'] < 25.5))
        
#
def extract_redshifts(PATH='./'):
    """ 
    
    The fitting code saves a FITS file with p(z) but doesn't print the 
    results.  Pull it out and make a catalog
    """
    import unicorn.interlace_test as test
    import unicorn.catalogs2 as cat2
    from threedhst import catIO
    
    c = catIO.Readfile('../Catalog/SN-MARSHALL-F160W.reform.cat')
    zsp = cat2.SpeczCatalog()
    dr, idx = zsp.match_list(c.x_world, c.y_world)
    c.z_spec = zsp.zspec[idx] #[dr < 0.3]
    c.z_spec[dr > 0.8] = -1
    
    ok = c.mag_auto < 26.5
    
    #original_path = os.getcwd()
    
    fp = open('marshall_redshifts.dat','w')
    fp.write('# id mag flag z_max z_peak l68 u68 l95 u95 z_spec z_source\n')
    fp.write('# %s' %(PATH))
    
    #os.chdir(PATH)
    
    for i in range(ok.sum()):
        id = c.number[ok][i]
        z_spec = c.z_spec[ok][i]
        if z_spec > 0:
            z_source = zsp.source[idx][ok][i]
        else:
            z_source = '-'
        #
        logstr = 'MARSHALL_%05d  %.2f ' %(id, c.mag_auto[ok][i])
        if os.path.exists('MARSHALL_%05d.new_zfit.pz.fits' %(id)):
            self = test.SimultaneousFit('MARSHALL_%05d' %(id))
            self.read_master_templates()
            status = self.new_load_fits()
            if status:
                logstr += ' 1 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(self.z_max_spec, self.z_peak_spec, self.c68[0], self.c68[1], self.c95[0], self.c95[1], z_spec, z_source)
            else:
                logstr += '-1 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(-1, -1, -1, -1, -1, -1, -1, '-')
        else:
            logstr += ' 0 %.4f %.4f %.4f %.4f %.4f %.4f  %.4f %s' %(-1, -1, -1, -1, -1, -1, -1, '-')
        #
        print logstr
        fp.write(logstr+'\n')
    
    fp.close()
    #os.chdir(original_path)
    
def check_redshifts():
    
    import astropy.cosmology as cc
    cosmo = cc.LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    
    c = catIO.Readfile('../Catalog/SN-MARSHALL-F160W.reform.cat')
    cok = c.number < -100
    zc = catIO.Readfile('marshall_redshifts.dat')
    for id in zc.id:
        idi = int(id[-5:])
        #print idi, cok.sum()
        cok = cok | (c.number == idi)
    #
    star = (c.mag_auto[cok] < 24) & (c.flux_radius[cok] < 2.6)
    
    cat, zout, fout = unicorn.analysis.read_catalogs('UDS')
    mat = catIO.CoordinateMatcher(cat)
    dr, idx = mat.match_list(c.x_world, c.y_world)
    logm, logm_z, z_phot = fout.lmass[idx][cok], fout.z[idx][cok], zout.z_peak[idx][cok]
    
    ok = (zc.z_max > 0) & ~star
    dz95 = (zc.u95-zc.l95)/(1+zc.z_max)
    dz68 = (zc.u68-zc.l68)/(1+zc.z_max)
    best = ok & (dz68 < 0.015)
    
    #plt.scatter(zc.z_max[best], zc.mag[best], alpha=0.5)
    #plt.scatter(zc.z_max[best], logm[best], alpha=0.5)
    
    zr = [0,4]
    zrange = np.log(1+np.array(zr))
    dz = 0.01
    nbins = zrange[1]/dz
    h = np.histogram(np.log(1+zc.z_max[best]), bins=nbins, range=zrange)
    h_idx = np.digitize(np.log(1+zc.z_max), h[1])
    z = np.exp(h[1][:-1]*0.5+h[1][1:]*0.5)-1
    cmv = h[1]*0.
    for i in range(len(h[1])):
        cmv[i] = cosmo.comoving_volume(h[1][i])/(4*np.pi*(360./2/np.pi)**2)
        
    survey_area = (2.2/60.)**2 #sq deg
    dcmv = np.diff(cmv)*survey_area
    
    nh = np.maximum(h[0], 0.01)
    plt.plot(z, nh, linestyle='steps')
        
    peak = np.abs(np.log(1+zc.z_max)-np.log(1+0.904)) < 0.005

    peak = np.abs(np.log(1+zc.z_max)-np.log(1+1.912)) < 0.005

    peak = np.abs(np.log(1+zc.z_max)-np.log(1+2.310)) < 0.005
    
    zbin = 0.904; bin_id = np.arange(len(h[1]))[z >= zbin][0]+1
    zbin = 1.912; bin_id = np.arange(len(h[1]))[z >= zbin][0]+1
    zbin = 1.522; bin_id = np.arange(len(h[1]))[z >= zbin][0]+1
    zbin = 2.310; bin_id = np.arange(len(h[1]))[z >= zbin][0]+1
    
    sel = (h_idx == bin_id) & best
    
    fp = open('/tmp/view_list','w')
    for id in zc.id[sel]:
        fp.write('%s.new_zfit.png\n' %(id))
        fp.write('%s.new_zfit.2D.png\n' %(id))
        fp.write('%s_stack.png\n' %(id))
    
    fp.close()
    os.system('open `cat /tmp/view_list`')
    
    plt.scatter(c.x_image, c.y_image, color='black', alpha=0.1)
    plt.scatter(c.x_image[cok][sel], c.y_image[cok][sel], color='red', alpha=0.8)
    
    
    
    
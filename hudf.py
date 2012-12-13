import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import pyfits

import threedhst
import threedhst.catIO as catIO
import unicorn

def make_ellis12_regions():
    
    id = ['UDF12-3954-6284', 'UDF12-4106-7304', 'UDF12-4265-7049', 'UDF12-3921-6322', 'UDF12-4344-6547', 'UDF12-3895-7114', 'UDF12-3947-8076', 'UDFj-38116243', 'UDFj-43696407', 'UDFj-35427336', 'UDFy-38135539', 'UDFy-37796000', 'UDFy-33436598']

    ra = ['3:32:39.54', '3:32:41.06', '3:32:42.65', '3:32:39.21', '3:32:43.44', '3:32:38.95', '3:32:39.47', '3:32:38.11', '3:32:43.69', '3:32:35.42', '3:32:38.13', '3:32:37.79', '3:32:33.43']
    
    dec = ['-27:46:28.4', '-27:47:30.4', '-27:47:04.9', '-27:46:32.2', '-27:46:54.7', '-27:47:11.4', '-27:48:07.6', '-27:46:24.3', '-27:46:40.7', '-27:47:33.6', '-27:45:53.9', '-27:46:00.0', '-27:46:59.8']
    
    redshift = [11.9, 9.5, 9.5, 8.8, 8.8, 8.6, 8.6, 0, 7.6, 7.9, 8.3, 8.1, 7.9]
    
    fp = open('ellis12.reg','w')
    fp.write('fk5\n')
    for i in range(len(id)):
        fp.write('circle(%s, %s, 1") # text={%s}\n' %(ra[i], dec[i], id[i]))
    
    fp.close()
    
    fp = open('ellis12.dat','w')
    fp.write('# id    rh     dh   ra    dec     z\n')
    for i in range(len(id)):
        fp.write('%15s  %s  %s  %.5f  %.5f  %.1f\n' %(id[i], ra[i], dec[i], threedhst.utils.DMS2decimal(ra[i], hours=True), threedhst.utils.DMS2decimal(dec[i], hours=False), redshift[i]))
    fp.close()
    
    
def extract_hudf12():
    """
    Make catalog and segmentation image from the 
    F160W image released by Ellis et al.
    """
    import threedhst
    import numpy as np
    
    ## Make catalog
    os.chdir("/research/HST/GRISM/3DHST/UDF/HUDF12")
    
    header = pyfits.getheader('hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits')
    ZP = -2.5*np.log10(header['PHOTFLAM']) - 21.10 - 5*np.log10(header['PHOTPLAM']) + 18.6921
         
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = 'hudf12.cat'
    se.options['CHECKIMAGE_NAME'] = 'hudf12_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_wht.fits'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_MINAREA']    = '10'
    se.options['DETECT_THRESH']    = '1.0'
    se.options['ANALYSIS_THRESH']  = '1.0'
    se.options['MAG_ZEROPOINT'] = '%.5f' %(ZP)
    se.options['DEBLEND_NTHRESH'] = '64'
    se.options['DEBLEND_MINCONT'] = '0.00002'
    se.options['PHOT_APERTURES'] = '6.6667,8.333333,16.667,8.335,20'
    se.options['GAIN']  = '1.0'
    
    status = se.sextractImage('hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits')
    
    threedhst.sex.sexcatRegions('hudf12.cat', 'hudf12.reg', format=2)
    
    #### Find matches to Ellis positions
    cat = threedhst.sex.mySexCat('hudf12.cat')
    cat.write('hudf12.reform.cat', reformat_header=True)
    cat = catIO.Readfile('hudf12.reform.cat')
    matcher = catIO.CoordinateMatcher(cat, USE_WORLD=True)
    
    lines = open('../ellis12.dat').readlines()
    fp = open('ellis12.matched.dat','w')
    fp.write(lines[0][:-1]+' my_id\n')
    for line in lines[1:]:
        sp = line.split()
        ra, dec = float(sp[3]), float(sp[4])
        dr, id = matcher.find_nearest(ra, dec, N=2)
        fp.write(line[:-1]+' %6d\n' %(cat.number[id][0]))
        print sp[0], dr, id
    
    fp.close()
    
def prepare():
    
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
    
    NGROW=125
    ROOT = 'PRIMO-1101'
    ROOT = 'GOODS-SOUTH-34'
    for p in [34, 36, 37, 38]:
        ROOT = 'GOODS-SOUTH-%d' %(p)
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = ROOT+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=ROOT+'-F140W', view=False, pad=60+200*(ROOT=='PRIMO-1026'), REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, auto_offsets=True, NSEGPIX=3)
        seg = pyfits.open(ROOT+'_inter_seg.fits', mode='update')
        seg[0].data[seg[0].data < 0] = 0
        seg.flush()
        #
        unicorn.reduce.interlace_combine(root=ROOT+'-G141', view=False, pad=60+200*(ROOT=='PRIMO-1026'),  NGROW=NGROW, auto_offsets=True)
        unicorn.reduce.interlace_combine(root=ROOT+'-F140W', view=False, pad=60+200*(ROOT=='PRIMO-1026'),  NGROW=NGROW, auto_offsets=True)
        #
        ref = pyfits.open(ROOT+'_ref_inter.fits', mode='update')
        ref[0].header['FILTER'] = 'F160W'
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
    
    if ROOT.startswith('PRIMO'):
        ref = pyfits.open(ROOT+'_ref_inter.fits')
        im140 = pyfits.open(ROOT+'-F140W_inter.fits', mode='update')
        im140[0].header['FILTER'] = 'F160W'
        im140[1].data = ref[1].data #/ 10**(-0.4*(26.46-25.96))
        #im140[2].data = im140[2].data # / 10**(-0.4*(26.46-25.96))
        im140.flush()
    
    #
    files = glob.glob('*ref_inter.fits')
    for file in files:
        ROOT=file.split('_ref_inter')[0]
        model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=25, REFINE_MAG_LIMIT=23, make_zeroth_model=False, grism='G141')
        model.make_wcs_region_file()
        
    #model = unicorn.reduce.GrismModel(root=ROOT, MAG_LIMIT=30, grism='G141')
    ### Force use F160W as detection image
    #if ROOT.startswith('PRIMO'):
    #model.direct[1].data = model.im[1].data*1./10**(-0.4*(26.46-25.96))
    
    model.get_corrected_wcs()
    model.make_wcs_region_file()
    
    #for p in [34,36,37,38]:
    #    ROOT = 'GOODS-SOUTH-%d' %(p)

def extract_dropouts():
    """
    Extract all targets from Ellis et al.
    """
    import unicorn.hudf as hudf
    cat = catIO.Readfile('../HUDF12/ellis12.matched.dat')
    for id in cat.my_id[0:]:
        hudf.extract_all(id, fit=False, miny=-200)
        
def extract_all(id=6818, fit=False, miny=-200):
    """
    Extract spectra from 3D-HST + CANDELS pointings
    """
    files = glob.glob('*ref_inter.fits')
    for file in files:
        ROOT=file.split('_ref_inter')[0]
        model = unicorn.reduce.process_GrismModel(root=ROOT, MAG_LIMIT=30, grism='G141')
        #model.get_corrected_wcs()
        #
        #id=6818  ### UDF dropout
        #id=6021  ### OIII emitter
        #id=7328
        #id=5836 # weaker line
        if id not in model.objects:
            continue
        #
        model.twod_spectrum(id, verbose=True, CONTAMINATING_MAGLIMIT=28, miny=miny)
        model.show_2d(savePNG=True, verbose=True)
        oned = unicorn.reduce.Interlace1D(ROOT+'_%05d.1D.fits' %(id), PNG=True, flux_units=True)
        if fit:
            gris = unicorn.interlace_fit.GrismSpectrumFit(ROOT+'_%05d' %(id))
            gris.fit_in_steps()
        
    stack(id, dy=20, save=True, inverse=True)
    
#
def stack(id=6818, dy=20, save=True, inverse=False, scale=[0.5,0.01], fcontam=0.):
    """
    Stack all UDF spectra for a given object ID
    """
    import scipy.ndimage as nd
    import glob
    
    if plt.cm.get_cmap().name != 'gray':
        plt.gray()
        
    files=glob.glob('[PG][RO]*%05d.2D.fits' %(id))

    #dy = 20
    NX = 270
    flux = np.zeros((len(files), 2*dy, NX))
    err = np.zeros((len(files), 2*dy, NX))
    contam = np.zeros((len(files), 2*dy, NX))
    model = np.zeros((len(files), 2*dy, NX))
    
    for i, file in enumerate(files):
        twod = unicorn.reduce.Interlace2D(file)
        im = twod.im
        #im = pyfits.open(file)
        #p = plt.plot(im['WAVE'].data-1.4e4, alpha=0.5, marker='o')
        y0 = np.interp(1.4e4, im['WAVE'].data, im['YTRACE'].data)
        sh = im['SCI'].data.shape
        x0 = np.interp(1.4e4, im['WAVE'].data, np.arange(sh[1]))
        #
        y0i = int(np.round(y0))
        xarr = np.arange(sh[1]*1.)
        #
        dx = int(np.round(149.68-x0))
        print y0, sh, x0, x0-149.68, dx
        if sh[1] < NX:
            err[i,:,:] = 10000.
            continue
        #
        flux[i,:,:] = nd.shift(im['SCI'].data[y0i-dy:y0i+dy,:], (0,dx), order=3)[:,:NX]
        contam[i,:,:] = nd.shift(im['CONTAM'].data[y0i-dy:y0i+dy,:], (0,dx), order=3)[:,:NX]
        err[i,:,:] = nd.shift(im['WHT'].data[y0i-dy:y0i+dy,:], (0,dx), order=3)[:,:NX]
        model[i,:,:] = nd.shift(im['MODEL'].data[y0i-dy:y0i+dy,:], (0,dx), order=3)[:,:NX]
        # p = plt.plot(nd.shift(xarr, -(149.68-x0), order=3), im['WAVE'].data-1.4e4, alpha=0.5, marker='o')
        wave = nd.shift(im['WAVE'].data, dx)[:NX]
        #p = plt.plot(wave-1.4e4, marker='o')
        #plt.xlim(148,152); plt.ylim(-40,40)

    NF = len(files)
    
    vmin, vmax = -1.e-4, 5.e-3
    if inverse:
        vmin, vmax = -vmax, -vmin
        flux = 0-flux
        contam = 0-contam
        model = 0-model
        
    weight = 1./((fcontam*contam)**2+err**2)
    weight[err <= 0] = 0
    weight[err < 1.e-6] = 0
    
    sum_weight = weight.sum(axis=0)
    sum_weight[sum_weight == 0] = 1.
    
    sum_flux = (flux*weight).sum(axis=0)/sum_weight
    sum_var = np.sqrt((err**2).sum(axis=0))
    sum_contam = (contam*weight).sum(axis=0)/sum_weight
    sum_model = (model*weight).sum(axis=0)/sum_weight
    
    diff = sum_flux-sum_contam
    if len(scale) == 0:
        scale = [scale, scale]
    vmin, vmax = np.percentile(diff.flatten(), scale[0]), np.percentile(diff.flatten(), 100-scale[1])
    if inverse:
        vmin, vmax = -vmax, -vmin
        
    #### Make the plot
    sh = sum_flux.shape
    aspect = (NF+2)*sh[0]*1./(3*sh[1])
    fig = unicorn.plotting.plot_init(square=True, left=0.01, right=0.01, top=0.01, bottom=0.01, hspace=0.0, wspace=0.0, aspect=aspect, xs=10, NO_GUI=False)
    
    counter = 4
    for i in range(NF):
        ax = fig.add_subplot(NF+2,3,counter); counter += 1
        a = ax.imshow(flux[i,:,:], interpolation='nearest', vmin=vmin, vmax=vmax)
        a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
        if i == 0: ax.set_title('Flux')
        #
        ax = fig.add_subplot(NF+2,3,counter); counter += 1
        a = ax.imshow(contam[i,:,:], interpolation='nearest', vmin=vmin, vmax=vmax)
        a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
        ax.text(0.95, 0.95, files[i].split('_')[0], ha='right', va='top', transform=ax.transAxes, color='white', size=8)
        if i == 0: 
            ax.set_title('Contam')
            ax.text(0,1.5,'UDF #%d' %(id), ha='left', va='bottom', transform=ax.transAxes, color='black', size=10)
        #
        ax = fig.add_subplot(NF+2,3,counter); counter += 1
        a = ax.imshow((flux-contam)[i,:,:], interpolation='nearest', vmin=vmin, vmax=vmax)
        a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
        if i == 0: ax.set_title(r'Flux $-$ Contam')
    
    #### Sum
    ax = fig.add_subplot(NF+2,3,counter); counter += 1
    a = ax.imshow(sum_flux, interpolation='nearest', vmin=vmin, vmax=vmax)
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    #
    ax = fig.add_subplot(NF+2,3,counter); counter += 1
    a = ax.imshow(sum_contam, interpolation='nearest', vmin=vmin, vmax=vmax)
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    ax.text(0.95, 0.95, 'Stack', ha='right', va='top', transform=ax.transAxes, color='white', size=8)
    #
    ax = fig.add_subplot(NF+2,3,counter); counter += 1
    a = ax.imshow((sum_flux-sum_contam), interpolation='nearest', vmin=vmin, vmax=vmax)
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    
    if save:
        outfile='UDF_%05d_stack.png' %(id)
        print outfile
        fig.savefig(outfile)
        plt.close()
        os.system('open %s' %(outfile))
    
    #
    fig = unicorn.plotting.plot_init(square=True, left=0.1, right=0.01, top=0.01, bottom=0.1, hspace=0.0, wspace=0.0, aspect=1, xs=5, NO_GUI=False)
    ax = fig.add_subplot(111)
    for i, file in enumerate(files):
        twod = unicorn.reduce.Interlace2D(file)
        im = twod.im
        sh = im['SCI'].data.shape
        x0 = int(np.round(np.interp(1.4e4, im['WAVE'].data, np.arange(sh[1]))))
        y0 = np.interp(1.4e4, im['WAVE'].data, im['YTRACE'].data)
        y0i = int(np.round(y0))
        profile = np.sum((1-2*inverse)*model[i,:,x0-10:x0+10], axis=1)
        #print y0i, x0, dy, model.shape, model[i,:,x0-10:x0+10].shape, profile.shape
        ax.plot(np.arange(2*dy), profile/np.sum(profile), label=file.split('2D')[0])
    
    ax.legend(prop=dict(size=8))
    
    fig.savefig('profile.pdf')
    plt.close()
    
    #
    udf = UDF(id=id, NPIX=dy, fcontam=fcontam)
    udf.go_stacks()
    udf.fix_oned_arrays()
    
    ### Ideas:
    
    ### Take a given f_lambda spectrum with emission lines and stack the
    ### computed models as with the object to put limit on emission line flux
    
    ### Should be able to sum spectra for any object like this....

class UDF():
    def __init__(self, id=6818, NPIX=20, verbose=True, fcontam=0.):
        """
        Read in the 2D FITS files
        """        
        self.verbose = verbose
        
        self.id = id
        self.files=glob.glob('[GP][OR]*%05d.2D.fits' %(id))
        self.NPIX = NPIX
        self.fcontam = fcontam
        
        self.N = len(self.files)
        self.nx = np.ones(self.N, dtype=int)
        self.ny = np.ones(self.N, dtype=int)
        self.y0 = np.ones(self.N)
        self.dx = np.ones(self.N)
        self.dxi = np.ones(self.N, dtype=int)
        ### Read in the spectra
        self.read_data()
        
    def go_stacks(self):
        ### Stack the observed spectra
        self.stack_spectra()
        
        ### Make the FITS files
        self.make_stacked_FITS(root='UDF')
        
        
    def read_data(self):
        """
        Read in 2D fits files and get information.
        """
        import scipy.ndimage as nd
    
        NPIX = self.NPIX
        
        self.twod = []
        for i in range(self.N):
            if self.verbose:
                print self.files[i]
            ### Size of FITS arrays
            self.twod.append(unicorn.reduce.Interlace2D(self.files[i]))
            im = self.twod[i].im
            self.ny[i], self.nx[i] = im['SCI'].data.shape
            ### Ycenter of trace
            self.y0[i] = np.interp(1.4e4, im['WAVE'].data, im['YTRACE'].data)
            ### Xoffset to register wavelength
            x0 = np.interp(1.4e4, im['WAVE'].data, np.arange(self.nx[i]))
            self.dx[i] = 150-x0
                
        #### Now store shifted flux arrays
        self.dx -= self.dx[-1]
        self.dxi[i] = np.cast[int](np.round(self.dx[i]))
        
        NX = self.nx.min()
        self.flux = np.zeros((self.N, 2*NPIX, NX))
        self.err = np.zeros((self.N, 2*NPIX, NX))
        self.contam = np.zeros((self.N, 2*NPIX, NX))
        self.model = np.zeros((self.N, 2*NPIX, NX))
        
        self.dsci = np.zeros((2*NPIX,2*NPIX))
        self.dseg = np.zeros((2*NPIX,2*NPIX))
        
        for i in range(self.N):
            y0i = int(np.round(self.y0[i]))
            dx = self.dxi[i]
            im = self.twod[i].im
            #
            self.flux[i,:,:] = nd.shift(im['SCI'].data[y0i-NPIX:y0i+NPIX,:], (0,dx), order=3)[:,:NX]
            self.contam[i,:,:] = nd.shift(im['CONTAM'].data[y0i-NPIX:y0i+NPIX,:], (0,dx), order=3)[:,:NX]
            self.err[i,:,:] = nd.shift(im['WHT'].data[y0i-NPIX:y0i+NPIX,:], (0,dx), order=3)[:,:NX]
            self.model[i,:,:] = nd.shift(im['MODEL'].data[y0i-NPIX:y0i+NPIX,:], (0,dx), order=3)[:,:NX]
            #
            self.dsci += im['DSCI'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX]/self.N
            self.dwht = im['DWHT'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX]
            self.dseg = ((self.dseg == self.id) | (im['DSEG'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX] == self.id))*int(self.id)
            
            # p = plt.plot(nd.shift(xarr, -(149.68-x0), order=3), im['WAVE'].data-1.4e4, alpha=0.5, marker='o')
        
        #print 'DX: ', dx
        self.wave = nd.shift(im['WAVE'].data, dx)[:NX]
        ### Fix case when you have zeros in the wavelength grid after shifting
        if (self.wave == 0).sum() > 0:
            i0 = np.where(self.wave > 0)[0][0]
            delta = np.median(np.diff(self.wave))
            self.wave = np.arange(NX)*delta+self.wave[i0]-(i0-1)*delta
        
        self.ytrace = nd.shift(im['YTRACE'].data, dx)[:NX]-y0i+NPIX
        self.sens = nd.shift(im['SENS'].data, dx)[:NX]

        ### Dummy 1D info
        # oned = self.twod[-1].oned
        # self.x_onedflux = nd.shift(oned.flux, dx)[:NX]
        # self.x_onedcont = nd.shift(oned.contam, dx)[:NX]

        ### Thumbnail
        self.dsci = im['DSCI'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX]
        self.dwht = im['DWHT'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX]
        self.dseg = im['DSEG'].data[y0i-NPIX:y0i+NPIX, y0i-NPIX:y0i+NPIX]

    def stack_spectra(self):
        """
        Stack the flux/contam/weight arrays of the available UDF spectra with
        inverse variance weighting.
        """
        weight = 1./((self.fcontam*self.contam)**2+self.err**2)
        weight[self.err <= 0] = 0
        weight[self.err < 1.e-6] = 0
        
        sum_weight = weight.sum(axis=0)
        sum_weight[sum_weight == 0] = 1.

        self.sum_flux = (self.flux*weight).sum(axis=0)/sum_weight
        self.sum_contam = (self.contam*weight).sum(axis=0)/sum_weight
        self.sum_model = (self.model*weight).sum(axis=0)/sum_weight
        self.sum_var = 1./sum_weight
        
    def make_stacked_FITS(self, root='UDF'):
        """
        Make a dummy 2D FITS file
        """
        import copy
        im = pyfits.open(self.files[-1])
        
        im['DSCI'].data = self.dsci
        im['DWHT'].data = self.dwht
        im['DSEG'].data = self.dseg
        
        im['SCI'].data = self.sum_flux
        im['WHT'].data = np.sqrt(self.sum_var)
        im['MODEL'].data = self.sum_model
        im['CONTAM'].data = self.sum_contam
        im['WAVE'].data = self.wave
        im['SENS'].data = self.sens
        im['YTRACE'].data = self.ytrace
        
        im.writeto('%s_%05d.2D.fits' %(root, self.id), clobber=True)
        
        #### 1D FITS file
        self.optimal_flux = self.wave*0.+1
        self.optimal_contam = self.wave*0.
        self.optimal_error = self.optimal_flux/10.
        self.trace_flux = self.wave*0.+1
        self.trace_err = self.wave*0.+1
        
        self.make_1D_spectrum()
        
    def make_1D_spectrum(self, root='UDF'):
        
        oned = pyfits.open(self.files[-1].replace('2D','1D'))
        
        c1 = pyfits.Column(name='wave', format='D', unit='ANGSTROMS', array=self.wave)
        c2 = pyfits.Column(name='flux', format='D', unit='ELECTRONS/S', array=self.optimal_flux)
        c3 = pyfits.Column(name='error', format='D', unit='ELECTRONS/S', array=self.optimal_error)
        c4 = pyfits.Column(name='contam', format='D', unit='ELECTRONS/S', array=self.optimal_contam)
        c5 = pyfits.Column(name='trace', format='D', unit='ELECTRONS/S', array=self.trace_flux)
        c6 = pyfits.Column(name='etrace', format='D', unit='ELECTRONS/S', array=self.trace_err)
        c7 = pyfits.Column(name='sensitivity', format='D', unit='E/S / 1E-17 CGS', array=self.sens*(self.twod[0].growx*self.twod[0].growy))
        #print 'MAX SENS: %.f' %(self.sens.max())
        
        coldefs = pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7])
        head = oned[1].header

        tbHDU = pyfits.new_table(coldefs, header=head)
        tbHDU.writeto('%s_%05d.1D.fits' %(root, self.id), clobber=True)
        
    def fix_oned_arrays(self, root='UDF'):
        """
        Optimal extraction for 1D spectra
        """
        two = unicorn.reduce.Interlace2D('%s_%05d.2D.fits' %(root, self.id))
        x0, self.optimal_contam = two.optimal_extract(two.im['CONTAM'].data)
        x0, self.optimal_flux = two.optimal_extract(two.im['SCI'].data)
        self.trace_flux = self.optimal_flux
        self.optimal_error = self.optimal_flux/10.
        self.make_1D_spectrum(root=root)
        
    def stack_model(self, lam_spec=None, flux_spec=None):
        """
        Take an input model (wavleength, f_lambda) and make
        a stacked 2D model.
        """
        import scipy.ndimage as nd

        NPIX = self.NPIX
        NX = self.nx.min()
        self.twod_model = np.zeros((self.N, 2*NPIX, NX))
        
        for i in range(self.N):
            self.twod[i].compute_model(lam_spec, flux_spec)
            y0i = int(np.round(self.y0[i]))
            dx = self.dxi[i]
            self.twod_model[i,:,:] = nd.shift(self.twod[i].model[y0i-NPIX:y0i+NPIX,:], (0,dx), order=3)[:,:NX]

        #
        weight = 1./self.err**2
        weight[self.err <= 0] = 0
        weight[self.err < 1.e-6] = 0
        
        sum_weight = weight.sum(axis=0)
        sum_weight[sum_weight == 0] = 1.

        self.sum_twod_model = (self.twod_model*weight).sum(axis=0)/sum_weight
        
    def show_stack(self):
        """
        Make diagnostic plot(s)
        """
        pass

def check_for_ivo():
    """
    Ivo has what looks like a massive z>4 galaxy that happens to be in the
    UDF.  Extract a stacked spectrum for him.
    """
    id_3dhst = 19356
    cat2.read_catalogs('GOODS-S')
    
    cat = catIO.Readfile('../HUDF12/hudf12.reform.cat')
    matcher = catIO.CoordinateMatcher(cat, USE_WORLD=True)
    ix = cat2.cat.id == id_3dhst
    dr, id = matcher.find_nearest(cat2.cat.ra[ix][0], cat2.cat.dec[ix][0], N=1)
    
    ###
    udf12_id = cat.number[id]
    hudf.extract_all(udf12_id, fit=False)
    
def check_lehnert():
    """
    Check for the z=8.6 line observed by Lehnert et al. on UDFy-38135539
    """
    id=7507
    
    twod = unicorn.reduce.Interlace2D('UDF_%05d.2D.fits' %(id))
    
    twod = unicorn.reduce.Interlace2D('PRIMO-1101_%05d.2D.fits' %(id))
    
    #### Make simple segmentation aperture with radius 4 pix, D~0.5"
    thumb = twod.im['DSCI'].data
    yi, xi = np.indices(thumb.shape)
    r = np.sqrt((xi-thumb.shape[1]/2.)**2+(yi-thumb.shape[0]/2.)**2)
    aper = r < 5
    #### Normalize the flux in the direct image
    twod.im['DSEG'].data = aper*id
    twod.im['DSCI'].data = twod.im['DSCI'].data / np.sum(twod.im['DSCI'].data*aper)
    #twod.flux *= aper
    twod.total_flux = np.sum(twod.flux*aper)
    twod.seg = twod.im['DSEG'].data*1
    twod.init_model()
    thumb = twod.flux*aper/np.sum(twod.flux*aper)
    
    #### Model the emission line, 1D + 2D models
    ztry = 8.5549
    
    ### Continuum from BC03
    lc, fc = np.loadtxt('/Users/gbrammer/research/drg/TEMPLATES/BC03_FROM_FAST/bc03_pr_ch_z02_ltau07.0_age08.0.dat', unpack=True)
    lcz = lc*(1+ztry)
    igm = threedhst.eazyPy.igm_factor(lcz, ztry)
    fc = fc*igm
    from scipy import interpolate
    interp_continuum = interpolate.InterpolatedUnivariateSpline(lcz,fc)
    
    ### Grid for G141 spectrum
    lx = np.arange(0.8e4,2.0e4,0.2)
    
    ### Continuum, normalized to F160W mag
    xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), 'F160W'), unpack=True)
    yf_int = np.interp(lcz, xf, yf, left=0., right=0.)
    
    fnu = fc*lcz**2/3.e18
    flux = np.trapz(fnu*yf_int, 1./lcz) / np.trapz(yf_int, 1./lcz)
    f160w_mag = 28.4
    f160w_fnu = 10**(-0.4*(f160w_mag+48.6))
    f160w_scale = f160w_fnu / flux
    
    continuum_flam = interp_continuum(lx)*f160w_scale
    
    ### Lehnert emission line
    fwhm = 9.
    sigma = fwhm / 2.35
    lam = 1216.*(1+ztry) #1.1615e4
    total_flux = 6e-18
    gauss = (1./np.sqrt(2*np.pi*sigma**2)*np.exp(-(lx-lam)**2/2/sigma**2)+1.e-5)*total_flux
    
    ### Check Ellis mags
    for f, ref in zip(['F125W','F140W', 'F160W'], [28.6, 28.5, 28.4]):
        xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), f), unpack=True)
        yf_int = np.interp(lx, xf, yf, left=0., right=0.)
        fnu = (continuum_flam + gauss)*lx**2/3.e18
        fnu2 = (continuum_flam + gauss*0)*lx**2/3.e18
        flux = np.trapz(fnu*yf_int, 1./lx) / np.trapz(yf_int, 1./lx)
        flux2 = np.trapz(fnu2*yf_int, 1./lx) / np.trapz(yf_int, 1./lx)
        mag = -2.5*np.log10(flux)-48.6
        print '%s: %.1f - %.3f  %.3f' %(f, ref, -2.5*np.log10(flux)-48.6, -2.5*np.log10(flux2)-48.6)
        
    #### Model 2D spectrum
    twod.compute_model(lx, (continuum_flam + gauss)/1.e-17/twod.total_flux)
    
    # new = unicorn.reduce.Interlace2D('PRIMO-1101_%05d.2D.fits' %(id))
    # new.im['SCI'].data += twod.model*6 - np.median(new.im['SCI'].data-new.im['CONTAM'].data) 
    # new.im.writeto('junk_07507.2D.fits', clobber=True)
    # os.system('cp PRIMO-1101_07507.1D.fits junk_07507.1D.fits')
    
    #### Test integrate model flux in 2d spectrum to compare to input
    sens = np.interp(lam, twod.oned.lam, twod.oned.sens)
    model_flux = np.sum(twod.model[:,43:54])/sens*1.e-17*np.median(np.diff(twod.im['WAVE'].data))
    
    #### Convolution
    NY, NX = twod.im['SCI'].data.shape
    flux = np.zeros((3*NY, NX+2*NY))
    contam = flux*0.
    model = flux*0.
    
    flux[NY:-NY, NY:-NY] = twod.im['SCI'].data
    contam[NY:-NY, NY:-NY] = twod.im['CONTAM'].data
    model[NY:-NY, NY:-NY] = twod.model
    
    import stsci.convolve.Convolve as sc
    
    cd = sc.correlate2d(flux-contam, thumb, mode='constant', cval=0.)
    cf = sc.correlate2d(flux, thumb, mode='constant', cval=0.)
    cc = sc.correlate2d(contam, thumb, mode='constant', cval=0.)
    cm = sc.correlate2d(model, thumb, mode='constant', cval=0.)
    ctotal = sc.correlate2d(flux-contam+model, thumb, mode='constant', cval=0.)
    
    ################# Comparison object for fluxes
    id=6752
    twod = unicorn.reduce.Interlace2D('PRIMO-1101_%05d.2D.fits' %(id))
    yi, xi = np.indices(twod.im['SCI'].data.shape)
    r = np.sqrt((xi-110)**2+(yi-26)**2)
    aper = r < 8
    
    r2 = np.sqrt((xi-130)**2+(yi-26)**2)
    aper2 = r2 < 8
    
    sens = np.interp(1.31e4, twod.oned.lam, twod.oned.sens)
    diff = twod.im['SCI'].data-twod.im['CONTAM'].data
    model_flux = np.sum(diff*aper-diff*aper2)/sens*1.e-17*np.median(np.diff(twod.im['WAVE'].data))
        
    lx = np.arange(0.9e4,1.8e4,0.2)
    fwhm = 9.
    sigma = fwhm / 2.35
    lam = 1.31e4
    total_flux = 4.5e-17
    gauss = (1./np.sqrt(2*np.pi*sigma**2)*np.exp(-(lx-lam)**2/2/sigma**2)+1.e-6)*total_flux
    twod.compute_model(lx, gauss/1.e-17/twod.total_flux)
    model_flux = np.sum(twod.model)/sens*1.e-17*np.median(np.diff(twod.im['WAVE'].data))

    
    
def check_z12_source(scale = 1, z=2.25, own_spec=1, observed_mag=29.3, save=False):
    import unicorn.hudf
    
    
    id=6001
    
    udf = unicorn.hudf.UDF(id=id, NPIX=20, fcontam=1)
    udf.go_stacks()
    udf.fix_oned_arrays()
    
    #unicorn.hudf.fix_2d_background(object = 'UDF_%05d' %(id), force=False)  
    #udf = unicorn.hudf.UDF(id=id, NPIX=100, fcontam=1)
    #udf.fix_oned_arrays()
    
    twod = unicorn.reduce.Interlace2D('UDF_%05d.2D.fits' %(id))
    # twod.flux /= twod.total_flux
    # twod.total_flux = 1.
    # twod.init_model()
    
    lam, flam, fnu = unicorn.hudf.generate_line_model(z=z, filter='F160W', observed_mag=observed_mag, own_spec=own_spec)
    twod.compute_model(lam, flam/1.e-17*scale/twod.total_flux)
    check = unicorn.hudf.check_mags(lam, fnu*scale)
    # 
    # ds9.frame(1)
    # ds9.view(twod.im['SCI'].data-twod.im['CONTAM'].data)
    # ds9.scale(-0.01, 0.01)
    # ds9.frame(2)
    # ds9.view(twod.model)
    # ds9.frame(3)
    # ds9.view(twod.im['SCI'].data-twod.im['CONTAM'].data+twod.model)
    # 
    #### Optimal extraction 
    w0, x0 = twod.optimal_extract(input=twod.im['SCI'].data-twod.im['CONTAM'].data+twod.model*0)
    w1, x1 = twod.optimal_extract(input=twod.im['SCI'].data-twod.im['CONTAM'].data+twod.model*1)
    wl, xl = twod.optimal_extract(twod.model*1)
    sens = twod.im['SENS'].data
    
    #### Make plot
    fig = unicorn.plotting.plot_init(square=True, left=0.11, top=0.01, right=0.01, bottom=0.07, aspect=6*twod.model.shape[0]*1./twod.model.shape[1], xs=6)
    
    ax = fig.add_subplot(611)
    a = ax.imshow(twod.im['SCI'].data-twod.im['CONTAM'].data, vmin=-0.005, vmax=0.005, aspect='auto', interpolation='nearest')
    xa = np.arange(twod.im['WAVE'].data.shape[0])
    ta = np.interp(np.arange(1.1e4,1.61e4,0.1e4), twod.im['WAVE'].data, xa)
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    a = ax.set_xticks(ta)
    a = ax.set_ylabel(r'Flux$-$Cont.')
    
    ax = fig.add_subplot(612)
    a = ax.imshow(twod.model, vmin=-0.005, vmax=0.005, aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    a = ax.set_xticks(ta)
    a = ax.set_ylabel(r'Model')
    
    ax = fig.add_subplot(613)
    a = ax.imshow(twod.im['SCI'].data-twod.im['CONTAM'].data+twod.model, vmin=-0.005, vmax=0.005, aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([])
    a = ax.set_xticks(ta)
    a = ax.set_ylabel(r'F$-$C+M')
    
    ax = fig.add_subplot(614)
    ellis = dict(F105W=-31.2, F140W=-30.5, F125W=-30.7, F160W=29.3)
    colors = dict(F105W='blue',F125W='green',F140W='yellow',F160W='orange')
    for f in ellis.keys():
        if ellis[f] < 0:
            m = 'v'
        else:
            m = 'o'
        #
        print f, check[f]
        if (check[f] > 34) | (not np.isfinite(check[f])):
            check[f] = 34.
            mc = 'v'
        else:
            mc = 'o'
        a = ax.plot(unicorn.reduce.PLAMs[f], np.abs(ellis[f]), marker=m, color=colors[f], markersize=10)
        a = ax.plot(unicorn.reduce.PLAMs[f], np.abs(check[f]), marker=mc, color='red', markersize=10, alpha=0.5)
    #
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[1], ylim[0])
    for f in ellis.keys():
        xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), f), unpack=True)
        r = np.random.rand()*0.2
        a = ax.plot(xf, ylim[1]- yf/yf.max()*(0.4+r)*(ylim[1]-ylim[0]), color='black', alpha=0.3)
        a = ax.fill_between(xf, ylim[1]- yf/yf.max()*(0.4+r)*(ylim[1]-ylim[0]), yf*0.+ylim[1], color=colors[f], alpha=0.3)
    #
    spec_ab = -2.5*np.log10(fnu*scale)-48.6
    ok = (lam > w0[0]) & (lam < w0[-1]) & (lam/(1+z) > 1210)
    #return lam, spec_ab, w0
    
    a = ax.plot(lam, spec_ab, color='red', alpha=0.1)
    a = ax.plot(lam[ok], spec_ab[ok], color='red')
    
    a = ax.set_xlim(w0[0], w0[-1])
    a = ax.set_xticklabels([])
    a = ax.set_ylabel(r'AB Mag.')
    
    ax = fig.add_subplot(313)
    
    a = ax.plot(w0, x0/sens*100, color='black', alpha=0.5)
    #a = ax.plot(w1, x1/sens*100, color='red', alpha=0.8)
    a = ax.fill_between(wl, xl*100-4*0, xl*0-4*0, color='red', alpha=0.4)
    a = ax.plot(lam, flam/1.e-19*scale-4, color='orange', alpha=0.8)
    
    a = ax.set_xlim(w0[0], w0[-1])
    a = ax.set_ylim(-4, 3.9)
    
    a = ax.text(0.5, 0.98, r'$z=%.3f$' %(z), ha='right', va='top', size=12, transform=ax.transAxes)
    a = ax.text(0.7, 0.98, r'$H_{160}=%.2f$' %(check['F160W']), ha='right', va='top', size=12, transform=ax.transAxes)
    continuum = np.interp(1.68e4, lam, flam)
    EW = np.trapz(((flam-continuum)/continuum)[ok], lam[ok])/(1+z)
    a = ax.text(0.95, 0.98, r'$EW_\mathrm{rest}=%d\,\AA$' %(EW), ha='right', va='top', size=12, transform=ax.transAxes)
    
    a = ax.set_xlabel(r'$\lambda$'); a = ax.set_ylabel(r'$f_\lambda\, (10^{-19}\,\mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}\,\mathrm{\AA}^{-1})$')
   
    if save is False:
        savefile='junk.png'
    else:
        savefile=save
    
    fig.savefig(savefile)
    plt.close()
    
    # ds9.frame(4)
    # rnd = np.random.normal(size=twod.model.shape)*twod.im['WHT'].data
    # ds9.view(rnd)
    
#
def go_check():
    comment = """
    Halpha and OIII+ HB at z=1.45 to match the potential observed lines.
    Need to bump up F160W flux to match the Halpha flux and then 
    get way to much blue flux from OIII for the bluer bands
    """
    hudf.check_z12_source(scale=1, z=1.44, own_spec=500, save='z1.4_29.3.pdf')
    
    comment = """
    Maybe a weak OIII line at z=2.195, but would make a F160W flux 3-4 times
    higher than observed (and a F140W detection with Hb)
    """
    hudf.check_z12_source(scale=1, z=2.20, own_spec=0, save='lens_spectrum_z2.2_29.3.pdf')
    hudf.check_z12_source(scale=3, z=2.20, own_spec=0, save='lens_spectrum_z2.2_29.3x3.pdf')

    hudf.check_z12_source(scale=1, z=2.20, own_spec=2000, observed_mag=29.3)
    hudf.check_z12_source(scale=1, z=2.20, own_spec=2000, observed_mag=28.7)

    comment = """
    OII, could maybe have an OIIx3727 line at z=3.3, 30% stronger than nominal
    to just make the upper limit of the F160W detection.  Would probably need
    pretty extreme equivalent width
    """
    hudf.check_z12_source(scale=1, z=3.30, own_spec=500, save='OII_z3.3_29.3.pdf')
    hudf.check_z12_source(scale=1/0.5, z=3.30, own_spec=100, save='OII_z3.3_29.3x2.pdf')
    
    comment = """
    Ly-alpha @ z=12.17???
    """
    hudf.check_z12_source(scale=1, z=12.17, own_spec=1000, save='z12_29.3.pdf')
    hudf.check_z12_source(scale=1, z=12.17, own_spec=1000, observed_mag=28.7, save='z12_28.7.pdf')
    
    #gris = unicorn.interlace_fit.GrismSpectrumFit('UDF_06001')
    #gris.fit_free_emlines(ztry=12.17)
    #gris.fit_free_emlines(ztry=1.44)
    
def generate_line_model(z=2.1, filter='F160W', observed_mag=29.3, own_spec=True):
    """
    Scale the UDS lens spectrum to a given observed magnitude and redshift
    
    ellis_m160 = 29.3
    ellis_m140 = 30.5
    ellis_m125 = 30.7
    ellis_m105 = 31.2
    
    """
    from scipy import interpolate
    
    #xf140, yf140 = np.loadtxt(os.getenv('iref')+'/F140W.dat', unpack=True)
    #xf160, yf160 = np.loadtxt(os.getenv('iref')+'/F160W.dat', unpack=True)
    #xf125, yf125 = np.loadtxt(os.getenv('iref')+'/F125W.dat', unpack=True)
    #xf105, yf105 = np.loadtxt(os.getenv('iref')+'/F105W.dat', unpack=True)
    xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), filter), unpack=True)

    #### Use 3D-HST lens as input
    lensx, lensy = np.loadtxt('lens_template.dat', unpack=True)

    #### Single set of lines
    if own_spec:
        waves = [1215.24, 2799.12, 3728., 4862.68, 4960.3, 5008.24, 6564.61]
        waves = [1216, 2799.12, 3728., 4862.68, 4960.3, 5008.24, 6564.61]
        fluxes = [1, 4, 4, 0.0, 1, 3., 0.4]
        dv = 60.
        lensx = np.arange(712.,1.e4,0.01)
        lensy = lensx*0.+1.e-6*own_spec
        for i in range(len(waves)):
            dlam = dv/3.e5*waves[i]
            lensy += 1./np.sqrt(2*np.pi*dlam**2)*fluxes[i]*np.exp(-(lensx-waves[i])**2/2/dlam**2)
    
    lensy[lensx < 1210] = 1.e-10
    lensy_fnu = lensy*lensx**2
    
    ### Range of filters    
    x0 = xf[yf > 0.05*yf.max()][[0,-1]]
    if ((x0[0]/5007.-1) > z) | ((x0[1]/5007.-1) < z):
        print 'OIII not in filter %s at z=%.3f' %(filter, z)
        
    #### Interpolate the filter curve and compute integrated flux
    s = interpolate.InterpolatedUnivariateSpline(xf,yf)
    lensx_z = (1+z)*lensx
    yf_int = s(lensx_z); yf_int[(lensx_z <= xf.min()) | (lensx_z >= xf.max())] = 0.
    
    filter_flux = np.trapz(lensy_fnu*yf_int, 3.e8/(lensx_z*1.e-10))/np.trapz(yf_int, 3.e8/(lensx_z*1.e-10))
    filter_mag = -2.5*np.log10(filter_flux)-48.6
    
    observed_fnu = 10**(-0.4*(observed_mag+48.6))
    scale_fnu = 10**(-0.4*(observed_mag-filter_mag))
    
    lensy_flam = lensy_fnu*scale_fnu*3.e18/lensx_z**2
    if observed_mag < 0:
        flam_flux = np.trapz(lensy_flam*yf_int, lensx_z)/np.trapz(yf_int, lensx_z)
        lensy_flam /= flam_flux
        scale_fnu /= flam_flux
        
    
    return lensx_z, lensy_flam, lensy_fnu*scale_fnu
    
def check_mags(wave, fnu):
    """
    Integrate a spectrum through WFC3/IR filters to get broad-band mags
    """
    from scipy import interpolate
    out = {}
    for filter in ['F105W','F125W','F140W','F160W']:
        xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), filter), unpack=True)
        s = interpolate.InterpolatedUnivariateSpline(xf,yf)
        yf_int = s(wave); yf_int[(wave <= xf.min()) | (wave >= xf.max())] = 0.
        filter_flux = np.trapz(fnu*yf_int, 3.e8/(wave*1.e-10))/np.trapz(yf_int, 3.e8/(wave*1.e-10))
        filter_mag = -2.5*np.log10(filter_flux)-48.6
        out[filter] = filter_mag
        
    return out
    
def rgb_browser():
    ### HJi
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(25.94-25.96))*1.5]
    
    ### HJY
    scales = [10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(26.2687-25.96))]
        
    rgb1 = 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f125w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f105w_v1.0_drz.fits[0]*%.3f' %(scales[0]*6, scales[1]*6*1.1, scales[2]*6*1.1)
    rgb2 = 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f125w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f105w_v1.0_drz.fits[0]*%.3f' %(scales[0]*20, scales[1]*20*1.1, scales[2]*20*1.1)
    rgb3 = 'hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f125w_v1.0_drz.fits[0]*%.3f, hlsp_hudf12_hst_wfc3ir_udfmain_f105w_v1.0_drz.fits[0]*%.3f' %(scales[0]*50, scales[1]*50*1.1, scales[2]*50*1.1)

    threedhst.gmap.makeImageMap([rgb3], aper_list=[15], tileroot=['HJY'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)

    threedhst.gmap.makeImageMap([rgb1,rgb2,rgb3], aper_list=[15], tileroot=['HJY','deep','deeper'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)

    threedhst.gmap.makeImageMap([rgb1,rgb2], aper_list=[15,16,17], tileroot=['HJY','deep'], extension=0, path='./HTML/', zmin=-0.05, zmax=1)
    
    
    #### Full PNG
    imr = pyfits.open('hlsp_hudf12_hst_wfc3ir_udfmain_f160w_v1.0_drz.fits')
    img = pyfits.open('hlsp_hudf12_hst_wfc3ir_udfmain_f125w_v1.0_drz.fits')
    imb = pyfits.open('hlsp_hudf12_hst_wfc3ir_udfmain_f105w_v1.0_drz.fits')
    
    f = 10
    scales = np.array([10**(-0.4*(25.96-25.96)), 10**(-0.4*(26.25-25.96)), 10**(-0.4*(26.2687-25.96))])*f
    #
    xc, yc, NX = 1800, 1800, 200
    xc, yc = 1570, 2358
    # sub_r = imr[0].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[0]
    # sub_g = img[0].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[1]*1.0
    # sub_b = imb[0].data[yc-NX:yc+NX, xc-NX:xc+NX]*scales[2]*1.0
    #
    sub_r = imr[0].data*scales[0]
    sub_g = img[0].data*scales[1]*1.0
    sub_b = imb[0].data*scales[2]*1.0
    #
    Q, alpha, m0 = 3., 8, -0.02
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='HUDF12_wfc3_%05.1f.png' %(f), shape=sub_r.shape)
    
def xcor_analysis():
    """
    Run a cross-correlation analysis to try to improve detection of the 
    z=12 line
    """
    
    id = 6001
    
    ### Read the spectrum
    ## Fix background
    # unicorn.hudf.fix_2d_background(object = 'UDF_%05d' %(id), force=False)
    twod = unicorn.reduce.Interlace2D('UDF_%05d.2D.fits' %(id))

    #### Xcor kernels
    thumb = twod.im['DSCI'].data*(twod.im['DSEG'].data == id)
    my_f160w_mag = -2.5*np.log10(thumb.sum())+unicorn.reduce.ZPs['F160W']
    kernel = thumb[22:-22,22:-22]
    my_f160w_mag = -2.5*np.log10(kernel.sum())+unicorn.reduce.ZPs['F160W']
    
    #### Get a line model, normalized to a specific F160W mag
    lam, flam, fnu = unicorn.hudf.generate_line_model(z=12.13, filter='F160W', observed_mag=my_f160w_mag, own_spec=1000)
    twod.compute_model(lam, flam/1.e-17/twod.total_flux)
    check = unicorn.hudf.check_mags(lam, fnu*scale)
    
    #### Get the integrated line flux and double check that it agrees with the
    #### flux in the 2D model
    ok = (lam > 1.59e4) & (lam < 1.6e4)
    line_flux = np.trapz(flam[ok], lam[ok])
    sens = np.interp(1.5965e4, twod.oned.lam, twod.oned.sens)
    model_flux = np.sum(twod.model)/sens*1.e-17*np.median(np.diff(twod.im['WAVE'].data))
    
    #### Normalize the cross-correlation kernels
    thumb /= np.sum(thumb)
    kernel /= np.sum(kernel)
    
    #### Set up full padded arrays
    NY, NX = twod.im['SCI'].data.shape
    flux = np.zeros((3*NY, NX+2*NY))
    contam = flux*0.
    model = flux*0.
    err = flux*0.
    
    flux[NY:-NY, NY:-NY] = twod.im['SCI'].data
    contam[NY:-NY, NY:-NY] = twod.im['CONTAM'].data
    model[NY:-NY, NY:-NY] = twod.model
    err[NY:-NY, NY:-NY] = twod.im['WHT'].data
    
    #### Run all of the cross-correlation
    import stsci.convolve.Convolve as sc
    
    cd = sc.correlate2d(flux-contam, thumb, mode='constant', cval=0.)
    cf = sc.correlate2d(flux, thumb, mode='constant', cval=0.)
    cc = sc.correlate2d(contam, thumb, mode='constant', cval=0.)
    cm = sc.correlate2d(model, thumb, mode='constant', cval=0.)
    cm2 = sc.correlate2d(model, kernel, mode='constant', cval=0.)
    noise = np.random.normal(size=model.shape)*err #+model
    crnd = sc.correlate2d(noise, thumb, mode='constant', cval=0.)
    
    #### Figure showing cross-correlation
    sh = twod.model.shape
    xv = np.arange(1.1,1.61,0.1)
    ta = np.interp(xv, twod.im['WAVE'].data/1.e4, np.arange(sh[1]))
    
    vm = (-0.01, 0.01)
        
    NF = 3
    i=0
    
    aspect = (NF)*sh[0]*1./(2*sh[1])
    fig = unicorn.plotting.plot_init(square=True, left=0.1, right=0.01, top=0.08, bottom=0.1, hspace=0.0, wspace=0.0, aspect=aspect, xs=10, NO_GUI=False)
        
    ## Data/flux
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(twod.im['SCI'].data, vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_ylabel('Observed')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)
    
    a = ax.set_title('Pixel values')
    
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(cf[NY:-NY,NY:-NY], vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)

    a = ax.set_title('Cross-correlation')
    
    ## Contam
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(twod.im['CONTAM'].data, vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_ylabel('Contam')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)
    
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(cc[NY:-NY,NY:-NY], vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)

    ## Cleaned
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(twod.im['SCI'].data-twod.im['CONTAM'].data, vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_ylabel('Cleaned')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)
    
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(cd[NY:-NY,NY:-NY], vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels([]); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)

    ## Line
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(twod.model, vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_ylabel('Line model')
    a = ax.set_xticklabels(xv); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)
    a = ax.set_xlabel(r'$\lambda\,[\mu\mathrm{m}]$')
    
    i += 1; ax = fig.add_subplot((NF)*100+20+i)
    a = ax.imshow(cm[NY:-NY,NY:-NY], vmin=vm[0], vmax=vm[1], aspect='auto', interpolation='nearest')
    a = ax.set_xticklabels(xv); a = ax.set_yticklabels([]); a = ax.set_xticks(ta)
    a = ax.set_xlabel(r'$\lambda\,[\mu\mathrm{m}]$')
    
    fig.savefig('UDF_%05d.xcor_2D.pdf' %(id))
    
    #############
    comment = """
    1) Define an aperture where the model flux is greater than some threshold 
       at a redshift consistent with the line
    
    2) Run a MC simulation varying the pixel fluxes according to the 2D 
       spectrum errors.  
    
    3) Extract an aperture "flux" for each simulation and compare the 
       distribution to the observed xcor aperture flux
    """

    ### Subregion for faster cross-correlation
    model_sub = model[60:115, 275:319]
    err_sub = err[60:115, 275:319]
    cd_sub = cd[60:115, 275:319]

    aper = model_sub > 1.e-4
    xcor_flux = np.sum(cd_sub*aper)
    
    NSIM = 10000
    xcor_rand = np.zeros(NSIM)
    for i in xrange(NSIM):
        print unicorn.noNewLine + '%d' %(i)
        noise = np.random.normal(size=err_sub.shape)*err_sub #+model
        crnd = sc.correlate2d(noise, kernel, mode='constant', cval=0.)
        xcor_rand[i] = np.sum(crnd*aper)
        
    test = xcor_rand < xcor_flux
    print 'f(simulated < observed) = %.4f (N=%d)' %(test.sum()*1. / NSIM, NSIM)
    
    ### Figure for aperture simulation
    fig = unicorn.plotting.plot_init(square=True, xs=6, aspect=1, left=0.14)
    ax = fig.add_subplot(111)
    h = ax.hist(xcor_rand, bins=100, alpha=0.5, color='black')
    a = ax.plot([xcor_flux, xcor_flux], [0,h[0].max()], color='red', alpha=0.5, linewidth=3)
    a = ax.text(0.05,0.95,'f(sim. < obs.) = %.4f' %(test.sum()*1. / NSIM), ha='left', va='top', transform=ax.transAxes)
    a = ax.text(0.05,0.9,'N = %d' %(NSIM), ha='left', va='top', transform=ax.transAxes)
    a = ax.set_xlabel('x-correlation aper flux')
    a = ax.set_ylabel('N')
    fig.savefig('UDF_%05d.xcor_aper.pdf' %(id))
    
def evaluate_2d_errors(id=6001):
    """
    Make a figure with histograms of the observed cleaned flux and 
    random pixel values with a distribution given by the error extension
    """
    
    twod = unicorn.reduce.Interlace2D('UDF_%05d.2D.fits' %(id))
    
    ran = (-0.015, 0.015)
    alpha = 0.6
    
    fig = unicorn.plotting.plot_init(square=True, xs=6, aspect=1, left=0.14)
    ax = fig.add_subplot(111)
    
    h = ax.hist(twod.im['SCI'].data.flatten(), bins=100, alpha=alpha, color='orange', range=ran, label='Raw', histtype='step')
    h = ax.hist(twod.im['SCI'].data.flatten() - twod.im['CONTAM'].data.flatten(), bins=100, alpha=alpha, color='blue', range=ran, label='Cleaned', histtype='stepfilled')
    
    noise = np.random.normal(size=twod.im['SCI'].data.shape) * twod.im['WHT'].data #*0.93  # to get about perfect
    h = ax.hist(noise.flatten(), bins=100, alpha=alpha, color='red', range=ran, label='Random', histtype='stepfilled')
    a = ax.set_xlabel('G141 pixel value (e/s)')
    a = ax.legend(prop=dict(size=10))
    
    fig.savefig('UDF_%05d.pixel_values.pdf' %(id))
    
def fix_2d_background(object = 'UDF_06001', force=False):
    import unicorn.hudf as hudf
    
    if (not os.path.exists('%s.bg.fits' %(object))) | force:
        hudf.get_2d_background(object=object)
        
    #
    twod = unicorn.reduce.Interlace2D('%s.2D.fits' %(object))
    chain = unicorn.interlace_fit.emceeChain(file='%s.bg.fits' %(object))
    
    flux = twod.im['SCI'].data
    ny, nx = flux.shape
    yi, xi = np.indices(flux.shape)
    indices, step_sig = [], []
    f = 0.0005
    indices.append(np.ones(flux.shape)); step_sig.append(f)
    indices.append(xi-nx/2.); step_sig.append(f/(nx/2.))
    indices.append(yi-ny/2.); step_sig.append(f/(ny/2.))
    
    obj_fun = hudf._objective_poly2d
    
    background = obj_fun(chain.median, indices, flux, flux, 1)
    twod.im['SCI'].data -= background
    
    if 'BG_C0' in twod.im[0].header.keys():
        if np.log10(np.abs(twod.im[0].header['BG_C0']/chain.median[0])) < 0.001:
            print 'BG seems to already be applied.'
            return False
            
    twod.im[0].header.update('BG_C0', chain.median[0])
    twod.im[0].header.update('BG_CX', chain.median[1])
    twod.im[0].header.update('BG_CY', chain.median[2])
        
    twod.im.writeto(twod.im.filename(), clobber='True')
    
def get_2d_background(object = 'UDF_06001'):
    """
    Fit a linear surface to a 2D image for background subtraction
    """
    import time
    import numpy as np
    import emcee

    import unicorn.hudf as hudf
    
    #twod = unicorn.reduce.Interlace2D('PRIMO-1101_06001.2D.fits')
    #twod = unicorn.reduce.Interlace2D('UDF_06001.2D.fits')
    twod = unicorn.reduce.Interlace2D('%s.2D.fits' %(object))
    
    flux = twod.im['SCI'].data
    err = twod.im['WHT'].data
    contam = twod.im['CONTAM'].data
    
    cleaned = flux-contam
    var = err**2+contam**2
    
    ### Fit arrays for polynomial background
    ny, nx = flux.shape
    yi, xi = np.indices(flux.shape)
    indices, step_sig = [], []
    f = 0.0005
    indices.append(np.ones(flux.shape)); step_sig.append(f)
    indices.append(xi-nx/2.); step_sig.append(f/(nx/2.))
    indices.append(yi-ny/2.); step_sig.append(f/(ny/2.))
    
    obj_fun = hudf._objective_poly2d
    obj_args = [indices, cleaned, var, 0]
    
    init = np.zeros(len(indices))
    
    ndim, nwalkers = len(init), 20
    p0 = [(init+np.random.normal(size=ndim)*step_sig) 
          for i in xrange(nwalkers)]
    
    NTHREADS, NSTEP = 1, 200
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, args = obj_args, 
                                    threads=NTHREADS)
    
    #
    t0 = time.time()
    result = sampler.run_mcmc(p0, NSTEP)
    t1 = time.time()
    print 'Sampler: %.1f s' %(t1-t0)
    
    chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names=['c0', 'cx', 'cy'])
    
    fig = unicorn.plotting.plot_init(square=True, xs=6, aspect=1)
    for ip, p in enumerate(chain.param_names):
        ax = fig.add_subplot(311+ip)
        if ip == 0:
            ax.text(0.5,0.95,r'$N_\mathrm{Walk}/N_\mathrm{Step} = %d/%d$; params $\times$ 1000' %(chain.nwalkers, chain.nstep), ha='center' ,va='top', transform=ax.transAxes)
        #
        chain.show_chain(p, ax=ax, scale=1000)
        ax.plot([0,chain.nstep],[0,0], color='red', alpha=0.4)
        if ip != 2:
            ax.set_xticklabels([])
            
    fig.savefig('%s.bg.png' %(object))
    plt.close()
    
    chain.save_fits('%s.bg.fits' %(object))
    fp = open('%s.bg.dat' %(object),'w')
    fp.write('# id  c0  cx  cy\n')
    fp.write('%s  %.2e  %.2e  %.2e\n' %(object, chain.stats['c0']['q50'], chain.stats['cx']['q50'], chain.stats['cy']['q50']))
    fp.close()
    
    #### Subtract from 2D spectrum file
    #background = obj_fun(chain.median, indices, cleaned, var, 1)
    
    #### Test
    # input = [0.001, -0.00002, -0.00008]
    # test = obj_fun(input, indices, cleaned, var, 1)
    # test += np.random.normal(size=err.shape)*err
    # obj_args = [indices, test, var, 0]
    # init[0] = np.median(test)
    # p0 = [(init+np.random.normal(size=ndim)*step_sig) 
    #       for i in xrange(nwalkers)]
    # #
    # NTHREADS, NSTEP = 1, 200
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, args = obj_args, 
    #                                 threads=NTHREADS)
    # 
    # #
    # t0 = time.time()
    # result = sampler.run_mcmc(p0, NSTEP)
    # t1 = time.time()
    # print 'Sampler: %.1f s' %(t1-t0)
    # 
    # chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names=['c0', 'x', 'y'])
    # background = obj_fun(chain.median, indices, test, var, 1)
    
def _objective_poly2d(params, indices, data, var, ret):
    """
    Objective fitting function for a 2D polynomical array with order=`order`
    and cross terms:
    
    xi = [1,x,y,x**2,y**2,x*y,...]
    """
    N = len(params)
    F = data*0.
    for i in range(N):
        F += params[i]*indices[i]
    
    if ret == 1:
        return F
        
    return -0.5*np.sum((data-F)**2/var)
    
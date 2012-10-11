"""

Run the full reduction scripts on the ACS parallels

"""
import threedhst
import threedhst.prep_flt_files

import unicorn
import glob
import os
import shutil
import numpy as np

USE_PLOT_GUI=False

from threedhst.prep_flt_files import process_3dhst_pair as pair
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

def go_all():
    import unicorn.go_acs
    unicorn.go_acs.test(root='jbhm32')
    unicorn.go_acs.test(root='jbhm51')
    unicorn.go_acs.test(root='jbhm54')

    #### For Danilo
    unicorn.go_acs.test(root='jbhm39')
    
def test(root='jbhm51'):
    """
    COSMOS-23 ACS overlaps with COSMOS-25 WFC3
    """
    #os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/PREP_FLT')
    
    #root='jbhm51'
    files=glob.glob('../RAW/%s*fits*' %(root))
    for file in files:
        os.system('cp %s .' %(file))
        os.system('gunzip %s' %(os.path.basename(file)))
        print file
    
    #### Destripe + CTE correction
    threedhst.prep_flt_files.prep_acs(force=True)
    
    os.system('rm *flt.fits')
    
    asn_direct_file = root+'010_asn.fits'
    asn_grism_file = root+'020_asn.fits'
    
    #### Copy corrected FLT files to . 
    asn = threedhst.utils.ASNFile(asn_direct_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    #
    asn = threedhst.utils.ASNFile(asn_grism_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    
    #### Flag CRs, subtract background and align WCS    
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    ALIGN = '/3DHST/Spectra/Work/COSMOS/MOSAIC/COSMOS-F140w_11-05-23_sci.fits'
    ALIGN = '/3DHST/Spectra/Work/COSMOS/MOSAIC/COSMOS-F140w_11-09-08_sci.fits'
    ALIGN_EXTENSION = 0
    
    #### UDS
    ALIGN = '/Users/gbrammer/CANDELS/UDS/PREP_FLT/UDS-F125W_drz.fits'
    ALIGN_EXTENSION = 1
    
    threedhst.shifts.run_tweakshifts(asn_direct_file, verbose=True)

    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)

    threedhst.shifts.refine_shifts(ROOT_DIRECT=asn_direct_file.split('_as')[0].upper(),
              ALIGN_IMAGE=ALIGN, 
              ALIGN_EXTENSION = ALIGN_EXTENSION,
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=False,
        updatewcs=False, clean=True, median=False)
        
    ### Grism
    threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)

    threedhst.prep_flt_files.startMultidrizzle(asn_grism_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=False, median=True)
    
    ## Check
    threedhst.gmap.makeImageMap(['JBHM54010_drz.fits[1]*4', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-14-F140W_drz.fits', 'JBHM54020_drz.fits[1]', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-14-G141_drz.fits'], aper_list=[16])
    
    ###################### Run the grism
    os.system('cp *shifts.txt *tweak.fits *asn.fits ../DATA')
    
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS')
    unicorn.go_3dhst.set_parameters(direct='F814W', LIMITING_MAGNITUDE=24)

    threedhst.process_grism.set_ACS_G800L()
    
    #threedhst.options['SKY_BACKGROUND'] = None
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/'+root.upper()+'010_drz.fits'
    #threedhst.options['PREFAB_GRISM_IMAGE'] = '../PREP_FLT/'+root.upper()+'020_drz.fits'
    
    ## UDS lens: threedhst.options['FORCE_CATALOG'] = 'lens_drz.cat'
    
    threedhst.process_grism.reduction_script(asn_grism_file=root+'020_asn.fits')
    
    ### SEDs unicorn.analysis.make_SED_plots(grism_root='jbhm51020')
    #os.chdir('../')
    grism_root=asn_grism_file.split('_asn')[0]
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root, uds=True)
    
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root.upper(), GRISM_NAME='G800L')
        
    print 'Matched catalog'
    unicorn.analysis.match_grism_to_phot(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat',
                  OUTPUT_DIRECTORY=OUTPUT_DIRECTORY,
                  MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE)
    
    ## make figures
    for id in grismCat.id:
        status = unicorn.analysis.specphot(id=id, grism_root=grism_root, SPC = SPC, 
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')
    
    
    unicorn.go_3dhst.clean_up()
    
    os.system('rsync -avz HTML/ ~/Sites_GLOBAL/P/GRISM_ACS/')
    
def testing_g800l_background(asn_file='jbhm39020_asn.fits'):
    """
    Divide by the aXe sky image, divide by the row to row variation
    and subtract the overall average in the ACS grism images.
    """
    import threedhst.dq
    import numpy as np
    import threedhst.prep_flt_files
    import pyfits
    import matplotlib.pyplot as plt

    asn = threedhst.utils.ASNFile(asn_file)
    
    #### ACS has entries in run file for each of two WFC chips
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    inst = flt[0].header.get('INSTRUME').strip()
    if inst == 'ACS':
        skip=2
    else:
        skip=1
        

    sky1 = pyfits.open('../CONF/ACS.WFC.CHIP1.msky.1.fits')
    sky2 = pyfits.open('../CONF/ACS.WFC.CHIP2.msky.1.fits')
    skies = [sky1, sky2]
    extensions = [1,4] ### SCI extensions
    
    run = threedhst.prep_flt_files.MultidrizzleRun((asn_file.split('_asn.fits')[0]).upper())
    
    for i in range(len(asn.exposures)):
        exp = asn.exposures[i]
        flt = pyfits.open(exp+'_flt.fits', mode='update')
        ### Loop through ACS chips
        for j in [0,1]:
            run.blot_back(ii=i*skip+j, copy_new=(i == 0) & (j == 0), ACS_CHIP = j+1)
            
            threedhst.prep_flt_files.make_segmap(run.flt[i*skip+j], IS_GRISM=True)
            seg = pyfits.open(exp+'_flt.seg.fits')
        
            ext = extensions[j]
            
            print unicorn.noNewLine+' 1) Mask'
            data = flt[ext].data/skies[j][0].data
            data[seg[0].data > 0] = np.nan
            
            ### Collapse along rows
            print unicorn.noNewLine+' 2) Profile'
            shp = data.shape
            avg = np.zeros(shp[0])
            for k in range(shp[0]):
                cut = data[k,:]
                avg[k] = threedhst.utils.biweight(cut[np.isfinite(cut)], mean=True)
        
            sig = 5
            xkern = np.arange(8*sig)-4*sig
            ykern = np.exp(-1*xkern**2/sig**2/2.)
            ykern /= np.trapz(ykern, xkern)
            
            #### Need to make a larger array for smoothing at the edges
            print unicorn.noNewLine+' 3) Correct'
            avg_grow = np.ones(shp[0]+8*sig)
            avg_grow[0:4*sig] *= avg[0]
            avg_grow[4*sig:-4*sig] = avg
            avg_grow[-4*sig:] *= avg[-1]
            smooth = np.convolve(avg_grow, ykern, mode='same')[4*sig:-4*sig]
            model = np.dot(smooth.reshape((shp[0],1)), np.ones((1,shp[1])))
            
            arr = data[seg[0].data == 0]
            stats = threedhst.utils.biweight(arr, both=True)
            
            flt[ext].data /= skies[j][0].data * model/stats[0]
            #flt[ext].data -= model
            
            data /= skies[j][0].data * model/stats[0]
            arr = data[seg[0].data == 0]
            stats2 = threedhst.utils.biweight(arr, both=True)
            flt[ext].data -= stats2[0]
            
            flt[ext].header.update('MDRIZSKY', stats2[0])
        
        
        print unicorn.noNewLine+' 4) Update FLT'
        
        flt.flush()
            
def testing_f814w_background(asn_file='jbhm39020_asn.fits'):
    """
    Divide by the aXe sky image, divide by the row to row variation
    and subtract the overall average in the ACS grism images.
    """
    import threedhst.dq
    import numpy as np
    import threedhst.prep_flt_files
    import pyfits
    import matplotlib.pyplot as plt

    asn = threedhst.utils.ASNFile(asn_file)
    
    #### ACS has entries in run file for each of two WFC chips
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    inst = flt[0].header.get('INSTRUME').strip()
    if inst == 'ACS':
        skip=2
    else:
        skip=1
    extensions = [1,4] ### SCI extensions
    
    run = threedhst.prep_flt_files.MultidrizzleRun((asn_file.split('_asn.fits')[0]).upper())
    
    for i in range(len(asn.exposures)):
        exp = asn.exposures[i]
        flt = pyfits.open(exp+'_flt.fits', mode='update')
        ### Loop through ACS chips
        for j in [0,1]:
            run.blot_back(ii=i*skip+j, copy_new=(i == 0) & (j == 0), ACS_CHIP = j+1)
            
            threedhst.prep_flt_files.make_segmap(run.flt[i*skip+j], IS_GRISM=True)
            seg = pyfits.open(exp+'_flt.seg.fits')
        
            ext = extensions[j]
            
            print unicorn.noNewLine+' 1) Mask'
            data = flt[ext].data/np.ones(flt[ext].data.shape)
            data[seg[0].data > 0] = np.nan
            
            print unicorn.noNewLine+' 2) Subtract Global Mean'            
            arr = data[seg[0].data == 0]
            stats = threedhst.utils.biweight(arr, both=True)
            
            flt[ext].data -= stats[0]
            data -= stats[0]
         
        flt.flush()    

def check_lens():
    """
    Look at the ACS spectrum of the lens arc in the UDS
    """
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/UDS/DATA')
    drz = pyfits.open('jbhm19020_drz.fits')
    cont = pyfits.open('jbhm19020CONT_drz.fits')
    
    twod = pyfits.open('../HTML/images/jbhm19020_00280_2D.fits.gz')
    thumb = pyfits.open('../HTML/images/jbhm19020_00280_thumb.fits.gz')
    
    ds9 = threedhst.dq.myDS9()
    
    ds9.frame(1)
    ds9.view(twod[1])
    ds9.scale(-10,200)
    
    ds9.frame(2)
    ds9.view(twod[1].data-twod[4].data*900, header=twod[1].header)
    ds9.scale(-10,200)
    
    ds9.frame(3)
    img = twod[1].data*1.
    shp = thumb[0].data.shape
    img[:,-shp[1]:] = thumb[0].data*100
    ds9.view(img, header=twod[1].header)
    ds9.scale(-10,200)
    
    ### line tracing out the arc in teh 280 thumbnail
    trace = np.array([179.6091, 49.73919, 178.48754, 43.850966, 176.43133, 39.738555, 174.93591, 37.68235, 169.98232, 32.167981, 164.74835, 28.803281, 151.56994, 24.036623, 139.32617, 23.849696])
    
    tracex = np.array(trace[0::2]) - (twod[1].data.shape[1] - twod[1].data.shape[0])-1
    tracey = np.array(trace[1::2])-1
    
    plt.plot(tracex, tracey, color='yellow', alpha=0.5, linewidth=2, marker='o', ms=10)

    from pylab import polyfit as sp_polyfit
    from pylab import polyval as sp_polyval
    
    coeff = sp_polyfit(tracey, tracex, 5)
    yy = np.arange(1000)/1000.*shp[1]
    plt.plot(sp_polyval(coeff, yy), yy, color='green', alpha=0.8)
    
    thumb_shifted = thumb[0].data*1.
    #thumb_shifted[:,10] = 0
    plt.imshow(thumb_shifted, vmin=-0.1, vmax=0.5, interpolation='nearest')
    
    twod_shifted = twod[1].data*1.

    shp = thumb_shifted.shape
    for yi in range(22,50):
        shift_i = shp[1]/2. - sp_polyval(coeff, yi)
        # print yi, shift_i
        thumb_shifted[yi,:] = xshift(thumb_shifted[yi,:], shift_i)
        twod_shifted[yi,:] = xshift(twod_shifted[yi,:], shift_i)
        #twod_shifted[yi-1,:] += xshift(twod_shifted[yi-1,:], shift_i)
        #twod_shifted[yi+1,:] += xshift(twod_shifted[yi+1,:], shift_i)
    #
    plt.imshow(thumb_shifted, vmin=-0.1, vmax=0.5, interpolation='nearest')
    
    ds9.frame(4)
    img = twod_shifted*1.
    shp = thumb[0].data.shape
    img[:,-shp[1]:] = thumb_shifted*100
    ds9.view(img, header=twod[1].header)
    ds9.scale(-10,200)
    
    
    ### OBJ = 257, 280
    
def xshift(array, shift, spline=True):
    """
    Shift an array by an arbitrary `shift` value, using spline interpolation
    and wrapping around array edges.
    """
    import scipy.interpolate as interp
    shp = array.shape
    xpix = np.arange(shp[0])
    if spline:
        shifted = interp.spline(xpix, array, (xpix-shift) % shp[0])
    else:
        shifted = np.interp((xpix-shift) % shp[0], xpix, array)
    #
    return shifted
    

def process_acs_pair(asn_direct_file='ib3706050_asn.fits',
                     asn_grism_file='ib3706060_asn.fits',
                       ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
                       ALIGN_EXTENSION=0,
                       SKIP_GRISM=False,
                       adjust_targname=True,
                       align_geometry='shift',
                       PATH_TO_RAW='../RAW',
                       get_shift=True,
                       TWEAKSHIFTS_ONLY=False): 

    """
    Does the basic processing for ACS F814W and G800L pointings: background subtraction, allignment and drizzlign.
    """

    import threedhst
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import make_targname_asn
    
    #### Copy corrected FLT files to . 
    asn = threedhst.utils.ASNFile(asn_direct_file)
    for exp in asn.exposures:
        print exp
        os.system('rm %s_flt.fits' %(exp))
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    #
    asn = threedhst.utils.ASNFile(asn_grism_file)
    for exp in asn.exposures:
        print exp
        os.system('rm %s_flt.fits' %(exp))
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))

    #DIRECT REDUCTION
    ROOT_DIRECT = asn_direct_file.split('_asn.fits')[0]
        
    from threedhst.prep_flt_files import make_targname_asn
    
    #this makes new asn.fits files but with ACS the names start with ANY
    #must add an optional tag to replace ANY with the field name
    if (asn_direct_file is not None) & adjust_targname:
        asn_direct_file = make_targname_asn(asn_direct_file,field='GOODS-S')
    
    if (asn_grism_file is not None) & adjust_targname:
        asn_grism_file = make_targname_asn(asn_grism_file,field='GOODS-S')

    #run = threedhst.prep_flt_files.MultidrizzleRun((asn_direct_file.split('_asn.fits')[0]).upper())
    threedhst.shifts.run_tweakshifts(asn_direct_file, verbose=True)
    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)

    for i,exp in enumerate(asn.exposures):
        asn_mask = asn.exposures[i]+'_flt.fits.mask.reg'
        print asn_mask
        if os.path.exists(asn_mask):
            threedhst.showMessage("Apply ASN mask: %s" %(asn_mask))
            threedhst.regions.apply_dq_mask(asn.exposures[i]+'_flt.fits', extension=3, 
                    mask_file = asn_mask)

    threedhst.shifts.refine_shifts(ROOT_DIRECT=asn_direct_file.split('_as')[0].upper(),
              ALIGN_IMAGE=ALIGN_IMAGE, 
              ALIGN_EXTENSION = ALIGN_EXTENSION,
              fitgeometry=align_geometry, clean=True)
    
    unicorn.go_acs.testing_f814w_background(asn_direct_file)

    SCALE = 0.06
    PIXFRAC=1.0

    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=False,
        updatewcs=True, clean=True, median=False)

    #GRISM REDUCTION

    threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)

    threedhst.prep_flt_files.startMultidrizzle(asn_grism_file, use_shiftfile=True,
        skysub=True,
        final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=True,
        updatewcs=True, clean=False, median=True)
        
    unicorn.go_acs.testing_g800l_background(asn_grism_file)

    threedhst.prep_flt_files.startMultidrizzle(asn_grism_file, use_shiftfile=True,
        skysub=True,
        final_scale=SCALE, pixfrac=PIXFRAC, driz_cr=True,
        updatewcs=True, clean=False, median=True)
    
def make_external_catalog(root='', master_segmentation='', master_catalog='', reference_image=''):

    """
    Given a master catalog, a master segmentation map and a reference image, make a catalog and a segmentation
    map which contain only the sources within the reference image.
    """

    import pyfits
    import numpy as np
    import pywcs
    
    if (root is None) or (master_segmentation is None) or (reference_image is None):
        print 'Missing input: provide root, master segmentation and master catalog.'
        return False
    
    if reference_image is None:
        reference_image = '../PREP_FLT/'+root+'_drz.fits'

    catalog_out = root+'.ext.cat'
    segmentation_out = root+'-G800L_seg.fits'
    
    print reference_image
    
    old_cat = threedhst.sex.mySexCat(master_catalog)
    ref = pyfits.open(reference_image)

    sw_seg = threedhst.sex.SWarp()
    sw_seg.swarpMatchImage(reference_image)
    sw_seg.options['IMAGEOUT_NAME'] = segmentation_out
    sw_seg.options['RESAMPLE'] = 'N'
    sw_seg.options['FSCALASTRO_TYPE']='NONE'
    if reference_image.find('GOODS-S-08') != -1:
        sw_seg.options['CENTER']='03:32:47.08, -27:53:59.57'
    sw_seg.swarpImage(master_segmentation)
    print sw_seg.options['CENTER']

    seg = pyfits.open(segmentation_out)

    inter_seg = np.array(seg[0].data)
    inter_seg[np.where(ref[1].data==0)]=0

    old_cat.change_MAG_AUTO_for_aXe(filter='F806W')

    objects_in_seg = np.unique(inter_seg)
    pop_id = []
    for id in old_cat.id[::-1]:
        if id not in objects_in_seg:
            #print 'Pop #%05d' %(id)
            pop_id.append(id)
    old_cat.popItem(np.array(pop_id))
    
    pop_mag = []
    for ii in range(old_cat.nrows):
        if (old_cat['MAG_F806W'][ii] == 99.0000) | (old_cat['MAG_F806W'][ii] == 0.00):
            print old_cat.MAG_F806W[ii]
            pop_mag.append(old_cat.id[ii])
    old_cat.popItem(np.array(pop_mag))
        
    tmp1 = np.array(old_cat['A_IMAGE'])*0.06/3600.
    old_cat.addColumn(data=tmp1, format='%f', name='A_WORLD', comment='Profile RMS along major axis [deg]', verbose=True)
    tmp2 = np.array(old_cat['B_IMAGE'])*0.06/3600.
    old_cat.addColumn(data=tmp2, format='%f', name='B_WORLD', comment='Profile RMS along minor axis [deg]', verbose=True)

    ref_wcs = pywcs.WCS(ref[1].header)
    new_x = np.array(old_cat['X_IMAGE'])
    new_y = np.array(old_cat['Y_IMAGE'])
    for ii in range(old_cat.nrows):
        new_x[ii], new_y[ii] = ref_wcs.wcs_sky2pix([[old_cat['X_WORLD'][ii],old_cat['Y_WORLD'][ii]]],1)[0]

    old_cat.renameColumn(original='X_IMAGE', new='X_OLD', verbose=True)
    old_cat.renameColumn(original='Y_IMAGE', new='Y_OLD', verbose=True)
    old_cat.addColumn(data=new_x, format='%f', name='X_IMAGE', comment = '', verbose=True)
    old_cat.addColumn(data=new_y, format='%f', name='Y_IMAGE', comment = '', verbose=True)

    old_cat.write(outfile=catalog_out)


def specphot_acs_and_wfc3(id=69, grism_root='ibhm45030',
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
    CACHE_FILE = 'Same', Verbose=False,
    SPC = None, cat=None, grismCat = None,
    zout = None, fout = None, OUT_PATH='/tmp/', OUT_FILE_FORMAT=True,
    OUT_FILE='junk.png', GET_SPEC_ONLY=False, GET_WFC3=False, WFC3_DIR='/3DHST/Spectra/Release/v2.0/GOODS-S'):
    """
    specphot_acs_and_wfc3(id)
    
    Get photometry/SED fit as well as WFC3 spectrum when available and overplot G141 spectrum.
    This is different from unicorn.analysis.specphot() which does not get the WFC3 spectrum.
    """
    
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import pyfits

    #### Get G141 spectrum
    if Verbose:
        print 'Read SPC'
    
    if SPC is None:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                    axe_drizzle_dir='DRIZZLE_G141')
                    
    spec = SPC.getSpec(id)
    if spec is False:
        return False
        
    xmin = 3000
    xmax = 2.4e4
    
    lam = spec.field('LAMBDA')
    flux = spec.field('FLUX')
    ffix = flux-spec.field('CONTAM')
    ferr = spec.field('FERROR') #*0.06/0.128254
        
    if Verbose:
        print 'Read grism catalog'
        
    #### Read the grism catalog and get coords of desired object
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    #### Source size
    R = np.sqrt(np.cast[float](grismCat.A_IMAGE)*np.cast[float](grismCat.B_IMAGE))
    grism_idx = np.where(grismCat.id == id)[0][0]
    
    Rmatch = R[grism_idx]*1.
    
    ra0 = grismCat.ra[grismCat.id == id][0]
    de0 = grismCat.dec[grismCat.id == id][0]
    
    #### Read EAZY outputs and get info for desired object
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
    
    
    photom_idx = np.where(dr == np.min(dr))[0][0]
    
    drMatch = dr[photom_idx]*1.
    
    if drMatch > 2:
        return False
        
    if Verbose:
        print 'Read zout'
    if zout is None:    
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
        
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    if Verbose:
        print 'Read binaries'
        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(photom_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE)
         
    try:
        lambdaz, temp_sed_sm = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
    except:
        temp_sed_sm = temp_sed*1.
 
    wfc3_exist = False
    
    if GET_WFC3:
        wfc3_file_path = WFC3_DIR+"/*/1D/FITS/*%05d.1D.fits" %(cat.id[photom_idx])
        wfc3_file = glob.glob(wfc3_file_path)
        if wfc3_file != []:
            wfc3_spec = pyfits.open(wfc3_file[0])
            wfc3_exist = True
        else: 
            print 'No WFC3 spectrum.'

    if Verbose: 
        print 'Normalize spectrum'
            
    q = np.where((lam > 0.55e4) & (lam < 1.0e4) & (flux > 0))[0]
        
    if len(q) == 0:
        return False

    yint = np.interp(lam[q], lambdaz, temp_sed_sm)
        
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    
    total_err = np.sqrt((ferr)**2+(1.0*spec.field('CONTAM'))**2)*anorm
    

    if GET_SPEC_ONLY:
        if drMatch > 1:
            return False
        else:
            return lam, ffix*anorm, total_err, lci, fobs, efobs, photom_idx
         
    if Verbose:
        print 'Start plot'
        
    #### Make the plot
    threedhst.plotting.defaultPlotParameters()
    
    xs=5.8
    ys = xs/4.8*3.2
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[xs,ys],dpi=100)
    else:
        fig = Figure(figsize=[xs,ys], dpi=100)
    
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.13*4.8/xs, bottom=0.15*4.8/xs,right=1.-0.02*4.8/xs,top=1-0.10*4.8/xs)
    
    ax = fig.add_subplot(111)
    
    ymax = np.max((ffix[q])*anorm)
    
    if Verbose:
        print 'Make the plot'
        
    ax.plot(lambdaz, temp_sed_sm, color='red')
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.2, linewidth=1)
    
    #### Show own extraction
    sp1d = threedhst.spec1d.extract1D(id, root=grism_root, path='./HTML', show=False, out2d=False)#, GRISM_NAME='G800L')
    lam = sp1d['lam']
    flux = sp1d['flux']
    ffix = sp1d['flux']-sp1d['contam'] 
    ferr = sp1d['error']
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.6, linewidth=1)
    
    #### Show photometry + eazy template
    ax.errorbar(lci, fobs, yerr=efobs, color='orange', marker='o', markersize=10, linestyle='None', alpha=0.4)
    ax.plot(lambdaz, temp_sed_sm, color='red', alpha=0.4)

    if wfc3_exist:
        q_wfc3 = np.where((wfc3_spec[1].data.wave > 1.08e4) & (wfc3_spec[1].data.wave < 1.68e4) & (wfc3_spec[1].data.flux > 0))[0]
        yint_wfc3 = np.interp(wfc3_spec[1].data.wave[q_wfc3], lambdaz, temp_sed_sm)
        spec_wfc3 = (wfc3_spec[1].data.flux-wfc3_spec[1].data.contam)/wfc3_spec[1].data.sensitivity
        anorm_wfc3 = np.sum(yint_wfc3*spec_wfc3[q_wfc3])/np.sum(spec_wfc3[q_wfc3]**2)
        print 'Scaling factors: ', anorm, anorm_wfc3
        ax.plot(wfc3_spec[1].data.wave[q_wfc3], spec_wfc3[q_wfc3]*anorm_wfc3, color='blue',alpha=0.6, linewidth=1)

    ax.set_ylabel(r'$f_{\lambda}$')
    
    if plt.rcParams['text.usetex']:
        ax.set_xlabel(r'$\lambda$ [\AA]')
        ax.set_title('%s: \#%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
    else:
        ax.set_xlabel(r'$\lambda$ [$\AA$]')
        ax.set_title('%s: #%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
        
    #kmag = 25-2.5*np.log10(cat.ktot[photom_idx])
    kmag = cat.kmag[photom_idx]
    
    ##### Labels
    label = 'ID='+r'%s   K=%4.1f  $\log M$=%4.1f' %(np.int(cat.id[photom_idx]),
        kmag, fout.field('lmass')[photom_idx])
        
    ax.text(5e3,1.08*ymax, label, horizontalalignment='left',
      verticalalignment='bottom')
    
    
    label = 'R=%4.1f"' %(drMatch)
    if drMatch > 1.1:
        label_color = 'red'
    else:
        label_color = 'black'
    ax.text(2.2e4,1.08*ymax, label, horizontalalignment='right',
      color=label_color, verticalalignment='bottom')
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.1*ymax,1.2*ymax)
    
    if Verbose:
        print 'Save the plot'
    
    if OUT_FILE_FORMAT:
        out_file = '%s_%05d_SED.png' %(grism_root, id)
    else:
        out_file = OUT_FILE
        
    if USE_PLOT_GUI:
        fig.savefig(OUT_PATH+'/'+out_file,dpi=100,transparent=False)
        plt.close()
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(OUT_PATH+'/'+out_file, dpi=100, transparent=False)
    
    print unicorn.noNewLine+OUT_PATH+'/'+out_file
    
    if Verbose:
        print 'Close the plot window'


def reduce_acs(root='',LIMITING_MAGNITUDE=20., match_wfc3 = False, WFC3_DIR='/3DHST/Spectra/Release/v2.0/GOODS-S'):
    
    """
    Set parameters for the aXe reduction, run aXe, read photometric catalog, find matches, make plots of spectra with photometry.
    """

    import numpy as np

    field = root[0:-3]
    if LIMITING_MAGNITUDE is None:
        LIMITING_MAGNITUDE=20.

    os_command = 'cp PREP_FLT/'+root+'*shifts.txt PREP_FLT/'+root+'*tweak.fits PREP_FLT/'+root+'*asn.fits DATA/' 
    os.system(os_command)
    
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/'+field)
    unicorn.go_3dhst.set_parameters(direct='F814W', LIMITING_MAGNITUDE=LIMITING_MAGNITUDE)

    threedhst.process_grism.set_ACS_G800L()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/'+root+'-F814W_drz.fits'
    threedhst.options['PREFAB_GRISM_IMAGE'] = root+'-G800L_drz.fits'
    if threedhst.options['PREFAB_GRISM_IMAGE'] != '': 
        print "Copy PREFAB_GRISM_IMAGE to DATA."
        os.system('cp PREP_FLT/'+threedhst.options['PREFAB_GRISM_IMAGE']+' DATA/')
        
    threedhst.options['FORCE_CATALOG']=root+'.ext.cat'
    ##### Note: segmentation image has to accompany the force_catalog!

    threedhst.process_grism.reduction_script(asn_grism_file=root+'-G800L_asn.fits')
    
    ### SEDs unicorn.analysis.make_SED_plots(grism_root='jbhm51020')
    #os.chdir('../')
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=field)
    cat.star_flag=np.ones(cat.size)

    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    ## read grism outputs
    if unicorn.hostname().startswith('uni'): 
        BASE_PATH='/Volumes/robot/3DHST/Spectra/Work/ACS_PARALLEL/'+field+'/'
    if unicorn.hostname().startswith('hyp'):
        BASE_PATH='/3DHST/Spectra/Work/ACS_PARALLEL/'+field+'/'
    print 'BASE PATH:', BASE_PATH
        
    grism_root = root+'-G800L'    
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root.upper(), BASE_PATH=BASE_PATH, GRISM_NAME='G800L')

    print 'Matched catalog'
    unicorn.analysis.match_grism_to_phot(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat',
                  OUTPUT_DIRECTORY=OUTPUT_DIRECTORY,
                  MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE)
    
    ## make figures
    if match_wfc3:
        for id in grismCat.id:
            status = unicorn.go_acs.specphot_acs_and_wfc3(id=id, grism_root=grism_root, SPC = SPC, 
                cat = cat,
                grismCat = grismCat, zout = zout, fout = fout, 
                OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
                MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
                OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
                CACHE_FILE = 'Same',GET_WFC3 = True, WFC3_DIR=WFC3_DIR)
    else:
        for id in grismCat.id:
            status = unicorn.analysis.specphot(id=id, grism_root=grism_root, SPC = SPC, 
                cat = cat,
                grismCat = grismCat, zout = zout, fout = fout, 
                OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
                MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
                OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
                CACHE_FILE = 'Same')
    
    unicorn.go_3dhst.clean_up()

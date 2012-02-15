"""

Run the full reduction scripts on the ACS parallels

"""
import threedhst
import threedhst.prep_flt_files

import unicorn
import glob
import os
import shutil

import pyraf
from pyraf import iraf

from threedhst.prep_flt_files import process_3dhst_pair as pair

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
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/PREP_FLT')
    
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
    Divide by the aXe sky image and subtract average values by row in the
    ACS grism images.
    """
    import threedhst.dq
    import numpy as np
    import threedhst.prep_flt_files
    import pyfits
    
    os.chdir('/3DHST/Spectra/Work/ACS_PARALLEL/COSMOS/PREP_FLT')
    
    #asn_file = 'jbhm39020_asn.fits'
    asn = threedhst.utils.ASNFile(asn_file)
    
    #### ACS has entries in run file for each of two WFC chips
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    inst = flt[0].header.get('INSTRUME').strip()
    if inst == 'ACS':
        skip=2
    else:
        skip=1
        
    #flt = pyfits.open('jbhm39avq_flt.fits')
    
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
            # ds9.v(flt[ext].data, vmin=100, vmax=200)
            # ds9.v(flt[ext].data/sky2[0].data, vmin=100, vmax=200)
            # ds9.v(data, vmin=100, vmax=200)
            
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
            #plt.hist(arr, range=(stats[0]-4*stats[1], stats[0]+4*stats[1]), bins=1000)
            
            flt[ext].data /= skies[j][0].data * model/stats[0]
            #flt[ext].data -= model
            
            data /= skies[j][0].data * model/stats[0]
            arr = data[seg[0].data == 0]
            stats2 = threedhst.utils.biweight(arr, both=True)
            flt[ext].data -= stats2[0]
            
            flt[ext].header.update('MDRIZSKY', stats2[0])
        
        #
        print unicorn.noNewLine+' 4) Update FLT'
        
        flt.flush()
    
    if (1 == 0):
        threedhst.prep_flt_files.startMultidrizzle(asn_file, use_shiftfile=True,
            skysub=False,
            final_scale=0.05, pixfrac=1, driz_cr=True,
            updatewcs=True, clean=True, median=True)
        
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
    
    
#!/usr/bin/env python
# encoding: utf-8
"""
reduce.py

Own reduction in distorted frame, interlacing.

Need to implement:

    o Image combination with interlacing
        
        - Consider running Multidrizzle with PIXFRAC=0 and blotting back 
          to distorted frame. Doing everything (CR-rejection, etc) on our own will be a
          pain.
          
        - make an inverse variance extension for SExtractor
        
    o Make a grism model from scratch with the geometry specified in the aXe conf 
      files.
      
        - Output spectra could have arbitrary scale (A/pix) to make the thumbnails
          line up with the 2D spectra in native pixels.
    
$URL: https://subversion.assembla.com/svn/threedhst_internal/trunk/reduce.py $
$Author: gbrammer $
$Date: 2011-05-22 02:01:43 -0400 (Sun, 22 May 2011) $


"""
__version__ = " $Rev: 5 $"

import os
import glob
import shutil
import time

try:
    import astropy.io.fits as pyfits
    if not hasattr(pyfits, '__version__'):
        pyfits.__version__ = 'astropy'
except:
    import pyfits

import numpy as np

import matplotlib.pyplot as plt

USE_PLOT_GUI = False

import time

import threedhst
import unicorn
import unicorn.utils_c as utils_c

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import scipy.ndimage as nd

try:
    ### Latest STSCI_PYTHON (Python2.7)
    import stsci.tools.wcsutil as wcsutil
except:
    try:
        ### Older STSCI_PYTHON (Python2.5.4)
        import pytools.wcsutil as wcsutil
    except:
        print "No WCSUtil, unicorn.reduce won't work."

###### Globals for compute_model
# conf_file = 'WFC3.IR.G141.V2.5.conf'
conf_file = 'G141.test27s.gbb.conf'
#conf_file = 'G141.test30.conf'
conf = threedhst.process_grism.Conf(conf_file, path=os.getenv('THREEDHST')+'/CONF/').params
conf_grism = 'G141'
sens_files = {}
for beam in ['A','B','C','D','E']:
    if 'SENSITIVITY_'+beam in conf.keys():
        sens = pyfits.open(os.getenv('THREEDHST') + '/CONF/'+conf['SENSITIVITY_'+beam])[1].data
        
        sens_wave = np.arange(sens['WAVELENGTH'].min() - float(conf['DLDP_%s_1' %(beam)].split()[0])*2, sens['WAVELENGTH'].max() + float(conf['DLDP_%s_1' %(beam)].split()[0])*2, np.diff(sens['WAVELENGTH'])[0]/2.)
        sens_sens = np.interp(sens_wave, sens['WAVELENGTH'], sens['SENSITIVITY'], left=1.e-10, right=1.e-10)
        sens_err = np.interp(sens_wave, sens['WAVELENGTH'], sens['ERROR'], left=0, right=0)
        
        sens_data = {'WAVELENGTH': sens_wave, 'SENSITIVITY': sens_sens, 'ERROR':sens_err}
        sens_files[beam] = sens_data
        

#### wavelength limits
grism_wlimit = {'G141':[1.05e4, 1.70e4, 22., 1.4e4], 'G102':[0.76e4, 1.17e4, 10., 1.05e4], 'G800L':[0.5e4, 1.05e4, 20., 0.75e4], 'GRS':[1.35e4, 1.95e4, 5., 1.65e4]}

ZPs = {'F105W':26.2687, 'F125W':26.25, 'F140W':26.46, 'F160W':25.96, 'F606W':26.486, 'F814W':25.937, 'F435W':25.65777, 'F110W':26.822, 'F098M':25.667, 'F555W':25.718, 'F475W':26.059, 'F625W':25.907, 'F775W':25.665, 'F850LP':24.842}

### STMag zeropoints
ZPsST = {'F105W':27.6933, 'F125W':28.0203, 'F140W':28.4790, 'F160W':28.1875, 'F606W':26.664, 'F814W':26.786, 'F435W':25.155, 'F110W':28.4401, 'F098M':29.9456}

PLAMs = {'F105W':1.0552e4, 'F125W':1.2486e4, 'F140W':1.3923e4, 'F160W': 1.5369e4, 'F606W':5917.678, 'F814W':8059.761, 'F435W':4350., 'F775W':7750., 'F850LP':9000, 'ch1':3.6e4, 'ch2':4.5e4, 'K':2.16e4, 'U':3828., 'G':4870., 'R':6245., 'I':7676., 'Z':8872., 'F110W':1.1534e4, 'F098M':9864.1}

# BWs = {}
# for filt in PLAMs.keys():
#     if filt.startswith('F1'):
#         xf, yf = np.loadtxt('%s/%s.dat' %(os.getenv('iref'), filt), unpack=True)
#         yf = yf / np.trapz(yf, xf)
#         BWs[filt] = 1./yf.max()
BWs = {'F105W': 2700.6102248267212,
 'F125W': 2845.3485021278179,
 'F140W': 3840.2311460823043,
 'F160W': 2682.9565027757253}

try:
    import astropy.wcs as pywcs
except:
    try:
        import pywcs
    except:
        print 'No pywcs/stwcs found.'
    
def set_grism_config(grism='G141', chip=1, use_new_config=True, force=False, use_config_file=None):
    import unicorn.reduce as red
    if (red.conf_grism == grism) & (not force):
        return None
    #
    config_file = {'G102':'WFC3.IR.G102.V2.0.conf', 'G141':'WFC3.IR.G141.V2.5.conf', 'G800L':'ACS.WFC.CHIP2.Cycle13.5.conf', 'GRS':'WFIRST.conf'}
    
    if use_new_config:
        config_file = {'G102':'G102.test27s.gbb.conf', 'G141':'G141.test27s.gbb.conf', 'G800L':'ACS.WFC.CHIP2.Cycle13.5.conf', 'GRS':'WFIRST.conf'}
        #config_file = {'G102':'G102.test27s.gbb.conf', 'G141':'G141.test30.conf', 'G800L':'ACS.WFC.CHIP2.Cycle13.5.conf', 'GRS':'WFIRST.conf'}
    
    BEAMS = ['A','B','C','D','E']
    if grism == 'G800L':
        config_file[grism] = 'ACS.WFC.CHIP%d.Cycle13.5.conf' %(chip)
        BEAMS = ['A','B','C','D','E','F','G']
        
    if use_config_file is None:
        use_config_file = config_file[grism]
        
    red.conf = threedhst.process_grism.Conf(use_config_file, path=os.getenv('THREEDHST')+'/CONF/').params
    red.conf_grism = grism
    red.conf_file = use_config_file
    red.sens_files = {}
    #
    
    for beam in BEAMS:
        if 'SENSITIVITY_'+beam in conf.keys():
            sens = pyfits.open(os.getenv('THREEDHST') + '/CONF/'+conf['SENSITIVITY_'+beam])[1].data
            
            sens_wave = np.arange(sens['WAVELENGTH'].min() - float(conf['DLDP_%s_1' %(beam)].split()[0])*2, sens['WAVELENGTH'].max() + float(conf['DLDP_%s_1' %(beam)].split()[0])*2, np.diff(sens['WAVELENGTH'])[0]/2.)
            
            #### ACS grism sensitivity cuts out where still no-zero at > 1µm
            if grism == 'G800L':
                w_extend = np.arange(sens['WAVELENGTH'].max()+1, 1.2e4)
                if len(w_extend) > 0:
                    width = 300
                    sens_extend = np.exp(-(w_extend-w_extend[0])**2/2/width**2) * sens['SENSITIVITY'][-1]
                    err_extend = w_extend*0.+sens['ERROR'][-1]
                    #plt.plot(w_extend, sens_extend*1.e-17)
                    sens_wave0 = np.append(sens['WAVELENGTH'], w_extend)
                    sens_sens0 = np.append(sens['SENSITIVITY'], sens_extend)
                    sens_err0 = np.append(sens['ERROR'], err_extend)
                    #
                    sens_wave = np.arange(sens_wave0.min() - float(conf['DLDP_%s_1' %(beam)].split()[0])*2, sens_wave0.max() + float(conf['DLDP_%s_1' %(beam)].split()[0])*2, np.diff(sens['WAVELENGTH'])[0]/2.)
                    
            else:
                sens_wave0 = sens['WAVELENGTH']
                sens_sens0 = sens['SENSITIVITY']
                sens_err0 = sens['ERROR']

            sens_sens = np.interp(sens_wave, sens_wave0, sens_sens0, left=1.e-10, right=1.e-10)
            sens_err = np.interp(sens_wave, sens_wave0, sens_err0, left=1.e-10, right=1.e-10)
                        
            sens_data = {'WAVELENGTH': sens_wave, 'SENSITIVITY': sens_sens, 'ERROR':sens_err}
            red.sens_files[beam] = sens_data
            
            #
            #print beam, np.diff(sens['WAVELENGTH'])[0], sens['WAVELENGTH'].min(), sens['WAVELENGTH'].max()
            #plt.plot(sens_data['WAVELENGTH'], sens_data['SENSITIVITY']*1.e-17, label=beam)
            
    # for beam in ['A','B','C','D','E','F','G']:
    #     if 'SENSITIVITY_'+beam in conf.keys():
    #         red.sens_files[beam] = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+conf['SENSITIVITY_'+beam])[1].data
    #         red.sens_files[beam].SENSITIVITY[-3:] *= 0.
    #         red.sens_files[beam].SENSITIVITY[:3] *= 0.
        
    print 'Set grism configuration files for %s: %s' %(red.conf_grism, use_config_file)
    
def go_all(clean_all=True, clean_spectra=True, make_images=True, make_model=True, fix_wcs=True, extract_limit=None, skip_completed_spectra=True, MAG_LIMIT=26, out_path='./'):
    """
    clean_all=True; clean_spectra=True; make_images=True; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
    """
    import unicorn
    
    #### Set up the interlaced images
    for field in ['AEGIS','COSMOS','GOODS-S','UDS'][1:]:
        
        os.chdir(unicorn.GRISM_HOME+field+'/PREP_FLT')
        files=glob.glob(field+'-[0-9]*G141_asn.fits')
        #
        ### First run
        clean_all=True; clean_spectra=True; make_images=True; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
        extract_limit=23.5
        
        ### reference mosaic + segm image
        clean_all=False; make_images=False
        
        ### Redo
        clean_all=False; clean_spectra=False; make_images=False; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
        extract_limit=23.5
        
        ### Fainter extraction limit, don't regenerate everything already done
        clean_all=False; clean_spectra=False; make_images=False; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
        extract_limit=25
        
        #### Regenerate 2D FITS files
        skip_completed=False
        
        for file in files:
            status = unicorn.reduce.reduce_pointing(file=file, clean_all=clean_all, clean_spectra=clean_spectra, make_images=make_images, make_model=make_model, fix_wcs=fix_wcs, extract_limit=extract_limit, skip_completed_spectra=skip_completed_spectra, MAG_LIMIT=MAG_LIMIT, out_path=out_path)


def reduce_pointing(file='AEGIS-1-G141_asn.fits', clean_all=True, clean_spectra=True, make_images=True, make_model=True, fix_wcs=True, extract_limit=None, skip_completed_spectra=True, MAG_LIMIT=26, out_path='./'):
    import unicorn

    print file
    root=file.split('-G141')[0].split('-F140W')[0]
    
    remove = []
    if clean_spectra | clean_all:
        #print 'Clean spectra'
        remove.extend(glob.glob(root+'_*1D.fits'))
        remove.extend(glob.glob(root+'_*1D.png'))
        remove.extend(glob.glob(root+'_*2D.fits'))
        remove.extend(glob.glob(root+'_*2D.png'))
        remove.extend(glob.glob(root+'_*2D.xxx'))
    
    if clean_all:
        #print 'Clean all'
        remove.extend(glob.glob(root+'_*inter*'))
        remove.extend(glob.glob(root+'_model.*'))
        remove.extend(glob.glob(root+'_seg.fits'))
        remove.extend(glob.glob(root+'-[FG]*inter.fits'))
        
    if remove != []:
        for ff in remove:
            print unicorn.noNewLine+'Remove %s' %(ff)
            os.remove(ff)
            
    if make_images | (not os.path.exists(root+'-F140W_inter.fits')):
        unicorn.reduce.interlace_combine(root+'-F140W', view=False, use_error=True, make_undistorted=False)

    if make_images | (not os.path.exists(root+'-G141_inter.fits')):
        unicorn.reduce.interlace_combine(root+'-G141', view=False, use_error=True, make_undistorted=False)
    
    if make_model:
        model = unicorn.reduce.process_GrismModel(root, MAG_LIMIT=MAG_LIMIT)
    else:
        return True
        
    if fix_wcs:
        #model.trim_edge_objects()
        #model.get_corrected_wcs()
        #model.make_wcs_region_file()
        pass
        
    if extract_limit is not None:
        model.extract_spectra_and_diagnostics(MAG_LIMIT=extract_limit, skip=skip_completed_spectra)
    
    del(model)
    
    return True

def combine_all(FORCE=False):
    import os
    import glob
    
    import unicorn
    
    files=glob.glob('[AUGC]*[0-9]-G141_asn.fits')
    for file in files:
        pointing=file.split('-G141')[0]
        if (not os.path.exists(pointing+'-F140W_inter.fits')) | FORCE:
            unicorn.reduce.interlace_combine(pointing+'-F140W', view=False, pad=60, NGROW=125)
        #
        if (not os.path.exists(pointing+'-G141_inter.fits')) | FORCE:
            unicorn.reduce.interlace_combine(pointing+'-G141', view=False, pad=60, NGROW=125)

def acs_interlace_combine(asn_file):
    # asn_file = 'jbhm19010_asn.fits'
    # asn_file = 'jbhm19020_asn.fits'

    asn = threedhst.utils.ASNFile(asn_file)
    
    xoff, yoff = unicorn.reduce.acs_interlace_offsets(asn_file, growx=1, growy=1)
    
    xoff = np.array([  0,  14,  -4, -18])
    yoff = np.array([ 0, 24, 30,  7])
    
    xoff[1] -= 2; yoff[1] += 2
    xoff[2] -= 2; yoff[2] += 1
    xoff[3] -= 1; yoff[3] -= 2
    
    xoff = np.append(xoff, xoff)
    yoff = np.append(yoff, yoff)
    
    im = pyfits.open(asn.exposures[0]+'_flt.fits')
    full = np.zeros((2200, 4200))
    narr = full*0.
    x0, y0 = 50, 50
    tile = im[4].data*0
    cr_mask = (im[6].data & 4096) == 0
    bg = np.median(im[4].data[cr_mask & (im[4].data != 0)])
    tile[cr_mask] = im[4].data[cr_mask]-bg
    full[y0+0:y0+0+2048, x0+0:x0+0+4096] += tile
    narr[y0+0:y0+0+2048, x0+0:x0+0+4096] += cr_mask*1.
    
    ds9.view(full/np.maximum(narr, 1))
    
    for i in np.arange(1,len(asn.exposures)):
        im = pyfits.open(asn.exposures[i]+'_flt.fits')
        tile = im[4].data*0
        cr_mask = (im[6].data & 4096) == 0
        bg = np.median(im[4].data[cr_mask & (im[4].data != 0)])
        tile[cr_mask] = im[4].data[cr_mask] - bg
        full[y0+yoff[i]:y0+yoff[i]+2048, x0+xoff[i]:x0+xoff[i]+4096] += tile
        narr[y0+yoff[i]:y0+yoff[i]+2048, x0+xoff[i]:x0+xoff[i]+4096] += cr_mask*1.
        #
        ds9.frame(i+1)
        ds9.view(full/np.maximum(narr, 1))
        
def acs_interlace_offsets(asn_file, growx=2, growy=2, path_to_flt='./', chip=1):
    """
    Get interlacing pixel offsets from postargs
    """
    if asn_file is None:
        Nexp = 4
    else:
        asn = threedhst.utils.ASNFile(asn_file)
        Nexp = len(asn.exposures)
    
    # xpos, ypos = [], []
    # for exp in asn.exposures:
    #     head = pyfits.getheader(os.path.join(path_to_flt, exp+'_flt.fits'))
    #     xpos.append(head['POSTARG1'])
    #     ypos.append(head['POSTARG2'])
    #
    ### WFC3 postargs
    # xpos = [0.0, 1.355, 0.881, -0.474]
    # ypos = [0.0, 0.424, 1.212,  0.788]
    # 
    # a10 = 0.00
    # a11 = 0.0494
    # b10 = 0.0494
    # b11 = 0.0040
    # 
    # acsang = 92.16 - -45.123
    # xpos_acs, ypos_acs = threedhst.utils.xyrot(np.array(xpos), np.array(ypos), acsang)
    # xoff = xpos_acs / a11
    # yoff = (ypos_acs-xoff*b11)/b10
    # 
    # # xinter = -np.cast[int](np.round(xoff*10)/10.*growx)
    # # yinter = -np.cast[int](np.round(yoff*10)/10.*growy)
    # # 
    # xinter = -np.cast[np.int32](np.round(xoff*growx))
    # yinter = -np.cast[np.int32](np.round(yoff*growy))
    
    if (Nexp % 4) == 0:
        if chip == 1:
            xinter = np.array([ 0,  12,  -6, -18])*growx
            yinter = np.array([ 0, 25, 31,  6])*growy
        else:
            xinter = np.array([  0,  12,  -5, -18])*growx
            yinter = np.array([0, 24, 29, 5])*growy
    else:
        if chip == 1:
            xinter = np.array([ 0,  12,  -6, -18, -3, -7])*growx
            yinter = np.array([ 0, 25, 31,  6, 12, 2])*growy
        else:
            xinter = np.array([  0,  12,  -5, -18, -2, -7])*growx
            yinter = np.array([0, 24, 29, 5, 12, 2])*growy
        
    return xinter, yinter
    
def wcs_interlace_offsets(asn_file, growx=2, growy=2, path_to_flt='./', verbose=False, ref_exp=0, raw=False, reference_pixel=None):
    """
    Get interlace offsets from header WCS itself rather than POS-TARGS
    """
    import drizzlepac.skytopix
    import drizzlepac.pixtosky
    import astropy.io.fits as pyfits
    
    asn = threedhst.utils.ASNFile(asn_file)
    
    im = pyfits.open('%s/%s_flt.fits' %(path_to_flt, asn.exposures[0]))
    if reference_pixel is None:
        r0, d0 = im[1].header['CRVAL1'], im[1].header['CRVAL2']
        x0, y0 = im[1].header['CRPIX1'], im[1].header['CRPIX2']
    else:
        x0, y0 = reference_pixel
        rd = drizzlepac.pixtosky.xy2rd('%s/%s_flt.fits[sci,1]' %(path_to_flt, asn.exposures[0]), x=x0, y=y0, verbose=False, hms=False)
        r0, d0 = rd[0][0], rd[1][0]
        
    N = len(asn.exposures)
    xoff, yoff = np.zeros(N), np.zeros(N)
    for i in range(N):
        xoff[i], yoff[i] = drizzlepac.skytopix.rd2xy('%s/%s_flt.fits[sci,1]' %(path_to_flt, asn.exposures[i]), ra=r0, dec=d0, verbose=False)
    
    xoff -= xoff[ref_exp]
    yoff -= yoff[ref_exp]
    
    if raw:
        return -xoff*growx, -yoff*growy
    
    #xinter = -np.cast[np.int32](np.round((xoff)*10)/10.*growx)
    #yinter = -np.cast[np.int32](np.round((yoff)*10)/10.*growy)
    xinter = -np.cast[np.int32](np.round(xoff*growx))
    yinter = -np.cast[np.int32](np.round(yoff*growy))
    return xinter, yinter
    
def get_interlace_offsets(asn_file, growx=2, growy=2, path_to_flt='./', verbose=False, first_zero=False, raw=False):
    """
    Compute the necessary interlace offsets for a set of 
    dithered exposures defined in an ASN table
    """
    asn = threedhst.utils.ASNFile(asn_file)
    xpos, ypos = [], []
    for exp in asn.exposures:
        fitsfile = threedhst.utils.find_fits_gz(os.path.join(path_to_flt, exp+'_flt.fits'))
        head = pyfits.getheader(fitsfile)
        xpos.append(head['POSTARG1'])
        ypos.append(head['POSTARG2'])
    
    # central arcsec / pix, from WFC3 instrument HB
    a11 = 0.1355
    b10 = 0.1211 
    
    xoff = np.array(xpos)/a11
    yoff = np.array(ypos)/b10
    
    if raw:
        return xoff, yoff
        
    ### round and reverse offsets for combining the images
    if verbose:
        print xoff*growx, yoff*growy
        #print '%.2f %.2f' %(xoff*growy, yoff*growy)
    
    if first_zero:
        x0, y0 = xoff[0], yoff[0]
    else:
        x0, y0 = 0., 0.
        
    xinter = -np.cast[np.int32](np.round((xoff-x0)*10)/10.*growx)
    yinter = -np.cast[np.int32](np.round((yoff-y0)*10)/10.*growy)
    plot = """
    xinter = -(np.round(xoff*10)/10.*growx)
    yinter = -(np.round(yoff*10)/10.*growy)
    plt.plot(xinter % growx, yinter % growy, marker='o', ms=10, alpha=0.5)
    plt.xlim(0,growx); plt.ylim(0,growy)
    """
    
    #print xinter, yinter
    return xinter-xinter[0], yinter-yinter[0]
    
def interlace_combine(root='COSMOS-1-F140W', view=True, use_error=True, make_undistorted=False, pad = 60, NGROWX=180, NGROWY=30, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=False, ref_exp=0, nhotpix=2, clip_negative=-3, min_neighbors=2, reference_pixel=None):
    
    # from pyraf import iraf
    # from iraf import dither
    
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import scipy.ndimage as nd
    
    if unicorn.hostname().startswith('uni'):
        view = False
    
    if view:
        import threedhst.dq
        ds9 = threedhst.dq.myDS9()

    #
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    
    im = pyfits.open(os.getenv('iref')+'ir_wfc3_map.fits')
    PAM = im[1].data
    im.close()
    
    if os.path.exists(root+'.run'):
        run = threedhst.prep_flt_files.MultidrizzleRun(root)
        #run.blot_back(ii=0, copy_new=True)
        scl = np.float(run.scl)
        xsh, ysh = threedhst.utils.xyrot(np.array(run.xsh)*scl, np.array(run.ysh)*scl, run.rot[0])
    else:
        xsh, ysh = np.zeros(len(asn.exposures)), np.zeros(len(asn.exposures))
        
    yi,xi = np.indices((1014,1014))
        
    N = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX), dtype=np.int)
        
    inter_sci = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX))
    inter_weight = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX))
    if use_error:
        inter_err = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX))
    
    xi+=pad/(2*growx)+NGROWX
    yi+=pad/(2*growy)+NGROWY
    
    #### From pieter
    dxs = np.array([0,-20,-13,7]) + np.int(np.round(xsh[0]))*0
    dys = np.array([0,-7,-20,-13]) + np.int(np.round(ysh[0]))*0
    
    #### GOODS-N field from Weiner program 11600
    if root.startswith('GOODS-N') | root.startswith('GNGRISM') | root.startswith('goodsn'):
        dxs = np.array([0,-9,-4,5]) + np.int(np.round(xsh[0]))*0
        dys = np.array([0,-3,-11,-8]) + np.int(np.round(ysh[0]))*0
    
    if root.startswith('TILE41'):
        dxs = np.array([-4, 9, -8, 5, 5, -8]) + np.int(np.round(xsh[0]))*0
        dys = np.array([-5, -4, 4, 5, 5, 4]) + np.int(np.round(ysh[0]))*0
        if len(asn.exposures) == 12:
            dxs = np.array([-4, 9, -8, 5, 5, -8, -4, 9, -8, 5, 5, -8]) + np.int(np.round(xsh[0]))*0
            dys = np.array([-5, -4, 4, 5, 5, 4, -5, -4, 4, 5, 5, 4]) + np.int(np.round(ysh[0]))*0
    
    if root.startswith('GEORGE'):
        dxs = np.array([0, 9, -4, 5, -8]) + np.int(np.round(xsh[0]))*0
        dys = np.array([0, -4, -5, 5, 4]) + np.int(np.round(ysh[0]))*0
    
    #### Rigby lens galaxy
    if root.startswith('RCS0327'):
        dxs = np.array([0, -26]) + np.int(np.round(xsh[0]))*0
        dys = np.array([0, -57]) + np.int(np.round(ysh[0]))*0
        
    #### Cooper program in AEGIS
    if root.startswith('EGS1'):
        auto_offsets=True
    
    if 'GDN' in root:
        if 'G1' in root:
            dxs = np.array([   0, -150,  -13,    7])
            dys = np.array([  0,  -7, -20, -13])
        else:
            dxs = np.array([-150, 7])
            dys = np.array([-7, -12])
            
    #### Erb quasar sightlines
    if flt[0].header['PROPOSID'] == 12471:
        dxs, dys = np.array([  0, -10])*growx, np.array([ 0, -7])
    
    if 'ERSII' in root:
        dxs = np.array([   0, -147, -148,-1])*growx
        dys = np.array([ 0,  0, 83, 83])*growy
       
    if auto_offsets:
        #xinter, yinter = red.wcs_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy, reference_pixel = reference_pixel)
        try:
            xinter, yinter = red.wcs_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy, reference_pixel = reference_pixel)
        except:
            xinter, yinter = red.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy)
            
        dxs = xinter #- xinter[ref_exp] #+ np.int(np.round(xsh[0]))*0
        dys = yinter #- yinter[ref_exp] #+ np.int(np.round(ysh[0]))*0
                
    dxs += ddx
    dys += ddy
    
    dxs -= dxs[ref_exp]
    dys -= dys[ref_exp]
    
    print 'Interlace offsets: ', dxs, dys
    
    dxi = np.cast[np.int32](np.ceil(dxs/growx))
    dyi = np.cast[np.int32](np.ceil(dys/growy))
    
    #### Find hot pixels, which are flagged as cosmic 
    #### rays at the same location in each of the flt files.  Ignore others.
    hot_pix = np.zeros((1014,1014),dtype='int')
    
    for flt in asn.exposures:
        im = pyfits.open(flt+'_flt.fits')
        hot_pix += (im[3].data & 4096) / 4096
        
    threedhst.showMessage('Flagged %d "hot pixels" marked as CRs in > %d images' %((hot_pix > nhotpix).sum(), nhotpix))
    hot_pix = hot_pix > nhotpix #2 # (len(asn.exposures)/2.)
    
    for i,flt in enumerate(asn.exposures):
        print flt
        #flt = run.flt[i]
        im = pyfits.open(flt+'_flt.fits')
        
        #### Use the pixel area map correction
        im[1].data *= PAM
        #### Divide by 4 to conserve surface brightness with smaller output pixels
        im[1].data /= 1.*(growx*growy)
        im[2].data /= 1.*(growx*growy)
        
        ### Mask cosmic rays
        if i == 0:
            h0 = im[0].header
            h1 = im[1].header
            header = red.scale_header_wcs(h1.copy(), factor=2, growx=growx, growy=growy, pad=pad,
                NGROWX=NGROWX, NGROWY=NGROWY)
            header['EXTNAME']='SCI'
            header['PAD']=pad
            header['GROWX']=growx
            header['GROWY']=growy
            header['NGROWX'] = NGROWX
            header['NGROWY'] = NGROWY
            header['REFIMAGE'] = ('','Source detection image')
            header_wht = header.copy()
            header_wht['EXTNAME']='ERR'
        
        dx = np.int(np.round((xsh[i]-xsh[0])*growx))
        dy = np.int(np.round((ysh[i]-ysh[0])*growy))
        #
        dx = dxs[i]
        dy = dys[i]
        #
        #use = ((im[3].data & 4096) == 0) & ((im[3].data & 4) == 0) #& (xi > np.abs(dx/2)) & (xi < (1014-np.abs(dx/2))) & (yi > np.abs(dy/2)) & (yi < (1014-np.abs(dy/2)))
        #        
        #use = ((im[3].data & (4+32+16+2048+4096)) == 0) & (~hot_pix)
        use = ((im[3].data & (4+16+2048+4096)) == 0) & (~hot_pix) & (im[2].data > 0) #& (im[2].data)
        
        #### Some pixels with discrepant small errors
        use = use & (im[2].data > 0.8*np.median(im[2].data[use])) & (im[2].data < 20000) & np.isfinite(im[2].data)
        
        #### Pixels more than -3 sigma negative
        neg = (im[1].data/im[2].data > clip_negative)
        print unicorn.noNewLine + '%s, Clip negative (< %.1f sigma): %4d / %.4f' %(flt, clip_negative, (use & ~neg).sum(), (use & ~neg).sum()*1./use.sum())
        use = use & neg
        
        ### debug
        # errs = im[2].data[use]
        # print im.filename(), errs.min(), np.median(errs), errs.max()
         
        if use_error:
            inter_sci[yi[use]*growy+dy,xi[use]*growx+dx] += im[1].data[use]/im[2].data[use]**2
            inter_weight[yi[use]*growy+dy,xi[use]*growx+dx] += 1./im[2].data[use]**2
        else:
            inter_sci[yi[use]*growy+dy,xi[use]*growx+dx] += im[1].data[use] 
            inter_weight[yi[use]*growy+dy,xi[use]*growx+dx] += 1
            
        N[yi[use]*growy+dy,xi[use]*growx+dx] += 1
        
        # if use_error:
        #     inter_err[yi[use]*growy+dy,xi[use]*growx+dx] += im[2].data[use]**2
        
        if view:
            ds9.view_array(inter_sci/np.maximum(N,1), header=header)
            ds9.scale(-0.1,5)
    
    #### Average for case when dither positions overlap, e.g. CANDELS SN fields
    #inter_sci /= np.maximum(N,1) 
    #inter_err = np.sqrt(inter_err) / np.maximum(N, 1)
    inter_weight[inter_weight == 0] = 1
    inter_sci = inter_sci / inter_weight
    inter_err = np.sqrt(1./inter_weight)
    inter_err[N == 0] = 0.
    #inter_err[~np.isfinite(inter_err)] = 0.
    
    if use_error:
        h0.update('WHTERROR',True,comment='WHT extension is FLT[err,1]')
    else:
        h0.update('WHTERROR',False,comment='WHT extension is 0/1 flagged bad pix')
    
    
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=np.cast[np.float32](inter_sci), header=header)
    if use_error:
        wht = pyfits.ImageHDU(data=np.cast[np.float32](inter_err), header=header_wht)
    else:
        wht = pyfits.ImageHDU(data=np.cast[np.int](N), header=header_wht)
            
    #### Mask singleton pixels 
    ok = (wht.data != 0)*1
    neighbor = nd.convolve(ok, np.ones((3,3))) - 1 # self
    mask = (wht.data != 0) & (neighbor < min_neighbors)
    wht.data[mask] = 0
    sci.data[mask] = 0
    mask_fraction = mask.sum()*1./ok.sum()
    threedhst.showMessage('Masked %d (%.3f) pixels with < %d valid neighbors.\n\nMost will be at image edge resulting from the dither extremes.' %(mask.sum(), mask_fraction, min_neighbors), warn=mask_fraction > 0.1)
    
    image = pyfits.HDUList([hdu,sci,wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.set('EXTEND', True, after='NAXIS')
    
    
    image.writeto(root+'_inter.fits', clobber=True)
    
    #pyfits.writeto('inter_N.fits', data=N, header=image[1].header, clobber=True)
    
    if make_undistorted:
        from pyraf import iraf
        from iraf import dither
        
        try:
            os.remove(root+'_inter_sci.fits')
            os.remove(root+'_inter_wht.fits')
        except:
            pass
        #
        new_coeffs_dat(input=asn.exposures[0]+'_flt_coeffs1.dat', output='scale_coeffs.dat')

        im = pyfits.open(root+'_inter.fits', mode='update')
        sf = threedhst.shifts.ShiftFile(root+'_shifts.txt')

        shift = threedhst.utils.xyrot(np.array([sf.xshift[0]]), np.array([sf.yshift[0]]), 180.-im[1].header['PA_APER'])

        im[1].header['CRPIX1'] = im[1].header['CRPIX1']+2*shift[0][0]
        im[1].header['CRPIX2'] = im[1].header['CRPIX2']+2*shift[1][0]
        im[2].header['CRPIX1'] = im[2].header['CRPIX1']+2*shift[0][0]
        im[2].header['CRPIX2'] = im[2].header['CRPIX2']+2*shift[1][0]

        im.flush()

        status = iraf.wdrizzle(data = root+'_inter.fits[1]', 
                               outdata = root+'_inter_sci.fits',
         outweig = root+'_inter_wht.fits', outcont = "", in_mask = root+'_inter.fits[2]', 
         wt_scl = 'exptime', 
         outnx = 2300, outny = 2300, geomode = 'wcs', kernel = 'square',
         pixfrac = 1.0, coeffs = "scale_coeffs.dat", lamb = 1392., xgeoim = "", 
         ygeoim = "", align = 'center', scale = 1.0, xsh = -10.0, ysh = 0.0, rot = 0.0,
         shft_un = 'output', shft_fr = 'output', outscl = 0.060,
         raref = header['CRVAL1'], decref = header['CRVAL2'], xrefpix = 1150, 
         yrefpix = 1150,
         orient = 0.0, dr2gpar = "", expkey = 'exptime', in_un = 'cps',
         out_un = 'cps', fillval = '0', mode = 'al')

def new_coeffs_dat(input='ibhm29wlq_flt_coeffs1.dat', output='scale_coeffs.dat', factor=2, pad=60):
    fp = open(input)
    lines = fp.readlines()
    fp.close()

    spl = np.cast[float](lines[2].split()[1:])
    lines[2] = 'refpix %f %f\n' %(spl[1]*factor+pad/2, spl[1]*factor+pad/2)
    
    ids = [4,8]
    for i in ids:
        spl = np.cast[float](lines[i].split())
        spl[1:3] /= factor
        spl[3:] /= factor**2
        lines[i] = ''
        for sp in spl:
            lines[i] += '%16e' %(sp)
        lines[i] += '\n'
    #
    ids = [5,9]
    for id in ids:
        spl = np.cast[float](lines[id].split())
        spl[0] /= factor**2
        spl[1:] /= factor**3
        lines[id] = ''
        for sp in spl:
            lines[id] += '%16e' %(sp)
        lines[id] += '\n'
    #
    ids = [6,10]
    for id in ids:
        spl = np.cast[float](lines[id].split())
        spl /= factor**4
        lines[id] = ''
        for sp in spl:
            lines[id] += '%16e' %(sp)
        lines[id] += '\n'
    
    open(output,'w').writelines(lines)
    
def scale_header_wcs(header, factor=2, growx=2, growy=2, pad=60, NGROWX=180, NGROWY=30):
    """
    Take an FLT header but scale the WCS keywords to account for the images 
    expanded/contracted by a factor of 'factor'.
    """
    import numpy as np
    
    header['NAXIS1'] = header['NAXIS1']*growx+NGROWX*growx*factor+pad
    header['NAXIS2'] = header['NAXIS2']*growy+NGROWY*growy*factor+pad

    header['CRPIX1'] = header['CRPIX1']*growx+NGROWX*growx+pad/2
    header['CRPIX2'] = header['CRPIX2']*growy+NGROWY*growy+pad/2
    
    ### SIP WCS keywords for distortion
    if ('A_ORDER' in header.keys()) & ('SIP' in header['CTYPE1']):
        a_order = header['A_ORDER']
        b_order = header['B_ORDER']
        for a in range(a_order):
            for b in range(b_order):
                for k in 'AB':
                    key = '%s_%d_%d' %(k, a, b)
                    if key in header.keys():
                        header.update(key, header[key]/growx**a/growy**b)
                        
    # keys = ['CD1_1','CD1_2','CD2_1','CD2_2']
    # for key in keys:
    #     header.update(key,header[key]/factor)
    
    header['CD1_1'] = header['CD1_1']/growx
    header['CD2_1'] = header['CD2_1']/growx
    header['CD1_2'] = header['CD1_2']/growy
    header['CD2_2'] = header['CD2_2']/growy
    
    header['IDCTAB'] = ''
    
    header['EXPTIME'] = 1
    
    return header
    
def grism_model(xc_full=244, yc_full=1244, lam_spec=None, flux_spec=None, grow_factor=2, growx=2, growy=2, pad = 60, ngrowx=0, ngrowy=0, BEAMS=['A','B','C','D'], dydx=True, grism='G141', smooth_binned_sigma=0.001, xmi=1000, xma=-1000, zeroth_position=False, conf=None, sens_files=None, finex=0, finey=0):
    """
    Main code for converting the grism dispersion calibration to a full 2D model spectrum
    
    The input lam_spec / flux_spec model is interpolated to the wavelength grid.  Need to 
    smooth it by a factor at least of order 1 pixel to get emission lines right 
    (smooth_binned_sigma)
    """
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    try:
        import threedhst.dq
    except:
        ### Probably no pysao installed
        pass
        
    import scipy.ndimage as nd
    
    if conf is None:
        conf = red.conf
    
    if sens_files is None:
        sens_files = red.sens_files
    
    #ds9 = threedhst.dq.myDS9()
    
    ## A star in COSMOS-15-F140W
    ##### For interlaced images
    
    # xc_full = 244.13
    # yc_full = 1244.323
    
    if zeroth_position:
        BEAMS = ['B']
        
    #### Coordinates in FLT frame
    ## xc_full[interlaced] = (x[FLT orig]+ngrow)*growx + pad/2
    xc = np.int(np.round((xc_full - pad/2. - ngrowx*growx)/growx))
    yc = np.int(np.round((yc_full - pad/2. - ngrowy*growy)/growy))
    
    if 'XOFF' in conf.keys():
        xoff = np.float(conf['XOFF'])
    else:
        xoff = 0.
    if 'YOFF' in conf.keys():
        yoff = np.float(conf['YOFF'])
    else:
        yoff = 0.
    
    bigX = xc - xoff
    bigY = yc - yoff
      
    #xmi, xma = -580, 730
    #xmi, xma = 1000,-1000
    # limit = {'G141': {'A':(10,213), 'B':(-210,-170), 'C':(207,464), 'D':(469,720), 'E':(-600,-400)},
    #          'G102': {'A':(38,248), 'B':(-280,-240), 'C':(330,670), 'D':(670,1014), 'E':(-740,-560)},
    #          'G800L': {'A':(-30,160), 'B':(-140,-80), 'C':(120,410), 'D':(260,660), 'E':(-590,-220), 'F':(-540,-300), 'G':(-980,-450)}}
    limit = {'G141': {'A':(10,213), 'B':(-210,-170), 'C':(207,464), 'D':(469,720), 'E':(-600,-400)},
             'G102': {'A':(38,248), 'B':(-280,-240), 'C':(330,670), 'D':(670,1014), 'E':(-740,-560)},
             'G800L': {'A':(-30,160), 'B':(-140,-80), 'C':(120,410), 'D':(260,660), 'E':(-590,-220), 'F':(-540,-300), 'G':(-980,-450)},
             'GRS': {'A':(-350,350)}}
    
    for BEAM in limit[grism].keys():
        if BEAM in BEAMS:
            lim = limit[grism][BEAM]
            xmi, xma = min(xmi,lim[0]), max(xma,lim[1])
            
    # if 'A' in BEAMS:
    #     xmi, xma = min(xmi,10), max(xma,213)
    # if 'B' in BEAMS:
    #     xmi, xma = min(xmi,-210), max(xma, -170)
    # if 'C' in BEAMS:
    #     xmi, xma = min(xmi,207), max(xma,454)
    # if 'D' in BEAMS:
    #     xmi, xma = min(xmi, 469), max(xma,668)
    # if 'E' in BEAMS:
    #     xmi, xma = min(xmi, -600), max(xma,-400)
        
    NX,NY = xma-xmi, 8

    NX,NY = xma-xmi, 15
    
    if grism == 'G800L':
        NY = 60
        
    xmi *= growx
    xma *= growx
    NX *= growx
    NY *= growy
    
    model = np.zeros((NY*2+1, NX), dtype=np.float)
    wavelength = model*0
    full_sens = model*0
    beam_index = np.zeros(NX, dtype=np.int)
    
    #print model.shape, NX, NY
    
    yi, xi = np.indices(model.shape)
    xi += xmi
    yi -= NY
    xarr, xpix, yarr = np.arange(xma-xmi)+xmi, np.arange(xma-xmi), np.arange(NY*2+1)-NY
    
    #
    # This wasn't needed and was giving dz/1+z = +0.0035 offset w.r.t spec-zs! 
    ### Still have offset, adding it back in gets wave/flux standards right
    ### with grow=1
    #xarr += 1
    #xarr += growx
            
    for ib, beam in enumerate(BEAMS):
        #print beam
        
        xmi_beam, xma_beam = np.cast[int](conf['BEAM'+beam].split())
        xmi_beam *= growx
        xma_beam *= growx
        
        dydx_order = np.int(conf['DYDX_ORDER_'+beam])
        dydx_0 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_0']) * growy
        dydx_1 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_1']) * growy / growx#/ grow_factor        
                    
        xoff_beam = field_dependent(bigX, bigY, conf['XOFF_'+beam]) * growx
        yoff_beam = field_dependent(bigX, bigY, conf['YOFF_'+beam]) * growy
        
        xoff_beam += finex
        yoff_beam += finey
        
        xmi_beam += int(np.round(xoff_beam))
        xma_beam += int(np.round(xoff_beam))
                
        disp_order = np.int(conf['DISP_ORDER_'+beam])
        dldp_0 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_0'])
        dldp_1 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_1']) / growx
        
        #### Improve alignment of zeroth order
        #if (beam == 'B') & (grism == 'G141'):
        #    dldp_0+=1500
            
        # if beam == 'Bx':
        #     #dydx_1 = 0.0
        #     dydx_1 = 0.1
        #     dydx_0 -= 2/growy
        #     dydx_0 += dydx_1 * 192 * growy
        #     
        #     dldp_x = dldp_1*1.
        #     f = 0.9*(1+np.abs(bigY-507)/507.*0.2)
        #     dldp_1 *= f
        #     dldp_0 += (1-f)*dldp_1*-213*2*(2-(1.02-0.02*np.abs(bigY-507)/507.))
            
        #print 'BEAM_%s: %5.2f %5.2f, %6.3f %5.3f, %9.2f %6.2f' %(beam, xoff_beam / grow_factor, yoff_beam / grow_factor, dydx_0, dydx_1*grow_factor, dldp_0, dldp_1*grow_factor)

        #### Wavelength
        dldp_x = (xarr-xoff_beam)
        dldp_coeffs = [dldp_0, dldp_1]
        
        ### Still have offset for 0th order
        if beam == 'B':
            dldp_x = (xarr+growx-xoff_beam)
        
        #if beam == 'A':
        #    print lam.min(), lam.max()
            
        if 'DLDP_'+beam+'_2' in conf.keys():
            dldp_2 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_2']) / growx**2
            #lam = dldp_0 + dldp_1*(xarr-xoff_beam) + dldp_2*(xarr-xoff_beam)**2
            dldp_coeffs.append(dldp_2)
            
        #
        #lam = dldp_0 + dldp_1*(xarr-xoff_beam)
        lam = xarr*0.
        for i in range(len(dldp_coeffs)):
            lam += dldp_coeffs[i]*dldp_x**i
        
        # if growx == 1:
        #     print 'XOFF!'
        #     lam += 0.5*(lam[1]-lam[0])
            
        if not dydx:
            dydx_0 = 0.
            dydx_1 = 0.
        
        #print 'DYDX: ', dydx, dydx_0, dydx_1
        
        #### Interpolate pixel at shifted ycenter along beam        
        ycenter = dydx_0 + dydx_1*(xarr-xoff_beam)
        
        if 'DYDX_'+beam+'_2' in conf.keys():
            dydx_2 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_2']) / growx**2
            ycenter = dydx_0 + dydx_1*(xarr-xoff_beam) + dydx_2*(xarr-xoff_beam)**2
        
        ### Offsets for G141 first order, 0.5 pix in y, 1.5 pix in lam (2x2)
        #dy_fudge = 0
        dy_fudge = yoff_beam
        
        # if (grism == 'G141') & (beam == 'A'):
        #     #dy_fudge = 0.5*growy/2
        #     #
        #     # dx_fudge = -1.*growx/2*dldp_1
        #     # lam += dx_fudge
        #     pass
            
        ycenter += dy_fudge
           
        ### Vertical offsets of the G800L beams
        # if grism == 'G800L':
        #     off = {'A':-1, 'B':0, 'C':-3, 'D':-3, 'E':-2, 'F':-1}
        #     if beam in off.keys():
        #         ycenter += off[beam]
        
        #### Get offset to zeroth order only
        if zeroth_position:
            lam_ref = grism_wlimit[grism][3]
            #print xoff_beam, lam_ref
            dx = np.interp(lam_ref, lam, xarr)
            dy = np.interp(dx, xarr, ycenter)
            return dx, dy
            
        stripe = model*0
        y0 = np.cast[int](np.floor(ycenter))
        f0 = ycenter-y0
        keep = (xarr >= xmi_beam) & (xarr <= xma_beam)
        if xarr[keep].size > 1:
            stripe[y0[keep]+NY,xpix[keep]]  = 1-f0[keep]
            stripe[y0[keep]+NY+1,xpix[keep]] = f0[keep]
            wavelength[y0[keep]+NY,xpix[keep]] += lam[keep]
            wavelength[y0[keep]+NY+1,xpix[keep]] += lam[keep]
            #
            beam_index[xpix[keep]] += 2**ib
            
        yoff_save = dydx_0+dydx_1*(xmi_beam+xma_beam)/2
        if beam == 'A':
            yoff_array = dydx_0+dydx_1*(xarr-xoff_beam)
            if 'DYDX_A_2' in conf.keys():
                yoff_array = dydx_0 + dydx_1*(xarr-xoff_beam) + dydx_2*(xarr-xoff_beam)**2
            
            yoff_array += dy_fudge
            
            #print 'Yoff_array', xarr.shape, yoff_array.shape
            
        #### Sensitivity
        #sens = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+conf['SENSITIVITY_'+beam])[1].data
        sens = sens_files[beam]
        #dl = np.diff(lam)[0]
        #if dl > 0:
        
        #sens_interp = np.interp(lam, sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)
        delta_wave = lam[1]-lam[0]
        if delta_wave > 0:
            sens_interp = unicorn.utils_c.interp_conserve_c(lam, sens['WAVELENGTH'], sens['SENSITIVITY']*1.e-17, left=0., right=0.)
        else:
            sens_interp = unicorn.utils_c.interp_conserve_c(lam[::-1], sens['WAVELENGTH'], sens['SENSITIVITY']*1.e-17, left=0., right=0.)[::-1]
            
        #else:
        #    sens_interp = np.interp(lam[::-1], sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)[::-1]
            
        #sens_interp = utils_c.interp_c(lam, np.array(sens.field('WAVELENGTH'), dtype=np.float64), np.array(sens.field('SENSITIVITY'), dtype=np.float64)*1.e-17, extrapolate=0.)
        
        #### Sensitivity curve is defined "per A" not "per pixel" so multiply by A / pixel and 
        #### the grow_factor term accounts for the fact that the pixels are smaller in the 
        #### spatial axis as well in the interlaced images
        #print '%s DLAM: %.3f %.3f' %(beam, np.median(np.diff(lam)), np.std(np.diff(lam)))
        sens_interp *= np.median(np.abs(np.diff(lam)))/(growx*growy)
        
        if xarr[keep].size > 1:
            full_sens[y0[keep]+NY,xpix[keep]] += sens_interp[keep]
            full_sens[y0[keep]+NY+1,xpix[keep]] += sens_interp[keep]

        if (lam_spec is not None) & (flux_spec is not None):
            #spec_interp = threedhst.utils.interp_conserve(lam, lam_spec, flux_spec, left=0., right=0.)            
            if delta_wave > 0:
                spec_interp = utils_c.interp_conserve_c(lam, lam_spec, flux_spec, left=0., right=0.)
            else:
                spec_interp = utils_c.interp_conserve_c(lam[::-1], lam_spec, flux_spec, left=0., right=0.)[::-1]
                
            #### Smooth input model
            # xg = np.arange(-5*smooth_binned_sigma, 5.1*smooth_binned_sigma, 1)
            # yg = np.exp(-xg**2/2/smooth_binned_sigma**2)
            # yg /= np.sum(yg)
            # spec_interp = nd.convolve1d(spec_interp, yg)

            ####### Debugging
            #print BEAM, spec_interp
            #spec_interp = np.interp(lam, lam_spec, flux_spec) 
            # plt.plot(lam, spec_interp, color='red')
            # plt.plot(lam_spec, flux_spec, color='blue')
            # print beam, dl, lam, spec_interp
            # np.savetxt('/tmp/lam',lam)
            # if beam == 'A':
            #     np.savetxt('/tmp/lam',np.array([lam, spec_interp]).T)
            #     np.savetxt('/tmp/lam_spec',np.array([lam_spec, flux_spec]).T)
            #     print 'xxx', lam, spec_interp, lam_spec, flux_spec
            #     
            #     if 0:
            #         lam, spec_interp = np.loadtxt('/tmp/lam', unpack=True)
            #         lam_spec, flux_spec = np.loadtxt('/tmp/lam_spec', unpack=True)
                    
            spec_interp[~np.isfinite(spec_interp)] = 0
            spec_interp[lam < np.min(lam_spec)] = 0
            
            sens_interp *= spec_interp
        
        #### Increase response in the zeroth order to match observations
        if (beam == 'B') & (grism == 'G141'):
            sens_interp *= 1.4 
        
        if (beam == 'B') & (grism == 'G102'):
            sens_interp *= 0.9 #3.6 * 10
        
        # if beam in ['C','D','E']:
        #     sens_interp *= 0.25
                
        stripe *= np.dot(np.ones((NY*2+1,1)), sens_interp.reshape(1,NX)) * (growx*growy)
        
        model += stripe
        #print beam, model.max()
        
    if 'A' not in BEAMS:
        yoff_array = np.array([0,0])
        
    return model, (xmi, xma, wavelength, full_sens, yoff_array, beam_index)
    
def process_GrismModel(root='GOODS-S-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=27, REFINE_MAG_LIMIT=23,  make_zeroth_model=True, use_segm=False, model_slope=0, model_list = None, normalize_model_list=True, direct='F140W', grism='G141', BEAMS=['A','B','C','D'], old_filenames=False, align_reference=True):
    import unicorn.reduce
    
    model = unicorn.reduce.GrismModel(root=root, grow_factor=grow_factor, growx=growx, growy=growy, MAG_LIMIT=MAG_LIMIT, use_segm=use_segm, direct=direct, grism=grism, old_filenames=old_filenames)
    
    if align_reference:
        if 'FINEX' not in model.im[0].header.keys():
            model.align_direct_reference(match_threshold=0.5, mag_limit=23., r_limit=8, apply=True)
        
    model.model_list = model_list
    
    if not use_segm:
        #model.get_corrected_wcs(verbose=True)
        pass
        
    status = model.load_model_spectra()
    if not status:
        if make_zeroth_model:
            ### "zeroth" iteration to flag 0th order contamination, no color / norm iteration
            model.compute_full_model(refine=False, MAG_LIMIT=MAG_LIMIT, save_pickle=True, BEAMS=['B'], normalize_model_list=normalize_model_list)   
            ### Reset the model
            model.init_object_spectra()
            model.model*=0
        
        ### First iteration with flat spectra and the object flux
        model.compute_full_model(refine=False, MAG_LIMIT=MAG_LIMIT, save_pickle=False, model_slope=model_slope, BEAMS=BEAMS, normalize_model_list=normalize_model_list)   
        ### For the brighter galaxies, refine the model with the observed spectrum         
        model.compute_full_model(refine=True, MAG_LIMIT=REFINE_MAG_LIMIT, save_pickle=True, model_slope=model_slope, BEAMS=BEAMS)
        
        ### Make zeroth-order region file
        model.make_zeroth_regions(MAG_LIMIT=23, scale_size=True, id_label=True)
        
    return model
    
class Interlace1D():
    def __init__(self, file='UDS-17_00319.1D.fits', PNG=True, flux_units=True, grism='G141'):
        """
        Wrapper for Interlaced 1D spectra.
        
        For convenience, FITS table columns are aliased to 
        
        self.lam = self.data['WAVE']
        self.flux = self.data['FLUX']
        self.error = self.data['ERROR']
        self.contam = self.data['CONTAM']
        self.sens = self.data['SENSITIVITY']
        
        The flux, error, and contam columns are divided by `sens` 
        automatically if `flux_units==True`.
        
        """

        self.file = file
        tab = pyfits.open(file)
        self.header = tab[1].header
        self.data = tab[1].data
        tab.close()
        
        self.sens = self.data['SENSITIVITY']
        self.N = self.sens.shape[0]
        
        if flux_units:
            convert = self.sens
        else:
            convert = np.ones(self.N)
            
        self.lam = self.data['WAVE']
        self.flux = self.data['FLUX']/convert
        self.error = self.data['ERROR']/convert
        self.contam = self.data['CONTAM']/convert
        self.grism = grism
                
        if PNG:
            self.show(savePNG=True)
    
    def find_em_lines(self, fp=None, verbose=True, ascii_file=True):
        """
        Search for emission lines with a wavelet search method modeled
        after SDSS.  Note line centers are just the center pixel of
        the line and aren't fit.  Signal-to-noise comes from the median
        flux error of the spectrum.

        Input variable `fp` is an optional file pointer to save the results.
        
        If fp is None, save to a separate ASCII file (in './'):
        
            os.path.basename(self.file).replace('.fits', '.wavelet.dat')
            
        """
        
        root = os.path.basename(self.file).split('.1D')[0]
        lines = threedhst.spec1d.findLines(self, idx=-1, show=False, verbose=False, trim_abs=True)
        line_info = '%s  %10.6f %10.6f %6.3f' %(root, self.header['RA'], self.header['Dec'], self.header['MAG'])
        
        if lines is None:
            lines = []
            
        result = []
        for line in lines:
            line_info += '   %8.1f %5.2f %5.2f %5.1f %4d %4.2f' %(line.wave, line.height, line.continuum, line.sn, line.waveMax-line.waveMin, line.fcontam)        
            result.extend([line.wave, line.height, line.continuum, line.sn, line.waveMax-line.waveMin, line.fcontam])
        
        if verbose:
            print line_info
            
        if fp:
            fp.write(line_info+'\n')
            
        if (fp is None) & ascii_file:
            fp = open(root+'.1D.wavelet.dat','w')
            fp.write('# id ra dec mag  (wave height cont. S/N  width  fcontam)\n')
            fp.write(line_info+'\n')
            fp.close()
            
        columns = ['wave','height','cont','sn','width','fcontam']
        result_dict = {}
        for i in range(len(columns)):
            result_dict[columns[i]] = np.array(result[i::len(columns)])
        
        return result_dict
        
    def show(self, savePNG=True):
        import unicorn
        
        #fig = unicorn.catalogs.plot_init(xs=4, left=0.1, right=0.02*1.618, bottom=0.08, top=0.02, aspect=1./1.618)
        fig = Figure(figsize=[4,4/1.618], dpi=100)

        fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.135,
                            bottom=0.15,right=0.97,top=1-0.02*1.618)
        
        ax = fig.add_subplot(111)
        
        if self.grism == 'G102':
            wavelength_region = (self.data.wave > 0.78e4) & (self.data.wave < 1.15e4)
        else:
            wavelength_region = (self.data.wave > 1.15e4) & (self.data.wave < 1.65e4)
        if np.sum(wavelength_region) > 0:
            ymax = ((self.data.flux-self.data.contam)/self.data.sensitivity)[wavelength_region].max()
            if ymax <= 0:
                ymax = (self.data.flux/self.data.sensitivity)[wavelength_region].max()
        else:
            ymax=0
            
        if ymax <= 0:
            ymax = (self.data.flux/self.data.sensitivity).max()
        if ymax <= 0:
            ymax = 1
        
        ax.errorbar(self.data.wave, (self.data.flux-self.data.contam)/self.data.sensitivity, self.data.error/self.data.sensitivity, color='black', ecolor='0.7', alpha=0.5, marker=',', ms=3, linestyle='None')
        ax.plot(self.data.wave, (self.data.flux-self.data.contam)/self.data.sensitivity, color='blue', alpha=0.5)
        ax.plot(self.data.wave, self.data.flux/self.data.sensitivity, color='red', alpha=0.2)
        ax.plot(self.data.wave, self.data.contam/self.data.sensitivity, color='red')
        
        ax.set_ylim(-0.1*ymax, 1.1*ymax)
        if self.grism == 'G102':
            ax.set_xlim(0.74e4, 1.17e4)
        else:
            ax.set_xlim(1.05e4, 1.72e4)
        ax.set_xlabel(r'$\lambda$ / $\AA$')
        ax.set_ylabel(r'$f_\lambda / 10^{-17}$ cgs')
        
        if savePNG:
            #fig.savefig(self.file.replace('fits','png'))
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(os.path.basename(self.file).replace('fits','png'), dpi=100, transparent=False)

def get_all_emission_lines(field='GOODS-N', force=True, skip=True):
    """
    Run the new wavelet line finder on *all* 1D spectra in the v2.0 release to look for
    emission line galaxies at mags fainter than the current extraction limit.
    
    (June 13, 2012)
    
    """
    import unicorn
    PWD = os.getcwd()
    
    BASE_PATH = unicorn.GRISM_HOME + '../Release/v2.0/%s' %(field)
    
    os.chdir(BASE_PATH)
    files = glob.glob(field+'-*')
    dirs = []
    for file in files:
        if os.path.isdir(file):
            dirs.append(file)
        
    for dir in dirs:
        if not os.path.exists(os.path.join(BASE_PATH, dir, '1D/FITS')):
            continue
        #
        try:
            os.mkdir(os.path.join(PWD, dir,'1D/Wavelet'))
        except:
            if force:
                ## Ignore couldn't make dir, e.g., dir already exists
                pass
            else:
                ## Don't force, e.g. maybe write permission problems in PWD
                continue
        #
        os.chdir(os.path.join(PWD, dir, '1D/Wavelet'))
        #
        oned_files = glob.glob(os.path.join(BASE_PATH, dir, '1D/FITS/*1D.fits'))
        for file in oned_files:
            if os.path.exists(os.path.join(PWD, dir, '1D/Wavelet', os.path.basename(file.replace('fits','wavelet.dat')))) & skip:
                continue
            #
            oned = unicorn.reduce.Interlace1D(file, PNG=False, flux_units=True)
            lines = oned.find_em_lines(verbose=True, ascii_file=True)
            
    #
    os.chdir(PWD)
           
class Interlace2D():
    def __init__(self, file='GOODS-S-34_00446.2D.fits', PNG=True, growx=2, growy=2, ngrowx=0, ngrowy=0):
        """
        
        Get 'NGROW' from the header if it exists, or use parameter.
        
        """
        self.id = int(file.split('.2D')[0].split('_')[-1])
        self.file = file
        self.im = pyfits.open(file)
            
        self.thumb = np.cast[np.double](self.im['DSCI'].data)
        self.thumb_weight = np.cast[np.double](self.im['DWHT'].data)
        self.seg = np.cast[np.uint](self.im['DSEG'].data)
        if 'GRISM' not in self.im[0].header.keys():
            self.grism_element = 'G141'
            self.direct_filter = 'F140W'
        else:
            self.grism_element = self.im[0].header['GRISM']
            self.direct_filter = self.im[0].header['FILTER']
        
        if self.grism_element == 'G800L':
            chip = self.im[0].header['CCDCHIP']
        else:
            chip=1
        
        if 'FINEX' in self.im[0].header.keys():
            self.fine_offset = [self.im[0].header['FINEX'], self.im[0].header['FINEY']]
        else:
            self.fine_offset = [0,0]
            
            
        self.profile2D = None
        
        #### Make sure to use the same configuration file that was used
        #### to generate the 2D fits file
        if 'GRISCONF' in self.im[0].header.keys():
            if unicorn.reduce.conf_file != self.im[0].header['GRISCONF']:
                unicorn.reduce.set_grism_config(grism=self.grism_element, force=True, use_config_file=self.im[0].header['GRISCONF'])
            
        unicorn.reduce.set_grism_config(self.grism_element, chip=chip)
        self.conf = unicorn.reduce.conf
        self.sens_files = unicorn.reduce.sens_files
        
        self.oned = unicorn.reduce.Interlace1D(file.replace('2D','1D'), PNG=False)
                    
        self.sh = self.im['SCI'].data.shape
        
        if 'WMIN' in self.im[0].header.keys():
            self.wlim = [self.im[0].header['WMIN'], self.im[0].header['WMAX'], self.im[0].header['DW'], self.im[0].header['W0']]
        else:
            self.wlim = grism_wlimit[self.grism_element]
            
        self.x_pix = self.im[0].header['X_PIX']
        self.y_pix = self.im[0].header['Y_PIX']
        self.pad = self.im[0].header['PAD']
        self.grow_factor = 2 ### add these parameters to 2D header
        try:
            self.growx = self.im[0].header['GROWX']
            self.growy = self.im[0].header['GROWY']
        except:
            self.growx = 2
            self.growy = 2
        
        if 'NGROWX' in self.im[0].header.keys():
            self.ngrowx = self.im[0].header['NGROWX']
        else:
            self.ngrowx = ngrowx
            
        if 'NGROWY' in self.im[0].header.keys():
            self.ngrowy = self.im[0].header['NGROWY']
        else:
            self.ngrowy = ngrowy
        
        ### Back-compatibility before NGROW split into X/Y
        if 'NGROW' in self.im[0].header.keys():
            self.ngrowx = self.ngrowy = self.im[0].header['NGROW']
            
        #
        #### XXX need to get the zeropoint for arbitrary filter
        self.flux = self.thumb * 10**(-0.4*(ZPs[self.direct_filter]+48.6))* 3.e18 / PLAMs[self.direct_filter]**2 / 1.e-17
        self.total_flux = np.sum(self.flux*(self.seg == self.id))
        #self.total_flux = np.sum(self.flux*(self.seg > 0))
        
        self.init_model()
        #self.model_quality(width=5)
        self.twod_mask = np.isfinite(self.im['SCI'].data) & (self.im['WHT'].data > 0)
        self.cleaned = self.im['SCI'].data-self.im['CONTAM'].data
        
    def init_model(self, lam_spec=None, flux_spec=None):
        """
        Initialize all of the junk needed to go from the pixels in the 
        direct thumbnail to the 2D model spectrum
        """
        orders, self.xi = unicorn.reduce.grism_model(self.x_pix, self.y_pix, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=['A'], grow_factor=self.grow_factor, growx=self.growx, growy=self.growy, pad=self.pad, ngrowx=self.ngrowx, ngrowy=self.ngrowy, grism=self.grism_element, conf=self.conf, sens_files=self.sens_files, finex=self.fine_offset[0], finey=self.fine_offset[1])
        
        yord, xord = np.indices(orders.shape)
        beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), self.xi[5].reshape((1,-1)))
        
        non_zero = orders > 0
        
        cast = np.int
        xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), self.xi[2][non_zero], self.xi[3][non_zero], beams[non_zero]
                
        ys = orders.shape
        xord += self.xi[0]
        yord -= (ys[0]-1)/2
        #        
        sh = self.thumb.shape
        NX = xord.max()+sh[0]
        ypad = 10
        self.object = np.zeros((sh[0]+2*ypad, NX), dtype=np.double)
        self.obj_flux = np.zeros((sh[0]+2*ypad, NX), dtype=np.double)
        self.obj_flux[ypad:-ypad,0:sh[0]] = self.flux
        self.obj_seg = np.zeros((sh[0]+2*ypad, NX), dtype=np.uint)
        self.obj_seg[ypad:-ypad,0:sh[0]] = self.seg
        self.ypad = ypad
        
        #### Don't consider pixels with little signal
        # self.obj_seg[self.obj_flux/self.obj_flux.max() < 0.01] = 0
        # self.total_flux = np.sum(self.obj_flux[self.obj_seg > 0])
        
        status = utils_c.disperse_grism_object(self.obj_flux, self.id, self.obj_seg, xord, yord, ford, self.object, xpix=[0, self.thumb.shape[0]], ypix=[self.ypad, self.thumb.shape[0]+self.ypad])
        
        xxi = np.int(np.round(sh[0]/2.))+xord
        use = (xxi >= 0) & (xxi < (NX))
        ### First order only
        #use = use & (xord > 10*self.grow_factor) & (xord < 213*self.grow_factor) 
        use_order = use 
        cut = np.zeros(NX)
        cut[xxi[use_order]] = word[use_order] 
        self.object_wave = cut.copy() #np.dot(np.ones((self.sh[0],1)),
        
        xarr = np.arange(NX)
        # wavelength_region = (self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)
        # if self.grism_element == 'G102':
        #     wavelength_region = (self.object_wave >= 0.73e4) & (self.object_wave <= 1.18e4)
        # #
        # if self.grism_element == 'G800L':
        #     wavelength_region = (self.object_wave >= 0.58e4) & (self.object_wave <= 0.92e4)
        #wlim = grism_wlimit[self.grism_element]
        wlim = self.wlim
        wavelength_region = (self.object_wave >= wlim[0]) & (self.object_wave <= wlim[1])
        
        limited = wavelength_region & ((self.object_wave-self.im['WAVE'].data.min()) > -1) & ((self.object_wave - self.im['WAVE'].data.max()) < 1)        
        
        xmin = xarr[limited].min()
        xmax = xarr[limited].max()
        ### need to account for [xmin:xmax] = (xmax-xmin)-1
        if limited.sum() < wavelength_region.sum():
            xmax += 1
        #
        self.xmin, self.xmax = xmin, xmax
        
        ### Do test better
        imw = self.im['WAVE'].data
        limited = (self.object_wave >= (imw[imw > 0].min()-1)) & (self.object_wave <= (imw[imw > 0].max()+1))
        xmin, xmax = xarr[limited][[0,-1]] #-1
        xmax += 1
        self.xmin, self.xmax = xmin, xmax
        #print 'XMINMAX: ', self.xmin, self.xmax
        
        #print 'Interlace2D xminmax: %d %d' %(xmin, xmax)
        
        # if 'DXMIN' in self.im[0].header.keys():
        #     print 'XXX DXMIN'
        #     self.xmin = self.im[0].header['DXMIN']
        #     self.xmax = self.im[0].header['DXMAX']

        # xoff_arr = np.arange(len(self.xi[4]), dtype=np.double)+self.xi[0]+sh[0]/2.
        # xoff_int = np.arange(xmin, xmax, dtype=np.double)
        # yoff_int = utils_c.interp_c(xoff_int, xoff_arr, self.xi[4])        
        #self.yc = np.int(np.round(yoff_int.mean()))
        self.yc = -self.im['SCI'].header['YTOFF']
        
        self.xord = xord
        self.yord = yord
        self.ford = ford
        self.non_zero = non_zero
        
        self.model = self.object[self.yc+ypad:self.yc-ypad,xmin:xmax]
    
    def compute_model(self, lam_spec=None, flux_spec=None, smooth_binned_sigma=0.8):
        """
        The nuts and bolts were determined in init_model.  Now just provide 
        and input spectrum and compute the model
        """
        ### Dispersion of single pixel        
        # if (lam_spec is not None):
        #     if (lam_spec.min() > 0.9e4) | (lam_spec.max() < 1.8e4):
        #         xint = np.arange(0.9e4,1.9e4,22)
        #         yint = xint*0.+1
        #         yint = unicorn.utils_c.interp_conserve_c(xint, lam_spec, flux_spec)
        #         lam_spec, flux_spec = xint, np.maximum(0,yint)
        #         print 'Reform', lam_spec, flux_spec
                
        #print lam_spec, flux_spec
        
        orders, self.xi = unicorn.reduce.grism_model(self.x_pix, self.y_pix, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=['A'], grow_factor=self.grow_factor, growx=self.growx, growy=self.growy, pad=self.pad, ngrowx=self.ngrowx, ngrowy=self.ngrowy, grism=self.grism_element, smooth_binned_sigma=smooth_binned_sigma, conf=self.conf, sens_files=self.sens_files, finex=self.fine_offset[0], finey=self.fine_offset[1])
        
        #print orders.shape
        
        if (lam_spec is not None):
            ford = np.array(orders[self.non_zero], dtype=np.float64)
        else:
            ford = self.ford
            
        ### Integrated over the thumbnail
        self.object *= 0.
        status = utils_c.disperse_grism_object(self.obj_flux, self.id, self.obj_seg, self.xord, self.yord, ford, self.object, xpix=[0, self.thumb.shape[0]], ypix=[self.ypad, self.thumb.shape[0]+self.ypad])
                
        ### Done!
        self.model = self.object[self.yc+self.ypad:self.yc-self.ypad,self.xmin:self.xmax]
    
    def optimal_extract(self, input=None, force_init=False, fcontam=0, get_var=False):
        
        if (self.profile2D is None) | (force_init):            
            #### Get extraction window where profile is > 10% of maximum 
            osh = self.im['MODEL'].data.shape
            xarr = np.arange(osh[1])
            # xint14 = int(np.interp(1.3e4, self.im['WAVE'].data, xarr))
            # if self.grism_element == 'G102':
            #     xint14 = int(np.interp(0.98e4, self.im['WAVE'].data, xarr))
            # #
            # if self.grism_element == 'G800L':
            #     xint14 = int(np.interp(0.75e4, self.im['WAVE'].data, xarr))
            # 
            #wlim = grism_wlimit[self.grism_element]
            wlim = self.wlim
            xint14 = int(np.interp(wlim[3], self.im['WAVE'].data, xarr))

            yprof = np.arange(osh[0])
            profile = np.sum(self.im['MODEL'].data[:,xint14-10:xint14+10], axis=1)
            profile = profile/profile.sum()

            #profile2D = np.dot(profile.reshape((-1,1)), np.ones((1,osh[1])))
            window = profile >= 0.1*profile.max()

            #### Use object as the profile since it contains the y shift along the trace
            prof_x = np.sum(self.im['MODEL'].data, axis=0)
            self.profile2D = self.im['MODEL'].data/np.dot(np.ones((osh[0],1)), prof_x.reshape((1,-1)))
        
        
        obj_cleaned = input*1
        #variance = self.im['WHT'].data**2
                        
        #### Optimal extraction
        opt_variance = self.im['WHT'].data**2 + (fcontam*self.im['CONTAM'].data)**2
        obj_cleaned[self.im['WHT'].data == 0] = 0
        opt_variance[self.im['WHT'].data == 0] = 1.e6
        
        num = np.sum(self.profile2D*obj_cleaned/opt_variance, axis=0)        
        denom = np.sum(self.profile2D**2 / opt_variance, axis=0)

        optimal_sum = num / denom
        
        if get_var:
            optimal_var = 1./np.sum(self.profile2D**2/opt_variance, axis=0)
            return self.im['WAVE'].data, optimal_sum, optimal_var
        else:   
            return self.im['WAVE'].data, optimal_sum
    
    def trace_extract(self, input=None, dy=0, width=0, get_apcorr=False, get_mask=False):
        """
        Extract pixel values along the trace
        
        With `get_apcorr`, also return an aperture correction computed from
        the ratio of the total model flux to the model flux along the trace
        as specified.
        
        """
        if input is None:
            obj_cleaned = self.im['SCI'].data*1
        else:
            obj_cleaned = input
                
        ytrace = np.cast[np.float64](self.im['YTRACE'].data)+dy
        trace_pix = np.cast[int](np.round(ytrace))
        #trace_spec = self.grism_sci[trace_pix,:]
        trace_spec = ytrace*0.
        
        #### Extract model in same way
        obj_model = self.im['MODEL'].data*1.
        obj_model[obj_cleaned == 0] = 0.
        trace_model = trace_spec*0.
        
        if width == 0:
            trace_lo = trace_pix
            trace_hi = trace_pix+1
        else:
            trace_lo = trace_pix-width
            trace_hi = trace_pix+width
        
        trace2D = obj_cleaned*0
        for i in range(len(ytrace)):
            trace2D[trace_lo[i]:trace_hi[i],i] = 1
            #trace_spec[i] = obj_cleaned[trace_pix[i],i]
            #trace_spec[i] = np.sum(obj_cleaned[trace_lo[i]:trace_hi[i],i])
        
        if get_mask:
            return trace2D
            
        trace_spec = np.sum(obj_cleaned*trace2D, axis=0)
        
                    
        if get_apcorr:
            #### Extract model in same way to get aperture correction
            obj_model = self.im['MODEL'].data*1.
            obj_model[obj_cleaned == 0] = 0.
            # trace_model = trace_spec*0.
            # for i in range(len(ytrace)):
            #     trace_model[i] = np.sum(obj_model[trace_lo[i]:trace_hi[i],i])
            trace_model = np.sum(obj_model*trace2D, axis=0)
            apcorr = obj_model.sum(axis=0) / trace_model
            return self.im['WAVE'].data, trace_spec, apcorr
        
        else:
            return self.im['WAVE'].data, trace_spec
            
    def oned_spectrum(self, verbose=True):
        """
        1D extraction after generating the 2D spectra.
        
        Take this out of the full model object and get everything from 
        the 2D spectrum itself
        """

        #### Paremters from the full model
        self.grism_sci = np.cast[np.float64](self.im['SCI'].data)
        self.full_model = np.cast[np.float64](self.im['CONTAM'].data + self.im['MODEL'].data)
        self.the_object = np.cast[np.float64](self.im['MODEL'].data)
        self.grism_wht = np.cast[np.float64](self.im['WHT'].data)
        self.wave = np.cast[np.float64](self.im['WAVE'].data)
        self.ytrace = np.cast[np.float64](self.im['YTRACE'].data)
        self.sens = np.cast[np.float64](self.im['SENS'].data)
        self.root = self.file.split('_%05d' %(self.id))[0]
        
        t0 = time.time()
        
        obj_cleaned = self.grism_sci - (self.full_model - self.the_object)
        contam_cleaned = self.grism_sci - self.the_object
        variance = self.grism_wht**2
            
        #### Get extraction window where profile is > 10% of maximum 
        osh = self.the_object.shape
        xarr = np.arange(osh[1])
        # xint14 = int(np.interp(1.3e4, self.wave, xarr))
        # if self.grism_element == 'G102':
        #     xint14 = int(np.interp(0.98e4, self.wave, xarr))
        # #
        # if self.grism_element == 'G800L':
        #     xint14 = int(np.interp(0.75e4, self.wave, xarr))
        #wlim = grism_wlimit[self.grism_element]
        wlim = self.wlim
        xint14 = int(np.interp(wlim[3], self.wave, xarr))
        
        yprof = np.arange(osh[0])
        profile = np.sum(self.the_object[:,xint14-10:xint14+10], axis=1)
        profile = profile/profile.sum()
        
        #profile2D = np.dot(profile.reshape((-1,1)), np.ones((1,osh[1])))
        window = profile >= 0.1*profile.max()
        
        #### Use object as the profile since it contains the y shift along the trace
        prof_x = np.sum(self.the_object, axis=0)
        profile2D = self.the_object/np.dot(np.ones((osh[0],1)), prof_x.reshape((1,-1)))
        
        if (1 == 0):
            plt.plot(yprof, profile, color='blue')
            plt.plot(yprof[window], profile[window], color='red', linewidth=3)
            print np.trapz(profile[window], xprof[window]) / np.trapz(profile, xprof)
        
        obj_cleaned[variance == 0] = 0
        contam_cleaned[variance == 0] = 0
        
        #### Simple sums, also weighted
        weights = 1./variance
        weights[variance == 0] = 0
        
        var_sum = np.sum(weights[window,:], axis=0)
        
        weighted_sum = np.sum((obj_cleaned*weights)[window,:], axis=0) / var_sum
        weighted_var = 1./np.sum(weights[window,:], axis=0)
        weighted_sig = np.sqrt(weighted_var)
        
        simple_sum = np.sum(obj_cleaned[window,:], axis=0)
        simple_var = np.sum(variance[window,:], axis=0)
        simple_sig = np.sqrt(simple_var)
        
        #### Optimal extraction
        opt_variance = variance.copy()
        opt_variance[opt_variance == 0] = 1.e6
        
        num = np.sum(profile2D*obj_cleaned/opt_variance, axis=0)
        #num_contam = np.sum(profile2D*contam_cleaned/opt_variance, axis=0)
        num_contam = np.sum(profile2D*(self.full_model - self.the_object)/opt_variance, axis=0)
        num_full = np.sum(profile2D*self.grism_sci/opt_variance, axis=0)
        
        denom = np.sum(profile2D**2 / opt_variance, axis=0)
        
        optimal_sum = num / denom
        optimal_sum_contam = num_contam / denom
        optimal_sum_full = num_full / denom
        
        optimal_var = 1./np.sum(profile2D**2/opt_variance, axis=0)
        optimal_sig = np.sqrt(optimal_var)
        
        trace_pix = np.cast[int](np.round(self.ytrace))
        #trace_spec = self.grism_sci[trace_pix,:]
        trace_spec = optimal_sum*0.
        trace_sig = optimal_sum*0.
        for i in range(osh[1]):
            trace_spec[i] = self.grism_sci[trace_pix[i],i]
            trace_sig[i] = self.grism_wht[trace_pix[i],i]
        
        scale_to_total = 1./np.max(profile)       
        if not np.isfinite(scale_to_total):
            scale_to_total=-1
            
        c1 = pyfits.Column(name='wave', format='D', unit='ANGSTROMS', array=self.wave)
        c2 = pyfits.Column(name='flux', format='D', unit='ELECTRONS/S', array=optimal_sum_full)
        c3 = pyfits.Column(name='error', format='D', unit='ELECTRONS/S', array=optimal_sig)
        c4 = pyfits.Column(name='contam', format='D', unit='ELECTRONS/S', array=optimal_sum_contam)
        c5 = pyfits.Column(name='trace', format='D', unit='ELECTRONS/S', array=trace_spec)
        c6 = pyfits.Column(name='etrace', format='D', unit='ELECTRONS/S', array=trace_sig)
        c7 = pyfits.Column(name='sensitivity', format='D', unit='E/S / 1E-17 CGS', array=self.sens*(self.growx*self.growy))
        #print 'MAX SENS: %.f' %(self.sens.max())
        
        coldefs = pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7])
        head = pyfits.Header()
        head.update('ATRACE', scale_to_total, comment='Factor to scale trace to total flux')
        
        #ii = np.where(np.cast[int](self.cat.id) == self.id)[0][0]
        try:
            head.update('RA', self.im[0].header['RA'], comment='Target R.A.')
            head.update('DEC', self.im[0].header['DEC'], comment='Target Dec.')
        except:
            head.update('RA', 0., comment='Target R.A.')
            head.update('DEC', 0., comment='Target Dec.')
        #
        head.update('X_PIX', self.x_pix, comment='X pixel in interlaced image')
        head.update('Y_PIX', self.y_pix, comment='Y pixel in interlaced image')
        head.update('MAG', self.im[0].header['MAG'], comment='MAG_AUTO from interlaced catalog')
        
        tbHDU = pyfits.new_table(coldefs, header=head)
        tbHDU.writeto(self.root+'_%05d.1D.fits' %(self.id), clobber=True)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print '1D spectrum (%.3f)' %(dt)
    
    def model_quality(self, dy=0, width=5):
        """
        Check distribution of flux - model - contam pixels to identify 
        bad contamination
        """
        yi, xi = np.indices(self.im['DSCI'].data.shape)
        xpix = xi[self.im['DSEG'].data == self.id]
        #width = (xpix.max()-xpix.min())/3
        
        ivar = 1./self.im['WHT'].data**2
        ivar[self.im['WHT'].data == 0] = 0
        
        trace_mask = (self.trace_extract(input=self.im['SCI'].data, dy=dy, width=width, get_apcorr=False, get_mask=True) > 0) & (self.im['SCI'].data != 0) & (self.im['WHT'].data != 0)
        if trace_mask.sum() < 10:
            self.resid_stats = np.ones(5)*99
        else:
            cleaned = self.im['SCI'].data - self.im['MODEL'].data - self.im['CONTAM'].data
            resid_sigma = cleaned / self.im['WHT'].data
            self.resid_stats = np.percentile(resid_sigma[trace_mask], [2.5,16,50,84,97.5])
          
    def thumb_image(self):
        ok = self.im['DSEG'].data == self.id
        max = self.im['DSCI'].data[ok].max()
        sh = ok.shape
        sh2 = self.im['SCI'].data.shape
        
        xarr, yarr = np.arange(sh[1]), np.arange(sh[0])
        xsec = (np.arange(sh[1]) - sh[1]/2. + ((sh[1]+1) % 2)*0.5)*0.128254/self.im[0].header['GROWX']
        ysec = (np.arange(sh[0]) - self.im['YTRACE'].data[sh2[1]/2.]+1)*0.128254/self.im[0].header['GROWY']
        
        NY, NX = int(np.max(np.abs(ysec)/0.5)), int(np.max(np.abs(xsec)/0.5))
        yint, xint = np.interp(np.arange(-NY*0.5,NY*0.5+0.1,0.5), ysec, yarr), np.interp(np.arange(-NX*0.5,NX*0.5+0.1,0.5), xsec, xarr)
        
        fig = unicorn.plotting.plot_init(aspect=1./2, left=0.01, right=0.01, bottom=0.01, top=0.01, xs=3, square=True, NO_GUI=True)

        ax = fig.add_axes((0.01, 0.01, 0.485, 0.98))
        ax.imshow(0-self.im['DSCI'].data, interpolation='Nearest', aspect='auto', vmin=-1.1*max, vmax=0.1*max)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xticks(xint); ax.set_yticks(yint)

        ax = fig.add_axes((0.495, 0.01, 0.485, 0.98))
        ax.imshow(0-ok*1., interpolation='Nearest', aspect='auto', vmin=-1, vmax=0)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xticks(xint); ax.set_yticks(yint)
        ax.text(0.5,0.5, self.id, ha='center', va='center', color='black', transform=ax.transAxes, size=8)
        ax.text(0.5,0.5, self.id, ha='center', va='center', color='white', alpha=0.9, transform=ax.transAxes, size=8)
        
        unicorn.plotting.savefig(fig, self.file.replace('2D.fits', 'thumb.png'))
        
    def init_fast_model(self):
        """
        Init of thumbnail-based approach building a model up summing
        the thumbnail
        """        
        thumb = self.im['DSCI'].data*1
        thumb[self.im['DSEG'].data != self.id] = 0
        thumb *= self.total_flux/np.sum(thumb)
        
        ### Set up
        sh = self.sh
        thumb_matrix = np.zeros((sh[1], sh[0], sh[1]))
        NX = sh[0]/2
        ytrace = self.im['YTRACE'].data-NX#-0.2
        #ytrace = np.cast[int](np.round(ytrace))
        
        ### For resampling pixels
        flux_scale = self.growx*self.growy
        
        for i in range(1, NX):
            #print i
            thumb_matrix[i,:,:NX+i] = nd.shift(thumb, (ytrace[i], 0))[:,NX-i:]*self.im['SENS'].data[i]*flux_scale

        for i in range(NX,sh[1]-NX):
            #print i
            thumb_matrix[i,:,i-NX:i+NX] = nd.shift(thumb, (ytrace[i], 0))*self.im['SENS'].data[i]*flux_scale

        #
        # for i in range(0,sh[1],2): 
        #     sh_thumb = nd.shift(thumb, (ytrace[i], 0),order=3)
        #     ds9.view(sh_thumb*100)
        #     print i, sh_thumb.max()/thumb.max(), sh_thumb.sum()/thumb.sum()
            
        for ix, i in enumerate(range(sh[1]-NX,sh[1])):
            #print i
            thumb_matrix[i,:,sh[1]-sh[0]+ix:] = nd.shift(thumb, (ytrace[i], 0))[:,:sh[0]-ix]*self.im['SENS'].data[i]*flux_scale
        
        self.thumb_matrix = thumb_matrix.reshape(sh[1],-1)
        self.wave = np.cast[float](self.im['WAVE'].data)
        
    def fast_compute_model(self, lam_spec=None, flux_spec=None):
        try:
            flx = self.wave*0
        except:
            self.init_fast_model()
        
        if lam_spec is None:
            flux_interp = self.wave*0.+1
        else:
            flux_interp = unicorn.utils_c.interp_conserve_c(self.wave, lam_spec, flux_spec)
            
        mflat = np.dot(flux_interp, self.thumb_matrix)
        self.model = mflat.reshape(self.sh)
    
    def simple_normalization(self, lam_spec=None, flux_spec=None):
        """
        The total flux in the spectrum should be similar to that of the 
        flat-spectrum model.
        """
        self.compute_model(lam_spec=lam_spec, flux_spec=flux_spec)
        clean = self.im['SCI'].data-self.im['CONTAM'].data
        mask = self.model > 0.05*self.model.max()
        norm =  clean[mask].sum() / self.model[mask].sum()
        return norm
        
    def timer(self):
        """
        Test how long it takes to execute.
        """
        import time
        t0 = time.time()
        N = 100
        m = self.model*1.
        for i in range(N):
            self.compute_model()
            m += self.model
        t1 = time.time()
        print '%d steps, no model: %.3f s' %(N, t1-t0)
        
        wave = np.exp(np.arange(7,10,0.1))
        flux = wave*0.+1

        wave = np.cast[np.double](self.oned.data.wave)
        flux = np.cast[np.double]((self.oned.data.flux-self.oned.data.contam) / self.oned.data.sensitivity)
        
        t0 = time.time()
        N = 100
        m = self.model*1.
        for i in range(N):
            self.compute_model(lam_spec = wave, flux_spec = flux)
            m += self.model
        t1 = time.time()
        print '%d steps, WITH model: %.3f s' %(N, t1-t0)
                
    def __add__(self, other, data=None, add_self=True):
        """
        Add the "model" attribute arrays of two `Interlace2D` objects with 
        the appropriate slices.  Requires the [X/Y]GRISM0 keywords in the 
        2D FITS file.
        """
        if not isinstance(other, self.__class__):
            threedhst.showMessage('Input is not a `%s` object (`%s`)' %(self.__class__, other.__class__), warn=True)
            return False
        
        if 'XGRISM0' not in self.im[0].header.keys():
            threedhst.showMessage('Self (%s) looks like an old \n2D FITS file without the necessary header keywords.  \n\nUpdate the module and extract again.' %(self.im.filename()), warn=True)
            return False
        
        if 'XGRISM0' not in other.im[0].header.keys():
            threedhst.showMessage('Self (%s) looks like an old \n2D FITS file without the necessary header keywords.  \n\nUpdate the module and extract again.' %(other.im.filename()), warn=True)
            return False
        
        n_self = self.model.shape[::-1]
        n_other = other.model.shape[::-1]
        
        xy_self = [self.im[0].header['XGRISM0'], self.im[0].header['YGRISM0']]
        xy_other = [other.im[0].header['XGRISM0'], other.im[0].header['YGRISM0']]
        
        xmin = np.min([xy_self[0], xy_other[0]])
        xmax = np.max([xy_self[0]+n_self[0], xy_other[0]+n_other[0]])
        
        ymin = np.min([xy_self[1], xy_other[1]])
        ymax = np.max([xy_self[1]+n_self[1], xy_other[1]+n_other[1]])
        
        out = np.zeros((ymax-ymin, xmax-xmin))
        if add_self:
            out[xy_self[1]-ymin:xy_self[1]-ymin+n_self[1], xy_self[0]-xmin:xy_self[0]-xmin+n_self[0]] += self.model

        if data is None:
            out[xy_other[1]-ymin:xy_other[1]-ymin+n_other[1], xy_other[0]-xmin:xy_other[0]-xmin+n_other[0]] += other.model
        else:
            out[xy_other[1]-ymin:xy_other[1]-ymin+n_other[1], xy_other[0]-xmin:xy_other[0]-xmin+n_other[0]] += data
            
        return out[xy_self[1]-ymin:xy_self[1]-ymin+n_self[1], xy_self[0]-xmin:xy_self[0]-xmin+n_self[0]]
    
    def get_sensitivity_error(self, scale=1./20):
        """
        Get a component to the counts/s uncertainty from the flux correction
        """
        sens = self.sens_files['A']
        s_max = sens['SENSITIVITY'].max()
        
        oflux = self.obj_flux*1.        
        self.obj_flux = np.abs(oflux)#/np.abs(self.oflux).sum()
        self.compute_model(sens['WAVELENGTH'], ((s_max/sens['SENSITIVITY'])**2-0)*scale)
        self.obj_flux = oflux*1.
        
        self.sens_var = self.model**2
    
class GrismModel():
    def __init__(self, root='GOODS-S-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, use_segm=False, grism='G141', direct='F140W', output_path='./', old_filenames=False):
        """
        Initialize: set padding, growth factor, read input files.
        """
        self.root=root
        if old_filenames:
            self.baseroot = os.path.join(output_path, os.path.basename(root))
        else:
            self.baseroot = os.path.join(output_path, os.path.basename(root)) + '-%s' %(grism)
            
        self.grow_factor = grow_factor
        self.growx = growx
        self.growy = growy
        self.use_segm = use_segm
        
        #### Read the direct interlaced image        
        self.direct = pyfits.open(self.root+'-%s_inter.fits' %(direct))
        self.direct[1].data = np.cast[np.double](self.direct[1].data)
        
        self.filter = self.direct[0].header['FILTER']
        if 'PAD' in self.direct[1].header.keys():
            self.pad = self.direct[1].header['PAD']
        else:
            self.pad = 0.
            
        self.REF_INTER = False
        
        #### If a "reference" interlaced image exists, use that
        if os.path.exists(self.baseroot+'_ref_inter.fits'):
            self.REF_INTER = True
            self.im = pyfits.open(self.baseroot+'_ref_inter.fits')
            #self.pad += 4*self.im[1].header['NGROW']
            self.filter = self.im[0].header['FILTER']
        else:
            self.im = self.direct
        
        self.fine_offset = [0, 0]
        if 'FINEX' in self.im[0].header.keys():
            self.fine_offset = [self.im[0].header['FINEX'], self.im[0].header['FINEY']]
            
        try:
            self.growx = self.im[0].header['GROWX']
            self.growy = self.im[0].header['GROWY']
        except:
            self.growx = growx
            self.growy = growy
        
        if 'NGROWX' not in self.direct[1].header.keys():
            self.ngrowx = (self.direct[1].data.shape[1] - 1014*self.growx - self.pad) / (self.growx*2) 
        else:
            self.ngrowx = self.direct[1].header['NGROWX']
            
        if 'NGROWY' not in self.direct[1].header.keys():
            self.ngrowy = (self.direct[1].data.shape[1] - 1014*self.growx - self.pad) / (self.growx*2) 
        else:
            self.ngrowy = self.direct[1].header['NGROWY']

        self.ZP = ZPs[self.filter]

        self.PLAM = PLAMs[self.filter]
                    
        self.make_flux()
        
        self.grism_element = grism
        self.direct_element = direct
        
        if self.grism_element == 'G800L':
            self.chip = self.direct[1].header['CCDCHIP']
        else:
            self.chip=1
            
        unicorn.reduce.set_grism_config(self.grism_element, chip=self.chip)
        #unicorn.reduce.set_grism_config(grism)
        
        self.gris = pyfits.open(self.root+'-%s_inter.fits' %(self.grism_element))
        if self.gris[1].data.shape != self.direct[1].data.shape:
            threedhst.showMessage('G141 and F140W images have different dimensions!\n')
        self.gris[1].data = np.array(self.gris[1].data, dtype=np.double)
        self.gris[2].data = np.array(self.gris[2].data, dtype=np.double)
        self.sh = self.im[1].data.shape
        
        if not os.path.exists(self.baseroot+'_inter_seg.fits'):
            threedhst.showMessage('Running SExtractor...')
            self.find_objects(MAG_LIMIT=MAG_LIMIT)
                
        threedhst.showMessage('Read FITS files and catalog')
        print 'Read files'
        self.read_files()
        
        print 'Init object spectra'
        self.init_object_spectra()
        
        self.make_total_flux()
        
        self.ra_wcs, self.dec_wcs = self.cat.ra, self.cat.dec
        self.model = np.zeros(self.sh, dtype=np.double)
        self.object = np.zeros(self.sh, dtype=np.double)
        #self.test_object = np.zeros(self.sh)
        self.yf, self.xf = np.indices(self.sh)
        
        #### Dict list of input spectra [ model_list[id] = (wave, flam) ]
        #### to use for the initial model
        self.model_list = None
        
        # if os.path.exists('%s_inter.cat.wcsfix' %(self.root)):
        #     self.get_corrected_wcs()
        
    def read_files(self):
        """
        Read FITS files, catalogs, and segmentation images for a pair of 
        interlaced exposures.
        """        
        self.cat = threedhst.sex.mySexCat('%s_inter.cat' %(self.baseroot))
        
        #self.trim_edge_objects(verbose=True)
        
        if 'X_FLT' in self.cat.column_names:
            self.cat.x_pix = np.cast[float](self.cat.X_FLT)
            self.cat.y_pix = np.cast[float](self.cat.Y_FLT)
        else:
            self.cat.x_pix = np.cast[float](self.cat.X_IMAGE)
            self.cat.y_pix = np.cast[float](self.cat.Y_IMAGE)
            
        self.cat.mag = np.cast[float](self.cat.MAG_AUTO)
        
        self.segm = pyfits.open('%s_inter_seg.fits' %(self.baseroot))
        self.segm[0].data[self.segm[0].data < 0] = 0
        self.segm[0].data[~np.isfinite(self.segm[0].data)] = 0
        self.segm[0].data = np.array(self.segm[0].data, dtype=np.uint)
        
        if self.use_segm:
            self.objects = np.unique(self.segm[0].data)[1:]
        else:
            self.objects = self.cat['NUMBER']
                 
    def make_flux(self):
        """
        Convert from e/s to physical fluxes using the ZP and pivot wavelength 
        for the direct image filter
        """
        self.flux = self.im[1].data * 10**(-0.4*(self.ZP+48.6))* 3.e18 / self.PLAM**2 / 1.e-17
        self.flux_fnu = self.im[1].data * 10**(-0.4*(self.ZP+48.6))
        
    def find_objects(self, MAG_LIMIT=23.5):
        #### Run SExtractor on the direct image, with the WHT 
        #### extension as a weight image
        se = threedhst.sex.SExtractor()
        
        ## Set the output parameters required for aXe 
        ## (stored in [threedhst source]/data/aXe.param) 
        se.aXeParams()
        
        ## XXX add test for user-defined .conv file
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = '%s_inter.cat' %(self.baseroot)
        se.options['CHECKIMAGE_NAME'] = '%s_inter_seg.fits' %(self.baseroot)
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_RMS'
        se.options['WEIGHT_IMAGE']    = '%s-%s_inter.fits[1]' %(self.root, self.direct_element)
        se.options['FILTER']    = 'Y'
        
        #### Detect thresholds (default = 1.5)
        se.options['DETECT_THRESH']    = '1' 
        se.options['ANALYSIS_THRESH']  = '1'
        se.options['MAG_ZEROPOINT'] = '%.2f' %(self.ZP)
        
        #### Run SExtractor
        status = se.sextractImage('%s-%s_inter.fits[0]' %(self.root, self.direct_element))
        self.cat = threedhst.sex.mySexCat('%s_inter.cat' %(self.baseroot))
        
        #### Trim faint sources
        mag = np.cast[float](self.cat.MAG_AUTO)
        q = np.where(mag > MAG_LIMIT)[0]
        if len(q) > 0:
            numbers = np.cast[int](self.cat.NUMBER)[q]
            self.cat.popItem(numbers)
        
        self.cat.write()
        
        #self.trim_edge_objects(verbose=True)
        
        #self.segm = pyfits.open(self.root+'_seg.fits')
        
        threedhst.sex.sexcatRegions('%s_inter.cat' %(self.baseroot), '%s_inter.reg' %(self.baseroot), format=1)
        
        # try:
        #     self.make_total_flux()
        # except:
        #     pass
        # self.init_object_spectra()
    
    def trim_edge_objects(self, verbose=True):
        """
        Trim objects from the SExtractor catalog that are in the ragged
        edge of the interlaced images.
        """
        
        x_edge = (self.cat.x_pix < self.pad/2+10) | (self.cat.x_pix > (self.sh[1]-self.pad/2-25))
        y_edge = (self.cat.y_pix < self.pad/2+10) | (self.cat.y_pix > (self.sh[0]-self.pad/2-25))
        
        q = x_edge | y_edge
        
        if q.sum() > 0:
            if verbose:
                print "Cutting %d likely fake and/or othwerwise problematic objects\n   at the image edges from the SExtractor catalog \n   (they'll still be in the segmentation image)." %(q.sum())
            numbers = np.cast[int](self.cat.NUMBER)[q]
            self.cat.popItem(numbers)
            self.cat.write()
            
    def get_corrected_wcs(self, verbose=True):
        """
        Use iraf.dither.tran to get the real, undistorted WCS coordinates of each object
        """
        from pyraf import iraf
        
        no = iraf.no
        yes = iraf.yes
        INDEF = iraf.INDEF
        
        from iraf import stsdas
        from iraf import dither
        import threedhst.catIO as catIO
        
        if os.path.exists('%s_inter.cat.wcsfix' %(self.baseroot)):
            wcs = catIO.Readfile('%s_inter.cat.wcsfix' %(self.baseroot))
            test = wcs.id == self.cat.id
            if np.sum(test) == len(self.cat.id):
                if verbose:
                    print 'Get corrected coordinates from %s_inter.cat.wcsfix' %(self.baseroot)
                self.ra_wcs = wcs.ra
                self.dec_wcs = wcs.dec
                return True
        
        if self.REF_INTER:
            self.ra_wcs, self.dec_wcs = self.cat.ra*1., self.cat.dec*1.       
            
            #### Write the wcsfix file
            fp = open('%s_inter.cat.wcsfix' %(self.baseroot),'w')
            fp.write('# id    mag   ra    dec\n')
            for i in range(self.cat.id.shape[0]):
                fp.write('%5d  %.3f  %.6f  %.6f\n' %(self.cat.id[i], self.cat.mag[i], self.ra_wcs[i], self.dec_wcs[i]))

            fp.close()

            if verbose:
                print 'Corrected coordinates: %s_inter.cat.wcsfix' %(self.baseroot)

            return True
        
        #### Else: use tran to get the pixel in the DRZ frame and then xy2rd to get coords        
        try:
            wcs = wcsutil.WCSObject(self.root+'-%s_drz.fits[1]' %(self.direct_element))
        except:
            wcs = wcsutil.WCSObject(self.root+'-%s_drz_sci.fits[0]' %(self.direct_element))
            
        NOBJ = len(self.cat.id)
        
        #### fix for different exposure structure of Barro G102 program
        if 'GDN' in self.root:
            filter = self.grism_element
        else:
            filter = self.direct_element
            
        asn = threedhst.utils.ASNFile(self.root+'-%s_asn.fits' %(filter))
        flt = asn.exposures[0]+'_flt.fits'
        fp = open('/tmp/%s.flt_xy' %(self.root),'w')
        for i in range(NOBJ):
            #1014*growy+pad+growy*2*NGROW
            xi = (self.cat.x_pix[i]-self.pad/2. - self.ngrowx*self.growx)/self.growx
            yi = (self.cat.y_pix[i]-self.pad/2. - self.ngrowy*self.growy)/self.growy
            fp.write('%.2f %.2f\n' %(xi,  yi))
        
        fp.close()
        
        print 'Running iraf.tran()'
        threedhst.process_grism.flprMulti()
        
        ### For some reason this dies occasionally "too many positional arguments for traxy"
        ### Try running a second time if it dies once
        try:
            status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-%s_drz.fits[1]' %(filter), direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
            os.remove("/tmp/%s.flt_xy" %(self.root))
        except:
            threedhst.process_grism.flprMulti()
            status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-%s_drz.fits[1]' %(filter), direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
            os.remove("/tmp/%s.flt_xy" %(self.root))
            
        self.ra_wcs, self.dec_wcs = np.zeros(NOBJ, dtype=np.double), np.zeros(NOBJ, dtype=np.double)
                    
        for i, line in enumerate(status[-NOBJ:]):
            if verbose:
                print unicorn.noNewLine+'Running iraf.xy2rd (%4d/%4d)' %(i+1, NOBJ)
            spl = np.cast[float](line.split())
            #
            if wcs is not None:
                self.ra_wcs[i], self.dec_wcs[i] = wcs.xy2rd((spl[2],spl[3]))
            else:
                status_xy2rd = iraf.xy2rd(infile=self.root+'-F140W_drz.fits[1]', x=spl[2], y=spl[3], hms=iraf.no, Stdout=1)
                self.ra_wcs[i], self.dec_wcs[i] = np.double(iraf.xy2rd.ra), np.double(iraf.xy2rd.dec)
            
        #### Write the wcsfix file
        fp = open('%s_inter.cat.wcsfix' %(self.baseroot),'w')
        fp.write('# id    mag   ra    dec\n')
        for i in range(self.cat.id.shape[0]):
            fp.write('%5d  %.3f  %.6f  %.6f\n' %(self.cat.id[i], self.cat.mag[i], self.ra_wcs[i], self.dec_wcs[i]))
            
        fp.close()
        
        if verbose:
            print 'Corrected coordinates: %s_inter.cat.wcsfix' %(self.baseroot)
            
        return True
        
    def make_wcs_region_file(self, filename=None):
        """
        Make a region file with the true WCS coordinates.  
        
        Default output `filename` is self.root+'_inter_wcs.reg'
        """
        if filename is None:
            filename = self.root+'_inter_wcs.reg'
        
        if self.ra_wcs is None:
            self.get_corrected_wcs()
        
        fp = open(filename,'w')
        fp.write('fk5\n')
        for i in range(len(self.cat.id)):
            fp.write('circle(%.6f, %.6f, 0.2") # text={%d}\n' %(self.ra_wcs[i], self.dec_wcs[i], self.cat.id[i]))
        fp.close()
        
        print 'WCS-correct region file: %s' %(filename)
        
    def init_object_spectra(self):
        
        self.lam_spec = []
        
        # if self.grism_element == 'G141':
        #     self.lam_spec = np.arange(0.95e4,1.8e4,22.)
        # 
        # if self.grism_element == 'G102':
        #     self.lam_spec = np.arange(0.73e4, 1.18e4, 10.)
        # 
        # #
        # if self.grism_element == 'G800L':
        #     self.lam_spec = np.arange(0.2e4, 1.2e4, 20.)
        wlim = grism_wlimit[self.grism_element]
        self.lam_spec = np.arange(wlim[0]*0.85, wlim[1]*1.15, wlim[2])
        
        self.obj_in_model = {}
        self.flux_specs = {}
        for obj in self.objects:
            #self.lam_specs[obj] = None
            self.flux_specs[obj] = None
            self.obj_in_model[obj] = False
            
    def make_total_flux(self):
        
        #### Total flux
        try:
            flux = self.flux
        except:
            self.make_flux()
            flux = self.flux
        #
        self.total_fluxes = {}
        print 'Total flux'
       
        # ## Test, compare c version
        # t0 = time.time()
        # for obj in self.cat['NUMBER']:
        #     print unicorn.noNewLine+'Object #%d' %(obj)
        #     self.total_fluxes[obj] = np.sum(flux[self.segm[0].data == obj])
        # 
        # told = time.time()
        #self.new_fluxes = {}

        fluxes = utils_c.total_flux(self.flux, self.segm[0].data)
        for obj in self.objects:
            #self.new_fluxes[obj] = fluxes[obj]
            self.total_fluxes[obj] = fluxes[obj]
                    
        # tnew = time.time()
        # print 'Compare total fluxes:  %.3f  %.3f' %(told-t0, tnew-told)
        # for obj in self.cat['NUMBER']:
        #     print ' %.7f %.7f' %(self.total_fluxes[obj], self.total_fluxes[obj]-self.new_fluxes[obj])
    
    def object_flux_thumb(self, id=0):
        ix = self.cat.id == id
        padx = self.ngrowx*self.growx + self.pad/2
        pady = self.ngrowy*self.growy + self.pad/2
        x_flt = (self.cat.x_pix[ix][0]-padx)/self.growx
        y_flt = (self.cat.y_pix[ix][0]-pady)/self.growy
        
        xc = int(np.round(self.cat.x_pix[ix][0]))
        yc = int(np.round(self.cat.y_pix[ix][0]))
        
        NX = np.cast[int](np.max([5*self.cat['FLUX_RADIUS'][ix][0], 10]))
        
        #NX = 17*np.maximum(self.growx, self.growy)
        thumb = self.flux[yc-NX:yc+NX, xc-NX:xc+NX]
        thumb_seg = self.segm[0].data[yc-NX:yc+NX, xc-NX:xc+NX] == id
        thumb *= thumb_seg
        return x_flt, y_flt, xc, yc, NX, thumb
        
    def fft_compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D'], verbose=False, normalize=True, place=True):
        import stsci.convolve as sc
        import mywfc3.grism
        conf = mywfc3.grism.aXeConf(os.getenv('THREEDHST') + '/CONF/G141.test30.conf')
        x_flt, y_flt, xc, yc, NX, thumb = self.object_flux_thumb(id=id)
        
        beam='A'
        
        beam_limit = conf.conf['BEAM%s' %(beam)]
        dx = np.arange(beam_limit[0]-1, beam_limit[1]+1, 1./self.growx)
        #dx = np.arange(beam_limit[0]-1, beam_limit[1]+1)
        dy, lam = conf.get_beam_trace(x=x_flt, y=y_flt, dx=dx, beam=beam)
        dy *= self.growy
        dx *= self.growx
        #
        delta_wave = lam[1]-lam[0]
        if delta_wave > 0:
            sens = unicorn.utils_c.interp_conserve_c(lam, unicorn.reduce.sens_files[beam]['WAVELENGTH'], unicorn.reduce.sens_files[beam]['SENSITIVITY']*1.e-17)*np.median(np.abs(np.diff(lam)))
        else:
            sens = unicorn.utils_c.interp_conserve_c(lam[::-1], unicorn.reduce.sens_files[beam]['WAVELENGTH'], unicorn.reduce.sens_files[beam]['SENSITIVITY']*1.e-17)[::-1]*np.median(np.abs(np.diff(lam)))
        
        if flux_spec is not None:
            sens *= unicorn.utils_c.interp_conserve_c(lam, lam_spec, flux_spec)
        
        trace = np.zeros((thumb.shape[0], dx.shape[0]))
        dyint = np.cast[int](dy)
        f = dy-dyint
        for i in range(dx.shape[0]):
            trace[NX+dyint[i]-1,i] = (1-f[i])*sens[i]
            trace[NX+dyint[i],i] = f[i]*sens[i]
        #
        obj = sc.convolve2d(trace, thumb, fft=1)
        if place:
            self.object *= 0
            self.object[yc-NX:yc+NX, xc+int(dx[0])-1:xc+int(dx[-1])] = obj
        else:
            return obj
            
    def fast_compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D'], verbose=False, normalize=True):
        import mywfc3.grism
        conf = mywfc3.grism.aXeConf(os.getenv('THREEDHST') + '/CONF/G141.test30.conf')

        x_flt, y_flt, xc, yc, NX, thumb = self.object_flux_thumb(id=id)
        
        self.object *= 0
        
        for beam in BEAMS:
            beam_limit = conf.conf['BEAM%s' %(beam)]
            dx = np.arange(beam_limit[0]-1, beam_limit[1]+1, 1./self.growx)
            #dx = np.arange(beam_limit[0]-1, beam_limit[1]+1)
            dy, lam = conf.get_beam_trace(x=x_flt, y=y_flt, dx=dx, beam=beam)
            dy *= self.growy
            dx *= self.growx
            #
            delta_wave = lam[1]-lam[0]
            if delta_wave > 0:
                sens = unicorn.utils_c.interp_conserve_c(lam, unicorn.reduce.sens_files[beam]['WAVELENGTH'], unicorn.reduce.sens_files[beam]['SENSITIVITY']*1.e-17)*np.median(np.abs(np.diff(lam)))
            else:
                sens = unicorn.utils_c.interp_conserve_c(lam[::-1], unicorn.reduce.sens_files[beam]['WAVELENGTH'], unicorn.reduce.sens_files[beam]['SENSITIVITY']*1.e-17)[::-1]*np.median(np.abs(np.diff(lam)))
            
            nonzero = sens > 0
            dy = dy[nonzero]
            dx = dx[nonzero]
            lam = lam[nonzero]
            sens = sens[nonzero]

            if flux_spec is not None:
                sens *= unicorn.utils_c.interp_conserve_c(lam, lam_spec, flux_spec)
                
            for i, di in enumerate(dx):
                #print i, di
                thumb_shift = nd.shift(thumb, (dy[i], 0))
                self.object[yc-NX:yc+NX, xc+di-NX:xc+di+NX] += thumb_shift*sens[i]
            
    def compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D'], verbose=False, normalize=True):
        import unicorn.reduce
        
        t0 = time.time()
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == id)[0][0]
        
        self.id = id
        
        xc = np.float(self.cat.x_pix[ii])
        yc = np.float(self.cat.y_pix[ii])
        
        #### Normalize the input spectrum to the direct image passband
        t0 = time.time()
        if (lam_spec is not None) & (normalize is True):            
            x_filt, y_filt = np.loadtxt(os.getenv('iref')+'/'+self.filter + '.dat', unpack=True)
            #y_filt_int = np.interp(lam_spec, x_filt, y_filt)
            y_filt_int = utils_c.interp_c(lam_spec, x_filt, y_filt)
            filt_norm = np.trapz(y_filt_int*flux_spec, lam_spec) / np.trapz(y_filt_int, lam_spec)
            if verbose:
                print 'Norm: %f' %(filt_norm)
            flux_spec /= filt_norm
            if verbose:
                t1 = time.time(); dt = t1-t0; t0=t1
                print 'Input spectrum is provided, normalize it for filter, %s (%.2f s)' %(self.filter, dt)

        #### Define the grism dispersion for the central pixel of an object.  
        #### Assume doesn't change across the object extent            
        t0 = time.time()
        orders, self.xi = unicorn.reduce.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=self.grow_factor, growx=self.growx, growy=self.growy, pad=self.pad, ngrowx=self.ngrowx, ngrowy =self.ngrowy, grism=self.grism_element, finex=self.fine_offset[0], finey=self.fine_offset[1]
        )
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Define grism dispersion (%.2f s).' %(dt)
            
        yord, xord = np.indices(orders.shape)
        beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), self.xi[5].reshape((1,-1)))
        
        non_zero = orders > 0
        #xord, yord, ford, word, sord = xord[non_zero], yord[non_zero], orders[non_zero], self.xi[2][non_zero], self.xi[3][non_zero]
        cast = np.int
        xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), self.xi[2][non_zero], self.xi[3][non_zero], beams[non_zero]
        
        #print non_zero.sum()
        
        ys = orders.shape
        xord += self.xi[0]
        yord -= (ys[0]-1)/2
        #
        self.object *= 0.
    
        mask = self.segm[0].data == id
        xpix = self.xf[mask]
        ypix = self.yf[mask]
        
        #### Loop through pixels defined within the segmentation region, summing
        #### dispersed flux
        # test_object = self.object*0.
        # for jj in range(xpix.size):
        #     x, y = xpix[jj], ypix[jj]
        #     xxi = x+xord
        #     yyi = y+yord
        #     use = (xxi >= 0) & (xxi < self.sh[1]) & (yyi >= 0) & (yyi < self.sh[0])
        #     test_object[yyi[use], xxi[use]] += ford[use]*self.flux[y,x] #*10
        # 
        # self.test_object = self.object*1.
        # self.object *= 0.
        # #self.object = object
        # if verbose:
        #     tx = time.time()
        #     print 'Loop through object pixels within the segmentation region. (%.3f)' %(tx-t0)
        #     t0 = tx
            
        #### Test cversion
        #print self.flux.dtype, self.segm[0].data.dtype, xord.dtype, yord.dtype, ford.dtype
        
        # yy, xx = np.indices(self.flux.shape)
        # ok = self.segm[0].data == id
        # xlim, ylim = [xx[ok].min(), xx[ok].max()], [yy[ok].min(), yy[ok].max()]

        ##### Only loop through segmentation pixels within 3.5*FLUX_RADIUS of
        ##### the central object pixel
        xlim, ylim = None, None
        # if ('FLUX_RADIUS' in self.cat.column_names):
        #     xci, yci = int(np.round(xc)), int(np.round(yc))
        #     ri = int(np.round(self.cat['FLUX_RADIUS'][ii]))
        #     NR = 3.5
        #     pix_radius = np.maximum(ri*NR, 10)
        #     xlim = np.clip([xci-pix_radius, xci+pix_radius], 0, self.sh[1])
        #     ylim = np.clip([yci-pix_radius, yci+pix_radius], 0, self.sh[0])
                
        status = utils_c.disperse_grism_object(self.flux, id, self.segm[0].data, xord, yord, ford, self.object, xpix=xlim, ypix=ylim)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Loop (c). (%.3f)' %(dt)
        
        #### Make full images that show the wavelength, sensitivity
        xxi = np.int(np.round(xc))+xord
        use = (xxi >= 0) & (xxi < self.sh[1])
        ### First order only
        #use = use & (xord > 10*self.grow_factor) & (xord < 213*self.grow_factor) 
        use_order = use & (bord == 1)
        #use_order = use & (bord & 1)
        cut = np.zeros(self.sh[1])
        cut[xxi[use_order]] = word[use_order] 
        self.object_wave = cut.copy() #np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
        self.full_wave = word[use]
        
        cut[xxi[use_order]] = sord[use_order] 
        self.object_sens = cut.copy() #np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Output wavelength and sensitivity arrays. (%.3f)' %(dt)
    
    def refine_model(self, id, BEAMS=['A','B','C','D'], view=False, verbose=False, MAG_FOR_SIMPLE=22):
        """
        MAG_FOR_SIMPLE:  use linear fit if mag > MAG_FOR_SIMPLE
        """
        from scipy import polyfit,polyval

        if id not in self.cat.id:
            print '#%d not in the catalog.' %(id)
            return False
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == id)[0][0]
                  
        #### First iteration
        #self.compute_object_model(id=id, BEAMS=['A'], lam_spec=self.lam_specs[id], flux_spec = self.flux_specs[id], verbose=verbose)
        self.compute_object_model(id=id, BEAMS=BEAMS, lam_spec=self.lam_spec, flux_spec = self.flux_specs[id], verbose=verbose, normalize=False)
        
        t0 = time.time()
        if view:
            first_model = self.object*1
        
        if self.obj_in_model[id]:
            self.model -= self.object
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Refine, first_model  (%.3f)' %(dt)
        
        #### Need to recompute with only the 
        if self.grism_element == 'G800L':
            self.compute_object_model(id=id, BEAMS=['A'], lam_spec=self.lam_spec, flux_spec = self.flux_specs[id], verbose=verbose, normalize=False)
            
        # whts = self.object**4
        # finite_mask = np.isfinite(self.model) & (self.gris[2].data != 0)
        # whts[~finite_mask] = 0
        # whts[self.gris[1].data == 0] = 0
        # 
        # ratio = (self.gris[1].data - (self.model)) / self.object
        # ratio[~finite_mask] = 0.
        # ratio[~np.isfinite(ratio)] = 0
        # whts[~np.isfinite(ratio)] = 0
        # 
        # ratio_extract = np.sum(ratio*whts, axis=0) / np.sum(whts, axis=0)
        # ratio_extract[np.sum(ratio*whts, axis=0) == 0] = 0.
        
        #ratio_extract = utils_c.get_model_ratio(self.object, self.model, self.gris[1].data)
        ### Try using only pixels above some fraction of the peak
        # clip = self.object*1
        # clip[clip/clip.max(axis=0) < 0.1] = 0
        # ratio_extract = utils_c.get_model_ratio_optimal(clip, self.model, self.gris[1].data, self.gris[2].data)
        ratio_extract = utils_c.get_model_ratio_optimal(self.object, self.model, self.gris[1].data, self.gris[2].data)
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Make ratio  (%.3f)' %(dt)
        
        #wave = np.mean(self.object_wave, axis=0)
        wave = self.object_wave
        xpix = np.arange(self.sh[1])
        
        if self.grism_element == 'G141':
            keep = (wave > 1.12e4) & (wave < 1.65e4) & (ratio_extract != 0)  # (xpix > (self.pad-22)) & (xpix < (self.sh[1]-self.pad-22))
        elif self.grism_element == 'G102':
            keep = (wave > 0.78e4) & (wave < 1.14e4) & (ratio_extract != 0)
        elif self.grism_element == 'GRS':
            keep = (wave > 1.35e4) & (wave < 1.95e4) & (ratio_extract != 0)
        else:
            keep = (wave > 0.55e4) & (wave < 1.14e4) & (ratio_extract != 0)
           
        #print len(wave), len(wave[wave > 0]), wave.max(), wave[2064], len(ratio_extract)
        
        if keep.sum() < 10:
            self.model += self.object
            self.obj_in_model[id] = True
            if verbose:
                print 'Skipping refine: only %d wavelength steps.' %(keep.sum())
            
            return True
            
        coeff = polyfit(wave[keep], ratio_extract[keep], 1)
        ratio_extract[~keep] = polyval(coeff, wave[~keep])
        #ratio_extract = threedhst.utils.medfilt(ratio_extract,21, AVERAGE=True)
        if self.cat.mag[ii] > MAG_FOR_SIMPLE:
            ratio_extract[keep] = polyval(coeff, wave[keep])
            
        #ratio_extract = np.maximum(ratio_extract, 0.33)
        
        ratio_extract = np.maximum(ratio_extract, 0.1)
        #print np.min(ratio_extract)
        
        #rint = utils_c.interp_c(self.lam_spec, wave[wave > 0], ratio_extract[wave > 0])
        rint = np.interp(self.lam_spec, wave[wave > 0], ratio_extract[wave > 0], left=0.0, right=0.0)
        
        #print np.min(rint)
        if self.flux_specs[id] is not None:
            self.flux_specs[id] *= rint
            #print np.min(self.flux_specs[id]), np.min(rint), self.flux_specs[id].shape, rint.shape
        else:
            #self.lam_specs[id] = wave[wave > 0]
            self.flux_specs[id] = rint
        
        #### Enforce minimum for "flux_specs" scaling to avoid divide-by-large-numbers
        #self.flux_specs[id] = np.maximum(self.flux_specs[id], 0.2)
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Extract ratio, save values  (%.3f)' %(dt)
            
        self.compute_object_model(id=id, BEAMS=BEAMS, lam_spec = self.lam_spec, flux_spec = self.flux_specs[id], normalize=False, verbose=verbose)
        if view:
            second_model = self.object.copy()
        
        #### View the results
        if view:
            print 'Display results for #%d to DS9' %(id)
            diff_first = self.gris[1].data - self.model - first_model
            diff_second = self.gris[1].data - self.model - second_model

            view.frame(1)
            view.v(self.gris[1].data, vmin=-0.1, vmax=5)
            view.frame(2)
            view.v(diff_first, vmin=-0.1, vmax=5)
            view.frame(3)
            view.v(diff_second, vmin=-0.1, vmax=5)

            idx = (self.cat['NUMBER'] == id)
            view.set('pan to %f %f' %(self.cat.x_pix[idx][0], self.cat.y_pix[idx][0]))
            view.set('zoom to 1.4')
            #view.set('cmap value 4.43936 0.462911')
            view.set('cmap HSV')
            view.set('cmap value 7.02222 0.0867647')
            view.set('regions file %s_inter.reg' %(self.root))
            view.match_all(alignType='image')
            
            plt.plot(self.lam_spec, self.flux_specs[id])
        
        ### Done, now update the total model
        if view:
            self.model += second_model
            self.obj_in_model[id] = True
        else:
            self.model += self.object
            self.obj_in_model[id] = True

    def show_ratio_spectrum(self, id, flux=False):
        if id not in self.cat['NUMBER']:
            print '#%d not in the catalog.' %(id)
            return False
        
        if self.flux_specs[id] is None:
            show = self.lam_spec*0.+1
        else:
            show = self.flux_specs[id].copy()
        
        if flux:
            show *= self.total_fluxes[id]
            
        plt.plot(self.lam_spec, show)
    
    def remove_object(self, id, BEAMS=['A','B','C','D'], verbose=False):
        
        if id not in self.cat.id:
            print '#%d not in the catalog.' %(id)
            return False
        
        self.compute_object_model(id=id, BEAMS=BEAMS, lam_spec=self.lam_spec, flux_spec = self.flux_specs[id], verbose=verbose, normalize=False)
                
        if self.obj_in_model[id]:
            self.model -= self.object
            self.obj_in_model[id] = False
            
    def compute_full_model(self, BEAMS=['A','B','C','D'], view=None, MAG_LIMIT=23., save_pickle=True, refine=True, model_slope=0, normalize_model_list=True):
        
        #self.model*=0.
        if refine:
            threedhst.showMessage('Making full object model (refined model)')
        else:
            threedhst.showMessage('Making full object model (flat spectra)')
            
        mag = self.cat['MAG_AUTO']
        so = np.argsort(mag)
        so = so[mag[so] <= MAG_LIMIT]
        
        if self.use_segm:
            so = np.arange(len(self.total_fluxes.keys()))
            
        N = len(so)
        for i in range(N):
            #for i, id in enumerate(self.cat['NUMBER'][so]):
            id = self.objects[so][i]
            # if np.abs(id-854) > 10:
            #     continue
            # 
            print unicorn.noNewLine+'Object #%d, m=%.2f (%d/%d)' %(id, mag[so][i], i+1, N)
            #print 'Object #%d, m=%.2f (%d/%d)' %(id, mag[so][i], i+1, N)
            if refine:
                self.refine_model(id, BEAMS=BEAMS, view=view)      
            else:
                # if self.grism_element == 'G141':
                #     lx = np.arange(0.9e4,1.8e4)
                #     #ly = 1+model_slope*(lx-1.4e4)/4.e3
                #     ly = (lx/1.4e4)**model_slope
                # 
                # if self.grism_element == 'G102':
                #     lx = np.arange(0.6e4,1.3e4)
                #     #ly = 1+model_slope*(lx-0.98e4)/2.8e3
                #     ly = (lx/0.98e4)**model_slope
                # 
                # #
                # if self.grism_element == 'G800L':
                #     lx = np.arange(0.3e4,1.2e4)
                #     #ly = 1+model_slope*(lx-0.75e4)/2.5e3
                #     ly = (lx/0.75e4)**model_slope                
                wlim = grism_wlimit[self.grism_element]
                lx = np.arange(wlim[0]*0.8,wlim[1]*1.2)
                ly = (lx/wlim[3])**model_slope
                
                #print lx.min(), lx.max()
                
                if self.model_list is not None:
                    if id in self.model_list.keys():
                        print unicorn.noNewLine+'Object #%d, m=%.2f (%d/%d) [model_list]' %(id, mag[so][i], i+1, N)
                        lx, ly = self.model_list[id][0]*1., self.model_list[id][1]*1.
                        self.flux_specs[id] = utils_c.interp_conserve_c(self.lam_spec, lx, ly, left=0., right=0.)
                else:
                    normalize_model_list=True
                    
                self.compute_object_model(id, BEAMS=BEAMS, lam_spec=lx, flux_spec=ly, normalize=normalize_model_list, verbose=False)      
                self.model += self.object
                self.obj_in_model[id] = True
                
                #### ly modified in previous step
                if self.model_list is not None:
                    if id in self.model_list.keys():
                        self.flux_specs[id] = utils_c.interp_conserve_c(self.lam_spec, lx, ly, left=0., right=0.)  
                                      
                if view:
                    try:
                        view.view(self.model)                    
                    except:
                        pass
                    #
                    # xx = raw_input('Continue?')
                    # if xx != '':
                    #     print this_breaks
                    
        if save_pickle:
            self.save_model_spectra(BEAMS=BEAMS)
            
    def save_model_spectra(self, BEAMS=['A','B','C','D']):
        """
        Save model "pickle" and the 2D model image.
        
        If just BEAM ["B"] is specified, save the 0th order contamination image.
        """
        import pickle
        import unicorn.reduce
        
        ####
        if BEAMS == ['B']:
            pyfits.writeto(self.root+'_inter_0th.fits', data=self.model, header=self.gris[1].header, clobber=True)
            return
            
        fp = open(self.baseroot+'_inter_model.pkl','wb')
        pickle.dump(self.cat.id, fp)
        pickle.dump(self.obj_in_model, fp)
        pickle.dump(self.lam_spec, fp)
        pickle.dump(self.flux_specs, fp)
        fp.close()
        
        header = self.gris[1].header.copy()
        #header['GRISCONF'] = (unicorn.reduce.conf_file, 'Grism configuration file')
        header.update('GRISCONF', unicorn.reduce.conf_file, comment='Grism configuration file')
        if 'VERSION' in unicorn.reduce.conf.keys():
            header.update('GRISCREV', unicorn.reduce.conf['VERSION'], comment='Grism conf. file revision')
        else:
            header.update('GRISCREV', "---", comment='Grism conf. file revision')
            
        pyfits.writeto(self.baseroot+'_inter_model.fits', data=self.model, header=header, clobber=True)
        
    def load_model_spectra(self, use_same_config=True):
        """
        Load model.pkl and model.fits files if they exist
        """
        import pickle
        import unicorn.reduce
        
        if ((not os.path.exists(self.baseroot+'_inter_model.pkl')) | 
            (not os.path.exists(self.baseroot+'_inter_model.fits'))):
            return False
            
        fp = open(self.baseroot+'_inter_model.pkl','rb')
        ids = pickle.load(fp)
        test = ids == self.cat.id
        if np.sum(test) == len(self.cat.id):
            print 'Load spectra from pickle %s_inter_model.pkl' %(self.baseroot)
            self.obj_in_model = pickle.load(fp)
            self.lam_spec = pickle.load(fp)
            self.flux_specs = pickle.load(fp)
            
            im = pyfits.open(self.baseroot+'_inter_model.fits')
            self.model = np.cast[np.double](im[0].data)
            
            #### Make sure we're using the same configuration file used 
            #### to generate the model
            if 'GRISCONF' in im[0].header.keys():
                if unicorn.reduce.conf_file != im[0].header['GRISCONF']:
                    unicorn.reduce.set_grism_config(grism=self.grism_element, force=True, use_config_file=im[0].header['GRISCONF'])
                
            im.close()
        else:
            print "ERROR: Object list in pickle doesn't match the current catalog."
            fp.close()
            return False
            
        fp.close()
        return True
        
    def twod_spectrum(self, id=328, grow=1, miny=18, maxy=None, CONTAMINATING_MAGLIMIT=23, refine=True, verbose=False, force_refine_nearby=False, USE_REFERENCE_THUMB=True, USE_FLUX_RADIUS_SCALE=3, BIG_THUMB=False, extract_1d=True, check_seg=True, wlim=None):
        """
        Extract a 2D spectrum from the interlaced image.
        
        `miny` is the minimum Y-size of the 2D spectrum.  If `miny` < 0: then 
         use that size explicitly rather than just setting a minimum.
         
        `CONTAMINATING_MAGLIMIT` is the limit for computing contaminating
        objects.
        
        `USE_REFERENCE_THUMB`: Use the reference image rather than the 
        nominal F140W image for the 2D thumbnail.
        """
        import unicorn.reduce
        
        if id not in self.cat.id:
            print 'Object #%d not in catalog.' %(id)
            return False
        
        unicorn.reduce.set_grism_config(grism=self.grism_element)            
        #
        ii = np.where(np.cast[int](self.cat.id) == id)[0][0]
        xc = np.int(np.round(self.cat.x_pix[ii]))-1
        yc = np.int(np.round(self.cat.y_pix[ii]))-1
        
        #count=0; print 'HERE%d' %(count)
        
        seg = np.cast[int](self.segm[0].data)
        
        #### Check that center within the image
        if (xc < 0) | (xc > (self.sh[1]-1)) | (yc < 0) | (yc > (self.sh[0]-1)):
            threedhst.showMessage('Object %d at (%d, %d) is off the image (%d,%d)' %(id, xc, yc, self.sh[1], self.sh[0]), warn=True)
            return False
        
        if (seg[yc,xc] != id) & (check_seg):
            threedhst.showMessage('Segmentation image at (%d, %d) not equal to object %d.' %(xc, yc, id), warn=True)
            return False
        
        # seg_mask = seg == id
        # ypix, xpix = self.yf, self.xf
        # 
        # NY = ypix[seg_mask].max()-ypix[seg_mask].min()
        # NX = xpix[seg_mask].max()-xpix[seg_mask].min()        
        if miny < 0:
            NT = -miny
        else:
            if (USE_FLUX_RADIUS_SCALE is not None) & ('FLUX_RADIUS' in self.cat.column_names):
                NT = np.max([USE_FLUX_RADIUS_SCALE*2*self.cat['FLUX_RADIUS'][ii], miny])
            else:
                #### Use only contiguous segmentation region around center pixel
                xmi, xma, ymi, yma = threedhst.utils.contiguous_extent(seg == id, xc, yc)
                NX = np.maximum(xc-xmi, xma-xc)*2
                NY = np.maximum(yc-ymi, yma-yc)*2
                NT = np.max([NY*grow, NX*grow, miny])
        
        #NT = np.min([self.pad/2, NT])
                                
        NT = np.min([NT, xc, yc, self.sh[1]-xc, self.sh[0]-yc])
        
        if (miny < 0) & (NT < np.abs(miny)):
            print "not enough y pixels"
            return False
                 
        #### Maximum size
        if maxy is not None:
            NT = np.min([NT, maxy])
        
        #print '!!! NT: %d' %(NT)
        
        #### Check that desired thumbnail actually within the image
        #print 'XXX', xc, NT, self.sh
        if ((xc-NT/2) < (0)) | ((xc+NT/2) > (self.sh[1])):
            threedhst.showMessage('Object %d (%d,%d) outside of image area (%d, %d) for desired thumbnail size (%d).' %(id, xc, yc, self.sh[1], self.sh[0], NT), warn=True)
            return False
        
        if verbose:
            print 'Generating direct image extensions'
            
        self.direct_thumb = self.direct[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_inter = self.direct_thumb*1.
        if len(self.direct) > 2:
            self.direct_wht = self.direct[2].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        else:
            self.direct_wht = self.direct_thumb*0+1.
            
        self.direct_seg = np.cast[np.int32](seg[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2])
                    
        #### Initialize the output FITS file
        prim = pyfits.PrimaryHDU()
        prim.header.update('POINTING', self.baseroot)
        prim.header.update('ID', id)
        prim.header.update('X_PIX', self.cat.x_pix[ii], comment='X pixel in interlaced image')
        prim.header.update('Y_PIX', self.cat.y_pix[ii], comment='Y pixel in interlaced image')
        prim.header.update('RA', self.ra_wcs[ii], comment='Target R.A.')
        prim.header.update('DEC', self.dec_wcs[ii], comment='Target Dec.')
        prim.header.update('MAG', self.cat.mag[ii], comment='MAG_AUTO from interlaced catalog')
        prim.header.update('FTOTAL', self.total_fluxes[id], comment='Total flux within segmentation image (1E-17)')
        prim.header.update('PAD', self.pad, comment='Padding at edge of interlaced image')
        prim.header.update('NGROWX', self.ngrowx, comment='Additional extra pixels at the image edges in X')
        prim.header.update('NGROWY', self.ngrowy, comment='Additional extra pixels at the image edges in Y')
        prim.header.update('GROW', self.grow_factor, comment='"grow_factor" from unicorn.reduce')
        prim.header.update('GROWX', self.growx, comment='"growx" from unicorn.reduce')
        prim.header.update('GROWY', self.growy, comment='"growy" from unicorn.reduce')
        prim.header.update('BUNIT', 'ELECTRONS/S')
        prim.header.update('EXPTIME', self.gris[0].header['EXPTIME'])
        prim.header.update('FILTER', self.filter)
        prim.header.update('GRISM', self.gris[0].header['FILTER'])
        prim.header.update('CCDCHIP', self.chip)
        
        prim.header.update('XDIRECT0', xc-NT/2, comment="Lower left pixel x in direct image")
        prim.header.update('YDIRECT0', yc-NT/2, comment="Lower left pixel y in direct image")
        
        prim.header.update('REFTHUMB', False, comment='Thumbnail comes from reference image')
        prim.header.update('GRISCONF', unicorn.reduce.conf_file, comment='Grism configuration file')
        if 'VERSION' in unicorn.reduce.conf.keys():
            prim.header.update('GRISCREV', unicorn.reduce.conf['VERSION'], comment='Grism conf. file revision')
        else:
            prim.header.update('GRISCREV', "---", comment='Grism conf. file revision')
        
        if '_ref_inter' in self.im.filename():
            ### Use reference image as thumbnail if it exists
            if (xc < ((self.pad + self.ngrowx*self.growx)/ 2)) | USE_REFERENCE_THUMB:
                prim.header.update('REFTHUMB', True)
                self.direct_thumb = self.im[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
                prim.header.update('FILTER', self.im[0].header['FILTER'])
        
        prim.header['FINEX'] = self.fine_offset[0]
        prim.header['FINEY'] = self.fine_offset[1]
                
        extensions = [prim]

        header = pyfits.ImageHDU().header
                
        #### Direct thumbnails
        header.update('EXTNAME','DSCI')
        extensions.append(pyfits.ImageHDU(self.direct_thumb, header))
        header.update('EXTNAME','DINTER')
        extensions.append(pyfits.ImageHDU(self.direct_inter, header))
        header.update('EXTNAME','DWHT')
        extensions.append(pyfits.ImageHDU(self.direct_wht, header))
        header.update('EXTNAME','DSEG')
        extensions.append(pyfits.ImageHDU(self.direct_seg, header))
                
        #### Find objects that could contaminate the current one
        radius = None
        if 'KRON_RADIUS' in self.cat.column_names:
            radius = self.cat['KRON_RADIUS']
        #
        if 'FLUX_RADIUS' in self.cat.column_names:
            radius = self.cat['FLUX_RADIUS']
        
        if radius is None:
            radius = 30.  ### stop-gap
                    
        dylo = self.cat.y_pix-self.cat.y_pix[ii]-radius*3
        dy0 = self.cat.y_pix-self.cat.y_pix[ii]
        dyhi = self.cat.y_pix-self.cat.y_pix[ii]+radius*3
        dy = np.minimum(np.abs(dylo), np.abs(dy0))
        dy = np.minimum(np.abs(dyhi), dy)
                
        dx = self.cat.x_pix-self.cat.x_pix[ii]
        nearby = (dy < (NT/2.+5)) & (dx < 420*self.growx) & (self.cat.mag < CONTAMINATING_MAGLIMIT)
        
        BEAMS=['A','B','C','D'] #,'E']
        view=False
        if nearby.sum() > 0:
            ids = self.cat.id[nearby][np.argsort(self.cat.mag[nearby])]
            for idi in ids:
                if verbose:
                    print 'Calculating contaminating object #%d' %(idi)
                #
                if refine:
                    if (not self.obj_in_model[idi]) | force_refine_nearby:
                        if verbose:
                            print unicorn.noNewLine + 'Calculating contaminating object #%d - ' %(idi)
                        self.refine_model(idi, BEAMS=BEAMS, view=view)      
                else:
                    if (not self.obj_in_model[idi]) | force_refine_nearby:
                        if verbose:
                            print unicorn.noNewLine + 'Calculating contaminating object #%d - ' %(idi)

                        self.compute_object_model(idi, BEAMS=BEAMS, lam_spec=None, flux_spec=None)      
                        self.model += self.object
                        self.obj_in_model[idi] = True
                        
        #### Generate model for current object
        if verbose:
            print 'Compute object model #%d' %(id)

        if refine:
            self.refine_model(id, BEAMS=BEAMS, view=view)      
        else:
            self.compute_object_model(id, BEAMS=BEAMS, lam_spec=self.lam_spec, flux_spec=self.flux_specs[id], normalize=False)      
            if not self.obj_in_model[id]:
                self.model += self.object
                self.obj_in_model[id] = True
        #
        self.compute_object_model(id, BEAMS=['A'], lam_spec=self.lam_spec, flux_spec=self.flux_specs[id], normalize=False)      
        
        #### Find pixels of the 1st order        
        xarr = np.arange(self.sh[1])
        # if self.grism_element == 'G141':
        #     wavelength_region = (self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)
        # #
        # if self.grism_element == 'G102':
        #     wavelength_region = (self.object_wave >= 0.73e4) & (self.object_wave <= 1.18e4)
        # #
        # if self.grism_element == 'G800L':
        #     wavelength_region = (self.object_wave >= 0.58e4) & (self.object_wave <= 0.92e4)
        #
        if wlim is None:
            wlim = grism_wlimit[self.grism_element]
        
        extensions[0].header['WMIN'] = wlim[0]
        extensions[0].header['WMAX'] = wlim[1]
        extensions[0].header['DW'] = wlim[2]
        extensions[0].header['W0'] = wlim[3]
        
        wavelength_region = (self.object_wave >= wlim[0]) & (self.object_wave <= wlim[1])
        
        if np.sum(wavelength_region) < 20:
            fp = open(self.baseroot+'_%05d.2D.xxx' %(id),'w')
            fp.write('%.2f %.2f  %d %d\n' %(xc, yc, -1, -1))
            fp.close()
            self.twod_status = False
            return False
        
        self.twod_status = True
        
        xmin = xarr[wavelength_region].min()
        xmax = xarr[wavelength_region].max()
        
        #print 'TWOD_SPECTRUM xminmax: %d %d' %(xmin, xmax)
        
        #### Try hard-coding the 1st order extraction region to keep constant
        #first_order = {'G141':[28,178], 'G102':[45, 227], 'G800L':[60, 192]}
        first_order = {'G141':[26,182], 'G102':[45, 225], 'G800L':[60, 192]}
        xx = first_order[self.grism_element]
        xmin, xmax = first_order[self.grism_element]
        xmin *= self.growx
        xmax *= self.growx
        # extensions[0].header['DXMIN'] = (xmin, 'Left edge of first order spectrum')
        # extensions[0].header['DXMAX'] = (xmax, 'Right edge of first order spectrum')
        xmin += xc
        xmax += xc
        #print 'FORCE: %d %d' %(xmin, xmax)
        
        if ((self.sh[1]-xmin) < 20) | ((xmax-xmin) < 20):
            fp = open(self.baseroot+'_%05d.2D.xxx' %(id),'w')
            fp.write('%.2f %.2f  %d %d\n' %(xc, yc, xmin, xmax))
            fp.close()
            return False
            
        #### dydx offset for grism trace
        xoff_arr = np.arange(len(self.xi[4]), dtype=np.double)+self.xi[0]+xc
        xoff_int = np.arange(xmin, xmax, dtype=np.double)
        yoff_int = utils_c.interp_c(xoff_int, xoff_arr, self.xi[4])

        yc_spec = np.int(np.round(self.cat.y_pix[ii]+yoff_int.mean()))-1
        spec_y_offset = yc-yc_spec
        
        self.ytrace = NT/2.+yoff_int-yoff_int.mean()
                        
        self.wave = self.object_wave[xmin:xmax]
        self.dwave = self.wave[1]-self.wave[0]
        
        #print '\nxxx: wave %f %f\nxxx\n' %(self.wave.min(), self.wave.max())
        
        if verbose:
            print 'Writing spectrum extensions.'
            
        ### WCS header
        header.update('CRPIX1',1)
        header.update('CD1_1',self.dwave)
        header.update('CD1_2',0)
        header.update('CRVAL1',self.wave[0])
        header.update('CUNIT1','Angstrom')
        header.update('CTYPE1','WAVE')
        
        # header.update('CRPIX2',NT/2.+1)
        # header.update('CRVAL2',0)
        # header.update('CDELT2',0.128254/self.growy)
        header.update('CRPIX2', NT/2+1)
        scl = 0.128254/self.growy
        header.update('CRVAL2', ((NT/2+1)-self.ytrace[0])*scl)
        header.update('CD2_1', -np.diff(self.ytrace)[0]*scl)
        header.update('CD2_2', scl)
        
        header.update('CUNIT2','arcsec')
        header.update('CTYPE2','LINEAR')
        
        header.update('YTOFF', spec_y_offset, comment='Offset of middle of spectrum trace')
        
        self.grism_sci = self.gris[1].data[yc_spec-NT/2:yc_spec+NT/2, xmin:xmax]
        
        header.update('EXTNAME','SCI')
        extensions.append(pyfits.ImageHDU(self.grism_sci, header))
        
        self.grism_wht = self.gris[2].data[yc_spec-NT/2:yc_spec+NT/2, xmin:xmax]
        header.update('EXTNAME','WHT')
        header.update('WHTERROR',self.gris[0].header['WHTERROR'])
        
        extensions.append(pyfits.ImageHDU(self.grism_wht, header))
        
        self.the_object = self.object[yc_spec-NT/2:yc_spec+NT/2, xmin:xmax]
        header.update('EXTNAME','MODEL')
        extensions.append(pyfits.ImageHDU(self.the_object, header))
        
        self.full_model = self.model[yc_spec-NT/2:yc_spec+NT/2, xmin:xmax]
        header.update('EXTNAME','CONTAM')
        extensions.append(pyfits.ImageHDU(self.full_model-self.the_object, header))
                
        header.update('EXTNAME','WAVE')
        extensions.append(pyfits.ImageHDU(self.wave, header))
        
        self.sens = self.object_sens[xmin:xmax]
        header.update('EXTNAME','SENS')
        header.update('BUNIT','(ELECTRONS/S) / (FLAM*1.E17)')
        extensions.append(pyfits.ImageHDU(self.sens, header))
        
        header.update('EXTNAME','YTRACE')
        header.update('BUNIT','PIXELS')
        extensions.append(pyfits.ImageHDU(self.ytrace, header))
        
        extensions[0].header.update('XGRISM0', xmin, comment="Lower left pixel x in grism image")
        extensions[0].header.update('YGRISM0', yc_spec-NT/2, comment="Lower left pixel y in grism image")
        
        hdu = pyfits.HDUList(extensions)
        if BIG_THUMB:
            hdu.writeto(self.baseroot+'-big_%05d.2D.fits' %(id), clobber=True)
        else:
            hdu.writeto(self.baseroot+'_%05d.2D.fits' %(id), clobber=True)
        
        self.twod = hdu
        
        if extract_1d:
            self.oned_spectrum(verbose=verbose)
        
        return True
        
    def show_2d(self, savePNG = False, verbose=False):
        import unicorn
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == self.id)[0][0]
        
        if verbose:
            print 'Setup figure'
            
        if savePNG:
            fig = Figure(figsize=[5,5*1.2], dpi=100)
            fig.subplots_adjust(wspace=0.2, hspace=0.08,left=0.02, bottom=0.02, right=0.98, top=0.98)
        else:
            fig = unicorn.catalogs.plot_init(xs=5, left=0.05, right=0.05, bottom=0.05, top=0.05, aspect=1./1.2)
            
        xx = np.arange(self.grism_sci.shape[1])
        ltick = np.array([1.2,1.3,1.4,1.5,1.6])
        xtick = np.interp(ltick*1.e4, self.wave, xx)
        
        if verbose:
            print unicorn.noNewLine+'Thumb panel'
        
        ax = fig.add_subplot(431)
        ax.imshow(self.direct_thumb, vmin=-0.5, vmax=2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        ax.text(1.1,0.8, self.baseroot, transform=ax.transAxes)
        ax.text(1.1,0.6, '#%d' %(self.id), transform=ax.transAxes)
        ax.text(1.1,0.4, r'$m_{140}=%.2f$' %(self.cat.mag[ii]), transform=ax.transAxes)
        ax.text(1.1,0.2, r'$(%.1f, %.1f)$' %(self.cat.x_pix[ii], self.cat.y_pix[ii]), transform=ax.transAxes)
        
        if verbose:
            print unicorn.noNewLine+'Plot flux_specs'
        
        if self.flux_specs[self.id] is not None:
            ax = fig.add_subplot(422)
            ax.plot(self.lam_spec, self.flux_specs[self.id]*self.total_fluxes[self.id])
            xx = ax.set_xticks(ltick*1.e4)
            ax.set_xticklabels(ltick)

            wlim = grism_wlimit[self.grism_element]
            ax.set_xlim(wlim[0], wlim[1])
            # ax.set_xlim(1.1e4,1.65e4)
            # if self.grism_element == 'G102':
            #     ax.set_xlim(0.73e4,1.18e4)
            # #
            # if self.grism_element == 'G800L':
            #     ax.set_xlim(0.58e4,0.92e4)
            
            y = self.flux_specs[self.id]*self.total_fluxes[self.id]
            ax.set_ylim(-0.1*y.max(), 1.1*y.max())
            ax.set_yticklabels([]);
        
        if verbose:
            print unicorn.noNewLine+'GRISM - science'
        
        ax = fig.add_subplot(412)
        ax.imshow(self.grism_sci, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        if verbose:
            print unicorn.noNewLine+'GRISM - model'

        ax = fig.add_subplot(413)
        ax.imshow(self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        if verbose:
            print unicorn.noNewLine+'GRISM - difference'
        
        ax = fig.add_subplot(414)
        ax.imshow(self.grism_sci-self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)
        
        if verbose:
            print unicorn.noNewLine+'Save the PNG'
        
        if savePNG:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(self.baseroot+'_%05d.2D.png' %(self.id), dpi=100, transparent=False)
                        
    def oned_spectrum(self, verbose=True):
        """
        1D extraction after generating the 2D spectra.
        
        Currently just testing
        """
        
        t0 = time.time()
        
        obj_cleaned = self.grism_sci - (self.full_model - self.the_object)
        contam_cleaned = self.grism_sci - self.the_object
        variance = self.grism_wht**2
            
        #### Get extraction window where profile is > 10% of maximum 
        osh = self.the_object.shape
        xarr = np.arange(osh[1])
        # xint14 = int(np.interp(1.3e4, self.wave, xarr))
        # if self.grism_element == 'G102':
        #     xint14 = int(np.interp(0.98e4, self.wave, xarr))
        # #
        # if self.grism_element == 'G800L':
        #     xint14 = int(np.interp(0.75e4, self.wave, xarr))
        wlim = grism_wlimit[self.grism_element]
        xint14 = int(np.interp(wlim[3], self.wave, xarr))
        
        yprof = np.arange(osh[0])
        profile = np.sum(self.the_object[:,xint14-10:xint14+10], axis=1)
        profile = profile/profile.sum()
        
        #profile2D = np.dot(profile.reshape((-1,1)), np.ones((1,osh[1])))
        window = profile >= 0.1*profile.max()
        
        #### Use object as the profile since it contains the y shift along the trace
        prof_x = np.sum(self.the_object, axis=0)
        profile2D = self.the_object/np.dot(np.ones((osh[0],1)), prof_x.reshape((1,-1)))
        
        if (1 == 0):
            plt.plot(yprof, profile, color='blue')
            plt.plot(yprof[window], profile[window], color='red', linewidth=3)
            print np.trapz(profile[window], xprof[window]) / np.trapz(profile, xprof)
        
        obj_cleaned[variance == 0] = 0
        contam_cleaned[variance == 0] = 0
        
        #### Simple sums, also weighted
        weights = 1./variance
        weights[variance == 0] = 0
        
        var_sum = np.sum(weights[window,:], axis=0)
        
        weighted_sum = np.sum((obj_cleaned*weights)[window,:], axis=0) / var_sum
        weighted_var = 1./np.sum(weights[window,:], axis=0)
        weighted_sig = np.sqrt(weighted_var)
        
        simple_sum = np.sum(obj_cleaned[window,:], axis=0)
        simple_var = np.sum(variance[window,:], axis=0)
        simple_sig = np.sqrt(simple_var)
        
        #### Optimal extraction
        opt_variance = variance.copy()
        opt_variance[opt_variance == 0] = 1.e6
        
        num = np.sum(profile2D*obj_cleaned/opt_variance, axis=0)
        #num_contam = np.sum(profile2D*contam_cleaned/opt_variance, axis=0)
        num_contam = np.sum(profile2D*(self.full_model - self.the_object)/opt_variance, axis=0)
        num_full = np.sum(profile2D*self.grism_sci/opt_variance, axis=0)
        
        denom = np.sum(profile2D**2 / opt_variance, axis=0)
        
        optimal_sum = num / denom
        optimal_sum_contam = num_contam / denom
        optimal_sum_full = num_full / denom
        
        optimal_var = 1./np.sum(profile2D**2/opt_variance, axis=0)
        optimal_sig = np.sqrt(optimal_var)
        
        trace_pix = np.cast[int](np.round(self.ytrace))
        #trace_spec = self.grism_sci[trace_pix,:]
        trace_spec = optimal_sum*0.
        trace_sig = optimal_sum*0.
        for i in range(osh[1]):
            trace_spec[i] = self.grism_sci[trace_pix[i],i]
            trace_sig[i] = self.grism_wht[trace_pix[i],i]
        
        scale_to_total = 1./np.max(profile)       
        if not np.isfinite(scale_to_total):
            scale_to_total=-1
            
        c1 = pyfits.Column(name='wave', format='D', unit='ANGSTROMS', array=self.wave)
        c2 = pyfits.Column(name='flux', format='D', unit='ELECTRONS/S', array=optimal_sum_full)
        c3 = pyfits.Column(name='error', format='D', unit='ELECTRONS/S', array=optimal_sig)
        c4 = pyfits.Column(name='contam', format='D', unit='ELECTRONS/S', array=optimal_sum_contam)
        c5 = pyfits.Column(name='trace', format='D', unit='ELECTRONS/S', array=trace_spec)
        c6 = pyfits.Column(name='etrace', format='D', unit='ELECTRONS/S', array=trace_sig)
        c7 = pyfits.Column(name='sensitivity', format='D', unit='E/S / 1E-17 CGS', array=self.sens*(self.growx*self.growy))
        #print 'MAX SENS: %.f' %(self.sens.max())
        
        coldefs = pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7])
        head = pyfits.Header()
        head.update('ATRACE', scale_to_total, comment='Factor to scale trace to total flux')
        
        ii = np.where(np.cast[int](self.cat.id) == self.id)[0][0]
        try:
            head.update('RA', self.ra_wcs[ii], comment='Target R.A.')
            head.update('DEC', self.dec_wcs[ii], comment='Target Dec.')
        except:
            head.update('RA', 0., comment='Target R.A.')
            head.update('DEC', 0., comment='Target R.A.')
        #
        head.update('X_PIX', self.cat.x_pix[ii], comment='X pixel in interlaced image')
        head.update('Y_PIX', self.cat.y_pix[ii], comment='Y pixel in interlaced image')
        head.update('MAG', self.cat.mag[ii], comment='MAG_AUTO from interlaced catalog')        
        
        tbHDU = pyfits.new_table(coldefs, header=head)
        tbHDU.writeto(self.baseroot+'_%05d.1D.fits' %(self.id), clobber=True)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print '1D spectrum (%.3f)' %(dt)
        
        #return self.wave, optimal_sum, optimal_sig, trace_spec, trace_sig, scale_to_total
        
        # #### Testing to see how errors look, fit a line and calculate chi**2
        # import gbb.pymc_line 
        # xarr, yarr, earr = self.wave*1., optimal_sum/self.sens, optimal_sig/self.sens
        # mc_line = gbb.pymc_line.init_model(xarr, yarr, earr, order=1)
        # R = pymc.MCMC(mc_line)
        # NSAMP, BURN = 10000, 3000
        # R.sample(NSAMP, burn=BURN)
        # 
        # plt.plot(xarr, yarr, 'o')
        # fit_trace = R.trace('eval_poly')[:]
        # for i in range(0,fit_trace.shape[0],fit_trace.shape[0]/100):
        #     p = plt.plot(xarr, fit_trace[i,:], color='red', alpha=0.05)
        # 
        # 
        # import scipy
        # stats = R.stats()
        # yfit = scipy.polyval([stats['a1']['mean'], stats['a0']['mean']], xarr)
        # dy = (yarr-yfit)/earr
        # plt.hist(dy, range=(-5,5), bins=100)
        # chi2 = np.sum((yarr-yfit)**2/earr**2)/(len(yarr)-1-3)
        
    def extract_spectra_and_diagnostics(self, list=None, skip=True, MAG_LIMIT=19, verbose=False, miny=26, USE_FLUX_RADIUS_SCALE=3, USE_REFERENCE_THUMB=True, largey=None):
        """
        Extract 1D and 2D spectra of objects with id in `list`.
        """
        import unicorn.reduce
        if list is None:
            list = self.cat.id[self.cat.mag < MAG_LIMIT]
        
        #fp = open(self.root+'_inter.failed','w')
        
        N = len(list)
        for i,id in enumerate(list):
            print unicorn.noNewLine+'Object #%-4d (%4d/%4d)' %(id,i+1,N)
            #
            if (os.path.exists(self.baseroot+'_%05d.2D.fits' %(id)) | os.path.exists(self.baseroot+'_%05d.2D.xxx' %(id))) & skip:
                continue
            
            if verbose:
                print '2D FITS'
            status = self.twod_spectrum(id=id, verbose=verbose, miny=miny, refine=True, USE_FLUX_RADIUS_SCALE=USE_FLUX_RADIUS_SCALE, USE_REFERENCE_THUMB=USE_REFERENCE_THUMB)
                        
            if status is False:
                if verbose:
                    print unicorn.noNewLine+'Object #%-4d (%4d/%4d) - No 2D' %(id,i+1,N)
                continue
            
            if largey is not None:
                if verbose:
                    print 'Making large image cutouts.'
                self.twod_spectrum(id=id, verbose=verbose, miny=largey, refine=True, USE_FLUX_RADIUS_SCALE=USE_FLUX_RADIUS_SCALE, USE_REFERENCE_THUMB=USE_REFERENCE_THUMB, BIG_THUMB=True, extract_1d=False)
                    
            if verbose:
                print unicorn.noNewLine+'2D FITS, 2D PNG, 1D FITS'
            self.show_2d(savePNG=True, verbose=verbose)

            if verbose:
                print unicorn.noNewLine+'2D FITS, 2D PNG, 1D FITS, 1D PNG'
            spec = unicorn.reduce.Interlace1D(self.baseroot+'_%05d.1D.fits' %(id), PNG=True)
        
    def make_zeroth_regions(self, MAG_LIMIT=22, id_label=False, scale_size=True, id_list=[]):
            """ 
            Make a regions file for just the zeroth orders. If list is supplied, only make regions for objects 
            in the list. Magnitude limit is applied 
            """
            fp = open(self.baseroot+'_inter_0th.reg','w')
            fp.write('image\n')

            N = len(self.cat.x_pix)

            for i in range(N):
                if self.cat.mag[i] < MAG_LIMIT:
                    if (id_list) and (self.cat['NUMBER'][i] not in id_list):
                        continue
                    else:
                        dx, dy = unicorn.reduce.grism_model(xc_full=self.cat.x_pix[i], yc_full=self.cat.y_pix[i], 
                            zeroth_position=True, growx=self.growx, growy=self.growy, pad=self.pad, 
                            ngrowx=self.ngrowx, ngrowy=self.ngrowy, grism=self.grism_element, finex=self.fine_offset[0], finey=self.fine_offset[1])
                        if scale_size:
                            a, b = self.cat['A_IMAGE'][i], self.cat['B_IMAGE'][i]
                            r = np.sqrt(a**2+b**2)
                        else:
                            r = 3
                        #theta = self.cat['THETA_IMAGE'][i]
                        region = "circle(%s, %s, %.2f)"  %(self.cat.x_pix[i]+dx, self.cat.y_pix[i]+dy, r)
                        if id_label:
                            region += '# text={%d}' %(self.cat['NUMBER'][i])
                    #
                        fp.write(region+'\n')
            fp.close()
    
    def mask_zeroth(self, MAG_LIMIT=22., min_radius=10., only_stars=True, addtocat=True, mask_grism=True):
        """
        Masks 0th order images of stars or of all objects brighter than MAG_LIMIT. 
        If mask_grsim=False, will only create the regions file, but no apply the mask. 
        Requires a catalog with STAR_FLAG to be available.
        """
       
        import pyregion
       
        cat, zout, fout = unicorn.analysis.read_catalogs(root=self.root)
        if (cat == None) and only_stars:
            print "CATALOG NOT FOUND. CANNOT MASK STARS."
            return
       
        star_flag = np.array([cat.star_flag[np.where(cat.id == id)[0][0]] for id in self.cat.id])
        print " {} objects with STAR_FLAG = 1.\n{} objects with STAR_FLAG = 2.".format(np.sum(star_flag==1), np.sum(star_flag==2))
       
        if addtocat:
            self.cat.addColumn(data=star_flag, name='STAR_FLAG')

        star_list = list(self.cat.id[star_flag == 1])
       
        if only_stars:
            self.make_zeroth_regions(MAG_LIMIT=MAG_LIMIT, id_label=False,
                scale_size=False, id_list=star_list)
        else:
            self.make_zeroth_regions(MAG_LIMIT=MAG_LIMIT, id_label=False, scale_size=True)
           
        if mask_grism:
            print "Masking 0th order."
            mask_reg = pyregion.open(self.baseroot+'_inter_0th.reg').as_imagecoord(header=self.gris['SCI'].header)
            print 'Masking {} regions.'.format(len(mask_reg))
       
            for reg in mask_reg:
                if (reg.params[-1] < min_radius) | (reg.coord_list[-1] < min_radius):
                    reg.params[-1] = pyregion.region_numbers.SimpleNumber(min_radius)
                    reg.coord_list[-1] = min_radius
       
            mask = mask_reg.get_mask(hdu=self.gris['SCI'])
            self.gris['SCI'].data[mask] = 0 
    
    
    def refine_mask_background(self, threshold=0.002, grow_mask=8, update=True, resid_threshold=4, clip_left=640, save_figure=True, interlace=True, column_average=False):
        """
        Use the computed model as an aggressive mask for refining the 
        background subtraction.  Also adjust the pixel errors in the second
        extension based on the observed masked pixel distribution.
        
        """
        import re
        
        #### Normalized residuals
        resid = self.gris[1].data / self.gris[2].data
        
        #### Mask on objects in the model and also left part of the image
        #### where model might not have objects
        yc, xc = np.indices(resid.shape)
        objects = nd.maximum_filter((self.model > threshold)*1, size=grow_mask) > 0
        mm = ~objects & (self.gris[2].data > 0) & (np.abs(resid) < resid_threshold) & (self.gris[2].data < 10)
        mm &= (xc > clip_left*(self.growx/2.)) 
        #flux_bg = np.median(self.gris[1].data[mm])
        #err_scale = np.std(resid[mm])
        #resid2 = (self.gris[1].data - flux_bg) / (self.gris[2].data*err_scale)
        #print 'BG: %.5f e/s, err_scale=%.3f' %(flux_bg, err_scale)
        
        #### Compute background for each input exposure
        flux_bgs = np.ones((self.growx, self.growy))
        resid2 = self.gris[1].data*1
        for i in range(self.growx):
            for j in range(self.growy):
                mm_ij = mm[j::self.growy,i::self.growx]
                #### Not enough pixels for this interlace position
                if mm_ij.sum() < 1000:
                    flux_bgs[i,j] = 0.
                else:
                    flux_bgs[i,j] = np.median(self.gris[1].data[j::self.growy,i::self.growx][mm_ij])
                #### Store array with new background subtracted
                resid2[j::self.growy,i::self.growx] -= flux_bgs[i,j]
        
        #### Scale factor to make the normalized residuals std=1
        err_scale = np.std((resid2/self.gris[2].data)[mm])
        
        ### Formatted string for multiple image residuals
        np.set_printoptions(precision=4)
        flux_bgs_str = re.sub('[\[\],\n]', '', (flux_bgs.__str__())).replace('  ',' ')
        np.set_printoptions(precision=8)
        
        print 'BG:%s e/s\nerr_scale=%.3f' %(flux_bgs_str, err_scale)
        
        xroot = os.path.basename(self.root)+'-'+self.grism_element
        
        ### Save information to a file
        fp = open(xroot+'_maskbg.dat','w')
        np.set_printoptions(precision=3)
        bg_head = ' '.join(['bg%d%d' %(i, j) for j in range(self.growx) for i in range(self.growy)])
        fp.write('# root %s err_scale\n%s %s %.5f\n' %(bg_head, xroot, flux_bgs_str, err_scale))
        np.set_printoptions(precision=8)
        fp.close()
        
        NX = mm.shape[1]
        column = np.zeros(NX)
        if column_average:
            """
            Refine subtracted column average background 
            """
            #msk = ~objects & (self.gris[2].data > 0) & (np.abs(resid) < resid_threshold) & (self.gris[2].data < 10)
            print 'Refine: subtract column average'
            cleaned = self.gris[1].data*mm
            for i in range(NX):
                subset = cleaned[:,np.maximum(i-20,0):i+20]
                ok = (subset != 0)
                if ok.sum() == 0:
                    continue
                #
                ok = ok & (subset > np.percentile(subset[ok], 16)) & (subset < np.percentile(subset[ok], 84))
                column[i] = np.mean(subset[ok])
            
            resid2 -= column
            
        ### Make plot
        if save_figure:
            fig = unicorn.plotting.plot_init(xs=8, aspect=0.5, left=0.01, right=0.01, top=0.01, bottom=0.1, square=True, NO_GUI=True)
            ax = fig.add_subplot(121)
            diff = self.gris[1].data*1.-column
            diff[~mm] = 0
            diff = diff[self.ngrowx*self.growx: -self.ngrowx*self.growx-1, self.ngrowy*self.growy:-self.ngrowy*self.growy-1]
            ax.imshow(nd.gaussian_filter(diff, 2), vmin=-0.02, vmax=0.02, interpolation='Nearest', aspect='auto')
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.set_xlabel('BG: %s e/s\n' %(flux_bgs_str)+r'$\sigma$ scale: %.3f' %(err_scale))
            if column_average:
                xx = self.ngrowx*self.growx
                ax.plot(column[xx:-xx-1]*80000 + diff.shape[0]/2., color='black')
            ax.set_xlim(0,diff.shape[1]); ax.set_ylim(0,diff.shape[0])
            #
            ax = fig.add_subplot(122)
            ax.hist((resid2/self.gris[2].data)[mm], range=(-3,3), bins=60, alpha=0.5, label='Observed', histtype='stepfilled', color='black')
            ax.hist(np.random.normal(size=mm.sum()), range=(-3,3), bins=60, alpha=0.8, color='red', histtype='step', linewidth=3)
            ax.set_ylabel(xroot)
            ax.set_yticklabels([])
            unicorn.plotting.savefig(fig, xroot+'_maskbg.png')
                
        ### Update the interlaced grism image in place
        if update:
            mx = pyfits.open(self.gris.filename(), mode='update')
            
            for i in range(self.growx):
                for j in range(self.growy):
                    key = 'FLUXBG%d%d' %(i,j)
                    if key in mx[0].header.keys():
                        mx[0].header[key] += flux_bgs[i,j]
                    else:
                        mx[0].header.update(key,flux_bgs[i,j])
            
            if 'ERRSCALE' in mx[0].header.keys():
                mx[0].header['ERRSCALE'] *= err_scale
            else:
                mx[0].header.update('ERRSCALE', err_scale)
            
            resid2[mx[2].data == 0] = 0
                            
            mx[1].data = resid2*1
            mx[2].data *= err_scale
            mx.flush()
            print 'Updated %s' %(self.gris.filename())
            
            self.gris = mx
            
    def run_lacosmic(self, update=True):
        """
        Run M. Tewes' Python version of LACosmic to flag bad remaining
        bad pixels.
        
        http://obswww.unige.ch/~tewes/cosmics_dot_py/
        
        """
        import scipy.ndimage as nd

        import lacosmic
        
        EXPTIME = self.gris[0].header['EXPTIME']
        SKY = self.gris[0].header['SKYSCALE']

        ### Fill zeros in interlaced grism
        orig = self.gris[1].data*1
        img = self.gris[1].data*1
        ma = nd.maximum_filter(img, size=3)
        mask = (self.gris['ERR'].data == 0)
        img[mask] = ma[mask]
        
        x = lacosmic.cosmicsimage(img*EXPTIME, gain=1, readnoise=20/self.growx/self.growy, pssl=SKY*EXPTIME, objlim=5.0, sigclip=5., sigfrac=1)
        
        x.run(verbose=True)
        self.cr_mask = x.mask
        
        if update:
            mx = pyfits.open(self.gris.filename(), mode='update')
            mx[0].header.update('LACOSMIC', 1, 'Additional CRs flagged with LACosmic')
            mx['SCI'].data[x.mask] = 0
            mx['ERR'].data[x.mask] = 0
            mx.flush()
    
    def align_direct_reference(self, match_threshold=0.5, mag_limit=23., r_limit=8, apply=False):
        """
        Cross-correlation alignment of blotted direct and nominal direct image
        """
        from threedhst import catIO
        
        #### Run SExtractor on the direct image, with the WHT 
        #### extension as a weight image
        se = threedhst.sex.SExtractor()
        
        ## Set the output parameters required for aXe 
        ## (stored in [threedhst source]/data/aXe.param) 
        se.aXeParams()
        
        ## XXX add test for user-defined .conv file
        se.copyConvFile()
        se.overwrite = True
        se.options['CATALOG_NAME']    = '%s_align0.cat' %(self.baseroot)
        se.options['CHECKIMAGE_TYPE'] = 'NONE'
        se.options['WEIGHT_TYPE']     = 'MAP_RMS'
        se.options['WEIGHT_IMAGE']    = '%s[1]' %(self.direct.filename())
        se.options['FILTER']    = 'Y'
        
        #### Detect thresholds (default = 1.5)
        se.options['THRESH_TYPE'] = 'ABSOLUTE'
        se.options['DETECT_THRESH']    = '0.2' 
        se.options['ANALYSIS_THRESH']  = '0.2'
        se.options['MAG_ZEROPOINT'] = '%.2f' %(self.ZP)
        
        #### Run SExtractor
        status = se.sextractImage('%s[0]' %(self.direct.filename()))

        se.options['CATALOG_NAME']    = '%s_align1.cat' %(self.baseroot)
        status = se.sextractImage('%s[0]' %(self.im.filename()))
        
        c0 = catIO.Table('%s_align0.cat' %(self.baseroot))
        c1 = catIO.Table('%s_align1.cat' %(self.baseroot))
        
        ok_obj = (c0['MAG_AUTO'] < mag_limit) & (c0['FLUX_RADIUS'] < r_limit)
        
        m = catIO.CoordinateMatcher(c0[ok_obj], 'X_IMAGE', 'Y_IMAGE', pixel_units=True)
        dr, idx = m.match_list(c1['X_IMAGE'], c1['Y_IMAGE'])
        ok = dr < match_threshold
        dx = c0['X_IMAGE'][ok_obj][idx][ok] - c1['X_IMAGE'][ok]
        dy = c0['Y_IMAGE'][ok_obj][idx][ok] - c1['Y_IMAGE'][ok]
        
        self.fine_offset = [np.median(dx), np.median(dy)]
        
        #### Make the alignment figure
        plt.ioff()
        fig = plt.figure(figsize=[5,5])
        ax = fig.add_subplot(111)

        ax.scatter(dx, dy, alpha=0.5)

        ax.plot([dx.min(), dx.max()], [0,0], color='black', linestyle=':', alpha=0.5)
        ax.plot([0,0], [dy.min(), dy.max()], color='black', linestyle=':', alpha=0.5)
        ax.plot([dx.min(), dx.max()], np.median(dy)*np.ones(2), color='red', linestyle='-', alpha=0.5)
        ax.plot(np.median(dx)*np.ones(2), [dy.min(), dy.max()], color='red', linestyle='-', alpha=0.5)

        ax.set_xlabel(r'$\Delta x$ (interlaced pixels)')
        ax.set_ylabel(r'$\Delta y$ (interlaced pixels)')

        ax.set_title('%s: (%.3f, %.3f)' %(self.baseroot, self.fine_offset[0], self.fine_offset[1]))

        fig.tight_layout()
        
        fig.savefig('%s_ref_align.png' %(self.baseroot))
        
        threedhst.showMessage('Align %s to %s:  \n\n   (dx, dy) = (%.3f, %.3f)  \n\n See %s_ref_align.png' %(self.im.filename(), self.direct.filename(), self.fine_offset[0], self.fine_offset[1], self.baseroot))
        
        #### Put alignment results in ref_inter FITS file
        if apply:
            im = pyfits.open(self.im.filename(), mode='update')
            im[0].header['FINEX'] = (self.fine_offset[0], 'Fine alignment x offset')
            im[0].header['FINEY'] = (self.fine_offset[1], 'Fine alignment y offset')
            im.flush()
            
            self.im = pyfits.open(self.im.filename())
            
        return self.fine_offset
        
    def get_eazy_templates(self, dr_min=0.2, MAG_LIMIT=23):
        """
        Get EAZY templates for objects in the field.  
        
        Returns a "model_list" dictionary suitable for input into 
        `processGrismModel`.
        """
        from threedhst import eazyPy as eazy
        from threedhst import catIO
        
        cat, zout, fout = unicorn.analysis.read_catalogs(root=self.root)
        #cat.kmag = 23.86-2.5*np.log10(cat.field(KTOT_COL))
        CAT_PATH = os.path.dirname(cat.filename)
        
        root = os.path.basename(zout.filename).split('.zout')[0]
        ZOUT_PATH = os.path.dirname(zout.filename)
        
        eazy_tempfilt, eazy_coeffs, eazy_temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE = root,                                                 OUTPUT_DIRECTORY=ZOUT_PATH, CACHE_FILE = 'Same')
        
        m = catIO.CoordinateMatcher(cat)
        dr, idx = m.match_list(self.cat.ra, self.cat.dec)
        ok = (dr < dr_min) & (self.cat.mag < MAG_LIMIT)
                            
        #### Read the templates
        nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed1_nolines.dat', unpack=True)
        NTEMP = eazy_coeffs['coeffs'].shape[0]    
        noline_temps = np.zeros((nlx.shape[0], NTEMP))
        noline_temps[:,0] = nly/eazy_coeffs['tnorm'][0]
        for i in range(2,7):
            nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed%d_nolines.dat' %(i), unpack=True)
            noline_temps[:,i-1] = nly/eazy_coeffs['tnorm'][i-1]
        
        #### The last v1.1 template was a BC03 model without lines, so just use it directly    
        i = 7
        lx, ly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.1_lines/eazy_v1.1_sed%d.dat' %(i), unpack=True)
        noline_temps[:,6] = np.interp(nlx, lx, ly)/eazy_coeffs['tnorm'][6]
        
        ### dusty old template
        if NTEMP == 8:
            lx, ly = np.loadtxt(unicorn.GRISM_HOME+'/templates/Ezgal/c09_del_8.6_z_0.019_chab_age09.40_av2.0.dat', unpack=True)
            noline_temps[:,7] = np.interp(nlx, lx, ly)/eazy_coeffs['tnorm'][7]
        
        param = eazy.EazyParam(zout.filename.replace('zout','param'))
        flam_factor = 10**(-0.4*(param['PRIOR_ABZP']+48.6))*3.e18/1.e-17
        
        ### Loop through objects
        model_list = {}
        for i in np.arange(len(idx))[ok]:
            print unicorn.noNewLine + 'ID: %d, mag: %.1f' %(self.cat.id[i], self.cat.mag[i])
            id = self.cat.id[i]
            best_fit_nolines = np.dot(noline_temps, eazy_coeffs['coeffs'][:,idx[i]])
            templam_nolines = nlx
            
            #### Convert to flux
            zi = eazy_tempfilt['zgrid'][eazy_coeffs['izbest'][idx[i]]]

            ###### Full template SED, observed frame
            lambdaz = nlx*(1+zi)
            temp_sed = np.dot(noline_temps, eazy_coeffs['coeffs'][:,idx[i]])
            temp_sed /= (1+zi)**2
            temp_sed *= (1/5500.)**2*flam_factor

            ###### IGM absorption
            lim1 = np.where(nlx < 912)
            lim2 = np.where((nlx >= 912) & (nlx < 1026))
            lim3 = np.where((nlx >= 1026) & (nlx < 1216))

            if lim1[0].size > 0: temp_sed[lim1] *= 0.
            if lim2[0].size > 0: temp_sed[lim2] *= 1.-eazy_temp_sed['db'][eazy_coeffs['izbest'][idx[i]]]
            if lim3[0].size > 0: temp_sed[lim3] *= 1.-eazy_temp_sed['da'][eazy_coeffs['izbest'][idx[i]]]
            
            model_list[id] = (lambdaz, temp_sed)

        return model_list
        
        if False:
            # Test
            model_list = model.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
            id = 335
            obj = model.object*0.
            N = len(model_list.keys())
            for i, id in enumerate(model_list.keys()):
                print '%d/%d %d' %(i, N, id)
                model.compute_object_model(id, lam_spec=model_list[id][0], flux_spec=model_list[id][1], verbose=False, normalize=True)
                obj += model.object
                
def field_dependent(xi, yi, str_coeffs):
    """ 
    Calculate field-dependent parameter for aXe conventions.
    
    a = a_0 + a_1 * xi + a_2 * yi + a_3 * xi**2 + a_4 * xi * yi + a_5 * yi**2
    
    """
    coeffs = np.cast[float](str_coeffs.split())
    xy = np.array([1,xi,yi,xi**2,xi*yi,yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3])
    a = np.sum(coeffs*xy[0:len(coeffs)])
    return a

def model_stripe():
    """
    Model a horizontal stripe across the image assuming uniform illumination
    to compare to the sky background images
    """
    import unicorn.reduce as red
    yran = range(800,900)
    yran = range(1014)
    
    model = np.zeros((1014,1014))
    lam_spec, flux_spec = None, None
    f125 = 0.68
    f160 = 0.6
    f160 *= 10**(-(25.96-26.25)/2.5)/(1.6/1.25)**2
    lam_spec = np.array([0.9e4,1.25e4,1.6e4,1.9e4])
    flux_spec = np.array([1.,1.,f160/f125,f160/f125])
    
    skip = 5
    
    BEAMS = ['A','B']
    BEAMS = ['A','B','C','D'] #,'E']
    BEAMS = ['A','B','C','D','E']
    
    model*=0
    
    orders, xi = red.grism_model(507,507, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0)
    yord, xord = np.indices(orders.shape)
    non_zero = orders > 0
    xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
        
    ys = orders.shape
    xord += xi[0]
    yord -= (ys[0]-1)/2
        
    for y in range(490, 510):
        print unicorn.noNewLine + 'Y: %d' %(y)
        #if (y % 20) == 0:
        #    pyfits.writeto('stripe.fits', model/model.max(), clobber=True)
        #
        # Range tuned to have the bands in the right place
        for x in range(-190,1014+85):
            if ((x % skip) + (y % skip)) == 0:
                orders, xi = red.grism_model(x+1, y+1, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0, dydx=False)
                yord, xord = np.indices(orders.shape)
                non_zero = orders > 0
                xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
                ys = orders.shape
                xord += xi[0]
                yord -= (ys[0]-1)/2
            #
            xxi = x+xord
            yyi = y+yord
            use = (xxi >= 0) & (xxi < 1014) & (yyi >= 0) & (yyi < 1014)
            model[yyi[use], xxi[use]] += ford[use]
        #
        #model += object #*norm
    #
    pyfits.writeto('stripe.fits', model/model.max(), clobber=True)
    #
    ####### Compare order contributions
    BEAMS = []
    flux_spec = (lam_spec/1.4e4)**-18
    for beam in ['B','C','D','E']:
        BEAMS = ['A']
        BEAMS.append(beam)
        orders, xi = red.grism_model(507,507, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0, xmi=-600, xma=700)
        orders_1d = np.sum(orders, axis=0)
        model_1d = orders_1d*0.
        for x in range(-190,1014+85):
            model_1d += nd.shift(orders_1d, xi[0]+x, mode='constant', cval=0)
        #
        plt.plot(model_1d[:1014], label=beam)
    
    plt.legend()
    
    ####### Compare vertical dependence
    BEAMS = ['A','B','C','D','E']
    for yi in range(0,1014,100):
        print yi
        orders, xi = red.grism_model(507, yi, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0, xmi=-600, xma=700)
        orders_1d = np.sum(orders, axis=0)
        model_1d = orders_1d*0.
        for x in range(-190,1014+85):
            model_1d += nd.shift(orders_1d, xi[0]+x, mode='constant', cval=0)
        #
        plt.plot(model_1d[:1014]/model_1d[800], label='y = %d' %(yi))
    
    ####### Compare background spectrum color
    BEAMS = ['A','B','C','D','E']
    for beta in range(-20,3,2):
        flux_spec = (lam_spec/1.4e4)**beta
        orders, xi = red.grism_model(507, 507, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0, xmi=-600, xma=700)
        orders_1d = np.sum(orders, axis=0)
        model_1d = orders_1d*0.
        for x in range(-190,1014+85):
            model_1d += nd.shift(orders_1d, xi[0]+x, mode='constant', cval=0)
        #
        plt.plot(model_1d[:1014]/model_1d[800], label='%d' %(beta))
    #
    ####### Compare background emission lines
    BEAMS = ['A','B','C','D'] #,'E']
    for l0 in np.arange(1.e4,1.61e4,0.1e4):
        flux_spec = (lam_spec/1.4e4)**-5*0.1+np.exp(-(lam_spec-l0)**2/2./100.**2)
        orders, xi = red.grism_model(507, 507, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, growx=1, growy=1, pad=0, xmi=-600, xma=700)
        orders_1d = np.sum(orders, axis=0)
        model_1d = orders_1d*0.
        for x in range(-190,1014+85):
            model_1d += nd.shift(orders_1d, xi[0]+x, mode='constant', cval=0)
        #
        plt.plot(model_1d[:1014]/model_1d[800], label='%f' %(l0/1.e4))
            
def interlace_combine_blot(root='COSMOS-19-F140W', view=True, pad=60, REF_ROOT = 'COSMOS_F160W', CATALOG='UCSC/catalogs/COSMOS_F160W_v1.cat',  NGROWX=180, NGROWY=30, verbose=True, growx=2, growy=2, auto_offsets=False, ref_exp=0, NSEGPIX=8, stop_early=False, grism='G141', old_filenames=False, reference_pixel=None):
    """
    Combine blotted image from the detection mosaic as if they were 
    interlaced FLT images
    
    xxxxx Still not modified for differing factors in x and y!
    """
    if not stop_early:
        from pyraf import iraf

    import threedhst.prep_flt_files
    from threedhst import catIO
    
    import unicorn
    
    if unicorn.hostname().startswith('uni'):
        view = False
    
    ### Strip off filter name    
    pointing = '-'.join(root.split('-')[:-1])
    
    if view:
        import threedhst.dq
        ds9 = threedhst.dq.myDS9()
    
    if not stop_early:
        run = threedhst.prep_flt_files.MultidrizzleRun(root)

        scl = np.float(run.scl)
        xsh, ysh = threedhst.utils.xyrot(np.array(run.xsh)*scl, np.array(run.ysh)*scl, run.rot[0])

    #### Pixel area map
    im = pyfits.open(os.getenv('iref')+'ir_wfc3_map.fits')
    PAM = im[1].data
    im.close()
    PAMy = np.median(PAM[:,50:-50], axis=1)
    PAM_grow = np.dot(PAMy.reshape((-1,1)), np.ones((1,1014+2*NGROWX)))
    PAM_grow = np.vstack([np.array([PAM_grow[0,:],]*NGROWY),PAM_grow, np.array([PAM_grow[-1,:],]*NGROWY)])

    yi,xi = np.indices((1014+2*NGROWY,1014+2*NGROWX))

    inter_sci = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX))
    inter_seg = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX), dtype=np.int32)
    inter_N = np.zeros((1014*growy+pad+growy*2*NGROWY, 1014*growx+pad+growx*2*NGROWX), dtype=np.int32)
    
    xi+=pad/(2*growx)
    yi+=pad/(2*growy)#+NGROW
    
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    
    #### From pieter
    dxs = np.array([0,-20,-13,7]) #+ np.int(np.round(xsh[0]))*0
    dys = np.array([0,-7,-20,-13]) #+ np.int(np.round(ysh[0]))*0
    
    #### GOODS-N field from Weiner program 11600
    if root.startswith('GOODS-N') | root.startswith('GNGRISM') | root.startswith('goodsn'):
        dxs = np.array([0,-9,-4,5]) #+ np.int(np.round(xsh[0]))*0
        dys = np.array([0,-3,-11,-8]) #+ np.int(np.round(ysh[0]))*0
    
    if root.startswith('TILE41'):
        dxs = np.array([-4, 9, -8, 5, 5, -8]) #+ np.int(np.round(xsh[0]))*0
        dys = np.array([-5, -4, 4, 5, 5, 4]) #+ np.int(np.round(ysh[0]))*0
        
    #### Erb QSO sightlines
    if flt[0].header['PROPOSID'] == 12471:
        dxs, dys = np.array([  0, -10])*growx, np.array([ 0, -7])
    
    if 'ERSII' in root:
        dxs = np.array([   0, -147, -148,-1])*growx
        dys = np.array([ 0,  0, 83, 83])*growy
        
    if auto_offsets:
        try:
            xinter, yinter = unicorn.reduce.wcs_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy, reference_pixel = reference_pixel)
        except:
            xinter, yinter = unicorn.reduce.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy)
        #xinter, yinter = unicorn.reduce.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy)
        dxs = xinter #+ np.int(np.round(xsh[0]))*0
        dys = yinter #+ np.int(np.round(ysh[0]))*0
        
        ## Check
        if False:
            xx, yy = unicorn.reduce.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy, raw=True)
            flts = []
            for exp in asn.exposures:
                flts.append(pyfits.open(exp+'_flt.fits')[1].data)
            #
            i=1
            xx, yy = xx-xx[0], yy-yy[0]
            ds9.view(flts[0]-nd.shift(flts[i], (-yy[i], -xx[i])))

    dxi = np.cast[np.int32](np.ceil(dxs/growx))
    dyi = np.cast[np.int32](np.ceil(dys/growy))
    
    dxi -= dxi[ref_exp]
    dyi -= dyi[ref_exp]

    dxs -= dxs[ref_exp]
    dys -= dys[ref_exp]

    print 'Interlace offsets: ', dxs, dys
    
    #### Find hot pixels, which are flagged as cosmic 
    #### rays at the same location in each of the flt files.  Ignore others.
    # hot_pix = np.zeros((1014,1014),dtype='int')
    # for flt in run.flt:
    #     im = pyfits.open(flt+'.fits')
    #     hot_pix += (im[3].data & 4096) / 4096
    #     
    # hot_pix = hot_pix > 2
        
    ### For COSMOS-2, somehow the order of exposures isn't sorted 
    ### alphabetically.  The multidrizzle ".run" file has the files in alpha
    ### order, so if you loop through them they won't correspond to the 
    ### correct dither positions.  Use the ASN file itself for the list
    ### of flt/blot files
    #for i,flt in enumerate(run.flt):
    #   flt = run.flt[i]
    for i in range(len(asn.exposures)):
        flt = asn.exposures[i]+'_flt'
        print flt
        im = pyfits.open(flt.replace('_flt','_blot')+'.fits')
        #im_wht = pyfits.open(flt.replace('_flt','_blot_wht')+'.fits')
        im_seg = pyfits.open(flt.replace('_flt','_seg')+'.fits')
        
        #### Use the pixel area map correction to preserve flux in FLT frame
        #print im[0].data.shape, PAM_grow.shape
        im[0].data *= PAM_grow
        
        #### Divide by 4 to conserve surface brightness with smaller output pixels
        im[0].data /= (growx*growy)
        #im_wht[0].data /= 4
        #im[2].data /= 4
        
        ### Mask cosmic rays
        if i == 0:
            im_flt = pyfits.open(flt+'.fits')
            h0 = im_flt[0].header
            h1 = im_flt[1].header
            h0['FILTER'] = im[0].header['FILTER']
            header = unicorn.reduce.scale_header_wcs(h1.copy(), factor=2, growx=growx, growy=growy, pad=pad, NGROWX=NGROWX, NGROWY=NGROWY)
            
            #header_wht = header.copy()
            #header_wht.update('EXTNAME','ERR')
            header['EXTNAME'] = 'SCI'
            header.update('REF_ROOT', REF_ROOT, comment='Source detection image')
            
        # dx = np.int(np.round((xsh[i]-xsh[0])*2))
        # dy = np.int(np.round((ysh[i]-ysh[0])*2))
        #
        dx = dxs[i]
        dy = dys[i]
        #
        #use = ((im[3].data & 4096) == 0) & ((im[3].data & 4) == 0) #& (xi > np.abs(dx/2)) & (xi < (1014-np.abs(dx/2))) & (yi > np.abs(dy/2)) & (yi < (1014-np.abs(dy/2)))
        #        
        #use = ((im[3].data & (4+32+16+2048+4096)) == 0) & (~hot_pix)
        print dy, dx
        inter_sci[yi*growy+dy,xi*growx+dx] += im[0].data
        #inter_err[yi*2+dy,xi*2+dx] += im_wht[0].data
        inter_seg[yi*growy+dy,xi*growx+dx] = im_seg[0].data
        inter_N[yi*growy+dy,xi*growx+dx] += 1
        
        #
        if view:
            ds9.view_array(inter_sci, header=header)
            ds9.scale(-0.1,5)
    
    inter_N[inter_N == 0] = 1
    
    #### Write interlaced sci/wht image
    #header.update('NGROW', NGROW, comment='Number of pixels added to X-axis (centered)')        
    #header.update('PAD', pad, comment='Additional padding around the edges') 
    
    header.update('PAD', pad, comment='Padding at edge of interlaced image')
    header.update('NGROWX', NGROWX, comment='Additional extra pixels at the image edges')
    header.update('NGROWY', NGROWY, comment='Additional extra pixels at the Y image edges')
    
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=inter_sci/inter_N, header=header)
    #wht = pyfits.ImageHDU(data=inter_err, header=header_wht)
            
    image = pyfits.HDUList([hdu, sci]) #, wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
    
    image.writeto(pointing+'-%s_ref_inter.fits' %(grism), clobber=True, output_verify="fix")
    
    #### Write interlaced segmentation image
    ## First clean up regions between adjacent objects that get weird segmentation values
    ## from blot
    # yh, xh = np.histogram(inter_seg, range=(0,inter_seg.max()), bins=inter_seg.max())
    # minseg = NSEGPIX
    # bad = (yh < minseg) & (yh > 0)
    # if verbose:
    #     print 'Clean %d objects with < %d segmentation pixels.  \nThese are probably artificial resulting from blotting the integer segmentation image.\n' %(bad.sum(), minseg)
    # for val in xh[:-1][bad]:
    #     if verbose:
    #         print unicorn.noNewLine+'%d' %(val)
    #     test = inter_seg == val
    #     if test.sum() > 0:
    #         inter_seg[test] = 0
        
    #### Clean up image, can get zeros and NaNs from above
    inter_seg[inter_seg < 0] = 0
    inter_seg[~np.isfinite(inter_seg)] = 0
        
    hdu = pyfits.PrimaryHDU(header=header, data=np.cast[np.int32](inter_seg))
    if old_filenames:
        seg_file = pointing + '_inter_seg.fits'
        cat_file = pointing + '_inter.cat'
    else:
        seg_file = pointing+'-%s_inter_seg.fits' %(grism)
        cat_file = pointing+'-%s_inter.cat' %(grism)
        
    hdu.writeto(seg_file, clobber=True)
    
    #### For use with astrodrizzle images: stop here and return offsets
    if stop_early:
        return dxs, dys
         
    if verbose:
        print 'Clean up segmentation image...'
        
    #### xxx not necessary now when using "nearest" with blot!!!
    #### Clean up overlap region of segmentation 
    unicorn.reduce.fill_inter_zero(seg_file)
    im = pyfits.open(seg_file) #, mode='update')
    inter_seg = im[0].data
    
    #s = np.ones((3,3))
    #labeled_array, num_features = nd.label(im[0].data, structure=s)
    
    #if verbose:
    #    print unicorn.noNewLine+'Clean up segmentation image...[1]'

    #### Maximum filter to flag overlap regions
    # max_filter = nd.maximum_filter(im[0].data, size=(3,3))
    # bad = (max_filter - im[0].data) > 2
    # inter_seg = max_filter*1
    # inter_seg[bad] = 0
    # 
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[2]'

    #### Orphan pixels in the max-filtered image
    # npix = nd.convolve((inter_seg > 0)*1, np.ones((3,3)))
    # bad = (inter_seg > 0) & (npix < 4)
    # inter_seg[bad] = 0
    # im[0].data = inter_seg
    # im.flush()
    
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[3]'

    #### Make a version of the catalog with transformed coordinates and only those objects
    #### that fall within the field    
    old_cat = threedhst.sex.mySexCat(CATALOG)
    objects_in_seg = np.unique(inter_seg)
    # pop_id = []
    # for id in old_cat.id[::-1]:
    #     if id not in objects_in_seg:
    #         #print 'Pop #%05d' %(id)
    #         pop_id.append(id)
    # 
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[3b]'
    
    ### Just rewrite a temporary catalog rather than using "popItem", which is very slow on large catalogs
    new_lines = []
    for i, id in enumerate(old_cat.id):
        if id in objects_in_seg:
            new_lines.append(old_cat.rowlines[i])
    
    fp = open('/tmp/%s.tmp.cat' %(root), 'w')
    fp.writelines(old_cat.headerlines)
    fp.writelines(new_lines)
    fp.close()
    
    old_cat = threedhst.sex.mySexCat('/tmp/%s.tmp.cat' %(root))
    
    #old_cat.popItem(np.array(pop_id), verbose=False)

    # import astropy.io.ascii as aa
    # old_cat = aa.SExtractor().read(CATALOG)
    # objects_in_seg = np.unique(inter_seg)
    # N = old_cat['NUMBER']
    
    # fp = open("/tmp/%s.drz_xy" %(root),'w')
    # #fpr = open("/tmp/%s.flt_xy.reg" %(root),'w')
    # xw, yw = old_cat['X_IMAGE'], old_cat['Y_IMAGE']
    # NOBJ = len(old_cat.id)
    # for i in range(NOBJ):
    #     fp.write('%.3f %.3f\n' %(xw[i], yw[i]))
    # fp.close()
    # 
    # if verbose:
    #     print 'Convert from reference coords to the distorted flt frame with iraf.wblot'
    
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[4]'

    #flt = run.flt[0]+'.fits'
    flt = asn.exposures[0]+'_flt.fits'
    
    threedhst.process_grism.flprMulti()
    #status = iraf.tran(origimage=flt+'[sci,1]', drizimage=root+'_drz_sci.fits', direction="backward", x=None, y=None, xylist="/tmp/%s.drz_xy" %(root), mode="h", Stdout=1)
    
    this_drz_head = pyfits.getheader(root+'_drz.fits', 1)
    # wcs_this = wcsutil.WCSObject('dummy', header=this_drz_head)
    print root, pywcs.__file__
    wcs_this = pywcs.WCS(header=this_drz_head)
    
    fp = open("/tmp/%s.drz_xy" %(root),'w')
    fpr = open("/tmp/%s.drz_xy.reg" %(root),'w')
    fpr2 = open("/tmp/%s.drz_rd.reg" %(root),'w')
    fpr.write('image\n')
    fpr2.write('fk5\n')
    NOBJ = len(old_cat.id)
    for i in range(NOBJ):
        #xw, yw = wcs_this.rd2xy((old_cat.ra[i], old_cat.dec[i]))
        try:
            xw, yw = wcs_this.wcs_world2pix([[old_cat.ra[i], old_cat.dec[i]]],1)[0]
        except:
            xw, yw = wcs_this.wcs_sky2pix([[old_cat.ra[i], old_cat.dec[i]]],1)[0]
            
        fp.write('%.2f %.2f\n' %(np.clip(xw, -1000, 4999), np.clip(yw, -1000, 4999)))
        fpr.write('circle(%.2f,%.2f,1")\n' %(np.clip(xw, -1000, 4999), np.clip(yw, -1000, 4999)))
        fpr2.write('circle(%.6f,%.6f,1") # color=magenta\n' %(old_cat.ra[i], old_cat.dec[i]))
    
    fp.close()
    fpr.close()
    fpr2.close()
    
    status = iraf.tran(origimage=flt+'[sci,1]', drizimage=root+'_drz.fits[1]', direction="backward", x=None, y=None, xylist="/tmp/%s.drz_xy" %(root), mode="h", Stdout=1)
    #### To avoid re-running multidrizzle, need to get xy coordinates in the original drz image,
    #### and *then* run tran
    
    ###xxx testing    
    if len(status) < NOBJ:
        threedhst.showMessage('Something went wrong with `iraf.tran`.', warn=True)
        for line in status:
            if len(line.split()) != 4:
                print line
        return False
        
    #os.remove("/tmp/%s.drz_xy" %(root))
    
    print unicorn.noNewLine+'Convert from reference coords to the distorted flt frame with iraf.wblot - Done.'
    
    #delta = catIO.Readfile('%s_delta.dat' %(root), save_fits=False)
    #delta_x, delta_y = delta.xoff[-1], delta.yoff[-1]
    delta_x, delta_y = 0, 0
    
    flt_x, flt_y, popID = [], [], []
    fpr = open("%s-%s_inter.reg" %(pointing, grism),'w')
    fpr.write("global color=green dashlist=8 3 width=1 font=\"helvetica 8 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\n")
    for i,line in enumerate(status[-NOBJ:]):
        #print line
        spl = np.cast[float](line.split())
        fxi, fyi = (spl[0]+NGROWX)*growx+pad/2-1, (spl[1]+NGROWY)*growy+pad/2-1
        if (fxi > pad/2) & (fxi < ((1014+2*NGROWX)*growx+pad/2.)) & (fyi > (pad/2.+NGROWY*growy)) & (fyi < (pad/2.+NGROWY*growy+1014*growy)):
            fxi += delta_x
            fyi += delta_y
            flt_x.append(fxi)
            flt_y.append(fyi)
            #
            fpr.write('ellipse(%.2f,%.2f,%.1f,%.1f) # text={%d}\n' %(fxi, fyi, 2.5*growx, 2.5*growy, old_cat.id[i]))
        else:
            popID.append(old_cat.id[i])
    
    if len(popID) > 0:
        old_cat.popItem(np.array(popID))               
    
    fpr.close()
    
    if verbose:
        print 'Make the catalog file: %s' %(cat_file)
        
    ### Put the x,y coordinates in the FLT frame into the SExtractor catalog
    status = old_cat.addColumn(np.array(flt_x), name='X_FLT', comment='X pixel in FLT frame', verbose=True)
    status = old_cat.addColumn(np.array(flt_y), name='Y_FLT', comment='Y pixel in FLT frame', verbose=True)
    old_cat.write(outfile=cat_file)
        
def prepare_blot_reference(REF_ROOT='COSMOS_F160W', filter='F160W', REFERENCE = 'UCSC/cosmos_sect11_wfc3ir_F160W_wfc3ir_drz_sci.fits', SEGM = 'UCSC/checkimages/COSMOS_F160W_v1.seg.fits', Force=False, sci_extension=0):
    """
    Need to strip some header keywords from Multidrizzled images and make sure 
    certain keywords exist (EXPTIME, CD1_2, etc.)
    """
    import shutil
    
    im_ref = pyfits.open(REFERENCE)
    
    new_head = unicorn.reduce.strip_header(im_ref[sci_extension].header, filter=filter)
    
    if (not os.path.exists(REF_ROOT+'_ref.fits')) | Force:
        ref = pyfits.PrimaryHDU(header=new_head, data=im_ref[sci_extension].data)
        ref.writeto(REF_ROOT+'_ref.fits', output_verify='fix', clobber=True)
        
    im_seg = pyfits.open(SEGM)
    
    if (not os.path.exists(REF_ROOT+'_seg.fits')) | Force:
        seg = pyfits.PrimaryHDU(header=new_head, data=np.cast[np.float32](im_seg[0].data))
        seg.writeto(REF_ROOT+'_seg.fits', output_verify='fix', clobber=True)

    #shutil.copy(SEGM, REF_ROOT+'_ones.fits')

    if (not os.path.exists(REF_ROOT+'_ones.fits')) | Force:
        test = im_seg[0].data > 0
        test *= 1.
        test = np.cast[np.float32](test)
        ones = pyfits.PrimaryHDU(data=test, header=new_head)
        ones.writeto(REF_ROOT+'_ones.fits', output_verify='fix', clobber=True)
        
    print '\n\n --- Can ignore "currupted HDU" warnings ---\n\n'
    
    fp = open(REF_ROOT+'.info','w')
    fp.write("%s_ref.fits  %s\n" %(REF_ROOT, REFERENCE))
    fp.write("%s_seg.fits  %s\n" %(REF_ROOT, SEGM))
    fp.close()
    
def adriz_blot_from_reference(pointing='cosmos-19-F140W', pad=60, NGROWX=180, NGROWY=30, growx=2, growy=2, auto_offsets=False, ref_exp=0, ref_image='Catalog/cosmos_3dhst.v4.0.IR_orig_sci.fits', ref_ext=0, ref_filter='F140W', seg_image='Catalog/cosmos_3dhst.v4.0.F160W_seg.fits', cat_file='Catalog/cosmos_3dhst.v4.0.IR_orig.cat', ACS=False, grism='G141', reference_pixel=None):
    """
    Use AstroDrizzle to blot reference and sci images and SExtractor catalogs
    to the FLT frame, assuming that the FLT headers have been TweakReg'ed 
    to the `ref_image` alignment image.
    
    pointing='cosmos-19-F140W'; pad=60; NGROW=125; growx=2; growy=2; auto_offsets=False;  ref_image='Catalog/cosmos_3dhst.v4.0.IR_orig_sci.fits'; ref_ext=0; seg_image='Catalog/cosmos_3dhst.v4.0.F160W_seg.fits'; cat_file='Catalog/cosmos_3dhst.v4.0.IR_orig.cat'; ref_filter='F140W'
    
    """
    
    import astropy.io.fits as pyfits
    import astropy.units as u
    
    from astropy.table import Table as table
    
    import drizzlepac
    import stwcs
    from drizzlepac import astrodrizzle, tweakreg, tweakback
    
    import threedhst
        
    pointing_root = '-'.join(pointing.split('-')[:-1])
    pointing_filt = pointing.split('-')[-1]
    #### Open files and get WCS
    asn = threedhst.utils.ASNFile('%s_asn.fits' %(pointing))
    
    print 'Read files...'
    ref = pyfits.open(ref_image)
    seg = pyfits.open(seg_image)
    seg_data = np.cast[np.float32](seg[0].data)
    seg_ones = np.cast[np.float32](seg_data > 0)
    
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=ref_ext)
    seg_wcs = stwcs.wcsutil.HSTWCS(seg, ext=0)
    
    if ACS:
        grism = 'G800L'
        fl_ext = 'flc'
        sci_ext = [4,1]
    else:
        fl_ext = 'flt'
        sci_ext = [1]
        
    #### Only use single reference ACS exposure since growx=growy=1
    ref_flt_exp = asn.exposures[ref_exp]
    if ACS:
        asn.exposures = [asn.exposures[ref_exp]]
        
    #### Loop through FLTs, blotting reference and segmentation
    for exp in asn.exposures:
        threedhst.showMessage('Blot from reference:\n%s -> %s_blot.fits\n%s -> %s_seg.fits' %(ref_image, exp, seg_image, exp))
        flt = pyfits.open('%s_%s.fits' %(exp, fl_ext))
        #
        print 'Blot image+segmentation: %s_blot/seg.fits' %(exp)
        if ACS:
            hdu_blot = [pyfits.PrimaryHDU()]
            hdu_seg = [pyfits.PrimaryHDU()]
        else:
            hdu_blot = []
            hdu_seg = []
            
        for ext in sci_ext:
            flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=ext)
            #
            flt_wcs.naxis1 = flt[ext].header['NAXIS1']+NGROWX*2
            flt_wcs.naxis2 = flt[ext].header['NAXIS2']+NGROWY*2
            flt_wcs.wcs.crpix[0] += NGROWX
            flt_wcs.sip.crpix[0] += NGROWX
            flt_wcs.wcs.crpix[1] += NGROWY
            flt_wcs.sip.crpix[1] += NGROWY           
            ### Put WCS info back in a header
            header = flt_wcs.wcs2header(sip2hdr=True)
            header['ORIGFILE'] = (ref_image, 'Reference image path')
            header['EXPTIME'] = ref[ref_ext].header['EXPTIME']
            header['FILTER'] = ref_filter
            header['PAD'] = (pad, 'Padding pixels to x and y')
            header['NGROWX'] = (NGROWX, 'Number of pixels added to X-axis (both sides)')
            header['NGROWY'] = (NGROWY, 'Number of pixels added to Y-axis (both sides)')
            #
            ### reference
            blotted_ref = astrodrizzle.ablot.do_blot(ref[ref_ext].data, ref_wcs, flt_wcs, 1, coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)
            hdu_blot.append(pyfits.ImageHDU(data=blotted_ref, header=header.copy()))

            ### Blot segmentation, need "ones" image for pixel areas
            blotted_seg = astrodrizzle.ablot.do_blot(seg_data, seg_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
            blotted_ones = astrodrizzle.ablot.do_blot(seg_ones, seg_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)

            blotted_ones[blotted_ones == 0] = 1
            ratio = np.round(blotted_seg/blotted_ones)
            grow = nd.maximum_filter(ratio, size=3, mode='constant', cval=0)
            ratio[ratio == 0] = grow[ratio == 0]
            header['ORIGFILE'] = (seg_image, 'Segmentation image path')
            hdu_seg.append(pyfits.ImageHDU(data=np.cast[np.int32](ratio), header=header.copy()))
        
        if ACS:
            hdulist_blot = pyfits.HDUList(hdu_blot)
            hdulist_blot.writeto('%s_blot.fits' %(exp), clobber=True)
            hdulist_seg = pyfits.HDUList(hdu_seg)
            hdulist_seg.writeto('%s_seg.fits' %(exp), clobber=True)
        else:
            pyfits.writeto('%s_blot.fits' %(exp), data=hdu_blot[0].data, header=hdu_blot[0].header, clobber=True)
            pyfits.writeto('%s_seg.fits' %(exp), data=hdu_seg[0].data, header=hdu_seg[0].header, clobber=True)

    #### Combine interlaced images, ignore NGROW
    if ACS:
        roots = [pointing_root+'-chip1', pointing_root+'-chip2']
        
        empty_sci = np.zeros((2048+pad, 4096+pad), dtype=np.float64)
        empty_seg = np.zeros((2048+pad, 4096+pad), dtype=np.int32)
        
        #### Make reference image/segimage for each of two ACS chips        
        for chip in [1,2]:
            hdu_blot = hdulist_blot[chip]
            hdu_seg = hdulist_seg[chip]
            
            #### Pixel area map to preserve photometry in FLT frame
            im = pyfits.open(os.getenv('jref')+'wfc%d_pam.fits' %(chip))
            PAM = im[0].data
            im.close()
            hdu_blot.data *= PAM
            
            for i in [1,2]:
                hdu_blot.header['CRPIX%d' %(i)] += pad/2
                hdu_seg.header['CRPIX%d' %(i)] += pad/2
                        
            empty_sci[pad/2:-pad/2, pad/2:-pad/2] = hdu_blot.data
            hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdu_blot.header), pyfits.ImageHDU(data=empty_sci, header=hdu_blot.header)])
            hdulist.writeto('%s-chip%d-%s_ref_inter.fits' %(pointing_root, chip, grism), clobber=True)

            empty_seg[pad/2:-pad/2, pad/2:-pad/2] = hdu_seg.data
            pyfits.writeto('%s-chip%d-%s_inter_seg.fits' %(pointing_root, chip, grism), data=empty_seg, header=hdu_seg.header, clobber=True)
        
        ### xx break out
        #return True
        dxs, dys = np.zeros(100), np.zeros(100)
        
    else:
        roots = [pointing_root]
        dxs, dys = unicorn.reduce.interlace_combine_blot(root=pointing, view=False, pad=pad, CATALOG='UCSC/catalogs/COSMOS_F160W_v1.cat',  NGROWX=NGROWX, NGROWY=NGROWY, verbose=True, growx=growx, growy=growy, auto_offsets=auto_offsets, NSEGPIX=8, stop_early=True, ref_exp=ref_exp, grism=grism, reference_pixel=reference_pixel)
        
    ### Make the pointing catalog
    for i, pointing_root in enumerate(roots):
        threedhst.showMessage('Make catalog: %s-%s_inter.cat' %(pointing_root, grism))
        #### Create blotted catalog of objects within blotted segm. image
        seg_inter = pyfits.open('%s-%s_inter_seg.fits' %(pointing_root, grism))
        objects = np.unique(seg_inter[0].data)[1:]
    
        cat = table.read(cat_file, format='ascii.sextractor')
        NOBJ = len(cat)
        keep = np.zeros(NOBJ) > 1
        for j in range(NOBJ):
            if cat['NUMBER'][j] in objects:
                keep[j] = True
        #
        sub_cat = cat[keep]
        np.savetxt('%s-%s_radec.dat' %(pointing_root, grism), np.array([sub_cat['X_WORLD'], sub_cat['Y_WORLD']]).T, fmt='%.7f')
     
        threedhst.showMessage('%s: %d objects from blotted catalog:\n  %s.' %(pointing_root, keep.sum(), cat_file))
    
        ### here's the line to translate the x,y coords to the flt distorted frame
        x_flt0, y_flt0 = drizzlepac.skytopix.rd2xy('%s_%s.fits[%d]' %(ref_flt_exp, fl_ext, sci_ext[i]), coordfile='%s-%s_radec.dat' %(pointing_root, grism), verbose=False)
        x_flt = (x_flt0+dxs[ref_exp]+NGROWX)*growx+pad/2#-1 
        y_flt = (y_flt0+dys[ref_exp]+NGROWY)*growy+pad/2#-1 
    
        if growx == 2:
            x_flt -= 1
        if growy == 2:
            y_flt -= 1
        
        #### Write catalog
        x_flt_column = table.Column(data=x_flt, name='X_FLT', description='Interlaced pixel coordinate, x', unit=u.Unit('pix'))
        y_flt_column = table.Column(data=y_flt, name='Y_FLT', description='Interlaced pixel coordinate, y', unit=u.Unit('pix'))
        sub_cat.add_column(x_flt_column)
        sub_cat.add_column(y_flt_column)
        sub_cat.write('%s-%s_inter.cat' %(pointing_root, grism), format='ascii.commented_header')
    
        cols = sub_cat.colnames
        cat_header = []
        for j in range(len(cols)):
            c = sub_cat[cols[j]]
            line = '# %3d %-24s %-59s' %(j+1, c.name, c.description)
            # if c.unit is not None:
            #     line += ' [%s]' %(c.unit.__str__())
            # #
            cat_header.append(line+'\n')
    
        lines = open('%s-%s_inter.cat' %(pointing_root, grism)).readlines()
        fp = open('%s-%s_inter.cat' %(pointing_root, grism), 'w')
        fp.writelines(cat_header)
        fp.writelines(lines[1:])
        fp.close()
    
        #### Make region file
        lines = ['global color=green dashlist=8 3 width=1 font="helvetica 8 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n', 'image\n']
        lines.extend(['circle(%.2f, %.2f, 5.0) # text={%d}\n' %(x_flt[j], y_flt[j], sub_cat['NUMBER'][j]) for j in range(keep.sum())])
        fp = open('%s-%s_inter.reg' %(pointing_root, grism),'w')
        fp.writelines(lines)
        fp.close()
    
    if not ACS:
        #### Diagnostic plot
        threedhst.showMessage('Make centering diagnostic plot: %s_XY_FLT.png\n\nunicorn.reduce.adjust_catalog_for_FLT(pointing=\'%s\', MAG_LIM=23.5, ref=True)' %(pointing_root, pointing_root))
        unicorn.reduce.adjust_catalog_for_FLT(pointing=pointing_root, MAG_LIM=23.5, ref=True, grism=grism)
            
def adjust_catalog_for_FLT(pointing='aegis-15', MAG_LIM=25, THUMB_SIZE=16, ref=True, filter='F140W', grism='G141'):
    """
    Check for residuals in the blotted catalog
    """ 
    from astropy.table import Table as table
    cat = table.read('%s-%s_inter.cat' %(pointing, grism), format='ascii.sextractor')
    if ref:
        im_ref = pyfits.open('%s-%s_ref_inter.fits' %(pointing, grism))
        ext='ref'
        qscl = 1./400
    else:
        im_ref = pyfits.open('%s-%s_inter.fits' %(pointing, filter))
        ext=filter
        qscl = 1./400
        
    im_seg = pyfits.open('%s-%s_inter_seg.fits' %(pointing, grism))
    #im = pyfits.open('aegis-15-F140W_inter.fits')
    
    ok = (cat['MAG_AUTO'] < MAG_LIM) & (cat['X_FLT'] > 2*THUMB_SIZE) & (cat['X_FLT'] < (im_ref[1].header['NAXIS1']-2*THUMB_SIZE)) & (cat['Y_FLT'] > 2*THUMB_SIZE) & (cat['Y_FLT'] < (im_ref[1].header['NAXIS1']-2*THUMB_SIZE))
    
    idx = np.arange(len(cat))[ok]
    dx = cat['X_FLT']*1.
    dy = cat['X_FLT']*1.
    
    #fp = open('/tmp/xx.reg','w')
    #fp.write('image\n')
    #THUMB_SIZE = 16
    iy, ix = np.indices((2*THUMB_SIZE, 2*THUMB_SIZE))

    for i in idx:
        xc, yc = int(np.round(cat['X_FLT'][i])), int(np.round(cat['Y_FLT'][i]))
        ddx, ddy = cat['X_FLT'][i]-xc, cat['Y_FLT'][i]-yc
        #
        xi, yi = ((cat['X_FLT'][i])), ((cat['Y_FLT'][i]))
        if ((yc-THUMB_SIZE >= 0) & (yc+THUMB_SIZE <= im_ref[1].data.shape[0]) & (xc-THUMB_SIZE >= 0) & (xc+THUMB_SIZE <= im_ref[1].data.shape[1])):
            thumb_ref = im_ref[1].data[yc-THUMB_SIZE:yc+THUMB_SIZE, xc-THUMB_SIZE:xc+THUMB_SIZE]
            thumb_seg = im_seg[0].data[yc-THUMB_SIZE:yc+THUMB_SIZE, xc-THUMB_SIZE:xc+THUMB_SIZE] == cat['NUMBER'][i]
            thumb_ref /= np.sum(thumb_ref[thumb_seg])
            thumb_ref[~thumb_seg] = 0
            dx[i] = np.sum(thumb_ref*(ix-(THUMB_SIZE-1)))-ddx
            dy[i] = np.sum(thumb_ref*(iy-(THUMB_SIZE-1)))-ddy
        
    #fp.close()
    
    fig = unicorn.plotting.plot_init(xs=10, aspect=1./3, square=True, left=0.1, bottom=0.08, top=0.05, NO_GUI=True)
    
    ### Histogram
    ax = fig.add_subplot(131)
    ax.hist(dx[ok], bins=50, range=(-2,2), alpha=0.5, label=r'$\Delta x$, %5.2f' %(np.median(dx[ok])), normed=True)
    ax.hist(dy[ok], bins=50, range=(-2,2), alpha=0.5, label=r'$\Delta y$, %5.2f' %(np.median(dy[ok])), normed=True)
    ax.legend(loc='upper right', prop={'size':8}, title=pointing)
    
    ### Offset quiver    
    ax = fig.add_subplot(132)
    ax.quiver(cat['X_FLT'][ok], cat['Y_FLT'][ok], np.clip(dx[ok], -3, 3), np.clip(dy[ok], -3, 3), scale=qscl, scale_units='xy', angles='xy', alpha=0.9)
    
    ax.quiver(1294, 2588-100, 1., 0, scale=qscl, scale_units='xy', angles='xy', alpha=1, color='red')
    ax.text(1294, 2508, r'$\Delta$ = 1 (interlaced) pixel', color='red', ha='center', va='bottom')
    ax.set_xlim(-200,2788); ax.set_ylim(-200,2788)
    ax.set_title('Offset')
    
    ### Residual quiver
    ax = fig.add_subplot(133)
    ax.quiver(cat['X_FLT'][ok], cat['Y_FLT'][ok], np.clip(dx[ok]-np.median(dx[ok]), -1, 1), np.clip(dy[ok]-np.median(dy[ok]), -1, 1), scale=qscl, scale_units='xy', angles='xy', alpha=0.9)
    
    ax.quiver(1294, 2588-100, 1., 0, scale=qscl, scale_units='xy', angles='xy', alpha=1, color='red')
    ax.text(1294, 2508, r'$\Delta$ = 1 (interlaced) pixel', color='red', ha='center', va='bottom')
    ax.set_xlim(-200,2788); ax.set_ylim(-200,2788)
    ax.set_title('Residual')
    
    unicorn.plotting.savefig(fig, '%s_%s_XY_FLT.png' %(pointing, ext))

def blot_from_reference(REF_ROOT = 'COSMOS_F160W', DRZ_ROOT = 'COSMOS-19-F140W', NGROW=125, verbose=True, run_multidrizzle=False, realign_blotted=True):
    """
    iraf.imcopy('COSMOS-full-F160W_drz_sci.fits[2262:6160,5497:8871]', 'COSMOS-full-F160W_drz_sci_subim.fits')
    hedit 'COSMOS-full-F160W_drz_sci_subim.fits' 'CD1_2' 0. add+ update+ verify-
    hedit 'COSMOS-full-F160W_drz_sci_subim.fits' 'CD2_1' 0. add+ update+ verify-

    iraf.imcopy('COSMOS-full-F160W_drz_weight.fits[2262:6160,5497:8871]', 'COSMOS-full-F160W_drz_wht_subim.fits')
    hedit 'COSMOS-full-F160W_drz_wht_subim.fits' 'CD1_2' 0. add+ update+ verify-
    hedit 'COSMOS-full-F160W_drz_wht_subim.fits' 'CD2_1' 0. add+ update+ verify-

    """
    from pyraf import iraf
    from iraf import dither
    import scipy.ndimage as nd
    
    import unicorn.reduce
    
    tmpname = threedhst.utils.gen_tempname()
            
    im_ref = pyfits.open(REF_ROOT+'_ref.fits')
    
    ref_pixscale = np.sqrt(im_ref[0].header['CD1_1']**2+im_ref[0].header['CD1_2']**2)*3600
    print '\n\nReference pixscale: %.3f\n\n' %(ref_pixscale)
    
    if not os.path.exists(DRZ_ROOT+'.run'):
        run_multidrizzle=True
    else:
        run = threedhst.prep_flt_files.MultidrizzleRun(DRZ_ROOT)
        #### Need multidrizzle ".run" file to tell how to blot the image to the reference
        #### The line below makes sure the ".run" parameters match the pixel grid of the reference image
        ####    !!! hard coded for WFC3/IR
        if ((run.outnx, run.outny) != (im_ref[0].header['NAXIS1'], im_ref[0].header['NAXIS2'])) | (int(ref_pixscale/float(run.scl)*1.e5) != 12825):
            run_multidrizzle=True
    
    if run_multidrizzle:
        if verbose:
            threedhst.showMessage('Get exact shifts needed to match %s to \n%s' %(DRZ_ROOT, REF_ROOT))
        #"refimage" wasn't giving me *exactly* the same image dimensions
        #threedhst.prep_flt_files.startMultidrizzle(DRZ_ROOT+'_asn.fits', use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False, refimage=REFERENCE, build_drz=False)
        if 'CD1_2' not in im_ref[0].header.keys():
            im_ref[0].header.update('CD1_2',0.)
        
        # ref_wcs = wcsutil.WCSObject('', header=im_ref[0].header)
        # ra0, dec0 = ref_wcs.xy2rd((im_ref[0].header['NAXIS1']/2., im_ref[0].header['NAXIS2']/2.))
        
        ref_wcs = pywcs.WCS(header=im_ref[0].header)
        try:
            ra0, dec0 = ref_wcs.wcs_pix2world([[im_ref[0].header['NAXIS1']/2., im_ref[0].header['NAXIS2']/2.]],1)[0]
        except:
            ra0, dec0 = ref_wcs.wcs_pix2sky([[im_ref[0].header['NAXIS1']/2., im_ref[0].header['NAXIS2']/2.]],1)[0]
            
        threedhst.prep_flt_files.startMultidrizzle(DRZ_ROOT+'_asn.fits', use_shiftfile=True, skysub=False, final_scale=ref_pixscale, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False, build_drz=False, ra=ra0, dec=dec0, final_outnx=im_ref[0].header['NAXIS1'], final_outny=im_ref[0].header['NAXIS2'], final_rot=-np.arctan(im_ref[0].header['CD1_2']/im_ref[0].header['CD1_1'])/2/np.pi*360, generate_run=True)
        
    if verbose:
        threedhst.showMessage('Prepare images for `BLOT`')
    
    run = threedhst.prep_flt_files.MultidrizzleRun(DRZ_ROOT)
    #im_drz = pyfits.open(DRZ_ROOT+'_drz_sci.fits')
    ## DRZ images probably not needed and can be deleted.  If so, they'll be generated
    ## each time "blot_from_reference" is run.
    #os.remove(DRZ_ROOT+'_drz_sci.fits')
    #os.remove(DRZ_ROOT+'_drz_weight.fits')
    
    threedhst.process_grism.flprMulti()
    
    #### Force reference image to have same header as newly-created DRZ image so that
    #### BLOT has all of the necessary keywords
    
    #im_wht = pyfits.open(WEIGHT)
    
    # ### xxx compare headers
    # for key in im_ref[0].header.keys():
    #     if key not in im_wht[0].header.keys():
    #         print key, im_ref[0].header[key]
            
    #### Segmentation image needs to be float, not int, for blot, also need binary segmentation 
    #### image with 1. for objects, 0 for nothing
    
    REF_FILTER = im_ref[0].header['FILTER']
    REF_EXPTIME = im_ref[0].header['EXPTIME']
    
    #im_ref[0].header = im_drz[0].header.copy()
    #im_ref.writeto('ref_sci.fits', clobber=True)
    #DATA = 'ref_sci.fits'
    
    #im_wht[0].header = im_drz[0].header.copy()
    #im_wht.writeto('ref_wht.fits', clobber=True)
    
    asn = threedhst.utils.ASNFile(DRZ_ROOT+'_asn.fits')
    for iflt, flt_root in enumerate(asn.exposures):
        
        if verbose:
            threedhst.showMessage('Blot sci/weight/seg images for %s' %(flt_root))
        
        #FLT = 'ibhm51x0q_flt.fits'
        FLT = flt_root+'_flt.fits'        
        im_flt = pyfits.open(FLT)

        #### Need to modify the coeffs file 
        fp = open(FLT.replace('.fits','_coeffs1.dat'))
        coeffs_lines = fp.readlines()
        fp.close()

        if NGROW > 0:
            for i, line in enumerate(coeffs_lines):
                if line.strip().startswith('refpix'):
                    ### Default to center pixel
                    coeffs_lines[i] = 'refpix %9.3f 507\n' %((NGROW*2+1014)*1./2)

        fp = open(tmpname+'_coeffs.dat','w')
        fp.writelines(coeffs_lines)
        fp.close()
        
        #### We've read in the input shifts from the Multidrizzle ".run" file.  Find 
        #### the index of the current FLT image in the .run arrays
        idx = -1
        for i in range(run.count):
            if run.flt[i].startswith(FLT.replace('.fits','')):
                idx = i
                break
        
        ##### Re-alignment of blotted images to the FLT frame.  The initial
        ##### shifts determined from the "prep" stage should have put everything
        ##### onto the same coordinate system, but there were apparent offsets between
        ##### the blotted and F140W interlaced images.
        #####
        ##### The blotted images **need** to be on the same frame as the FLT images because
        ##### they determine the reference position of the grism spectra
        if iflt == 0:
          if realign_blotted:
            fp = open(DRZ_ROOT+'_delta.dat','w')
            fp.write('# xoff  yoff  roff dx  dy  xrms  yrms\n')
            #
            xoff, yoff, roff = 0., 0., 0.
            #for i in range(3):
            for geom in ['rotate', 'shift', 'rotate', 'shift', 'shift', 'shift']:
              #
                threedhst.process_grism.flprMulti()
                if os.path.exists('align_blot.fits'):
                    os.remove('align_blot.fits')
                #
                iraf.wblot.setParam('data', REF_ROOT+'_ref.fits')
                iraf.wblot.setParam('outdata', 'align_blot.fits')
                iraf.wblot.setParam('outnx', 1014)
                iraf.wblot.setParam('outny', 1014)
                iraf.wblot.setParam('geomode', 'user')
                iraf.wblot.setParam('interpol', 'poly5')
                iraf.wblot.setParam('sinscl', 1.0)
                iraf.wblot.setParam('coeffs', run.flt[idx]+'_coeffs1.dat')
                iraf.wblot.setParam('lambda', 555.)
                iraf.wblot.setParam('xgeoim', '')
                iraf.wblot.setParam('ygeoim', '')
                iraf.wblot.setParam('align', 'center')
                iraf.wblot.setParam('scale', float(run.scl))
                iraf.wblot.setParam('xsh', run.xsh[idx]+xoff)
                iraf.wblot.setParam('ysh', run.ysh[idx]+yoff)
                iraf.wblot.setParam('rot', run.rot[idx]+roff)
                iraf.wblot.setParam('shft_un', 'input')
                iraf.wblot.setParam('shft_fr', 'input')
                iraf.wblot.setParam('refimage', '')
                iraf.wblot.setParam('outscl', 0.128254)
                iraf.wblot.setParam('raref', im_flt[1].header['CRVAL1'])
                iraf.wblot.setParam('decref', im_flt[1].header['CRVAL2'])
                iraf.wblot.setParam('xrefpix', im_flt[1].header['CRPIX1'])
                iraf.wblot.setParam('yrefpix', im_flt[1].header['CRPIX2'])
                iraf.wblot.setParam('orient', im_flt[1].header['ORIENTAT'])
                iraf.wblot.setParam('dr2gpar', '')
                iraf.wblot.setParam('expkey', 'exptime')
                iraf.wblot.setParam('expout', 'input')
                iraf.wblot.setParam('in_un', 'cps')
                iraf.wblot.setParam('out_un', 'cps')
                iraf.wblot.setParam('fillval', 0.0)
                iraf.wblot.setParam('mode', 'al')
                #
                iraf.wblot.setParam('mode', 'h')
                iraf.wblot()
                
                # status = iraf.wblot( data = REF_ROOT+'_ref.fits', outdata = 'align_blot.fits', 
                #    outnx = 1014, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
                #    coeffs = run.flt[idx]+'_coeffs1.dat', xgeoim = '', ygeoim = '', 
                #    align = 'center', scale = float(run.scl), xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
                #    rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
                #    refimage = '', outscl = 0.128254, raref = im_flt[1].header['CRVAL1'], 
                #    decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1'], 
                #    yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
                #    dr2gpar = '', expkey = 'exptime', expout = 'input', 
                #    in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=0)
                # print status
                
                # import drizzlepac
                # from drizzlepac import astrodrizzle
                # import stwcs
                # ref = pyfits.open(REF_ROOT+'_ref.fits')
                # ref_wcs = stwcs.wcsutil.HSTWCS(ref)
                # flt_wcs = stwcs.wcsutil.HSTWCS(im_flt, ext=1)
                # ref_wcs.wcs.crpix += np.array([xoff, yoff])
                # ref_wcs.rotateCD(roff/360.*2*np.pi)
                # 
                # blotted = astrodrizzle.ablot.do_blot(ref[0].data, ref_wcs, flt_wcs, im_flt[0].header['EXPTIME'], coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)/im_flt[0].header['EXPTIME']
                # pyfits.writeto('align_blot.fits', data=blotted, header=im_flt[1].header)
                
                #
                dx, dy, rot, xrms, yrms = unicorn.reduce.realign_blotted(flt=FLT, blotted='align_blot.fits', fitgeometry=geom)
                #### Transform the shifts/rotation in the FLT frame to the reference image frame
                alpha = np.pi+(run.rot[idx]+roff)/360*2*np.pi
                dx_flt = np.cos(alpha)*dx-np.sin(alpha)*dy
                dy_flt = np.sin(alpha)*dx+np.cos(alpha)*dy
                #
                if geom == 'shift':
                    xoff += dx_flt/np.float(run.scl)
                    yoff += dy_flt/np.float(run.scl)
                #
                roff += rot
                fp.write('%.3f %.3f %.4f  %.3f %.3f  %.2f %.2f\n' %(xoff, yoff, roff, dx, dy, xrms, yrms))
                
                #threedhst.sex.sexcatRegions('direct.cat','direct.reg')
                #threedhst.sex.sexcatRegions('align.cat','align.reg')
                
            fp.close()
          else:
             fp = open(DRZ_ROOT+'_delta.dat','w')
             fp.write('# xoff  yoff  roff dx  dy  xrms  yrms\n')
             xoff, yoff, roff, dx, dy, xrms, yrms = 0., 0., 0., 0., 0., 0., 0.
             fp.write('%.3f %.3f %.4f  %.3f %.3f  %.2f %.2f\n' %(xoff, yoff, roff, dx, dy, xrms, yrms))
             fp.close()
             
        ##################################
        #### Run BLOT for SCI extension
        ##################################
        if verbose:
            print unicorn.noNewLine+'sci'
            
        if os.path.exists(FLT.replace('_flt','_blot')):
            os.remove(FLT.replace('_flt','_blot'))
        
        threedhst.process_grism.flprMulti()
        # status = iraf.wblot( data = REF_ROOT+'_ref.fits', outdata = FLT.replace('_flt','_blot'), 
        #    outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
        #    coeffs = tmpname+'_coeffs.dat', xgeoim = '', ygeoim = '', 
        #    align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
        #    rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
        #    refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
        #    decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
        #    yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
        #    dr2gpar = '', expkey = 'exptime', expout = 'input', 
        #    in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)

        iraf.wblot.setParam('data', REF_ROOT+'_ref.fits')
        iraf.wblot.setParam('outdata', FLT.replace('_flt','_blot'))
        iraf.wblot.setParam('outnx', 1014+2*NGROW)
        iraf.wblot.setParam('outny', 1014)
        iraf.wblot.setParam('geomode', 'user')
        iraf.wblot.setParam('interpol', 'poly5')
        iraf.wblot.setParam('sinscl', 1.0)
        iraf.wblot.setParam('coeffs', tmpname+'_coeffs.dat')
        iraf.wblot.setParam('lambda', 555.)
        iraf.wblot.setParam('xgeoim', '')
        iraf.wblot.setParam('ygeoim', '')
        iraf.wblot.setParam('align', 'center')
        iraf.wblot.setParam('scale', float(run.scl))
        iraf.wblot.setParam('xsh', run.xsh[idx]+xoff)
        iraf.wblot.setParam('ysh', run.ysh[idx]+yoff)
        iraf.wblot.setParam('rot', run.rot[idx]+roff)
        iraf.wblot.setParam('shft_un', 'input')
        iraf.wblot.setParam('shft_fr', 'input')
        iraf.wblot.setParam('refimage', '')
        iraf.wblot.setParam('outscl', 0.128254)
        iraf.wblot.setParam('raref', im_flt[1].header['CRVAL1'])
        iraf.wblot.setParam('decref', im_flt[1].header['CRVAL2'])
        iraf.wblot.setParam('xrefpix', im_flt[1].header['CRPIX1']+NGROW)
        iraf.wblot.setParam('yrefpix', im_flt[1].header['CRPIX2'])
        iraf.wblot.setParam('orient', im_flt[1].header['ORIENTAT'])
        iraf.wblot.setParam('dr2gpar', '')
        iraf.wblot.setParam('expkey', 'exptime')
        iraf.wblot.setParam('expout', 'input')
        iraf.wblot.setParam('in_un', 'cps')
        iraf.wblot.setParam('out_un', 'cps')
        iraf.wblot.setParam('fillval', 0.0)
        iraf.wblot.setParam('mode', 'al')

        iraf.wblot.setParam('mode', 'h')
        iraf.wblot()


        # import drizzlepac
        # from drizzlepac import astrodrizzle
        # import stwcs
        # ref = pyfits.open(REF_ROOT+'_ref.fits')
        # ref_wcs = stwcs.wcsutil.HSTWCS(ref)
        # flt_wcs = stwcs.wcsutil.HSTWCS(im_flt, ext=1)
        # ref_wcs.wcs.crpix += np.array([xoff, yoff])
        # ref_wcs.rotateCD(roff/360.*2*np.pi)
        # 
        # blotted = astrodrizzle.ablot.do_blot(ref[0].data, ref_wcs, flt_wcs, im_flt[0].header['EXPTIME'], coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)/im_flt[0].header['EXPTIME']
        # pyfits.writeto(FLT.replace('_flt','_blot'), data=blotted, header=im_flt[1].header)
                
        #### Add NGROW to header
        im_blt = pyfits.open(FLT.replace('_flt','_blot'), mode='update')
        im_blt[0].header.update('NGROW',NGROW, comment='Number of pixels added to X-axis (centered)')
        im_blt[0].header.update('FILTER', REF_FILTER)
        im_blt[0].header.update('EXPTIME', REF_EXPTIME)
        im_blt.flush()
        
        ##################################
        #### Blot for WHT extension
        ##################################
        # if verbose:
        #     print unicorn.noNewLine+'sci, weight'
        # 
        # if os.path.exists(FLT.replace('_flt','_blot_wht')):
        #     os.remove(FLT.replace('_flt','_blot_wht'))
        # 
        # threedhst.process_grism.flprMulti()
        # status = iraf.wblot( data = WEIGHT, outdata = FLT.replace('_flt','_blot_wht'), 
        #    outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
        #    coeffs = tmpname+'_coeffs.dat', lambd = 1392.0, xgeoim = '', ygeoim = '', 
        #    align = 'center', scale = run.scl, xsh = run.xsh[idx], ysh = run.ysh[idx], 
        #    rot = run.rot[idx], shft_un = 'input', shft_fr = 'input', 
        #    refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
        #    decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
        #    yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
        #    dr2gpar = '', expkey = 'exptime', expout = 'input', 
        #    in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)
        # 
        # #### Add NGROW to header
        # im_blt = pyfits.open(FLT.replace('_flt','_blot_wht'), mode='update')
        # 
        # #### Weight is inverse variance.  Interlace assumes sigma
        # im_blt[0].data = 1/np.sqrt(im_blt[0].data)
        # im_blt[0].header.update('NGROW',NGROW, comment='Number of pixels added to X-axis (centered)')
        # im_blt[0].header.update('FILTER', REF_FILTER)
        # im_blt[0].header.update('EXPTIME', REF_EXPTIME)
        # im_blt.flush()
        
        #################################
        #### Segmentation image
        #################################
        if verbose:
            print unicorn.noNewLine+'sci, weight, seg0'

        if os.path.exists(FLT.replace('_flt','_seg')):
            os.remove(FLT.replace('_flt','_seg'))
        
        # threedhst.process_grism.flprMulti()
        # status = iraf.wblot( data = REF_ROOT+'_seg.fits', outdata = FLT.replace('_flt','_seg'), 
        #    outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'nearest', sinscl = 1.0, 
        #    coeffs = tmpname+'_coeffs.dat', lambd = 1392.0, xgeoim = '', ygeoim = '', 
        #    align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
        #    rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
        #    refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
        #    decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
        #    yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
        #    dr2gpar = '', expkey = 'exptime', expout = 'input', 
        #    in_un = 'counts', out_un = 'counts', fillval = 0.0, mode = 'al', Stdout=1)
        
        iraf.wblot.setParam('data', REF_ROOT+'_seg.fits')
        iraf.wblot.setParam('outdata', FLT.replace('_flt','_seg'))
        iraf.wblot.setParam('outnx', 1014+2*NGROW)
        iraf.wblot.setParam('outny', 1014)
        iraf.wblot.setParam('geomode', 'user')
        iraf.wblot.setParam('interpol', 'nearest')
        iraf.wblot.setParam('sinscl', 1.0)
        iraf.wblot.setParam('coeffs', tmpname+'_coeffs.dat')
        iraf.wblot.setParam('lambda', 1392.)
        iraf.wblot.setParam('xgeoim', '')
        iraf.wblot.setParam('ygeoim', '')
        iraf.wblot.setParam('align', 'center')
        iraf.wblot.setParam('scale', float(run.scl))
        iraf.wblot.setParam('xsh', run.xsh[idx]+xoff)
        iraf.wblot.setParam('ysh', run.ysh[idx]+yoff)
        iraf.wblot.setParam('rot', run.rot[idx]+roff)
        iraf.wblot.setParam('shft_un', 'input')
        iraf.wblot.setParam('shft_fr', 'input')
        iraf.wblot.setParam('refimage', '')
        iraf.wblot.setParam('outscl', 0.128254)
        iraf.wblot.setParam('raref', im_flt[1].header['CRVAL1'])
        iraf.wblot.setParam('decref', im_flt[1].header['CRVAL2'])
        iraf.wblot.setParam('xrefpix', im_flt[1].header['CRPIX1']+NGROW)
        iraf.wblot.setParam('yrefpix', im_flt[1].header['CRPIX2'])
        iraf.wblot.setParam('orient', im_flt[1].header['ORIENTAT'])
        iraf.wblot.setParam('dr2gpar', '')
        iraf.wblot.setParam('expkey', 'exptime')
        iraf.wblot.setParam('expout', 'input')
        iraf.wblot.setParam('in_un', 'cps')
        iraf.wblot.setParam('out_un', 'cps')
        iraf.wblot.setParam('fillval', 0.0)
        iraf.wblot.setParam('mode', 'al')

        iraf.wblot.setParam('mode', 'h')
        iraf.wblot()
        
        #print status
        
        # seg = pyfits.open(REF_ROOT+'_seg.fits')
        # blotted = astrodrizzle.ablot.do_blot(seg[0].data, ref_wcs, flt_wcs, im_flt[0].header['EXPTIME'], coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)/im_flt[0].header['EXPTIME']
        # pyfits.writeto(FLT.replace('_flt','_seg'), data=blotted, header=im_flt[1].header)
        
        #
        if verbose:
            print unicorn.noNewLine+'sci, weight, seg0, seg1'
        
        if os.path.exists(FLT.replace('_flt','_ones')):
            os.remove(FLT.replace('_flt','_ones'))
        
        threedhst.process_grism.flprMulti()
        # status = iraf.wblot( data = REF_ROOT+'_ones.fits', outdata = FLT.replace('_flt','_ones'), 
        #    outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'nearest', sinscl = 1.0, 
        #    coeffs = tmpname+'_coeffs.dat', lambd = 1392.0, xgeoim = '', ygeoim = '', 
        #    align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
        #    rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
        #    refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
        #    decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
        #    yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
        #    dr2gpar = '', expkey = 'exptime', expout = 'input', 
        #    in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)

        iraf.wblot.setParam('data', REF_ROOT+'_ones.fits')
        iraf.wblot.setParam('outdata', FLT.replace('_flt','_ones'))
        iraf.wblot.setParam('outnx', 1014+2*NGROW)
        iraf.wblot.setParam('outny', 1014)
        iraf.wblot.setParam('geomode', 'user')
        iraf.wblot.setParam('interpol', 'nearest')
        iraf.wblot.setParam('sinscl', 1.0)
        iraf.wblot.setParam('coeffs', tmpname+'_coeffs.dat')
        iraf.wblot.setParam('lambda', 1392.)
        iraf.wblot.setParam('xgeoim', '')
        iraf.wblot.setParam('ygeoim', '')
        iraf.wblot.setParam('align', 'center')
        iraf.wblot.setParam('scale', float(run.scl))
        iraf.wblot.setParam('xsh', run.xsh[idx]+xoff)
        iraf.wblot.setParam('ysh', run.ysh[idx]+yoff)
        iraf.wblot.setParam('rot', run.rot[idx]+roff)
        iraf.wblot.setParam('shft_un', 'input')
        iraf.wblot.setParam('shft_fr', 'input')
        iraf.wblot.setParam('refimage', '')
        iraf.wblot.setParam('outscl', 0.128254)
        iraf.wblot.setParam('raref', im_flt[1].header['CRVAL1'])
        iraf.wblot.setParam('decref', im_flt[1].header['CRVAL2'])
        iraf.wblot.setParam('xrefpix', im_flt[1].header['CRPIX1']+NGROW)
        iraf.wblot.setParam('yrefpix', im_flt[1].header['CRPIX2'])
        iraf.wblot.setParam('orient', im_flt[1].header['ORIENTAT'])
        iraf.wblot.setParam('dr2gpar', '')
        iraf.wblot.setParam('expkey', 'exptime')
        iraf.wblot.setParam('expout', 'input')
        iraf.wblot.setParam('in_un', 'cps')
        iraf.wblot.setParam('out_un', 'cps')
        iraf.wblot.setParam('fillval', 0.0)
        iraf.wblot.setParam('mode', 'al')

        iraf.wblot.setParam('mode', 'h')
        iraf.wblot()
        
        # seg = pyfits.open(REF_ROOT+'_ones.fits')
        # blotted = astrodrizzle.ablot.do_blot(seg[0].data, ref_wcs, flt_wcs, im_flt[0].header['EXPTIME'], coeffs=True, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)/im_flt[0].header['EXPTIME']
        # pyfits.writeto(FLT.replace('_flt','_ones'), data=blotted, header=im_flt[1].header)
        
        #### Add NGROW to header
        im_seg = pyfits.open(FLT.replace('_flt','_seg'), mode='update')
        ones = pyfits.open(FLT.replace('_flt','_ones'))
        # yh, xh = np.histogram(ones[0].data.flatten(), range=(0.1,ones[0].data.max()), bins=100)
        # keep = ones[0].data > (0.5*xh[:-1][yh == yh.max()])
        # 
        ones[0].data[ones[0].data == 0] = 1
        ratio = im_seg[0].data / ones[0].data
        #test = (np.abs(np.log(ratio/np.round(ratio))) < 1.e-5) & keep
        
        # test = (np.abs(ratio-np.round(ratio)) < 0.01) & keep
        # ratio[test] = np.round(ratio[test])
        # ratio[~test | (ones[0].data == 0)] = 0
        
        # im_seg[0].data[keep] /= ones[0].data[keep]
        # im_seg[0].data[~keep] = 0
        #im_seg[0].data = np.cast[np.int32](np.round(im_seg[0].data))
        seg = np.cast[np.int32](np.round(ratio))
        ### Grow regions a bit, but not where already defined at adjacent segments
        grow = nd.maximum_filter(seg, size=3, mode='constant', cval=0)
        grow_pix = (seg == 0) & (grow > 0)
        seg[grow_pix] = grow[grow_pix]
        im_seg[0].data = seg
        
         
        im_seg[0].header.update('NGROW',NGROW, comment='Number of pixels added to X-axis (centered)')
        im_seg[0].header.update('FILTER', REF_FILTER)
        im_seg[0].header.update('EXPTIME', REF_EXPTIME)
        
        if not pyfits.__version__.startswith('3'):
            pyfits.writeto(im_seg.filename(), data=im_seg[0].data, header=im_seg[0].header, clobber='True')
        else:
            im_seg.flush()
    
        # os.remove(FLT.replace('_flt','_ones'))
    
        #### 
        # import threedhst.dq
        # ds9 = threedhst.dq.myDS9()
        # 
        # ds9.frame(1)
        # im = pyfits.open(FLT)
        # ds9.view(im[1].data)
        # ds9.scale(-0.05, 2)
        # 
        # ds9.frame(2)
        # im = pyfits.open(FLT.replace('flt','blot'))
        # ds9.view(im[0].data[:,NGROW:-NGROW])
        # ds9.scale(-0.05, 2)
        # 
        # ds9.frame(3)
        # im = pyfits.open(FLT.replace('flt','seg'))
        # ds9.view(im[0].data[:,NGROW:-NGROW])
        # 
        # ds9.frame(4)
        # im = pyfits.open(FLT.replace('flt','blot_wht'))
        # ds9.view(im[0].data[:,NGROW:-NGROW])
        # ds9.scale(-0.05, 2)
        
        
    ##### Clean up
    tmp_files = glob.glob(tmpname+'*')
    print 'Clean up:\n'
    iline = ''
    for ifile in tmp_files:
        iline += ', '+ifile
        print unicorn.noNewLine+iline
        os.remove(ifile)
#
def fill_inter_zero(image='TILE41-132-F160W_inter.fits', fill_error=True):
    from scipy.signal import convolve2d

    im = pyfits.open(image, mode='update')
    
    HAS_ERR = True
    x0, x1 = 'SCI', 'ERR'
    minpix = 0
    
    if len(im) == 1:
        HAS_ERR = False
        x0, x1 = 0, 0
        minpix = 4
        
    if len(im) == 2:
        HAS_ERR = False
        x0, x1 = 1, 1
        
    kernel = np.ones((3,3))
    wht = im[x1].data != 0
    bad = (wht == 0)
    if wht.sum() < 10:
        return None
        
    npix = convolve2d(wht, kernel, boundary='fill', fillvalue=0, mode='same')
    sum = convolve2d(im[x0].data, kernel, boundary='fill', fillvalue=0, mode='same')
    
    fill_pix = bad & (npix > minpix)
    im[x0].data[fill_pix] = (sum/npix)[fill_pix]*1
    
    if HAS_ERR & fill_error:
        sumwht = convolve2d(im[x1].data**2, kernel, boundary='fill', fillvalue=0, mode='same')
        im[x1].data[fill_pix] = np.sqrt((sumwht/npix)[fill_pix]*1.)
    
    im.flush()
    
def strip_header(header=None, filter='F140W'):
    """
    BLOT tends to die with a segmentation fault for some reason if I don't
    strip out all but the WCS keywords from the image header.
    """
    new_header = pyfits.Header()
    #new_header.update('NAXIS',2)
    #new_header.update('SIMPLE',True)
    
    if 'EXPTIME' in header.keys():
        new_header.update('EXPTIME', header.get('EXPTIME'))
    else:
        new_header.update('EXPTIME', 1.)
    
    if 'FILTER' in header.keys():
        new_header.update('FILTER', header.get('FILTER'))
    else:
        new_header.update('FILTER', filter)
            
    new_header.update('CDELT1', header.get('CD1_1'))
    new_header.update('CDELT2', header.get('CD2_2'))
    new_header.update('LTV1', 0.)
    new_header.update('LTV2', 0.)

    copy_keys = ['WCSAXES', 'CTYPE1', 'CTYPE2','CRVAL1', 'CRVAL2', 'CRPIX1','CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'LTM1_1', 'LTM2_2', 'PA_APER', 'ORIENTAT', 'RA_APER', 'DEC_APER', 'VAFACTOR']
    
    for key in copy_keys:
        #print 'Key: [%s]' %(key)
        if key in header.keys():
            new_header.update(key,  header.get(key))
        else:
            new_header.update(key, 0.)
            
    return new_header

def realign_blotted(flt='ibhj34h6q_flt.fits', blotted='align_blot.fits', fitgeometry='shift'):
    """
    Blot wasn't getting the images perfectly aligned, so do a catalog matching alignment
    on them to get some small additional shifts.
    """
    import threedhst
    import numpy as np
    from pyraf import iraf

    import os
    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    verbose = False
    #fitgeometry='shift'
    
    #flt = 'ibhj34h6q_flt.fits'
    #blt = 'ibhj34h6q_blot.fits'
    toler = 3
    
    try:
        os.remove('SCI.fits')
        os.remove('WHT.fits')
    except:
        print 'No SCI/WHT'
        
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'NONE'
    se.options['WEIGHT_TYPE']     = 'MAP_RMS'
    se.options['WEIGHT_IMAGE']    = 'WHT.fits'
    se.options['FILTER']    = 'Y'
    ## Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '3' 
    se.options['ANALYSIS_THRESH']  = '3' 
    se.options['MAG_ZEROPOINT'] = '26.46'

    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['CATALOG_NAME']    = 'direct.cat'
    iraf.imcopy('%s[1]' %(flt),"SCI.fits")
    iraf.imcopy('%s[2]' %(flt),"WHT.fits")
    
    ### Clean up bad pixels and cosmic rays
    flt_im = pyfits.open(flt)
    bad = (flt_im['SCI'].data > 1.e5) | ((flt_im['DQ'].data & 4096) > 0)
    sci = pyfits.open('SCI.fits', mode='update')
    wht = pyfits.open('WHT.fits', mode='update')
    sci[0].data[bad] = 0
    wht[0].data[bad] = 1.e6
    sci.flush()
    wht.flush()
    
    status = se.sextractImage('SCI.fits')

    ## alignment image
    se.options['CATALOG_NAME']    = 'align.cat'
    status = se.sextractImage(blotted)

    ## Read the catalogs
    directCat = threedhst.sex.mySexCat('direct.cat')
    alignCat = threedhst.sex.mySexCat('align.cat')
        
    xshift = 0
    yshift = 0
    rot = 0
    scale = 1.
    
    xrms = 2
    yrms = 2
    
    NITER = 5
    IT = 0
    while (IT < NITER):
        IT = IT+1
        
        #### Get x,y coordinates of detected objects
        ## direct image
        fp = open('direct.xy','w')
        for i in range(len(directCat.X_IMAGE)):
            fp.write('%s  %s\n' %(directCat.X_IMAGE[i],directCat.Y_IMAGE[i]))
        fp.close()

        ## alignment image
        fp = open('align.xy','w')
        for i in range(len(alignCat.X_IMAGE)):
            fp.write('%s  %s\n' %(np.float(alignCat.X_IMAGE[i])+xshift,
                       np.float(alignCat.Y_IMAGE[i])+yshift))
        fp.close()

        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        #### iraf.xyxymatch to find matches between the two catalogs
        pow = toler*1.
        try:
            os.remove('align.match')
        except:
            pass
        
        status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                       output="align.match",
                       tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
        
        while status1[-1].startswith('0'):
            pow+=1
            os.remove('align.match')
            status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                           output="align.match",
                           tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
            
        if verbose:
            for line in status1:
                print line
                
        #### Compute shifts with iraf.geomap
        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        try:
            os.remove("align.map")
        except:
            pass
            
        status2 = iraf.geomap(input="align.match", database="align.map",
                    fitgeometry=fitgeometry, interactive=no, 
                    xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF,
                    maxiter = 10, reject = 2.0, Stdout=1)
        if verbose:
            for line in status2:
                print line
        
        #fp = open(root+'.iraf.log','a')
        #fp.writelines(status1)
        #fp.writelines(status2)
        #fp.close()
                
        #### Parse geomap.output 
        fp = open("align.map","r")
        for line in fp.readlines():
            spl = line.split()
            if spl[0].startswith('xshift'):
                xshift += float(spl[1])    
            if spl[0].startswith('yshift'):
                yshift += float(spl[1])    
            if spl[0].startswith('xrotation'):
                rot = float(spl[1])    
            if spl[0].startswith('xmag'):
                scale = float(spl[1])    
            if spl[0].startswith('xrms'):
                xrms = float(spl[1])    
            if spl[0].startswith('yrms'):
                yrms = float(spl[1])    
            
        fp.close()
        
        #os.system('wc align.match')
        print 'Shift iteration #%d, xshift=%f, yshift=%f, rot=%f, scl=%f (rms: %5.2f,%5.2f)' %(IT, xshift, yshift, rot, scale, xrms, yrms)
        
        #clean_files = ['align.cat', 'align.map','align.match', 'align.xy','direct.xy', 'sex_stderr', 'SCI.fits','WHT.fits','direct.cat','threedhst_auto.param','threedhst_auto.sex', 'align_blot.fits','default.conv','default.nnw']
        #for file in clean_files:
        #    os.remove(file)
            
    return xshift, yshift, rot, xrms, yrms
    
def plot_defaults(point=True):
    plt.rcParams['patch.edgecolor'] = 'None'
    plt.rcParams['font.size'] = 10

    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['image.interpolation'] = 'nearest'

    if point:
        plt.rcParams['lines.linestyle'] = 'None'
        plt.rcParams['lines.marker'] = 'o'
    else:
        plt.rcParams['lines.linestyle'] = '-'
        plt.rcParams['lines.marker'] = 'None'
        
def all_deep():
    files = glob.glob('*G141_inter.fits')
    for file in files: 
        root=file.split('-G141')[0]
        deep_model(root=root, MAG_LIMIT=28)

def check_kluge_shifts():
    files=glob.glob('AEG*G141_inter.fits')

    for file in files:
        pass
        #unicorn.reduce.model_kluge(file)
    #
    files=glob.glob('AEG*G141_inter.fits')
    for file in files:
        print file
        g141 = pyfits.open(file, mode='update')    
        model = pyfits.open(file.replace('-G141_inter', '_inter_model'))
        #pyfits.writeto(file.replace('inter','diff'), g141[1].data-model[0].data, header=g141[1].header, clobber=True)
        gshift = g141[1].data
        gshift[:,:-1] = gshift[:,1:]
        pyfits.writeto(file.replace('inter','diff_shift1'), gshift-model[0].data, header=g141[1].header, clobber=True)
        gshift[:,:-1] = gshift[:,1:]
        pyfits.writeto(file.replace('inter','diff_shift2'), gshift-model[0].data, header=g141[1].header, clobber=True)
        
def model_kluge(inter_grism='GOODS-S-34-G141_inter.fits', xshift=1):
    """
    Some indexing error somewhere results in redshifts that are too high by ~1 pixel.  Shift the 
    interlaced grism images "left" by this amount to take it out by hand....
    """
    
    import scipy.ndimage as nd
    import stsci.convolve
    import unicorn
    
    #inter_grism = 'COSMOS-10-G141_inter.fits'
    #inter_grism='GOODS-S-34-G141_inter.fits'
    print inter_grism
    g141 = pyfits.open(inter_grism, mode='update')
    unicorn.reduce.fill_inter_zero(inter_grism, fill_error=False)
    
    model = pyfits.open(inter_grism.replace('-G141_inter', '_inter_model'))
    
    ### Taper for cross correlation
    print unicorn.noNewLine+'%s  Taper' %(inter_grism)
    taper = nd.gaussian_filter((model[0].data > 0.03)*1., 8)
    NGROW = g141[1].header['NGROW']
    pad = 300
    g141_tapered = (g141[1].data*taper)[2*NGROW+pad*2:-2*NGROW-pad*2, 2*NGROW+pad*2:-2*NGROW-pad*2]
    model_tapered = (model[0].data*taper)[2*NGROW+pad*2:-2*NGROW-pad*2, 2*NGROW+pad*2:-2*NGROW-pad*2]
    
    print unicorn.noNewLine+'%s  Corr1' %(inter_grism)
    cross_corr = stsci.convolve.correlate2d(g141_tapered, model_tapered, fft=1)
    print unicorn.noNewLine+'%s  Corr2' %(inter_grism)
    cross_corr_ref = stsci.convolve.correlate2d(g141_tapered, g141_tapered, fft=1)
    sh = cross_corr.shape
    yi, xi = np.indices(cross_corr[sh[1]/2-20:sh[1]/2+20, sh[0]/2-20:sh[0]/2+20].shape)
    xc, yc = xi.flatten()[np.argmax(cross_corr[sh[1]/2-20:sh[1]/2+20, sh[0]/2-20:sh[0]/2+20])], yi.flatten()[np.argmax(cross_corr[sh[1]/2-20:sh[1]/2+20, sh[0]/2-20:sh[0]/2+20])]
    xc_ref, yc_ref = xi.flatten()[np.argmax(cross_corr_ref[sh[1]/2-20:sh[1]/2+20, sh[0]/2-20:sh[0]/2+20])], yi.flatten()[np.argmax(cross_corr_ref[sh[1]/2-20:sh[1]/2+20, sh[0]/2-20:sh[0]/2+20])]
    dx, dy = xc-xc_ref, yc-yc_ref
    
    #g141[1].data[:,:-4] = g141[1].data[:,4:]
    #g141[2].data[:,:-4] = g141[2].data[:,4:]
    
    print unicorn.noNewLine+'%s  %d  %d' %(inter_grism, dx, dy)
    
    #ds9.view(g141[1].data-model[0].data)
    
def check_stars(pointing='GOODS-S-34'):
    os.chdir('/Users/brammer/3DHST/Spectra/Work/GOODS-S/PREP-FLT')
    
    model = unicorn.reduce.GrismModel(pointing)
    
    plt.scatter(model.cat['MAG_AUTO'], model.cat['FLUX_RADIUS'], alpha=0.5)
    
    stars = (model.cat['FLUX_RADIUS'] < 4) & (model.cat['MAG_AUTO'] < 23)
    for id in model.cat['NUMBER'][stars]:
        model.twod_spectrum(id, miny=-23)
        model.show_2d(savePNG=True)
        
    for id in model.cat['NUMBER'][stars]:
        twod = unicorn.reduce.Interlace2D('%s_%05d.2D.fits' %(pointing, id))
        wave, flux = twod.optimal_extract(twod.im['SCI'].data-twod.im['CONTAM'].data)
        sens = twod.im['SENS'].data
        plt.plot(wave, np.roll(flux, 0)/sens, color='black')
        plt.plot(wave, np.roll(flux, 2)/sens, color='blue')
        plt.plot(wave, np.roll(flux, 4)/sens, color='red')
        plt.plot(wave, np.roll(flux, -2)/sens, color='green')
        
def deep_model(root='COSMOS-19', MAG_LIMIT=28):
    """ 
    Use the model generator to make a full model of all objects in a pointing down to 
    very faint limits to test the grism background subtraction.
    """
    import matplotlib 
    
    model = unicorn.reduce.GrismModel(root=root, grow_factor=2, growx=2, growy=2, MAG_LIMIT=MAG_LIMIT)
    
    if model.cat.mag.max() < (MAG_LIMIT-0.3):
        model.find_objects(MAG_LIMIT=MAG_LIMIT)
        model = unicorn.reduce.GrismModel(root=root, grow_factor=2, growx=2, growy=2, MAG_LIMIT=MAG_LIMIT)
        
    model.get_corrected_wcs(verbose=True)
    model.load_model_spectra()
    
    ### First iteration with flat spectra and the object flux
    model.compute_full_model(refine=False, MAG_LIMIT=28, save_pickle=False)   
    ### For the brighter galaxies, refine the model with the observed spectrum         
    model.compute_full_model(refine=True, MAG_LIMIT=24, save_pickle=True)
    
    err = np.random.normal(size=model.model.shape)*model.gris[2].data
    
    factors = [1,0.1,0.01,0.001,0.000001]
    
    yi, xi = np.indices(model.model.shape)
    
    fig = unicorn.catalogs.plot_init(xs=6, aspect=1, left=0.18)
    ax = fig.add_subplot(211)
    #ax.plot([0,0])
    
    fp = open(root+'_inter_deep_residuals.dat','w')
    fp.write('# pointing  sn_factor   mean  sigma\n')
    for factor in factors:
        mask = (model.model < (model.gris[2].data*factor)) & (model.gris[2].data != 0) & (xi > 350) & (yi > 120)
        #
        yh, xh = np.histogram(model.gris[1].data[mask].flatten(), range=(-0.05,0.05), bins=200, normed=True)
        stats = threedhst.utils.biweight(model.gris[1].data[mask].flatten(), both=True)
        a = ax.plot(xh[1:]*4, yh, linestyle='steps-', alpha=0.8, label=r'%5.3f %7.5f$\pm$%.3f' %(factor, stats[0]*4, stats[1]*4))
        fp.write('%s  %.2e  %.3e  %.3e\n' %(root, factor, stats[0]*4, stats[1]*4))
    #
    fp.close()
    
    ax.legend(prop=matplotlib.font_manager.FontProperties(size=10))
    ax.set_xlabel(r'$\delta$ electrons / s')
    ax.set_ylabel('N')
    
    # import gbb.pymc_gauss
    # gfit = gbb.pymc_gauss.init_model(xh[1:], yh, yh*0.+0.1)
    # gfit.sample(5000,1000)
    # plt.plot(xh[1:], yh-gfit.stats()['eval_gaussian']['mean'])
    
    ax = fig.add_subplot(212)
    ### See the flux level that corresponds to the over-subtraction
    wave = unicorn.reduce.sens_files['A']['WAVELENGTH']
    flam = np.abs(stats[0])*4/(unicorn.reduce.sens_files['A']['SENSITIVITY']*46.5)
    fnu = flam*wave**2/3.e18
    mag = -2.5*np.log10(fnu)-48.6
    ax.plot(wave, mag)
    ax.set_xlim(1.05e4,1.7e4)
    ran = (wave > 1.05e4) & (wave < 1.7e4)
    ax.set_ylim(mag[ran].min()-0.2, mag[ran].max()+0.2)
    ax.set_xlabel(r'$\lambda$') 
    ax.set_ylabel(r'mag of offset')
    
    unicorn.catalogs.savefig(fig, root+'_inter_deep_residuals.png')
    
##### Show parameter variation across the field
def show_conf_params():
    conf = unicorn.reduce.conf
    beam='A'
    bigX, bigY = np.indices((1014,1014))
    
    dydx_order = np.int(conf['DYDX_ORDER_'+beam])
    dydx_0 = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['DYDX_'+beam+'_0']) 
    dydx_1 = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['DYDX_'+beam+'_1'])
    
    plt.imshow(dydx_0 - dydx_0[507, 507], vmin=-1, vmax=1, origin='lower')
    cb = plt.colorbar()
    cb.set_label(r'$\Delta y$ (FLT pixels)')
    plt.scatter([507,507+125], [507,507+125], marker='o', color='black')
    plt.text(507+125, 507+125+20, r'%.2f pix' %(dydx_0[507+125, 507+125]-dydx_0[507,507]), ha='center', va='bottom')
    plt.savefig('/tmp/dy_FLT.pdf')
    plt.close()
    
    xoff_beam = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['XOFF_'+beam])
    yoff_beam = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['YOFF_'+beam])
    
    disp_order = np.int(conf['DISP_ORDER_'+beam])
    dldp_0 = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['DLDP_'+beam+'_0'])
    dldp_1 = unicorn.reduce.field_dependent_multi(bigX, bigY, conf['DLDP_'+beam+'_1'])

    plt.imshow(dldp_0 - dldp_0[507, 507], vmin=-50, vmax=50, origin='lower')
    cb = plt.colorbar()
    cb.set_label(r'$\Delta \lambda (\mathrm{\AA})$')
    plt.scatter([507,507+125], [507,507+125], marker='o', color='black')
    plt.text(507+125, 507+125+20, r'%.1f $\mathrm{\AA}$' %(dldp_0[507+125, 507+125]-dldp_0[507,507]), ha='center', va='bottom')
    plt.savefig('/tmp/dlam0_FLT.pdf')
    plt.close()
    
    plt.imshow(dldp_1 / dldp_1[507, 507], vmin=0.95, vmax=1.05, origin='lower')
    cb = plt.colorbar()
    cb.set_label(r'$\Delta$ dispersion (%)')
    plt.scatter([507,507+125], [507,507+125], marker='o', color='black')
    plt.text(507+125, 507+125+20, r'%.2f' %((dldp_1[507+125, 507+125]/dldp_1[507,507]-1)*100)+'%', ha='center', va='bottom')
    plt.savefig('/tmp/dlam1_FLT.pdf')
    plt.close()
    
def field_dependent_multi(xi, yi, str_coeffs):
    """ 
    Calculate field-dependent parameter for aXe conventions.
    
    a = a_0 + a_1 * xi + a_2 * yi + a_3 * xi**2 + a_4 * xi * yi + a_5 * yi**2
    
    """
    coeffs = np.cast[float](str_coeffs.split())
    xy = np.array([xi*0+1,xi,yi,xi**2,xi*yi,yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3])
    a = np.sum(coeffs*xy[0:len(coeffs),:].T, axis=-1)
    return a

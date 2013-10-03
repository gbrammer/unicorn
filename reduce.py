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

import numpy as np
import matplotlib.pyplot as plt

USE_PLOT_GUI = False

import pyfits

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
conf = threedhst.process_grism.Conf('WFC3.IR.G141.V2.5.conf', path=os.getenv('THREEDHST')+'/CONF/').params
conf_grism = 'G141'
sens_files = {}
for beam in ['A','B','C','D','E']:
    sens_files[beam] = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+conf['SENSITIVITY_'+beam])[1].data

ZPs = {'F105W':26.2687, 'F125W':26.25, 'F140W':26.46, 'F160W':25.96, 'F606W':26.486, 'F814W':25.937}
PLAMs = {'F105W':1.0552e4, 'F125W':1.2486e4, 'F140W':1.3923e4, 'F160W': 1.5369e4, 'F606W':5917.678, 'F814W':8059.761, 'F435W':4350., 'F775W':7750., 'F850LP':9000, 'ch1':3.6e4, 'ch2':4.5e4, 'K':2.16e4, 'U':3828., 'G':4870., 'R':6245., 'I':7676., 'Z':8872.}

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
    import pywcs
except:
    print 'No pywcs found.'
    
def set_grism_config(grism='G141'):
    import unicorn.reduce as red
    if red.conf_grism == grism:
        return None
    #
    red.conf = threedhst.process_grism.Conf('WFC3.IR.%s.V2.0.conf' %(grism), path=os.getenv('THREEDHST')+'/CONF/').params
    red.conf_grism = grism
    red.sens_files = {}
    for beam in ['A','B','C','D','E']:
        red.sens_files[beam] = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+conf['SENSITIVITY_'+beam])[1].data
    
    print 'Set grism configuration files: %s' %(red.conf_grism)
    
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
        model.get_corrected_wcs()
        model.make_wcs_region_file()
        
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
    asn_file = 'jbhm19010_asn.fits'
    asn_file = 'jbhm19020_asn.fits'

    asn = threedhst.utils.ASNFile(asn_file)
    
    xoff, yoff = unicorn.reduce.acs_interlace_offsets(asn_file, growx=1, growy=1)
    
    xoff = np.array([  0,  14,  -4, -18])
    yoff = np.array([ 0, 24, 30,  7])
    
    xoff[1] -= 2; yoff[1] += 2
    xoff[2] -= 2; yoff[2] += 1
    xoff[3] -= 1; yoff[3] -= 2
    
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
        
def acs_interlace_offsets(asn_file, growx=2, growy=2, path_to_flt='./'):
    """
    Get interlacing pixel offsets from postargs
    """
    # asn = threedhst.utils.ASNFile(asn_file)
    # xpos, ypos = [], []
    # for exp in asn.exposures:
    #     head = pyfits.getheader(os.path.join(path_to_flt, exp+'_flt.fits'))
    #     xpos.append(head['POSTARG1'])
    #     ypos.append(head['POSTARG2'])
    #
    ### WFC3 postargs
    xpos = [0, 1.355, 0.881, -0.474]
    ypos = [0.0, 0.424, 1.212, 0.788]
    
    a10 = 0.00
    a11 = 0.0494
    b10 = 0.0494
    b11 = 0.0040

    acsang = 92.16 - -45.123
    xpos_acs, ypos_acs = threedhst.utils.xyrot(np.array(xpos), np.array(ypos), acsang)
    xoff = xpos_acs / a11
    yoff = (ypos_acs-xoff*b11)/b10

    # xinter = -np.cast[int](np.round(xoff*10)/10.*growx)
    # yinter = -np.cast[int](np.round(yoff*10)/10.*growy)
    # 
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
    
def interlace_combine(root='COSMOS-1-F140W', view=True, use_error=True, make_undistorted=False, pad = 60, NGROW=0, ddx=0, ddy=0, growx=2, growy=2, auto_offsets=False):
    
    from pyraf import iraf
    from iraf import dither
    
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    
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
    
    #pad += NGROW*4
    
    N = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW), dtype=np.int)
        
    inter_sci = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW))
    inter_weight = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW))
    if use_error:
        inter_err = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW))
    
    xi+=pad/(2*growx)+NGROW
    yi+=pad/(2*growy)+NGROW
    
    #### From pieter
    dxs = np.array([0,-20,-13,7]) + np.int(np.round(xsh[0]))*0
    dys = np.array([0,-7,-20,-13]) + np.int(np.round(ysh[0]))*0
    
    #### GOODS-N field from Weiner program 11600
    if root.startswith('GOODS-N') | root.startswith('GNGRISM'):
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
    
    #### Erb quasar sightlines
    if flt[0].header['PROPOSID'] == 12471:
        dxs, dys = np.array([  0, -10])*growx, np.array([ 0, -7])
        
    if auto_offsets:
        xinter, yinter = red.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy)
        dxs = xinter + np.int(np.round(xsh[0]))*0
        dys = yinter + np.int(np.round(ysh[0]))*0
        
    dxs += ddx
    dys += ddy
    
    dxi = np.cast[np.int32](np.ceil(dxs/growx))
    dyi = np.cast[np.int32](np.ceil(dys/growy))
    
    #### Find hot pixels, which are flagged as cosmic 
    #### rays at the same location in each of the flt files.  Ignore others.
    hot_pix = np.zeros((1014,1014),dtype='int')
    
    for flt in asn.exposures:
        im = pyfits.open(flt+'_flt.fits')
        hot_pix += (im[3].data & 4096) / 4096
        
    hot_pix = hot_pix > 2
        
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
            header = red.scale_header_wcs(h1.copy(), factor=2, growx=growx, growy=growy, pad=pad, NGROW=NGROW)
            header.update('EXTNAME','SCI')
            header.update('PAD',pad)
            header.update('GROWX', growx)
            header.update('GROWY', growy)
            header.update('NGROW', NGROW)
            header.update('REFIMAGE', '', comment='Source detection image')
            header_wht = header.copy()
            header_wht.update('EXTNAME','ERR')
        
        dx = np.int(np.round((xsh[i]-xsh[0])*growx))
        dy = np.int(np.round((ysh[i]-ysh[0])*growy))
        #
        dx = dxs[i]
        dy = dys[i]
        #
        #use = ((im[3].data & 4096) == 0) & ((im[3].data & 4) == 0) #& (xi > np.abs(dx/2)) & (xi < (1014-np.abs(dx/2))) & (yi > np.abs(dy/2)) & (yi < (1014-np.abs(dy/2)))
        #        
        use = ((im[3].data & (4+32+16+2048+4096)) == 0) & (~hot_pix)
        
        ### debug
        # errs = im[2].data[use]
        # print im.filename(), errs.min(), np.median(errs), errs.max()
         
        inter_sci[yi[use]*growy+dy,xi[use]*growx+dx] += im[1].data[use]/im[2].data[use]**2
        inter_weight[yi[use]*growy+dy,xi[use]*growx+dx] += 1./im[2].data[use]**2
        
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
        
    image = pyfits.HDUList([hdu,sci,wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
    
    image.writeto(root+'_inter.fits', clobber=True)
    
    pyfits.writeto('inter_N.fits', data=N, header=image[1].header, clobber=True)
    
    if make_undistorted:
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

        im[1].header.update('CRPIX1', im[1].header['CRPIX1']+2*shift[0][0])
        im[1].header.update('CRPIX2', im[1].header['CRPIX2']+2*shift[1][0])
        im[2].header.update('CRPIX1', im[2].header['CRPIX1']+2*shift[0][0])
        im[2].header.update('CRPIX2', im[2].header['CRPIX2']+2*shift[1][0])

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
    
def scale_header_wcs(header, factor=2, growx=2, growy=2, pad=60, NGROW=0):
    """
    Take an FLT header but scale the WCS keywords to account for the images 
    expanded/contracted by a factor of 'factor'.
    """
    import numpy as np
    
    header.update('NAXIS1',header['NAXIS1']*growx+NGROW*growx*factor+pad)
    header.update('NAXIS2',header['NAXIS2']*growy+NGROW*growy*factor+pad)

    header.update('CRPIX1',header['CRPIX1']*growx+NGROW*growx+pad/2)
    header.update('CRPIX2',header['CRPIX2']*growy+NGROW*growy+pad/2)
    
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
    
    header.update('CD1_1', header['CD1_1']/growx)
    header.update('CD2_1', header['CD2_1']/growx)
    header.update('CD1_2', header['CD1_2']/growy)
    header.update('CD2_2', header['CD2_2']/growy)
    
    header.update('IDCTAB','')
    
    header.update('EXPTIME',1)
    
    return header
    
def grism_model(xc_full=244, yc_full=1244, lam_spec=None, flux_spec=None, grow_factor=2, growx=2, growy=2, pad = 60, BEAMS=['A','B','C','D'], dydx=True, grism='G141', smooth_binned_sigma=0.001, xmi=1000, xma=-1000):
    """
    Main code for converting the grism dispersion calibration to a full 2D model spectrum
    
    The input lam_spec / flux_spec model is interpolated to the wavelength grid.  Need to 
    smooth it by a factor at least of order 1 pixel to get emission lines right 
    (smooth_binned_sigma)
    """
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    import scipy.ndimage as nd
    #ds9 = threedhst.dq.myDS9()
    
    ## A star in COSMOS-15-F140W
    ##### For interlaced images
    
    # xc_full = 244.13
    # yc_full = 1244.323
    
    xc = np.int(np.round((xc_full - pad/2)/growx))
    yc = np.int(np.round((yc_full - pad/2)/growy))
    
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
    if 'A' in BEAMS:
        xmi, xma = min(xmi,10), max(xma,213)
    if 'B' in BEAMS:
        xmi, xma = min(xmi,-210), max(xma, -170)
    if 'C' in BEAMS:
        xmi, xma = min(xmi,207), max(xma,454)
    if 'D' in BEAMS:
        xmi, xma = min(xmi, 469), max(xma,668)
    if 'E' in BEAMS:
        xmi, xma = min(xmi, -600), max(xma,-400)
        
    NX,NY = xma-xmi, 8

    NX,NY = xma-xmi, 15
    
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
    
    # This wasn't needed and was giving dz/1+z = +0.0035 offset w.r.t spec-zs! 
    # xarr += 1
    
    for ib, beam in enumerate(BEAMS):
       
        xmi_beam, xma_beam = tuple(np.cast[int](conf['BEAM'+beam].split())) 
        xmi_beam *= growx
        xma_beam *= growx
        
        dydx_order = np.int(conf['DYDX_ORDER_'+beam])
        dydx_0 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_0']) * growy
        dydx_1 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_1']) * growy / growx#/ grow_factor        
                    
        xoff_beam = field_dependent(bigX, bigY, conf['XOFF_'+beam]) * growx
        yoff_beam = field_dependent(bigX, bigY, conf['YOFF_'+beam]) * growy
                
        disp_order = np.int(conf['DISP_ORDER_'+beam])
        dldp_0 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_0'])
        dldp_1 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_1']) / growx
        
        #### Improve alignment of zeroth order
        if (beam == 'B') & (grism == 'G141'):
            dldp_0+=1500
            
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
        lam = dldp_0 + dldp_1*(xarr-xoff_beam)
                
        if not dydx:
            dydx_0 = 0.
            dydx_1 = 0.
        
        #print 'DYDX: ', dydx, dydx_0, dydx_1
        
        #### Interpolate pixel at shifted ycenter along beam        
        ycenter = dydx_0 + dydx_1*(xarr-xoff_beam)
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
            #print 'Yoff_array', xarr.shape, yoff_array.shape
            
        #### Sensitivity
        #sens = pyfits.open(os.getenv('THREEDHST')+'/CONF/'+conf['SENSITIVITY_'+beam])[1].data
        sens = sens_files[beam]
        #dl = np.diff(lam)[0]
        #if dl > 0:
        sens_interp = np.interp(lam, sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)
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
            dl = lam[1]-lam[0]
            if dl > 0:
                spec_interp = utils_c.interp_conserve_c(lam, lam_spec, flux_spec, left=0., right=0.)
            else:
                spec_interp = utils_c.interp_conserve_c(lam[::-1], lam_spec, flux_spec, left=0., right=0.)[::-1]
                
            #### Smooth input model
            # xg = np.arange(-5*smooth_binned_sigma, 5.1*smooth_binned_sigma, 1)
            # yg = np.exp(-xg**2/2/smooth_binned_sigma**2)
            # yg /= np.sum(yg)
            # spec_interp = nd.convolve1d(spec_interp, yg)

            spec_interp[~np.isfinite(spec_interp)] = 0
            spec_interp[lam < np.min(lam_spec)] = 0
            
            sens_interp *= spec_interp
        
        ##### Needed these factors when I was having unit problems and not including
        ##### the "sensitivity per A" but instead assuming "sensitivity per pixel"
        #if beam == 'A':
        #    sens_interp *= 0.5
        
        #### Increase response in the zeroth order to match observations
        if (beam == 'B') & (grism == 'G141'):
            sens_interp *= 1.5 #3.6 * 10
        
        # if beam in ['C','D','E']:
        #     sens_interp *= 0.25
                
        stripe *= np.dot(np.ones((NY*2+1,1)), sens_interp.reshape(1,NX)) * (growx*growy)
        
        model += stripe
        #print beam, model.max()
        
    if 'A' not in BEAMS:
        yoff_array = np.array([0,0])
        
    return model, (xmi, xma, wavelength, full_sens, yoff_array, beam_index)
    
def process_GrismModel(root='GOODS-S-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=27, REFINE_MAG_LIMIT=23,  make_zeroth_model=True, use_segm=False, model_slope=0, model_list = None, grism='G141'):
    import unicorn.reduce
    
    model = unicorn.reduce.GrismModel(root=root, grow_factor=grow_factor, growx=growx, growy=growy, MAG_LIMIT=MAG_LIMIT, use_segm=use_segm, grism=grism)
    
    model.model_list = model_list
    
    if not use_segm:
        model.get_corrected_wcs(verbose=True)
    
    test = model.load_model_spectra()
    if not test:
        if make_zeroth_model:
            ### "zeroth" iteration to flag 0th order contamination, no color / norm iteration
            model.compute_full_model(refine=False, MAG_LIMIT=MAG_LIMIT, save_pickle=True, BEAMS=['B'])   
            ### Reset the model
            model.init_object_spectra()
            model.model*=0
        
        ### First iteration with flat spectra and the object flux
        model.compute_full_model(refine=False, MAG_LIMIT=MAG_LIMIT, save_pickle=False, model_slope=model_slope)   
        ### For the brighter galaxies, refine the model with the observed spectrum         
        model.compute_full_model(refine=True, MAG_LIMIT=REFINE_MAG_LIMIT, save_pickle=True, model_slope=model_slope)
        
    return model
    
class Interlace1D():
    def __init__(self, file='UDS-17_00319.1D.fits', PNG=True, flux_units=True):
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
    def __init__(self, file='GOODS-S-34_00446.2D.fits', PNG=True, growx=2, growy=2):
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
        
        unicorn.reduce.set_grism_config(self.grism_element)
        
        self.oned = unicorn.reduce.Interlace1D(file.replace('2D','1D'), PNG=False)
            
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
        #
        #### XXX need to get the zeropoint for arbitrary filter
        self.flux = self.thumb * 10**(-0.4*(ZPs[self.direct_filter]+48.6))* 3.e18 / PLAMs[self.direct_filter]**2 / 1.e-17
        self.total_flux = np.sum(self.flux*(self.seg == self.id))
        #self.total_flux = np.sum(self.flux*(self.seg > 0))
        
        self.init_model()
        
    def init_model(self, lam_spec=None, flux_spec=None):
        """
        Initialize all of the junk needed to go from the pixels in the 
        direct thumbnail to the 2D model spectrum
        """
        orders, self.xi = unicorn.reduce.grism_model(self.x_pix, self.y_pix, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=['A'], grow_factor=self.grow_factor, growx = self.growx, growy=self.growy, pad=self.pad)
        
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
        
        status = utils_c.disperse_grism_object(self.obj_flux, self.id, self.obj_seg, xord, yord, ford, self.object)
        
        xxi = np.int(np.round(sh[0]/2.))+xord
        use = (xxi >= 0) & (xxi < (NX))
        ### First order only
        #use = use & (xord > 10*self.grow_factor) & (xord < 213*self.grow_factor) 
        use_order = use 
        cut = np.zeros(NX)
        cut[xxi[use_order]] = word[use_order] 
        self.object_wave = cut.copy() #np.dot(np.ones((self.sh[0],1)),
        
        xarr = np.arange(NX)
        wavelength_region = (self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)
        if self.grism_element == 'G102':
            wavelength_region = (self.object_wave >= 0.75e4) & (self.object_wave <= 1.15e4)
            
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
        
        # xoff_arr = np.arange(len(self.xi[4]), dtype=np.double)+self.xi[0]+sh[0]/2.
        # xoff_int = np.arange(xmin, xmax, dtype=np.double)
        # yoff_int = utils_c.interp_c(xoff_int, xoff_arr, self.xi[4])        
        #self.yc = np.int(np.round(yoff_int.mean()))
        self.yc = -self.im['SCI'].header['YTOFF']
        self.ypad = ypad
        
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
        
        orders, self.xi = unicorn.reduce.grism_model(self.x_pix, self.y_pix, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=['A'], grow_factor=self.grow_factor, growx=self.growx, growy=self.growy, pad=self.pad, grism=self.grism_element, smooth_binned_sigma=smooth_binned_sigma)
        
        #print orders.shape
        
        if (lam_spec is not None):
            ford = np.array(orders[self.non_zero], dtype=np.float64)
        else:
            ford = self.ford
            
        ### Integrated over the thumbnail
        self.object *= 0.
        status = utils_c.disperse_grism_object(self.obj_flux, self.id, self.obj_seg, self.xord, self.yord, ford, self.object)
        
        ### Done!
        self.model = self.object[self.yc+self.ypad:self.yc-self.ypad,self.xmin:self.xmax]
    
    def optimal_extract(self, input=None):
        
        obj_cleaned = input
        variance = self.im['WHT'].data**2
            
        #### Get extraction window where profile is > 10% of maximum 
        osh = self.im['MODEL'].data.shape
        xarr = np.arange(osh[1])
        xint14 = int(np.interp(1.3e4, self.im['WAVE'].data, xarr))
        if self.grism_element == 'G102':
            xint14 = int(np.interp(0.98e4, self.im['WAVE'].data, xarr))
        
        yprof = np.arange(osh[0])
        profile = np.sum(self.im['MODEL'].data[:,xint14-10:xint14+10], axis=1)
        profile = profile/profile.sum()
        
        #profile2D = np.dot(profile.reshape((-1,1)), np.ones((1,osh[1])))
        window = profile >= 0.1*profile.max()
        
        #### Use object as the profile since it contains the y shift along the trace
        prof_x = np.sum(self.im['MODEL'].data, axis=0)
        profile2D = self.im['MODEL'].data/np.dot(np.ones((osh[0],1)), prof_x.reshape((1,-1)))
                
        obj_cleaned[variance == 0] = 0
        
        #### Optimal extraction
        opt_variance = variance.copy()
        opt_variance[opt_variance == 0] = 1.e6
        
        num = np.sum(profile2D*obj_cleaned/opt_variance, axis=0)        
        denom = np.sum(profile2D**2 / opt_variance, axis=0)
        
        optimal_sum = num / denom
        return self.im['WAVE'].data, optimal_sum
    
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
        xint14 = int(np.interp(1.3e4, self.wave, xarr))
        if self.grism_element == 'G102':
            xint14 = int(np.interp(0.98e4, self.wave, xarr))
            
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
                
        
class GrismModel():
    def __init__(self, root='GOODS-S-24', grow_factor=2, growx=2, growy=2, MAG_LIMIT=24, use_segm=False, grism='G141'):
        """
        Initialize: set padding, growth factor, read input files.
        """
        self.root=root
        self.grow_factor = grow_factor
        self.growx = growx
        self.growy = growy
        self.use_segm = use_segm
        
        #### Read the direct interlaced image        
        self.direct = pyfits.open(self.root+'-F140W_inter.fits')
        self.direct[1].data = np.cast[np.double](self.direct[1].data)
        
        self.filter = self.direct[0].header['FILTER']
        if 'PAD' in self.direct[1].header.keys():
            self.pad = self.direct[1].header['PAD']
        else:
            self.pad = 0.
            
        self.REF_INTER = False
        
        #### If a "reference" interlaced image exists, use that
        if os.path.exists(self.root+'_ref_inter.fits'):
            self.REF_INTER = True
            self.im = pyfits.open(self.root+'_ref_inter.fits')
            #self.pad += 4*self.im[1].header['NGROW']
            self.filter = self.im[0].header['FILTER']
        else:
            self.im = self.direct
        #
        try:
            self.growx = self.im[0].header['GROWX']
            self.growy = self.im[0].header['GROWY']
        except:
            self.growx = growx
            self.growy = growy
        
        if 'NGROW' not in self.direct[1].header.keys():
            self.ngrow = (self.direct[1].data.shape[1] - 1014*self.growx - self.pad) / (self.growx*2) 
        else:
            self.ngrow = self.direct[1].header['NGROW']
            
        self.ZP = ZPs[self.filter]

        self.PLAM = PLAMs[self.filter]
                    
        self.make_flux()
        
        self.grism_element = grism
        unicorn.reduce.set_grism_config(grism)
        
        self.gris = pyfits.open(self.root+'-%s_inter.fits' %(self.grism_element))
        if self.gris[1].data.shape != self.direct[1].data.shape:
            threedhst.showMessage('G141 and F140W images have different dimensions!\n')
        self.gris[1].data = np.array(self.gris[1].data, dtype=np.double)
        self.gris[2].data = np.array(self.gris[2].data, dtype=np.double)
        self.sh = self.im[1].data.shape
        
        if not os.path.exists(root+'_inter_seg.fits'):
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
        self.cat = threedhst.sex.mySexCat(self.root+'_inter.cat')
        
        #self.trim_edge_objects(verbose=True)
        
        if 'X_FLT' in self.cat.column_names:
            self.cat.x_pix = np.cast[float](self.cat.X_FLT)
            self.cat.y_pix = np.cast[float](self.cat.Y_FLT)
        else:
            self.cat.x_pix = np.cast[float](self.cat.X_IMAGE)
            self.cat.y_pix = np.cast[float](self.cat.Y_IMAGE)
            
        self.cat.mag = np.cast[float](self.cat.MAG_AUTO)
        
        self.segm = pyfits.open(self.root+'_inter_seg.fits')
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
        se.options['CATALOG_NAME']    = self.root+'_inter.cat'
        se.options['CHECKIMAGE_NAME'] = self.root+'_inter_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = self.root+'-F140W_inter.fits[1]'
        se.options['FILTER']    = 'Y'
        
        #### Detect thresholds (default = 1.5)
        se.options['DETECT_THRESH']    = '1' 
        se.options['ANALYSIS_THRESH']  = '1'
        se.options['MAG_ZEROPOINT'] = '%.2f' %(self.ZP)
        
        #### Run SExtractor
        status = se.sextractImage(self.root+'-F140W_inter.fits[0]')
        self.cat = threedhst.sex.mySexCat(self.root+'_inter.cat')
        
        #### Trim faint sources
        mag = np.cast[float](self.cat.MAG_AUTO)
        q = np.where(mag > MAG_LIMIT)[0]
        if len(q) > 0:
            numbers = np.cast[int](self.cat.NUMBER)[q]
            self.cat.popItem(numbers)
        
        self.cat.write()
        
        #self.trim_edge_objects(verbose=True)
        
        #self.segm = pyfits.open(self.root+'_seg.fits')
        
        threedhst.sex.sexcatRegions(self.root+'_inter.cat', self.root+'_inter.reg', format=1)
        
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
        
        if os.path.exists('%s_inter.cat.wcsfix' %(self.root)):
            wcs = catIO.Readfile('%s_inter.cat.wcsfix' %(self.root))
            test = wcs.id == self.cat.id
            if np.sum(test) == len(self.cat.id):
                if verbose:
                    print 'Get corrected coordinates from %s_inter.cat.wcsfix' %(self.root)
                self.ra_wcs = wcs.ra
                self.dec_wcs = wcs.dec
                return True
        
        if self.REF_INTER:
            self.ra_wcs, self.dec_wcs = self.cat.ra*1., self.cat.dec*1.       
            
            #### Write the wcsfix file
            fp = open('%s_inter.cat.wcsfix' %(self.root),'w')
            fp.write('# id    mag   ra    dec\n')
            for i in range(self.cat.id.shape[0]):
                fp.write('%5d  %.3f  %.6f  %.6f\n' %(self.cat.id[i], self.cat.mag[i], self.ra_wcs[i], self.dec_wcs[i]))

            fp.close()

            if verbose:
                print 'Corrected coordinates: %s_inter.cat.wcsfix' %(self.root)

            return True
        
        #### Else: use tran to get the pixel in the DRZ frame and then xy2rd to get coords        
        wcs = wcsutil.WCSObject(self.root+'-F140W_drz.fits[1]')
            
        NOBJ = len(self.cat.id)
        
        asn = threedhst.utils.ASNFile(self.root+'-F140W_asn.fits')
        flt = asn.exposures[0]+'_flt.fits'
        fp = open('/tmp/%s.flt_xy' %(self.root),'w')
        for i in range(NOBJ):
            #1014*growy+pad+growy*2*NGROW
            xi = (self.cat.x_pix[i]-self.pad/2. - self.ngrow*self.growx)/self.growx
            yi = (self.cat.y_pix[i]-self.pad/2. - self.ngrow*self.growy)/self.growy
            fp.write('%.2f %.2f\n' %(xi,  yi))
        
        fp.close()
        
        print 'Running iraf.tran()'
        threedhst.process_grism.flprMulti()
        
        ### For some reason this dies occasionally "too many positional arguments for traxy"
        ### Try running a second time if it dies once
        try:
            status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-F140W_drz.fits[1]', direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
            os.remove("/tmp/%s.flt_xy" %(self.root))
        except:
            threedhst.process_grism.flprMulti()
            status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-F140W_drz.fits[1]', direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
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
        fp = open('%s_inter.cat.wcsfix' %(self.root),'w')
        fp.write('# id    mag   ra    dec\n')
        for i in range(self.cat.id.shape[0]):
            fp.write('%5d  %.3f  %.6f  %.6f\n' %(self.cat.id[i], self.cat.mag[i], self.ra_wcs[i], self.dec_wcs[i]))
            
        fp.close()
        
        if verbose:
            print 'Corrected coordinates: %s_inter.cat.wcsfix' %(self.root)
            
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
        
        if self.grism_element == 'G141':
            self.lam_spec = np.arange(0.95e4,1.8e4,22.)
        
        if self.grism_element == 'G102':
            self.lam_spec = np.arange(0.7e4, 1.25e4, 12.)
        
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
        orders, self.xi = unicorn.reduce.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=self.grow_factor, growx=self.growx, growy=self.growy, pad=self.pad, grism=self.grism_element)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Define grism dispersion (%.2f s).' %(dt)
            
        yord, xord = np.indices(orders.shape)
        beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), self.xi[5].reshape((1,-1)))
        
        non_zero = orders > 0
        #xord, yord, ford, word, sord = xord[non_zero], yord[non_zero], orders[non_zero], self.xi[2][non_zero], self.xi[3][non_zero]
        cast = np.int
        xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), self.xi[2][non_zero], self.xi[3][non_zero], beams[non_zero]
                
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
        
        status = utils_c.disperse_grism_object(self.flux, id, self.segm[0].data, xord, yord, ford, self.object)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Loop (c). (%.3f)' %(dt)
        
        #### Make full images that show the wavelength, sensitivity
        xxi = np.int(np.round(xc))+xord
        use = (xxi >= 0) & (xxi < self.sh[1])
        ### First order only
        #use = use & (xord > 10*self.grow_factor) & (xord < 213*self.grow_factor) 
        use_order = use &  (bord == 1)
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
        if view.__class__ == threedhst.dq.myDS9:
            first_model = self.object*1
        
        if self.obj_in_model[id]:
            self.model -= self.object
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Refine, first_model  (%.3f)' %(dt)
            
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
        ratio_extract = utils_c.get_model_ratio_optimal(self.object, self.model, self.gris[1].data, self.gris[2].data)
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Make ratio  (%.3f)' %(dt)
        
        #wave = np.mean(self.object_wave, axis=0)
        wave = self.object_wave
        xpix = np.arange(self.sh[1])
        
        if self.grism_element == 'G141':
            keep = (wave > 1.12e4) & (wave < 1.65e4) & (ratio_extract != 0)  # (xpix > (self.pad-22)) & (xpix < (self.sh[1]-self.pad-22))
        else:
            keep = (wave > 0.8e4) & (wave < 1.10e4) & (ratio_extract != 0)
            
        #print len(wave), len(wave[wave > 0]), wave.max(), wave[2064], len(ratio_extract)
        
        if keep.sum() < 10:
            self.model += self.object
            self.obj_in_model[id] = True
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
        if view.__class__ == threedhst.dq.myDS9:
            second_model = self.object.copy()
        
        #### View the results
        if view.__class__ == threedhst.dq.myDS9:
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
        if view.__class__ == threedhst.dq.myDS9:
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
            
    def compute_full_model(self, BEAMS=['A','B','C','D'], view=None, MAG_LIMIT=23., save_pickle=True, refine=True, model_slope=0):
        
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
            if refine:
                self.refine_model(id, BEAMS=BEAMS, view=view)      
            else:
                if self.grism_element == 'G141':
                    lx = np.arange(0.9e4,1.8e4)
                    ly = 1+model_slope*(lx-1.4e4)/4.e3
                
                if self.grism_element == 'G102':
                    lx = np.arange(0.7e4,1.2e4)
                    ly = 1+model_slope*(lx-0.98e4)/2.8e3
                
                if self.model_list is not None:
                    if id in self.model_list.keys():
                        print unicorn.noNewLine+'Object #%d, m=%.2f (%d/%d) [model_list]' %(id, mag[so][i], i+1, N)
                        lx, ly = self.model_list[id]
                        self.flux_specs[id] = utils_c.interp_conserve_c(self.lam_spec, lx, ly, left=0., right=0.)
                        
                self.compute_object_model(id, BEAMS=BEAMS, lam_spec=lx, flux_spec=ly, normalize=True)      
                self.model += self.object
                self.obj_in_model[id] = True
                
        if save_pickle:
            self.save_model_spectra(BEAMS=BEAMS)
            
    def save_model_spectra(self, BEAMS=['A','B','C','D']):
        """
        Save model "pickle" and the 2D model image.
        
        If just BEAM ["B"] is specified, save the 0th order contamination image.
        """
        import pickle
        
        ####
        if BEAMS == ['B']:
            pyfits.writeto(self.root+'_inter_0th.fits', data=self.model, header=self.gris[1].header, clobber=True)
            return
            
        fp = open(self.root+'_inter_model.pkl','wb')
        pickle.dump(self.cat.id, fp)
        pickle.dump(self.obj_in_model, fp)
        pickle.dump(self.lam_spec, fp)
        pickle.dump(self.flux_specs, fp)
        fp.close()
        
        pyfits.writeto(self.root+'_inter_model.fits', data=self.model, header=self.gris[1].header, clobber=True)
        
    def load_model_spectra(self):
        import pickle
        if ((not os.path.exists(self.root+'_inter_model.pkl')) | 
            (not os.path.exists(self.root+'_inter_model.fits'))):
            return False
            
        fp = open(self.root+'_inter_model.pkl','rb')
        ids = pickle.load(fp)
        test = ids == self.cat.id
        if np.sum(test) == len(self.cat.id):
            print 'Load spectra from pickle'
            self.obj_in_model = pickle.load(fp)
            self.lam_spec = pickle.load(fp)
            self.flux_specs = pickle.load(fp)
            
            im = pyfits.open(self.root+'_inter_model.fits')
            self.model = np.cast[np.double](im[0].data)
            im.close()
        else:
            print "ERROR: Object list in pickle doesn't match the current catalog."
            fp.close()
            return False
            
        fp.close()
        return True
        
    def twod_spectrum(self, id=328, grow=1, miny=18, maxy=None, CONTAMINATING_MAGLIMIT=23, refine=True, verbose=False, force_refine_nearby=False, USE_REFERENCE_THUMB=False, USE_FLUX_RADIUS_SCALE=3):
        """
        Extract a 2D spectrum from the interlaced image.
        
        `miny` is the minimum Y-size of the 2D spectrum.  If `miny` < 0: then 
         use that size explicitly rather than just setting a minimum.
         
        `CONTAMINATING_MAGLIMIT` is the limit for computing contaminating
        objets.
        
        `USE_REFERENCE_THUMB`: Use the reference image rather than the 
        nominal F140W image for the 2D thumbnail.
        """
        if id not in self.cat.id:
            print 'Object #%d not in catalog.' %(id)
            return False
            
        #
        ii = np.where(np.cast[int](self.cat.id) == id)[0][0]
        xc = np.int(np.round(self.cat.x_pix[ii]))#-1
        yc = np.int(np.round(self.cat.y_pix[ii]))#-1
        
        seg = np.cast[int](self.segm[0].data)
        # seg_mask = seg == id
        # ypix, xpix = self.yf, self.xf
        # 
        # NY = ypix[seg_mask].max()-ypix[seg_mask].min()
        # NX = xpix[seg_mask].max()-xpix[seg_mask].min()
        #### Use only contiguous segmentation region around center pixel
        xmi, xma, ymi, yma = threedhst.utils.contiguous_extent(seg == id, xc, yc)
        NX = np.maximum(xc-xmi, xma-xc)*2
        NY = np.maximum(yc-ymi, yma-yc)*2
        
        NT = np.max([NY*grow, NX*grow, miny])
        if (USE_FLUX_RADIUS_SCALE is not None) & ('FLUX_RADIUS' in self.cat.column_names):
            NT = np.max([USE_FLUX_RADIUS_SCALE*2*self.cat['FLUX_RADIUS'][ii], miny])
        
        if miny < 0:
            NT = -miny
        
        #NT = np.min([self.pad/2, NT])
                                
        NT = np.min([NT, xc, yc, self.sh[1]-xc, self.sh[0]-yc])
        
        #### Maximum size
        if maxy is not None:
            NT = np.min([NT, maxy])
        
        #print '!!! NT: %d' %(NT)
        
        if verbose:
            print 'Generating direct image extensions'
            
        self.direct_thumb = self.direct[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_wht = self.direct[2].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_seg = np.cast[np.int32](seg[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2])
                    
        #### Initialize the output FITS file
        prim = pyfits.PrimaryHDU()
        prim.header.update('POINTING', self.root)
        prim.header.update('ID', id)
        prim.header.update('X_PIX', self.cat.x_pix[ii], comment='X pixel in interlaced image')
        prim.header.update('Y_PIX', self.cat.y_pix[ii], comment='Y pixel in interlaced image')
        prim.header.update('RA', self.ra_wcs[ii], comment='Target R.A.')
        prim.header.update('DEC', self.dec_wcs[ii], comment='Target Dec.')
        prim.header.update('MAG', self.cat.mag[ii], comment='MAG_AUTO from interlaced catalog')
        prim.header.update('FTOTAL', self.total_fluxes[id], comment='Total flux within segmentation image (1E-17)')
        prim.header.update('PAD', self.pad, comment='Padding at edge of interlaced image')
        prim.header.update('GROW', self.grow_factor, comment='"grow_factor" from unicorn.reduce')
        prim.header.update('GROWX', self.growx, comment='"growx" from unicorn.reduce')
        prim.header.update('GROWY', self.growy, comment='"growy" from unicorn.reduce')
        prim.header.update('BUNIT', 'ELECTRONS/S')
        prim.header.update('EXPTIME', self.gris[0].header['EXPTIME'])
        prim.header.update('FILTER', self.filter)
        prim.header.update('GRISM', self.gris[0].header['FILTER'])
        
        prim.header.update('REFTHUMB', False, comment='Thumbnail comes from reference image')
        
        if '_ref_inter' in self.im.filename():
            ### Use reference image as thumbnail if it exists
            if (xc < (self.pad / 2)) | USE_REFERENCE_THUMB:
                prim.header.update('REFTHUMB', True)
                self.direct_thumb = self.im[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
                prim.header.update('FILTER', self.im[0].header['FILTER'])
                
        extensions = [prim]

        header = pyfits.ImageHDU().header
        
        #### Direct thumbnails
        header.update('EXTNAME','DSCI')
        extensions.append(pyfits.ImageHDU(self.direct_thumb, header))
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
        
        BEAMS=['A','B','C','D']
        view=False
        if nearby.sum() > 0:
            ids = self.cat.id[nearby][np.argsort(self.cat.mag[nearby])]
            for idi in ids:
                if verbose:
                    print 'Calculating contaminating object #%d' %(idi)
                #
                if refine:
                    if (not self.obj_in_model[idi]) | force_refine_nearby:
                        self.refine_model(idi, BEAMS=BEAMS, view=view)      
                else:
                    if (not self.obj_in_model[idi]) | force_refine_nearby:
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
        
        #### Find pixels of the 1st order        
        xarr = np.arange(self.sh[0])
        if self.grism_element == 'G141':
            wavelength_region = (self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)
        #
        if self.grism_element == 'G102':
            wavelength_region = (self.object_wave >= 0.75e4) & (self.object_wave <= 1.15e4)
        
        if np.sum(wavelength_region) < 20:
            fp = open(self.root+'_%05d.2D.xxx' %(id),'w')
            fp.write('%.2f %.2f  %d %d\n' %(xc, yc, -1, -1))
            fp.close()
            self.twod_status = False
            return False
        
        self.twod_status = True
        
        xmin = xarr[wavelength_region].min()
        xmax = xarr[wavelength_region].max()
        
        if ((self.sh[1]-xmin) < 20) | ((xmax-xmin) < 20):
            fp = open(self.root+'_%05d.2D.xxx' %(id),'w')
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
        
        if verbose:
            print 'Writing spectrum extensions.'
            
        ### WCS header
        header.update('CRPIX1',1)
        header.update('CDELT1',self.dwave)
        header.update('CRVAL1',self.wave[0])
        header.update('CUNIT1','Angstrom')
        header.update('CTYPE1','WAVE')
        
        header.update('CRPIX2',NT/2.+1)
        header.update('CRVAL2',0)
        header.update('CDELT2',0.128254/2)
        header.update('CUNIT2','arcsec')
        header.update('CTYPE2','CRDIST')
        
        header.update('YTOFF',spec_y_offset,comment='Offset of middle of spectrum trace')
        
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
        
        hdu = pyfits.HDUList(extensions)
        hdu.writeto(self.root+'_%05d.2D.fits' %(id), clobber=True)
        
        self.twod = hdu
        
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
        ax.text(1.1,0.8, self.root, transform=ax.transAxes)
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
            ax.set_xlim(1.1e4,1.65e4)
            if self.grism_element == 'G102':
                ax.set_xlim(0.75e4,1.12e4)
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
            canvas.print_figure(self.root+'_%05d.2D.png' %(self.id), dpi=100, transparent=False)
                        
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
        xint14 = int(np.interp(1.3e4, self.wave, xarr))
        if self.grism_element == 'G102':
            xint14 = int(np.interp(0.98e4, self.wave, xarr))
            
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
        tbHDU.writeto(self.root+'_%05d.1D.fits' %(self.id), clobber=True)
        
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
        
    def extract_spectra_and_diagnostics(self, list=None, skip=True, MAG_LIMIT=19, verbose=False, miny=26, USE_FLUX_RADIUS_SCALE=3):
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
            if (os.path.exists(self.root+'_%05d.2D.fits' %(id)) | os.path.exists(self.root+'_%05d.2D.xxx' %(id))) & skip:
                continue
            
            if verbose:
                print '2D FITS'
            status = self.twod_spectrum(id=id, verbose=verbose, miny=miny, refine=True, USE_FLUX_RADIUS_SCALE=USE_FLUX_RADIUS_SCALE)
            
            if status is False:
                if verbose:
                    print unicorn.noNewLine+'Object #%-4d (%4d/%4d) - No 2D' %(id,i+1,N)
                continue
                    
            if verbose:
                print unicorn.noNewLine+'2D FITS, 2D PNG, 1D FITS'
            self.show_2d(savePNG=True, verbose=verbose)

            if verbose:
                print unicorn.noNewLine+'2D FITS, 2D PNG, 1D FITS, 1D PNG'
            spec = unicorn.reduce.Interlace1D(self.root+'_%05d.1D.fits' %(id), PNG=True)
                
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
    
def interlace_goodsn():
    """
    Reduce the GOODS-N pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    
    #os.chdir(unicorn.GRISM_HOME+'GOODS-N/Interlace_GBB')
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/INTERLACE')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    #unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-N_F140W_v1', filter='F140W', REFERENCE = 'goodsn_for_arjen/goodsn_f140w_sci_sub.fits', SEGM = 'goodsn_for_arjen/goodsn_f140w_v1.seg.fits')
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-N_F140W', filter='F140W', REFERENCE = '/3DHST/Photometry/Release/GOODS-N/v1.8/detection/goodsn_f140w_sci_sub.fits', SEGM = '/3DHST/Photometry/Release/GOODS-N/v1.8/detection/goodsn_f140w_v1.seg.fits')
    
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Release/GOODS-N/v1.8/detection/goodsn_f140w_v1_det.cat'
    
    #direct=glob.glob('ib*50_asn.fits')
    direct = glob.glob('GOODS-N-*-F140W_asn.fits')

    extract_limit = 24
    skip_completed=True
    #REF_ROOT='GOODS-N_F140W_v1'
    REF_ROOT='GOODS-N_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        #
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)
    
    ##### Generate the spectral model and extract spectra
    inter = glob.glob('GOODS-N-1[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if not (os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=27.1)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=27.1)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-N-11')
        
    skip_completed = False
    
    #Fit the redshift for objects brighter than F140W=24.
    models = glob.glob('GOODS-N-3[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))# & (model.cat.mag > 17.)) #& (model.mag.id != 1876))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root,p_flat = 0.0)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('GOODS-N-4[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)


    ### Get some quality flags on the reduced spectra
    root='GOODS-N'
    
    files = glob.glob('GOODS-N-4*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])
            
    fp.close()
    fp2.close()
    
    ###################### Analysis plots ################
    ## 
    os.system('paste %s.zfit.dat %s.dqflag.dat > %s.specz.dat' %(root, root, root))
    
    zfit = catIO.Readfile('%s.specz.dat' %(root))
    
    #### v2.0 release
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-N')
    zfit = catIO.Readfile('GOODS-N.zfit.linematched.dat')
    dq = catIO.Readfile('GOODS-N.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm

    #### dq has some ****** IDL format values in max_contam
    lines = pyfits.open('GOODS-N.linefit.fits')
    cat = catIO.Readfile('goodsn_f140w_v1_det_reform.cat')
    root = 'goodsn'
    
    ####  COSMOS
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/COSMOS')
    zfit = catIO.Readfile('COSMOS.zfit.linematched.dat')
    dq = catIO.Readfile('COSMOS.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    lines = pyfits.open('COSMOS.linefit.fits')
    cat = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Catalog/3dhst.cosmos.v2.0.cat')
    cat.x_image, cat.y_image = cat.x, cat.y
    root = 'cosmos'
    fout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Fast/3dhst.cosmos.v2.0.fout')
    zout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Eazy/3dhst.cosmos.v2.0.zout')
    
    ####  GOODS-S
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-S')
    zfit = catIO.Readfile('GOODS-S.zfit.linematched.dat')
    dq = catIO.Readfile('GOODS-S.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    lines = pyfits.open('GOODS-S.linefit.fits')
    cat = catIO.Readfile('GOODS-S_v2.0_PHOTOMETRY/catalog/GOODS-S_v2.0.fullz_wzp.cat')
    cat.x_image, cat.y_image = cat.ximage, cat.yimage
    root = 'goodss'
    
    zfit.mag = dq.mag
    
    dz_spec = (zfit.z_max_spec-zfit.z_spec)/(1+zfit.z_spec)
    dz_peak = (zfit.z_peak_spec-zfit.z_spec)/(1+zfit.z_spec)
    sz = np.maximum((np.abs(dz_spec)/0.001)**(0.8),8)
    
    dz_phot = (zfit.z_peak_phot - zfit.z_spec)/(1+zfit.z_spec)
    
    fig = unicorn.plotting.plot_init(square=True, xs=4, left=0.13)
    ax = fig.add_subplot(111)
    ax.scatter(zfit.z_spec, zfit.z_max_spec, color='blue', alpha=0.2, s=2)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 4)
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{gris}$')
    unicorn.plotting.savefig(fig, root+'_zphot_zspec.png')
    
    z0 = (zfit.z_spec > 0) & (zfit.z_spec < 10)
    plt.scatter(zfit.z_spec, dz_spec, color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_spec[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_z_dz_max.png')
    
    plt.scatter(zfit.z_spec, dz_peak, color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_peak[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{peak}$')
    plt.savefig(root+'_z_dz_peak.png')

    plt.scatter(zfit.z_spec, dz_phot, color='orange', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_phot[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{phot}$')
    plt.savefig(root+'_z_dz_phot.png')
    
    zp6 = zfit.z_spec > 0.6
    keep = zp6
    
    plt.plot(dq.mag[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(18,24.5)
    plt.plot(dq.mag[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.1,0.1); plt.xlim(18,24.5)
    xm, ym, ys, N = threedhst.utils.runmed(dq.mag[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(dq.mag[keep], dz_phot[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    plt.xlabel(r'dqflag.mag'); plt.ylabel(r'$\Delta z_\mathrm{phot}$')
    plt.savefig(root+'_mag_dz_gris.png')
    keep = keep & (dq.mag < 24)
    
    plt.plot(dq.max_contam[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.02,0.02); plt.xlim(0.002, 30)
    xm, ym, ys, N = threedhst.utils.runmed(dq.max_contam[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.semilogx()
    plt.xlabel(r'dqflag.max_contam'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_max_contam_dz_gris.png')
    keep = keep & (dq.max_contam < 1)

    plt.plot(dq.int_contam[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.008,2)
    xm, ym, ys, N = threedhst.utils.runmed(dq.int_contam[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.semilogx()
    plt.xlabel(r'dqflag.int_contam'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_int_contam_dz_gris.png')
    
    keep = keep & (dq.int_contam < 0.6)
    
    plt.plot(dq.q_z[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.1,0.1); plt.xlim(0.001,100) ; plt.semilogx()
    plt.plot(dq.q_z[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.001,100) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.q_z[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    keep = keep & (zfit.q_z < 5)
    
    plt.plot(dq.f_cover[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_cover[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.01,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_cover[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.xlabel(r'dqflag.f_cover'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_cover_dz_gris.png')
    
    keep = keep & (dq.f_cover > 0.5)
    
    plt.plot(dq.f_flagged[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_flagged[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.01,0.01); plt.xlim(0.005,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_flagged[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')

    plt.xlabel(r'dqflag.f_flagged'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_flagged_dz_gris.png')

    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_negative = np.maximum(dq.f_negative, 1.e-2)*norm
    
    plt.plot(dq.f_negative[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_negative[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.01,0.01); plt.xlim(0.005,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_negative[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')

    plt.xlabel(r'dqflag.f_negative'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_negative_dz_gris.png')
    
    keep = keep & (dq.f_flagged < 0.55)
    nopartial = zfit.f_flagged < 0.2
    
    #### Ha eqw
    z15 = (zfit.z_spec > 0.7) & (zfit.z_spec < 1.5)
    eqw_rf = np.maximum(lines[1].data['HALPHA_EQW'] / (1+lines[1].data['Z']), 0.2)
    plt.plot(eqw_rf[z15], dz_spec[z15], color='blue', alpha=0.4); plt.ylim(-0.05,0.05); plt.xlim(0.1,5000) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(eqw_rf[z15], dz_spec[z15], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(eqw_rf[z15], dz_phot[z15], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    plt.xlabel(r'HALPHA_EQW / (1+z)'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_ha_eqw_dz_gris_all.png')
    
    plt.xlim(10,1000)
    plt.ylim(-0.01,0.01)
    plt.savefig(root+'_ha_eqw_dz_gris.png')
        
    d99 = ((zout.u99-zout.l99)/(1+zout.z_peak))[zfit.phot_id-1]
    plt.plot(d99[~keep], dz_spec[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(d99[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(d99[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    #### Do the optimal selection on the full catalog
    fig = unicorn.plotting.plot_init(square=True, xs=6, aspect=0.66, left=0.08)
    ax = fig.add_subplot(111)
    has_gris = zfit.z_max_spec > 0.01
    full_selection = has_gris & (dq.mag < 24) & (dq.max_contam < 1) & (dq.f_cover > 0.5) & (dq.f_flagged < 0.55) & (dq.f_negative < 0.55)
    yh, xh = np.histogram(np.log(1+zfit.z_max_spec[full_selection]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', label='grism', color='blue')
    yh, xh = np.histogram(np.log(1+zfit.z_max_spec[has_gris]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', label='full', color='blue', alpha=0.3)
    yh, xh = np.histogram(np.log(1+zfit.z_spec[full_selection]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, 0-yh, linestyle='steps', marker='None', color='green', label='spec')
    yh, xh = np.histogram(np.log(1+zfit.z_spec[has_gris]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, 0-yh, linestyle='steps', marker='None', color='green', label='spec/full', alpha=0.3)
    ax.legend()

    ax.set_ylim(-50,50)
    ax.set_xlim(0,3)
    
    if root=='cosmos':
        ax.set_ylim(-80,80)
        ax.set_xlim(0,3)
        
    ax.set_xlabel('z')
    unicorn.plotting.savefig(fig, root+'_full_zdist_compare_zspec.png')
    
    #### Figure looking for overdensities
    zmax = 3
    r0 = np.median(cat.ra)
    
    width = 0.005
    ztry = 0.850

    zr = (0.7,2.4)
    for count, lz in enumerate(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), width)):
        print count
        ztry = np.exp(lz)-1
        #
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        #
        fig = unicorn.plotting.plot_init(square=True, xs=5, aspect=1/0.7, NO_GUI=False)
        ax = fig.add_axes((0.05, 0.8, 0.90, 0.19))
        #
        yh, xh = np.histogram(np.log(1+zfit.z_max_spec[has_gris]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None')
        yh, xh = np.histogram(np.log(1+zfit.z_max_spec[bin]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', color='red', linewidth=2)
        a = ax.set_yticklabels([])
        #
        ax = fig.add_axes((0.05, 0.05, 0.90, 0.70))
        a = ax.plot(r0-cat.ra, cat.dec, ',', color='0.8', alpha=0.1)
        a = ax.plot(r0-cat.ra[bin], cat.dec[bin], 'o', alpha=0.4, color='red', ms=5)
        a = ax.set_yticklabels([])
        a = ax.set_xticklabels([])
        #
        unicorn.plotting.savefig(fig, 'goodsn_zhist_%04d.png' %(count+1))
    
    ### nearest neighbor density plot
    width = 0.005
    NX, NY = cat.x_image.max(), cat.y_image.max()
    fact = 100
    yi, xi = np.indices((int(NY/fact),int(NX/fact)))
        
    import scipy.spatial

    bin = cat.x_image > 0
    xy = np.array([cat.x_image[bin]/fact, cat.y_image[bin]/fact])
    tree = scipy.spatial.cKDTree(xy.T, 10)
    #
    N = 7
    mask = xi*0.
    for x in xrange(int(NX/fact)):
        #print unicorn.noNewLine+'%d' %(x)
        for y in xrange(int(NY/fact)):
            mask[y,x] = tree.query([x,y], k=N)[0][-1]
    
    mask = mask < 3
    
    ### which redshift to use
    zuse = zfit.z_max_spec
    has_gris = zuse > 0.01
    
    zuse = zout.z_peak
    has_gris = (dq.mag < 25) & (cat.use == 1)
    
    aspect, xs = 1/0.6, 5
    max_scale = np.log(0.2)
            
    if root == 'cosmos':
        aspect = 1./0.4
        xs = 5
        max_scale = np.log(0.1)
        
    zr = (0.65, 2.4)
    
    for count, lz in enumerate(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), width)):
        print unicorn.noNewLine+'%d %f' %(count, lz)
        ztry = np.exp(lz)-1
        #
        bin = has_gris & (np.abs(zuse-ztry)/(1+ztry) < width)
        #
        fig = unicorn.plotting.plot_init(square=True, xs=xs, aspect=aspect, NO_GUI=True)
        ax = fig.add_axes((0.05, 0.8, 0.90, 0.19))
        #
        yh, xh = np.histogram(np.log(1+zuse[has_gris]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None')
        yh, xh = np.histogram(np.log(1+zuse[bin]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', color='red', linewidth=2)
        a = ax.set_yticklabels([])
        a = ax.text(0.95, 0.8, 'z=%.3f' %(ztry), transform=ax.transAxes, horizontalalignment='right', fontsize=12)
        #
        ax = fig.add_axes((0.05, 0.05, 0.90, 0.70))
        #
        xy = np.array([cat.x_image[bin]/fact, cat.y_image[bin]/fact])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        #
        R = xi*0.
        for x in xrange(int(NX/fact)):
            #print unicorn.noNewLine+'%d' %(x)
            for y in xrange(int(NY/fact)):
                R[y,x] = tree.query([x,y], k=N)[0][-1]
        #
        density = N/np.pi/R**2
        #a = ax.imshow(0-np.log(density*mask), interpolation='Nearest', vmin=max_scale, vmax=6.5)
        #### Uncertainty
        # bin2 = has_gris & (np.abs(zuse-ztry)/(1+ztry) < width)
        # xy = np.array([cat.x_image[bin2]/fact, cat.y_image[bin2]/fact])
        # tree = scipy.spatial.cKDTree(xy.T, 10)
        # Nbin = bin2.sum()
        # points = np.zeros(Nbin)
        # for i in range(Nbin):
        #     points[i] = tree.query([cat.x_image[bin2][i]/fact, cat.y_image[bin2][i]/fact], k=N+1)[0][-1]
        # #
        # sigma = np.std(N/np.pi/points**2) 
        density[~mask] = -10*sigma
        #ai = ax.imshow(0-density/sigma, interpolation='Nearest', vmin=-10, vmax=0)
        #ai = ax.imshow(density/sigma, interpolation='Nearest', vmin=0, vmax=10)
        ai = ax.imshow(np.log(density/sigma), interpolation='Nearest', vmin=-6, vmax=max_scale)
        a = ax.scatter(xy[0,:], xy[1,:], alpha=0.4, color='white', s=5)
        a = ax.set_yticklabels([])
        a = ax.set_xticklabels([])
        a = ax.set_xlim(0, R.shape[1])#
        a = ax.set_ylim(0, R.shape[0])
        # ac = fig.add_axes((0.03, 0.2, 0.07, 0.5))
        # cb = plt.colorbar(ai, cax=ac)
        # cb.set_label(r'$\Sigma / \sigma$')
        #
        unicorn.plotting.savefig(fig, root+'_nn_%5.3f.png' %(ztry))
        
    #### Add zFOURGE
    centers = [[150.133, 2.302], ['10:00:15.769','+02:15:39.52'], ['10:00:18.394','+02:14:58.78'], ['10:00:23.562','+02:14:34.10']]
    xy = np.array([cat.ra, cat.dec])
    tree = scipy.spatial.cKDTree(xy.T, 10)
    for center in centers:
        if isinstance(center[0],str):
            ra = threedhst.utils.DMS2decimal(center[0], hours=True)
            dec = threedhst.utils.DMS2decimal(center[1], hours=False)
        else:
            ra, dec = center[0], center[1]
        #
        idx = tree.query([ra, dec], k=2)[1][0]
        a = ax.plot(cat.x_image[idx]/fact, cat.y_image[idx]/fact, marker='o', color='white', alpha=0.2, ms=20)
    #
    unicorn.plotting.savefig(fig, root+'_zFOURGE_%5.3f.png' %(ztry))
    
    ### 3D plot, doesn't really work
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    zsel = (zfit.z_max_spec > 0.7) & (zfit.z_max_spec < 1.5)
    ax.scatter(cat.ra[has_gris & zsel], cat.dec[has_gris & zsel], zfit.z_max_spec[has_gris & zsel])
        
    rf = catIO.Readfile('/Users/gbrammer/research/drg/PHOTZ/EAZY/GOODS_F140W/HIGHRES/OUTPUT/goodsn1.7.153-155.rf')
    
    rf_uv = catIO.Readfile('/research/HST/GRISM/3DHST/GOODS-S/FIREWORKS/fireworks.153-155.rf')
    rf_vj = catIO.Readfile('/research/HST/GRISM/3DHST/GOODS-S/FIREWORKS/fireworks.155-161.rf')
    
    uv = (-2.5*np.log10(rf_uv.l153/rf_uv.l155))[zfit.phot_id-1]
    vj = (-2.5*np.log10(rf_vj.l155/rf_vj.l161))[zfit.phot_id-1]
    
    plt.plot(uv[~keep], dz[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,2.5)
    plt.plot(uv[keep], dz[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.01,2.5) 
    xm, ym, ys, N = threedhst.utils.runmed(uv[keep], dz[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(uv[keep], dzphot[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    
    plt.scatter(vj[keep], uv[keep], marker='o', s=sz[keep], alpha=0.4)
    
    plt.scatter(zfit.z_max_spec[keep], uv[keep], marker='o', s=sz[keep], alpha=0.4)

    plt.scatter(zfit.z_max_spec[keep], zfit.mag[keep], marker='o', s=sz[keep], alpha=0.4)

    plt.scatter(zfit.z_max_spec[keep], cat['[24um]_totf'][zfit.phot_id-1][keep], marker='o', s=sz[keep], alpha=0.4)
    plt.ylim(0.01,1000); plt.semilogy()

    plt.scatter(zfit.z_max_spec[keep], cat.Kr50[keep], marker='o', s=sz[keep], alpha=0.4)
    
    plt.plot(np.abs(dz[~keep]), np.abs(dzphot[~keep]), color='red', alpha=0.2)
    plt.plot(np.abs(dz[keep]), np.abs(dzphot[keep]), color='blue', alpha=0.4)
    plt.plot([1.e-6,10],[1e-6,10], color='black', alpha=0.4, linewidth=2, marker='', linestyle='-')
    plt.plot([1.e-6,10],[1e-5,100], color='black', alpha=0.4, linewidth=2, marker='', linestyle='--')
    plt.loglog()
    plt.xlim(1.e-5,8); plt.ylim(1.e-5,8)
    
    delta = np.abs(zfit.z_max_spec - zfit.z_peak_phot)/(1+zfit.z_peak_phot)
    plt.plot(delta[~keep], dz_spec[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.001,1.5) ; plt.semilogx()
    plt.plot(delta[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.001,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(delta[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    keep = keep & (delta < 0.1)
    
    yh, xh, qq = plt.hist(dz_spec[keep], range=(-0.01, 0.01), bins=40, alpha=0.5)
    xx = xh[:-1]/2.+xh[1:]/2.
    import gbb.pymc_gauss as gaussfit
    fit = gaussfit.init_model(xx, yh, np.maximum(np.sqrt(yh),1))
    NSAMP,NBURN=11000,1000
    fit.sample(NSAMP, NBURN)
    trace = fit.trace('eval_gaussian')[::NSAMP/25,:]
    plt.errorbar(xx, yh, np.sqrt(yh))
    for i in range(trace.shape[0]):
        plt.plot(xx, trace[i,:], color='red', alpha=0.2, linestyle='-', marker='')
        
    ###
    test = keep
    test = keep & nopartial
    plt.plot(zfit.z_spec[test], dz[test], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[test], dz[test], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[test], dzphot[test], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    
    bad = keep & nopartial & (np.abs(dz) > 0.02)
    fp = open('/tmp/bad','w')
    for id in zfit.id[bad]:
        fp.write('/tmp/%s.zfit.png\n' %(id))
    fp.close()
    
    #### Local environmental density
    width = 0.02
    N = 7
    r0, d0 = np.median(cat.ra), np.median(cat.dec)
    dra = (cat.ra-r0)*3600.*np.cos(d0/360*2*np.pi)
    ddec = (cat.dec-d0)*3600.
    
    idx = np.arange(len(cat.ra))
    
    nearest = np.zeros((cat.ra.shape[0], N))
    for i in idx[has_gris]:
        ztry = zfit.z_max_spec[i]
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        xy = np.array([dra[bin], ddec[bin]])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        dist, did = tree.query([dra[i], ddec[i]], k=N+1)
        nearest[i,:] = dist[1:]
    
    zbin = (zfit.z_max_spec > 0.9) & (zfit.z_max_spec < 1.4) & (fout.lmass > 10.5)
    halpha = lines[1].data['HALPHA_FLUX']*1.
    halpha[halpha == 0] += 0.01 * (1+0.05*np.random.normal(size=(halpha == 0).sum()))
    plt.plot(1./nearest[zbin, 0]**2, halpha[zbin], alpha=0.2)
    
    #### Pairs
    min_dist = 1.5
    pairs = (zfit.z_max_spec > 0.8) & (zfit.z_max_spec < 2.5) & (fout.lmass > 10.) & (nearest[:,0] < min_dist)
    detect = pyfits.open('COSMOS_v2.0_PHOTOMETRY/Detection/COSMOS_F125W-F140W-F160W_detect.fits')
    
    ii = idx[pairs]
    ii = ii[np.argsort(zfit.z_max_spec[ii])]
    
    DX = 5/0.06
    for i in ii:
        ztry = zfit.z_max_spec[i]
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        xy = np.array([dra[bin], ddec[bin]])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        dist, did = tree.query([dra[i], ddec[i]], k=N+1)
        neighbors = did[(dist > 0) & (dist < min_dist)]
        #
        xc, yc = int(np.round(cat.x_image[i])), int(np.round(cat.y_image[i]))
        subim = detect[0].data[yc-DX:yc+DX, xc-DX:xc+DX]
        plt.imshow(0-subim, vmin=-10, vmax=0.03, interpolation='Nearest')
        plt.text(DX, DX+5, '%5.3f' %(ztry), horizontalalignment='center', verticalalignment='bottom', color='red')
        for neighbor in neighbors:
            plt.plot(cat.x_image[bin][neighbor]-xc+DX, cat.y_image[bin][neighbor]-yc+DX, marker='o', ms=12, color='yellow', alpha=0.3)
            plt.text(cat.x_image[bin][neighbor]-xc+DX, cat.y_image[bin][neighbor]-yc+DX+5, '%5.3f' %(zfit.z_max_spec[bin][neighbor]), horizontalalignment='center', verticalalignment='bottom', color='green')
        #
        plt.xlim(0,2*DX)
        plt.ylim(0,2*DX)
        
def interlace_goodss():
    """
    Reduce the GOODS-S pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F140W', REFERENCE = '/3DHST/Photometry/Work/GOODS-S/v2/images/GOODS-S_F140W_sci_sub.fits', SEGM = '/3DHST/Photometry/Work/GOODS-S/v2/sextr/checkimages/GOODS-S_125_140_160_v2.seg.fits')
         
    NGROW=125
    pad=60
    #CATALOG='phot/GOODS-S_F160W_v1.cat'
    CATALOG='/3DHST/Photometry/Work/GOODS-S/v2/sextr/catalogs/GOODS-S_F140W_conv_v2.cat'

    direct=glob.glob('GOODS-S-*-F140W_asn.fits')
    direct.remove('GOODS-S-1-F140W_asn.fits')
    direct.remove('GOODS-S-28-F140W_asn.fits')
        
    extract_limit = 24
    skip_completed=False
    REF_ROOT='GOODS-S_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        #
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-1[0-9]-G141_inter.fits')
    inter.remove('GOODS-S-14-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-2')

    skip_completed = True

    models = glob.glob('GOODS-S-3[0-9]_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where(model.cat.mag < 24.)
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    models = glob.glob('GOODS-S-3[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

    ### Get some quality flags on the reduced spectra
    root='GOODS-S'

    files = glob.glob('GOODS-S-*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	  
    

def interlace_hudf():
    """
    Reduce the GOODS-S pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    #os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE')
    os.chdir(unicorn.GRISM_HOME+'HUDF/INTERLACE')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    #unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F160W', REFERENCE = 'phot/GOODS-S_F160W_sci.fits', SEGM = 'phot/GOODS-S_F160W_v1.seg.fits')
    #unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F160W', REFERENCE = 'phot/HUDF_F160W_detection.fits', SEGM = 'phot/HUDF_F160W_seg.fits')
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F140W', REFERENCE = 'phot/UDF-F140W_sci_drz.fits', SEGM = 'phot/HUDF_F160W_seg.fits')
        
    NGROW=125
    pad=60
    #CATALOG='phot/GOODS-S_F160W_v1.cat'
    CATALOG='phot/convolved_h_nomean.cat'

    direct=glob.glob('GOODS-S-3[4678]-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='GOODS-S_F140W'

#    ##### Generate the interlaced images, including the "blotted" detection image
#    for i in range(len(direct)):
#        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
#        #
#        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
#        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
#        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
#        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-3[4678]*-G141_inter.fits')
    redo = False

    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]

        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            #model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
            #model.extract_spectra_and_diagnostics(MAG_LIMIT=24)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-34')

    skip_completed = True

    models = glob.glob('GOODS-S-36_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        #
        ii = np.where(model.cat.mag < 24.)

        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    models = glob.glob('GOODS-S-34_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

    ### Get some quality flags on the reduced spectra
    root='HUDF'

    files = glob.glob('GOODS-S-3[4678]*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	
            
def interlace_uds18():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    os.chdir(unicorn.GRISM_HOME+'UDS/UDS-18')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='UDS-18_F140W', filter='F140W', REFERENCE = '/3DHST/Photometry/Work/UDS/UDS-18/psfmatch/run/UDS-18_F140W_sci_new/UDS_18_F140W_sci_conv.fits', SEGM = '/3DHST/Photometry/Work/UDS/UDS-18/sextr/checkimages/UDS-18_F140W.seg.fits')

        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = '/3DHST/Photometry/Work/UDS/UDS-18/sextr/catalogs/UDS-18_F140W.cat'
    direct=glob.glob('UDS-18-F140W_asn.fits')

    extract_limit = 25.8
    skip_completed=False
    REF_ROOT='UDS-18_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('UDS-18-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-18')

    skip_completed = True

    models = glob.glob('UDS-18_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('UDS-18_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)	

    ### Get some quality flags on the reduced spectra
    root='UDS-18'

    files = glob.glob('UDS-18*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	  	

def interlace_sntile41():
    
    #
    ### SWarp for color image
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F160W.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_wht.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F160W_wht.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F125W.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/COSMOS/ACS/*sci.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F814W.fits')

    ### F160W
    sub_r = pyfits.open('match_F160W.fits')[0].data*10**(-0.4*(25.96-25.96))
    ### F125W
    sub_g = pyfits.open('match_F125W.fits')[0].data*10**(-0.4*(26.25-25.96))
    ### F814W
    sub_b = pyfits.open('match_F814W.fits')[0].data*10**(-0.4*(25.94-25.96))*1.5

    Q, alpha, m0 = 4, 12, -0.005
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='TILE41.png', shape=sub_r.shape)
    Q, alpha, m0 = 3.5, 5, -0.01
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='TILE41_2.png', shape=sub_r.shape)
    
    unicorn.reduce.interlace_combine('TILE41-132-G141', NGROW=125)
    unicorn.reduce.interlace_combine('TILE41-132-F160W', NGROW=125, auto_offsets=True)
    
    model = unicorn.reduce.process_GrismModel('TILE41-132', MAG_LIMIT=25)
        
    se = threedhst.sex.SExtractor()
    
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    
    ## XXX add test for user-defined .conv file
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = 'test.cat'
    se.options['CHECKIMAGE_NAME'] = 'test_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'match_F160W_wht.fits[0]'
    se.options['FILTER']    = 'Y'
    threedhst.sex.USE_CONVFILE =  'gauss_2.0_5x5.conv'
    threedhst.sex.USE_CONVFILE =  'gauss_4.0_7x7.conv'
    se.copyConvFile()
    se.options['FILTER_NAME'] = threedhst.sex.USE_CONVFILE
    
    #### Detect thresholds (default = 1.5)
    se.options['DEBLEND_NTHRESH']    = '16' 
    se.options['DEBLEND_MINCONT']    = '0.001' 
    se.options['DETECT_MINAREA']    = '36' 
    se.options['DETECT_THRESH']    = '1.5' 
    se.options['ANALYSIS_THRESH']  = '1.5'
    se.options['MAG_ZEROPOINT'] = '%.2f' %(model.ZP)
    
    #### Run SExtractor
    status = se.sextractImage('match_F160W.fits[0]')
    threedhst.sex.sexcatRegions('test.cat', 'test.reg', format=2)
    
    unicorn.reduce.prepare_blot_reference(REF_ROOT='TILE41_F160W', filter='F160W', REFERENCE = 'match_F160W.fits', SEGM = 'test_seg.fits', Force=False, sci_extension=0)
    
    unicorn.reduce.blot_from_reference(REF_ROOT='TILE41_F160W', DRZ_ROOT = 'TILE41-132-F160W', NGROW=125, verbose=True)
    unicorn.reduce.interlace_combine_blot(root='TILE41-132-F160W', view=True, pad=60, REF_ROOT = 'TILE41_F160W', CATALOG='test.cat',  NGROW=125, verbose=True, growx=2, growy=2, auto_offsets=True, NSEGPIX=0)
    #unicorn.reduce.fill_inter_zero('TILE41-132_inter_seg.fits')
    unicorn.reduce.fill_inter_zero('TILE41-132_ref_inter.fits')
    
    unicorn.reduce.interlace_combine('TILE41-132-G141', NGROW=125, auto_offsets=True)
    unicorn.reduce.interlace_combine('TILE41-132-F160W', NGROW=125, auto_offsets=True)
    unicorn.reduce.fill_inter_zero('TILE41-132-F160W_inter.fits')
    
    model = unicorn.reduce.process_GrismModel('TILE41-132', MAG_LIMIT=25, REFINE_MAG_LIMIT=23, make_zeroth_model=False)
    
    unicorn.reduce.blot_from_reference(REF_ROOT='TILE41_F160W', DRZ_ROOT = 'TILE41-152-F140W', NGROW=125, verbose=True)
    unicorn.reduce.interlace_combine_blot(root='TILE41-152-F140W', view=True, pad=60, REF_ROOT = 'TILE41_F160W', CATALOG='test.cat',  NGROW=125, verbose=True, growx=2, growy=2, auto_offsets=True, NSEGPIX=0)
    
    unicorn.reduce.fill_inter_zero('TILE41-152_ref_inter.fits')
    
    unicorn.reduce.interlace_combine('TILE41-152-G141', NGROW=125, auto_offsets=True)
    unicorn.reduce.interlace_combine('TILE41-152-F140W', NGROW=125, auto_offsets=True)
    unicorn.reduce.fill_inter_zero('TILE41-152-F140W_inter.fits')
    
    model = unicorn.reduce.process_GrismModel('TILE41-152', MAG_LIMIT=25, REFINE_MAG_LIMIT=23)
    
    #### Seems OK, but haven't been able to check against z_specs
    
    ### Extract all objects > z>0.6, m140 > X
    cat, zout, fout = unicorn.analysis.read_catalogs('TILE41')
    matcher = catIO.CoordinateMatcher(cat)
    
    dr, idx = matcher.match_list(ra=model.ra_wcs, dec=model.dec_wcs)
    
    ok = (dr < 0.8) & ((zout.z_peak[idx] > 0.65) & (model.cat.mag < 23)) | ((zout.z_peak[idx] > 1.2) & (model.cat.mag < 24))
    
    for id in model.objects[ok]:
        try:
            model.twod_spectrum(id)
            gris = unicorn.interlace_fit.GrismSpectrumFit('TILE41-152_%05d' %(id))
            gris.fit_in_steps(zrfirst=(0.5, 4))
            #os.system('open TILE41-132_%05d*zfit*png' %(id))
        except:
            pass
            
def interlace_combine_blot(root='COSMOS-19-F140W', view=True, pad=60, REF_ROOT = 'COSMOS_F160W', CATALOG='UCSC/catalogs/COSMOS_F160W_v1.cat',  NGROW=125, verbose=True, growx=2, growy=2, auto_offsets=False, NSEGPIX=8):
    """
    Combine blotted image from the detection mosaic as if they were 
    interlaced FLT images
    
    xxxxx Still not modified for differing factors in x and y!
    """
    from pyraf import iraf

    import threedhst.prep_flt_files
    import unicorn
    
    if unicorn.hostname().startswith('uni'):
        view = False
        
    pointing = root.split('-F1')[0]
    
    #unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = root, NGROW=NGROW, verbose=verbose)
    
    if view:
        import threedhst.dq
        ds9 = threedhst.dq.myDS9()

    run = threedhst.prep_flt_files.MultidrizzleRun(root)
    #run.blot_back(ii=0, copy_new=True)
    
    #im = pyfits.open(os.getenv('iref')+'ir_wfc3_map.fits')
    #PAM = im[1].data
    #im.close()
    
    scl = np.float(run.scl)
    xsh, ysh = threedhst.utils.xyrot(np.array(run.xsh)*scl, np.array(run.ysh)*scl, run.rot[0])
    yi,xi = np.indices((1014,1014+2*NGROW))

    #N = np.zeros((2028+pad+2*2*NGROW, 2028+pad+2*2*NGROW))
    inter_sci = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW))
    #inter_err = np.zeros((2028+pad+2*2*NGROW, 2028+pad+2*2*NGROW))
    inter_seg = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW), dtype=np.int32)
    inter_N = np.zeros((1014*growy+pad+growy*2*NGROW, 1014*growx+pad+growx*2*NGROW), dtype=np.int32)
    
    xi+=pad/(2*growx)
    yi+=pad/(2*growy)+NGROW
    
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    
    #### These were needed for COSMOS-19
    #xi -= 1
    #yi += 1
    
    #### From pieter
    dxs = np.array([0,-20,-13,7]) + np.int(np.round(xsh[0]))*0
    dys = np.array([0,-7,-20,-13]) + np.int(np.round(ysh[0]))*0
    
    #### GOODS-N field from Weiner program 11600
    if root.startswith('GOODS-N') | root.startswith('GNGRISM'):
        dxs = np.array([0,-9,-4,5]) + np.int(np.round(xsh[0]))*0
        dys = np.array([0,-3,-11,-8]) + np.int(np.round(ysh[0]))*0
    
    if root.startswith('TILE41'):
        dxs = np.array([-4, 9, -8, 5, 5, -8]) + np.int(np.round(xsh[0]))*0
        dys = np.array([-5, -4, 4, 5, 5, 4]) + np.int(np.round(ysh[0]))*0
        
    #### Erb QSO sightlines
    if flt[0].header['PROPOSID'] == 12471:
        dxs, dys = np.array([  0, -10])*growx, np.array([ 0, -7])
    
    if auto_offsets:
        xinter, yinter = unicorn.reduce.get_interlace_offsets(root+'_asn.fits', growx=growx, growy=growy)
        dxs = xinter + np.int(np.round(xsh[0]))*0
        dys = yinter + np.int(np.round(ysh[0]))*0
        print 'Auto offsets: ', dxs, dys
        
    dxi = np.cast[np.int32](np.ceil(dxs/growx))
    dyi = np.cast[np.int32](np.ceil(dys/growy))
    
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
        
        #### Use the pixel area map correction
        #im[1].data *= PAM
        #### Divide by 4 to conserve surface brightness with smaller output pixels
        im[0].data /= (growx*growy)
        #im_wht[0].data /= 4
        #im[2].data /= 4
        
        ### Mask cosmic rays
        if i == 0:
            im_flt = pyfits.open(flt+'.fits')
            h0 = im_flt[0].header
            h1 = im_flt[1].header
            header = unicorn.reduce.scale_header_wcs(h1.copy(), factor=2, growx=growx, growy=growy, pad=pad, NGROW=NGROW)
            
            #header_wht = header.copy()
            #header_wht.update('EXTNAME','ERR')
            header.update('EXTNAME','SCI')
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
        
        inter_sci[yi*growy+dy,xi*growx+dx] += im[0].data
        #inter_err[yi*2+dy,xi*2+dx] += im_wht[0].data
        inter_seg[yi*growy+dy,xi*growx+dx] += im_seg[0].data
        inter_N[yi*growy+dy,xi*growx+dx] += 1
        
        #
        if view:
            ds9.view_array(inter_sci, header=header)
            ds9.scale(-0.1,5)
    
    inter_N[inter_N == 0] = 1
    
    #### Write interlaced sci/wht image
    header.update('NGROW', NGROW, comment='Number of pixels added to X-axis (centered)')        
    header.update('PAD', pad, comment='Additional padding around the edges') 
    
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=inter_sci/inter_N, header=header)
    #wht = pyfits.ImageHDU(data=inter_err, header=header_wht)
            
    image = pyfits.HDUList([hdu, sci]) #, wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
    
    image.writeto(pointing+'_ref_inter.fits', clobber=True, output_verify="fix")
    
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
        
    hdu = pyfits.PrimaryHDU(header=header, data=np.cast[np.int32](inter_seg/inter_N))
    hdu.writeto(pointing+'_inter_seg.fits', clobber=True)
    
    if verbose:
        print 'Clean up segmentation image...'
        
    #### Clean up overlap region of segmentation 
    unicorn.reduce.fill_inter_zero(pointing+'_inter_seg.fits')
    im = pyfits.open(pointing+'_inter_seg.fits', mode='update')
    #s = np.ones((3,3))
    #labeled_array, num_features = nd.label(im[0].data, structure=s)
    
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[1]'

    #### Maximum filter to flag overlap regions
    max_filter = nd.maximum_filter(im[0].data, size=(3,3))
    bad = (max_filter - im[0].data) > 2
    inter_seg = max_filter*1
    inter_seg[bad] = 0
    
    if verbose:
        print unicorn.noNewLine+'Clean up segmentation image...[2]'

    #### Orphan pixels in the max-filtered image
    npix = nd.convolve((inter_seg > 0)*1, np.ones((3,3)))
    bad = (inter_seg > 0) & (npix < 4)
    inter_seg[bad] = 0
    im[0].data = inter_seg
    im.flush()
    
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
    wcs_this = pywcs.WCS(header=this_drz_head)
    
    fp = open("/tmp/%s.drz_xy" %(root),'w')
    #fpr = open("/tmp/%s.flt_xy.reg" %(root),'w')
    NOBJ = len(old_cat.id)
    for i in range(NOBJ):
        #xw, yw = wcs_this.rd2xy((old_cat.ra[i], old_cat.dec[i]))
        xw, yw = wcs_this.wcs_sky2pix([[old_cat.ra[i], old_cat.dec[i]]],1)[0]
        fp.write('%.2f %.2f\n' %(np.clip(xw, -200, 4999), np.clip(yw, -200, 4999)))
    
    fp.close()
    
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
        
    os.remove("/tmp/%s.drz_xy" %(root))
    
    print unicorn.noNewLine+'Convert from reference coords to the distorted flt frame with iraf.wblot - Done.'
    
    flt_x, flt_y, popID = [], [], []
    fpr = open("%s_inter.reg" %(pointing),'w')
    fpr.write("global color=green dashlist=8 3 width=1 font=\"helvetica 8 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nimage\n")
    for i,line in enumerate(status[-NOBJ:]):
        spl = np.cast[float](line.split())
        fxi, fyi = (spl[0]+NGROW)*growx+pad/2-1, (spl[1]+NGROW)*growy+pad/2-1
        if (fxi > pad/2) & (fxi < ((1014+2*NGROW)*growx+pad/2.)) & (fyi > (pad/2.+NGROW*2)) & (fyi < (pad/2.+NGROW*growy+1014*growy)):
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
        print 'Make the catalog file: %s' %(pointing+'_inter.cat')
        
    ### Put the x,y coordinates in the FLT frame into the SExtractor catalog
    status = old_cat.addColumn(np.array(flt_x), name='X_FLT', comment='X pixel in FLT frame', verbose=True)
    status = old_cat.addColumn(np.array(flt_y), name='Y_FLT', comment='Y pixel in FLT frame', verbose=True)
    old_cat.write(outfile=pointing+'_inter.cat')
        
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
    
    import unicorn.reduce
    
    tmpname = threedhst.utils.gen_tempname()
            
    im_ref = pyfits.open(REF_ROOT+'_ref.fits')
    
    if not os.path.exists(DRZ_ROOT+'.run'):
        run_multidrizzle=True
    else:
        run = threedhst.prep_flt_files.MultidrizzleRun(DRZ_ROOT)
        if (run.outnx, run.outny) != (im_ref[0].header['NAXIS1'], im_ref[0].header['NAXIS2']):
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
        ra0, dec0 = ref_wcs.wcs_pix2sky([[im_ref[0].header['NAXIS1']/2., im_ref[0].header['NAXIS2']/2.]],1)[0]
        
        threedhst.prep_flt_files.startMultidrizzle(DRZ_ROOT+'_asn.fits', use_shiftfile=True, skysub=False, final_scale=0.06, pixfrac=0.8, driz_cr=False, updatewcs=False, clean=True, median=False, build_drz=False, ra=ra0, dec=dec0, final_outnx=im_ref[0].header['NAXIS1'], final_outny=im_ref[0].header['NAXIS2'], final_rot=-np.arctan(im_ref[0].header['CD1_2']/im_ref[0].header['CD1_1'])/2/np.pi*360, generate_run=True)
        
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
            
            xoff, yoff, roff = 0., 0., 0.
            #for i in range(3):
            for geom in ['rotate', 'shift', 'rotate', 'shift', 'shift', 'shift']:
              #
                threedhst.process_grism.flprMulti()
                if os.path.exists('align_blot.fits'):
                    os.remove('align_blot.fits')
                #
                status = iraf.wblot( data = REF_ROOT+'_ref.fits', outdata = 'align_blot.fits', 
                   outnx = 1014, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
                   coeffs = run.flt[idx]+'_coeffs1.dat', xgeoim = '', ygeoim = '', 
                   align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
                   rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
                   refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
                   decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1'], 
                   yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
                   dr2gpar = '', expkey = 'exptime', expout = 'input', 
                   in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)
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
        status = iraf.wblot( data = REF_ROOT+'_ref.fits', outdata = FLT.replace('_flt','_blot'), 
           outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
           coeffs = tmpname+'_coeffs.dat', xgeoim = '', ygeoim = '', 
           align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
           rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
           refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
           decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
           yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
           dr2gpar = '', expkey = 'exptime', expout = 'input', 
           in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)
                
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
        
        threedhst.process_grism.flprMulti()
        status = iraf.wblot( data = REF_ROOT+'_seg.fits', outdata = FLT.replace('_flt','_seg'), 
           outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
           coeffs = tmpname+'_coeffs.dat', lambd = 1392.0, xgeoim = '', ygeoim = '', 
           align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
           rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
           refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
           decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
           yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
           dr2gpar = '', expkey = 'exptime', expout = 'input', 
           in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)
        
        #
        if verbose:
            print unicorn.noNewLine+'sci, weight, seg0, seg1'

        if os.path.exists(FLT.replace('_flt','_ones')):
            os.remove(FLT.replace('_flt','_ones'))
        
        threedhst.process_grism.flprMulti()
        status = iraf.wblot( data = REF_ROOT+'_ones.fits', outdata = FLT.replace('_flt','_ones'), 
           outnx = 1014+2*NGROW, outny = 1014, geomode = 'user', interpol = 'poly5', sinscl = 1.0, 
           coeffs = tmpname+'_coeffs.dat', lambd = 1392.0, xgeoim = '', ygeoim = '', 
           align = 'center', scale = run.scl, xsh = run.xsh[idx]+xoff, ysh = run.ysh[idx]+yoff, 
           rot = run.rot[idx]+roff, shft_un = 'input', shft_fr = 'input', 
           refimage = '', outscl = 0.128, raref = im_flt[1].header['CRVAL1'], 
           decref = im_flt[1].header['CRVAL2'], xrefpix = im_flt[1].header['CRPIX1']+NGROW, 
           yrefpix = im_flt[1].header['CRPIX2'], orient = im_flt[1].header['ORIENTAT'], 
           dr2gpar = '', expkey = 'exptime', expout = 'input', 
           in_un = 'cps', out_un = 'cps', fillval = 0.0, mode = 'al', Stdout=1)
        
        #### Add NGROW to header
        im_seg = pyfits.open(FLT.replace('_flt','_seg'), mode='update')
        ones = pyfits.open(FLT.replace('_flt','_ones'))
        yh, xh = np.histogram(ones[0].data.flatten(), range=(0.1,ones[0].data.max()), bins=100)
        keep = ones[0].data > (0.5*xh[:-1][yh == yh.max()])
        
        ratio = im_seg[0].data / ones[0].data
        #test = (np.abs(np.log(ratio/np.round(ratio))) < 1.e-5) & keep
        
        # test = (np.abs(ratio-np.round(ratio)) < 0.01) & keep
        # ratio[test] = np.round(ratio[test])
        # ratio[~test | (ones[0].data == 0)] = 0
        
        # im_seg[0].data[keep] /= ones[0].data[keep]
        # im_seg[0].data[~keep] = 0
        # im_seg[0].data = np.cast[int](np.round(im_seg[0].data))
        im_seg[0].data = np.cast[np.int32](np.round(ratio))
        
        im_seg[0].header.update('NGROW',NGROW, comment='Number of pixels added to X-axis (centered)')
        im_seg[0].header.update('FILTER', REF_FILTER)
        im_seg[0].header.update('EXPTIME', REF_EXPTIME)
        
        if not pyfits.__version__.startswith('3'):
            pyfits.writeto(im_seg.filename(), data=im_seg[0].data, header=im_seg[0].header, clobber='True')
        else:
            im_seg.flush()
    
        os.remove(FLT.replace('_flt','_ones'))
    
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
def fill_inter_zero(image='TILE41-132-F160W_inter.fits'):
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
    
    if HAS_ERR:
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
    wave = unicorn.reduce.sens_files['A'].field('WAVELENGTH')
    flam = np.abs(stats[0])*4/(unicorn.reduce.sens_files['A'].field('SENSITIVITY')*46.5)
    fnu = flam*wave**2/3.e18
    mag = -2.5*np.log10(fnu)-48.6
    ax.plot(wave, mag)
    ax.set_xlim(1.05e4,1.7e4)
    ran = (wave > 1.05e4) & (wave < 1.7e4)
    ax.set_ylim(mag[ran].min()-0.2, mag[ran].max()+0.2)
    ax.set_xlabel(r'$\lambda$') 
    ax.set_ylabel(r'mag of offset')
    
    unicorn.catalogs.savefig(fig, root+'_inter_deep_residuals.png')
    
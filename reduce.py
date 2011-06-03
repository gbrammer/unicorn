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

import numpy as np
import matplotlib.pyplot as plt
USE_PLOT_GUI = False

import pyfits
import pyraf
from pyraf import iraf

import threedhst

def prep():
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    
    os.chdir("/research/HST/GRISM/3DHST/COSMOS/OWN_REDUCTION")
    
    ALIGN = '../NMBS/COSMOS-1.V4.K_nosky.fits'
    
    threedhst.prep_flt_files.process_3dhst_pair('ibhm46030_asn.fits', 'ibhm46040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], GET_SHIFT=True, SKIP_DIRECT=True)
    
    # threedhst.prep_flt_files.startMultidrizzle('COSMOS-1-F140W_asn.fits',
    #     use_shiftfile=True, skysub=False,
    #     final_scale=0.128254, pixfrac=1.0, driz_cr=True,
    #     updatewcs=True, clean=True, median=True)
    
    # threedhst.prep_flt_files.startMultidrizzle('COSMOS-15-F140W_asn.fits',
    #     use_shiftfile=True, skysub=False,
    #     final_scale=0.06, pixfrac=0.8, driz_cr=False,
    #     updatewcs=False, clean=True, median=False)

#
def combine(root='COSMOS-1-F140W'):
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    ds9 = threedhst.dq.myDS9()

    #os.chdir("/research/HST/GRISM/3DHST/COSMOS/OWN_REDUCTION")

    run = threedhst.prep_flt_files.MultidrizzleRun(root)
    #run.blot_back(ii=0, copy_new=True)
    
    scl = np.float(run.scl)
    xsh, ysh = threedhst.utils.xyrot(np.array(run.xsh)*scl, np.array(run.ysh)*scl, run.rot[0])
    yi,xi = np.indices((1014,1014))
    pad = 60
    N = np.zeros((2028+pad,2028+pad))
    inter = np.zeros((2028+pad,2028+pad))
    xi+=pad/4
    yi+=pad/4
    
    #### From pieter
    dxs = np.array([0,-20,-13,7]) + np.int(np.round(xsh[0]))*0
    dys = np.array([0,-7,-20,-13]) + np.int(np.round(ysh[0]))*0
    
    dxi = np.cast[int](np.ceil(dxs/2))
    dyi = np.cast[int](np.ceil(dys/2))
    
    #### Find hot pixels, which are flagged as cosmic 
    #### rays at the same location in each of the flt files.  Ignore others.
    hot_pix = np.zeros((1014,1014),dtype='int')
    for flt in run.flt:
        im = pyfits.open(flt+'.fits')
        hot_pix += (im[3].data & 4096) / 4096
    
    hot_pix = hot_pix > 2
    
    # #### Background
    # bgsum = np.zeros((1014+pad/2,1014+pad/2))
    # bgN = bgsum*0.
    # pfl = pyfits.open('/research/HST/GRISM/IREF/uc721143i_pfl.fits')
    # 
    # for i,flt in enumerate(run.flt):
    #     print flt
    #     flt = run.flt[i]
    #     im = pyfits.open(flt+'.fits')
    #     #im[1].data /= pfl[1].data[5:-5,5:-5]
    #     use = ((im[3].data & 4) == 0) & (~hot_pix) & ((im[3].data & 32) == 0)
    #     bgsum[yi[use]+dyi[i],xi[use]+dxi[i]] += im[1].data[use]
    #     bgN[yi[use]+dyi[i],xi[use]+dxi[i]] += 1
    # 
    # bgsum /= bgN
    # 
    # i=0
    # diff = bgsum
    # flt = run.flt[i]
    # im = pyfits.open(flt+'.fits','update')
    # #im[1].data /= pfl[1].data[5:-5,5:-5]
    # use = ((im[3].data & 4) == 0) & (~hot_pix) & ((im[3].data & 32) == 0)
    # im[1].data[use] -= diff[yi[use]+dyi[i],xi[use]+dxi[i]]
    # im[1].data[~use] = -1
    # sub = np.abs(im[1].data) < 0.2
    # im[1].data -= threedhst.utils.biweight(im[1].data[sub], mean=True)
    # im.flush()
    # ds9.view(im[1].data)
    # ds9.scale(-0.1,1)
    
    for i,flt in enumerate(run.flt):
        print flt
        flt = run.flt[i]
        im = pyfits.open(flt+'.fits')
        if i == 0:
            h0 = im[0].header
            h1 = im[1].header
            header = red.scale_header_wcs(h1.copy(), factor=2, pad=pad)
        dx = np.int(np.round((xsh[i]-xsh[0])*2))
        dy = np.int(np.round((ysh[i]-ysh[0])*2))
        #
        dx = dxs[i]
        dy = dys[i]
        #
        #use = ((im[3].data & 4096) == 0) & ((im[3].data & 4) == 0) #& (xi > np.abs(dx/2)) & (xi < (1014-np.abs(dx/2))) & (yi > np.abs(dy/2)) & (yi < (1014-np.abs(dy/2)))
        #
        use = ((im[3].data & 4) == 0) & (~hot_pix) & ((im[3].data & 32) == 0)
        inter[yi[use]*2+dy,xi[use]*2+dx] += im[1].data[use]
        N[yi[use]*2+dy,xi[use]*2+dx] += 1
        #
        ds9.view_array(inter, header=header)
        ds9.scale(-0.1,5)
            
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=inter, header=header)
    wht = pyfits.ImageHDU(data=N, header=header)
    pyfits.HDUList([hdu,sci,wht]).writeto(root+'_inter.fits', clobber=True)
    
    try:
        os.remove('sci.fits')
        os.remove('wht.fits')
        os.remove('scix.fits')
    except:
        pass
    #
    new_coeffs_dat(input=run.flt[0]+'_coeffs1.dat', output='scale_coeffs.dat')
    
    status = iraf.wdrizzle(data = root+'_inter.fits[1]', outdata = root+'_inter_sci.fits', \
     outweig = root+'_inter_wht.fits', outcont = "", in_mask = root+'_inter.fits[2]', 
     wt_scl = 'exptime', \
     outnx = 2300, outny = 2300, geomode = 'wcs', kernel = 'square', \
     pixfrac = 1.0, coeffs = "scale_coeffs.dat", lamb = 1392., xgeoim = "", ygeoim = "", \
     align = 'center', scale = 1.0, xsh = 0.0, ysh = 0.0, rot = 0.0, \
     shft_un = 'output', shft_fr = 'output', outscl = 0.060, \
     raref = header['CRVAL1'], decref = header['CRVAL2'], xrefpix = 1150, 
     yrefpix = 1150,
     orient = 0.0, dr2gpar = "", expkey = 'exptime', in_un = 'cps', \
     out_un = 'cps', fillval = '0', mode = 'al')

def new_coeffs_dat(input='ibhm29wlq_flt_coeffs1.dat',output='scale_coeffs.dat', factor=2, pad=60):
    lines = open(input).readlines()

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
    
def scale_header_wcs(header, factor=2, pad=60):
    """
    Take an FLT header but scale the WCS keywords to account for the images 
    expanded/contracted by a factor of 'factor'.
    """
    import numpy as np
    
    header.update('NAXIS1',header['NAXIS1']*factor+pad)
    header.update('NAXIS2',header['NAXIS2']*factor+pad)

    header.update('CRPIX1',header['CRPIX1']*factor+pad/2)
    header.update('CRPIX2',header['CRPIX2']*factor+pad/2)
    
    keys = ['CD1_1','CD1_2','CD2_1','CD2_2']
    
    # header.update('ONAXIS1',header['ONAXIS1']*factor+pad)
    # header.update('ONAXIS2',header['ONAXIS2']*factor+pad)    
    # header.update('OCRPIX1',header['OCRPIX1']*factor+pad/2)
    # header.update('OCRPIX2',header['OCRPIX2']*factor+pad/2)
    # keys.extend(['IDCSCALE','OCX10','OCX11','OCY10','OCY11', 'OCD1_1', 'OCD1_2', 'OCD2_1', 'OCD2_2'])
    
    for key in keys:
        header.update(key,header[key]/factor)
    
    header.update('IDCTAB','')
    
    #### Scale distortion keywords, since X = (x-CRPIX1)/factor
    # for key in header.keys():
    #     if key.startswith('A_') | key.startswith('B_'):
    #         try:
    #             coeffs = np.cast[int](key.split('_')[1:])
    #             header.update(key, header[key] / factor**(np.sum(coeffs)))
    #             #print key, np.sum(coeffs)
    #         except:
    #             pass
                
    header.update('EXPTIME',1)
    
    return header

def grism_model(xc_full=244, yc_full=1244, lam_spec=None, flux_spec=None, grow_factor=2, pad = 60, BEAMS=['A','B','C','D','E'], dydx=True):
    
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    #ds9 = threedhst.dq.myDS9()

    #os.chdir("/research/HST/GRISM/3DHST/COSMOS/OWN_REDUCTION")

    conf = threedhst.process_grism.Conf('WFC3.IR.G141.V2.0.conf', path='../CONF/').params
    
    ## A star in COSMOS-15-F140W
    ##### For interlaced images
    
    # xc_full = 244.13
    # yc_full = 1244.323
    
    xc = np.int(np.round((xc_full - pad/2)/grow_factor))
    yc = np.int(np.round((yc_full - pad/2)/grow_factor))
    
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
    
    #### Can do everything by pixel, perhaps make an init loop that calculates a 
    #### spectrum for a normalized flux at each pixel and stores in a lookup table
    
    #### Will need to interpolate for the DYDX offset
        
    ##### Configuration parameters into arrays, eventually loop for all BEAMs
    xmi, xma = -580, 730
    NX,NY = xma-xmi, 20
    
    xmi *= grow_factor
    xma *= grow_factor
    NX *= grow_factor
    NY *= grow_factor
    
    model = np.zeros((NY*2+1, NX), dtype=np.float)
    wavelength = model*0
    full_sens = model*0
    
    yi, xi = np.indices(model.shape)
    xi += xmi
    yi -= NY
    xarr, xpix, yarr = np.arange(xma-xmi)+xmi, np.arange(xma-xmi), np.arange(NY*2+1)-NY
    
    for beam in BEAMS:
       
        xmi_beam, xma_beam = tuple(np.cast[int](conf['BEAM'+beam].split())) 
        xmi_beam *= grow_factor
        xma_beam *= grow_factor
        
        dydx_order = np.int(conf['DYDX_ORDER_'+beam])
        dydx_0 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_0']) * grow_factor
        dydx_1 = field_dependent(bigX, bigY, conf['DYDX_'+beam+'_1']) #/ grow_factor        
                    
        xoff_beam = field_dependent(bigX, bigY, conf['XOFF_'+beam]) * grow_factor
        yoff_beam = field_dependent(bigX, bigY, conf['YOFF_'+beam]) * grow_factor
                
        disp_order = np.int(conf['DISP_ORDER_'+beam])
        dldp_0 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_0'])
        dldp_1 = field_dependent(bigX, bigY, conf['DLDP_'+beam+'_1']) / grow_factor
        
        #### Improve alignment of zeroth order
        if beam == 'BXr':
            #dydx_1 = 0.0
            dydx_1 = 0.1
            dydx_0 -= 2/grow_factor
            dydx_0 += dydx_1 * 192 * grow_factor
            
            dldp_x = dldp_1*1.
            f = 0.9*(1+np.abs(bigY-507)/507.*0.2)
            dldp_1 *= f
            dldp_0 += (1-f)*dldp_1*-213*2*(2-(1.02-0.02*np.abs(bigY-507)/507.))
            
        #print 'BEAM_%s: %5.2f %5.2f, %6.3f %5.3f, %9.2f %6.2f' %(beam, xoff_beam / grow_factor, yoff_beam / grow_factor, dydx_0, dydx_1*grow_factor, dldp_0, dldp_1*grow_factor)

        #### Wavelength
        lam = dldp_0 + dldp_1*(xarr-xoff_beam)
                
        #### Interpolate pixel at shifted ycenter along beam
        if not dydx:
            dydx_0 = 0.
            dydx_1 = 0.
        
        #print 'DYDX: ', dydx, dydx_0, dydx_1
        
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
            
        #### Sensitivity
        sens = pyfits.open('../CONF/'+conf['SENSITIVITY_'+beam])[1].data
        sens_interp = np.interp(lam, sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)
        
        if xarr[keep].size > 1:
            full_sens[y0[keep]+NY,xpix[keep]] += sens_interp[keep]
            full_sens[y0[keep]+NY+1,xpix[keep]] += sens_interp[keep]
        
        if (lam_spec is not None) & (flux_spec is not None):
            spec_interp = np.interp(lam, lam_spec, flux_spec, left=0., right=0.)
            sens_interp *= spec_interp
        #
        if beam == 'A':
            sens_interp *= 0.5
        
        if beam == 'B':
            sens_interp *= 3.6 * 10
        
        if beam in ['C','D','E']:
            sens_interp *= 0.25
                
        stripe *= np.dot(np.ones((NY*2+1,1)), sens_interp.reshape(1,NX)) * grow_factor**2
        
        model += stripe
        
    return model, (xmi, xma, wavelength, full_sens)
    
def synthesize():
    import unicorn.reduce as red
    import time
    
    im = pyfits.open('COSMOS-15-F140W_inter.fits')
    gris = pyfits.open('COSMOS-15-G141_inter.fits')
    flux = im[1].data * 10**(-0.4*(26.46+48.6))*3.e18/1.392e4**2/1.e-17
    sh = im[1].data.shape
        
    ## stars
    xc, yc = 1184, 853
    dd = 50
    lam_spec = np.arange(1.e4,1.7e4)
    flux_spec = (lam_spec-1.4e4)*-0.00015 + 1

    # xc, yc = 1756, 1015
    # dd = 50
    # lam_spec = np.arange(1.e4,1.7e4)
    # flux_spec = (lam_spec-1.4e4)*-0.00001 + 1
    
    ## galaxy
    xc, yc = 513, 622
    dd = 50
    lam_spec = np.arange(1.e4,1.7e4)
    flux_spec = (lam_spec-1.4e4)*-0.00013 + 1
    
    orders, xi = red.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec)
    yord, xord = np.indices(orders.shape)
    
    non_zero = orders > 0
    xord, yord, ford = xord[non_zero], yord[non_zero], orders[non_zero]
    
    ys = orders.shape
    xord += xi[0]
    yord -= (ys[0]-1)/2
    
    # lam_spec = None
    # flux_spec = None
    # yi, xi = np.indices(sh)
    # pix = im[1].data > 0.5
    # flux_pix = flux[pix]
    
    model = im[1].data*0.

    t0 = time.time()
    
    for y in range(yc-dd, yc+dd):
        yyi = y+yord
        yuse = (yyi >= 0) & (yyi < sh[0])
        for x in range(xc-dd, xc+dd):
            xxi = x+xord
            use = (xxi >= 0) & (xxi < sh[1]) & yuse
            model[yyi[use], xxi[use]] += ford[use]*flux[y,x]*10
            # xxs = np.append(xxs, xxi[use])
            # yys = np.append(yys, yyi[use])
            # #print x,y
            # x0i, x1i = x+xi[0]+1, x+xi[1]+1
            # y0i, y1i = y-(ys[0]-1)/2+1, y+(ys[0]-1)/2+1+1
            # #
            # x0m, x1m = 0, ys[1]
            # y0m, y1m = 0, ys[0]
            # if x0i < 0:
            #     x0m = -x0i
            #     x0i = 0
            # if y0i < 0:
            #     y0m = -y0i
            #     y0i = 0
            # if x1i >= sh[1]:
            #     x1m = sh[1]-x1i
            #     x1i = sh[1]
            # if y1i >= sh[0]:
            #     y1m = sh[0]-y1i
            #     y1i = sh[0]
            #
            # print orders[y0m:y1m,x0m:x1m].shape, model[y0i:y1i,x0i:x1i].shape
            #
            ### speed test
            #model[y0i:y1i,x0i:x1i] += orders[y0m:y1m,x0m:x1m]*flux_xy*10
                
    t1 = time.time()
    print 'Run: %10.2f s.' %(t1-t0)
    
    #pyfits.writeto('test.fits', model, header=im[1].header, clobber=True)
    
    ds9 = threedhst.dq.myDS9()
    
    ds9.frame(1)
    ds9.view_fits(im[1])
    ds9.scale(-0.1*4.5, 5*4.5)
    
    ds9.frame(2)
    ds9.view_fits(gris[1])
    ds9.scale(-0.1,5)
    
    use = model > (0.1*model.max())
    norm = np.sum(gris[1].data[use]*model[use])/np.sum(model[use]**2)
    
    ds9.frame(3)
    ds9.view_array(model, header=gris[1].header)
    ds9.scale(-0.1/norm,5/norm)
    
    ds9.frame(4)
    ds9.view_array(gris[1].data-model*norm, header=gris[1].header)
    ds9.scale(-0.1,5)

    ds9.frame(3)
    ds9.xpa.set('pan to 702.2 627.2') # galaxy
    #ds9.xpa.set('pan to 1378.2 1022.2') # star2
    #ds9.xpa.set('pan to 1397.2 852.2') # star1
    ds9.xpa.set('zoom to 0.833333')
    ds9.xpa.set('match frames image')

class GrismModel():
    def __init__(self, root='GOODS-S-24'):
        im = pyfits.open(root+'-F140W_inter.fits')
        gris = pyfits.open(root+'-G141_inter.fits')
        flux = im[1].data * 10**(-0.4*(26.46+48.6))*3.e18/1.392e4**2/1.e-17
        sh = im[1].data.shape

        cat = threedhst.sex.mySexCat(root+'_inter.cat')
        seg = pyfits.open(root+'_seg.fits')
        full_model = pyfits.open(root+'_model0.fits')
        
          
def twod_spectrum(root='GOODS-S-24', id=328):
    import unicorn.reduce as red
    
    im = pyfits.open(root+'-F140W_inter.fits')
    gris = pyfits.open(root+'-G141_inter.fits')
    flux = im[1].data * 10**(-0.4*(26.46+48.6))*3.e18/1.392e4**2/1.e-17
    sh = im[1].data.shape
    
    cat = threedhst.sex.mySexCat(root+'_inter.cat')
    seg = pyfits.open(root+'_seg.fits')
    full_model = pyfits.open(root+'_model0.fits')
    
    ii = np.where(np.cast[int](cat.NUMBER) == id)[0][0]
    
    xc = np.float(cat.X_IMAGE[ii])
    yc = np.float(cat.Y_IMAGE[ii])
    
    #### Generate the model for a given object
    orders, xi = red.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=['A'])
    
    yord, xord = np.indices(orders.shape)
    # non_zero = orders > 0
    xord, yord, ford, word = xord, yord, orders, xi[2]
    
    ys = orders.shape
    xord += xi[0]
    yord -= (ys[0]-1)/2
    
    mask = seg[0].data == id
    xpix = xf[mask]
    ypix = yf[mask]
    object = model*0.
    
    for jj in range(xpix.size):
        x, y = xpix[jj], ypix[jj]
        xxi = x+xord
        yyi = y+yord
        use = (xxi >= 0) & (xxi < sh[1]) & (yyi >= 0) & (yyi < sh[0])
        object[yyi[use], xxi[use]] += ford[use]*flux[y,x]*10
    
    wave_2d = xi[2]
    sha = wave_2d.shape
    wave = np.zeros(sha[1])
    for i in range(sha[1]):
        wave[i] = wave_2d[:,i].max()
    
    xxi = np.round(xc)+xord
    
    dlam = (wave[wave > 0][1:]-wave[wave > 0][0:-1]).mean()
    ypix, xpix = np.indices(seg[0].data.shape)
    NY = ypix[seg[0].data == id].max()-ypix[seg[0].data == id].min()
    NX = xpix[seg[0].data == id].max()-xpix[seg[0].data == id].min()
    
    
def model_full_image(root='COSMOS-18'):
    import unicorn.reduce as red
    import time
    
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()

    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()

    ## XXX add test for user-defined .conv file
    se.copyConvFile()

    se.overwrite = True
    se.options['CATALOG_NAME']    = root+'_inter.cat'
    se.options['CHECKIMAGE_NAME'] = root+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = root+'-F140W_inter.fits[1]'
    se.options['FILTER']    = 'Y'

    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '1' 
    se.options['ANALYSIS_THRESH']  = '1'
    se.options['MAG_ZEROPOINT'] = '26.46'

    #### Run SExtractor
    status = se.sextractImage(root+'-F140W_inter.fits[0]')
        
    im = pyfits.open(root+'-F140W_inter.fits')
    gris = pyfits.open(root+'-G141_inter.fits')
    flux = im[1].data * 10**(-0.4*(26.46+48.6))*3.e18/1.392e4**2/1.e-17
    sh = im[1].data.shape
    
    t0 = time.time()
    
    seg = pyfits.open(root+'_seg.fits')[0].data
        
    cat = threedhst.sex.mySexCat(root+'_inter.cat')
    #### Trim faint sources
    mag = np.cast[float](cat.MAG_AUTO)
    q = np.where(mag > 23.5)[0]
    if len(q) > 0:
        threedhst.showMessage('Trimming objects with direct M < %6.2f.' %(24))
        numbers = np.cast[int](cat.NUMBER)[q]
        cat.popItem(numbers)
    
    cat.write()
    print len(seg[seg > 0])
    
    yf, xf = np.indices(im[1].data.shape)
    
    lam_spec = None
    flux_spec = None
    
    #### Try a default blue spectrum
    # lam_spec = np.arange(1.e4,1.7e4)
    # flux_spec = (lam_spec-1.4e4)*-0.0001 + 1
    
    noNewLine = '\x1b[1A\x1b[1M'

    model = im[1].data*0.
    
    BEAMS = ['A','B','C','D','E']
    
    for ii,id in enumerate(cat.id):
        print noNewLine + 'ID: %d' %(id)
        xc = np.float(cat.X_IMAGE[ii])
        yc = np.float(cat.Y_IMAGE[ii])
        #
        orders, xi = red.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS)
        yord, xord = np.indices(orders.shape)
        non_zero = orders > 0
        xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
        
        ys = orders.shape
        xord += xi[0]
        yord -= (ys[0]-1)/2
        #
        mask = seg == id
        xpix = xf[mask]
        ypix = yf[mask]
        object = model*0.
        #
        for jj in range(xpix.size):
            x, y = xpix[jj], ypix[jj]
            xxi = x+xord
            yyi = y+yord
            use = (xxi >= 0) & (xxi < sh[1]) & (yyi >= 0) & (yyi < sh[0])
            object[yyi[use], xxi[use]] += ford[use]*flux[y,x]*10
        #
        model += object #*norm
        
    #### Second loop, subtract all other object and renormalize for a given object
    t1 = time.time()
    print '\nFirst pass: %.1fs.\n\nSecond pass...\n' %(t1-t0)
    model0 = model*1.
    model1 = model0*0.
    norm_coeffs = cat.id*0.
    
    pyfits.writeto(root+'_model0.fits', model0, header=gris[1].header, clobber=True)
    
    so = np.argsort(np.cast[float](cat.MAG_AUTO))
    for count,ii in enumerate(so):
        id = cat.id[ii]
        print noNewLine + 'ID: %d, %3d' %(id, count*1./(len(so)-1)*100) + '%'
        xc = np.float(cat.X_IMAGE[ii])
        yc = np.float(cat.Y_IMAGE[ii])
        #
        orders, xi = red.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec)
        yord, xord = np.indices(orders.shape)
        non_zero = orders > 0
        xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
        ys = orders.shape
        xord += xi[0]
        yord -= (ys[0]-1)/2
        #
        mask = seg == id
        xpix = xf[mask]
        ypix = yf[mask]
        object = model*0.
        # wave = model*0
        # wave_weight = model*0
        # sens = model*0
        #
        for jj in range(xpix.size):
            x, y = xpix[jj], ypix[jj]
            xxi = x+xord
            yyi = y+yord
            use = (xxi >= 0) & (xxi < sh[1]) & (yyi >= 0) & (yyi < sh[0])
            object[yyi[use], xxi[use]] += ford[use]*flux[y,x]*10
            # wave[yyi[use], xxi[use]] += word[use]*flux[y,x]
            # wave_weight[yyi[use], xxi[use]] += flux[y,x]
            # sens[yyi[use], xxi[use]] += ford[use]
        #
        #### Get slope
        # data = gris[1].data-(model0 - object)
        # wok = (wave_weight > 0) & (np.isfinite(data))
        # wave[wok] /= wave_weight[wok]
        # data[wok == False] = 0
        # wht = data*1.
        # yd, xd = np.indices(data.shape)
        # yd = np.unique(yd[wok])
        # for jj in yd:
        #     wline = wave[jj,:]
        #     dline = data[jj,:]
        #     sline = sens[jj,:]/np.max(sens[jj,:])
        #     dok = (dline > 0.2*dline.max()) & (np.abs(wline-1.4e4) < 300)
        #     dline /= sline
        #     dline -= np.median(dline[dok])
        #     data[jj,:] = dline  
        # #
        # wok = (wave_weight > 0) & (np.isfinite(data)) & (np.abs(wave-1.4e4) > 600) & (wht > 0.4*wht.max())
        # data[wok == False] = 0
        # wht[wok == False] = 0
        # sens[wok == False] = 0
        # a = (data)/((wave-1.4e4))
        # aout = np.sum(a*wht)/np.sum(wht)
        # lam_spec = np.arange(1.e4,1.7e4)
        # flux_spec = (lam_spec-1.4e4)*aout*2 + 1
        
        overlap = (object > 0.1*object.max()) & (model0 < 2*object.max()) & (gris[1].data != 0)
        if len(object[overlap]) > 1:
            cleaned = gris[1].data - (model0 - object)
            norm = np.sum(object[overlap]*cleaned[overlap])/np.sum(object[overlap]**2)
        else:
            norm = 1.
        if (not np.isfinite(norm)) | (norm > 2.5):
            norm = 1.
        #
        norm_coeffs[ii] = norm
        model1 += object*norm
        
    #
    t2 = time.time()
    print '\n\nSecond pass: %.1fs.\n' %(t2-t1)
    
    pyfits.writeto(root+'_model1.fits', model1, header=gris[1].header, clobber=True)
    #pyfits.writeto('junk_ratio.fits', gris[1].data/model, header=gris[1].header, clobber=True)
    
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
    
    BEAMS = ['A','B']
    BEAMS = ['A','B','C','D'] #,'E']
    noNewLine = '\x1b[1A\x1b[1M'

    # orders, xi = red.grism_model(507, 850, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, pad=0)
    # yord, xord = np.indices(orders.shape)
    # non_zero = orders > 0
    # xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
    # 
    # ys = orders.shape
    # xord += xi[0]
    # yord -= (ys[0]-1)/2
    
    orders, xi = red.grism_model(507,507, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, pad=0)
    yord, xord = np.indices(orders.shape)
    non_zero = orders > 0
    xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]

    ys = orders.shape
    xord += xi[0]
    yord -= (ys[0]-1)/2
    
    skip = 5
    
    for y in range(460, 540):
        print noNewLine + 'Y: %d' %(y)
        if (y % 20) == 0:
            pyfits.writeto('stripe.fits', model/model.max(), clobber=True)
        
        for x in range(-190,1014+85):
            if ((x % skip) + (y % skip)) == 0:
                orders, xi = red.grism_model(x+1, y+1, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=1, pad=0, dydx=False)
                yord, xord = np.indices(orders.shape)
                non_zero = orders > 0
                xord, yord, ford, word = xord[non_zero], yord[non_zero], orders[non_zero], xi[2][non_zero]
                ys = orders.shape
                xord += xi[0]
                yord -= (ys[0]-1)/2
        
            xxi = x+xord
            yyi = y+yord
            use = (xxi >= 0) & (xxi < 1014) & (yyi >= 0) & (yyi < 1014)
            model[yyi[use], xxi[use]] += ford[use]
        #
        #model += object #*norm
    #
    pyfits.writeto('stripe.fits', model/model.max(), clobber=True)
    
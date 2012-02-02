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

noNewLine = '\x1b[1A\x1b[1M'

import time

import threedhst

def interlace_combine(root='COSMOS-1-F140W', view=True):
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    
    if view:
        ds9 = threedhst.dq.myDS9()

    run = threedhst.prep_flt_files.MultidrizzleRun(root)
    #run.blot_back(ii=0, copy_new=True)
    
    im = pyfits.open(os.getenv('iref')+'ir_wfc3_map.fits')
    PAM = im[1].data
    im.close()
    
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
        
    for i,flt in enumerate(run.flt):
        print flt
        flt = run.flt[i]
        im = pyfits.open(flt+'.fits')
        
        #### Use the pixel area map correction
        im[1].data *= PAM
        #### Divide by 4 to conserve surface brightness with smaller output pixels
        im[1].data /= 4
        
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
        if view:
            ds9.view_array(inter, header=header)
            ds9.scale(-0.1,5)
            
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=inter, header=header)
    wht = pyfits.ImageHDU(data=N, header=header)
    pyfits.HDUList([hdu,sci,wht]).writeto(root+'_inter.fits', clobber=True)
    
    try:
        os.remove(root+'_inter_sci.fits')
        os.remove(root+'_inter_wht.fits')
    except:
        pass
    #
    new_coeffs_dat(input=run.flt[0]+'_coeffs1.dat', output='scale_coeffs.dat')
    
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

def new_coeffs_dat(input='ibhm29wlq_flt_coeffs1.dat',output='scale_coeffs.dat', factor=2, pad=60):
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
    
    for key in keys:
        header.update(key,header[key]/factor)
    
    header.update('IDCTAB','')
    
    header.update('EXPTIME',1)
    
    return header

def grism_model(xc_full=244, yc_full=1244, lam_spec=None, flux_spec=None, grow_factor=2, pad = 60, BEAMS=['A','B','C','D','E'], dydx=True):
    
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    import threedhst.dq
    #ds9 = threedhst.dq.myDS9()

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
      
    #xmi, xma = -580, 730
    xmi, xma = 1000,-1000
    if 'A' in BEAMS:
        xmi, xma = min(xmi,10), max(xma,213)
    if 'B' in BEAMS:
        xmi, xma = min(xmi,-210), max(xma, -170)
    if 'C' in BEAMS:
        xmi, xma = min(xmi,207), max(xma,454)
    if 'D' in BEAMS:
        xmi, xma = min(xmi, 469), max(xma,668)
    if 'E' in BEAMS:
        xmi, xma = min(xmi, -550), max(xma,-400)
        
    NX,NY = xma-xmi, 8
    
    xmi *= grow_factor
    xma *= grow_factor
    NX *= grow_factor
    NY *= grow_factor
    
    model = np.zeros((NY*2+1, NX), dtype=np.float)
    wavelength = model*0
    full_sens = model*0
    
    #print model.shape, NX, NY
    
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
        if beam == 'B':
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
        
        yoff_save = dydx_0+dydx_1*(xmi_beam+xma_beam)/2
              
        #### Sensitivity
        sens = pyfits.open('../CONF/'+conf['SENSITIVITY_'+beam])[1].data
        sens_interp = np.interp(lam, sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)
        
        if xarr[keep].size > 1:
            full_sens[y0[keep]+NY,xpix[keep]] += sens_interp[keep]
            full_sens[y0[keep]+NY+1,xpix[keep]] += sens_interp[keep]
        
        if (lam_spec is not None) & (flux_spec is not None):
            #spec_interp = threedhst.utils.interp_conserve(lam, lam_spec, flux_spec, left=0., right=0.)            
            spec_interp = threedhst.utils.interp_conserve_c(lam, lam_spec, flux_spec, left=0., right=0.)
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
        
    return model, (xmi, xma, wavelength, full_sens, yoff_save)


class GrismModel():
    def __init__(self, root='GOODS-S-24', grow_factor=2, pad=60, LIMITING_MAGNITUDE=24):
        """
        Initialize: set padding, growth factor, read input files.
        """
        self.root=root
        
        header = pyfits.getheader(root+'-F140W_inter.fits',0)
        self.filter = header['FILTER']
        ZPs = {'F125W':26.25, 'F140W':26.46, 'F160W':25.96}
        self.ZP = ZPs[self.filter]

        PLAMs = {'F125W':1.2486e4, 'F140W':1.3923e4, 'F160W': 1.5369e4}
        self.PLAM = PLAMs[self.filter]
        
        if not os.path.exists(root+'_seg.fits'):
            threedhst.showMessage('Running SExtractor...')
            self.find_objects(LIMITING_MAGNITUDE=LIMITING_MAGNITUDE)
        
        self.grow_factor = grow_factor
        self.pad = pad
        
        threedhst.showMessage('Reading FITS files and catalog...')
        self.read_files()
        self.init_object_spectra()
        try:
            test = self.total_fluxes
        except:
            self.make_total_flux()
            
    def read_files(self):
        """
        Read FITS files, catalogs, and segmentation images for a pair of 
        interlaced exposures.
        """
        self.im = pyfits.open(self.root+'-F140W_inter.fits')
        
        self.make_flux()
                
        self.gris = pyfits.open(self.root+'-G141_inter.fits')
        self.sh = self.im[1].data.shape
        self.cat = threedhst.sex.mySexCat(self.root+'_inter.cat')
        self.segm = pyfits.open(self.root+'_seg.fits')
        self.model = self.im[1].data*0.
        self.yf, self.xf = np.indices(self.im[1].data.shape)
        
    def make_flux(self):
        """
        Convert from e/s to physical fluxes using the ZP and pivot wavelength 
        for the direct image filter
        """
        self.flux = self.im[1].data * 10**(-0.4*(self.ZP+48.6))* 3.e18 / self.PLAM**2 / 1.e-17
        self.flux_fnu = self.im[1].data * 10**(-0.4*(self.ZP+48.6))
        
    def find_objects(self, LIMITING_MAGNITUDE=23.5):
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
        se.options['CHECKIMAGE_NAME'] = self.root+'_seg.fits'
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
        q = np.where(mag > LIMITING_MAGNITUDE)[0]
        if len(q) > 0:
            numbers = np.cast[int](self.cat.NUMBER)[q]
            self.cat.popItem(numbers)

        self.cat.write()
        self.segm = pyfits.open(self.root+'_seg.fits')
        
        threedhst.sex.sexcatRegions(self.root+'_inter.cat', self.root+'_inter.reg', format=1)
        
        self.make_total_flux()
        self.init_object_spectra()
        
    def init_object_spectra(self):
        
        self.lam_specs = {}
        self.flux_specs = {}
        for obj in self.cat['NUMBER']:
            self.lam_specs[obj] = None
            self.flux_specs[obj] = None
        
    def make_total_flux(self):
        
        #### Total flux
        try:
            flux = self.flux
        except:
            self.make_flux()
            flux = self.flux
        #
        self.total_fluxes = {}
        print 'Total flux...\n\n'
        for obj in self.cat['NUMBER']:
            print noNewLine+'Object #%d' %(obj)
            self.total_fluxes[obj] = np.sum(flux[self.segm[0].data == obj])
            
    def compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D','E'], verbose=False, normalize=True):
        import unicorn.reduce as red
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == id)[0][0]
        
        xc = np.float(self.cat.X_IMAGE[ii])
        yc = np.float(self.cat.Y_IMAGE[ii])
        
        #### Normalize the input spectrum to the direct image passband
        t0 = time.time()
        if (lam_spec is not None) & (normalize is True):
            if verbose:
                t1 = time.time(); dt = t1-t0; t0=t1
                print 'Input spectrum is provided, normalize it for filter, %s (%.2f s)' %(self.filter, dt)
            
            x_filt, y_filt = np.loadtxt(os.getenv('iref')+'/'+self.filter + '.dat', unpack=True)
            y_filt_int = np.interp(lam_spec, x_filt, y_filt)
            filt_norm = np.trapz(y_filt_int*flux_spec, lam_spec) / np.trapz(y_filt_int, lam_spec)
            flux_spec /= filt_norm
            
        #### Define the grism dispersion for the central pixel of an object.  
        #### Assume doesn't change across the object extent
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Define grism dispersion (%.2f s).' %(dt)
            
        orders, self.xi = red.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=self.grow_factor, pad=self.pad)
        yord, xord = np.indices(orders.shape)
        non_zero = orders > 0
        xord, yord, ford, word, sord = xord[non_zero], yord[non_zero], orders[non_zero], self.xi[2][non_zero], self.xi[3][non_zero]
        
        ys = orders.shape
        xord += self.xi[0]
        yord -= (ys[0]-1)/2
        #
        mask = self.segm[0].data == id
        xpix = self.xf[mask]
        ypix = self.yf[mask]
        object = self.model*0.
        
        #### Loop through pixels defined within the segmentation region, summing
        #### dispersed flux
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Loop through object pixels within the segmentation region (%.2f s).' %(dt)
        for jj in range(xpix.size):
            x, y = xpix[jj], ypix[jj]
            xxi = x+xord
            yyi = y+yord
            use = (xxi >= 0) & (xxi < self.sh[1]) & (yyi >= 0) & (yyi < self.sh[0])
            object[yyi[use], xxi[use]] += ford[use]*self.flux[y,x]*10

        self.object = object
        
        #### Make full images that show the wavelength, sensitivity
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Make full image arrays that show the wavelength, sensitivity (%.2f s)' %(dt)
        
        xxi = np.int(np.round(xc))+xord
        use = (xxi >= 0) & (xxi < self.sh[1])
        cut = np.zeros(self.sh[1])
        cut[xxi[use]] = word[use] 
        self.object_wave = np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
        
        cut[xxi[use]] = sord[use] 
        self.object_sens = np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
    
    def refine_model(self, ID, BEAMS=['A','B','C','D','E'], view=False, model_has_object=True):
        from scipy import polyfit,polyval

        if ID not in self.cat['NUMBER']:
            print '#%d not in the catalog.' %(ID)
            return False
        
        #### First iteration
        self.compute_object_model(id=ID, BEAMS=['A'], lam_spec=self.lam_specs[ID], flux_spec = self.flux_specs[ID])
        first_model = self.object*1.
        
        if model_has_object:
            self.model -= first_model
            
        whts = self.object**4
        finite_mask = np.isfinite(self.model) 
        whts[~finite_mask] = 0
        whts[self.gris[1].data == 0] = 0

        ratio = (self.gris[1].data - self.model) / self.object
        ratio[~finite_mask] = 0.
        ratio[~np.isfinite(ratio)] = 0
        whts[~np.isfinite(ratio)] = 0

        ratio_extract = np.sum(ratio*whts, axis=0) / np.sum(whts, axis=0)
        ratio_extract[np.sum(ratio*whts, axis=0) == 0] = 0.

        wave = np.mean(self.object_wave, axis=0)
        
        keep = (wave > 1.1e4) & (wave < 1.65e4)
        if len(keep[keep]) < 10:
            self.model += self.object
            return True
            
        coeff = polyfit(wave[keep], ratio_extract[keep], 1)
        ratio_extract[~keep] = polyval(coeff, wave[~keep])
        
        if self.flux_specs[ID] is not None:
            self.flux_specs[ID] *= ratio_extract[wave > 0]
        else:
            self.lam_specs[ID] = wave[wave > 0]
            self.flux_specs[ID] = ratio_extract[wave > 0]
            
        self.compute_object_model(id=ID, BEAMS=BEAMS, lam_spec = self.lam_specs[ID], flux_spec = self.flux_specs[ID], normalize=False)
        second_model = self.object*1.
        
        #### View the results
        if view.__class__ == threedhst.dq.myDS9:
            print 'Display results for #%d to DS9' %(ID)
            diff_first = self.gris[1].data - self.model - first_model
            diff_second = self.gris[1].data - self.model - second_model

            view.frame(1)
            view.v(self.gris[1].data, vmin=-0.1, vmax=5)
            view.frame(2)
            view.v(diff_first, vmin=-0.1, vmax=5)
            view.frame(3)
            view.v(diff_second, vmin=-0.1, vmax=5)

            idx = (self.cat['NUMBER'] == ID)
            view.set('pan to %f %f' %(self.cat['X_IMAGE'][idx][0], self.cat['Y_IMAGE'][idx][0]))
            view.set('zoom to 1.4')
            view.set('cmap value 4.43936 0.462911')
            view.set('cmap HSV')
            view.set('cmap value 7.20824 0.48544')
            view.set('regions file %s_inter.reg' %(self.root))
            view.match_all(alignType='image')
            
        self.model += second_model
    
    def show_ratio_spectrum(self, ID):
        if ID not in self.cat['NUMBER']:
            print '#%d not in the catalog.' %(ID)
            return False
        
        plt.plot(self.lam_specs[ID], self.flux_specs[ID]*self.total_fluxes[ID])
        
    def full_model(self, BEAMS=['A','B','C','D','E'], save_ext='_model0', iteration_number=0):
        
        self.model*=0.
        threedhst.showMessage('Making full object model...')
        
        so = np.argsort(self.cat['MAG_AUTO'])
        N = len(so)
        for i, id in enumerate(self.cat['NUMBER'][so]):
            print noNewLine+'Full model, object #%d (%d/%d)' %(id, i+1, N)
            self.refine_model(id, BEAMS=BEAMS, view=ds9, model_has_object=(iteration_number > 0))      
        if save_ext:
            pyfits.writeto(self.root+save_ext+'.fits', self.model, clobber=True)
    
    def twod_spectrum(self, id=328, grow=1, miny=80):
        """
        Extract 2D direct + grism thumbnails.
        """
        self.compute_object_model(id=id, BEAMS=['A'])
        #self.model = pyfits.open('GOODS-S-24_model0.fits')[0].data
        
        seg = self.segm[0].data
        ypix, xpix = np.indices(seg.shape)
        NY = ypix[seg == id].max()-ypix[seg == id].min()
        NX = xpix[seg == id].max()-xpix[seg == id].min()
        NT = np.max([NY*grow, NX*grow, miny])
        #print NT, NX, NY
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == id)[0][0]
        
        xc = np.int(np.round(np.float(self.cat.X_IMAGE[ii])))-1
        yc = np.int(np.round(np.float(self.cat.Y_IMAGE[ii])))-1
        
        self.direct_thumb = self.im[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_wht = self.im[2].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        
        #### dydx offset for grism trace
        yc = np.int(np.round(np.float(self.cat.Y_IMAGE[ii])+self.xi[4]))-1
        
        #pyfits.writeto('junk.fits', self.direct_thumb, clobber=True)
        
        xmin = xpix[self.object_wave > 0].min()
        xmax = xpix[self.object_wave > 0].max()
                
        prim = pyfits.PrimaryHDU()
        prim.header.update('POINTING', self.root)
        prim.header.update('ID', id)
        prim.header.update('BUNIT', 'ELECTRONS/S')
        prim.header.update('EXPTIME', self.gris[0].header['EXPTIME'])
        prim.header.update('FILTER', self.gris[0].header['FILTER'])
        
        extensions = [prim]
        header = pyfits.ImageHDU().header
        
        #### Direct thumbnails
        header.update('EXTNAME','DSCI')
        extensions.append(pyfits.ImageHDU(self.direct_thumb, header))
        header.update('EXTNAME','DWHT')
        extensions.append(pyfits.ImageHDU(np.cast[int](self.direct_wht), header))
        
        self.wave = self.object_wave[507, xmin:xmax]
        self.dwave = self.wave[1]-self.wave[0]
        
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
        
        self.grism_sci = self.gris[1].data[yc-NT/2:yc+NT/2, xmin:xmax]
        header.update('EXTNAME','SCI')
        extensions.append(pyfits.ImageHDU(self.grism_sci, header))
        
        self.grism_wht = np.cast[int](self.gris[2].data[yc-NT/2:yc+NT/2, xmin:xmax])
        header.update('EXTNAME','WHT')
        extensions.append(pyfits.ImageHDU(self.grism_wht, header))
        
        self.the_object = self.object[yc-NT/2:yc+NT/2, xmin:xmax]
        header.update('EXTNAME','MODEL')
        extensions.append(pyfits.ImageHDU(self.the_object, header))
        
        self.full_model = self.model[yc-NT/2:yc+NT/2, xmin:xmax]
        header.update('EXTNAME','CONTAM')
        extensions.append(pyfits.ImageHDU(self.full_model-self.the_object, header))
                
        header.update('EXTNAME','WAVE')
        extensions.append(pyfits.ImageHDU(self.wave, header))
        
        self.sens = self.object_sens[507, xmin:xmax]
        header.update('EXTNAME','SENS')
        header.update('BUNIT','(ELECTRONS/S) / (FLAM*1.E17)')
        extensions.append(pyfits.ImageHDU(self.sens, header))

        hdu = pyfits.HDUList(extensions)
        hdu.writeto(self.root+'_%05d.fits' %(id), clobber=True)
        
        self.twod = hdu
    
    def extract_1d(self, id=284):
        """
        Extract a 1D spectrum after cutting out the 2D spectrum using `twod_spectrum`.
        
        Trying to implement smoothing of the sensitivity curve, but still doesn't
        look right.  Should the whole curve be smoothed, or just the wings
        extended?  In principle, it should work if you provide the correct input 
        spectrum to the model.
        
        Still need to implement computing the quantitative contamination.
        """
        from scipy import polyfit, polyval
        
        if not os.path.exists(self.root+'_%05d.fits' %(id)):
            self.twod_spectrum(id=id)
        
        im = pyfits.open(self.root+'_%05d.fits' %(id))
        shape = im['SCI'].data.shape
        
        profile = np.sum(im['MODEL'].data, axis=1)
        profile /= np.sum(profile)        
        profile = np.dot(profile.reshape(-1,1), np.ones((1,shape[1])))
        
        flux_dn = im['SCI'].data*im[0].header['EXPTIME']
        
        wht = im['WHT'].data
        
        x,y = np.where(wht == 0)
        
        pad = 1
        for xi, yi in zip(x,y):
            sub = flux_dn[xi-pad:xi+pad+1,yi-pad:yi+pad+1]
            wht_i = wht[xi-pad:xi+pad+1,yi-pad:yi+pad+1]
            if (np.sum(wht_i) != 0.0):
                flux_dn[xi,yi] = np.mean(sub[np.isfinite(sub)])
                wht[xi,yi] = 1
                
        sigma = np.sqrt(flux_dn + 21**2)
        sigma[flux_dn < 0] = 21
        
        sum_dn = np.sum(flux_dn/im[0].header['EXPTIME']*profile/sigma**2*wht, axis=0)/np.sum(profile**2/sigma**2*wht)
        
        obj_profile = np.sum(im['DSCI'].data, axis=0)
        obj_profile /= np.sum(obj_profile)
        
        wave = im['WAVE'].data
        sens = im['SENS'].data
        
        sens_smooth = np.convolve(sens, obj_profile, mode='same')
        
        ## lo = 1.1157e4, 1.6536
        yint = np.interp(np.array([1.1157e4, 1.6536e4]), wave, sens)
        yint_smooth = np.interp(np.array([1.1157e4, 1.6536e4]), wave, sens_smooth)
        
        sens_corrected = sens*1.
        sens_corrected[wave < 1.1157e4] = sens_smooth[wave < 1.1157e4]+yint[0]-yint_smooth[0]
        sens_corrected[wave > 1.6536e4] = sens_smooth[wave > 1.6536e4]+yint[1]-yint_smooth[1]
        
        self.oned_wave = wave
        self.oned_flux = sum_dn / sens_corrected
        
        coeffs = polyfit(wave[(wave > 1.15e4) & (wave < 1.64e4)], self.oned_flux[(wave > 1.15e4) & (wave < 1.64e4)], 1)
        self.oned_flux_linear = polyval(coeffs, wave)/polyval(coeffs, 1.4e4)
        
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
    
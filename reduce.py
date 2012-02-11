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
import pyraf
from pyraf import iraf

noNewLine = '\x1b[1A\x1b[1M'

import time

import threedhst

import unicorn.utils_c as utils_c

def go_all():
    
    for field in ['AEGIS','COSMOS','GOODS-S','UDS'][1:]:
        os.chdir(unicorn.GRISM_HOME+field+'/PREP_FLT')
        files=glob.glob(field+'-[0-9]*asn.fits')
        for file in files:
            print file
            #unicorn.reduce.interlace_combine(file.split('_asn')[0], view=False, use_error=True, make_undistorted=False)
            if 'G141' in file:
                model = unicorn.reduce.GrismModel(file.split('-G141')[0], LIMITING_MAGNITUDE=26)
                model.get_corrected_wcs()
                model.make_wcs_region_file()

def interlace_combine(root='COSMOS-1-F140W', view=True, use_error=True, make_undistorted=False):
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
    inter_sci = np.zeros((2028+pad,2028+pad))
    if use_error:
        inter_err = np.zeros((2028+pad,2028+pad))
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
        im[2].data /= 4
        
        if i == 0:
            h0 = im[0].header
            h1 = im[1].header
            header = red.scale_header_wcs(h1.copy(), factor=2, pad=pad)
            header_wht = header.copy(); header_wht.update('EXTNAME','ERR')
        
        dx = np.int(np.round((xsh[i]-xsh[0])*2))
        dy = np.int(np.round((ysh[i]-ysh[0])*2))
        #
        dx = dxs[i]
        dy = dys[i]
        #
        #use = ((im[3].data & 4096) == 0) & ((im[3].data & 4) == 0) #& (xi > np.abs(dx/2)) & (xi < (1014-np.abs(dx/2))) & (yi > np.abs(dy/2)) & (yi < (1014-np.abs(dy/2)))
        #
        use = ((im[3].data & 4) == 0) & (~hot_pix) & ((im[3].data & 32) == 0)
        inter_sci[yi[use]*2+dy,xi[use]*2+dx] += im[1].data[use]
        if use_error:
            inter_err[yi[use]*2+dy,xi[use]*2+dx] += im[2].data[use]
        
        N[yi[use]*2+dy,xi[use]*2+dx] += 1
        #
        if view:
            ds9.view_array(inter_sci, header=header)
            ds9.scale(-0.1,5)
    
    if use_error:
        h0.update('WHTERROR',True,comment='WHT extension is FLT[err,1]')
    else:
        h0.update('WHTERROR',True,comment='WHT extension is FLT[err,1]')
        
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=inter_sci, header=header)
    if use_error:
        wht = pyfits.ImageHDU(data=inter_err*N, header=header)
    else:
        wht = pyfits.ImageHDU(data=N, header=header)
        
    image = pyfits.HDUList([hdu,sci,wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
    
    image.writeto(root+'_inter.fits', clobber=True)
    
    if make_undistorted:
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
    beam_index = np.zeros(NX, dtype=np.int)
    
    #print model.shape, NX, NY
    
    yi, xi = np.indices(model.shape)
    xi += xmi
    yi -= NY
    xarr, xpix, yarr = np.arange(xma-xmi)+xmi, np.arange(xma-xmi), np.arange(NY*2+1)-NY
    
    for ib,beam in enumerate(BEAMS):
       
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
            #
            beam_index[xpix[keep]] += 2**ib
            
        yoff_save = dydx_0+dydx_1*(xmi_beam+xma_beam)/2
        if beam == 'A':
            yoff_array = dydx_0+dydx_1*(xarr-xoff_beam)
            #print 'Yoff_array', xarr.shape, yoff_array.shape
            
        #### Sensitivity
        sens = pyfits.open('../CONF/'+conf['SENSITIVITY_'+beam])[1].data
        #sens_interp = np.interp(lam, sens.field('WAVELENGTH'), sens.field('SENSITIVITY')*1.e-17, left=0., right=0.)
        sens_interp = utils_c.interp_c(lam, np.array(sens.field('WAVELENGTH'), dtype=np.float64), np.array(sens.field('SENSITIVITY'), dtype=np.float64)*1.e-17, extrapolate=0.)
        
        if xarr[keep].size > 1:
            full_sens[y0[keep]+NY,xpix[keep]] += sens_interp[keep]
            full_sens[y0[keep]+NY+1,xpix[keep]] += sens_interp[keep]
        
        if (lam_spec is not None) & (flux_spec is not None):
            #spec_interp = threedhst.utils.interp_conserve(lam, lam_spec, flux_spec, left=0., right=0.)            
            spec_interp = utils_c.interp_conserve_c(lam, lam_spec, flux_spec, left=0., right=0.)
            spec_interp[~np.isfinite(spec_interp)] = 0
            spec_interp[lam < np.min(lam_spec)] = 0
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
        #print beam, model.max()
        
    return model, (xmi, xma, wavelength, full_sens, yoff_array, beam_index)


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
        print 'Read files...'
        self.read_files()
        print 'Init object spectra...'
        self.init_object_spectra()
        try:
            test = self.total_fluxes
        except:
            self.make_total_flux()
        
        self.ra_wcs, self.dec_wcs = None, None
        
    def read_files(self):
        """
        Read FITS files, catalogs, and segmentation images for a pair of 
        interlaced exposures.
        """
        self.im = pyfits.open(self.root+'-F140W_inter.fits')
        self.make_flux()
                
        self.gris = pyfits.open(self.root+'-G141_inter.fits')
        self.gris[1].data = np.array(self.gris[1].data, dtype=np.double)
        self.sh = self.im[1].data.shape
        
        self.cat = threedhst.sex.mySexCat(self.root+'_inter.cat')
        self.cat.x_pix = np.cast[float](self.cat.X_IMAGE)
        self.cat.y_pix = np.cast[float](self.cat.Y_IMAGE)
        self.cat.mag = np.cast[float](self.cat.MAG_AUTO)
        
        self.segm = pyfits.open(self.root+'_seg.fits')
        self.segm[0].data = np.array(self.segm[0].data, dtype=np.uint)
        
        self.model = np.zeros(self.sh, dtype=np.double)
        self.object = np.zeros(self.sh, dtype=np.double)
        self.test_object = np.zeros(self.sh)
        self.yf, self.xf = np.indices(self.sh)
        
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
        
        try:
            self.make_total_flux()
        except:
            pass
        self.init_object_spectra()
        
    def get_corrected_wcs(self, verbose=True):
        """
        Use iraf.dither.tran to get the real, undistorted WCS coordinates of each object
        """
        import iraf
        from iraf import stsdas
        from iraf import dither
        
        try:
            import stsci.tools.wcsutil as wcsutil
            wcs = wcsutil.WCSObject(self.root+'-F140W_drz.fits[1]')
        except:
            print 'Failed: import stsci.tools.wcsutil'
            wcs = None
            
        NOBJ = len(self.cat.id)
        
        asn = threedhst.utils.ASNFile(self.root+'-F140W_asn.fits')
        flt = asn.exposures[0]+'_flt.fits'
        fp = open('/tmp/%s.flt_xy' %(self.root),'w')
        for i in range(NOBJ):
            xi = (self.cat['X_IMAGE'][i]-self.pad/2.)/self.grow_factor
            yi = (self.cat['Y_IMAGE'][i]-self.pad/2.)/self.grow_factor
            fp.write('%.2f %.2f\n' %(xi,  yi))
        
        fp.close()
        
        print 'Running iraf.tran()'
        threedhst.process_grism.flprMulti()
        
        status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-F140W_drz.fits[1]', direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
        
        self.ra_wcs, self.dec_wcs = np.zeros(NOBJ, dtype=np.double), np.zeros(NOBJ, dtype=np.double)
                    
        fp = open('%s_inter.cat.wcsfix' %(self.root),'w')
        fp.write('# id    mag   ra    dec\n')
        for i, line in enumerate(status[-NOBJ:]):
            if verbose:
                print noNewLine+'Running iraf.xy2rd (%4d/%4d)' %(i+1, NOBJ)
            spl = np.cast[float](line.split())
            #
            if wcs is not None:
                self.ra_wcs[i], self.dec_wcs[i] = wcs.xy2rd((spl[2],spl[3]))
            else:
                status_xy2rd = iraf.xy2rd(infile=self.root+'-F140W_drz.fits[1]', x=spl[2], y=spl[3], hms=iraf.no, Stdout=1)
                self.ra_wcs[i], self.dec_wcs[i] = np.double(iraf.xy2rd.ra), np.double(iraf.xy2rd.dec)
            
            fp.write('%5d  %.3f  %.6f  %.6f\n' %(self.cat.id[i], self.cat.mag[i], self.ra_wcs[i], self.dec_wcs[i]))
            
        fp.close()
        if verbose:
            print 'Corrected coordinates: %s_inter.cat.wcsfix' %(self.root)
            
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
        print 'Total flux...'
       
        # ## Test, compare c version
        # t0 = time.time()
        # for obj in self.cat['NUMBER']:
        #     print noNewLine+'Object #%d' %(obj)
        #     self.total_fluxes[obj] = np.sum(flux[self.segm[0].data == obj])
        # 
        # told = time.time()
        #self.new_fluxes = {}

        fluxes = utils_c.total_flux(self.flux, self.segm[0].data)
        for obj in self.cat['NUMBER']:
            #self.new_fluxes[obj] = fluxes[obj]
            self.total_fluxes[obj] = fluxes[obj]
                    
        # tnew = time.time()
        # print 'Compare total fluxes:  %.3f  %.3f' %(told-t0, tnew-told)
        # for obj in self.cat['NUMBER']:
        #     print ' %.7f %.7f' %(self.total_fluxes[obj], self.total_fluxes[obj]-self.new_fluxes[obj])
        
            
    def compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D','E'], verbose=False, normalize=True):
        import unicorn.reduce
        
        t0 = time.time()
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == id)[0][0]
        
        self.id = id
        
        xc = np.float(self.cat.X_IMAGE[ii])
        yc = np.float(self.cat.Y_IMAGE[ii])
        
        #### Normalize the input spectrum to the direct image passband
        t0 = time.time()
        if (lam_spec is not None) & (normalize is True):            
            x_filt, y_filt = np.loadtxt(os.getenv('iref')+'/'+self.filter + '.dat', unpack=True)
            #y_filt_int = np.interp(lam_spec, x_filt, y_filt)
            y_filt_int = utils_c.interp_c(lam_spec, x_filt, y_filt)
            filt_norm = np.trapz(y_filt_int*flux_spec, lam_spec) / np.trapz(y_filt_int, lam_spec)
            flux_spec /= filt_norm
            if verbose:
                t1 = time.time(); dt = t1-t0; t0=t1
                print 'Input spectrum is provided, normalize it for filter, %s (%.2f s)' %(self.filter, dt)

        #### Define the grism dispersion for the central pixel of an object.  
        #### Assume doesn't change across the object extent            
        t0 = time.time()
        orders, self.xi = unicorn.reduce.grism_model(xc, yc, lam_spec=lam_spec, flux_spec=flux_spec, BEAMS=BEAMS, grow_factor=self.grow_factor, pad=self.pad)
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Define grism dispersion (%.2f s).' %(dt)
            
        yord, xord = np.indices(orders.shape)
        beams = np.dot(np.ones((orders.shape[0],1), dtype=np.int), self.xi[5].reshape((1,-1)))
        
        non_zero = orders > 0
        #xord, yord, ford, word, sord = xord[non_zero], yord[non_zero], orders[non_zero], self.xi[2][non_zero], self.xi[3][non_zero]
        cast = np.uint
        xord, yord, ford, word, sord, bord = np.array(xord[non_zero], dtype=cast), np.array(yord[non_zero], dtype=cast), np.array(orders[non_zero], dtype=np.float64), self.xi[2][non_zero], self.xi[3][non_zero], beams[non_zero]
                
        ys = orders.shape
        xord += self.xi[0]
        yord -= (ys[0]-1)/2
        #
        self.object *= 0.
    
        # mask = self.segm[0].data == id
        # xpix = self.xf[mask]
        # ypix = self.yf[mask]
        # 
        # #### Loop through pixels defined within the segmentation region, summing
        # #### dispersed flux
        # for jj in range(xpix.size):
        #     x, y = xpix[jj], ypix[jj]
        #     xxi = x+xord
        #     yyi = y+yord
        #     use = (xxi >= 0) & (xxi < self.sh[1]) & (yyi >= 0) & (yyi < self.sh[0])
        #     self.object[yyi[use], xxi[use]] += ford[use]*self.flux[y,x]*10
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
        use = use &  (bord == 1)
        cut = np.zeros(self.sh[1])
        cut[xxi[use]] = word[use] 
        self.object_wave = cut*1. #np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
        
        cut[xxi[use]] = sord[use] 
        self.object_sens = cut*1. #np.dot(np.ones((self.sh[0],1)), cut.reshape(1,self.sh[1]))
        
        if verbose:
            t1 = time.time(); dt = t1-t0; t0=t1
            print 'Output wavelength and sensitivity arrays. (%.3f)' %(dt)
        
    def refine_model(self, ID, BEAMS=['A','B','C','D','E'], view=False, model_has_object=True, verbose=False):
        from scipy import polyfit,polyval

        if ID not in self.cat.id:
            print '#%d not in the catalog.' %(ID)
            return False
                    
        #### First iteration
        #self.compute_object_model(id=ID, BEAMS=['A'], lam_spec=self.lam_specs[ID], flux_spec = self.flux_specs[ID], verbose=verbose)
        self.compute_object_model(id=ID, BEAMS=BEAMS, lam_spec=self.lam_specs[ID], flux_spec = self.flux_specs[ID], verbose=verbose, normalize=False)
        
        t0 = time.time()
        if view.__class__ == threedhst.dq.myDS9:
            first_model = self.object*1
        
        if model_has_object:
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
        ratio_extract = utils_c.get_model_ratio(self.object, self.model, self.gris[1].data)
        
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Make ratio  (%.3f)' %(dt)
        
        #wave = np.mean(self.object_wave, axis=0)
        wave = self.object_wave
        
        keep = (wave > 1.12e4) & (wave < 1.65e4) & (ratio_extract != 0)
        #print len(wave), len(wave[wave > 0]), wave.max(), wave[2064], len(ratio_extract)
        
        if len(keep[keep]) < 10:
            self.model += self.object
            return True
            
        coeff = polyfit(wave[keep], ratio_extract[keep], 1)
        ratio_extract[~keep] = polyval(coeff, wave[~keep])
        
        if self.flux_specs[ID] is not None:
            rint = utils_c.interp_c(self.lam_specs[ID], wave[wave > 0], ratio_extract[wave > 0])
            self.flux_specs[ID] *= rint
        else:
            self.lam_specs[ID] = wave[wave > 0]
            self.flux_specs[ID] = ratio_extract[wave > 0]
        #
        if verbose:
            t1 = time.time(); dt=t1-t0; t0=t1
            print 'Extract ratio, save values  (%.3f)' %(dt)
            
        self.compute_object_model(id=ID, BEAMS=BEAMS, lam_spec = self.lam_specs[ID], flux_spec = self.flux_specs[ID], normalize=False, verbose=verbose)
        if view.__class__ == threedhst.dq.myDS9:
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
            #view.set('cmap value 4.43936 0.462911')
            view.set('cmap HSV')
            view.set('cmap value 7.02222 0.0867647')
            view.set('regions file %s_inter.reg' %(self.root))
            view.match_all(alignType='image')
            
            plt.plot(self.lam_specs[ID], self.flux_specs[ID])
        
        ### Done, now update the total model
        if view.__class__ == threedhst.dq.myDS9:
            self.model += second_model
        else:
            self.model += self.object

    def show_ratio_spectrum(self, ID):
        if ID not in self.cat['NUMBER']:
            print '#%d not in the catalog.' %(ID)
            return False
        
        plt.plot(self.lam_specs[ID], self.flux_specs[ID]*self.total_fluxes[ID])
        
    def compute_full_model(self, BEAMS=['A','B','C','D','E'], save_ext='_model0', iteration_number=0, view=None, MAG_LIMIT=23.5):
        
        self.model*=0.
        threedhst.showMessage('Making full object model...')
        
        mag = self.cat['MAG_AUTO']
        so = np.argsort(mag)
        so = so[mag[so] <= MAG_LIMIT]
        N = len(so)
        
        for i in range(N):
            #for i, id in enumerate(self.cat['NUMBER'][so]):
            id = self.cat['NUMBER'][so][i]
            print noNewLine+'Object #%d, m=%.2f (%d/%d)' %(id, mag[so][i], i+1, N)
            self.refine_model(id, BEAMS=BEAMS, view=view, model_has_object=(iteration_number > 0))      
        
        if save_ext is not None:
            pyfits.writeto(self.root+save_ext+'.fits', self.model, clobber=True)
    
    def twod_spectrum(self, id=328, grow=1, miny=50, CONTAMINATING_MAGLIMIT=28, refine=True, verbose=False):
        
        seg = self.segm[0].data
        seg_mask = seg == id
        ypix, xpix = self.yf, self.xf
    
        NY = ypix[seg_mask].max()-ypix[seg_mask].min()
        NX = xpix[seg_mask].max()-xpix[seg_mask].min()
        NT = np.max([NY*grow, NX*grow, miny])
        
        ii = np.where(np.cast[int](self.cat.id) == id)[0][0]
        
        xc = np.int(np.round(self.cat.x_pix[ii]))-1
        yc = np.int(np.round(self.cat.y_pix[ii]))-1
        
        if verbose:
            print 'Generating direct image extensions'
            
        self.direct_thumb = self.im[1].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_wht = self.im[2].data[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2]
        self.direct_seg = np.cast[np.int16](seg[yc-NT/2:yc+NT/2, xc-NT/2:xc+NT/2])
        
        #### Initialize the output FITS file
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
        extensions.append(pyfits.ImageHDU(self.direct_wht, header))
        header.update('EXTNAME','DSEG')
        extensions.append(pyfits.ImageHDU(self.direct_seg, header))
        
        #### Find objects that could contaminate the current one
        dylo = self.cat.y_pix-self.cat.y_pix[ii]-self.cat['FLUX_RADIUS']*3
        dy0 = self.cat.y_pix-self.cat.y_pix[ii]
        dyhi = self.cat.y_pix-self.cat.y_pix[ii]+self.cat['FLUX_RADIUS']*3
        dy = np.minimum(np.abs(dylo), np.abs(dy0))
        dy = np.minimum(np.abs(dyhi), dy)
        
        dx = self.cat.x_pix-self.cat.x_pix[ii]
        nearby = (dy < (NT/2.+5)) & (dx < 420*self.grow_factor) & (self.cat.mag < CONTAMINATING_MAGLIMIT)
        
        BEAMS=['A','B','C','D','E']
        view=False
        if nearby.sum() > 0:
            ids = self.cat.id[nearby][np.argsort(self.cat.mag[nearby])]
            for idi in ids:
                if verbose:
                    print 'Calculating contaminating object #%d' %(idi)
                #
                if refine:
                    self.refine_model(idi, BEAMS=BEAMS, view=view, model_has_object=(self.lam_specs[idi] is not None))      
                else:
                    self.compute_object_model(idi, BEAMS=BEAMS, lam_spec=self.lam_specs[idi], flux_spec=self.flux_specs[idi], normalize=False)      
                                
        #### Generate model for current object
        if verbose:
            print 'Compute object model #%d' %(id)

        if refine:
            self.refine_model(id, BEAMS=BEAMS, view=view, model_has_object=(self.lam_specs[id] is not None))      
        else:
            self.compute_object_model(id, BEAMS=BEAMS, lam_spec=self.lam_specs[id], flux_spec=self.flux_specs[id], normalize=False)      
        
        #### Find pixels of the 1st order        
        xarr = np.arange(self.sh[0])
        xmin = xarr[(self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)].min()
        xmax = xarr[(self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)].max()
        
        #### dydx offset for grism trace
        xoff_arr = np.arange(len(self.xi[4]), dtype=np.double)+self.xi[0]+xc
        xoff_int = np.arange(xmin, xmax, dtype=np.double)
        yoff_int = utils_c.interp_c(xoff_int, xoff_arr, self.xi[4])

        yc = np.int(np.round(self.cat.y_pix[ii]+yoff_int.mean()))-1
        #yc_off = np.cast[int](np.round(self.cat.y_pix[ii]+yoff_int))-1
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
        
        self.grism_sci = self.gris[1].data[yc-NT/2:yc+NT/2, xmin:xmax]
        
        #### testing shifting the pixels:  better to just provide the (float) trace
        ####  than shifting by integer pixels
        
        # self.grism_sci_yoff = self.grism_sci*0.
        # for i in range(len(yoff_int)):
        #     self.grism_sci_yoff[:,i] = self.gris[1].data[yc_off[i]-NT/2:yc_off[i]+NT/2, xmin+i]

        header.update('EXTNAME','SCI')
        extensions.append(pyfits.ImageHDU(self.grism_sci, header))
        
        self.grism_wht = self.gris[2].data[yc-NT/2:yc+NT/2, xmin:xmax]
        header.update('EXTNAME','WHT')
        header.update('WHTERROR',self.gris[0].header['WHTERROR'])
        
        extensions.append(pyfits.ImageHDU(self.grism_wht, header))
        
        self.the_object = self.object[yc-NT/2:yc+NT/2, xmin:xmax]
        header.update('EXTNAME','MODEL')
        extensions.append(pyfits.ImageHDU(self.the_object, header))
        
        self.full_model = self.model[yc-NT/2:yc+NT/2, xmin:xmax]
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
        hdu.writeto(self.root+'_%05d.fits' %(id), clobber=True)
        
        self.twod = hdu
    
    def show_2d(self, savePNG = False):
        import unicorn
        
        ii = np.where(np.cast[int](self.cat.NUMBER) == self.id)[0][0]
        
        fig = unicorn.catalogs.plot_init(xs=5, left=0.05, right=0.05, bottom=0.05, top=0.05, aspect=1./1.2)
        
        xx = np.arange(self.grism_sci.shape[1])
        ltick = np.array([1.2,1.3,1.4,1.5,1.6])
        xtick = np.interp(ltick*1.e4, self.wave, xx)
        
        ax = fig.add_subplot(431)
        ax.imshow(self.direct_thumb, vmin=-0.5, vmax=2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        ax.text(1.1,0.8, self.root+r'_%05d' %(self.id), transform=ax.transAxes)
        ax.text(1.1,0.6, r'$m_{140}=%.2f$' %(self.cat.mag[ii]), transform=ax.transAxes)
        
        if self.lam_specs[self.id] is not None:
            ax = fig.add_subplot(422)
            ax.plot(self.lam_specs[self.id], self.flux_specs[self.id]*self.total_fluxes[self.id])
            xx = ax.set_xticks(ltick*1.e4)
            ax.set_xticklabels(ltick)
            ax.set_xlim(1.1e4,1.65e4)
            y = self.flux_specs[self.id]*self.total_fluxes[self.id]
            ax.set_ylim(-0.1*y.max(), 1.1*y.max())
            ax.set_yticklabels([]);
        
        ax = fig.add_subplot(412)
        ax.imshow(self.grism_sci, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        ax = fig.add_subplot(413)
        ax.imshow(self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        ax = fig.add_subplot(414)
        ax.imshow(self.grism_sci-self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)
        
        if savePNG:
            fig.savefig(self.root+'_%05d.png' %(self.id))
            plt.close()
            
    def twod_spectrum_old(self, id=328, grow=1, miny=80):
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
        header.update('WHTERROR',self.gris[0].header['WHTERROR'])
        
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
        yint = utils_c.interp_c(np.array([1.1157e4, 1.6536e4]), wave, sens)
        yint_smooth = utils_c.interp_c(np.array([1.1157e4, 1.6536e4]), wave, sens_smooth)
        
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
    
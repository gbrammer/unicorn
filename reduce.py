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
plt.rcParams['patch.edgecolor'] = 'None'
plt.rcParams['font.size'] = 10

plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.interpolation'] = 'nearest'

USE_PLOT_GUI = False

import pyfits
import pyraf
from pyraf import iraf

noNewLine = '\x1b[1A\x1b[1M'

import time

import threedhst

import unicorn.utils_c as utils_c

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

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
        ### Redo
        clean_all=False; clean_spectra=False; make_images=False; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
        extract_limit=23.5
        
        ### Fainter extraction limit, don't regenerate everything already done
        clean_all=False; clean_spectra=False; make_images=False; make_model=True; fix_wcs=True; extract_limit=None; skip_completed_spectra=True; MAG_LIMIT=26; out_path='./'
        extract_limit=25
        
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
            print noNewLine+'Remove %s' %(ff)
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
        model.trim_edge_objects()
        model.get_corrected_wcs()
        model.make_wcs_region_file()
        
    if extract_limit is not None:
        model.extract_spectra_and_diagnostics(MAG_LIMIT=extract_limit, skip=skip_completed_spectra)
    
    del(model)
    
    return True
                        
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
        
        ### Mask cosmic rays
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
        use = ((im[3].data & (4+32+16+2048+4096)) == 0) & (~hot_pix)
        
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
    
def grism_model(xc_full=244, yc_full=1244, lam_spec=None, flux_spec=None, grow_factor=2, pad = 60, BEAMS=['A','B','C','D'], dydx=True):
    
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
    
def process_GrismModel(root='GOODS-S-24', grow_factor=2, pad=60, MAG_LIMIT=24):
    import unicorn.reduce
    
    model = unicorn.reduce.GrismModel(root=root, grow_factor=grow_factor, pad=pad, MAG_LIMIT=MAG_LIMIT)
    
    model.get_corrected_wcs(verbose=True)
    
    test = model.load_model_spectra()
    if not test:
        ### First iteration with flat spectra and the object flux
        model.compute_full_model(refine=False, MAG_LIMIT=model.cat.mag.max(), save_pickle=False)   
        ### For the brighter galaxies, refine the model with the observed spectrum         
        model.compute_full_model(refine=True, MAG_LIMIT=21, save_pickle=True)
        
    return model
    
class Interlace1D():
    def __init__(self, file='UDS-17_00319.1D.fits', PNG=True):
        self.file = file
        tab = pyfits.open(file)
        self.header = tab[1].header
        self.data = tab[1].data
        tab.close()
        
        if PNG:
            self.show(savePNG=True)
            
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
            canvas.print_figure(self.file.replace('fits','png'), dpi=100, transparent=False)
            
class GrismModel():
    def __init__(self, root='GOODS-S-24', grow_factor=2, pad=60, MAG_LIMIT=24):
        """
        Initialize: set padding, growth factor, read input files.
        """
        self.root=root
        self.grow_factor = grow_factor
        self.pad = pad
        
        header = pyfits.getheader(root+'-F140W_inter.fits',0)
        self.filter = header['FILTER']
        ZPs = {'F125W':26.25, 'F140W':26.46, 'F160W':25.96}
        self.ZP = ZPs[self.filter]

        PLAMs = {'F125W':1.2486e4, 'F140W':1.3923e4, 'F160W': 1.5369e4}
        self.PLAM = PLAMs[self.filter]
        
        self.im = pyfits.open(self.root+'-F140W_inter.fits')
        self.make_flux()
                
        self.gris = pyfits.open(self.root+'-G141_inter.fits')
        self.gris[1].data = np.array(self.gris[1].data, dtype=np.double)
        self.gris[2].data = np.array(self.gris[2].data, dtype=np.double)
        self.sh = self.im[1].data.shape
        
        if not os.path.exists(root+'_seg.fits'):
            threedhst.showMessage('Running SExtractor...')
            self.find_objects(MAG_LIMIT=MAG_LIMIT)
                
        threedhst.showMessage('Read FITS files and catalog')
        print 'Read files'
        self.read_files()
        
        print 'Init object spectra'
        self.init_object_spectra()
        
        self.make_total_flux()
        
        self.ra_wcs, self.dec_wcs = None, None
        self.model = np.zeros(self.sh, dtype=np.double)
        self.object = np.zeros(self.sh, dtype=np.double)
        #self.test_object = np.zeros(self.sh)
        self.yf, self.xf = np.indices(self.sh)
                
    def read_files(self):
        """
        Read FITS files, catalogs, and segmentation images for a pair of 
        interlaced exposures.
        """        
        self.cat = threedhst.sex.mySexCat(self.root+'_inter.cat')
        
        self.trim_edge_objects(verbose=True)
        
        self.cat.x_pix = np.cast[float](self.cat.X_IMAGE)
        self.cat.y_pix = np.cast[float](self.cat.Y_IMAGE)
        self.cat.mag = np.cast[float](self.cat.MAG_AUTO)
        
        self.segm = pyfits.open(self.root+'_seg.fits')
        self.segm[0].data = np.array(self.segm[0].data, dtype=np.uint)
                
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
        
        x_edge = (self.cat['X_IMAGE'] < self.pad/2+10) | (self.cat['X_IMAGE'] > (self.sh[1]-self.pad/2-25))
        y_edge = (self.cat['Y_IMAGE'] < self.pad/2+10) | (self.cat['Y_IMAGE'] > (self.sh[0]-self.pad/2-25))
        
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
        import iraf
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
        
        ### For some reason this dies occasionally "too many positional arguments for traxy"
        ### Try running a second time if it dies once
        try:
            status = iraf.tran(origimage=flt+'[sci,1]', drizimage=self.root+'-F140W_drz.fits[1]', direction="forward", x=None, y=None, xylist="/tmp/%s.flt_xy" %(self.root), mode="h", Stdout=1)
        except:
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
        
        self.lam_spec = np.arange(0.95e4,1.8e4,22.)
        
        self.obj_in_model = {}
        self.flux_specs = {}
        for obj in self.cat['NUMBER']:
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
        
            
    def compute_object_model(self, id=328, lam_spec=None, flux_spec=None, BEAMS=['A','B','C','D'], verbose=False, normalize=True):
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
        cast = np.int
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
        keep = (wave > 1.12e4) & (wave < 1.65e4) & (ratio_extract != 0) & (xpix > (self.pad-22)) & (xpix < (self.sh[1]-self.pad-22))
        #print len(wave), len(wave[wave > 0]), wave.max(), wave[2064], len(ratio_extract)
        
        if keep.sum() < 10:
            self.model += self.object
            self.obj_in_model[id] = True
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
            view.set('pan to %f %f' %(self.cat['X_IMAGE'][idx][0], self.cat['Y_IMAGE'][idx][0]))
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

    def show_ratio_spectrum(self, id, flux=True):
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
            
    def compute_full_model(self, BEAMS=['A','B','C','D'], view=None, MAG_LIMIT=21., save_pickle=True, refine=True):
        
        #self.model*=0.
        if refine:
            threedhst.showMessage('Making full object model (refined model)')
        else:
            threedhst.showMessage('Making full object model (flat spectra)')
            
        mag = self.cat['MAG_AUTO']
        so = np.argsort(mag)
        so = so[mag[so] <= MAG_LIMIT]
        N = len(so)
                
        for i in range(N):
            #for i, id in enumerate(self.cat['NUMBER'][so]):
            id = self.cat['NUMBER'][so][i]
            print noNewLine+'Object #%d, m=%.2f (%d/%d)' %(id, mag[so][i], i+1, N)
            if refine:
                self.refine_model(id, BEAMS=BEAMS, view=view)      
            else:
                self.compute_object_model(id, BEAMS=BEAMS, lam_spec=None, flux_spec=None)      
                self.model += self.object
                self.obj_in_model[id] = True
                
        if save_pickle:
            self.save_model_spectra()
            
    def save_model_spectra(self):
        import pickle
        
        fp = open(self.root+'_model.pkl','wb')
        pickle.dump(self.cat.id, fp)
        pickle.dump(self.obj_in_model, fp)
        pickle.dump(self.lam_spec, fp)
        pickle.dump(self.flux_specs, fp)
        fp.close()
        
        pyfits.writeto(self.root+'_model.fits', data=self.model, header=self.gris[1].header, clobber=True)
        
    def load_model_spectra(self):
        import pickle
        if not os.path.exists(self.root+'_model.pkl'):
            return False
            
        fp = open(self.root+'_model.pkl','rb')
        ids = pickle.load(fp)
        test = ids == self.cat.id
        if np.sum(test) == len(self.cat.id):
            print 'Load spectra from pickle'
            self.obj_in_model = pickle.load(fp)
            self.lam_spec = pickle.load(fp)
            self.flux_specs = pickle.load(fp)
            
            im = pyfits.open(self.root+'_model.fits')
            self.model = np.cast[np.double](im[0].data)
            im.close()
        else:
            print "ERROR: Object list in pickle doesn't match the current catalog."
            fp.close()
            return False
            
        fp.close()
        return True
        
    def twod_spectrum(self, id=328, grow=1, miny=50, CONTAMINATING_MAGLIMIT=23, refine=True, verbose=False, force_refine_nearby=False):
        
        seg = self.segm[0].data
        seg_mask = seg == id
        ypix, xpix = self.yf, self.xf
    
        NY = ypix[seg_mask].max()-ypix[seg_mask].min()
        NX = xpix[seg_mask].max()-xpix[seg_mask].min()
        NT = np.max([NY*grow, NX*grow, miny])
        
        #NT = np.min([self.pad/2, NT])
        
        ii = np.where(np.cast[int](self.cat.id) == id)[0][0]
        
        xc = np.int(np.round(self.cat.x_pix[ii]))-1
        yc = np.int(np.round(self.cat.y_pix[ii]))-1
        
        NT = np.min([NT, xc, yc, self.sh[1]-xc, self.sh[0]-yc])
        
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
            self.model += self.object
            self.obj_in_model[idi] = True
        
        #### Find pixels of the 1st order        
        xarr = np.arange(self.sh[0])
        wavelength_region = (self.object_wave >= 1.05e4) & (self.object_wave <= 1.70e4)
        
        if np.sum(wavelength_region) < 20:
            fp = open(self.root+'_%05d.2D.xxx' %(id),'w')
            fp.write('%.2f %.2f  %d %d\n' %(xc, yc, -1, -1))
            fp.close()
            return False
        
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
            print noNewLine+'Thumb panel'
        
        ax = fig.add_subplot(431)
        ax.imshow(self.direct_thumb, vmin=-0.5, vmax=2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        ax.text(1.1,0.8, self.root, transform=ax.transAxes)
        ax.text(1.1,0.6, '#%d' %(self.id), transform=ax.transAxes)
        ax.text(1.1,0.4, r'$m_{140}=%.2f$' %(self.cat.mag[ii]), transform=ax.transAxes)
        ax.text(1.1,0.2, r'$(%.1f, %.1f)$' %(self.cat.x_pix[ii], self.cat.y_pix[ii]), transform=ax.transAxes)
        
        if verbose:
            print noNewLine+'Plot flux_specs'
        
        if self.flux_specs[self.id] is not None:
            ax = fig.add_subplot(422)
            ax.plot(self.lam_spec, self.flux_specs[self.id]*self.total_fluxes[self.id])
            xx = ax.set_xticks(ltick*1.e4)
            ax.set_xticklabels(ltick)
            ax.set_xlim(1.1e4,1.65e4)
            y = self.flux_specs[self.id]*self.total_fluxes[self.id]
            ax.set_ylim(-0.1*y.max(), 1.1*y.max())
            ax.set_yticklabels([]);
        
        if verbose:
            print noNewLine+'GRISM - science'
        
        ax = fig.add_subplot(412)
        ax.imshow(self.grism_sci, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        if verbose:
            print noNewLine+'GRISM - model'

        ax = fig.add_subplot(413)
        ax.imshow(self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)

        if verbose:
            print noNewLine+'GRISM - difference'
        
        ax = fig.add_subplot(414)
        ax.imshow(self.grism_sci-self.full_model, vmin=-0.05, vmax=0.2, interpolation='nearest')
        ax.set_yticklabels([]); ax.set_xticklabels([]); 
        xx = ax.set_xticks(xtick)
        
        if verbose:
            print noNewLine+'Save the PNG'
        
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
        c7 = pyfits.Column(name='sensitivity', format='D', unit='COUNTS / 1E-17 CGS', array=self.sens)
        
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
        
    def extract_spectra_and_diagnostics(self, list=None, skip=True, MAG_LIMIT=19, verbose=False):
        """
        Extract 1D and 2D spectra of objects with id in `list`.
        """
        import unicorn.reduce
        if list is None:
            list = self.cat.id[self.cat.mag < MAG_LIMIT]
        
        #fp = open(self.root+'_inter.failed','w')
        
        N = len(list)
        for i,id in enumerate(list):
            print noNewLine+'Object #%-4d (%4d/%4d)' %(id,i+1,N)
            #
            if (os.path.exists(self.root+'_%05d.2D.fits' %(id)) | os.path.exists(self.root+'_%05d.2D.xxx' %(id))) & skip:
                continue
            
            if verbose:
                print '2D FITS'
            status = self.twod_spectrum(id=id, verbose=verbose, miny=30, refine=True)
            
            if status is False:
                if verbose:
                    print noNewLine+'Object #%-4d (%4d/%4d) - No 2D' %(id,i+1,N)
                continue
                    
            if verbose:
                print noNewLine+'2D FITS, 2D PNG, 1D FITS'
            self.show_2d(savePNG=True, verbose=verbose)

            if verbose:
                print noNewLine+'2D FITS, 2D PNG, 1D FITS, 1D PNG'
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
    
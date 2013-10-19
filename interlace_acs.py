"""
Combine dithered ACS exposures for use with 3D-HST interlacing software.

Model and fits look pretty good, the former with some artificial offsets to some of the grism orders.  Some hints
that redshifts *underestimated* from the ACS fits.

"""
import unicorn
import pyfits
import os
import matplotlib.pyplot
import threedhst.dq
import numpy as np

def interlace_combine_acs(root='jbhm19010', view=True, use_error=True, make_undistorted=False, pad = 120, NGROW=0, ddx=0, ddy=0, growx=1, growy=1, auto_offsets=False, chip=1, filter='F814W', outroot='UDS-19'):
    """
    Combine four dithered ACS exposures in an interlaced image, but use the same native ACS pixel grid
    since the dither offsets don't evenly sample the ACS pixel.  This also simplifies the image distortions
    to following the native definitions.
    
    Input is a visit ID which will read the attached ASN table.  FLT images should be destriped and CTE corrected
    and run through [Multi/Astro]Drizzle to flag CRs (DQ=4096).
    
    Provide `filter` and `outroot` keywords to rename the visit name to something more meaningful.
    
    Make separate images for each chip, note [SCI,1] is CCDCHIP=2, [SCI,2] is CCDCHIP=1.
    """
    # from pyraf import iraf
    # #from iraf import iraf
    # from iraf import dither
    # 
    import threedhst.prep_flt_files
    import unicorn.reduce as red
    
    chip_ext = {1:2,2:1} ### Chip1 is SCI,2 and Chip2 is SCI,1 extensions
    
    if unicorn.hostname().startswith('uni'):
        view = False
    
    if view:
        import threedhst.dq
        ds9 = threedhst.dq.myDS9()

    #
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    
    ### xx update for ACS PAM
    im = pyfits.open(os.getenv('jref')+'wfc%d_pam.fits' %(chip))
    PAM = im[0].data
    im.close()
    
    #PAM = flt[1].data*0.+1
    
    if os.path.exists(root+'.run'):
        run = threedhst.prep_flt_files.MultidrizzleRun(root)
        #run.blot_back(ii=0, copy_new=True)
        scl = np.float(run.scl)
        xsh, ysh = threedhst.utils.xyrot(np.array(run.xsh)*scl, np.array(run.ysh)*scl, run.rot[0])
    else:
        xsh, ysh = np.zeros(len(asn.exposures)), np.zeros(len(asn.exposures))
    
    NX, NY = 4096, 2048
    yi,xi = np.indices((NY, NX))
    
    #pad += NGROW*4
    
    N = np.zeros((NY*growy+pad+growy*2*NGROW, NX*growx+pad+growx*2*NGROW), dtype=np.int)
        
    inter_sci = np.zeros((NY*growy+pad+growy*2*NGROW, NX*growx+pad+growx*2*NGROW))
    inter_weight = np.zeros((NY*growy+pad+growy*2*NGROW, NX*growx+pad+growx*2*NGROW))
    # if use_error:
    #     inter_err = np.zeros((NY*growy+pad+growy*2*NGROW, NX*growx+pad+growx*2*NGROW))
    
    xi+=pad/(2*growx)+NGROW
    yi+=pad/(2*growy)+NGROW
    
    #### Interlaced offsets
    dxs, dys = unicorn.reduce.acs_interlace_offsets(None, growx=growx, growy=growy, path_to_flt='./')
    
    dxs += ddx
    dys += ddy
    
    #### Find hot pixels, which are flagged as cosmic 
    #### rays at the same location in each of the flt files.  Ignore others.
    hot_pix = np.zeros((NY, NX),dtype='int')
    
    for flt in asn.exposures:
        im = pyfits.open(flt+'_flt.fits')
        hot_pix += (im['DQ',chip_ext[chip]].data & 4096) / 4096
        
    hot_pix = hot_pix > 2
        
    for i,flt in enumerate(asn.exposures):
        print flt
        #flt = run.flt[i]
        im = pyfits.open(flt+'_flt.fits')
        #ds9.frame(i+1); ds9.view(im[1].data); ds9.scale(0,5)
        
        #### Use the pixel area map correction
        im['SCI',chip_ext[chip]].data *= PAM
        
        #### Divide by 4 to conserve surface brightness with smaller output pixels
        im['SCI',chip_ext[chip]].data /= 1.*(growx*growy)
        im['ERR',chip_ext[chip]].data /= 1.*(growx*growy)
        
        ### Mask cosmic rays
        if i == 0:
            h0 = im[0].header
            h1 = im['SCI',chip_ext[chip]].header
            #header = red.scale_header_wcs(h1.copy(), factor=2, growx=growx, growy=growy, pad=pad)
            header = h1.copy()
            header.update('EXTNAME','SCI')
            #header.update('PAD',pad)
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
        use = ((im['DQ',chip_ext[chip]].data & (4+32+16+2048+4096)) == 0) & (~hot_pix)
        
        #### preparation step doesn't subtract background
        im['SCI',chip_ext[chip]].data -= np.median(im['SCI',chip_ext[chip]].data[use])
        im['SCI',chip_ext[chip]].data /= im[0].header['EXPTIME']
        im['ERR',chip_ext[chip]].data /= im[0].header['EXPTIME']
        
        #inter_sci[yi[use]*growy+dy,xi[use]*growx+dx] += im['SCI',chip_ext[chip]].data[use]
        inter_sci[yi[use]*growy+dy,xi[use]*growx+dx] += im['SCI',chip_ext[chip]].data[use]/im['ERR',chip_ext[chip]].data[use]**2
        inter_weight[yi[use]*growy+dy,xi[use]*growx+dx] += 1./im['ERR',chip_ext[chip]].data[use]**2
        
        N[yi[use]*growy+dy,xi[use]*growx+dx] += 1
        
        # if use_error:
        #     inter_err[yi[use]*growy+dy,xi[use]*growx+dx] += im['ERR',chip_ext[chip]].data[use]**2
        
        if view:
            ds9.view_array(inter_sci/np.maximum(N,1), header=header)
            ds9.scale(-0.1,5)
    
    #### Average for case when dither positions overlap, e.g. CANDELS SN fields
    # inter_sci /= np.maximum(N,1) 
    # inter_err = np.sqrt(inter_err) / np.maximum(N, 1)
    # inter_err[N == 0] = 0.
    inter_weight[inter_weight == 0] = 1
    inter_sci = inter_sci / inter_weight
    inter_err = np.sqrt(1./inter_weight)
    inter_err[N == 0] = 0.
    
    if use_error:
        h0.update('WHTERROR',True,comment='WHT extension is FLT[err,1]')
    else:
        h0.update('WHTERROR',False,comment='WHT extension is 0/1 flagged bad pix')
    
    ### Clip back to native size
    inter_sci = inter_sci[pad/2:-pad/2, pad/2:-pad/2]
    inter_err = inter_err[pad/2:-pad/2, pad/2:-pad/2]
    N = N[pad/2:-pad/2, pad/2:-pad/2]
    #print pad, inter_sci.shape
    
    hdu = pyfits.PrimaryHDU(header=h0)
    sci = pyfits.ImageHDU(data=np.cast[np.float32](inter_sci), header=header)
    if use_error:
        wht = pyfits.ImageHDU(data=np.cast[np.float32](inter_err), header=header_wht)
    else:
        wht = pyfits.ImageHDU(data=np.cast[np.int](N), header=header_wht)
        
    image = pyfits.HDUList([hdu,sci,wht])
    if 'EXTEND' not in hdu.header.keys():
        hdu.header.update('EXTEND', True, after='NAXIS')
    
    image[0].header.update('FILTER', filter)
    if outroot is None:
        outroot = root
        
    image.writeto('%s-chip%d-%s_inter.fits' %(outroot, chip, filter), clobber=True)
    
    if 0:
        flt = pyfits.open('jbhj07yfq_flt.fits')
        inter = pyfits.open('GOODS-S-7-chip2-F814W_inter.fits')
        ds9.view(flt[1].data/flt[0].header['EXPTIME']-inter[1].data[pad:-pad,pad:-pad])
        
def subtract_grism_background(inter_image='jbhm19020_chip1_inter.fits', ORDER=3, SIGCLIP=3, verbose=True):
    """
    Fit a 2D polynomial to the grism background image since this isn't necessarily subtracted well
    in the preparation steps.
    """
    import scipy.linalg
    
    im = pyfits.open(inter_image, mode='update')
    sh = im[1].data.shape
    
    fit = threedhst.prep_flt_files.fit_2D_background(ORDER=ORDER, IMAGES=[], NX=sh[1], NY=sh[0])#, x0=507, y0=507)
    #fit.fit_image(im[1].data, A=fit.A, show=False, overwrite=False, save_fit=False)
    
    IMG = im[1].data
    IMGout = IMG*0
    
    #### First iteration
    if verbose:
        print 'First iteration.  No bad pix.'
        
    ok = (IMG != 0) & (IMG < 100)
    Aq = np.transpose(fit.A[:,ok])
    IMGq = IMG[ok].flatten()
    p, resid, rank, s = scipy.linalg.lstsq(Aq,IMGq)
    
    #### Second iteration, sigma clipped
    IMGout = IMG*0.
    for i in range(fit.NPARAM):
        IMGout += fit.A[i,:,:]*p[i]
    
    #
    if verbose:
        print 'Second iteration.  No bad pix, %f-sigma clipped.' %(SIGCLIP)
    
    ok = (IMG != 0) & (np.abs(IMG-IMGout)/im[2].data < 3)
    Aq = np.transpose(fit.A[:,ok])
    IMGq = IMG[ok].flatten()
    p, resid, rank, s = scipy.linalg.lstsq(Aq,IMGq)
    
    IMGout = IMG*0.
    for i in range(fit.NPARAM):
        IMGout += fit.A[i,:,:]*p[i]
    
    mask = IMG != 0
    IMG[mask] -= IMGout[mask]
    im.flush()

def test_model():
    """
    Test dataset in a GOODS-S pointing
    """    
    import unicorn.interlace_acs
    ### GOODS-S test
    root, new = 'jbhj07', 'GOODS-S-7'
    
    #### combine interlaced images in both chips
    for ext, filt in zip(['010','020'], ['F814W', 'G800L']):
        for chip in [1]:
            unicorn.interlace_acs.interlace_combine_acs('%s%s' %(root, ext), outroot=new, filter=filt, chip=chip)
            if filt == 'G800L':
                unicorn.interlace_acs.subtract_grism_background(inter_image='%s-chip%d-%s_inter.fits' %(new, chip, filt))
    
    #### Try blotting the reference segmentation image with the new astrodrizzle tools
    #### !!! segmentation images need to be "cleaned" images from unicorn.reduce.prepare_blot_reference
    
    threedhst.prep_flt_astrodrizzle.ablot_segmentation(segimage='GOODS-S_F125W_F140W_F160W.seg.fits', refimage='GOODS-S-7-chip1-F814W_inter.fits', output='blot_flt_seg.fits', use_coeffs=True)
    threedhst.prep_flt_astrodrizzle.ablot_segmentation(segimage='test_seg.fits', refimage='jbhj07yfq_flt.fits', output='blot_flt_seg.fits', use_coeffs=True)
    threedhst.prep_flt_astrodrizzle.ablot_segmentation(segimage='XXX_seg.fits', refimage='GOODS-S-7-chip1-F814W_inter.fits', output='blot_flt_seg.fits', use_coeffs=True)
    
    #### Make full multibeam models (not using reference segmentation images from above)
    for chip in [1,2]:    
        model = unicorn.reduce.GrismModel('%s-chip%d' %(new, chip), grow_factor=1, growx=1, growy=1, MAG_LIMIT=28, use_segm=False, grism='G800L', direct='F814W')    
        model.compute_full_model(MAG_LIMIT=26, refine=False, BEAMS=['A','B','C','D','E','F'])
    
    #### Extract a spectrum
    id=170; model.twod_spectrum(id); gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(model.root, id), dr_match=2, lowz_thresh=0.1)
    
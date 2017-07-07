"""
Try "drizzling" a region of grism FLTs at a specific wavelength to an output stamp
"""

import os
import glob

import numpy as np
import astropy.io.fits as pyfits

from drizzlepac import astrodrizzle
import stwcs

def diff(xarr):
    """
    Like np diff but make same size as input array filling first element
    with diff[0]
    """    
    import numpy as np
    d = np.diff(xarr)
    return np.append(d[0], d)


def get_drizzled_cutouts():
    
    ### Optinally display products to ds9 with, e.g., pysao.ds9
    try:
        import pysao
        ds9 = pysao.ds9()
    
        ds9.set('scale limits -0.08 8')
        ds9.set('lock colorbar')
        ds9.set('frame lock wcs')
        ds9.set('tile')
        for f in range(6):
            ds9.frame(f+1)

    except:
        ds9 = None
        
    
    ### Full mosaic, also provides a reference WCS
    im_mosaic = pyfits.open('../../MACS1149/Catalog/MACS1149-F160W_drz_sci.fits')
    wcs_mosaic = stwcs.wcsutil.HSTWCS(im_mosaic, ext=0)
    
    ### Get cutout of NX pixels around the center position from the reference mosaic
    NX, pix_scale = 50, 0.065
    wcs_mosaic.updatePscale(pix_scale)

    line_wavelengths = {'O2':3727, 'Ne3':3869, 'Hb':4861, 'O3':5007, 'Ha':6563., 'S2':6724}
    
    ### Object    
    stats = {'id':404, 'z':1.40707, 'lines':{'G141':['O3','Ha','S2'], 'G102':['O2']}}
    stats = {'id':1422, 'z':2.27763, 'lines':{'G141':['O2','Hb','O3']}}
    stats = {'id':1917, 'z':1.89105, 'lines':{'G141':['Hb','O3'], 'G102':['O2']}}

    stats = {'id':2315, 'z':1.8936, 'lines':{'G141':['Hb','O3'], 'G102':['O2']}}
    stats = {'id':2389, 'z':1.8936, 'lines':{'G141':['Hb','O3'], 'G102':['O2']}}

    stats = {'id':3746, 'z':1.2477, 'lines':{'G141':['O3','Ha','S2'], 'G102':['O2','Hb','O3']}}
    
    import collections
    drizzled = collections.OrderedDict()
    for grism in stats['lines'].keys():
        files=glob.glob('MACS1149-???-%s_%05d*2D.fits' %(grism, stats['id']))
        for line in stats['lines'][grism]:
            lam = line_wavelengths[line]*(1+stats['z'])
            drizzled['%s_%s' %(grism, line)] = driz_from_twod(files, lam=lam, pixfrac=0.5, wcs_mosaic=wcs_mosaic, ds9=ds9)
    
    ### Show all
    keys = drizzled.keys()
    for i, key in enumerate(keys):
        print 'Frame %d: %s' %(i+1, key)
        ds9.frame(i+1)
        ds9.view(drizzled[key][1], header=drizzled[key][2])
        
def driz_from_twod(files, lam=1.4e4, pixfrac=0.5, wcs_mosaic=None, ds9=None):
    """
    Drizzle from 2D cutout spectra
    """
    import os
            
    ### Drizzle 2D spectra!
    for i, file in enumerate(files):
        
        ### Open the 2D spectrum
        twod = pyfits.open(file)
        shg = twod['SCI'].data.shape
        shd = twod['DSCI'].data.shape
        
        ### If "unicorn" scripts available for computing 2D models
        ### will still oversubtract continuum if there are bright lines
        try:
            import unicorn
            gris = unicorn.reduce.Interlace2D(file)
            gris.compute_mode()
            twod['MODEL'].data = gris.model*1.
        except:
            pass
            
        #### (distorted) WCS of the cutout
        wcs_grism = stwcs.wcsutil.HSTWCS(twod, ext=1)
        
        ### Initialize if first spectrum
        if i==0:
            
            ### Output WCS
            xy = np.array(wcs_mosaic.all_world2pix([twod[0].header['RA']], [twod[0].header['DEC']], 0)).flatten()
            xy0 = np.cast[int](xy)
            
            if wcs_mosaic is None:
                wcs_mosaic = wcs_grism.copy()
                
            slx_mosaic = slice(xy0[0]-NX,xy0[0]+NX)
            sly_mosaic = slice(xy0[1]-NX,xy0[1]+NX)
            
            wcs_mosaic_slice = wcs_mosaic.slice((sly_mosaic, slx_mosaic))
            
            out_head = wcs_mosaic_slice.to_header()
            for ii in [1,2]:
                for jj in [1,2]:
                    if 'PC%d_%d' %(ii,jj) in out_head:
                        out_head['CD%d_%d' %(ii,jj)] = out_head['PC%d_%d' %(ii,jj)]
                        out_head.remove('PC%d_%d' %(ii,jj))
            
            ### Initialize drizzle products
            dsci = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            dwht = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            dcon = np.zeros((2*NX, 2*NX), dtype=np.int32)        

            gsci = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            gwht = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            gcon = np.zeros((2*NX, 2*NX), dtype=np.int32)        
                    
        ### Drizzle direct thumbnail
        astrodrizzle.adrizzle.do_driz(twod['DINTER'].data, wcs_grism, twod['DWHT'].data, wcs_mosaic_slice, dsci, dwht, dcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=np.nan, stepsize=10, wcsmap=None)    
        
        #############
        ### Drizzle 2D spectrum
        
        ### Center pixel along the trace at the desired wavelength
        x0 = np.interp(lam, twod['WAVE'].data, np.arange(shg[1])) - shd[1]/2
        y0 = np.interp(lam, twod['WAVE'].data, twod['YTRACE'].data) - shd[0]/2
        
        ### Drizzle the contamination & model-subtracted spectrum.  Here, the ['MODEL'] extension is likely only 
        ### a rough approxtwodation, since it won't include strong emission lines
        clean = twod['SCI'].data-twod['CONTAM'].data-twod['MODEL'].data
        
        ### Units of line flux (1e-17 erg/s/cm2)
        clean_flux = clean / twod['SENS'].data * diff(twod['WAVE'].data)
        
        ### Shift the 2D spectrum so that it should line up on top of the direct thumbnail, where the WCS is defined
        gris_thumb = np.roll(np.roll(clean_flux, -int(np.floor(x0)), axis=1), -int(np.floor(y0)), axis=0)
        
        ### Need the shifted uncertainty extension as well, which will indicate the masked pixels
        wht_thumb = np.cast[np.float32](1/np.roll(np.roll(twod['WHT'].data, -int(np.floor(x0)), axis=1), -int(np.floor(y0)), axis=0)**2)
        wht_thumb[~np.isfinite(wht_thumb)] = 0
        
        ## Slightly shifted WCS to account for pixel centering
        hdu = wcs_grism.to_fits(relax=True)
        hdu[0].data = gris_thumb*1
        hdu[0].header['CRPIX1'] += x0-np.floor(x0)
        hdu[0].header['CRPIX2'] += y0-np.floor(y0)
        wcs_grism_spec = stwcs.wcsutil.HSTWCS(hdu, ext=0)
        
        ## Drizzle it
        astrodrizzle.adrizzle.do_driz(gris_thumb, wcs_grism_spec, wht_thumb, wcs_mosaic_slice, gsci, gwht, gcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=np.nan, stepsize=10, wcsmap=None)    
        
        if ds9:
            ds9.frame(i+1)
            ds9.view(twod['DINTER'].data, header=twod['DINTER'].header)
            
            ds9.frame(5)
            ds9.view(dsci/wcs_grism.pscale**2, header=out_head)
            
            ds9.frame(6)
            ds9.view(gsci/wcs_grism_spec.pscale**2, header=out_head)
            
        #time.sleep(1)
    
    dsci /= wcs_grism.pscale**2
    gsci /= wcs_grism_spec.pscale**2
    
    return dsci, gsci, out_head
    
    #ds9.frame(6)
    #ds9.view(im_mosaic[0].data[sly_mosaic, slx_mosaic], header=out_h)
    
    
########## This all testing below

def test_driz_from_twod():
    
    import os
    import astropy.io.fits as pyfits
    
    from drizzlepac import astrodrizzle
    import stwcs
    
    ### Checking
    flt = pyfits.open('../PREP_Apr15/ica521naq_flt.fits')
    inter = pyfits.open('../PREP_Apr15/MACS1149-032-F160W_inter.fits')

    h_grow = unicorn.reduce.scale_header_wcs(flt[1].header.copy(), factor=2, growx=2, growy=2, pad=60, NGROWX=200, NGROWY=50)
    out = inter[1].data*1
    
    # ### Try dummy interlaced
    # flt_grow = np.zeros((2028,2028))
    # for i in [0,1]:
    #     for j in [0,1]:
    #         flt_grow[i::2,j::2] += flt[1].data
    # 
    # flt_wcs = stwcs.wcsutil.HSTWCS(flt, 1)
    # flt_wcs_grow = stwcs.wcsutil.HSTWCS(flt, 1)
    # h_grow = unicorn.reduce.scale_header_wcs(flt[1].header.copy(), factor=1, growx=2, growy=2, pad=0, NGROWX=0, NGROWY=0)
    ### try just CDELT
    #h_grow = flt[1].header.copy()
    #h_grow['CDELT1'] /= 2
    #h_grow['CDELT2'] /= 2
    
    hdu_grow = pyfits.HDUList([pyfits.ImageHDU(data=out, header=h_grow)])
    grow_wcs = stwcs.wcsutil.HSTWCS(hdu_grow, ext=0)
    
    sk = 50
    yp, xp = np.indices((1014/sk, 1014/sk))
    rd = flt_wcs.all_pix2world(xp.flatten()*sk, yp.flatten()*sk, 0)
    xp_grow, yp_grow = grow_wcs.all_world2pix(rd[0], rd[1], 0)

    ds9.frame(1)
    ds9.view(flt[1].data, header=flt[1].header)

    ds9.frame(2)
    ds9.view(inter[1].data, header=inter[1].header)
    
    
    #ds9.frame(4)
    ds9.view(inter[1].data, header=header)
    
    files=glob.glob('*model.fits')
    models = []
    for file in files:
        im = pyfits.open(file.replace('_model',''))
        grism = im[0].header['FILTER']
        direct = {'G102':'F105W', 'G141':'F160W'}[grism]
        model = unicorn.reduce.process_GrismModel(root=file.split('-%s' %(grism))[0], grism=grism, direct=direct)
        
        #### Fix headers with updated "scale_header_wcs"
        root = file.split('_')[0]
        asn = threedhst.utils.ASNFile('../PREP_Apr15/%s_asn.fits' %(root))
        flt = pyfits.open('../PREP_Apr15/%s_flt.fits' %(asn.exposures[0]))
        h = model.gris[1].header
        growx, growy, pad = h['GROWX'], h['GROWY'], h['PAD']
        ngrowx, ngrowy = h['NGROWX'], h['NGROWY']
        
        header = unicorn.reduce.scale_header_wcs(flt[1].header.copy(), factor=1, growx=growx, growy=growy, pad=pad, NGROWX=ngrowx, NGROWY=ngrowy)
        model.gris[1].header.update(header)
        model.gris_wcs = stwcs.wcsutil.HSTWCS(model.gris, ext=1)
        
        models.append(model)
        
    # id=3746
    # lam = 1.1275e4
    # 
    # im_mosaic = pyfits.open('../MACS1149/Catalog/MACS1149-F160W_drz_sci.fits')
    # wcs_mosaic = stwcs.wcsutil.HSTWCS(im_mosaic, ext=0)
    # 
    # N = 10
    # 
    # ### re-extract 2D
    if len(glob.glob('*%05d.2D.fits' %(id))) < 100:
        for model in models:
            model.twod_spectrum(id, miny=-80, refine=False)
        
    files=glob.glob('MACS1149-???-%s_%05d*2D.fits' %(grism, id))

    #files=glob.glob('MACS1149-???-G102_%05d*2D.fits' %(id))

    for i, file in enumerate(files):
        im = pyfits.open(file)
        shg = im['SCI'].data.shape
        shd = im['DSCI'].data.shape
        
        if i==0:
            NX = 50
            xy = np.array(wcs_mosaic.all_world2pix([im[0].header['RA']], [im[0].header['DEC']], 0)).flatten()
            xy0 = np.cast[int](xy)
            
            slx_mosaic = slice(xy0[0]-NX,xy0[0]+NX)
            sly_mosaic = slice(xy0[1]-NX,xy0[1]+NX)
            
            wcs_mosaic_slice = wcs_mosaic.slice((sly_mosaic, slx_mosaic))
            out_h = wcs_mosaic_slice.to_header()
            
            outsci = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            outwht = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            outcon = np.zeros((2*NX, 2*NX), dtype=np.int32)        

            gsci = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            gwht = np.zeros((2*NX, 2*NX), dtype=np.float32)        
            gcon = np.zeros((2*NX, 2*NX), dtype=np.int32)        
            
        root = file.split('_')[0]
        inter_grism = pyfits.open('../PREP_Apr15/%s_inter.fits' %(root))
        inter_ref_grism = pyfits.open('../PREP_Apr15/%s_ref_inter.fits' %(root))
        #wcs_grism = stwcs.wcsutil.HSTWCS(inter_grism, ext=1)
        
        #### Drizzle individual FLTs
        #asn = threedhst.utils.ASNFile('../PREP_Apr15/%s-F160W_asn.fits' %(root.split('-G1')[0]))
        # for exp in asn.exposures:
        #     print exp
        #     flt = pyfits.open('../PREP_Apr15/%s_flt.fits' %(exp))
        #     wcs_flt = stwcs.wcsutil.HSTWCS(flt, ext=1)
        #     wht = 1/flt['ERR'].data**2
        #     wht[flt['DQ'].data > 0] = 0
        #     wht[flt['SCI'].data/flt['ERR'].data < -3] = 0
        #     
        #     astrodrizzle.adrizzle.do_driz(flt['SCI'].data, wcs_flt, wht, wcs_mosaic_slice, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=0.5, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
        #     ds9.view(outsci)
            
        #### Get fixed interlaced header
        asn = threedhst.utils.ASNFile('../PREP_Apr15/%s_asn.fits' %(root))
        flt = pyfits.open('../PREP_Apr15/%s_flt.fits' %(asn.exposures[0]))
        print asn.exposures[0]
        header = unicorn.reduce.scale_header_wcs(flt[1].header.copy(), factor=1, growx=2, growy=2, pad=60, NGROWX=200, NGROWY=50)
        inter_grism[1].header = header.copy()
        wcs_grism = stwcs.wcsutil.HSTWCS(inter_grism, ext=1)
        
        h = im[0].header
        slx = slice(h['XDIRECT0'], h['XDIRECT0']+shd[1])
        sly = slice(h['YDIRECT0'], h['YDIRECT0']+shd[0])
        
        wcs_grism_slice = wcs_grism.slice((sly, slx))
        grism_cutout = inter_ref_grism[1].data[sly, slx]*1.
        #ds9.view(grism_cutout, header=wcs_grism_slice.to_header(relax=True))
        ds9.frame(i+1)
        ds9.view(im['DINTER'].data, header=wcs_grism_slice.to_header(relax=True))
        
        #wcs_grism_slice.wcs.crpix += dx #np.array(shd)/2.-dx
        wcs_grism_slice.naxis1 = im['DINTER'].header['NAXIS1']
        wcs_grism_slice.naxis2 = im['DINTER'].header['NAXIS2']
        
        hdu = wcs_grism_slice.to_fits(relax=True)
        hdu[0].data = im['DINTER'].data*1
        wcs_grism_slice2 = stwcs.wcsutil.HSTWCS(hdu, ext=0)
        
        #astrodrizzle.adrizzle.do_driz(im['DINTER'].data, wcs_grism_slice, im['DWHT'].data, wcs_mosaic_slice, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=0.1, kernel='square', fillval=np.nan, stepsize=10, wcsmap=wcs_functions.WCSMap)    
        astrodrizzle.adrizzle.do_driz(hdu[0].data, wcs_grism_slice2, im['DWHT'].data, wcs_mosaic_slice, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=0.1, kernel='square', fillval=np.nan, stepsize=10, wcsmap=wcs_functions.WCSMap)    
        #astrodrizzle.adrizzle.do_driz(inter_ref_grism[1].data*1, wcs_grism, np.ones_like(inter_ref_grism[1].data, dtype=np.float32), wcs_mosaic_slice, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=0., kernel='square', fillval=np.nan, stepsize=10, wcsmap=None)    
        
        ### Drizzle the spectrum
        x0 = np.interp(lam, im['WAVE'].data, np.arange(shg[1]))-shd[1]/2
        y0 = np.interp(lam, im['WAVE'].data, im['YTRACE'].data)-shd[0]/2-dy
        
        ### contam and model subtracted
        clean = im['SCI'].data-im['CONTAM'].data-im['MODEL'].data
        clean_flux = clean / im['SENS'].data * threedhst.utils.diff(im['WAVE'].data)
        
        gris_thumb = np.roll(np.roll(clean_flux, -int(np.floor(x0)), axis=1), -int(np.floor(y0)), axis=0)
        wht_thumb = np.cast[np.float32](1/np.roll(np.roll(im['WHT'].data, -int(np.floor(x0)), axis=1), -int(np.floor(y0)), axis=0)**2)
        wht_thumb[~np.isfinite(wht_thumb)] = 0
        
        hdu = wcs_grism_slice.to_fits(relax=True)
        hdu[0].data = gris_thumb*1
        hdu[0].header['CRPIX1'] += x0-np.floor(x0)
        hdu[0].header['CRPIX2'] += y0-np.floor(y0)
        wcs_grism_slice3 = stwcs.wcsutil.HSTWCS(hdu, ext=0)
        
        astrodrizzle.adrizzle.do_driz(gris_thumb, wcs_grism_slice3, wht_thumb, wcs_mosaic_slice, gsci, gwht, gcon, 1., 'cps', 1, wcslin_pscale=1, uniqid=1, pixfrac=0.5, kernel='square', fillval=np.nan, stepsize=10, wcsmap=wcs_functions.WCSMap)    
        
        # out = wcs_functions.WCSMap(wcs_grism_slice, wcs_mosaic_slice)
        # xpo, ypo = out.forward(xp, yp)
        # plt.scatter(xpo, ypo, alpha=0.5)
        
        ds9.frame(5)
        ds9.view(gsci/wcs_grism_slice3.pscale**2, header=out_h)
        #time.sleep(1)
    
    gsci /= wcs_grism_slice3.pscale**2
    outsci /= wcs_grism_slice2.pscale**2
    
    ds9.frame(6)
    ds9.view(im_mosaic[0].data[sly_mosaic, slx_mosaic], header=out_h)

def test():
    
    import numpy as np
    
    from drizzlepac import astrodrizzle
    import stwcs
    
    from threedhst import catIO
    import research.hawki_driz
    import unicorn
    import unicorn.interlace_test
    
    import os
    import astropy.io.fits as pyfits
    
    grism, filter = 'G102', 'F105W'
    
    main_root = 'FIGS-GS1-G102'
    
    #flt_files = glob.glob('icoi1*flt.fits')
    
    id=28445 # two knots
    id = 29649 # edge-on
    id = 30691 # fuzzy
    
    lrest = 6563.
    z = None
    
    size=3
    pixscale=0.05
    pixfrac = 0.4
    
    #twod_files = glob.glob('FIGS-GS1-*-G102_%05d.2D.fits' %(id))
    #twod = unicorn.reduce.Interlace2D(twod_files[0])
    #wave = 10625.5
    gris = unicorn.interlace_test.SimultaneousFit('%s_%05d' %(main_root, id), fast=True)
    if z is None:
        wave = lrest*(1+gris.z_max_spec)
    else:
        wave = lrest*(1+z)
    
    print 'Continue? (wave=%.1f, lrest=%.1f)' %(wave, lrest)
    var = raw_input()
    if var != '':
        print breakme
    else:
        print 'OK!'
            
    ra, dec = gris.twod.im[0].header['RA'], gris.twod.im[0].header['DEC']
    
    #ra, dec = 53.158163486, -27.7810907
    
    hout, wcs_out = research.hawki_driz.get_wcsheader(ra=ra, dec=dec, size=size, pixscale=pixscale, get_hdu=False, theta=0)
    sh = (int(hout['NAXIS2']), int(hout['NAXIS1']))
    
    files = catIO.Table('flt_files')
    flt_files = files['FILE'][files['FILTER'] == filter]
        
    outsci = np.zeros(sh, dtype=np.float32)
    outwht = np.zeros(sh, dtype=np.float32)
    outcon = np.zeros(sh, dtype=np.float32)
    
    ds9.frame(1)
    
    for file in flt_files:
        print file
        flt = pyfits.open(file)
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, 1)
        
        # xflt, yflt = flt_wcs.all_world2pix([ra], [dec], 0)
        # x0, y0 = int(xflt[0]), int(yflt[0])
        # sub_flt = flt[1].data[y0-11:y0+10, x0-30:x0+30]
        # ds9.view(sub_flt)
        
        sci = flt[1].data
        wht = 1/flt['ERR'].data**2
        wht = np.ones((1014,1014), dtype=np.float32)*np.median(wht)
        dq = flt['DQ'].data.copy()
        for bit in [64,512]:
            dq -= dq & bit
            
        msk = (dq > 0) | (~np.isfinite(wht))
        wht[msk] = 0.
        sci[msk] = 0.
        
        astrodrizzle.adrizzle.do_driz(sci, flt_wcs, wht, wcs_out, outsci, outwht, outcon, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
        ds9.view(outsci, header=hout)
        
    ## Try drizzling spectra     
    flt_files = files['FILE'][files['FILTER'] == grism]
     
    import mywfc3.grism
    conf = mywfc3.grism.aXeConf('%s/CONF/%s.test27s.gbb.conf' %(os.getenv('THREEDHST'), grism))
     
    goutsci = np.zeros(sh, dtype=np.float32)
    goutwht = np.zeros(sh, dtype=np.float32)
    goutcon = np.zeros(sh, dtype=np.float32)
    
    xpix = np.arange(0,300)
    
    ds9.frame(2)
    
    for file in flt_files:
        flt = pyfits.open(file)
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, 1)
    
        sci = flt[1].data
        wht = 1/flt['ERR'].data**2
        wht = np.ones((1014,1014), dtype=np.float32)*np.median(wht)
        dq = flt['DQ'].data.copy()
        for bit in [64,512]:
         dq -= dq & bit
        #
        msk = (dq > 0) | (~np.isfinite(wht))
        wht[msk] = 0.
        sci[msk] = 0.
    
        ### offset for wavelength
        xflt, yflt = flt_wcs.all_world2pix([ra], [dec], 0)
        ytrace, lam = conf.get_beam_trace(x=xflt[0], y=yflt[0], dx=xpix, beam='A')
        dx = np.interp(wave, lam, xpix)
        dy = np.interp(wave, lam, ytrace)
        print '%s (%.1f,%.1f), (%.2f,%.2f)' %(file, xflt[0], yflt[0], dx, dy)
        flt_wcs.wcs.crpix += np.array([dx, dy])
        try:
            flt_wcs.sip.crpix += np.array([dx, dy])
        except:
            pass
            
        x0, y0 = int(xflt[0]+dx), int(yflt[0]+dy)
        sub_flt = flt[1].data[y0-30:y0+30, x0-50:x0+30]
        prof = np.median(sub_flt, axis=1)
        yprof = np.zeros(1014)
        yprof[y0-30:y0+30] += prof
        
        sci = (sci.T-yprof).T
        #ds9.view(sub_flt)
        
        #dx = int(np.round(dx))
        #dy = int(np.round(dy))
        #sci = 
    
        astrodrizzle.adrizzle.do_driz(sci, flt_wcs, wht, wcs_out, goutsci, goutwht, goutcon, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
        ds9.view(goutsci, header=hout)
    
    ### Get continuum extraction from the 2D files
    import mywfc3.grism
    conf = mywfc3.grism.aXeConf('%s/CONF/%s.test27s.gbb.conf' %(os.getenv('THREEDHST'), grism))
     
    xpix = np.arange(0,300)
    
    ds9.frame(3)
    
    im = pyfits.open(os.getenv('iref')+'ir_wfc3_map.fits')
    PAM = im[1].data
    im.close()
    
    drizzled_sci = {}
    drizzled_wht = {}
    drizzled_ctx = {}
    drizzled_mod = {}
    drizzled_contam = {}
    
    twod_files=glob.glob('*-???-%s_%05d.2D.fits' %(grism, id))
    for twod_file in twod_files:
        twod = unicorn.reduce.Interlace2D(twod_file)
        twod.compute_model()
        
        root = twod_file.split('_')[0]
        print 'Interlaced: %s' %(root)
        inter = pyfits.open('%s_inter.fits' %(root))
        
        sh2 = twod.model.shape
        h = twod.im[0].header
        sly = slice(h['YGRISM0'], h['YGRISM0']+sh2[0])
        slx = slice(h['XGRISM0'], h['XGRISM0']+sh2[1])
        
        ### continuum subtraction
        clean = inter[1].data*0
        clean[sly, slx] += twod.model #+ twod.im['CONTAM'].data

        contam = inter[1].data*0
        contam[sly, slx] += twod.im['CONTAM'].data
        
        ### Get FLT components and "de-interlace" model
        asn = threedhst.utils.ASNFile('%s_asn.fits' %(root))
        
        if root.split('-G10')[0] in grism_ref_exp.keys():
            ref_exp = grism_ref_exp[root.split('-G10')[0]]
        else:
            ref_exp = 0
            
        dxs, dys = unicorn.reduce.wcs_interlace_offsets('%s_asn.fits' %(root), growx=2, growy=2, path_to_flt='./', verbose=False, ref_exp=ref_exp, raw=False, reference_pixel=None)
        
        hi = inter[1].header
        x0 = hi['PAD']/2 + hi['NGROWX']*hi['GROWX']
        y0 = hi['PAD']/2 + hi['NGROWY']*hi['GROWY']
        
        foutsci = np.zeros(sh, dtype=np.float32)
        foutwht = np.zeros(sh, dtype=np.float32)
        foutctx = np.zeros(sh, dtype=np.float32)

        coutsci = np.zeros(sh, dtype=np.float32)
        coutwht = np.zeros(sh, dtype=np.float32)
        coutctx = np.zeros(sh, dtype=np.float32)

        moutsci = np.zeros(sh, dtype=np.float32)
        moutwht = np.zeros(sh, dtype=np.float32)
        moutctx = np.zeros(sh, dtype=np.float32)
            
        for i,exp in enumerate(asn.exposures):
            #print exp
            flt = pyfits.open(exp+'_flt.fits')
            islx = slice(x0+dxs[i], x0+dxs[i]+1014*hi['GROWX'], hi['GROWX'])
            isly = slice(y0+dys[i], y0+dys[i]+1014*hi['GROWY'], hi['GROWY'])
        
            flt_wcs = stwcs.wcsutil.HSTWCS(flt, 1)

            sci = flt[1].data #- cutout
            cutout = clean[isly, islx]*(hi['GROWX']*hi['GROWY'])/PAM
            contam_cutout = contam[isly, islx]*(hi['GROWX']*hi['GROWY'])

            #wht = 1/flt['ERR'].data**2
            #wht = np.ones((1014,1014), dtype=np.float32)*np.median(wht)
            err = (np.median(flt['ERR'].data)**2 + contam_cutout**2*2)
            wht = 1/err
            dq = flt['DQ'].data.copy()
            for bit in [64,512]:
             dq -= dq & bit
            #
            msk = (dq > 0) | (~np.isfinite(wht))
            wht[msk] = 0.
            sci[msk] = 0.

            ### offset for wavelength
            xflt, yflt = flt_wcs.all_world2pix([ra], [dec], 0)
            ytrace, lam = conf.get_beam_trace(x=xflt[0], y=yflt[0], dx=xpix, beam='A')
            dx = np.interp(wave, lam, xpix)
            dy = np.interp(wave, lam, ytrace)
            print '%s (%.1f,%.1f), (%.2f,%.2f)' %(exp, xflt[0], yflt[0], dx, dy)
            flt_wcs.wcs.crpix += np.array([dx, dy])
            try:
                flt_wcs.sip.crpix += np.array([dx, dy])
            except:
                pass

            astrodrizzle.adrizzle.do_driz(sci, flt_wcs, wht, wcs_out, foutsci, foutwht, foutctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
            astrodrizzle.adrizzle.do_driz(contam_cutout, flt_wcs, wht, wcs_out, coutsci, coutwht, coutctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
            astrodrizzle.adrizzle.do_driz(cutout, flt_wcs, wht, wcs_out, moutsci, moutwht, moutctx, 1., 'cps', 1, wcslin_pscale=1.0, uniqid=1, pixfrac=pixfrac, kernel='square', fillval=0, stepsize=10, wcsmap=None)    
            ds9.frame(3)
            ds9.view(foutsci, header=hout)
            ds9.frame(4)
            ds9.view(foutsci-moutsci-coutsci, header=hout)
            
        drizzled_sci[root] = foutsci
        drizzled_wht[root] = foutwht
        drizzled_ctx[root] = foutctx
        drizzled_mod[root] = moutsci
        drizzled_contam[root] = coutsci
            
    sci = foutsci*0
    wht = sci*0
    csci = sci*0
    for root in drizzled_sci.keys():
        print root
        wht += drizzled_wht[root]
        sci += drizzled_sci[root]*drizzled_wht[root]
        csci += (drizzled_mod[root]+drizzled_contam[root])*drizzled_wht[root]
        
    sci /= wht
    csci /= wht
    ds9.frame(3)
    ds9.view(sci, header=hout)
    ds9.frame(4)
    ds9.view(sci-csci, header=hout)
    
    #wave = 10625.5
    
         
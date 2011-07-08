#!/usr/bin/env python
# encoding: utf-8
"""
galfit.py

Wrapper around GALFIT to fit the 3D-HST thumbnails.

$URL: https://subversion.assembla.com/svn/threedhst_internal/trunk/reduce.py $
$Author: gbrammer $
$Date: 2011-05-22 02:01:43 -0400 (Sun, 22 May 2011) $

"""
__version__ = " $Rev: 5 $"

import sys, os
import shutil
import glob
import numpy as np

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import pyfits

import threedhst.catIO as catIO

def fit_cosmos():
    import unicorn.galfit
    import threedhst
    
    sexCat = threedhst.sex.mySexCat('orient1_drz.cat')
    mag = np.cast[float](sexCat.MAG_F1392W)
    use = mag < 23.5
    
    root='orient1'
    
    #### matched catalog, do all objects that have NMBS matches to some K limit
    match = catIO.Readfile('../../../../COSMOS/HTML_v1.0/SED/'+root+'_match.cat')
    q = []
    for id in match.id_f140w:
        mat = np.where(sexCat.id == id)[0][0]
        q.append(mat)
    
    q = np.array(q)
    use = q[(match.mag_ktot < 23) & (match.rmatch < 1) & (match.logm >= 10.9)]
    
    for id in sexCat.id[use]:
        try:
            #
            status = unicorn.galfit.go_fit(id=id, fix_n=False, fit_sky=False,
                    tar_root=root,
                    PSF_IMAGE='star_'+root+'_thumb.fits')
            # status = unicorn.galfit.go_fit(id=id, fix_n=False,
            #       tar_root='orient2',
            #       PSF_IMAGE='star_orient2_01732_thumb.fits')
        except:
            pass
    
def fit_goodsn():
    import unicorn.galfit
    import threedhst
    
    tar_path = '/research/HST/GRISM/3DHST/GOODS-N/HTML_TAR/'
    
    catalogs = glob.glob('/research/HST/GRISM/3DHST/GOODS-N/HTML_v1.0/*drz.cat')
    psfs = glob.glob('star_ib*fits')
    matches = glob.glob('/research/HST/GRISM/3DHST/GOODS-N/HTML_v1.0/SED/*match.cat')
    for i,catalog in enumerate(catalogs):
        cat = threedhst.sex.mySexCat(catalog)
        mag = np.cast[float](cat.MAG_F1392W)
        
        match = catIO.Readfile(matches[i])
        
        q = []
        for id in match.id_f140w:
            mat = np.where(cat.id == id)[0][0]
            q.append(mat)
        
        q = np.array(q)
        use = q[(match.logm > 10.6) & (match.rmatch < 1)]
        
        root=os.path.basename(catalog).split('_drz')[0]
        for id in cat.id[use]:
            status = unicorn.galfit.go_fit(id=id, fix_n=False, tar_root=root,
                            tar_path = tar_path,
                            PSF_IMAGE=psfs[i])

def fit_uds_marshall():
    import unicorn.galfit
    import threedhst
    
    tar_path = '/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/images/'
    catalogs = glob.glob('/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/*drz.cat')
    PSF_IMAGE = 'star_UDS_direct1_00360.thumb.fits'
    matches = glob.glob('/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/SED/*match.cat')
    
    for i,catalog in enumerate(catalogs):
        cat = threedhst.sex.mySexCat(catalog)
        mag = np.cast[float](cat.MAG_F1249W)
        use = mag < 24.5
        root=os.path.basename(catalog).split('_drz')[0]
        
        match = catIO.Readfile(matches[i])
        q = []
        for id in match.id_f140w:
            mat = np.where(cat.id == id)[0][0]
            q.append(mat)
        
        q = np.array(q)
        use = q[(match.logm >= 10.9) & (match.rmatch < 1)]
        
        for id in cat.id[use]:
            if not os.path.exists('%s_%05d_galfit.png' %(root,id)):
                status = unicorn.galfit.go_fit(id=id, fix_n=False,
                    tar_root=root, fit_sky=False,
                    tar_path = tar_path, PSF_IMAGE=PSF_IMAGE)
    

def make_cat(symbol_color='red', masslim=(11,13), CIRCULARIZE=True):
    import unicorn.galfit
    import threedhst
    import numpy as np
    
    sexCat = threedhst.sex.mySexCat('orient1_drz.cat')
    sexCat.re = sexCat.id*0.-99
    sexCat.n = sexCat.id*0.-99
    sexCat.bovera = sexCat.id*0.-99
    sexCat.chi2 = sexCat.id*0.+100
    
    for i,id in enumerate(sexCat.id):
        #print id
        try:
            root='orient1_%05d_galfit' %(id)
            params = unicorn.galfit.read_log(root+'.log')
            x0, y0, mag, re, n, bovera, chi2 = params
            sexCat.re[i] = np.float(re)
            sexCat.n[i] = np.float(n)
            sexCat.bovera[i] = np.float(bovera)
            sexCat.chi2[i] = np.float(chi2)
            #print re
        except:
            pass
        #
        try:
            root='orient2_%05d_galfit' %(id)
            params = unicorn.galfit.read_log(root+'.log')
            x0, y0, mag, re, n, bovera, chi2 = params
            sexCat.re[i] = np.float(re)
            sexCat.n[i] = np.float(n)
            sexCat.bovera[i] = np.float(bovera)
            sexCat.chi2[i] = np.float(chi2)
            #print re
        except:
            pass
        
    ##### Make plot of sizes vs. z    
    match = catIO.Readfile('../../../../COSMOS/HTML_v1.0/SED/orient1_match.cat')

    match1 = catIO.Readfile('../../../../COSMOS/HTML_v1.0/SED/orient1_match.cat')
    q = []
    for id in match1.id_f140w:
        mat = np.where(sexCat.id == id)[0][0]
        q.append(mat)
    #
    match2 = catIO.Readfile('../../../../COSMOS/HTML_v1.0/SED/orient2_match.cat')
    for id in match2.id_f140w:
        mat = np.where(sexCat.id == id)[0][0]
        q.append(mat)
    
    for col in match.columns:
        match[col] = np.append(match1[col],match2[col])
        
    q = np.array(q)
    
    fp = open('orient1_galfit.cat','w')
    fp.write('# id_f140w re n bovera chi2 id_nmbs46 z_peak logm star_flag r_match  \n')
    for i,qi in enumerate(q):
        fp.write('%-5d %7.2f %7.2f %7.2f %9.1e  %-5d %8.3f %6.2f %1d  %5.2f\n' %(sexCat.id[qi], sexCat.re[qi], sexCat.n[qi], sexCat.bovera[qi], sexCat.chi2[qi], match.id_phot[i], match.z_peak[i], match.logm[i], match.star_flag[i], match.rmatch[i]))
    
    fp.close()
    
    zgrid = np.array([0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
    scale = np.array([3.268, 4.421, 5.343, 5.733,  6.082, 6.394, 6.673, 6.922, 7.144, 7.518, 7.812, 8.041, 8.216, 8.346, 8.439, 8.502, 8.539, 8.556, 8.555, 8.540, 8.512, 8.475, 8.430, 8.377, 8.320, 8.257, 8.192])
    
    use = (match.star_flag == 0) & (match.rmatch < 1)
    use = (match.logm > masslim[0]) & (match.logm < masslim[1])
    use = use & (sexCat.chi2[q] < 2)
    
    sint = np.interp(match.z_peak[use], zgrid, scale)
    
    xvd = np.array([0, 0.6, 1.1, 1.6, 2.0])
    yvd = np.array([12.4, 8.0, 5.3, 4.1, 3.0])
    
    if CIRCULARIZE:
        sexCat.re *= np.sqrt(1./np.abs(sexCat.bovera))
    
    plt.semilogy(match.z_peak[use], sexCat.re[q][use]*sint*0.06, marker='o', linestyle='None', color=symbol_color, markersize=12, alpha=0.7)
    print symbol_color
    plt.plot(xvd, yvd, marker='None', linewidth=10, alpha=0.3, color='black')
    
    #### Add GOODS-N
    # catalogs = glob.glob('/research/HST/GRISM/3DHST/GOODS-N/HTML_v1.0/*drz.cat')
    # tar_path = '/research/HST/GRISM/3DHST/GOODS-N/HTML_TAR/'
    # for catalog in catalogs:
    #     print catalog
    #     cat = threedhst.sex.mySexCat(catalog)
    #     mag = np.cast[float](cat.MAG_F1392W)
    #     use = mag < 23.5
    #     root=os.path.basename(catalog).split('_drz')[0]
    #     cat.re = cat.id*0.-99
    #     cat.n = cat.id*0.-99
    #     cat.bovera = cat.id*0.-99
    #     cat.chi2 = cat.id*0.+100
    #     for i,id in enumerate(cat.id):
    #         #print id
    #         try:
    #             rooti=root+'_%05d_galfit' %(id)
    #             params = unicorn.galfit.read_log(rooti+'.log')
    #             x0, y0, mag, re, n, bovera, chi2 = params
    #             cat.re[i] = np.float(re)
    #             cat.n[i] = np.float(n)
    #             cat.bovera[i] = np.float(bovera)
    #             cat.chi2[i] = np.float(chi2)
    #             #print re
    #         except:
    #             pass
    #     #
    #     match = catIO.Readfile('/research/HST/GRISM/3DHST/GOODS-N/HTML_v1.0/SED/'+root+'_match.cat')
    #     q = []
    #     for id in match.id_f140w:
    #         mat = np.where(cat.id == id)[0][0]
    #         q.append(mat)
    #     #
    #     q = np.array(q)
    #     #
    #     use = (match.rmatch < 1)
    #     use = (match.logm >= masslim[0]) & (match.logm <= masslim[1])
    #     use = use & (cat.chi2[q] < 2)
    #     sint = np.interp(match.z_peak[use], zgrid, scale)
    #     #
    #     if CIRCULARIZE:
    #         cat.re *= np.sqrt(1./np.abs(cat.bovera))
    #     
    #     print match.z_peak[use], cat.re[q][use]*sint*0.06
    #     
    #     plt.semilogy(match.z_peak[use], cat.re[q][use]*sint*0.06, marker='o', linestyle='None', color='green', markersize=10, alpha=0.6)
    
    # #### Add UDS
    catalogs = glob.glob('/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/*drz.cat')
    tar_path = '/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/images/'
    for catalog in catalogs[0:1]:
        cat = threedhst.sex.mySexCat(catalog)
        mag = np.cast[float](cat.MAG_F1249W)
        use = mag < 23.5
        root=os.path.basename(catalog).split('_drz')[0]
        cat.re = cat.id*0.-99
        cat.n = cat.id*0.-99
        cat.bovera = cat.id*0.-99
        cat.chi2 = cat.id*0.+100
        for i,id in enumerate(cat.id):
            #print id
            try:
                rooti=root+'_%05d_galfit' %(id)
                params = unicorn.galfit.read_log(rooti+'.log')
                x0, y0, mag, re, n, bovera, chi2 = params
                cat.re[i] = np.float(re)
                cat.n[i] = np.float(n)
                cat.bovera[i] = np.float(bovera)
                cat.chi2[i] = np.float(chi2)
                #print re
            except:
                pass
        #
        match = catIO.Readfile('/research/HST/GRISM/3DHST/SN-MARSHALL/HTML/SED/'+root+'_match.cat')
        q = []
        for id in match.id_f140w:
            mat = np.where(cat.id == id)[0][0]
            q.append(mat)
        #
        q = np.array(q)
        #
        use = (match.star_flag == 0) & (match.rmatch < 1)
        use = (match.logm >= masslim[0]) & (match.logm <= masslim[1])
        use = use & (cat.chi2[q] < 2)
        sint = np.interp(match.z_peak[use], zgrid, scale)
        
        print match.z_peak[use]
        if CIRCULARIZE:
            cat.re *= np.sqrt(1./np.abs(cat.bovera))        
        #
        plt.plot(match.z_peak[use], cat.re[q][use]*sint*0.06, marker='o', linestyle='None', color='orange', markersize=10, alpha=0.6)
       
        
    plt.ylim(0.5,30)
    #plt.ylim(-5,25)
    
    plt.xlim(0,3.5)
    plt.xlabel(r'$z_\mathrm{phot}$')
    plt.ylabel(r'$r_e$ [kpc]')
    plt.text(2.0, 13, r'F140W, $%4.1f\ <\ \log\ M/M_\odot\ <\ %4.1f$' %(masslim[0], masslim[1]))
    plt.savefig('3dhst_cosmos_sizes.png')

def z2_galaxies(zrange=(1.6,2.4)):
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/GALFIT/TEST/COSMOS/Z2_COSMOS')
    
    root_path = '/research/HST/GRISM/3DHST/COSMOS/HTML_v1.0/'
    root = 'orient1'
    
    co = catIO.Readfile('../orient1_galfit.cat')
    co.re *= np.sqrt(1./co.bovera)
    
    #### Select z~2
    use = (co.z_peak >= zrange[0]) & (co.z_peak <= zrange[1]) & (co.logm > 10.9) & (co.r_match < 1) & (co.re > 0)
    idx = np.arange(len(use))
    use = idx[use]
    use = use[np.argsort(co.re[use])]
    
    #### Angular scale
    zgrid = np.array([0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
    scale = np.array([3.268, 4.421, 5.343, 5.733,  6.082, 6.394, 6.673, 6.922, 7.144, 7.518, 7.812, 8.041, 8.216, 8.346, 8.439, 8.502, 8.539, 8.556, 8.555, 8.540, 8.512, 8.475, 8.430, 8.377, 8.320, 8.257, 8.192])
    sint = np.interp(co.z_peak, zgrid, scale)
    
    co.re_kpc = co.re*0.06*sint
    
    fp = open('z2.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="http://localhost/~gbrammer/COSMOS/scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[4,4]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> z </th>
        <th> logM </th>
        <th> r_e </th>
        <th> Thumb </th>
        <th> SED </th>
    </thead>
    <tbody>
    """)
    
    old = 0
    for i in use:
        if co.id_f140w[i] == old:
            root='orient2'
        else:
            root='orient1'
        old=co.id_f140w[i]
        
        file="%s_%05d" %(root, co.id_f140w[i])

        if not os.path.exists('%s/SED/%s_SED.png' %(root_path, file)):
            root='orient2'
            
        file="%s_%05d" %(root, co.id_f140w[i])
        
        fp.write("""
    <tr>
        <td> %s </td>
        <td> %5.2f </td>
        <td> %5.2f </td>
        <td> %5.2f </td>
        <td> <img src=../%s_galfit.png height=180px> </td>
        <td> <img src=%s/SED/%s_SED.png height=180px> </td>
    </tr>
    """ %(file, co.z_peak[i], co.logm[i], co.re_kpc[i],
          file,
          root_path, file))
    
    fp.write("</tbody></table></body></html>")
    fp.close()

def get_thumb(id=None, root='orient1',
    path='/research/HST/GRISM/3DHST/COSMOS/HTML_TAR/v1.0/', tarfile=True):
    
    thumb_file = root+'_%05d_thumb.fits' %(id)
    try:
        os.remove(thumb_file)
    except:
        pass
    
    if tarfile:
        status = os.system('tar xzvf '+path+'/'+root+'_thumbs.tar.gz %s.gz' %(thumb_file))
        if status != 0:
            return False
    else:
        shutil.copy(path+'/'+thumb_file+'.gz','.')
    
    return thumb_file
    
    #os.system('gunzip %s.gz' %(thumb_file))

def sync():
    import os
    os.system('rsync -avz /3DHST/Spectra/Work/ANALYSIS/GALFIT/*galfit.png /Users/gbrammer/Sites_GLOBAL/P/GRISM/GALFIT/')
    
def fit_3dhst_object(object='COSMOS-15-G141_00388', fit_sky=True, open=False, PSF_IMAGE='star_PSF.fits'):
    
    import unicorn.analysis
    import unicorn.galfit
    import os
    
    PATH = unicorn.analysis.get_grism_path(root=object.split('G141')[0]) 
    thumb_path = PATH+'/HTML/images/%s_thumb.fits.gz' %(object)
    print thumb_path
    
    if os.path.exists(thumb_path):
        #### First try, fit sersic + sky
        result = unicorn.galfit.go_fit(thumb_file=thumb_path, fix_n=False, fit_sky=fit_sky, PSF_IMAGE = PSF_IMAGE)
        if result:
            log = unicorn.galfit.GalfitLogfile(object+'_galfit.log')
            chi_s = log.chi2
            r_e_fit = log.list[0].re.value
            n_fit = log.list[0].n.value
        else:
            chi_s = 1.e5
            r_e_fit = 0.
            n_fit = 100
            
        # #### Fit is high-n, high-re, try adding a disk
        # if (r_e_fit*0.06 > 4) & (n_fit > 5):
        #     result = unicorn.galfit.go_fit(thumb_file=thumb_path, fix_n=False, add_disk=True)
        #     if result:
        #         log = unicorn.galfit.GalfitLogfile(object+'_galfit.log')
        #         chi_s = log.chi2
        #     else:
        #         chi_s = 1.e5
                
        #### If died or high chi2, fit PSF + sky
        if (result is False) | (chi_s > 50) | (r_e_fit == 0.01):
            result = unicorn.galfit.go_fit(thumb_file=thumb_path, psf_only=True, fit_sky=fit_sky, PSF_IMAGE=PSF_IMAGE)
            
            #### If died again, try turning off sky fit
            if not result:
                result = unicorn.galfit.go_fit(thumb_file=thumb_path, psf_only=True, fit_sky=False, PSF_IMAGE=PSF_IMAGE)
            
            if result:
                log = unicorn.galfit.GalfitLogfile(object+'_galfit.log')
                chi_p = log.chi2
            else:
                chi_p = chi_s*10
                
            if chi_p > chi_s:
                ### Revert to sersic fit
                result = unicorn.galfit.go_fit(thumb_file=thumb_path, fix_n=False, fit_sky=fit_sky)
            
def go_fit(thumb_file=None, mask=True, add_disk=False, fix_n=4,
           PSF_IMAGE='star_PSF.fits',
           fit_sky=True,
           psf_only=False,
           extension='_galfit'):
    
    import unicorn.galfit
    
    try:
        os.remove('galfit.01')
    except:
        pass

    root = os.path.basename(thumb_file).split('_thumb')[0]
    
    imgal = pyfits.open(thumb_file)
    #### Direct image exposure time (not right in full orient1 image)
    TEXP = 808
    #### Convert to ADU
    imgal[0].header.update('EXPTIME',TEXP)
    imgal[0].header.update('GAIN',1.0)
    imgal[0].data *= TEXP
    imgal.writeto('gal.fits', clobber=True)
    
    if mask:
        unicorn.galfit.make_segmap_mask(root)
    
    #### Star PSF
    star = pyfits.open(PSF_IMAGE)
    star[0].header.update('EXPTIME',TEXP)
    star[0].data *= TEXP
    star.writeto('psf.fits', clobber=True)
    
    unicorn.galfit.make_galfit_param(imgal, mask=mask, add_disk=add_disk, fix_n=fix_n, psf_only=psf_only, fit_sky=fit_sky)
    
    os.system('galfit galfit.feedme')
        
    try:
        shutil.move('fit.log',root+extension+'.log')
        shutil.move('imgblock.fits', root+extension+'.fits')
    except:
        return False
        
    os.system('ls galfit.[0-9]*')
    #logs = glob.glob('galfit.[0-9]*')
    #print logs
    #print logs[-1]
    shutil.move('galfit.01', root+extension+'.param')
    #os.system('mv galfit.01 '+root+'.param')
    unicorn.galfit.make_plot(root+extension)
    #os.system('ls galfit.[0-9]*')
    
    # if open:
    #     os.system('open '+root+extension+'.png')
        
    return True
    
def make_segmap_mask(root, extension='_galfit'):
    import threedhst
    
    se = threedhst.sex.SExtractor()
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    ## XXX add test for user-defined .conv file
    se.copyConvFile()

    se.overwrite = True
    se.options['CATALOG_NAME']    = root+extension+'.cat'
    se.options['CHECKIMAGE_NAME'] = root+extension+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'NONE'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '1.5'
    se.options['ANALYSIS_THRESH']  = '1.' 
    se.options['MAG_ZEROPOINT'] = '26.46'
    status = se.sextractImage('gal.fits')
    
    im = pyfits.open(root+extension+'_seg.fits')
    NY, NX = im[0].data.shape
    self_id = im[0].data[NY/2., NX/2]
    other = (im[0].data > 0) & (im[0].data != self_id)
    im[0].data *= 0
    im[0].data[other] = 1
    im.writeto('mask.fits', clobber=True)
    
    th = pyfits.open('gal.fits')
    seg = pyfits.open(root+extension+'_seg.fits')
    mask = seg[0].data == 0
    
    ### Subtract background
    th[0].data -= np.median(th[0].data[mask])
    th.writeto('gal.fits', clobber=True)
    
    ### Make rms image
    rms = np.std(th[0].data[mask])
    #th[0].data = th[0].data*0.+rms #*np.random.randn(NY, NX)
    GAIN = 2.5
    th[0].data = np.sqrt(np.abs(th[0].data)*GAIN/np.sqrt(4.)+rms**2)
    th.writeto('rms.fits', clobber=True)
    
def make_galfit_param(img, mask=False, fit_sky=True, add_disk=False, fix_n=4,
    psf_only=False):
   
    NY, NX = img[0].data.shape
        
    if mask:
        mask_file = 'mask.fits'
    else:
        mask_file = 'none'
    
    if fit_sky:
        fit_sky = 1
    else:
        fit_sky = 0
        
    if fix_n is False:
        nval = 4.
        fix_n = 1
    else:
        nval = fix_n
        fix_n = 0
    
    header = """ ================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) gal.fits            # Input data image (FITS file)
B) imgblock.fits       # Output data image block
C) none            # Sigma image name (made from data if blank or "none") 
D) psf.fits   #        # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) %s                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1    %d   1    %d   # Image region to fit (xmin xmax ymin ymax)
I) %d    %d          # Size of the convolution box (x y)
J) 26.46              # Magnitude photometric zeropoint 
K) 0.06  0.06        # Plate scale (dx dy)    [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
""" %(mask_file, NX, NY, int(NX*3./4), int(NY*3/4))
    
    #### Sky
    sky_model = """
# Sky
0) sky                    #  object type
1) 0.0         %d          #  sky background at center of fitting region [ADUs]
2) 0.0000      0          #  dsky/dx (sky gradient in x)
3) 0.0000      0          #  dsky/dy (sky gradient in y)
Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
""" %(fit_sky)
    
    #### Sersic, fix or free n
    sersic_model = """
# Object number: 1 (bulge)
 0) sersic                 #  object type
 1) %f %f  1 1  #  position x, y
 3) 21.31     1          #  Integrated magnitude	
 4) 4.0         1          #  R_e (half-light radius)   [pix]
 5) %f         %d          #  Sersic index n (de Vaucouleurs n=4) 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 1.0         1          #  axis ratio (b/a)  
10) 45.3690     1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
    """ %(NX/2.+0.5, NY/2.+0.5, nval, fix_n)
    
    #### Disk component, n=1 fixed
    disk_model = """
 0) sersic                 #  object type
 1) %f %f  1 1  #  position x, y
 3) 23.31     1          #  Integrated magnitude	
 4) 4.0         1          #  R_e (half-light radius)   [pix]
 5) 1.0         0          #  Sersic index n (de Vaucouleurs n=4) 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.2570      1          #  axis ratio (b/a)  
10) 45.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
    """ %(NX/2.+0.5, NY/2.+0.5)
    
    #### PSF component
    psf_model = """
 0) psf                 #  object type
 1) %f %f  1 1  #  position x, y
 3) 23.31     1          #  Integrated magnitude	
 4) 4.0         1          #  R_e (half-light radius)   [pix]
 5) 1.0         0          #  Sersic index n (de Vaucouleurs n=4) 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.2570      0          #  axis ratio (b/a)  
10) 45.3690    0          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
        """ %(NX/2.+0.5, NY/2.+0.5)
    
    #### Make the parameter file
    fp = open('galfit.feedme','w')
    fp.write(header)
    if psf_only:
        fp.write(psf_model)
    else:
        fp.write(sersic_model)
        if add_disk:
            fp.write(disk_model)
    
    fp.write(sky_model)
    fp.close()

def axis_ticks(im, ax):
    ## Axis ticks in units of arcsec
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    size = im[0].header.get('NAXIS1')
    pixel_scale = 0.06
    asec_pix = 1./pixel_scale
    nasec = int(size/asec_pix/2)
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size/2)
    ax.set_xticklabels([])
    ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size/2)
    
def make_plot(root):
    import matplotlib.pyplot as plt
    import unicorn.galfit
    
    plt.rcParams['image.cmap'] = 'gray'
    im = pyfits.open(root+'.fits')
    seg = pyfits.open(root+'_seg.fits')
    
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[6,2.3],dpi=100)
    else:
        fig = Figure(figsize=[6,2.3], dpi=100)
    
    fig.subplots_adjust(wspace=0.0,hspace=0.0,left=0.02,
                        bottom=0.02,right=0.98,top=0.87)
    
    #
    im_max = im[1].data.max()
    
    im[1].data = 0-im[1].data
    im[2].data = 0-im[2].data
    im[3].data = 0-im[3].data
    #im_max *= -1
    
    scl = 0.8
    
    ax = fig.add_subplot(131)    
    ax.imshow(im[1].data, vmin=-1*im_max*scl, vmax=0.1*im_max*scl, interpolation='nearest')
    axis_ticks(im, ax)

    ax = fig.add_subplot(132)    
    ax.imshow(im[2].data, vmin=-1*im_max*scl, vmax=0.1*im_max*scl, interpolation='nearest')
    axis_ticks(im, ax)
    
    scl = scl*0.2
    ax = fig.add_subplot(133)    
    ax.imshow(im[3].data, vmin=-1*im_max*scl, vmax=0.1*im_max*scl, interpolation='nearest')
    axis_ticks(im, ax)
    
    #### Segmentation mask
    # ax = fig.add_subplot(133)
    # #obj = seg[0].data > 0
    # #seg[0].data[obj] = -1
    # mask = 0-seg[0].data*1./seg[0].data.max()  
    # plt.imshow(mask, vmin=-1, vmax=0, interpolation='nearest')
        
    id = root.split('_')[1]
    id = root.replace('-G141','').replace('_galfit','')
    
    x0, y0, mag, re, n, bovera, chi2 = unicorn.galfit.read_log(root+'.log')
    re = '%5.2f' %(np.float(re.replace('*',''))*0.06)
    label = '#'+id+r'  $r_e$='+re+r'$^{\prime\prime}$  $n$='+n+r'  $b/a$='+bovera+r'  $\chi^2_\nu$='+chi2
    
    #### Smarter label
    log = unicorn.galfit.GalfitLogfile(root+'.log')
    label = id +r'  $\log\ \chi^2_\nu $=%.2f' %(np.log10(log.chi2))
    NCOMP = 0

    ax.text(-2,1.08, label, horizontalalignment='left',
      verticalalignment='center',transform = ax.transAxes)

    y0 = 1.08
    for comp in log.list:
        label = ''
        if comp.type == 'sersic':
            label = r'$r_e=%5.2f^{\prime\prime}$  $n$=%4.2f  $b/a$=%4.2f' %(comp.re.value*0.06, comp.n.value, comp.ba.value)
            NCOMP += 1
        if comp.type == 'psf':
            label = 'PSF '
            NCOMP += 1
        
        if NCOMP > 1:
            y0-=0.1
            
        if label:
            ax.text(1,y0, label, horizontalalignment='right',
                verticalalignment='center',transform = ax.transAxes)
        
    if USE_PLOT_GUI:
        fig.savefig(root+'.png',dpi=100,transparent=False)
        plt.close()
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(root+'.png', dpi=100, transparent=False)
    

def process_parameters(parline, errline):
    for char in ',:()':
        parline = parline.replace(char,'')   
        errline = errline.replace(char,'')   
    
    parline = parline.split()
    errline = errline.split()
     
    return parline, errline

class GalfitParameter():
    def __init__(self, value, error):
        self.fixed=False
        self.warning=False
        if isinstance(value,'string'.__class__):
            if value.find('*') >= 0:
                self.warning=True
                value = value.replace('*','')
                error = error.replace('*','')
        #print value
        
        try:
            self.value = float(value)
            self.error = float(error)
            #print self.value, self.error
        except:
            self.value = np.nan
            self.error = np.nan
            
class GalfitComponent():

    def __getitem__(self, param):
        if param not in self.params:
            print 'Parameter, %s, not found for component "%s"\n"' %(param, self.type)
            return np.nan
        
        stri = 'val = self.'+param+'.value'
        exec(stri)
        stri = 'err = self.'+param+'.error'
        exec(stri)
        return val, err

    def __init__(self, parline, errline):
        for char in ',:()[]':
            parline = parline.replace(char,'')   
            errline = errline.replace(char,'')   

        parline = parline.split()
        errline = errline.split()
        
        self.type = parline[0]
        if self.type == 'sersic':
            self.x0 = GalfitParameter(parline[1], errline[0])
            self.y0 = GalfitParameter(parline[2], errline[1])
            self.mag = GalfitParameter(parline[3], errline[2])
            self.re = GalfitParameter(parline[4], errline[3])
            self.n = GalfitParameter(parline[5], errline[4])
            self.ba = GalfitParameter(parline[6], errline[5])
            self.PA = GalfitParameter(parline[7], errline[6])
            self.params = ['x0','y0','mag','re','n','ba','PA']
        
        if self.type == 'psf':
            self.x0 = GalfitParameter(parline[1], errline[0])
            self.y0 = GalfitParameter(parline[2], errline[1])
            self.mag = GalfitParameter(parline[3], errline[2])
            self.params = ['x0','y0','mag']
        
        if self.type == 'sky':
            self.sky = GalfitParameter(parline[3], errline[0])
            self.dx = GalfitParameter(parline[4], errline[1])
            self.dy = GalfitParameter(parline[5], errline[2])
            self.params = ['sky','dx','dy']
        
        #### Info string
        info = 'Type: %s\n' %self.type
        for param in self.params:
            info+='%6s' %param + '=%g (%g)\n' %(self.__getitem__(param))
        self.info = info
        
class GalfitLogfile():
    def __init__(self, logfile):
        loglines = open(logfile).readlines()
        components = []
        for i,line in enumerate(loglines):
            if line.strip().startswith('sersic'):
                components.append(GalfitComponent(loglines[i], loglines[i+1]))
            
            if line.strip().startswith('psf'):
                components.append(GalfitComponent(loglines[i], loglines[i+1]))
                
            if line.strip().startswith('sky'):
                components.append(GalfitComponent(loglines[i], loglines[i+1]))
            
            if line.strip().startswith('Chi^2/nu'):
                chi2 = np.float(line.split()[2])
                self.chi2 = chi2
        
        self.list = components
        self.NCOMP = len(components)
        self.components = []
        for comp in components:
            self.components.append(comp.type)
            
    def __getitem__(self, idx):
        if isinstance(idx,'string'.__class__):
            for i,comp in enumerate(self.components):
                if idx == comp:
                    idx=i
                    break
        
        return self.list[idx]
        
def read_log(logfile):
    loglines = open(logfile).readlines()
    
    for i,line in enumerate(loglines):
        if line.strip().startswith('sersic'):
            parline, errline = process_parameters(line, loglines[i+1])
        #
        if line.strip().startswith('psf'):
            parline, errline = process_parameters(line, loglines[i+1])
        
        if line.strip().startswith('Chi^2/nu'):
            chi2 = np.float(line.split()[2])
            chi2 = '%4.2f' %(chi2)
        
    if parline[0] == 'sersic':
        x0 = parline[1]
        y0 = parline[2]
        mag = parline[3]
        re = parline[4]
        n = parline[5]
        bovera = parline[6]
    
    if parline[0] == 'psf':
        x0 = parline[1]
        y0 = parline[2]
        mag = parline[3]
        re = '-1'
        n = '-1'
        bovera = '-1'
        
    return x0, y0, mag, re, n, bovera, chi2
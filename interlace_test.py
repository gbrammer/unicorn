"""
Tests on the UDF interlaced data
"""

import unicorn
import threedhst
from threedhst.prep_flt_files import process_3dhst_pair as pair
import threedhst.prep_flt_files
import glob
import os

def prepare():
    os.chdir(unicorn.GRISM_HOME+'UDF/PREP_FLT')

    ALIGN = '/research/HST/GRISM/3DHST/GOODS-S/UDF/hlsp_hudf09_hst_wfc3ir_hudf09_F160W_v1_sci.fits'

    direct_files = glob.glob('i*30_asn.fits')
    grism_files = glob.glob('i*40_asn.fits')
    for direct, grism in zip(direct_files[1:], grism_files[1:]):
        pair(direct, grism, ALIGN_IMAGE = ALIGN, ALIGN_EXTENSION=0, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')

def reduce_interlace():
    
    import unicorn
    
    clean_all=True
    clean_spectra=False
    make_images=True
    make_model=True
    fix_wcs=True
    skip_completed_spectra=True
    MAG_LIMIT=26
    out_path='./'
    extract_limit=25
    
    files=glob.glob('GOODS-S-3*G141_asn.fits')
    for file in files[1:]:
        status = unicorn.reduce.reduce_pointing(file=file, clean_all=clean_all, clean_spectra=clean_spectra, make_images=make_images, make_model=make_model, fix_wcs=fix_wcs, extract_limit=extract_limit, skip_completed_spectra=skip_completed_spectra, MAG_LIMIT=MAG_LIMIT, out_path=out_path)
    
def compare_pieter_2D():
    
    import threedhst.dq
    import pyfits
    
    input = pyfits.open('GOODS-S-34-G141_inter.fits')
    me = pyfits.open('GOODS-S-34_model.fits')
    pvd1 = pyfits.open('Pieter/2D_out/GOODS-S/34/GOODS-S-34_G141_mod1.fits')
    pvd2 = pyfits.open('Pieter/2D_out/GOODS-S/34/GOODS-S-34_G141_mod2.fits')
    
    xi, yi = 400, 400
    dx, dy=10, 10
    
    NX, NY = 400, 400
    inp_sub = input[1].data[yi+dy:yi+dy+NY, xi+dx:xi+dx+NX]
    me_sub = me[0].data[yi+dy:yi+dy+NY, xi+dx:xi+dx+NX]
    pvd1_sub = pvd1[0].data[yi:yi+NY, xi:xi+NX]
    pvd2_sub = pvd2[0].data[yi:yi+NY, xi:xi+NX]
    
    ds9.frame(1)
    ds9.v(inp_sub, vmin=-0.05, vmax=0.2)
    ds9.frame(2)
    ds9.v(inp_sub-me_sub, vmin=-0.05, vmax=0.2)
    ds9.frame(3)
    ds9.v(inp_sub-pvd1_sub/4., vmin=-0.05, vmax=0.2)
    ds9.frame(4)
    ds9.v(inp_sub-pvd2_sub/4., vmin=-0.05, vmax=0.2)
    
def compare_pieter_1D():
    import threedhst
    import threedhst.catIO as catIO
    import unicorn
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    theirs = threedhst.sex.mySexCat('Pieter/1D/UDF-34-G141_drz.cat')
    idx = np.arange(theirs.id.shape[0])
    
    print 'Read SPC'
    SPC = threedhst.plotting.SPCFile('Pieter/1D/GOODS-S-34-G141_opt.SPC.fits', axe_drizzle_dir='./')
    print 'Done.'
    
    mine = catIO.Readfile('GOODS-S-34_inter.cat.wcsfix')
    
    mag = np.cast[float](theirs.MAG_F1392W)
    use = (mag > 17) & (mag < 23)
    
    for id in idx[use][np.argsort(mag[use])]:
        
        fig = unicorn.catalogs.plot_init(square=True, aspect=1./1.6, xs=5, left=0.08, right=0.05, bottom=0.08, top=0.08)
        
        spec = SPC.getSpec(theirs.id[id])
        if (spec is False) | (spec == -1):
            plt.close()
            continue
            
        dr = np.sqrt((mine.ra-theirs.ra[id])**2*np.cos(theirs.dec[id]/360.*2*np.pi)**2 + (mine.dec-theirs.dec[id])**2)*3600.
        match_id = mine.id[dr == dr.min()][0]
        my_spec = unicorn.reduce.Interlace1D('GOODS-S-34_%05d.1D.fits' %(match_id), PNG=False)
        
        plt.plot(spec['LAMBDA'], spec['FLUX']-spec['CONTAM'], color='blue')
        keep = (spec['LAMBDA'] > 1.2e4) & (spec['LAMBDA'] < 1.5e4) & np.isfinite(spec['FERROR'])
        yint = np.interp(spec['LAMBDA'], my_spec.data.wave, (my_spec.data.flux-my_spec.data.contam)/my_spec.data.sensitivity*1.e-17)
        anorm = np.sum(yint[keep]*(spec['FLUX']-spec['CONTAM'])[keep])/np.sum(yint[keep]**2)
        
        plt.plot(my_spec.data.wave, (my_spec.data.flux-my_spec.data.contam)/my_spec.data.sensitivity*1.e-17*anorm, color='red')
        
        plt.plot(spec['LAMBDA'], spec['FERROR'], color='blue', alpha=0.5)
        plt.plot(my_spec.data.wave, (my_spec.data.error)/my_spec.data.sensitivity*1.e-17*2.8, color='red', alpha=0.5)
        
        ymax = ((my_spec.data.flux-my_spec.data.contam)/my_spec.data.sensitivity*1.e-17*anorm)[(my_spec.data.wave > 1.1e4) & (my_spec.data.wave < 1.6e4)].max()
        plt.ylim(-0.05*ymax, 1.3*ymax)
        plt.xlim(1.05e4, 1.7e4)
        
        plt.title('GOODS-S-34_%05d - %.1f - %.1f' %(theirs.id[id], mag[id], anorm))
        
        plt.savefig('GOODS-S-34_%05d.compare.png' %(theirs.id[id]))
        plt.close()
        
def compare_multiple_spectra():
    import threedhst
    import threedhst.catIO as catIO
    import unicorn
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    cats = []    
    for i in [34,36,37,38]:
        cats.append(catIO.Readfile('GOODS-S-%d_inter.cat.wcsfix' %i))
        
    id34=466
    
    sp = unicorn.reduce.Interlace1D('GOODS-S-34_%05d.1D.fits' %(id34), PNG=False)
    wave = sp.data.wave*1.
    flux = sp.data.flux*1.
    err = sp.data.error*1.
    contam = sp.data.contam*1.
    sens = sp.data.sensitivity*1.
    pid = sp.data.wave*0.+34
    ix = np.arange(cats[0].id.shape[0])[cats[0].id == id34][0]
    
    pointings = [34,36,37,38]
    for i in range(1,4):
        dr = np.sqrt((cats[i].ra-cats[0].ra[ix])**2*np.cos(cats[i].ra/360.*2*np.pi)**2 + (cats[i].dec-cats[0].dec[ix])**2)*3600.
        match = dr == dr.min()
        if dr.min() > 1:
            continue
        #
        sp = unicorn.reduce.Interlace1D('GOODS-S-%d_%05d.1D.fits' %(pointings[i], cats[i].id[match][0]), PNG=False)
        wave = np.append(wave, sp.data.wave*1.)
        flux = np.append(flux, sp.data.flux*1.)
        err = np.append(err, sp.data.error*1.)
        contam = np.append(contam, sp.data.contam*1.)
        sens = np.append(sens, sp.data.sensitivity*1.)
        pid = np.append(pid, sp.data.contam*0+pointings[i])
        ix = np.arange(cats[0].id.shape[0])[cats[0].id == id34][0]
        
    #### Show results
    for p in pointings:
        mat = pid == p
        pl = plt.plot(wave[mat], (flux[mat]-contam[mat])/sens[mat], alpha=0.3)
        
    xwave = np.arange(1.1e4,1.705e4,22)
    xint, yint, yerr, narr = interlace_test.rebin(xwave, wave, (flux-contam)/sens, err/sens)
    #pl = plt.errorbar(xint, yint, yerr, marker='o', linestyle='None', alpha=0.5, ecolor='black', color='black')
    pl = plt.fill_between(xint, yint+yerr, yint-yerr, alpha=0.3, color='black')
    pl = plt.plot(xint, yint, alpha=0.5, color='black', linewidth=2)
    
    ymax = yint.max()
    plt.ylim(-0.05*ymax, 1.2*ymax)
    plt.xlim(1.05e4, 1.7e4)
    
    
    
def rebin(xint, x, y, err):
    import numpy as np
    xarr = xint*0.
    yint = xint*0.
    yerr = xint*0.
    narr = xint*0.
    var = err**2
    N = len(xint)
    for i in range(N-1):
        mat = (x >= xint[i]) & (x < xint[i+1])
        NMAT = mat.sum()
        if NMAT > 0:
            xint[i] = x[mat].sum()/NMAT
            yint[i] = y[mat].sum()/NMAT
            yerr[i] = np.sqrt(var[mat].sum())/NMAT
            narr[i] = NMAT
            
    return xint, yint, yerr, narr
    
    
    
    
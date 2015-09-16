"""
Tests on the UDF interlaced data
"""

import numpy as np
import matplotlib.pyplot as plt
import pyfits

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
    """
    Reduce an interlaced dataset
    """
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
    
def check_refit_redshifts():
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    
    os.system('echo "# id     mag     z_spec   z_old    z_new    z_source" > refit_redshifts.dat')    
    os.system('cat *new_zfit*dat |grep -v "id" >> refit_redshifts.dat')
     
    z = catIO.Readfile('refit_redshifts.dat')
    dz_old = (z.z_old-z.z_spec)/(1+z.z_spec)
    dz_new = (z.z_new-z.z_spec)/(1+z.z_spec)
    
    ok = (z.mag < 24) & (np.abs(dz_old) < 0.1)
    plt.scatter(z.z_spec[ok], dz_old[ok], marker='o', color='blue', alpha=0.5)
    xm, ym, ys, yn = threedhst.utils.runmed(z.z_spec[ok], dz_old[ok], NBIN=20, use_median=True, use_nmad=True)
    plt.plot(xm, ym, color='blue')
    plt.plot(xm, ys, color='red')
    plt.ylim(-0.1, 0.1); plt.xlim(0, 4)
    plt.savefig('refit_old.png'); plt.close()
    
    ok = (z.mag < 24) & (np.abs(dz_new) < 0.1)
    plt.scatter(z.z_spec[ok], dz_new[ok], marker='o', color='blue', alpha=0.5)
    xm, ym, ys, yn = threedhst.utils.runmed(z.z_spec[ok], dz_new[ok], NBIN=20, use_median=True, use_nmad=True)
    plt.plot(xm, ym, color='blue')
    plt.plot(xm, ys, color='red')
    plt.ylim(-0.1, 0.1); plt.xlim(0, 4)
    plt.savefig('refit_new.png'); plt.close()
    
    
###### New simultaneous redshift fitting
def refit_specz_objects():
    import unicorn.catalogs2 as cat2
    import unicorn.interlace_test as test
    
    zsp = cat2.SpeczCatalog()
    
    models = glob.glob('*model.pkl')#[::-1]
    for model_file in models:
        model = unicorn.reduce.process_GrismModel(model_file.split('_inter')[0])
    
        dr, idx = zsp.match_list(model.ra_wcs, model.dec_wcs)
        ok = dr < 0.3
        ids = model.objects[ok]
        mags = model.cat.mag[ok]
        zspec = zsp.zspec[idx][ok]

        for i in range(ok.sum()):
            id = model.objects[ok][i]
            mag = model.cat.mag[ok][i]
            zspec = zsp.zspec[idx][ok][i]
            zsource = zsp.source[idx][ok][i]
            #
            if os.path.exists('%s_%05d.new_zfit.png' %(model.root, id)):
                print '%s_%05d.new_zfit.png' %(model.root, id)
                continue
            #
            model.twod_spectrum(id, miny=12, maxy=20, USE_REFERENCE_THUMB=True)
            model.show_2d(savePNG=True)
            #
            self = test.SimultaneousFit('%s_%05d' %(model.root, id), RELEASE=False, p_flat=1.e-8, lowz_thresh=0.)
            if (self.status is False):
                continue
            #### Check for valid pixels
            ngood = np.sum(self.twod.im['WHT'].data != 0, axis=0)
            okpix = (ngood > 0) & (self.oned.lam > 1.1e4) & (self.oned.lam < 1.65e4)
            if okpix.sum() < 10:
                continue
            #
            print ' '
            fp = open('%s_%05d.new_zfit.dat' %(model.root, id),'w')
            header = 'id     mag     z_spec   z_old    z_new    z_source\n'
            fp.write(header)
            #
            self.zout.z_spec[self.ix] = zspec*1. 
            #### New fit, spectrum+photometry
            status = self.read_master_templates() # ; self.get_spectrum_tilt()
            if self.phot_use.sum() < 4:
                fp.close()
                continue
            #self.get_spectrum_tilt(fit_for_lines=True)
            #self.new_fit_in_steps(dzfirst=0.01, ignore_spectrum=False)
            self.new_fit_constrained(dzfirst=0.005, ignore_spectrum=False)
            z_new = self.z_max_spec*1.
            #### EMCEE fit
            self.emcee_fit(NWALKERS=20)
            #### Old fit, spectrum separate
            self.fit_in_steps()
            z_old = self.z_max_spec*1.
            #### Store results
            log = '%05d  %.3f  %.5f  %.5f  %.5f  %s\n' %(id, mag, zspec, z_old, z_new, zsource)
            print header+log
            fp.write(log)
            fp.close()

    #zz = zsp.zspec[idx][model.objects == id]

class MultiFit():
    def __init__(self, files=[], fcontam=2, min_err=0.05, z=None, znorm='z_peak', norm_maglimit=23.5, ignore_photometry=False, clip_edges=False, fitbg=False, smooth_oned=1, outroot=None, flatten_thumb=False):
        import threedhst.eazyPy as eazy
        
        if len(files) == 0:
            return False
        
        self.files = files
        self.fcontam = fcontam
        self.smooth_oned = smooth_oned
        
        #### Use GrismSpectrumFit as helper for reading EAZY products
        status = False
        for file in files:
            self.gris = unicorn.interlace_fit.GrismSpectrumFit(file.split('.2D.fits')[0].split('_2d.fits')[0], lowz_thresh=0)
            if self.gris.status:
                status = True
                break
                
        self.status = status
        if not status:
            threedhst.showMessage("Couldn't initialize GrismSpectrumFit", warn=True)
            return
                    
        self.lambdaz, self.temp_sed, self.lci, self.obs_sed, self.fobs, self.efobs = self.gris.eazy_fit
        self.efobs = np.sqrt(self.efobs**2+(min_err*self.fobs)**2)
        
        param = eazy.EazyParam(self.gris.zout.filename.replace('.zout', '.param'), READ_FILTERS=True)
        self.eazy_param = param

        self.phot_lines = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='full', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.gris.zout)
        # for i in range(self.phot_lines.NTEMP):
        #     plt.plot(self.phot_lines.templam, self.phot_lines.temp_seds[:,i]/np.interp(7000, self.phot_lines.templam, self.phot_lines.temp_seds[:,i]), label='%d' %(i), alpha=0.5)
            
        self.phot_mask = (self.fobs > -99)
        if ignore_photometry:
            self.efobs *= 100000
            #self.fobs *= 0
            
            self.gris.ix = 0
            self.gris.z_peak = -1
            
        self.Nphot = self.phot_mask.sum()
        
        self.Ndata = self.Nphot
        #### Read spectra
        self.twod = []
        for file in files:
            t = unicorn.reduce.Interlace2D(file, flatten_thumb=flatten_thumb)
            if not t:
                continue
            
            t.init_fast_model()
            
            t.bg, t.scl = 0., 1
            t.cwht = 1/np.sqrt(t.im['WHT'].data**2 + (self.fcontam*t.im['CONTAM'].data)**2)
            t.wht = 1/(t.im['WHT'].data**2 + (self.fcontam*t.im['CONTAM'].data)**2)
            wht_mask = (t.im['WHT'].data == 0) | ~np.isfinite(t.cwht) | (t.im['WHT'].data > 100)
            t.fast_compute_model()
            wht_mask |= (t.model < 0.05*t.model.max())
            t.wht[wht_mask] = 0.
            t.cwht[wht_mask] = 0.
            t.wave = np.cast[float](t.im['WAVE'].data)
            t.sens = np.cast[float](t.im['SENS'].data*(t.im[0].header['GROWX']*t.im[0].header['GROWY']))
                
            t.ok_mask = ~wht_mask
            
            ### don't use edges of the arrays where uncertain for large extended galaxies
            if clip_edges:
                #print 'Clip edges'
                #print t.ok_mask.sum()
                if t.grism_element == 'G102':
                    t.ok_mask &= (t.wave > 0.84e4) & (t.wave < 1.1e4)
                if t.grism_element == 'G141':
                    t.ok_mask &= (t.wave > 1.13e4) & (t.wave < 1.63e4)
                #
                if t.grism_element == 'G800L':
                    t.ok_mask &= (t.wave > 6000.) & (t.wave < 9500)
                
                t.wht[~t.ok_mask] = 0
                
                #print t.ok_mask.sum()
                
            t.NOK = t.ok_mask.sum()
            if t.NOK == 0:
                continue
                        
            self.twod.append(t)
            self.Ndata += t.NOK
            
        #
        if len(self.twod) == 0:
            self.status = False
            return
            
        if outroot is None:
            self.root = self.gris.root.split('-G1')[0]
        else:
            self.root = outroot
            
        self.id = self.twod[0].id
        
        # self.norm_to_template(); plt.ylim(-0.1, 1)
        
        #### Set up arrays for fits
        self.data_sci = np.zeros(self.Ndata)
        self.data_wht = np.zeros(self.Ndata)
        
        self.data_sci[:self.Nphot] = self.fobs[self.phot_mask]
        self.data_wht[:self.Nphot] = 1/self.efobs[self.phot_mask]**2
        
        ix = self.Nphot
        for t in self.twod:
            self.data_sci[ix:ix+t.NOK] = t.cleaned[t.ok_mask].flatten()
            self.data_wht[ix:ix+t.NOK] = t.wht[t.ok_mask].flatten()#*100
            ix += t.NOK
        
        ### Don't renorm if very faint
        self.mag = self.twod[-1].im[0].header['MAG']
        if self.mag > norm_maglimit:
            znorm = None
            
        self.znorm = znorm
        self.amin = 0.3
        self.amax = 5
        self.bmax = 0.08
        self.fitbg = fitbg
        
        if znorm is not None:
            if znorm == 'z_peak':
                znorm = self.gris.z_peak
            else:
                znorm = znorm

            self.fit_at_z(z=znorm, renorm=True, show=False)
        
        self.z_1 = None
        self.z_max_spec = self.gris.z_peak
            
        if z is not None:
            self.z = z
        
        self.status = True
    
    def mask_photometry(self, idx, z=None):
        self.efobs[idx] = 10000
        for t in self.twod:
            t.scl = 1.
            t.bg = 0
        
        self.update_data()
        if z is None:
            if self.gris.z_peak > 0:
                z = self.gris.z_peak*1
        
        if (z is not None) & (z > 0):        
            self.fit_at_z(z, renorm=True)
    
    def update_data(self):
        print 'update'
        
        self.data_sci[:self.Nphot] = self.fobs[self.phot_mask]
        self.data_wht[:self.Nphot] = 1/self.efobs[self.phot_mask]**2
        
        ix = self.Nphot
        for t in self.twod:
            self.data_sci[ix:ix+t.NOK] = t.cleaned[t.ok_mask].flatten()/t.scl - t.bg
            self.data_wht[ix:ix+t.NOK] = t.wht[t.ok_mask].flatten()/t.scl**2
            ix += t.NOK
    
    def ln_zgrid(self, zrange=[0,4], dz=0.001):
        zgrid = np.exp(np.arange(np.log(1.0+zrange[0]), np.log(1.+zrange[1]), dz))-1
        return zgrid
            
    def __call__(self, verbose=True, nmf_toler=1.e-5):
        self.chi2 = np.zeros(len(self.z))
        for i in range(len(self.z)):
            self.chi2[i] = self.fit_at_z(z=self.z[i], renorm=False, show=False, nmf_toler=nmf_toler)
            if verbose:
                print unicorn.noNewLine+'%.3f %10.2f' %(self.z[i], self.chi2[i])
                
        return self.chi2
                
    def norm_to_template(self, temp=None, eazy=True, show=True, simple=True):
        """
        Normalize using the EAZY template
        """
        
        flat = False
        if temp:
            xspec = temp[0]
            yspec = temp[1]
        elif eazy:
            xspec=np.cast[float](self.gris.eazy_fit[0])
            yspec=np.cast[float](self.gris.eazy_fit[1])
        else:
            xspec=np.cast[float](self.gris.eazy_fit[0])
            yspec=xspec*0.+1
            flat = True
        
        if show:
            fig = plt.figure(figsize=[8,4])
            ax = fig.add_subplot(111)
        
        for t in self.twod:  
            if flat:
                norm = 1
            else:
                norm = t.total_flux
                                          
            self.get_bg(t, xspec=xspec, yspec=yspec/norm, show=show, simple=simple)
            #print t.scl, t.bg
        
        self.update_data()
            
        if show:
            # ax2 = ax.twinx()
            # zgrid = np.arange(0.3,3.3)
            # for line in [3727, 5007, 6563.]:
            #     ax2.plot(line*(1+zgrid), zgrid)
            
            return fig
        
    def get_bg(self, t, xspec=None, yspec=None, show=False, simple=False, slope=False):
        """
        Fit for a small background component
        """
        import scipy.optimize
        
        width=6
        if yspec is not None:
            t.fast_compute_model(xspec, yspec)
        else:
            t.fast_compute_model(xspec, yspec)
            
        mask = (t.model > (0.02*t.model.max())) & t.ok_mask  #(t.wht > 0)
        
        if simple:
            #t.scl = 1./((t.cleaned*t.model*t.wht)[mask].sum()/(t.model**2*t.wht)[mask].sum())
            
            b = 0
            a0 = (t.cleaned*t.model*t.wht)[mask].sum()
            a1 = (t.model*t.wht)[mask].sum()
            a2 = (t.model**2*t.wht)[mask].sum()
            
            b0 = (t.cleaned*t.wht)[mask].sum()
            b1 = (t.model*t.wht)[mask].sum()
            b2 = t.wht[mask].sum()
            
            for i in range(10):
                a = (a0 - b*a1) / a2
                if self.fitbg:
                    b = (b0 - a*b1) / b2
                #print '%.4f %.4f' %(a, b)
                            
            #print 'x %.2f %.4f' %(a, b)
            #t.scl = t.cleaned[mask].sum()/t.model[mask].sum()
            if (a < self.amin) | (a > self.amax) | (np.abs(b) > self.bmax):
                a = 1.
                b = 0.
                            
            t.scl = a
            t.bg = b
            
            if slope:
                # y = (beta*(L-L0)+gamma)*model+b = a*model + b
                pass
                
        else:    
            wave, yspec_data, apcorr = t.trace_extract(t.cleaned-t.bg, width=width, dy=0, get_apcorr=True)
            #print xspec, yspec
            wave, yspec_model = t.trace_extract(t.model, width=width, dy=0)
            wave, yspec_wht = t.trace_extract(t.wht, width=width, dy=0)
            wht = 1/yspec_wht
            wht[yspec_wht == 0] = 0.
        
            components = np.array([yspec_model*t.scl, yspec_model*0.+1.])
            p0 = np.array([1, 0.])
            
            tol=1.49e-8  # not sure what this controls
            popt = scipy.optimize.leastsq(self.obj_lstsq, p0, args=(wave, yspec_data-t.bg*(2*width+1), components, wht), full_output=True, ftol=tol/10000., xtol=tol/10000.)
            xcoeff = popt[0]
            #xcoeff = [1,0]
            t.scl *= xcoeff[0]
            t.bg += xcoeff[1]/(2*width+1)
        
        if show:
            wave, yspec_data, apcorr = t.trace_extract(t.cleaned-t.bg, width=width, dy=0, get_apcorr=True)
            wave, yspec_model = t.trace_extract(t.model, width=width, dy=0)  
            
            ax = plt.gca()
            xok = t.im['SENS'].data > 0.15*t.im['SENS'].data.max()
            ax.plot(wave[xok], ((yspec_data-t.bg)/t.im['SENS'].data)[xok]/t.im[0].header['GROW']**2/t.scl*apcorr[xok], color='black', alpha=0.7)
            ax.plot(wave[xok], (yspec_model/t.im['SENS'].data)[xok]/t.im[0].header['GROW']**2*apcorr[xok], color='red', alpha=0.4, linewidth=2)
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.gris.eazy_fit
            ok = (fobs > -10) & (efobs > 0)
            ax.errorbar(lci[ok], fobs[ok], efobs[ok], marker='o', linestyle='None', alpha=0.5, color='blue', markersize=10)
            
    def fit_at_z(self, z=0, renorm=False, show=False, nmf_toler=1.e-5, norm_phot_only=True):    
        self.templates = np.zeros((self.Ndata, self.phot_lines.NTEMP))
        
        #### Photometry
        self.template_photometry = self.phot_lines.interpolate_photometry(z)
        self.templates[:self.Nphot, :] = self.template_photometry[self.phot_mask, :]
        
        self.z_fit = z
        
        ix = self.Nphot
        for t in self.twod:
            for i in range(self.phot_lines.NTEMP):
                t.fast_compute_model(self.phot_lines.templam*(1+z), self.phot_lines.temp_seds[:,i]/(1+z)**2/t.total_flux)
                self.templates[ix:ix+t.NOK, i] = t.model[t.ok_mask].flatten()
            
            ix += t.NOK
        
        init_coeffs = np.ones(10)
        if norm_phot_only:
            init_coeffs[7:9] = 0 #1.e-16
            ix_max = self.Nphot
        else:
            ix_max = -1
            
        ####
        var = 1/self.data_wht
        
        amatrix = unicorn.utils_c.prepare_nmf_amatrix(var[:ix_max], self.templates[:ix_max,:].T)
        self.coeffs = unicorn.utils_c.run_nmf(self.data_sci[:ix_max], var[:ix_max], self.templates[:ix_max,:].T, amatrix, toler=nmf_toler, init_coeffs = init_coeffs, verbose=False)
        
        ### Normalize spectra to best-fit 
        best_temp = np.dot(self.phot_lines.temp_seds, self.coeffs)
        self.best_temp = best_temp*1.
        if show:
            plt.plot(self.phot_lines.templam*(1+z), best_temp/(1+z)**2, color='purple', alpha=0.8)
        
        if renorm:
            self.norm_to_template(temp=[self.phot_lines.templam*(1+z), self.best_temp/(1+z)**2], eazy=False, show=False, simple=True)
            for t in self.twod:
                print 'iter1 %s %.2f %.4f' %(t.file, t.scl, t.bg)
                     
        self.coeffs[7:9] = 1.e-14
        
        #self.data_wht[:self.Nphot] *= 100
        for iter in range(1):
            var = 1/self.data_wht
            amatrix = unicorn.utils_c.prepare_nmf_amatrix(var, self.templates.T)
            self.coeffs = unicorn.utils_c.run_nmf(self.data_sci, var, self.templates.T, amatrix, toler=nmf_toler, MAXITER=1e6, init_coeffs = self.coeffs, verbose=False)
        
            best_temp = np.dot(self.phot_lines.temp_seds, self.coeffs)
            self.best_temp = best_temp*1.
        
            # if renorm:
            #     self.norm_to_template(temp=[self.phot_lines.templam*(1+z), self.best_temp/(1+z)**2], eazy=False, show=False, simple=True)
            #     for t in self.twod:
            #         print 'iter2 %s %.2f %.4f' %(t.file, t.scl, t.bg)
        
        # self.data_wht[:self.Nphot] /= 100.
        # amatrix = unicorn.utils_c.prepare_nmf_amatrix(var, self.templates.T)
        # self.coeffs = unicorn.utils_c.run_nmf(self.data_sci, var, self.templates.T, amatrix, toler=nmf_toler, MAXITER=1e6, init_coeffs = self.coeffs, verbose=False)
        #     
        # best_temp = np.dot(self.phot_lines.temp_seds, self.coeffs)
        # self.best_temp = best_temp*1.
        
        ### Check
        if False:
            for t in self.twod:
                m = t.model*0.
                for i in range(self.phot_lines.NTEMP):
                    t.fast_compute_model(self.phot_lines.templam*(1+z), self.phot_lines.temp_seds[:,i]/(1+z)**2/t.total_flux*self.coeffs[i])
                    m += t.model
            
                t.fast_compute_model(self.phot_lines.templam*(1+z), self.best_temp/(1+z)**2/t.total_flux)
                resid = (t.cleaned - m*t.scl - t.bg)#*np.sqrt(t.wht)
                    
        if show:
            best_phot = np.dot(self.template_photometry, self.coeffs)
            best_temp = np.dot(self.phot_lines.temp_seds, self.coeffs)
            plt.plot(self.phot_lines.templam*(1+z), best_temp/(1+z)**2, color='blue', alpha=0.8)
            plt.scatter(self.phot_lines.lc, best_phot, color='red')
            plt.errorbar(self.lci[self.phot_mask], self.fobs[self.phot_mask], self.efobs[self.phot_mask], marker='o', linestyle='None', color='black')
            for t in self.twod:
                wave, yspec_data = t.optimal_extract((t.cleaned-t.bg)/t.scl)
                xok = t.im['SENS'].data > 0.15*t.im['SENS'].data.max()
                plt.plot(wave[xok], ((yspec_data)/t.im['SENS'].data)[xok]/t.im[0].header['GROW']**2, color='black', alpha=0.7)
                
                t.fast_compute_model(self.phot_lines.templam*(1+z), best_temp/(1+z)**2/t.total_flux)
                wave, yspec_model = t.optimal_extract(t.model) 
                plt.plot(wave[xok], (yspec_model/t.im['SENS'].data)[xok]/t.im[0].header['GROW']**2, color='orange', linewidth=2, alpha=0.7)
    
            plt.xlim(3000,3.e4)
        
        full_model = np.dot(self.coeffs, self.templates.T)
        chi2 = np.sum((self.data_sci-full_model)**2*self.data_wht)
        return chi2
    
    def fit_in_steps(self, zrfirst=[0,3.8], dzfirst=0.01, dzsecond=0.0003, verbose=True, minwidth=0.02, save_dat=True, make_plot=True):
        """
        Fit grism redshifts
        """
        self.z_1 = self.ln_zgrid(zrange=zrfirst, dz=dzfirst)
        self.chi_1 = self.z_1*0+1e10
        
        print '\n'
        for i in range(len(self.z_1)):
            self.chi_1[i] = self.fit_at_z(z=self.z_1[i], renorm=False, show=False, nmf_toler=1.e-4)
            if verbose:
                print unicorn.noNewLine+'%.3f %10.2f (%.4f)' %(self.z_1[i], self.chi_1[i], self.z_1[np.argmin(self.chi_1)])
        
        self.chi2scale_1 = self.chi_1.min()/self.Ndata
        self.lnp_1 = np.exp(-0.5/self.chi2scale_1*(self.chi_1-self.chi_1.min()))
        self.lnp_1 /= np.trapz(self.lnp_1, self.z_1)
        
        ### Get interval for zoom-in
        pzcum = np.cumsum(self.lnp_1*dzfirst*(1+self.z_1))
        self.ci_1 = np.interp([0.01, 0.5, 0.99], pzcum/pzcum[-1], self.z_1)
        self.lim_1 = np.array([-1,1])*np.maximum(minwidth,  np.diff(self.ci_1)/(1+self.ci_1[1]))*(1+self.ci_1[1]) + self.ci_1[1]
        self.lim_1 = np.clip(self.lim_1, 0, 13)
        
        print '\nFirst %.4f (%.4f, %.4f)\n\n' %(self.ci_1[1], self.lim_1[0], self.lim_1[1])
        
        self.z_2 = self.ln_zgrid(zrange=self.lim_1, dz=dzsecond)
        self.chi_2 = self.z_2*0.+1e10
        
        for i in range(len(self.z_2)):
            self.chi_2[i] = self.fit_at_z(z=self.z_2[i], renorm=False, show=False, nmf_toler=1.e-4)
            if verbose:
                print unicorn.noNewLine+'%.3f %10.2f (%.4f)' %(self.z_2[i], self.chi_2[i], self.z_2[np.argmin(self.chi_2)])
        
        self.chi2scale_2 = self.chi_2.min()/self.Ndata
        self.lnp_2 = np.exp(-0.5/self.chi2scale_2*(self.chi_2-self.chi_2.min()))
        self.lnp_2 /= np.trapz(self.lnp_2, self.z_2)
        
        pzcum = np.cumsum(self.lnp_2*dzfirst*(1+self.z_2))
        self.ci_2 = np.interp([0.01, 0.025, 0.158, 0.5, 0.842, 0.975, 0.99], pzcum/pzcum[-1], self.z_2)
        print '\nSecond %.4f (%.4f, %.4f)\n' %(self.ci_2[3], self.ci_2[1], self.ci_2[-2])
        
        self.z_max_spec = self.ci_2[3]
        
        if save_dat:
            fp = open('%s_%05d.multifit.dat' %(self.root, self.id),'w')
            names = []
            scls = []
            bgs = []
            for t in self.twod:
                names.append(t.file)
                scls.append('%16.3f' %(t.scl))
                bgs.append('%16.5f' %(t.bg))

            fp.write(' '.join(names)+'\n')
            fp.write('scale %05d' %(self.id) + ' '.join(scls)+'\n')
            fp.write('bg    %05d' %(self.id) + ' '.join(bgs)+'\n')
            ### 2.5, 50, 97.5 percentiles
            fp.write('%s %d z_gris %.4f %.4f %.4f' %(self.root, self.id, self.ci_2[1], self.ci_2[3], self.ci_2[5]))        
            fp.close()
            
        if make_plot:
            self.diagnostic_figure()
            
    def diagnostic_figure(self, xlim=[0.4,2.3]):
        """
        Make figure with results of MultiFit
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import scipy.ndimage as nd
        
        gs1 = gridspec.GridSpec(3,4)
        gs1.update(left=0.01, right=0.99, bottom=0.75, top=0.99, wspace=0.01, hspace=0.01)
        
        chi_best = self.fit_at_z(z=self.z_max_spec, renorm=(self.znorm is not None), show=False, nmf_toler=1.e-4)

        best_phot = np.dot(self.template_photometry, self.coeffs)
        best_temp = np.dot(self.phot_lines.temp_seds, self.coeffs)
        c_nolines = self.coeffs*1.
        c_nolines[7:9] = 0
        best_temp_nolines = np.dot(self.phot_lines.temp_seds, c_nolines)
        
        plt.rcParams['font.size'] = 10
        fig = plt.figure(figsize=[10,4])
        plt.set_cmap('gray_r')
        lmin, lmax = 1e10, 0
        for i, t in enumerate(self.twod):
            lmin = np.minimum(lmin, t.wave.min())
            lmax = np.maximum(lmax, t.wave.max())
            scl = np.clip(t.model.max(), 0.03, 1)
            scl = np.clip(np.percentile((t.cleaned/t.scl-t.bg), 95)*1.1, 0.03, 1)
            #print scl
                            
            ax = fig.add_subplot(gs1[i*2])
            ax.imshow((t.cleaned/t.scl-t.bg)/scl, origin='lower', vmin=-0.1, vmax=1, interpolation='nearest', aspect='auto')
            self.axis_ticks(ax, t)

            ### Show continuum residuals
            ax = fig.add_subplot(gs1[i*2+1])
            t.fast_compute_model(self.phot_lines.templam*(1+self.z_fit), best_temp_nolines/(1+self.z_fit)**2/t.total_flux)            
            show = (t.cleaned/t.scl-t.bg - t.model)/scl*(t.wht != 0)
            show = nd.gaussian_filter(show, 1.5)*3
            ax.imshow(show, origin='lower', vmin=-0.1, vmax=1, interpolation='nearest', aspect='auto')
            ax.text(0.99, 0.95, t.file.split('.2D.fits')[0], ha='right', va='top', transform=ax.transAxes, fontsize=6, color='red', backgroundcolor='white')
            
            t.fast_compute_model(self.phot_lines.templam*(1+self.z_fit), best_temp/(1+self.z_fit)**2/t.total_flux)
            
            #ax.imshow(t.model/scl*(t.wht != 0), origin='lower', vmin=-0.02, vmax=1, interpolation='Nearest', aspect='auto')
            self.axis_ticks(ax, t)
        
        #
        gs2 = gridspec.GridSpec(1,1)
        gs2.update(left=0.06, right=0.5, bottom=0.08, top=0.74, wspace=0.1, hspace=0.01)

        ### SED + spec
        ax = fig.add_subplot(gs2[0])
        
        ax.plot(self.phot_lines.templam*(1+self.z_fit)/1.e4, best_temp/(1+self.z_fit)**2, color='black', alpha=0.2)
        ax.scatter(self.lci/1.e4, best_phot, color='red', alpha=0.5)
        ax.errorbar(self.lci[self.phot_mask]/1.e4, self.fobs[self.phot_mask], self.efobs[self.phot_mask], marker='o', linestyle='None', color='black', alpha=0.5)
        
        ymax = 0
        ymin = 1e10
        
        for t in self.twod:
            wave, yspec_data = t.optimal_extract((t.cleaned-t.bg)/t.scl)
            xok = t.im['SENS'].data > 0.15*t.im['SENS'].data.max()
            t.wopt, t.fopt = wave[xok]/1.e4, ((yspec_data)/t.im['SENS'].data)[xok]/t.im[0].header['GROWX']/t.im[0].header['GROWY']
            if self.smooth_oned > 0:
                ax.plot(t.wopt, nd.gaussian_filter(t.fopt, self.smooth_oned), alpha=0.5, linestyle='steps-mid')
            else:
                ax.plot(t.wopt, t.fopt, alpha=0.5, linestyle='steps-mid')
                
            #unicorn.plotting.fill_between_steps(t.wopt, nd.gaussian_filter(t.fopt, 1.4), t.fopt*0, alpha=0.1, color='black')
            
            t.fast_compute_model(self.phot_lines.templam*(1+self.z_fit), best_temp/(1+self.z_fit)**2/t.total_flux)
            wave, yspec_model = t.optimal_extract(t.model) 
            t.mopt = (yspec_model/t.im['SENS'].data)[xok]/t.im[0].header['GROWX']/t.im[0].header['GROWY']
            ax.plot(wave[xok]/1.e4, t.mopt, color='black', linewidth=2, alpha=0.5)
            ymax = np.maximum(ymax, t.mopt.max())
            ymin = np.minimum(ymin, t.mopt.min())
            
        ax.set_xlim(lmin/1.e4, lmax/1.e4)
        ax.set_ylim(0.6*ymin, ymax*1.2)
        ax.yaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune='both'))
        
        ### Overview
        gs3 = gridspec.GridSpec(1,2, width_ratios=[1,1])
        gs3.update(left=0.52, right=0.99, bottom=0.07, top=0.45, wspace=0.05, hspace=0.01)
        
        ax = fig.add_subplot(gs3[0])
        
        ax.plot(self.phot_lines.templam*(1+self.z_fit)/1.e4, best_temp/(1+self.z_fit)**2, color='black', alpha=0.6, zorder=4)
        ax.scatter(self.lci/1.e4, best_phot, color='red', alpha=0.5)
        ax.errorbar(self.lci[self.phot_mask]/1.e4, self.fobs[self.phot_mask], self.efobs[self.phot_mask], marker='o', linestyle='None', color='black', alpha=0.5)
        
        ax.semilogx()
        ax.set_yticklabels([])
        for t in self.twod:
            ax.plot(t.wopt, nd.gaussian_filter(t.fopt, 3), alpha=0.5, linestyle='steps-mid', zorder=1)
        
        ax.set_xlim(xlim)
        ymax = self.fobs[self.phot_mask].max()
        ax.set_ylim(-0.1*ymax, 1.5*ymax)
        
        #### Label
        ax.text(0.1, 1.6, r'%s %05d  $m_\mathrm{%s}=%0.2f$ (ID = %6s, dr = %6.2f")' %(self.root, self.id, self.twod[-1].im[0].header['FILTER'], self.mag, str(self.gris.zout.id[self.gris.ix]), self.gris.dr), transform=ax.transAxes, ha='left', va='top')
        ax.text(0.1, 1.4, 'z_spec = %6.3f   z_peak = %6.3f   z_grism = %6.3f' %(self.gris.zout.z_spec[self.gris.ix], self.gris.z_peak, self.z_max_spec), transform=ax.transAxes, ha='left', va='top')
        zsp = self.gris.zout.z_spec[self.gris.ix]
        if zsp > 0:
            ax.text(0.1, 1.2, r'$\Delta z$=%.4f' %((self.z_max_spec-zsp)/(1+zsp)), transform=ax.transAxes, ha='left', va='top')

        ### p(z)
        ax = fig.add_subplot(gs3[1])
        phot_pz = np.exp(self.gris.phot_lnprob)
        ax.plot(self.gris.phot_zgrid, phot_pz/phot_pz.max(), color='green', alpha=0.5, linewidth=2)
        
        ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune='both'))

        if self.z_1 is not None:
            ax.plot(self.z_1, self.lnp_1/self.lnp_1.max(), color='blue', alpha=0.2, linewidth=1)
            ax.plot(self.z_2, self.lnp_2/self.lnp_2.max(), color='blue', alpha=0.6, linewidth=2)
            ax.set_xlim(self.z_2.min(), self.z_2.max())
        
        ax.semilogy()
        ax.set_ylim(1.e-3, 1)
        ax.set_yticklabels([])
        
        fig.savefig('%s_%05d.multifit.png' %(self.root, self.id))
        plt.close()
        plt.rcParams['font.size'] = 12
        
    def axis_ticks(self, ax, t):
        if t.grism_element == 'G141':
            xticks = np.arange(1.1,1.75,0.1)   
        else:
            xticks = np.arange(0.8,1.25,0.1)
                
        sh = t.cleaned.shape
        tick_xarr = np.arange(sh[1])
        tick_int = np.interp(xticks*1.e4, t.wave, tick_xarr)
        ax.set_xticks(tick_int); ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        wlim = unicorn.reduce.grism_wlimit[t.grism_element]
        xlim_int = np.interp(np.array(wlim[:2]), t.wave, tick_xarr)
        ax.set_xlim(xlim_int)
        ax.set_ylim([0, sh[0]-1])
        
class SimultaneousFit(unicorn.interlace_fit.GrismSpectrumFit):
    """
    Extension of the `unicorn.interlace_fit.GrismSpectrumFit` class. 
    
    Now allows simultaneous fits of spectra + photometry
    
    """
    def __init__(self, *args, **kwargs):
        unicorn.interlace_fit.GrismSpectrumFit.__init__(self, *args, **kwargs)
        self.read_master_templates()
        self.new_load_fits()
        
    def read_master_templates(self):
        """
        Read the eazy-interpolated continuum + emission line templates
        """
        import scipy.ndimage as nd
        import threedhst.eazyPy as eazy
        
        #### Read template error function
        self.te_phot = eazy.TemplateError(unicorn.GRISM_HOME+'templates/TEMPLATE_ERROR.eazy_v1.0')
        self.te_spec = eazy.TemplateError(unicorn.GRISM_HOME+'templates/TEMPLATE_ERROR.v1.0.emlines')
        
        #### Read template photometry
        #self.phot = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='eazy_v1.1', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.zout)
        self.phot_lines = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='full', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.zout)
        self.phot_eazy = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='eazy_v1.1', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.zout)
        self.phot = self.phot_lines
        
        # if True:
        #     threedhst.showMessage('Conroy model', warn=True)
        #     cvd12 = np.loadtxt(unicorn.GRISM_HOME+'/templates/cvd12_t11_solar_Chabrier.dat')
        #     self.phot.temp_seds[:,0] = np.interp(self.phot.templam, cvd12[:,0], cvd12[:,1])
            
        #### SPECTRUM
        self.poly_tilt_coeffs = [1]
        self.beta_tilt_coeffs = [0]
        self.tilt_function = self.beta_tilt
        
        tilt_file = '%s.new_zfit_tilt.dat' %(self.root)
        if os.path.exists(tilt_file):
            lines = open(tilt_file).readlines()
            i=1
            self.poly_tilt_coeffs = []
            while 'beta' not in lines[i]:
                self.poly_tilt_coeffs.append(float(lines[i]))
                i += 1
            #
            self.beta_tilt_coeffs = []
            i += 1
            for line in lines[i:]:
                self.beta_tilt_coeffs.append(float(lines[i]))
                i += 1
                
        self.shape2D = self.twod.im['SCI'].data.shape
        
        ## Flatten the 2D spectroscopic flux and variance arrays
        ## Variance
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        #var += ((0.02*self.twod.im['SCI'].data)**2).flatten()
        #var += ((0.05*self.twod.im['SCI'].data)**2).flatten()
        
        #### Smooth contam model along dispersion axis
        sm_kern = np.ones((3,30*2./self.twod.growx))
        sm_err = nd.convolve(self.twod.im['CONTAM'].data, sm_kern/sm_kern.sum())
        var += (self.contam_ferror*sm_err**2).flatten()
        
        var[var == 0] = 1.e6
        #var += ((0.5*self.twod.im['CONTAM'].data)**2).flatten()
        var += (self.contam_ferror*(self.twod.im['CONTAM'].data)**2).flatten()
        
        #### Include sensitivity curve uncertainty
        # status = self.compute_model_function(np.cast[float](unicorn.reduce.sens_files['A']['WAVELENGTH']), np.cast[float](unicorn.reduce.sens_files['A']['WAVELENGTH'])*0 + 1.)
        # sens_err_flat = self.twod.model*1
        # 
        # status = self.compute_model_function(np.cast[float](unicorn.reduce.sens_files['A']['WAVELENGTH']), np.cast[float](unicorn.reduce.sens_files['A']['ERROR']/unicorn.reduce.sens_files['A']['SENSITIVITY'])**2)
        # sens_err_twod = np.sqrt(self.twod.model / sens_err_flat)*2
        # sens_err_twod[sens_err_flat == 0] = 0
        # 
        # #### Flat model
        # status = self.compute_model_function(np.cast[float](unicorn.reduce.sens_files['A']['WAVELENGTH']), np.cast[float](1./unicorn.reduce.sens_files['A']['SENSITIVITY']/1.e-17/6.))
        # var += ((self.twod.model*sens_err_twod)**2).flatten()
        
        ### inverse sensitivity, minimum error ~5%
        MIN_ERR = 0.05
        inv = 1./unicorn.reduce.sens_files['A']['SENSITIVITY']
        inv *= 1./np.min(inv)
        status = self.compute_model_function(np.cast[float](unicorn.reduce.sens_files['A']['WAVELENGTH']), np.cast[float](inv)**2*MIN_ERR) #/self.twod.total_flux*0.05)
        var += (self.twod.model**2).flatten()
        
        
        #### Downweight edges of grism sensitivity
        x = self.twod.im['WAVE'].data
        dx = 200
        #### G141
        xlimit = [1.1e4, 1.64e4]
        #xlimit = [1.2e4, 1.6e4]
        if self.grism_element == 'G102':
            xlimit = [0.8e4, 1.1e4]
        
        sens_weight = x*0.+0.
        sens_weight += np.minimum(np.exp(-(xlimit[0]-x)**2/2./dx**2), 1); sens_weight[x < xlimit[0]] = 1.
        sens_weight += np.minimum(np.exp(-(xlimit[1]-x)**2/2./dx**2), 1); sens_weight[x > xlimit[1]] = 1.
        var *= (np.ones((self.shape2D[0],1))*(1+sens_weight*2.)).flatten()
        
        self.spec_var = var
        
        #### Bad pixel outliers
        # import scipy.ndimage as nd
        # resid = (self.twod.im['SCI'].data - self.twod.im['CONTAM'].data)#/self.twod.im['WHT'].data
        # resid[self.twod.im['WHT'].data == 0] = 0
        # max = nd.maximum_filter(resid, size=3)#-resid
        # med = nd.median_filter(resid, size=3)#-resid
        
        # 
        # avg[avg == 0] = 0
        # flag = ((resid-avg) > 20).flatten()
        # self.spec_flux[flag] = 0
        #self.twod.im['WHT'].data[flag] = 0
        #self.twod.im['SCI'].data[flag] = 0
        #print 'Flagged %d outlier pixels.' %(flag.sum())
        
        ## Flux
        self.spec_flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        self.spec_use = np.isfinite(self.spec_flux) & (self.twod.im['WHT'].data.flatten() != 0)
        
        ### Use only pixels with flux in the model
        self.compute_model_function()
        xmax = np.max(self.twod.model, axis=0)
        self.spec_use &= (self.twod.model/xmax > 0.1).flatten()
        
        #### PHOTOMETRY
        ## best-fit EAZY SED
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
        self.lc = lci*1
        self.phot_flam = fobs*1
        self.phot_eflam = efobs*1
        #self.phot_fnu = fobs*(lci/5500.)**2
        #self.phot_efnu = efobs*(lci/5500.)**2
        self.phot_use = (fobs > -99) & (efobs > 0)

        #### Masks for fitting the photometry and/or spectra
        self.fitting_mask = {}
        self.fitting_mask['phot_only'] = self.phot_use
        
        self.fitting_mask['spec_only'] = np.zeros(self.phot.NFILT + np.size(self.spec_flux), dtype=np.bool)
        self.fitting_mask['spec_only'][self.phot.NFILT:] |= self.spec_use
        
        self.fitting_mask['both'] = np.zeros(self.phot.NFILT + np.size(self.spec_flux), dtype=np.bool) | self.fitting_mask['spec_only']
        self.fitting_mask['both'][:self.phot.NFILT] |= self.phot_use
         
        self.initialize_prior()          
        #
        #self.coeffs = np.ones(self.phot.NTEMP)
        self.coeffs = None
        self.zgrid = None
        self.lnprob_spec = None
         
    def initialize_prior(self):
        """
        Get EAZY apparent mag prior
        """
        from scipy import interpolate
        
        prior = open('%s/data/prior_K_zmax7_coeff.dat' %(os.path.dirname(threedhst.__file__))).readlines()
        prior_mag = np.cast[float](prior[0].split()[2:])
        prior_z0 = np.cast[float](prior[1].split()[1:])
        prior_gamma = np.cast[float](prior[2].split()[1:])
        
        z0 = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_z0, left=prior_z0[0], right=prior_z0[-1])
        gamma = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_gamma, left=prior_gamma[0], right=prior_gamma[-1])
        
        zgrid = self.ln_zgrid(zrange=[0.0001,10], dz=0.01)
        #threedhst.showMessage('z0=%f; gamma=%f; mag=%f' %(z0, gamma, self.twod.im[0].header['MAG']), warn=True)
        prior_pz = np.log(zgrid**gamma * np.exp(-(zgrid / z0)**gamma))
        #prior_pz[0] = prior_pz[1]-10
        ok = np.isfinite(prior_pz)
        prior_pz[~ok] = np.min(prior_pz[ok])
        prior_pz -= np.log(np.trapz(np.exp(prior_pz), zgrid)) #prior_pz.max()
        
        self.get_prior = interpolate.InterpolatedUnivariateSpline(zgrid, prior_pz, k=1)
        #self.get_prior = interpolate.LinearNDInterpolator(zgrid, prior_pz)
        
        #return prior_pz
            
    def poly_tilt(self, x):
        """
        Polynomial tilt scaling
        """
        from scipy import polyval
        return polyval(self.poly_tilt_coeffs, x)
    
    def beta_tilt(self, x):
        """
        c = self.beta_tilt_coeffs
        beta, alpha = c[0], np.exp(c[1])
        Scaling = alpha * (x/1.e4)**beta
        """
        from scipy import polyval
        return np.exp(polyval(self.beta_tilt_coeffs, np.log(x/1.4e4)))
        
    def get_spectrum_tilt(self, z=None, make_figure=True, order=6, faint_limit=24, fit_for_lines=True, NO_GUI=True):
        from scipy import polyfit, polyval
        
        if z is None:
            z = self.z_peak
                
        if fit_for_lines:
            self.new_fit_zgrid(dz=0.01, ignore_spectrum=False, ignore_photometry=True, zrange=[1.19, 2.35], nmf_toler=1.e-7)
            z = self.z_max_spec
            
        #### xxx new idea:  fit redshift grids separately to get good continua
        # self.new_fit_zgrid(dz=0.02, ignore_spectrum=False, ignore_photometry=True, zrange=[0.7, 1.5])
        # self.z_show = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        # self.fit_combined(self.z_show, nmf_toler=1.e-7, te_scale = 0.5, ignore_photometry=True)
        # 
        # model_spec = np.dot(self.coeffs, self.phot.temp_seds.T)
        # model_spec_t = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        # model_spec_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        # xf, yf = self.twod.optimal_extract(f2d)
        # xspec, yspec = self.twod.optimal_extract(model_spec_2D)
        # 
        # #### xxx Doesn't work because don't have spectrum templates
        # self.new_fit_zgrid(dz=0.02, ignore_spectrum=True, ignore_photometry=False, zrange=[0.7, 1.5])
        # self.z_show = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        # self.fit_combined(self.z_show, nmf_toler=1.e-7, te_scale = 0.5, ignore_spectrum=True)
        
        igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
        self.poly_tilt_coeffs = [1]
        self.beta_tilt_coeffs = [0]
        
        #### Only fit photometric filters around the spectrum
        orig_errors = self.phot_eflam*1.
        keep = (self.phot.lc > 0.7e4) & (self.phot.lc < 2.7e4)
        self.phot_eflam[~keep] = self.phot_flam[~keep]*100
        
        self.fit_combined(z, nmf_toler=1.e-8, te_scale = 0.5, ignore_spectrum=True)
        model_phot = np.dot(self.coeffs, self.phot.temp_seds.T)
        model_phot_t = np.dot(self.coeffs, self.templates)
        phot_coeffs = self.coeffs*1.
        
        #### restore photoemtric errors
        self.phot_eflam = orig_errors*1.
        
        #model_phot_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(selfse2D)
        
        self.fit_combined(z, nmf_toler=1.e-9, te_scale = 0.5, ignore_photometry=True)
        model_spec = np.dot(self.coeffs, self.phot.temp_seds.T)
        model_spec_t = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        model_spec_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        model_phot_2D = np.dot(phot_coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
                
        ### 
        # f2d = self.spec_flux.reshape(self.shape2D)
        # xf, yf = self.twod.optimal_extract(f2d)
        # xspec, yspec = self.twod.optimal_extract(model_spec_2D)
        # xphot, yphot = self.twod.optimal_extract(model_phot_2D)
        # xmod, ymod = self.twod.optimal_extract(model_both_2D)
        
        #### Compute the offset scaling
        subregion = (igmz > 1.1e4) & (igmz < 1.6e4)
        if self.grism_element == 'G102':
            subregion = (igmz > 0.83e4) & (igmz < 1.1e4)
            
        subregion = subregion[1:] & (np.diff(igmz) > np.percentile(np.diff(igmz[subregion]), 80.)-5) & (np.diff(igmz) > 5) 
        x, y = igmz[1:][subregion], (model_spec/model_phot)[1:][subregion]
        dy = np.append(0, np.diff(y))
        sub2 = (dy >= np.percentile(dy, 2)) & (dy <= np.percentile(dy, 98))
        x, y = x[sub2], y[sub2]
        self.poly_tilt_coeffs = polyfit(x, y, order)
        self.beta_tilt_coeffs = polyfit(np.log(x/1.4e4), np.log(y), 1)
        #ybeta = np.exp(polyval(self.beta_tilt_coeffs, np.log(x/1.4e4)))
        if np.isfinite(self.beta_tilt_coeffs).sum() != 2:
            self.beta_tilt_coeffs = [0,0]
        
        #### Skip for faint objects
        if self.twod.im[0].header['MAG'] > faint_limit:
            self.poly_tilt_coeffs = [1]
            self.beta_tilt_coeffs = [0]
            
        # f = polyval(self.tilt_coeffs, x)
        # plt.plot(x, f)
        # plt.plot(x, y)
        # plt.ylim(0.1, 2)
        # plt.semilogy()
        # 
        
        if make_figure:
            self.fit_combined(z, nmf_toler=1.e-6, te_scale = 0.5, ignore_photometry=False)
            model_both = np.dot(self.coeffs, self.phot.temp_seds.T)
            model_both_t = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
            model_both_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
            
            f2d = self.spec_flux.reshape(self.shape2D)
            xf, yf = self.twod.optimal_extract(f2d)
            xspec, yspec = self.twod.optimal_extract(model_spec_2D)
            xphot, yphot = self.twod.optimal_extract(model_phot_2D)
            xmod, ymod = self.twod.optimal_extract(model_both_2D)
            
            fig = unicorn.plotting.plot_init(xs=8, aspect=0.4, square=True, NO_GUI=NO_GUI, left=0.18)
            
            ax = fig.add_subplot(121)
            f = self.tilt_function(x)
            s = np.interp(1.4e4, x, f)
            if self.grism_element == 'G102':
                s = np.interp(1.e4, x, f)

            ax.plot(x/1.e4, y, color='black', linewidth=4)
            ax.plot(x/1.e4, f, color='white', linewidth=2, alpha=0.5)
            ax.plot(x/1.e4, f, color='red', linewidth=1)
            #ax.set_ylim(0.1, 2)
            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel('(spectrum / photometry)')
            ax.text(0.95, 0.95, r'z=%.3f, 1.4$\mu$m=%.3f' %(z, s), va='top', ha='right', transform=ax.transAxes)
            
            #for tab, xrange in zip([222, 212], [[1.0e4,1.75e4], [0.3e4, 8.e4]]):
            ax = fig.add_subplot(122)
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
            
            ax.plot(lambdaz, temp_sed, color='blue')
                            
            #### Observed photometry
            ax.errorbar(self.phot.lc, self.phot_flam, self.phot_eflam, marker='o', color='purple', markersize=8, linestyle='None', zorder=100, alpha=0.8)
            
            #### Observed spectrum
            ax.plot(xf, yf/self.oned.sens, color='red')
            ax.plot(igmz, model_spec/(1+z)**2, color='red')
            
            ax.scatter(self.phot.lc, model_spec_t, color='red', marker='s', s=50, alpha=0.5)
            w2d, f2d = self.twod.optimal_extract(model_spec_2D)
            ax.plot(w2d, f2d/self.oned.data.sensitivity, color='white', linewidth=2); ax.plot(w2d, f2d/self.oned.data.sensitivity, color='red', linewidth=1)
            
            #### Fit to photometry
            ax.plot(igmz, model_phot/(1+z)**2, color='green')
            ax.scatter(self.phot.lc, model_phot_t, color='green', marker='s', s=50, alpha=0.5)

            #### Fit to both
            ax.plot(igmz, model_both/(1+z)**2, color='black', alpha=0.4)
            ax.scatter(self.phot.lc, model_both_t, color='black', marker='s', s=50, alpha=0.5)
            ax.plot(xf, yf/self.oned.sens/self.tilt_function(self.oned.lam), color='orange')
            w2d, f2d = self.twod.optimal_extract(model_both_2D)
            f2d /= self.oned.data.sensitivity*self.tilt_function(w2d)
            ax.plot(w2d, f2d, color='white', linewidth=2); ax.plot(w2d, f2d, color='orange', linewidth=1)
            
            ax.set_xlim(0.6e4, 1.9e4)
                        
            #ymax = np.array([self.phot_flam.max(), model_phot_t.max(), self.oned.flux.max()]).max()
            #ymin = np.array([self.phot_flam[self.phot_flam > 0].min()]).min()
            ok = np.isfinite(f2d)
            ymin, ymax = f2d[ok].min(), f2d[ok].max()
            ax.set_ylim(0.5*ymin,1.5*ymax)
            ax.text(0.95, 0.95, self.root, va='top', ha='right', transform=ax.transAxes)
            #ax.semilogy()
            
            #ax.set_ylim(-0.1*ymax,2*ymax)
            #ax.semilogx()
            
            #unicorn.plotting.savefig(fig, self.root+'.zfit.tilt.png')
            unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit_tilt.%s' %(self.FIGURE_FORMAT))
            fp = open(self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit_tilt.dat','w')
            fp.write('poly\n')
            np.savetxt(fp, self.poly_tilt_coeffs)
            fp.write('beta\n')
            np.savetxt(fp, self.beta_tilt_coeffs)
            fp.close()
            
        
        return True
        
        ##### xxx all testing below

        #plt.plot(self.oned.lam, self.oned.flux/polyval(self.tilt_coeffs, self.oned.lam))
        plt.plot(igmz, model_both/(1+z)**2, color='black', alpha=0.4)
        plt.scatter(self.phot.lc, model_both_t, color='black', marker='s', s=50, alpha=0.5)
        w2d, f2d = self.twod.optimal_extract(model_both_2D)
        f2d /= self.oned.data.sensitivity*self.tilt_function(w2d)
        plt.plot(w2d, f2d, color='white', linewidth=2); plt.plot(w2d, f2d, color='orange', linewidth=1)
        
        ### Make diagnostic plot, separate function
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
        plt.plot(lambdaz, temp_sed, color='blue')
        plt.plot(self.oned.lam, self.oned.flux/self.tilt_function(self.oned.lam), color='orange')
        plt.plot(self.oned.lam, self.oned.flux, color='red')
        plt.plot(igmz, model_spec/(1+z)**2, color='red')
        w2d, f2d = self.twod.optimal_extract(model_spec_2D)
        f2d /= self.oned.data.sensitivity
        plt.plot(w2d, f2d, color='white', linewidth=2); plt.plot(w2d, f2d, color='red', linewidth=1)
        plt.scatter(self.phot.lc, model_spec_t, color='red', marker='s', s=50, alpha=0.5)
        plt.plot(igmz, model_phot/(1+z)**2, color='green')
        plt.scatter(self.phot.lc, model_phot_t, color='green', marker='s', s=50, alpha=0.5)
        plt.scatter(self.phot.lc, self.phot_flam, marker='o', color='purple', s=50)
        plt.loglog()
        plt.xlim(0.3e4, 6.e4)
        plt.ylim(0.01,10)
    
    def ln_zgrid(self, zrange=[0,4], dz=0.001):
        zgrid = np.exp(np.arange(np.log(1.0+zrange[0]), np.log(1.+zrange[1]), dz))-1
        return zgrid
        
    def new_fit_zgrid(self, zrange=[0,4], dz=0.001, ignore_photometry=False, ignore_spectrum=False, is_grid=False, apply_prior=True, nmf_toler=1.e-6, refit_norm=False):
        
        if is_grid:
            zgrid = zrange
        else:
            zgrid = self.ln_zgrid(zrange, dz)
        
        chi2 = zgrid*0.+1.e6
        for i in range(len(zgrid)):
            chi2[i] = self.fit_combined(zgrid[i], nmf_toler=nmf_toler, te_scale = 0.5, get_chi2=True, ignore_photometry=ignore_photometry, ignore_spectrum=ignore_spectrum, refit_norm=refit_norm)
            print unicorn.noNewLine+'%.4f  %.4e  [%.4f]' %(zgrid[i], chi2[i], zgrid[np.argmin(chi2)])
            #chi2[i] = self.fit_at_z(zgrid[i], nmf_toler=1.e-5, te_scale = 0.5)
                        
        pz = -0.5*chi2
        if apply_prior:
            pz += self.get_prior(zgrid)
            
        pz -= pz.max()
        
        self.zgrid = zgrid*1.
        self.lnprob_spec = pz*1.
        self.z_max_spec = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        
        # plt.plot(zgrid, pz-pz.max())
        # 
        # plt.plot(self.phot_zgrid, self.phot_lnprob-self.phot_lnprob.max(), color='green')
        # plt.plot(zgrid, prior_pz-prior_pz.max(), color='orange')
        # plt.ylim(-20,0.5)
        # zsp = self.zout.z_spec[self.ix]
        # if zsp > 0:
        #     plt.fill_between([zsp-0.003*(1+zsp), zsp+0.003*(1+zsp)], [-20,-20], [0.5,0.5], color='red', alpha=0.5)
        #     plt.plot([zsp,zsp], [-20, 0.5], color='red')
        #     plt.xlim(np.maximum(0, np.minimum(zrange[0], zsp-0.1)), np.maximum(zsp+0.1, zrange[1]))

    def new_fit_constrained(self, zrfirst=[0.,4.0], dzfirst=0.005, dzsecond=0.0005, make_plot=True, ignore_photometry=False, ignore_spectrum=False, refit_norm=False, get_tilt=True, faint_limit=24):
        import scipy.ndimage as nd
        
        print '\n First run: photometry alone \n'
        #### First fit photometry only with the new line templates
        #self.phot = self.phot_eazy
        self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=True, ignore_photometry=False, apply_prior=False, nmf_toler=1.e-5)
        self.zgrid_first = self.zgrid*1.
        self.lnprob_first_photom = self.lnprob_spec*1. #/1
        self.phot_zgrid = self.zgrid_first
        self.phot_lnprob = self.lnprob_spec*1.
        
        if self.dr > 0.5:
            self.phot_lnprob*=0
            self.lnprob_first_photom*=0
            
        #### Now fit just the spectrum
        
        self.phot = self.phot_lines
        print '\n Second run: spectrum alone \n'
        #self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=False, ignore_photometry=True, apply_prior=False, nmf_toler=1.e-5)
        #self.lnprob_first_spec = self.lnprob_spec*1.
        keep = self.phot_lnprob > -200
        self.new_fit_zgrid(zrange=self.zgrid_first[keep], is_grid=True, ignore_spectrum=False, ignore_photometry=False, apply_prior=False, nmf_toler=1.e-5)
        self.lnprob_first_spec = self.phot_lnprob*0.-200
        self.lnprob_first_spec[keep] = self.lnprob_spec*1.
        
        #### Total probability with the prior
        lnprior = self.get_prior(self.zgrid_first)
        self.lnprob_first_total = self.lnprob_first_photom + self.lnprob_first_spec + lnprior
        self.lnprob_first_total -= self.lnprob_first_total.max()
        
        #### Smooth it
        zsecond = self.ln_zgrid(zrfirst, dzsecond) #np.arange(zrfirst[0], zrfirst[1], dzsecond)
        pzint_photom = np.interp(zsecond, self.zgrid_first, self.lnprob_first_photom)
        pzint_photom -= pzint_photom.max()
        pzint_total = np.interp(zsecond, self.zgrid_first, self.lnprob_first_total)
        pzint_total -= pzint_total.max()
        
        z_max = self.zgrid_first[np.argmax(self.lnprob_first_total)]
        if get_tilt:
            self.get_spectrum_tilt(z_max, fit_for_lines=False, NO_GUI=True, faint_limit=faint_limit)
        
        min_width = 0.003*(1+z_max)
        
        #### Smooth it with a gaussian
        wz = 0.015
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm_photom = nd.convolve1d(np.exp(pzint_photom), yg, mode='constant', cval=0.)
        sm_photom = np.log(sm_photom/sm_photom.max())
        sm_total = nd.convolve1d(np.exp(pzint_total), yg, mode='constant', cval=0.)
        sm_total = np.log(sm_total/sm_total.max())
        
        zsub = sm_total > np.log(5.e-4)
        if zsub.sum() == 0:
            zsub = zsecond < zrfirst[1] #0.1
        
        #print zsecond[zsub]
        #### Fit spectrum again on finer grid
        print '\n Last run: photometry+spectrum, high z resolution (%.3f) [%.3f, %.3f]\n\n' %(z_max, zsecond[zsub].min(), zsecond[zsub].max())
        self.new_fit_zgrid(zrange=zsecond[zsub], dz=dzsecond, is_grid=True, ignore_spectrum=False, ignore_photometry=(self.dr > 0.5), apply_prior=False, refit_norm=False)
        self.zgrid_second = self.zgrid
        self.lnprob_second_spec = self.lnprob_spec*1.

        #### Smooth it with a gaussian
        wz = 0.0005
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm_spec = nd.convolve1d(np.exp(self.lnprob_second_spec*1.), yg, mode='constant', cval=0.)
        self.lnprob_second_spec = np.log(sm_spec/sm_spec.max())

        #### Total probability with the prior
        lnprior = self.get_prior(self.zgrid_second)
        self.lnprob_second_photom = sm_photom[zsub]
        self.lnprob_second_total = self.lnprob_second_photom*0 + self.lnprob_second_spec + lnprior
        self.lnprob_second_total -= self.lnprob_second_total.max()

        self.get_redshift_stats()
        
        self.z_max_spec = self.zgrid[np.argmax(self.lnprob_second_total)]

        self.norm_prob = np.exp(self.lnprob_second_total) / np.trapz(np.exp(self.lnprob_second_total), self.zgrid_second)
        self.z_peak_spec = np.trapz(self.zgrid_second*self.norm_prob, self.zgrid_second)
        self.pzsum = np.cumsum(self.norm_prob[1:]*np.diff(self.zgrid_second))
        self.c68 = np.interp([0.16,0.84], self.pzsum, self.zgrid_second[1:])
        self.c95 = np.interp([0.025,0.975], self.pzsum, self.zgrid_second[1:])

        ###### Save fit information
        fp = open(self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit_tilt.dat','w')
        fp.write('poly\n')
        np.savetxt(fp, self.poly_tilt_coeffs)
        fp.write('beta\n')
        np.savetxt(fp, self.beta_tilt_coeffs)
        fp.close()

        self.lnprob_spec = self.lnprob_second_total*1.
        
        self.new_save_results()
        
        if make_plot:
            self.make_new_fit_figure(ignore_photometry=(self.dr > 0.5), NO_GUI=True, log_pz=True)
            self.new_save_fits()
            
    def new_save_fits(self):
        import pyfits
        header = pyfits.Header()
        if self.skip_photometric:
            header.update('PHOTID', -1)
        else:
            header.update('PHOTID', self.cat.id[self.ix], comment='ID in %s' %(os.path.basename(self.cat.filename)))
        header.update('GRISID', self.grism_id, comment='Grism ID')
        
        hdu = [pyfits.PrimaryHDU(header=header)]

        #hdu.append(pyfits.ImageHDU(data=self.phot_zgrid, name='PZGRID'))
        #hdu.append(pyfits.ImageHDU(data=self.phot_lnprob, name='LN_PZ_PROB'))
                
        hdu.append(pyfits.ImageHDU(data=self.zgrid_first, name='ZGRID0'))
        hdu.append(pyfits.ImageHDU(data=self.lnprob_first_photom, name='LNP_PHOT_0'))
        hdu.append(pyfits.ImageHDU(data=self.lnprob_first_spec, name='LNP_SPEC_0'))
        hdu.append(pyfits.ImageHDU(data=self.zgrid_second, name='ZGRID1'))
        hdu.append(pyfits.ImageHDU(data=self.lnprob_second_spec, name='LNP_BOTH_1'))
        hduList = pyfits.HDUList(hdu)
        hduList.writeto(self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit.pz.fits', clobber=True, output_verify='silentfix')
    
    def new_load_fits(self, PATH='./'):
        #pzfits = '%s/%s.new_zfit.pz.fits' %(PATH, self.grism_id)
        pzfits = self.twod.file.replace('2D','new_zfit.pz')
        if not os.path.exists(pzfits):
            return False
        
        im = pyfits.open(pzfits)
        self.zgrid_first = im['ZGRID0'].data*1
        self.lnprob_first_phot = im['LNP_PHOT_0'].data*1
        self.lnprob_first_spec = im['LNP_SPEC_0'].data*1
        self.zgrid_second = im['ZGRID1'].data*1
        self.lnprob_second_spec = im['LNP_BOTH_1'].data*1
        
        self.lnprior = self.get_prior(self.zgrid_second)
        self.lnprob_second_total = self.lnprob_second_spec + self.lnprior
        self.lnprob_second_total -= self.lnprob_second_total.max()
        ### 
        self.zgrid = self.zgrid_second*1
        self.lnprob_spec = self.lnprob_second_spec*1
                
        self.get_redshift_stats()         
        
        ### for compatibility with interlace_fit.new_fit_free_emlines
        self.zgrid1 = self.zgrid_second*1
        self.full_prob1 = self.lnprob_second_total*1
        
        self.fit_combined(self.z_max_spec, nmf_toler=1.e-6, te_scale = 0.5, ignore_photometry=(self.dr > 0.5), ignore_spectrum=False, background=False)
        
        igmz, igm_factor = self.phot.get_IGM(self.z_max_spec, matrix=False)
        self.best_spec = np.dot(self.coeffs, self.phot.temp_seds.T)/(1+self.z_max_spec)**2*igm_factor
        self.best_photom = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        self.best_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        flux_2D = self.spec_flux.reshape(self.shape2D)
        xflux_1D, yflux_1D = self.twod.optimal_extract(flux_2D)
        
        self.oned_wave, self.best_1D = self.twod.optimal_extract(self.best_2D)
        self.model_1D = self.best_1D
                
        print 'Read p(z) from %s: z_max = %.4f [%.3f, %.3f]' %(pzfits, self.z_max_spec, self.c68[0], self.c68[1])
        
        linefits = self.twod.file.replace('2D.f','linefit.f')
        if os.path.exists(linefits):
            print 'Read %s' %(linefits)
            self.chain = unicorn.interlace_fit.emceeChain(file=linefits)
        
        return True
        
    def get_redshift_stats(self):
        """
        Compute best redshift and confidence intervals from p(z)
        """
        from scipy import interpolate

        self.z_max_spec = self.zgrid_second[np.argmax(self.lnprob_second_total)]
        self.norm_prob = np.exp(self.lnprob_second_total) / np.trapz(np.exp(self.lnprob_second_total), self.zgrid_second)
        self.z_peak_spec = np.trapz(self.zgrid_second*self.norm_prob, self.zgrid_second)
        
        self.pzsum = np.cumsum(self.norm_prob[1:]*np.diff(self.zgrid_second))
        self.pzsum /= self.pzsum[-1]
        self.c68 = np.interp([0.16,0.84], self.pzsum, self.zgrid_second[1:])
        self.c95 = np.interp([0.025,0.975], self.pzsum, self.zgrid_second[1:])
        self.width = [np.diff(self.c68), np.diff(self.c95)]
        
    def new_fit_in_steps(self, zrfirst=[0.,4.0], dzfirst=0.003, dzsecond=0.0005, make_plot=True, ignore_photometry=False, ignore_spectrum=False):
        import scipy.ndimage as nd
        
        # self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=True, ignore_photometry=False)
        # self.new_z_phot = self.fit_zgrid[np.argmax(self.lnprob_spec+self.get_prior(self.fit_zgrid))]
        # self.get_spectrum_tilt(self.new_z_phot)
        
        #### First loop through a coarser grid
        self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=ignore_spectrum, ignore_photometry=ignore_photometry)
        self.zgrid_first = self.zgrid*1.
        self.lnprob_first = self.lnprob_spec
        self.lnprob_first[0] = self.lnprob_first[1]
        
        zsecond = self.ln_zgrid(zrfirst, dzsecond) #np.arange(zrfirst[0], zrfirst[1], dzsecond)
        pzint = np.interp(zsecond, self.zgrid_first, self.lnprob_first)
        pzint -= pzint.max()
        
        z_max = self.zgrid_first[np.argmax(self.lnprob_first)]
        min_width = 0.003*(1+z_max)
        
        #### Smooth it with a gaussian
        wz = 0.008
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm = nd.convolve1d(np.exp(pzint), yg/yg.max(), mode='constant', cval=0.)
        
        #zsub = pzint > np.log(1.e-5)
        zsub = (sm/sm.max()) > 1.e-3 ### 3.7 sigma 
        zsub = (sm/sm.max()) > 1.e-4 ### xxx sigma 
        if zsub.sum() == 0:
            threedhst.showMessage('Something went wrong with the redshift grid...', warn=True)
            print pzint.max(), pzint.min(), self.lnprob_first.max()
            return False
        #
        if (zsecond[zsub].max() - zsecond[zsub].min()) < min_width*2:
            zrange = self.ln_zgrid([z_max-min_width, z_max + min_width], dzsecond)
        else:
            zrange = zsecond[zsub]
        
        print 'Finer grid step: %.4f (%.4f,%.4f)\n\n' %(z_max, zrange.min(), zrange.max())
        
        #### Second loop over a finer grid
        #self.new_fit_zgrid(zrange=zrange, is_grid=True, ignore_spectrum=ignore_spectrum, ignore_photometry=ignore_photometry)
        self.new_fit_zgrid(zrange=zrange, is_grid=True, ignore_spectrum=ignore_spectrum, ignore_photometry=True)
        
        wz = 0.0005
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm = np.log(nd.convolve1d(np.exp(self.lnprob_spec), yg/yg.max(), mode='constant', cval=0.))
        self.lnprob_spec = sm-sm.max()
        
        self.zgrid_second = self.zgrid*1.
        self.lnprob_second = self.lnprob_spec
        self.lnprob_second[0] = self.lnprob_second[1]
        
        if make_plot:
            self.make_new_fit_figure(ignore_photometry=True)
        
        
    def make_new_fit_figure(self, z_show=None, force_refit=True, NO_GUI=True, ignore_photometry=False, ignore_spectrum=False, log_pz=True, OUT_ROOT='new_zfit', background=0):
        
        from scipy import polyval
        
        if z_show is None:
            self.z_show = self.z_max_spec
        else:
            self.z_show = z_show
            
        if (self.coeffs is None) | force_refit:
            self.fit_combined(self.z_show, nmf_toler=1.e-6, te_scale = 0.5, ignore_photometry=ignore_photometry, ignore_spectrum=ignore_spectrum, background=background)
                         
        #
        igmz, igm_factor = self.phot.get_IGM(self.z_show, matrix=False)
        self.best_spec = np.dot(self.coeffs, self.phot.temp_seds.T)/(1+self.z_show)**2*igm_factor
        self.best_photom = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        self.best_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        flux_2D = self.spec_flux.reshape(self.shape2D)
        xflux_1D, yflux_1D = self.twod.optimal_extract(flux_2D)
         
        self.oned_wave, self.best_1D = self.twod.optimal_extract(self.best_2D)
        self.model_1D = self.best_1D
        
        #### Initialize the figure
        fig = unicorn.catalogs.plot_init(xs=10,aspect=1./3.8, left=0.1, right=0.02, bottom=0.09, top=0.08, NO_GUI=NO_GUI)

        show = self.oned.data.flux != 0.0
        #### Spectrum in e/s
        ax = fig.add_subplot(141)
        
        wuse = (self.oned.data.wave > 1.15e4) & (self.oned.data.wave < 1.6e4)
        
        if self.grism_element == 'G102':
            wuse = (self.oned.data.wave > 0.78e4) & (self.oned.data.wave < 1.15e4)
        #
        if self.grism_element == 'G800L':
            wuse = (self.oned.data.wave > 0.58e4) & (self.oned.data.wave < 0.92e4)
            #wuse = (self.oned.data.wave > 0.5e4) & (self.oned.data.wave < 1.05e4)
        
        yflux, ycont = self.oned.data.flux, self.oned.data.contam
        wuse = wuse & np.isfinite(yflux)
                
        y = yflux-ycont
        #y = yflux_1D
        
        yerr = self.oned.data.error #[show]
        ax.fill_between(self.oned.data.wave[show]/1.e4, (y+yerr-background)[show], (y-yerr-background)[show], color='blue', alpha=0.1)
        
        ax.plot(self.oned.data.wave[show]/1.e4, yflux[show]-background, color='black', alpha=0.1)
        ax.plot(self.oned.data.wave[show]/1.e4, (yflux-ycont)[show]-background, color='black')
        
        ax.plot(self.oned_wave[show]/1.e4, self.best_1D[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave/1.e4, self.best_1D, color='red', alpha=0.08, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'e$^-$ / s')
        if wuse.sum() > 5:
            ymax = yflux[wuse].max(); 
        else:
            ymax = yflux.max()
                
        ax.set_ylim(-0.05*ymax, 1.1*ymax) 
        ax.set_xlim(1.0, 1.73)

        if self.grism_element == 'G102':
            ax.set_xlim(0.74, 1.17)
        #
        if self.grism_element == 'G800L':
            ax.set_xlim(0.5, 1.05)
            
        #### Spectrum in f_lambda
        #self.oned.data.sensitivity /= 100
         
        ax = fig.add_subplot(142)
        ax.plot(self.oned.data.wave[show]/1.e4, yflux[show]/self.oned.data.sensitivity[show]*100, color='black', alpha=0.1)
        
        show_flux = (yflux[show]-ycont[show])/self.oned.data.sensitivity[show]*100
        show_err = self.oned.data.error[show]/self.oned.data.sensitivity[show]*100
        ax.plot(self.oned.data.wave[show]/1.e4, show_flux, color='black')
        #ax.fill_between(self.oned.data.wave[show]/1.e4, show_flux+show_err, show_flux-show_err, color='0.5', alpha=0.2)
        
        ax.plot(self.oned_wave/1.e4, self.best_1D/self.oned.data.sensitivity*100, color='red', alpha=0.08, linewidth=2)
        ax.plot(self.oned_wave[show]/1.e4, self.best_1D[show]/self.oned.data.sensitivity[show]*100, color='red', alpha=0.5, linewidth=2)
        #ax.plot(self.oned_wave[show]/1.e4, self.slope_1D[show]/self.oned.data.sensitivity[show], color='orange', alpha=0.2, linewidth=1)
        
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'$f_\lambda$')
        if wuse.sum() > 5:
            ymax = (yflux/self.oned.data.sensitivity*100)[wuse].max()
        else:
            ymax = (yflux/self.oned.data.sensitivity*100).max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax)
        ax.set_xlim(1.0, 1.73)

        if self.grism_element == 'G102':
            ax.set_xlim(0.74, 1.17)
        #
        if self.grism_element == 'G800L':
            #ax.set_xlim(0.58, 0.92)
            ax.set_xlim(0.5, 1.05)

        
        #### p(z)
        ax = fig.add_subplot(143)
        
        ax.plot(self.phot_zgrid, np.exp(self.phot_lnprob-self.phot_lnprob.max()), color='green')

        if self.zgrid is not None:
            ax.plot(self.zgrid, np.exp(self.lnprob_spec), color='blue', alpha=0.8)        
            #ax.fill_between(self.zgrid, np.exp(self.lnprob_spec), np.exp(self.lnprob_spec)*1.e-5, color='blue', alpha=0.2)
            
        
        zsp = -1
        if not self.skip_photometric:
            zsp = self.zout.z_spec[self.ix]
            if zsp > 0:
                ax.plot(zsp*np.array([1,1]), [1.e-8,1], color='red', alpha=0.7, linewidth=1)
        
        self.z_spec = zsp
        
        if self.dr < 1:
            ax.set_xlim(np.min([self.zgrid.min(), self.zout.l99[self.ix]]), np.max([self.zgrid.max(), self.zout.u99[self.ix]]))
        
        if log_pz:
            ax.set_ylim(1.e-3, 1)
            ax.semilogy()
            
        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'$p(z)$')
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune='both'))
         
        #### Make title text
        if not self.skip_photometric:
            deltaz = ''
            if zsp > 0:
                deltaz = '$\Delta z$ = %.4f' %(-(self.zout.z_spec[self.ix]-self.z_show)/(1+self.zout.z_spec[self.ix]))
            #
            ax.text(-0.05, 1.1, r'%s  %s  $H_{140}=$%.2f $z_\mathrm{spec}$=%.3f  $z_\mathrm{phot}$=%.3f  $z_\mathrm{gris}$=%.3f  %s' %(self.grism_id, self.zout.id[self.ix], self.twod.im[0].header['MAG'], self.zout.z_spec[self.ix], self.zout.z_peak[self.ix], self.z_show, deltaz), transform=ax.transAxes, horizontalalignment='center')
        else:
            ax.text(-0.05, 1.1, r'%s  $H_{140}=$%.2f  $z_\mathrm{gris}$=%.3f' %(self.grism_id, self.twod.im[0].header['MAG'], self.z_show), transform=ax.transAxes, horizontalalignment='center')
            unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.zfit.%s' %(self.FIGURE_FORMAT))
            return True

        ####  Show the (photometric) SED with the spectrum overplotted
        ax = fig.add_subplot(144)

        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
        
        #### Use new fit
        lambdaz, temp_sed = self.phot.templam*(1+self.z_show), self.best_spec
        obs_sed = self.best_photom
        
        temp_sed_int = np.interp(self.oned.data.wave, lambdaz, temp_sed)
        keep = (self.oned.data.wave > 1.2e4) & (self.oned.data.wave < 1.5e4)
        if self.grism_element == 'G102':
            keep = (self.oned.data.wave > 0.85) & (self.oned.data.wave < 1.05e4)
        #
        if self.grism_element == 'G800L':
            keep = (self.oned.data.wave > 0.65) & (self.oned.data.wave < 0.85e4)
        
        #flux_spec = (self.oned.data.flux-self.oned.data.contam)/self.oned.data.sensitivity*100
        flux_spec = (yflux_1D-background)/self.oned.sens*100
        
        ### factor of 100 to convert from 1.e-17 to 1.e-19 flux units
        #anorm = np.sum(temp_sed_int[keep]*flux_spec[keep])/np.sum(flux_spec[keep]**2)
        #scale = 100.
        anorm = 1
        scale = 100
        
        ax.plot(lambdaz, temp_sed*scale, color='blue', alpha=0.5, zorder=-1)
        ax.scatter(self.phot.lc, obs_sed*scale, color='purple', alpha=0.5, zorder=100)
        
        ax.errorbar(lci, fobs*scale, efobs*scale, color='black', marker='o', ms=7, alpha=0.7, linestyle='None')
       # ax.plot(self.oned.data.wave, flux_spec*anorm/polyval(self.tilt_coeffs, self.oned.data.wave)*s, color='red', alpha=0.3)
        ax.plot(self.oned.data.wave, (flux_spec)*anorm/self.tilt_function(self.oned.data.wave), color='red', alpha=0.3)
        bin = 4
        binned = unicorn.utils_c.interp_conserve(self.oned.data.wave[::4], self.oned.data.wave, flux_spec)
        ax.plot(self.oned.data.wave[::4], binned/self.tilt_function(self.oned.data.wave[::4]), color='red', alpha=0.7)
        
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$f_\lambda$')
        
        good = efobs > 0
        if good.sum() > 0:
            ymax = fobs[good].max()
        else:
            ymax = fobs.max()
        
        ymax *= scale
        
        ax.semilogx(); ax.set_xlim(3000.,8.e4); ax.set_ylim(-0.05*ymax, 1.5*ymax)
        #ax.set_xlim(0.8e4,2.e4)
        
        #### Save the result to a file
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.%s.%s' %(OUT_ROOT, self.FIGURE_FORMAT))
        
        ### Make 2D spectrum figure
        continuum_coeffs = self.coeffs*1
        continuum_coeffs[-3:-1] = 0.
        self.cont_model = np.dot(continuum_coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        self.flux_model = self.best_2D
        self.twod_figure(base='new_zfit')
        
        #### For new_free_emlines
        self.zgrid1 = self.zgrid_second
        self.full_prob1 = self.lnprob_second_total
        
    def emcee_fit(self, NWALKERS=20, NSTEP=100, verbose=True, refit=False, fit_background=False):
        """
        Fit redshift / tilt with emcee
        """
        import time
        import emcee
        
        #### Get a probability distribution for the redshift
        try:
            dummy = self.zgrid_second
        except:
            refit = True
        
        if refit:
            self.new_fit_constrained(dzfirst=0.01, dzsecond=0.001, make_plot=False)
        
        ### Draw z from p(z)
        soften = 3.
        pz = np.exp(self.lnprob_second_total/soften)
        pzsum = np.cumsum(pz[1:]*np.diff(self.zgrid_second))/np.trapz(pz, self.zgrid_second)
        z_draw = [np.interp(np.random.rand(), pzsum, self.zgrid_second[1:]) for i in range(NWALKERS)]
        
        #### normalization
        b0 = [self.beta_tilt_coeffs[0]+np.random.normal()*1 for i in range(NWALKERS)]
        b1 = [self.beta_tilt_coeffs[1]+np.random.normal()*0.2 for i in range(NWALKERS)]
        if fit_background:
            back = [np.random.normal()*0.02 for i in range(NWALKERS)]
        else:
            back = [0 for i in range(NWALKERS)]
            
        p0 = np.array([z_draw,b0,b1,back]).T
        
        obj_fun = self._objective_combined
        obj_args = [self, False, False] ## self, ignore_spectrum, ignore_photometry
        # init_pz = np.zeros(NWALKERS)
        # for i in range(NWALKERS):
        #     init_pz[i] = obj_fun(p0[i,:], *obj_args)
            
        NTHREADS=1
        ndim = 4
        self.emcee_step = 0
        
        if verbose:
            print 'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d)' %(NWALKERS, NSTEP)
            t0 = time.time()
        
        self.sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
        result = self.sampler.run_mcmc(p0, NSTEP)
        
        param_names = ['z', 'b0', 'b1', 'bkg']
        self.chain = unicorn.interlace_fit.emceeChain(chain=self.sampler.chain, param_names = param_names)
        
        if verbose:
            print unicorn.noNewLine+'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d) -> Done (%.1f s).  Make plot...' %(NWALKERS, NSTEP, time.time()-t0)
        
        yh, xh = np.histogram(self.chain['z'].flatten(), bins=80)
        self.zgrid = 0.5*(xh[1:]+xh[:-1])
        self.lnprob_spec = np.log(np.maximum(yh*1./yh.max(), 0.1/yh.max()))
        
        self.make_new_fit_figure(self.chain.median[0], ignore_photometry=False, NO_GUI=True, log_pz=True, OUT_ROOT='new_zfit_chain', background=self.chain.median[3])
        self.chain.save_fits(file='%s.new_zfit_chain.fits' %(self.root))
        
        ### test
        # p = self.chain.draw_random(N=100)
        # x = np.arange(1.1e4, 1.7e4, 100)
        # for i in range(N):
        #     self.beta_tilt_coeffs = p[i,1:3]
        #     plt.plot(x, self.tilt_function(x))
            
    @staticmethod
    def _objective_combined(params, self, ignore_spectrum, ignore_photometry):
        """
        Objective for emcee fit
        """
        z = params[0]
        if z <= 0.01:
            return -1.e10
        
        self.emcee_step += 1
        
        self.beta_tilt_coeffs = params[1:3]
        print unicorn.noNewLine, params
        chi2 = self.fit_combined(z, nmf_toler=1.e-4, te_scale=0.5, ignore_spectrum=ignore_spectrum, ignore_photometry=ignore_photometry, get_chi2=True, refit_norm=False, init_coeffs = None, background=params[3])
        print unicorn.noNewLine, self.emcee_step, params, -0.5*chi2
        
        return -0.5*chi2+self.get_prior(z)#[0]
        
    def fit_combined(self, z, nmf_toler=1.e-4, te_scale = 0.5, ignore_spectrum=False, ignore_photometry=False, get_chi2=False, refit_norm=False, init_coeffs=None, background=0):
        """
        
        Fit the spectrum AND the photometry at redshift, z
        
        Default is to store the template fit coefficients in 
        self.coeffs and the derived templates in self.templates
        
        If get_chi2==True, generate the model and compute chi2
        
        """
        from scipy import polyfit, polyval
        from unicorn import utils_c
        
        #### Masks and full template array.  If only want to fit photometry,
        #### don't make the large output arrays to save some time
        if ignore_spectrum:
            ### Use photometric template error function
            self.te = self.te_phot         
            ### Mask spectrum
            mask = self.fitting_mask['phot_only']
            ### Photometry templates
            self.templates = np.zeros((self.phot.NTEMP, self.phot.NFILT))
            full_variance = np.zeros(self.phot.NFILT)
            full_flux = np.zeros(self.phot.NFILT)   
        else:
            ### Use spectroscopic template error function
            #self.te = self.te_spec
            self.te = self.te_phot
            
            if ignore_photometry:
                mask = self.fitting_mask['spec_only']
            else:
                mask = self.fitting_mask['both']
            
            ### Photometry + spectrum templates
            self.templates = np.zeros((self.phot.NTEMP, self.phot.NFILT+np.product(self.shape2D)))
            full_variance = np.zeros(self.phot.NFILT+np.product(self.shape2D))
            full_flux = np.zeros(self.phot.NFILT+np.product(self.shape2D))
        
        #print np.sum(mask)
        
        ###### PHOTOMETRY
        ##
        #### Template error function and minimum photometric error
        #MIN_ERR = 0.025
        MIN_ERR = 0.05
        
        te_func = self.te.interpolate(self.phot.lc, z)*te_scale
        phot_var = self.phot_eflam**2+(self.phot_flam*np.maximum(te_func, MIN_ERR))**2
        
        full_variance[:self.phot.NFILT] += phot_var
        full_flux[:self.phot.NFILT] += self.phot_flam
        #### Interpolated template photometry
        self.templates[:, :self.phot.NFILT] = self.phot.interpolate_photometry(z).reshape((self.phot.NFILT, self.phot.NTEMP)).T
        
        #### Get normalization at grid point
        if refit_norm:
            self.get_spectrum_tilt(z, fit_for_lines=False, NO_GUI=False, make_figure=False)
            
        ###### SPECTRUM
        ##
        #### IGM factor
        if not ignore_spectrum:
            igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
            #### Include both IGM and computed tilt
            scale = igm_factor*self.tilt_function(igmz)/(1+z)**2
            #### Get full 2D template models
            for i in range(self.phot.NTEMP):
                self.compute_model_function(igmz, self.phot.temp_seds[:,i]*scale/self.twod.total_flux)
                self.templates[i,self.phot.NFILT:] = self.twod.model.flatten()
            #
            full_variance[self.phot.NFILT:] += self.spec_var
            full_flux[self.phot.NFILT:] += self.spec_flux-background
                            
        #### Fit Non-negative coefficients
        #### Need to seed the dusty template coefficient with a small number
        if init_coeffs is None:
            init_coeffs = np.ones(self.phot.NTEMP)
            if self.phot.NTEMP == 10:
                init_coeffs[-3:] *= 1.e-12
                
        amatrix = utils_c.prepare_nmf_amatrix(full_variance[mask], self.templates[:,mask])
        self.coeffs = utils_c.run_nmf(full_flux[mask], full_variance[mask], self.templates[:,mask], amatrix, toler=nmf_toler, MAXITER=1e6, init_coeffs=init_coeffs)
        
        #### Chi2 of the fit
        if get_chi2:
            self.model = np.dot(self.coeffs, self.templates)        
            chi2 = np.sum(((full_flux-self.model)**2/full_variance)[mask])
            return chi2
                    

    def new_save_results(self, verbose=True):
        """
        Print the results of the fit to a file, like 
        goodss-34_00622.zfit.dat.

        """
        
        self.get_redshift_stats()
        
        fp = open(self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit.dat','w')
        fp.write('#  spec_id   phot_id   mag dr   z_spec  z_peak_phot  z_max_grism z_peak_grism l95 l68 u68 u95\n')
        
        if self.skip_photometric:
            fp.write('#  Phot: None\n')
            file_string = '%s %d %.2f %.3f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n' %(self.grism_id, -1, -1, -1, -1, -1, self.z_max_spec, self.z_peak_spec, self.c95[0], self.c68[0], self.c68[1], self.c95[1])
        else:
            fp.write('#  Phot: %s\n' %(self.cat.filename))
            file_string = '%s %s %.2f %.3f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n' %(self.grism_id, self.cat.id[self.ix], self.twod.im[0].header['MAG'], self.dr, self.zout.z_spec[self.ix], self.zout.z_peak[self.ix], self.z_max_spec, self.z_peak_spec, self.c95[0], self.c68[0], self.c68[1], self.c95[1])
            
        fp.write(file_string)
        fp.close()

        if verbose:
            print file_string

    def make_2d_model(self, base='new_zfit', write_fits=True):
        
        if (self.coeffs is None):
            self.fit_combined(self.z_max_spec, nmf_toler=1.e-6, te_scale = 0.5, ignore_photometry=False, ignore_spectrum=False, background=0.)
            print 'Refitting to get coefficients.'
        
        continuum_coeffs = self.coeffs*1
        continuum_coeffs[-3:-1] = 0.
        self.cont_model = np.dot(continuum_coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        line_coeffs = self.coeffs*1.
        msk = np.array(line_coeffs, dtype=bool)
        msk[-3:-1] = False
        line_coeffs[msk] = 0.
        self.line_model = np.dot(line_coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        oned_wave_x, self.cont_1D = self.twod.optimal_extract(self.cont_model)
        oned_wave_x, self.line_1D = self.twod.optimal_extract(self.line_model)
    
        if write_fits:
            self.make_spectrum_fits(base=base)
        
        
"""
Tests on the UDF interlaced data
"""

import numpy as np
import matplotlib.pyplot as plt

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
    
    
###### New simultaneous redshift fitting

class SimultaneousFit(unicorn.interlace_fit.GrismSpectrumFit):
    def read_master_templates(self):
        """
        Read the eazy-interpolated continuum + emission line templates
        """
        import threedhst.eazyPy as eazy
        
        #### Read template error function
        self.te = eazy.TemplateError(unicorn.GRISM_HOME+'templates/TEMPLATE_ERROR.eazy_v1.0')
        
        #### Read template photometry
        #self.phot = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='eazy_v1.1', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.zout)
        self.phot = eazy.TemplateInterpolator(None, MAIN_OUTPUT_FILE='full', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/templates/FULL_EAZY/OUTPUT', zout=self.zout)
        
        #### SPECTRUM
        self.tilt_coeffs = [1]
        self.shape2D = self.twod.im['SCI'].data.shape
        
        ## Flatten the 2D spectroscopic flux and variance arrays
        ## Variance
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        var[var == 0] = 1.e6
        var += (0.1*self.twod.im['CONTAM'].data**2).flatten()
        self.spec_var = var
        
        ## Flux
        self.spec_flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        self.spec_use = np.isfinite(self.spec_flux) & (self.twod.im['WHT'] != 0)
                
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
                   
        #
        self.coeffs = np.ones(self.phot.NTEMP)
        self.zgrid = None
        self.lnprob_spec = None
         
    def get_prior(self, zgrid):
        """
        Get EAZY apparent mag prior
        """
        prior = open('%s/data/prior_K_zmax7_coeff.dat' %(os.path.dirname(threedhst.__file__))).readlines()
        prior_mag = np.cast[float](prior[0].split()[2:])
        prior_z0 = np.cast[float](prior[1].split()[1:])
        prior_gamma = np.cast[float](prior[2].split()[1:])
        z0 = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_z0, left=prior_z0[0], right=prior_z0[-1])
        gamma = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_gamma, left=prior_gamma[0], right=prior_gamma[-1])
        
        prior_pz = np.log(zgrid**gamma * np.exp(-(zgrid / z0)**gamma))
        prior_pz[~np.isfinite(prior_pz)] = -300
        prior_pz -= prior_pz.max()
        
        return prior_pz
            
    def get_spectrum_tilt(self, z=None, make_figure=True, order=6):
        from scipy import polyfit, polyval
        
        if z is None:
            z = self.z_peak
        
        #### xxx new idea:  fit redshift grids separately to get good continua
        self.new_fit_zgrid(dz=0.02, ignore_spectrum=False, ignore_photometry=True, zrange=[0.7, 1.5])
        self.z_show = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        self.fit_combined(self.z_show, nmf_toler=1.e-7, te_scale = 0.5, ignore_photometry=True)
        
        model_spec = np.dot(self.coeffs, self.phot.temp_seds.T)
        model_spec_t = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        model_spec_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        xf, yf = self.twod.optimal_extract(f2d)
        xspec, yspec = self.twod.optimal_extract(model_spec_2D)
        
        #### xxx Doesn't work because don't have spectrum templates
        self.new_fit_zgrid(dz=0.02, ignore_spectrum=True, ignore_photometry=False, zrange=[0.7, 1.5])
        self.z_show = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        self.fit_combined(self.z_show, nmf_toler=1.e-7, te_scale = 0.5, ignore_spectrum=True)
        
        igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
        self.tilt_coeffs = [1]
        
        #### Only fit photometric filters around the spectrum
        orig_errors = self.phot_eflam*1.
        keep = (self.phot.lc > 1.e4) & (self.phot.lc < 1.7e4)
        self.phot_eflam[~keep] = self.phot_flam[~keep]*100
        
        self.fit_combined(z, nmf_toler=1.e-6, te_scale = 0.5, ignore_spectrum=True)
        model_phot = np.dot(self.coeffs, self.phot.temp_seds.T)
        model_phot_t = np.dot(self.coeffs, self.templates)
        phot_coeffs = self.coeffs*1.
        
        #### restore photoemtric errors
        self.phot_eflam = orig_errors*1.
        
        #model_phot_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        self.fit_combined(z, nmf_toler=1.e-7, te_scale = 0.5, ignore_photometry=True)
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
        subregion = (igmz > 1.e4) & (igmz < 1.8e4)
        subregion = subregion[1:] & (np.diff(igmz) > np.max(np.diff(igmz[subregion]))-5) 
        x, y = igmz[1:][subregion], (model_spec/model_phot)[1:][subregion]
        dy = np.append(0, np.diff(y))
        sub2 = (dy >= np.percentile(dy, 2)) & (dy <= np.percentile(dy, 98))
        x, y = x[sub2], y[sub2]
        self.tilt_coeffs = polyfit(x, y, order)
        
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
            
            fig = unicorn.plotting.plot_init(xs=6, aspect=1, square=True)
            
            ax = fig.add_subplot(221)
            f = polyval(self.tilt_coeffs, x)
            s = np.interp(1.4e4, x, f)
            ax.plot(x, y, color='black', linewidth=4)
            ax.plot(x, f, color='white', linewidth=2, alpha=0.5)
            ax.plot(x, f, color='red', linewidth=1)
            #ax.set_ylim(0.1, 2)
            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel('Template fit (spectrum / photometry)')
            ax.text(0.95, 0.95, r'z=%.3f, 1.4$\mu$m=%.3f' %(z, s), va='top', ha='right', transform=ax.transAxes)
            
            for tab, xrange in zip([222, 212], [[1.0e4,1.75e4], [0.3e4, 8.e4]]):
                ax = fig.add_subplot(tab)
                lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
                
                ax.plot(lambdaz, temp_sed, color='blue')
                                
                #### Observed photometry
                ax.errorbar(self.phot.lc, self.phot_flam, self.phot_eflam, marker='o', color='purple', markersize=8, linestyle='None')
                
                #### Observed spectrum
                ax.plot(xf, yf/self.oned.sens, color='red')
                ax.plot(igmz, model_spec/(1+z)**2, color='red')
                
                ax.scatter(self.phot.lc, model_spec_t, color='red', marker='s', s=50, alpha=0.5)
                w2d, f2d = self.twod.optimal_extract(model_spec_2D)
                ax.plot(w2d, f2d/self.oned.data.sensitivity, color='white', linewidth=2); plt.plot(w2d, f2d/self.oned.data.sensitivity, color='red', linewidth=1)
                
                #### Fit to photometry
                ax.plot(igmz, model_phot/(1+z)**2, color='green')
                ax.scatter(self.phot.lc, model_phot_t, color='green', marker='s', s=50, alpha=0.5)

                #### Fit to both
                ax.plot(igmz, model_both/(1+z)**2, color='black', alpha=0.4)
                ax.scatter(self.phot.lc, model_both_t, color='black', marker='s', s=50, alpha=0.5)
                ax.plot(xf, yf/self.oned.sens/polyval(self.tilt_coeffs, self.oned.lam), color='orange')
                w2d, f2d = self.twod.optimal_extract(model_both_2D)
                f2d /= self.oned.data.sensitivity*polyval(self.tilt_coeffs, w2d)
                ax.plot(w2d, f2d, color='white', linewidth=2); plt.plot(w2d, f2d, color='orange', linewidth=1)
                
                ax.set_xlim(xrange)
                ymax = np.array([self.phot_flam.max(), model_phot_t.max(), self.oned.flux.max()]).max()
                ymin = np.array([self.phot_flam[self.phot_flam > 0].min()]).min()
                ax.set_ylim(0.5*ymin,2*ymax)
                ax.text(0.95, 0.95, self.root, va='top', ha='right', transform=ax.transAxes)
                ax.semilogy()
            
            ax.set_ylim(-0.1*ymax,2*ymax)
            ax.semilogx()
            
            #unicorn.plotting.savefig(fig, self.root+'.zfit.tilt.png')
            unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.zfit_tilt.%s' %(self.FIGURE_FORMAT))
            
        
        return True
        

        #plt.plot(self.oned.lam, self.oned.flux/polyval(self.tilt_coeffs, self.oned.lam))
        plt.plot(igmz, model_both/(1+z)**2, color='black', alpha=0.4)
        plt.scatter(self.phot.lc, model_both_t, color='black', marker='s', s=50, alpha=0.5)
        w2d, f2d = self.twod.optimal_extract(model_both_2D)
        f2d /= self.oned.data.sensitivity*polyval(self.tilt_coeffs, w2d)
        plt.plot(w2d, f2d, color='white', linewidth=2); plt.plot(w2d, f2d, color='orange', linewidth=1)
        
        ### Make diagnostic plot, separate function
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
        plt.plot(lambdaz, temp_sed, color='blue')
        plt.plot(self.oned.lam, self.oned.flux/polyval(self.tilt_coeffs, self.oned.lam), color='orange')
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
        
    def new_fit_zgrid(self, zrange=[0,4], dz=0.001, ignore_photometry=False, ignore_spectrum=False, is_grid=False):
        
        if is_grid:
            zgrid = zrange
        else:
            zgrid = self.ln_zgrid(zrange, dz)
        
        chi2 = zgrid*0.
        for i in range(len(zgrid)):
            chi2[i] = self.fit_combined(zgrid[i], nmf_toler=1.e-6, te_scale = 0.5, get_chi2=True, ignore_photometry=ignore_photometry, ignore_spectrum=ignore_spectrum)
            print unicorn.noNewLine+'%.4f  %.4e' %(zgrid[i], chi2[i])
            #chi2[i] = self.fit_at_z(zgrid[i], nmf_toler=1.e-5, te_scale = 0.5)
                        
        pz = -0.5*chi2
        pz -= pz.max()
        
        self.zgrid = zgrid*1.
        self.lnprob_spec = pz*1.
        
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
    
    def new_fit_in_steps(self, zrfirst=[0.,3.8], dzfirst=0.003, dzsecond=0.0005, make_plot=True, ignore_photometry=False, ignore_spectrum=False):
        import scipy.ndimage as nd
        
        # self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=True, ignore_photometry=False)
        # self.new_z_phot = self.fit_zgrid[np.argmax(self.lnprob_spec+self.get_prior(self.fit_zgrid))]
        # self.get_spectrum_tilt(self.new_z_phot)
        
        #### First loop through a coarser grid
        self.new_fit_zgrid(zrange=zrfirst, dz=dzfirst, ignore_spectrum=ignore_spectrum, ignore_photometry=ignore_photometry)
        self.zgrid_first = self.zgrid*1.
        self.fit_lnprob_first = self.lnprob_spec+self.get_prior(self.zgrid_first)
        self.fit_lnprob_first[0] = self.fit_lnprob_first[1]
        
        zsecond = self.ln_zgrid(zrfirst, dzsecond) #np.arange(zrfirst[0], zrfirst[1], dzsecond)
        pzint = np.interp(zsecond, self.zgrid_first, self.fit_lnprob_first)
        pzint -= pzint.max()
        
        z_max = self.zgrid_first[np.argmax(self.fit_lnprob_first)]
        min_width = 0.003*(1+z_max)
        
        #### Smooth it with a gaussian
        wz = 0.008
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm = nd.convolve1d(np.exp(pzint), yg/yg.max(), mode='constant', cval=0.)
        
        #zsub = pzint > np.log(1.e-5)
        zsub = (sm/sm.max()) > 1.e-3 ### 3.7 sigma 
        if zsub.sum() == 0:
            threedhst.showMessage('Something went wrong with the redshift grid...', warn=True)
            print pzint.max(), pzint.min(), self.fit_lnprob_first.max()
            return False
        #
        if (zsecond[zsub].max() - zsecond[zsub].min()) < min_width*2:
            zrange = self.ln_zgrid([z_max-min_width, z_max + min_width], dzsecond)
        else:
            zrange = zsecond[zsub]
        
        print 'Finer grid step: %.4f (%.4f,%.4f)\n\n' %(z_max, zrange.min(), zrange.max())
        
        #### Second loop over a finer grid
        self.new_fit_zgrid(zrange=zrange, is_grid=True, ignore_spectrum=ignore_spectrum, ignore_photometry=ignore_photometry)
        
        wz = 0.0005
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2/wz**2)
        sm = np.log(nd.convolve1d(np.exp(self.lnprob_spec), yg/yg.max(), mode='constant', cval=0.))
        self.lnprob_spec = sm-sm.max()
        
        self.zgrid_second = self.zgrid*1.
        self.fit_lnprob_second = self.lnprob_spec + self.get_prior(self.zgrid_second)
        self.fit_lnprob_second[0] = self.fit_lnprob_second[1]
        
        if make_plot:
            self.make_new_fit_figure()
        
        
    def make_new_fit_figure(self, z_show=None):
        
        from scipy import polyval
        
        if z_show is None:
            self.z_show = self.zgrid[self.lnprob_spec == self.lnprob_spec.max()][0]
        else:
            self.z_show = z_show
            
        if self.coeffs is None:
            self.fit_combined(self.z_show, nmf_toler=1.e-6, te_scale = 0.5, ignore_photometry=False)
        
        self.best_spec = np.dot(self.coeffs, self.phot.temp_seds.T)/(1+self.z_show)**2
        self.best_photom = np.dot(self.coeffs, self.templates)[:self.phot.NFILT]
        self.best_2D = np.dot(self.coeffs, self.templates)[self.phot.NFILT:].reshape(self.shape2D)
        
        self.oned_wave, self.best_1D = self.twod.optimal_extract(self.best_2D)
        
        #### Initialize the figure
        fig = unicorn.catalogs.plot_init(xs=10,aspect=1./3.8, left=0.1, right=0.02, bottom=0.09, top=0.08, NO_GUI=False)

        show = self.oned.data.flux != 0.0
        #### Spectrum in e/s
        ax = fig.add_subplot(141)
        
        wuse = (self.oned.data.wave > 1.15e4) & (self.oned.data.wave < 1.6e4)
        
        if self.grism_element == 'G102':
            wuse = (self.oned.data.wave > 0.78e4) & (self.oned.data.wave < 1.15e4)
        
        yflux, ycont = self.oned.data.flux, self.oned.data.contam
        y = yflux-ycont
        
        yerr = self.oned.data.error #[show]
        ax.fill_between(self.oned.data.wave[show]/1.e4, (y+yerr)[show], (y-yerr)[show], color='blue', alpha=0.1)
        
        ax.plot(self.oned.data.wave[show]/1.e4, yflux[show], color='black', alpha=0.1)
        ax.plot(self.oned.data.wave[show]/1.e4, (yflux-ycont)[show], color='black')
        
        ax.plot(self.oned_wave[show]/1.e4, self.best_1D[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave/1.e4, self.best_1D, color='red', alpha=0.08, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'e$^-$ / s')
        if wuse.sum() > 5:
            ymax = yflux[wuse].max(); 
        else:
            ymax = yflux.max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax) 
        if self.grism_element == 'G102':
            ax.set_xlim(0.7, 1.15)
        else:
            ax.set_xlim(1.0, 1.73)
            
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
        if self.grism_element == 'G102':
            ax.set_xlim(0.7, 1.15)
        else:
            ax.set_xlim(1.0, 1.73)

        
        #### p(z)
        ax = fig.add_subplot(143)
        ax.plot(self.phot_zgrid, np.exp(self.phot_lnprob-self.phot_lnprob.max()), color='green')
        if self.zgrid is not None:
            ax.plot(self.zgrid, np.exp(self.lnprob_spec), color='blue', alpha=0.4)
            ax.fill_between(self.zgrid, np.exp(self.lnprob_spec), self.zgrid*0., color='blue', alpha=0.2)
            
        #ax.plot(zgrid1, np.exp(full_prob1), color='blue')

        # ax.plot(self.phot_zgrid, (self.phot_lnprob-self.phot_lnprob.max()), color='green')
        # ax.plot(zgrid0, (full_prob0), color='blue', alpha=0.4)
        # ax.plot(zgrid1, (full_prob1), color='blue')
        
        zsp = -1
        if not self.skip_photometric:
            zsp = self.zout.z_spec[self.ix]
            if zsp > 0:
                ax.plot(zsp*np.array([1,1]), [0,1], color='red', alpha=0.9, linewidth=2)

        if self.dr < 1:
            ax.set_xlim(np.min([self.zgrid.min(), self.zout.l99[self.ix]]), np.max([self.zgrid.max(), self.zout.u99[self.ix]]))

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
        
        flux_spec = (self.oned.data.flux-self.oned.data.contam)/self.oned.data.sensitivity*100
        
        ### factor of 100 to convert from 1.e-17 to 1.e-19 flux units
        #anorm = np.sum(temp_sed_int[keep]*flux_spec[keep])/np.sum(flux_spec[keep]**2)
        #scale = 100.
        anorm = 1
        scale = 100
        
        ax.plot(lambdaz, temp_sed*scale, color='blue', alpha=0.5, zorder=-1)
        ax.scatter(self.phot.lc, obs_sed*scale, color='purple', alpha=0.5, zorder=100)
        
        ax.errorbar(lci, fobs*scale, efobs*scale, color='black', marker='o', ms=7, alpha=0.7, linestyle='None')
        s = polyval(self.tilt_coeffs, 1.4e4)
        ax.plot(self.oned.data.wave, flux_spec*anorm/polyval(self.tilt_coeffs, self.oned.data.wave)*s, color='red', alpha=0.3)
        bin = 4
        binned = unicorn.utils_c.interp_conserve(self.oned.data.wave[::4], self.oned.data.wave, flux_spec)
        ax.plot(self.oned.data.wave[::4], binned/polyval(self.tilt_coeffs, self.oned.data.wave[::4])*s, color='red', alpha=0.7)
        
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
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.new_zfit.%s' %(self.FIGURE_FORMAT))
        
    def fit_combined(self, z, nmf_toler=1.e-4, te_scale = 0.5, ignore_spectrum=False, ignore_photometry=False, get_chi2=False):
        
        from scipy import polyfit, polyval
        from unicorn import utils_c
        
        #### Masks and full template array.  If only want to fit photometry,
        #### don't make the large output arrays to save some time
        if ignore_spectrum:
            mask = self.fitting_mask['phot_only']
            self.templates = np.zeros((self.phot.NTEMP, self.phot.NFILT))
            full_variance = np.zeros(self.phot.NFILT)
            full_flux = np.zeros(self.phot.NFILT)            
        else:
            if ignore_photometry:
                mask = self.fitting_mask['spec_only']
            else:
                mask = self.fitting_mask['both']
            
            self.templates = np.zeros((self.phot.NTEMP, self.phot.NFILT+np.product(self.shape2D)))
            full_variance = np.zeros(self.phot.NFILT+np.product(self.shape2D))
            full_flux = np.zeros(self.phot.NFILT+np.product(self.shape2D))
             
        print np.sum(mask)
        
        ###### PHOTOMETRY
        ##
        #### Template error function
        te_func = self.te.interpolate(self.phot.lc, z)*te_scale
        phot_var = self.phot_eflam**2+(self.phot_flam*np.maximum(te_func, 0.025))**2
        full_variance[:self.phot.NFILT] += phot_var
        full_flux[:self.phot.NFILT] += self.phot_flam
        #### Interpolated template photometry
        self.templates[:, :self.phot.NFILT] = self.phot.interpolate_photometry(z).reshape((self.phot.NFILT, self.phot.NTEMP)).T
        
        ###### SPECTRUM
        ##
        #### IGM factor
        if not ignore_spectrum:
            igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
            #### Include both IGM and computed tilt
            scale = igm_factor*polyval(self.tilt_coeffs, igmz)/(1+z)**2
            #### Get full 2D template models
            for i in range(self.phot.NTEMP):
                self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*scale/self.twod.total_flux)
                self.templates[i,self.phot.NFILT:] = self.twod.model.flatten()
            #
            full_variance[self.phot.NFILT:] += self.spec_var
            full_flux[self.phot.NFILT:] += self.spec_flux
            
        #### Fit Non-negative coefficients
        amatrix = utils_c.prepare_nmf_amatrix(full_variance[mask], self.templates[:,mask])
        self.coeffs = utils_c.run_nmf(full_flux[mask], full_variance[mask], self.templates[:,mask], amatrix, toler=nmf_toler, MAXITER=1e6)
        
        #### Chi2 of the fit
        if get_chi2:
            self.model = np.dot(self.coeffs, self.templates)        
            chi2 = np.sum(((full_flux-self.model)**2/full_variance)[mask])
            return chi2
            
    # def fit_spec(self, z, nmf_toler=1.e-5, te_scale = 0.5, tilt_coeffs = [1]):
    #     
    #     ####
    #     ######## Get template spectra at redshift, z
    #     ####
    # 
    #     #### Full flattened template array
    #     template_spec = np.zeros((self.phot.NTEMP, self.flux.shape[0]))
    #     #### IGM factor
    #     igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
    #     #### Include both IGM and computed tilt
    #     scale = igm_factor*polyval(tilt_coeffs, igmz)
    #     #### Get full 2D template models
    #     for i in range(self.phot.NTEMP):
    #         self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*scale)
    #         template_spec[i,:] = self.twod.model.flatten()
    #             
    #     #### Fit Non-negative coefficients
    #     amatrix = utils_c.prepare_nmf_amatrix(self.spec_var[self.spec_use], template_spec[:,self.spec_use])
    #     coeffs = utils_c.run_nmf(flux[self.spec_use], self.spec_var[self.spec_use], template_spec[:,self.spec_use], amatrix, toler=nmf_toler, MAXITER=1e6)
    #     
    #     #### Store results
    #     self.spec_coeffs = coeffs
    #     self.model_spec_full = np.dot(coeffs, self.phot.temp_seds.T)*igm_factor
    #     self.model_spec = np.dot(coeffs, template_spec)        
    #     
    #     #### Chi2 of the fit
    #     self.chi2_spec = np.sum(((flux-model_spec)**2/var)[self.spec_use])
    # 
    #     subregion = (np.diff(igmz) > 1) & (igmz[1:] > 1.e4) & (igmz[1:] < 1.8e4)
    #     
    #     x, y = igmz[1:][subregion], (model_spec_full/model_phot_full)[1:][subregion]
    #     
    #     from scipy import polyfit, polyval
    #     p = polyfit(x, y, 6)
    #     f = polyval(p, x)
    #     plt.plot(x, f)
    #     plt.plot(x, y)
    #     
    #     for i in range(self.phot.NTEMP):
    #         self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*igm_factor*polyval(p, igmz))
    #         template_spec[i,:] = self.twod.model.flatten()
    #                                                  
    #     #
    #     amatrix = utils_c.prepare_nmf_amatrix(var[use], template_spec[:,use])
    #     coeffs = utils_c.run_nmf(flux[use], var[use], template_spec[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)
    # 
    #     model_spec2 = np.dot(coeffs, template_spec)        
    #     
    #     # model_spec = np.dot(coeffs, template_spec).reshape(self.twod.im['SCI'].data.shape)
    #     # v2d = var.reshape(model_spec.shape)
    #     # u2d = use.reshape(v2d.shape)
    #     # f2d = flux.reshape(v2d.shape)
    #     # 
    #     # shifts = np.arange(-1,1,0.1)
    #     # chi2 = shifts*0.
    #     # for i in range(len(shifts)):
    #     #     sh = nd.shift(model_spec, (shifts[i], 0))
    #     #     chi2[i] = np.sum(((f2d-sh)**2/v2d)[u2d])
    #     # 
    #     # best_yshift = shifts[chi2 == chi2.min()][0]
    #     # sh = nd.shift(model_spec, (best_yshift, 0))
    #     # 
    #     # flux = np.cast[float](nd.shift(self.twod.im['SCI'].data, (-best_yshift, 0))-self.twod.im['CONTAM'].data).flatten()
    #     # coeffs = utils_c.run_nmf(flux[use], var[use], template_spec[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)
    #     # model_spec2 = np.dot(coeffs, template_spec).reshape(self.twod.im['SCI'].data.shape)
    #     
    #     
    # def fit_at_z(self, z, nmf_toler=1.e-5, te_scale = 0.5):
    #     
    #     #### Template error
    #     te_func = self.te.interpolate(lci, z)*te_scale
    #     phot_var = self.efobs_fnu**2+(self.fobs_fnu*te_func)**2
    #     
    #     template_phot = self.phot.interpolate_photometry(z).reshape((self.phot.NFILT, self.phot.NTEMP)).T #[:-3,:]
    #     
    #     #tfull = self.phot.interpolate_photometry(zgrid)
    #     #i = 148
    #     #template_phot = tfull[:,:,i].reshape((self.phot.NFILT, self.phot.NTEMP)).T
    #     
    #     #### Fit the template to photometry
    #     amatrix = utils_c.prepare_nmf_amatrix(var[use], template_phot[:,use])
    #     coeffs = utils_c.run_nmf(fobs_fnu[use], var[use], template_phot[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)
    #     model_phot = np.dot(coeffs, template_phot)
    #     igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
    #     model_phot_full = np.dot(coeffs, self.phot.temp_seds.T)*igm_factor
    #     
    #     chi2 = np.sum((model_phot[use]-fobs_fnu[use])**2/var[use])
    #     chi2_obs = np.sum((model_phot[use]-obs_sed[use])**2/var[use])
    #     resid =  (model_phot[use]-fobs_fnu[use])/np.sqrt(var[use])
    #     #print resid
    #     return chi2
    #     
    #     #plt.errorbar(lci, fobs_fnu/(lci/5500.)**2, efobs_fnu/(lci/5500.)**2, color='blue')
    #     # plt.errorbar(lci, model_phot-fobs_fnu, efobs_fnu, color='red', marker='o', alpha=0.5)
    #     # plt.errorbar(lci, obs_sed*(lci/5500.)**2-fobs_fnu, efobs_fnu, color='green', marker='o', alpha=0.5)
    #     plt.errorbar(lci, fobs_fnu, efobs_fnu, color='blue', marker='o', alpha=0.5)
    #     plt.scatter(lci, obs_sed*(lci/5500.)**2, color='green', marker='o', alpha=0.5)
    #     plt.scatter(lci, model_phot, color='red', marker='o', alpha=0.5)
        

    
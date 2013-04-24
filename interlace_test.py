"""
Tests on the UDF interlaced data
"""

import numpy as np

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
        self.phot_fnu = fobs*(lci/5500.)**2
        self.phot_efnu = efobs*(lci/5500.)**2
        self.phot_use = (fobs > -99) & (efobs > 0)
        
        #### Masks for fitting the photometry and/or spectra
        self.fitting_mask = {}
        self.fitting_mask['phot_only'] = self.phot_use
        
        self.fitting_mask['spec_only'] = np.zeros(self.phot.NFILT + np.size(self.spec_flux), dtype=np.bool)
        self.fitting_mask['spec_only'][self.phot.NFILT:] |= self.spec_use
        
        self.fitting_mask['both'] = np.zeros(self.phot.NFILT + np.size(self.spec_flux), dtype=np.bool) | self.fitting_mask['spec_only']
        self.fitting_mask['both'][:self.phot.NFILT] |= self.phot_use
                    
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
        prior_pz -= prior_pz.max()
        return prior_pz
            
    def get_spectrum_tilt(self, z=None):
        from scipy import polyfit, polyval
        
        if z is None:
            z = self.z_peak
        
        igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
        self.tilt_coeffs = [1]
         
        self.fit_combined(z, nmf_toler=1.e-4, te_scale = 0.5, ignore_spectrum=True)
        model_phot = np.dot(self.coeffs, self.phot.temp_seds.T)
        
        self.fit_combined(z, nmf_toler=1.e-4, te_scale = 0.5, ignore_photometry=True)
        model_spec = np.dot(self.coeffs, self.phot.temp_seds.T)
        
        #### Compute the offset scaling
        subregion = (np.diff(igmz) > 1) & (igmz[1:] > 1.e4) & (igmz[1:] < 1.8e4)
        x, y = igmz[1:][subregion], 1./(model_spec/model_phot)[1:][subregion]
        
        self.tilt_coeffs = polyfit(x, y, 6)
        f = polyval(self.tilt_coeffs, x)
        plt.plot(x, f)
        plt.plot(x, y)
        
        return True
        
        self.fit_combined(z, nmf_toler=1.e-4, te_scale = 0.5, ignore_photometry=False)
        model_both = np.dot(self.coeffs, self.phot.temp_seds.T)
        plt.plot(self.oned.lam, self.oned.flux)
        #plt.plot(self.oned.lam, self.oned.flux/polyval(self.tilt_coeffs, self.oned.lam))
        plt.plot(igmz, model_both)
        plt.plot(self.phot.lc, fobs, marker='o')
        plt.plot(lambdaz, temp_sed)
        
        ### Make diagnostic plot, separate function
        plt.plot(lambdaz, temp_sed)
        plt.plot(self.oned.lam, self.oned.flux/polyval(self.tilt_coeffs, self.oned.lam))
        plt.plot(self.oned.lam, self.oned.flux)
        plt.plot(igmz, model_spec/(1+z)**2)
        plt.plot(igmz, model_phot/(1+z)**2)
        plt.plot(self.phot.lc, fobs, linestyle='None', marker='o')
        plt.loglog()
        plt.xlim(0.8e4, 2.e4)
        plt.ylim(0.1,10)
    
    def fit_photz(self, zrange=[0,4], dz=0.005):
        zgrid = np.exp(np.arange(np.log(1.0+zrange[0]), np.log(1.+zrange[1]), dz))-1
        chi2 = zgrid*0.
        for i in range(len(zgrid)):
            chi2[i] = self.fit_combined(zgrid[i], nmf_toler=1.e-4, te_scale = 0.5, get_chi2=True, ignore_photometry=True)
            print zgrid[i], chi2[i]
            #chi2[i] = self.fit_at_z(zgrid[i], nmf_toler=1.e-5, te_scale = 0.5)
        
        prior_pz = self.get_prior(zgrid)
                
        pz = -0.5*chi2+prior_pz
        plt.plot(self.phot_zgrid, self.phot_lnprob-self.phot_lnprob.max())
        plt.plot(zgrid, pz-pz.max())
        plt.plot(zgrid, prior_pz-prior_pz.max())
        plt.ylim(-20,0.5)
    
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
             
        
        ###### PHOTOMETRY
        ##
        #### Template error function
        te_func = self.te.interpolate(self.phot.lc, z)*te_scale
        phot_var = self.phot_efnu**2+(self.phot_fnu*te_func)**2
        full_variance[:self.phot.NFILT] += phot_var
        full_flux[:self.phot.NFILT] += self.phot_fnu
        #### Interpolated template photometry
        self.templates[:, :self.phot.NFILT] = self.phot.interpolate_photometry(z).reshape((self.phot.NFILT, self.phot.NTEMP)).T
        
        ###### SPECTRUM
        ##
        #### IGM factor
        if not ignore_spectrum:
            igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
            #### Include both IGM and computed tilt
            scale = igm_factor*polyval(self.tilt_coeffs, igmz)
            #### Get full 2D template models
            for i in range(self.phot.NTEMP):
                self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*scale)
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
            
    def fit_spec(self, z, nmf_toler=1.e-5, te_scale = 0.5, tilt_coeffs = [1]):
        
        ####
        ######## Get template spectra at redshift, z
        ####

        #### Full flattened template array
        template_spec = np.zeros((self.phot.NTEMP, self.flux.shape[0]))
        #### IGM factor
        igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
        #### Include both IGM and computed tilt
        scale = igm_factor*polyval(tilt_coeffs, igmz)
        #### Get full 2D template models
        for i in range(self.phot.NTEMP):
            self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*scale)
            template_spec[i,:] = self.twod.model.flatten()
                
        #### Fit Non-negative coefficients
        amatrix = utils_c.prepare_nmf_amatrix(self.spec_var[self.spec_use], template_spec[:,self.spec_use])
        coeffs = utils_c.run_nmf(flux[self.spec_use], self.spec_var[self.spec_use], template_spec[:,self.spec_use], amatrix, toler=nmf_toler, MAXITER=1e6)
        
        #### Store results
        self.spec_coeffs = coeffs
        self.model_spec_full = np.dot(coeffs, self.phot.temp_seds.T)*igm_factor
        self.model_spec = np.dot(coeffs, template_spec)        
        
        #### Chi2 of the fit
        self.chi2_spec = np.sum(((flux-model_spec)**2/var)[self.spec_use])

        subregion = (np.diff(igmz) > 1) & (igmz[1:] > 1.e4) & (igmz[1:] < 1.8e4)
        
        x, y = igmz[1:][subregion], (model_spec_full/model_phot_full)[1:][subregion]
        
        from scipy import polyfit, polyval
        p = polyfit(x, y, 6)
        f = polyval(p, x)
        plt.plot(x, f)
        plt.plot(x, y)
        
        for i in range(self.phot.NTEMP):
            self.twod.compute_model(igmz, self.phot.temp_seds[:,i]*igm_factor*polyval(p, igmz))
            template_spec[i,:] = self.twod.model.flatten()
                                                     
        #
        amatrix = utils_c.prepare_nmf_amatrix(var[use], template_spec[:,use])
        coeffs = utils_c.run_nmf(flux[use], var[use], template_spec[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)

        model_spec2 = np.dot(coeffs, template_spec)        
        
        # model_spec = np.dot(coeffs, template_spec).reshape(self.twod.im['SCI'].data.shape)
        # v2d = var.reshape(model_spec.shape)
        # u2d = use.reshape(v2d.shape)
        # f2d = flux.reshape(v2d.shape)
        # 
        # shifts = np.arange(-1,1,0.1)
        # chi2 = shifts*0.
        # for i in range(len(shifts)):
        #     sh = nd.shift(model_spec, (shifts[i], 0))
        #     chi2[i] = np.sum(((f2d-sh)**2/v2d)[u2d])
        # 
        # best_yshift = shifts[chi2 == chi2.min()][0]
        # sh = nd.shift(model_spec, (best_yshift, 0))
        # 
        # flux = np.cast[float](nd.shift(self.twod.im['SCI'].data, (-best_yshift, 0))-self.twod.im['CONTAM'].data).flatten()
        # coeffs = utils_c.run_nmf(flux[use], var[use], template_spec[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)
        # model_spec2 = np.dot(coeffs, template_spec).reshape(self.twod.im['SCI'].data.shape)
        
        
    def fit_at_z(self, z, nmf_toler=1.e-5, te_scale = 0.5):
        
        #### Template error
        te_func = self.te.interpolate(lci, z)*te_scale
        phot_var = self.efobs_fnu**2+(self.fobs_fnu*te_func)**2
        
        template_phot = self.phot.interpolate_photometry(z).reshape((self.phot.NFILT, self.phot.NTEMP)).T #[:-3,:]
        
        #tfull = self.phot.interpolate_photometry(zgrid)
        #i = 148
        #template_phot = tfull[:,:,i].reshape((self.phot.NFILT, self.phot.NTEMP)).T
        
        #### Fit the template to photometry
        amatrix = utils_c.prepare_nmf_amatrix(var[use], template_phot[:,use])
        coeffs = utils_c.run_nmf(fobs_fnu[use], var[use], template_phot[:,use], amatrix, toler=nmf_toler, MAXITER=1e6)
        model_phot = np.dot(coeffs, template_phot)
        igmz, igm_factor = self.phot.get_IGM(z, matrix=False)
        model_phot_full = np.dot(coeffs, self.phot.temp_seds.T)*igm_factor
        
        chi2 = np.sum((model_phot[use]-fobs_fnu[use])**2/var[use])
        chi2_obs = np.sum((model_phot[use]-obs_sed[use])**2/var[use])
        resid =  (model_phot[use]-fobs_fnu[use])/np.sqrt(var[use])
        #print resid
        return chi2
        
        #plt.errorbar(lci, fobs_fnu/(lci/5500.)**2, efobs_fnu/(lci/5500.)**2, color='blue')
        # plt.errorbar(lci, model_phot-fobs_fnu, efobs_fnu, color='red', marker='o', alpha=0.5)
        # plt.errorbar(lci, obs_sed*(lci/5500.)**2-fobs_fnu, efobs_fnu, color='green', marker='o', alpha=0.5)
        plt.errorbar(lci, fobs_fnu, efobs_fnu, color='blue', marker='o', alpha=0.5)
        plt.scatter(lci, obs_sed*(lci/5500.)**2, color='green', marker='o', alpha=0.5)
        plt.scatter(lci, model_phot, color='red', marker='o', alpha=0.5)
        

    
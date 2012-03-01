"""
Fit templates to the interlaced grism spectra to measure redshifts, and....

"""
import os

import pyfits

##### First option is STSCI_PYTHON v2.12 / Python2.7.  Second is older Python2.5.4
try:
    import stsci.convolve as convolve
except:
    import convolve
    
import numpy as np
import scipy.stats as stats

import matplotlib.pyplot as plt

# import emcee

import unicorn
import threedhst
import threedhst.catIO as catIO

import unicorn.utils_c as utils_c
import threedhst.eazyPy as eazy

def go_bright(skip_completed=True):
    """
    Fit spectra of bright galaxies in the test UDF pointing
    """
    import unicorn.interlace_fit
    
    os.chdir('/research/HST/GRISM/3DHST/UDF/PREP_FLT')
    model = unicorn.reduce.GrismModel('GOODS-S-34')
    os.chdir('mcz')
    
    os.chdir('/research/HST/GRISM/3DHST/UDS/PREP_FLT')

    os.chdir(unicorn.GRISM_HOME+'/GOODS-S/PREP_FLT')
    
    pointings = glob.glob('*_model.fits')
    skip_completed=True
    
    for point in pointings:
        pointing = point.split('_model')[0]
        #pointing = 'UDS-18'
        model = unicorn.reduce.GrismModel(pointing)
        bright = model.cat.mag < 23.5
        ids = model.cat.id[bright]
        for id in ids:     
            root='%s_%05d' %(pointing, id)
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #gris = unicorn.interlace_fit.GrismSpectrumFit(root='../GOODS-S-34_%05d' %(id))
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            if gris.dr > 1:
                continue
            #
            #gris.fit_in_steps(dzfirst=0.01, dzsecond=0.0002)
            gris.make_figure()
        
class GrismSpectrumFit():
    """
    Functions for fitting (redshifts for now) the interlaced grism spectra
    """
    def __init__(self, root='../GOODS-S-34_00280', FIGURE_FORMAT='png', verbose=True):
        """
        Read the 1D/2D spectra and get the photometric constraints
        necessary for the spectrum fits.
        """
        self.FIGURE_FORMAT=FIGURE_FORMAT
        
        #### Get the 1D/2D spectra
        self.twod = unicorn.reduce.Interlace2D(root+'.2D.fits', PNG=False)
        self.oned = unicorn.reduce.Interlace1D(root+'.1D.fits', PNG=False)
        self.grism_id = os.path.basename(self.twod.file.split('.2D.fits')[0])

        #### Get the photometric match.
        ra, dec = self.oned.header['RA'], self.oned.header['DEC']
        self.get_photometric_constraints(ra=ra, dec=dec, verbose=verbose)

        #### If match distance > 1", prior is flat
        if self.dr > 1:
            self.phot_lnprob *= 0

        #### Initialize the continuum and emission line templates    
        self.line_free_template()
        self.linex, self.liney = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
        
        #### Try to read the previous p(z) from the pickle file
        if self.read_pickle():
            if verbose:
                print 'Read p(z) from %s.zfit.pkl' %(self.grism_id)
                self.get_best_fit()
            
    def fit_zgrid(self, zrange = (0.01,6), dz = 0.02, get_model_at_z=None, verbose=True):
        """
        Take the best-fit template from the photo-z fit and just march through in redshift

        The redshift grid is set with the input parameters zrange and dz.  

        If `get_model_at_z` is set, just compute and return the best-fit 2D spectrum at the
        specified redshift.

        """
        #### Initialize redshift grid, steps are dz*(1+z)
        zgrid = np.exp(np.arange(np.log(1+zrange[0]), np.log(1+zrange[1]), dz))-1

        if get_model_at_z is not None:
            get_model = True
            zgrid = np.array([get_model_at_z])
        else:
            get_model = False

        NZ = len(zgrid)
        spec_lnprob = zgrid*0

        #### Observed flux and variance images of the 2D spectrum
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        var[var == 0] = 1.e6
        var += (0.1*self.twod.im['CONTAM'].data**2).flatten()

        flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        use = np.isfinite(flux)

        show = False

        #### Templates to allow blue/red tilt to the spectrum with respect to the EAZY templates.
        #### These are simple lines normalized at lam=1.4e4.  That it works appears a bit magical, but
        #### it does seem to.  Only need to compute them once since they don't change with the 
        #### redshift grid.
        tilt_blue = -(self.templam_nolines*(1+0.5)-1.4e4)/4000.+1
        self.twod.compute_model(self.templam_nolines*(1+0.5), tilt_blue)
        tilt_blue_model = self.twod.model*1.
        
        tilt_red = (self.templam_nolines*(1+0.5)-1.4e4)/4000.+1
        self.twod.compute_model(self.templam_nolines*(1+0.5), tilt_red)
        tilt_red_model = self.twod.model*1.

        #### Loop through redshift grid
        for i in range(NZ):
            if verbose:
                print unicorn.noNewLine+'z: %.2f' %(zgrid[i])
            
            z_test = zgrid[i]

            #### Generate the 2D spectrum model for the continuum and emission line templates
            if self.best_fit_nolines.sum() == 0.:
                self.best_fit_nolines = self.best_fit_nolines*0.+1
            self.twod.compute_model(self.templam_nolines*(1+z_test), self.best_fit_nolines); continuum_model = self.twod.model*1.
                
            self.twod.compute_model(self.linex*(1+z_test), self.liney); line_model = self.twod.model*1.

            #### Initialize if on the first grid point
            if i == 0:
                templates = np.ones((4,line_model.size))
                templates[2,:] = tilt_blue_model.flatten()
                templates[3,:] = tilt_red_model.flatten()

            #### Put in the spectral templates
            templates[0,:] = continuum_model.flatten()
            templates[1,:] = line_model.flatten()

            #### Compute the non-negative normalizations of the template components
            amatrix = utils_c.prepare_nmf_amatrix(var[use], templates[:,use])
            coeffs = utils_c.run_nmf(flux[use], var[use], templates[:,use], amatrix, toler=1.e-5)
            flux_fit = np.dot(coeffs.reshape((1,-1)), templates) #.reshape(line_model.shape)

            #### If `get_model_at_z`, return the model at that redshift
            if (get_model):
                flux_model = flux_fit.reshape(line_model.shape)
                oned_wave, model_oned = self.twod.optimal_extract(flux_model)
                return flux_model, oned_wave, model_oned

            #### Log likelihood at redshift zgrid[i] (-0.5*chi2)
            spec_lnprob[i] = -0.5*np.sum((flux-flux_fit)**2/var)

        #### Done!
        return zgrid, spec_lnprob
    
    def fit_in_steps(self, dzfirst=0.02, dzsecond=0.0001, save=True, make_plot=True):
        """
        Do two fit iterations, the first on a coarse redshift grid over the full z=(0.01,6)
        and the second refined grid around the peak found in the first iteration
        
        If `save` is set, save the results to the ascii and pickle files (.dat, .pkl).
        
        If `make_plot` is set, make the diagnostic plot.
        """
        zout = self.zout

        #################
        #### First iteration, fit full range at low dz resolution
        #################
        zgrid0, spec_lnprob0 = self.fit_zgrid(zrange=(0.01,6), dz=dzfirst)

        #### Interpolate the photometric p(z) to apply it as a prior
        phot_int0 = np.interp(zgrid0, self.phot_zgrid,  self.phot_lnprob)
        full_prob0 = phot_int0 + spec_lnprob0 - spec_lnprob0.max()
        full_prob0 -= full_prob0.max()  ### normalize p_max = 1

        z_max = zgrid0[full_prob0 == full_prob0.max()][0]
        zsub = full_prob0 > np.log(1.e-10)

        #### Define the range and step size for the next iteration
        width = 0.02*(1+z_max)
        zrange, dz = (z_max-width, z_max+width), dzsecond

        print '\nStep 1: z_max = %.3f (%.3f,%.3f)\n' %(z_max, zrange[0], zrange[1])

        #################
        #### Second iteration, small dz steps around the probability peak
        #################
        zgrid1, spec_lnprob1 = self.fit_zgrid(zrange=zrange, dz=dz)

        #### Interpolate the photometric p(z) to apply it as a prior
        phot_int1 = np.interp(zgrid1, self.phot_zgrid,  self.phot_lnprob)
        full_prob1 = phot_int1 + spec_lnprob1 - spec_lnprob1.max()
        full_prob1 -= full_prob1.max()

        self.z_max_spec = zgrid1[full_prob1 == full_prob1.max()][0]
        
        #### Save results and make the diagnostic plot
        self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1 = zgrid0, full_prob0, zgrid1, full_prob1
        
        ### Get the best 2D/1D fit
        self.get_best_fit()
        
        if save:
            ### Write the ascii file
            self.save_results()
            
            ### Save p(z) to a pickle
            self.write_pickle()    
        
        if make_plot:
            ### Make diagnostic figure
            self.make_figure()
    
    def save_results(self, verbose=True):
        """
        Print the results of the fit to a file, like 
        GOODS-S-34_00622.zfit.dat.

        """
        zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
        
        z_max_spec = zgrid1[full_prob1 == full_prob1.max()][0]
        z_peak_spec = np.trapz(zgrid1*np.exp(full_prob1), zgrid1)/np.trapz(np.exp(full_prob1), zgrid1)

        fp = open(self.grism_id+'.zfit.dat','w')
        fp.write('#  spec_id   phot_id   dr   z_spec  z_peak_phot  z_max_spec z_peak_spec\n')
        fp.write('#  Phot: %s\n' %(self.cat.filename))
        file_string = '%s %d  %.3f  %.5f  %.5f  %.5f  %.5f\n' %(self.grism_id, self.cat.id[self.ix], self.dr, self.zout.z_spec[self.ix], self.zout.z_peak[self.ix], z_max_spec, z_peak_spec)
        fp.write(file_string)
        fp.close()

        if verbose:
            print file_string
    
    def write_pickle(self):
        """

        Save the redshift/probability grid to a pickle file, like
        GOODS-S-34_00622.zfit.pkl

        """
        import pickle

        fp = open(self.grism_id+'.zfit.pkl','wb')
        pickle.dump(self.zgrid0, fp)
        pickle.dump(self.full_prob0, fp)
        pickle.dump(self.zgrid1, fp)
        pickle.dump(self.full_prob1, fp)
        fp.close()
    
    def read_pickle(self):
        """

        Save the redshift/probability grid to a pickle file, like
        GOODS-S-34_00622.zfit.pkl

        """
        import pickle
        
        if not os.path.exists(self.grism_id+'.zfit.pkl'):
            return False
            
        fp = open(self.grism_id+'.zfit.pkl','rb')
        self.zgrid0 = pickle.load(fp)
        self.full_prob0 = pickle.load(fp)
        self.zgrid1 = pickle.load(fp)
        self.full_prob1 = pickle.load(fp)
        fp.close()
        
        self.z_max_spec = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        return True
        
    def get_best_fit(self):
        """
        Get the 2D fit at the redshift of peak probability
        
        """
        z_max_spec = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]

        self.flux_model, self.oned_wave, self.model_oned = self.fit_zgrid(get_model_at_z=z_max_spec, verbose=False)
        
    def make_figure(self):
        """
        Make a plot showing the fit in e/s and f_lambda units, along with p(z) and the
        photometric SED

        """
        zout = self.zout
        
        zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
                
        #### Initialize the figure
        fig = unicorn.catalogs.plot_init(xs=10,aspect=1./3.8, left=0.1, right=0.02, bottom=0.09, top=0.08, NO_GUI=True)

        #### Spectrum in e/s
        ax = fig.add_subplot(141)
        wuse = (self.oned.data.wave > 1.15e4) & (self.oned.data.wave < 1.6e4)
        ax.plot(self.oned.data.wave/1.e4, self.oned.data.flux, color='black', alpha=0.1)
        ax.plot(self.oned.data.wave/1.e4, self.oned.data.flux-self.oned.data.contam, color='black')
        ax.plot(self.oned_wave/1.e4, self.model_oned, color='red', alpha=0.5, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'e$^-$ / s')
        if wuse.sum() > 5:
            ymax = self.oned.data.flux[wuse].max(); 
        else:
            ymax = self.oned.data.flux.max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax) ; ax.set_xlim(1.0,1.73)

        #### Spectrum in f_lambda
        ax = fig.add_subplot(142)
        ax.plot(self.oned.data.wave/1.e4, self.oned.data.flux/self.oned.data.sensitivity, color='black', alpha=0.1)
        ax.plot(self.oned.data.wave/1.e4, (self.oned.data.flux-self.oned.data.contam)/self.oned.data.sensitivity, color='black')
        ax.plot(self.oned_wave/1.e4, self.model_oned/self.oned.data.sensitivity, color='red', alpha=0.5, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'$f_\lambda$')
        if wuse.sum() > 5:
            ymax = (self.oned.data.flux/self.oned.data.sensitivity)[wuse].max()
        else:
            ymax = (self.oned.data.flux/self.oned.data.sensitivity).max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax) ; ax.set_xlim(1.0,1.73)

        #### p(z)
        ax = fig.add_subplot(143)
        ax.plot(self.phot_zgrid, np.exp(self.phot_lnprob), color='green')
        ax.plot(zgrid0, np.exp(full_prob0), color='blue', alpha=0.4)
        ax.plot(zgrid1, np.exp(full_prob1), color='blue')
        if zout.z_spec[self.ix] > 0:
            ax.plot(zout.z_spec[self.ix]*np.array([1,1]), [0,1], color='red')

        if self.dr < 1:
            ax.set_xlim(zout.l99[self.ix], zout.u99[self.ix])

        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'$p(z)$')
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune='both'))

        #### Make title text
        if zout.z_spec[self.ix] > 0:
            deltaz = '$\Delta z$ = %.4f' %((zout.z_spec[self.ix]-self.z_max_spec)/(1+zout.z_spec[self.ix]))
        else:
            deltaz = ''

        ax.text(-0.05, 1.1, r'%s  %d  $H_{140}=$%.2f $z_\mathrm{spec}$=%.3f  $z_\mathrm{phot}$=%.3f  $z_\mathrm{gris}$=%.3f  %s' %(self.grism_id, self.zout.id[self.ix], self.twod.im[0].header['MAG'], zout.z_spec[self.ix], zout.z_peak[self.ix], self.z_max_spec, deltaz), transform=ax.transAxes, horizontalalignment='center')

        ####  Show the (photometric) SED with the spectrum overplotted
        ax = fig.add_subplot(144)

        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit

        temp_sed_int = np.interp(self.oned.data.wave, lambdaz, temp_sed)
        keep = (self.oned.data.wave > 1.2e4) & (self.oned.data.wave < 1.5e4)
        flux_spec = (self.oned.data.flux-self.oned.data.contam)/self.oned.data.sensitivity
        anorm = np.sum(temp_sed_int[keep]*flux_spec[keep])/np.sum(flux_spec[keep]**2)

        ax.plot(lambdaz, temp_sed, color='blue', alpha=0.5)
        ax.plot(lci, fobs, color='black', marker='o', ms=7, alpha=0.7, linestyle='None')
        ax.plot(self.oned.data.wave, flux_spec*anorm, color='red', alpha=0.5)
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$f_\lambda$')
        
        good = efobs > 0
        if good.sum() > 0:
            ymax = fobs[good].max()
        else:
            ymax = fobs.max()
            
        ax.semilogx(); ax.set_xlim(3000.,8.e4); ax.set_ylim(-0.05*ymax, 1.1*ymax)

        #### Save the result to a file
        unicorn.catalogs.savefig(fig, self.grism_id+'.zfit.%s' %(self.FIGURE_FORMAT))

        #z_peaks = zgrid[1:-1][(full_prob[1:-1] > np.log(0.05)) & (np.diff(full_prob,2) < 0)]
        #zrange, dz = (z_peaks[0]-0.02*(1+z_peaks[0]), z_peaks[0]+0.02*(1+z_peaks[0])), 0.0002
    
    def get_photometric_constraints(self, ra=0., dec=0., verbose=True):
        """ 
        Read the overlapping photometric catalog and retrieve the photometry and EAZY fit.

        """
        # GRISM_PATH=unicorn.GRISM_HOME+'GOODS-S/'
        # CAT_PATH = GRISM_PATH+'FIREWORKS/'
        # if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
        #     CAT_PATH = '/3DHST/Ancillary/GOODS-S/FIREWORKS/FAST/'
        # #
        # CAT_FILE = CAT_PATH+'fireworks.cat'
        # ZOUT_FILE = CAT_PATH+'fireworks.zout'
        # FOUT_FILE = CAT_PATH + 'fireworks.fout'
        # KTOT_COL = 'Ks_totf'

        #### Read in EAZY files
        cat, zout, fout = unicorn.analysis.read_catalogs(root=self.grism_id)
        #cat.kmag = 23.86-2.5*np.log10(cat.field(KTOT_COL))
        CAT_PATH = os.path.dirname(cat.filename)
        
        root = os.path.basename(zout.filename).split('.zout')[0]
        ZOUT_PATH = os.path.dirname(zout.filename)
        
        tempfilt, self.eazy_coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE = root,                                                 OUTPUT_DIRECTORY=ZOUT_PATH, CACHE_FILE = 'Same')

        # eazyParam = eazy.EazyParam(ZOUT_PATH+'/'+root+'.param')

        # id = 6187
        # ix = np.arange(cat.id.shape[0])[cat.id == id][0]

        dr = np.sqrt((cat.ra-ra)**2*np.cos(dec/360.*2*np.pi)**2+(cat.dec-dec)**2)*3600.
        ix = np.arange(cat.id.shape[0])[dr == dr.min()][0]

        self.eazy_fit = eazy.getEazySED(ix, MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=ZOUT_PATH, CACHE_FILE = 'Same')

        self.ix = ix
        self.zout = zout
        self.cat = cat
        self.z_peak = zout.z_peak[ix]
        self.dz = (zout.u68[ix]-zout.l68[ix])/2.

        self.dr = dr.min()
        
        if verbose:
            print 'Match in %s\n #%d  dr=%.2f"  z_spec=%7.3f  z_peak=%7.3f' %(cat.filename, cat.id[ix], self.dr, zout.z_spec[ix], self.z_peak)
            
        self.best_fit = np.dot(temp_sed['temp_seds'], self.eazy_coeffs['coeffs'][:,ix])
        self.templam = temp_sed['templam']*1.

        # self.phot_fnu = tempfilt['fnu'][:,ix]
        # self.phot_efnu = tempfilt['efnu'][:,ix]
        # self.phot_lc = tempfilt['lc']
        # 
        # self.tempfilt = tempfilt
        # self.tempgrid = tempfilt['tempfilt']

        self.phot_zgrid = tempfilt['zgrid']
        self.phot_lnprob = -0.5*pz['chi2fit'][:,ix]
        self.phot_lnprob -= self.phot_lnprob.max()

        # self.temp_err_x, self.temp_err_y = np.loadtxt(os.getenv('TEMPLATE_DIR')+'/'+eazyParam.params['TEMP_ERR_FILE'], unpack=True)
        # self.temp_err_y *= eazyParam.params['TEMP_ERR_A2']

    def line_free_template(self):
        """
        Assuming that the v1.1 line templates were used in the fit, generate the best-fit template
        from the version *without* emission lines.
        """

        #### Read in the line-free templates and scale them with "tnorm" to match those stored
        #### in the temp_sed file
        nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed1_nolines.dat', unpack=True)
        noline_temps = np.zeros((nlx.shape[0],7))
        noline_temps[:,0] = nly/self.eazy_coeffs['tnorm'][0]
        for i in range(2,7):
            nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed%d_nolines.dat' %(i), unpack=True)
            noline_temps[:,i-1] = nly/self.eazy_coeffs['tnorm'][i-1]

        #### The last v1.1 template was a BC03 model without lines, so just use it directly    
        i = 7
        lx, ly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.1_lines/eazy_v1.1_sed%d.dat' %(i), unpack=True)
        noline_temps[:,6] = np.interp(nlx, lx, ly)/self.eazy_coeffs['tnorm'][6]

        self.best_fit_nolines = np.dot(noline_temps, self.eazy_coeffs['coeffs'][:,self.ix])
        self.templam_nolines = nlx
    
    def flag_contamination(self, FEXCESS=2.5):
        # import threedhst.dq
        # 
        # root='UDS-18_00180'        
        # gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
        # 
        # gris.twod.compute_model()
        # ds9 = threedhst.dq.myDS9()
        # ds9.frame(1)
        # ds9.v(gris.twod.im['SCI'].data, vmin=-0.05, vmax=0.2)
        # ds9.frame(2)
        # ds9.v(gris.twod.model, vmin=-0.05, vmax=0.2)
        # ds9.frame(3)
        # ds9.v(gris.twod.im['SCI'].data-gris.twod.model, vmin=-0.05, vmax=0.2)
        # 
        # flag = (gris.twod.im['SCI'].data-(gris.twod.im['CONTAM'].data+gris.twod.model)) / gris.twod.im['WHT'].data
        
        #### Find wavelengths counting in from the edges where the observed extracted
        #### spectrum is a factor FEXCESS greater than the model (with contamination).
        #### This is intended to find cases at the left edge of the frame where objects
        #### fall out of the direct image but still have a 1st order spectrum in the grism image.
        ####
        #### Note can't apply blindly because emission lines will also satisfy such a cut.
        self.twod.compute_model()
        wave_model, flux_model = self.twod.optimal_extract(self.twod.model+self.twod.im['CONTAM'].data)
        wave_obs, flux_obs = self.twod.optimal_extract(self.twod.im['SCI'].data)
        # FEXCESS = 2.5
        test = flux_obs/flux_model < FEXCESS
        if test.sum() > 0:
            first = wave_model[test][0]
            self.twod.im['CONTAM'].data[:,self.twod.im['WAVE'].data < first] += (self.twod.im['SCI'].data-self.twod.model)[:,self.twod.im['WAVE'].data < first]
            print 'Unflagged contamination: lam < %.2f' %(first)

        
        
def go_MCMC_fit():
    """
    Try a MCMC fit, TESTING 
    """
    #### EMCEE
    import emcee
    init = np.array([unicorn.interlace_fit.z_peak, 0])
    step_sig = np.array([unicorn.interlace_fit.dz, 0.03])
    obj_fun = unicorn.interlace_fit.objective_twod
    
    ndim, nwalkers = len(init), 100
    p0 = [(init+np.random.normal(size=2)*step_sig) for i in xrange(nwalkers)]
    
    NTHREADS, NSTEP = 1, 100
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, threads=NTHREADS)
    result = sampler.run_mcmc(p0, NSTEP)
    chain = sampler.flatchain
    
    z_fit = np.median(chain[:,0])
    slope = np.median(chain[:,1])
    best_model = unicorn.interlace_fit.objective_twod([z_fit, slope], return_model=True)
    
    ds9.view(twod.im['SCI'].data-twod.im['CONTAM'].data)
    
def objective_twod(params, return_model=False):
    """
    Objective function for the MCMC redshift fit
    
    params[0] = z
    params[1] = continuum slope
    
    """
    import unicorn.interlace_fit
    
    continuum = params[1]*(unicorn.interlace_fit.linex-1.4e4)/1000.+1
    unicorn.interlace_fit.twod.compute_model(unicorn.interlace_fit.linex*(1+params[0]), continuum); continuum_model = unicorn.interlace_fit.twod.model*1.
    unicorn.interlace_fit.twod.compute_model(unicorn.interlace_fit.linex*(1+params[0]), unicorn.interlace_fit.liney); line_model = unicorn.interlace_fit.twod.model*1.
    
    templates = np.ones((2,line_model.size))
    templates[0,:] = continuum_model.flatten()
    templates[1,:] = line_model.flatten()
    
    var = np.cast[float](twod.im['WHT'].data**2).flatten()
    var[var == 0] = 1.e6
    flux = np.cast[float](twod.im['SCI'].data-twod.im['CONTAM'].data).flatten()
    
    amatrix = utils_c.prepare_nmf_amatrix(var, templates)
    coeffs = utils_c.run_nmf(flux, var, templates, amatrix, toler=1.e-7)
    flux_fit = np.dot(coeffs.reshape((1,-1)), templates) #.reshape(line_model.shape)
    if return_model:
        flux_model = flux_fit.reshape(line_model.shape)
        return flux_model
        
    lnprob = -0.5*np.sum((flux-flux_fit)**2/var)
    prior_z = np.interp(params[0], unicorn.interlace_fit.phot_zgrid, unicorn.interlace_fit.phot_lnprob)
    
    #print "%.3f %.5f %.3f %.3f" %(params[0], params[1], lnprob, prior_z)
    return lnprob+prior_z

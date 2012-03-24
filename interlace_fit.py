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
        bright = model.cat.mag < 23
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
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
            #gris.make_figure()
        
class GrismSpectrumFit():
    """
    Functions for fitting (redshifts for now) the interlaced grism spectra
    """
    def __init__(self, root='../GOODS-S-34_00280', FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.55):
        """
        Read the 1D/2D spectra and get the photometric constraints
        necessary for the spectrum fits.
        """
        self.FIGURE_FORMAT=FIGURE_FORMAT
        self.root = os.path.basename(root)

        plt.rcParams['lines.marker'] = ''
        plt.rcParams['lines.linestyle'] = '-'        
        
        self.use_lines = None
        
        #### Get the 1D/2D spectra
        self.twod = unicorn.reduce.Interlace2D(root+'.2D.fits', PNG=False)
        self.status = True
        if self.twod.im['SCI'].data.max() <= 0:
            if verbose:
                threedhst.showMessage('%s: \nNo non-zero pixels in the 2D spectrum.' %(root), warn=True)
            self.status = False
            return None
            
        self.oned = unicorn.reduce.Interlace1D(root+'.1D.fits', PNG=False)
        if self.oned.data.flux.max() <= 0:
            print '%s: No valid pixels in 1D spectrum.' %(root)
            self.status = False
            return None
            
        #### Convert to 10**-17 ergs / s / cm**2 / A
        #self.oned.data.sensitivity *= np.diff(self.oned.data.wave)[0]
        
        self.grism_id = os.path.basename(self.twod.file.split('.2D.fits')[0])

        #### Get the photometric match.
        ra, dec = self.oned.header['RA'], self.oned.header['DEC']
        self.get_photometric_constraints(ra=ra, dec=dec, verbose=verbose)

        #### If match distance > 1", prior is flat
        if self.dr > 1:
            self.phot_lnprob *= 0
        else:
            #### Put less weight on spectra when z_phot < `lowz_thresh` (0.55) and there aren't 
            #### many features expected in the spectrum.
            zgt = self.phot_zgrid > lowz_thresh
            zprob = np.trapz(self.phot_linear[zgt], self.phot_zgrid[zgt])
            if zprob < 0.3:
                if verbose:
                    print '\n!! p(z < %.2f | phot) = %.3f, I\'ll decrease spectrum weight in fit.\n' %(lowz_thresh, 1-zprob)
                self.twod.im['WHT'].data *= 3.
                #self.twod.im['WHT'].data = np.maximum(self.twod.im['WHT'].data, self.twod.im['SCI'].data/5.)
                
        #### Initialize the continuum and emission line templates    
        self.line_free_template()
        self.linex, self.liney = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
        
        #### Try to read the previous p(z) from the pickle file
        if self.read_pickle():
            if verbose:
                print 'Read p(z) from %s.zfit.pkl' %(self.grism_id)
                self.get_best_fit()
    
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
        
        ### normalized probability
        self.phot_linear = np.exp(self.phot_lnprob)
        self.phot_linear /= np.trapz(self.phot_linear, self.phot_zgrid)
        
        # self.temp_err_x, self.temp_err_y = np.loadtxt(os.getenv('TEMPLATE_DIR')+'/'+eazyParam.params['TEMP_ERR_FILE'], unpack=True)
        # self.temp_err_y *= eazyParam.params['TEMP_ERR_A2']
                    
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
            
            ##### Probably no flux in the direct image
            if templates.max() == 0:
                self.status = False
                return False
                
            #### Compute the non-negative normalizations of the template components
            amatrix = utils_c.prepare_nmf_amatrix(var[use], templates[:,use])
            coeffs = utils_c.run_nmf(flux[use], var[use], templates[:,use], amatrix, toler=1.e-5)
            flux_fit = np.dot(coeffs.reshape((1,-1)), templates) #.reshape(line_model.shape)

            #### If `get_model_at_z`, return the model at that redshift
            if (get_model):
                ixc = np.array([0,2,3])
                self.cont_model = np.dot(coeffs[ixc].reshape((1,-1)), templates[ixc,:]).reshape(line_model.shape)
                self.slope_model = np.dot(coeffs[2:].reshape((1,-1)), templates[2:,:]).reshape(line_model.shape)
                self.line_model = (coeffs[1]*templates[1,:]).reshape(line_model.shape)
                self.flux_model = flux_fit.reshape(line_model.shape)
                self.oned_wave, self.model_1D = self.twod.optimal_extract(self.flux_model)
                oned_wave_x, self.cont_1D = self.twod.optimal_extract(self.cont_model)
                oned_wave_x, self.slope_1D = self.twod.optimal_extract(self.slope_model)
                oned_wave_x, self.line_1D = self.twod.optimal_extract(self.line_model)
                
                return True
                
                #return flux_model, cont_model, line_model, oned_wave, model_1D, cont_1D, line_1D, slope_1D

            #### Log likelihood at redshift zgrid[i] (-0.5*chi2)
            spec_lnprob[i] = -0.5*np.sum((flux-flux_fit)**2/var)

        #### Done!
        return zgrid, spec_lnprob
    
    def fit_in_steps(self, dzfirst=0.01, zrfirst = (0.01,6), dzsecond=0.0002, save=True, make_plot=True):
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
        result = self.fit_zgrid(zrange=zrfirst, dz=dzfirst)
        if result is False:
            self.status = False
            return False
        
        zgrid0, spec_lnprob0 = result
        
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
            self.twod_figure()
            self.make_fits()
            
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
        
        #self.flux_model, self.cont_model, self.line_model, self.oned_wave, self.model_1D, self.cont_1D, self.line_1D, self.slope_1D = self.fit_zgrid(get_model_at_z=z_max_spec, verbose=False)
        status = self.fit_zgrid(get_model_at_z=z_max_spec, verbose=False)
    
    def make_fits(self):
        """
        Make an output FITS file containing the continuum and line models, 1D and 2D
        """
        zout = self.zout
        zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
        
        prim = pyfits.PrimaryHDU()
        prim.header = self.twod.im[0].header.copy()
        hdus = [prim]
        
        twod_header = self.twod.im['SCI'].header.copy()
        
        #self.model_1D
        #self.oned.data.sensitivity
        #self.cont_model
        #self.line_model
        
        hdus.append(pyfits.ImageHDU(data=self.cont_model, header=twod_header, name='CONT2D'))
        hdus.append(pyfits.ImageHDU(data=self.line_model, header=twod_header, name='LINE2D'))
        hdus.append(pyfits.ImageHDU(data=self.cont_1D, name='CONT1D'))
        hdus.append(pyfits.ImageHDU(data=self.line_1D, name='LINE1D'))
        
        pyfits.HDUList(hdus).writeto(self.grism_id+'.zfit.fits', clobber=True)
        
    def make_figure(self):
        """
        Make a plot showing the fit in e/s and f_lambda units, along with p(z) and the
        photometric SED
        
        xxx To do: store continuum and 1D fits, fit emission lines fixed at z_grism
        """
        zout = self.zout
        
        zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
                
        #### Initialize the figure
        fig = unicorn.catalogs.plot_init(xs=10,aspect=1./3.8, left=0.1, right=0.02, bottom=0.09, top=0.08, NO_GUI=True)

        show = self.oned.data.flux != 0.0
        #### Spectrum in e/s
        ax = fig.add_subplot(141)
        wuse = (self.oned.data.wave > 1.15e4) & (self.oned.data.wave < 1.6e4)

        ax.plot(self.oned.data.wave[show]/1.e4, self.oned.data.flux[show], color='black', alpha=0.1)
        ax.plot(self.oned.data.wave[show]/1.e4, self.oned.data.flux[show]-self.oned.data.contam[show], color='black')
        ax.plot(self.oned_wave[show]/1.e4, self.model_1D[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave/1.e4, self.model_1D, color='red', alpha=0.08, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'e$^-$ / s')
        if wuse.sum() > 5:
            ymax = self.oned.data.flux[wuse].max(); 
        else:
            ymax = self.oned.data.flux.max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax) ; ax.set_xlim(1.0,1.73)

        #### Spectrum in f_lambda
        ax = fig.add_subplot(142)
        ax.plot(self.oned.data.wave[show]/1.e4, self.oned.data.flux[show]/self.oned.data.sensitivity[show], color='black', alpha=0.1)
        
        show_flux = (self.oned.data.flux[show]-self.oned.data.contam[show])/self.oned.data.sensitivity[show]
        show_err = self.oned.data.error[show]/self.oned.data.sensitivity[show]
        ax.plot(self.oned.data.wave[show]/1.e4, show_flux, color='black')
        #ax.fill_between(self.oned.data.wave[show]/1.e4, show_flux+show_err, show_flux-show_err, color='0.5', alpha=0.2)
        
        ax.plot(self.oned_wave/1.e4, self.model_1D/self.oned.data.sensitivity, color='red', alpha=0.08, linewidth=2)
        ax.plot(self.oned_wave[show]/1.e4, self.model_1D[show]/self.oned.data.sensitivity[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave[show]/1.e4, self.slope_1D[show]/self.oned.data.sensitivity[show], color='orange', alpha=0.2, linewidth=1)
        
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
            deltaz = '$\Delta z$ = %.4f' %(-(zout.z_spec[self.ix]-self.z_max_spec)/(1+zout.z_spec[self.ix]))
        else:
            deltaz = ''

        ax.text(-0.05, 1.1, r'%s  %d  $H_{140}=$%.2f $z_\mathrm{spec}$=%.3f  $z_\mathrm{phot}$=%.3f  $z_\mathrm{gris}$=%.3f  %s' %(self.grism_id, self.zout.id[self.ix], self.twod.im[0].header['MAG'], zout.z_spec[self.ix], zout.z_peak[self.ix], self.z_max_spec, deltaz), transform=ax.transAxes, horizontalalignment='center')

        ####  Show the (photometric) SED with the spectrum overplotted
        ax = fig.add_subplot(144)

        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit

        temp_sed_int = np.interp(self.oned.data.wave, lambdaz, temp_sed)
        keep = (self.oned.data.wave > 1.2e4) & (self.oned.data.wave < 1.5e4)
        flux_spec = (self.oned.data.flux-self.oned.data.contam-self.slope_1D)/self.oned.data.sensitivity
        anorm = np.sum(temp_sed_int[keep]*flux_spec[keep])/np.sum(flux_spec[keep]**2)
        
        ### factor of 100 to convert from 1.e-17 to 1.e-19 flux units
        scale = 100
        
        ax.plot(lambdaz, temp_sed*scale, color='blue', alpha=0.5)
        ax.errorbar(lci, fobs*scale, efobs*scale, color='black', marker='o', ms=7, alpha=0.7, linestyle='None')
        ax.plot(self.oned.data.wave, flux_spec*anorm*100, color='red', alpha=0.5)
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$f_\lambda$')
        
        good = efobs > 0
        if good.sum() > 0:
            ymax = fobs[good].max()
        else:
            ymax = fobs.max()
        
        ymax *= scale
        
        ax.semilogx(); ax.set_xlim(3000.,8.e4); ax.set_ylim(-0.05*ymax, 1.1*ymax)
        #ax.set_xlim(0.8e4,2.e4)
        
        #### Save the result to a file
        unicorn.catalogs.savefig(fig, self.grism_id+'.zfit.%s' %(self.FIGURE_FORMAT))

        #z_peaks = zgrid[1:-1][(full_prob[1:-1] > np.log(0.05)) & (np.diff(full_prob,2) < 0)]
        #zrange, dz = (z_peaks[0]-0.02*(1+z_peaks[0]), z_peaks[0]+0.02*(1+z_peaks[0])), 0.0002
    
    def twod_figure(self, vmax=None):
        """
        Make a figure showing the raw twod-spectrum, plus contamination- and 
        continuum-subtracted versions.
        """
        
        sh = self.twod.im[4].data.shape
        top = 0.055
        bottom = 0.07
        left = 0.06
        aspect = 3.*sh[0]/sh[1]/(1-(top+bottom-left))
        
        wave = self.twod.im[8].data
        xint = [1.1,1.2,1.3,1.4,1.5,1.6]
        ax_int = np.interp(np.array(xint)*1.e4, wave, np.arange(wave.shape[0]))
        
        fig = unicorn.catalogs.plot_init(xs=5,aspect=aspect, left=left, right=0.02, bottom=bottom, top=top, NO_GUI=True)
        fig.subplots_adjust(hspace=0.001)
        
        if vmax==None:
            #values = self.twod.im[4].data.flatten()
            #vmax = values[np.argsort(values)][-10]
            vmax = self.flux_model.max()
            
        #### Raw spectrum
        ax = fig.add_subplot(311)
        ax.imshow(0-self.twod.im[4].data, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks(ax_int)
        ax.set_ylabel('Raw')
        #ax.text(0.95,1+0.1*aspect,self.grism_id, transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')
        ax.set_title(self.grism_id)
        
        #### Contam-subtracted spectrum
        ax = fig.add_subplot(312)
        ax.imshow(0-self.twod.im[4].data+self.twod.im[7].data, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks(ax_int)
        ax.set_ylabel('-Contam.')

        #### Continuum-subtracted spectrum
        ax = fig.add_subplot(313)
        ax.imshow(0-self.twod.im[4].data+self.twod.im[7].data+self.cont_model, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto')
        ax.set_yticklabels([])
        ax.set_xticklabels(xint)
        ax.set_xticks(ax_int)
        ax.set_ylabel('-Continuum')
        ax.set_xlabel(r'$\lambda\ (\mu\mathrm{m})$')
        
        unicorn.catalogs.savefig(fig, self.grism_id+'.zfit.2D.%s' %(self.FIGURE_FORMAT))
        
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
    
    def fit_free_emlines(self, ztry=None, verbose=True, NTHREADS=1, NWALKERS=100, NSTEP=100, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100):
        import emcee

        if ztry is None:
            ztry = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        if verbose:
            print 'Fit lines: z=%.4f' %(ztry)
            
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        var[var == 0] = 1.e6
        var += (0.1*self.twod.im['CONTAM'].data**2).flatten()

        flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        use = np.isfinite(flux)
        
        (use_lines, fancy), templates = self.get_emline_templates(ztry=ztry, line_width=line_width0)
        
        #### Initialize fit
        amatrix = utils_c.prepare_nmf_amatrix(var[use], templates[:,use])
        coeffs = utils_c.run_nmf(flux[use], var[use], templates[:,use], amatrix, toler=1.e-5)
        flux_fit = np.dot(coeffs.reshape((1,-1)), templates).reshape((self.twod.im['SCI'].data.shape))
        
        coeffs = np.maximum(coeffs, 1.e-5)
        #### MCMC fit with EMCEE sampler
        init = np.log(coeffs)
        step_sig = np.ones(len(coeffs))*np.log(1.05)
        
        #### Fit just emission lines, z fixed to z_grism
        obj_fun = unicorn.interlace_fit._objective_lineonly
        obj_args = [flux, var, templates]
        
        #### Try to fit the redshift, maybe will take forever
        if FIT_REDSHIFT | FIT_WIDTH:
            obj_fun = unicorn.interlace_fit._objective_z_line
            obj_args = [flux, var, self]
            ### Add redshift as parameter
            init = np.append(init, ztry)
            step_sig = np.append(step_sig, 0.003*FIT_REDSHIFT)
            ### Add (log) line width as parameter
            init = np.append(init, np.log(line_width0))
            step_sig = np.append(step_sig, 0.4*FIT_WIDTH)
        
        #### Run the Markov chain
        ndim = len(init)
        p0 = [(init+np.random.normal(size=len(init))*step_sig) for i in xrange(NWALKERS)]

        # NTHREADS, NSTEP = 1, 100
        if verbose:
            print 'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d)' %(NWALKERS, NSTEP)
            
        self.sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
        result = self.sampler.run_mcmc(p0, NSTEP)
        chain = np.exp(self.sampler.flatchain)
        
        if verbose:
            print unicorn.noNewLine+'emcee MCMC fit: (nwalkers x NSTEPS) = (%d x %d) ==> Done.\nMake figure.' %(NWALKERS, NSTEP)
        
        #### Store the results
        fig = unicorn.catalogs.plot_init(xs=5,aspect=1./1.65, left=0.15, right=0.02, bottom=0.09, top=0.02, NO_GUI=True)
        ax = fig.add_subplot(111)
        
        wok = (self.oned_wave > 1.e4)
        ax.plot(self.oned_wave[wok]/1.e4, self.oned.data.flux[wok]-self.oned.data.contam[wok], color='black', alpha=0.8)
        
        temps = []
        
        ntemp = templates.shape[0]
        for i in range(ntemp):
            t2d = templates[i,:].reshape(self.flux_model.shape)
            wave_model, flux_model = self.twod.optimal_extract(t2d)
            temps.append(flux_model)
        
        #### Line equivalent widths
        eqw = np.ones(ntemp-1)
        eqw_err = np.ones(ntemp-1)
        #halpha_eqw = -np.trapz((-halpha/continuum)[1:-1], temp_seds['templam'][1:-1]
        ok = temps[0] > 0
        for i in range(ntemp-1):
            eqw_post = -np.trapz(-temps[i+1][ok]/temps[0][ok], wave_model[ok])*chain[:,i+1]/chain[:,0]
            eqw[i] = np.median(eqw_post)
            eqw_err[i] = threedhst.utils.biweight(eqw_post)
        
        NSHOW = 100; ix_show = np.cast[int](np.random.rand(NSHOW)*NSTEP*NWALKERS)
        for i in ix_show:
            mc_model = temps[0]*chain[i,0]
            for j in range(1,ntemp):
                mc_model += temps[j]*chain[i,j]
            #
            aa = ax.plot(self.oned_wave[wok]/1.e4, mc_model[wok], color='red', alpha=0.05)
        
        ax.plot(self.oned_wave[wok]/1.e4, self.model_1D[wok], color='green', alpha=0.5, linewidth=2, zorder=1)
        
        fp = open(self.grism_id+'.linefit.dat','w')
        fp.write('# line  flux  error  EQW_obs EQW_obs_err \n# z=%.5f\n# flux: 10**-17 ergs / s / cm**2\n' %(ztry))
        
        print use_lines
        for i, line in enumerate(use_lines):
            #hist = plt.hist(chain[:,i+1], bins=50)
            stats = threedhst.utils.biweight(chain[:,i+1], both=True)
            #median = np.median(chain[:,i+1])
            #range = np.percentile(chain[:,i+1], [15.8655,100-15.8655])
            fp.write('%4s  %6.2f  %.2f  %6.2f %6.2f\n' %(line, stats[0], stats[1], eqw[i], eqw_err[i]))
            ax.fill_between([0.03,0.53],np.array([0.95,0.95])-0.07*i, np.array([0.88, 0.88])-0.07*i, color='white', alpha=0.8, transform=ax.transAxes, zorder=19)
            
            ax.text(0.05, 0.95-0.07*i, '%4s  %6.1f$\pm$%.1f  %4.1f  %6.1f\n' %(fancy[line], stats[0], stats[1], stats[0]/stats[1], eqw[i]), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, zorder=20)
            
        ax.fill_between([0.61,0.96],np.array([0.95,0.95]), np.array([0.88, 0.88]), color='white', alpha=0.8, transform=ax.transAxes, zorder=19)
        ax.text(0.95, 0.95, self.grism_id, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, zorder=20)
        ax.fill_between([0.8,0.96],np.array([0.85,0.85]), np.array([0.78, 0.78]), color='white', alpha=0.8, transform=ax.transAxes, zorder=19)
        ax.text(0.95, 0.85, r'$z=%.4f$' %(ztry), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, zorder=20)
        
        ax.set_xlim(1.0,1.73)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'Flux (e - / s)')
        ymax = self.model_1D.max()
        ax.set_ylim(-0.1*ymax, 1.5*ymax)
        
        unicorn.catalogs.savefig(fig, self.grism_id+'.linefit.'+self.FIGURE_FORMAT)
        
        fp.close()
        
        return True
        ###
        files=glob.glob('GOODS-S-*2D.fits')
        for file in files:
            root = file.replace('.2D.fits','')
            if os.path.exists(root+'.linefit.png'):
                continue
            self = unicorn.interlace_fit.GrismSpectrumFit(root)
            if self.status:
                self.fit_free_emlines()
    
    def get_emline_templates(self, ztry=0., use_determined_lines=True, line_width=100):
        
        ### Continuum at z=ztry
        status = self.fit_zgrid(get_model_at_z=ztry, verbose=False)
        
        line_wavelengths = {} ; line_ratios = {}
        line_wavelengths['Ha'] = [6564.61]; line_ratios['Ha'] = [1.]
        line_wavelengths['Hb'] = [4862.68]; line_ratios['Hb'] = [1.]
        line_wavelengths['Hg'] = [4341.68]; line_ratios['Hg'] = [1.]
        line_wavelengths['OIII'] = [5008.240, 4960.295]; line_ratios['OIII'] = [2.98, 1]
        line_wavelengths['OII'] = [3729.875]; line_ratios['OII'] = [1]
        line_wavelengths['SII'] = [6718.29, 6732.67]; line_ratios['SII'] = [1, 1]
        line_wavelengths['SIII'] = [9068.6, 9530.6]; line_ratios['SIII'] = [1, 2.44]
        line_wavelengths['HeI'] = [3889.0]; line_ratios['HeI'] = [1.]
        line_wavelengths['MgII'] = [2799.117]; line_ratios['MgII'] = [1.]
        line_wavelengths['CIV'] = [1549.480]; line_ratios['CIV'] = [1.]
        line_wavelengths['Lya'] = [1215.4]; line_ratios['Lya'] = [1.]
        
        fancy = {}
        fancy['Ha'] = r'H$\alpha$'
        fancy['Hb'] = r'H$\beta$'
        fancy['Hg'] = r'H$\gamma$'
        fancy['OIII'] =  'O III'
        fancy['OII']  =  'O II' 
        fancy['SII']  =  'S II' 
        fancy['SIII'] =  'S III'
        fancy['HeI']  =  'He I' 
        fancy['MgII'] =  'Mg II'
        fancy['CIV']  =  'C IV' 
        fancy['Lya'] = r'Ly$\alpha$'
        
        if use_determined_lines & (self.use_lines is not None):
            use_lines = self.use_lines
        else:
            use_lines = []
            for line in line_wavelengths.keys():
                lam = line_wavelengths[line][0]
                if (lam*(1+ztry) > self.oned_wave[self.oned_wave > 0.9e4].min()) & (lam*(1+ztry) < self.oned_wave.max()):
                    use_lines.append(line)
        
        templates = np.ones((1+len(use_lines), self.twod.im['SCI'].data.size))
        templates[0,:] = self.cont_model.flatten()
        
        #line_width = 100 # km/s
        xline = np.exp(np.arange(np.log(900), np.log(1.e4), line_width/3.e5/3.))
        sens_int = {}
        
        has_line_flux = np.ones(1+len(use_lines)) > 0
        
        for i,line in enumerate(use_lines):
            yline = xline*0.
            nmult = len(line_wavelengths[line])
            for j in range(nmult):
                lam = line_wavelengths[line][j]
                norm = line_ratios[line][j]
                dlam = lam*line_width/3.e5
                yline += 1./np.sqrt(np.pi*2*dlam**2)*np.exp(-(xline-lam)**2/2/dlam**2)*norm
            # 
            sens_int[line] = np.interp(lam, self.oned_wave, self.oned.data.sensitivity)
            #print lam
            #
            self.twod.compute_model(xline*(1+ztry), yline)
            self.twod.model = np.maximum(self.twod.model, 0.)
            #
            #### templates normalized to  total flux 1.e-17 erg / s / cm2 / A
            wave_model, flux_model = self.twod.optimal_extract(self.twod.model)
            cgs_model = flux_model/self.oned.data.sensitivity
            bad = ~np.isfinite(cgs_model)
            if bad.sum() > 0:
                cgs_model[bad] = 0.
            #
            if cgs_model.max() == 0.:
                has_line_flux[i+1] = False
            #
            norm_flux = np.trapz(cgs_model, self.oned_wave)
            templates[i+1,:] = self.twod.model.flatten() / norm_flux
        
        #### Some lines might not actually have any flux, so remove them
        if ~use_determined_lines | (self.use_lines is None):
            for i, line in enumerate(use_lines):
                if not has_line_flux[-(i+1)]:
                    out = use_lines.pop(-(i+1))
        
        templates = templates[has_line_flux,:]
        
        self.use_lines = use_lines
        
        return (use_lines, fancy), templates
            
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

    #
    def stats(self, return_string=False):
        """
        Get some information about the spectrum: contamination, fcover, etc.
        """
        
        if self.twod.im['SCI'].data.max() <= 0:
            print '%s: No valid pixels in 2D spectrum.' %(self.twod.im.filename())
            self.status = False
            return False
        
        #
        self.oned_wave = self.oned.data.wave
        
        lower, upper = 1.15e4, 1.6e4
        wuse = (self.oned_wave > lower) & (self.oned_wave < upper)
        non_zero = self.oned.data.flux != 0
        
        if (wuse & non_zero).sum() < 3:
            print '%s: No valid pixels in 2D spectrum.' %(self.twod.im.filename())
            self.status = False
            return False
        
        ### Collapse along wave axis for fcover
        spec = np.sum(self.twod.im['SCI'].data, axis=0)
        wmin = np.maximum(self.oned_wave[spec != 0].min(), lower) 
        wmax = np.minimum(self.oned_wave[spec != 0].max(), upper)
        
        #### Fraction of bad pixels in 2D spec within wmin < wave < wmax
        #### -- for the high sky-background GOODS-N images, there might only be one of the
        ####    four dithered exposures that contributes to the spectrum.
        wregion = (self.oned_wave >= wmin) & (self.oned_wave <= wmax)
        test = self.twod.im['SCI'].data[:, wregion] == 0.
        F_FLAGGED = test.sum()*1./test.size
        
        #### Wavelength coverage with respect to the nominal range `lower` > wave > `upper`
        F_COVER = (wmax-wmin)/(upper-lower)
        
        #### Maximum contamination within the nominal wavelength range
        MAX_CONTAM = (self.oned.data.contam/self.oned.data.flux)[wuse & non_zero].max()
        
        #### Average contamination over the spectrum = sum(contam) / sum(flux)
        sum_flux = np.trapz(self.oned.data.flux[wuse & non_zero], self.oned_wave[wuse & non_zero])
        sum_contam = np.trapz(self.oned.data.contam[wuse & non_zero], self.oned_wave[wuse & non_zero])
        
        INT_CONTAM = sum_contam / sum_flux
        
        #### Number of negative contamination-corrected pixels as a proxy for a bad contamination
        #### correction
        is_negative = (self.oned.data.flux-self.oned.data.contam) < 0
        F_NEGATIVE = (is_negative & wuse & non_zero).sum()*1. / (wuse & non_zero).sum()
        
        #### Q_z of the photo-z fit
        Q_Z = self.zout.q_z[self.ix]
        
        #### MAG_AUTO from the SExtractor catalog
        DIRECT_MAG = self.twod.im[0].header['MAG']
        
        if return_string:
            header = '# id       mag   q_z     f_cover f_flagged max_contam int_contam f_negative\n'
            params = '%s %.3f  %.2e  %.2f  %.2f %6.2f %6.2f  %.2f\n' %(self.root, DIRECT_MAG, Q_Z, F_COVER, F_FLAGGED, MAX_CONTAM, INT_CONTAM, F_NEGATIVE)
            return [header, params]
        else:    
            return DIRECT_MAG, Q_Z, F_COVER, F_FLAGGED, MAX_CONTAM, INT_CONTAM, F_NEGATIVE
#
def _objective_lineonly(coeffs, observed, var, templates):
    ### The "minimum" function limits the exponent to acceptable float values
    flux_fit = np.dot(np.exp(np.minimum(coeffs,345)).reshape((1,-1)), templates)
    lnprob = -0.5*np.sum((observed-flux_fit)**2/var)
    return lnprob
#
def _objective_z_line(params, observed, var, GrisModel):
    
    #### Get the lines at the new redshift    
    (use_lines, fancy), templates = GrisModel.get_emline_templates(ztry=params[-2], line_width=np.exp(params[-1]), use_determined_lines=True)
    
    ### The "minimum" function limits the exponent to acceptable float values
    flux_fit = np.dot(np.exp(np.minimum(params[:-2],345)).reshape((1,-1)), templates)
    lnprob = -0.5*np.sum((observed-flux_fit)**2/var)
    
    print 'z_try: %.4f, dv = %.1f, lnprob: %.2f' %(params[-2], np.exp(params[-1]), lnprob)
    
    return lnprob
        
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
    
#
from multiprocessing import Pool, Process
from threading import Thread

def test_threading(nproc=4, N=5):
    """
    Testing parallel threading, nothing promising because it seems like the threading/multiprocessing bookkeeping overhead is significant compared to the "model.compute_model()" execution time.
    """
    import stsci.convolve
    
    ext = unicorn.reduce.Interlace2D('COSMOS-12_00038.2D.fits')
    #point = unicorn.reduce.Interlace2D('zSPEC/COSMOS-12_00348.2D.fits')
    
    #N = 200
    
    import time
    t0 = time.time()
    models = []
    for i in range(N):
        for j in range(5):
            ext.compute_model()
            models.append(ext.model.copy())
    
    t1 = time.time()
    print t1-t0
    
    #### threading
    for i in range(N):
        models = []
        for j in range(5):
            current = threadedModel(ext)
            models.append(current)
            current.start()
        
        for model in models:
            model.join()
        
    print np.std(model.model)
        
    t2 = time.time()
    print t2-t1
    
    #### Multi Process
    for i in range(N):
        models = []
        for j in range(5):
            current = multiModel(ext)
            models.append(current)
            current.start()
        
        for model in models:
            model.join()
        
    print np.std(model.model)
        
    t3 = time.time()
    print t3-t2
    
    #### Multiprocessing
    for i in range(N):
        pool = Pool(processes=nproc)
        results = []
        for j in range(5):
            result = pool.apply_async(poolModel, [ext,j])
            out = result.get()
        pool.terminate()
        
        # #
        # for result in results:
        #     out = result.get(timeout=10)
    
    #out = results[-1].get(timeout=10)
    print np.std(out)
    
    t4 = time.time()
    print t4-t3
    
    print t1-t0, t2-t1, t3-t2
    
    
def poolModel(twod, f):
    twod.compute_model()
    return twod.model.copy()
#
class multiModel(Process):
    def __init__(self, spec_2d):
        Process.__init__(self)
        self.spec_2d = spec_2d
        self.status = False
    #
    def run(self):
        self.spec_2d.compute_model()
        self.model = self.spec_2d.model.copy()
        self.status = True
    
class threadedModel(Thread):
    def __init__(self, spec_2d):
        Thread.__init__(self)
        self.spec_2d = spec_2d
        self.status = False
    #
    def run(self):
        self.spec_2d.compute_model()
        self.model = self.spec_2d.model.copy()
        self.status = True

    
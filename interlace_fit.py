"""
Fit templates to the interlaced grism spectra to measure redshifts, and....

"""
import os

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

##### First option is STSCI_PYTHON v2.12 / Python2.7.  Second is older Python2.5.4
try:
    import stsci.convolve as convolve
except:
    try:
        import convolve
    except:
        print 'No stsci.convolve found.'
        
import numpy as np
import scipy.stats as stats
from scipy import polyval

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
        pointing = point.split('_inter_model')[0]
        #pointing = 'UDS-18'
        model = unicorn.reduce.process_GrismModel(pointing)
        bright = model.cat.mag < 23
        ids = model.cat.id[bright]
        for id in ids:     
            root='%s_%05d' %(pointing, id)
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #gris = unicorn.interlace_fit.GrismSpectrumFit(root='../GOODS-S-34_%05d' %(id))
            if not os.path.exists(root+'.2D.fits'):
                model.twod_spectrum(id)
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            #if gris.dr > 1:
            #    continue
            #
            print '\n'
            status = gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
            if status:
                gris.new_fit_free_emlines()
            #gris.make_figure()
        
class GrismSpectrumFit():
    """
    Functions for fitting (redshifts for now) the interlaced grism spectra
    
    gris = unicorn.interlace_fit.GrismSpectrumFit('GOODS-S-34_00280')
    gris.fit_in_steps()     ## fit redshift
    gris.fit_free_emlines() ## fit emission lines
    
    """
    def __init__(self, root='GOODS-S-34_00280', FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.55, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=1.e-4, use_mag_prior=True, dr_match=1., fast=False, contam_ferror=0.1, flatten_thumb=False):
        """
        Read the 1D/2D spectra and get the photometric constraints
        necessary for the spectrum fits.
        
        The "fix_direct_thumbnail" option runs the `interpolate_direct_thumb`
        method do interpolate over missing pixels in the direct thumbnail. 
        """
        self.FIGURE_FORMAT=FIGURE_FORMAT
        self.root = os.path.basename(root)

        plt.rcParams['lines.marker'] = ''
        plt.rcParams['lines.linestyle'] = '-'        
        
        self.use_lines = None

        self.pointing = root.split('_')[0]
        self.field = '-'.join(root.split('-')[:-1])
        self.contam_ferror = contam_ferror
        
        if RELEASE:
            self.OUTPUT_PATH='refit'
        else:
            self.OUTPUT_PATH=OUTPUT_PATH
        
        #### Get the 1D/2D spectra
        if RELEASE:
            # BASE_PATH = '%s/../Release/v2.0/%s' %(unicorn.GRISM_HOME, self.field)            
            #BASE_PATH = '/3DHST/Spectra/Release/v2.1/%s/%s-WFC3_v2.1_SPECTRA/' %(self.field, self.field)
            
            #self.twod = unicorn.reduce.Interlace2D('%s/%s/2D/FITS/%s.2D.fits' %(BASE_PATH, self.pointing, root), PNG=False)
            BASE_PATH = '/Volumes/3DHST_Gabe/RELEASE_v4.0/Spectra/%s/WFC3/' %(self.field)
            BASE_PATH = '/Volumes/3DHST_Gabe/RELEASE_v4.0/Spectra/%s/WFC3/' %(self.field)
            BASE_PATH = '/Volumes/3DHST_Gabe/3DHST/v4.1.4/%s-wfc3-spectra_v4.1.4' %(self.field.split('-')[0])
            
            self.twod = unicorn.reduce.Interlace2D('%s/%s/2D/FITS/%s.2D.fits' %(BASE_PATH, self.pointing.split('-G141')[0], root), PNG=False, flatten_thumb=flatten_thumb)

        else:
            self.twod = unicorn.reduce.Interlace2D(root+'.2D.fits', PNG=False, flatten_thumb=flatten_thumb)            
        #
        if 'GRISM' not in self.twod.im[0].header.keys():
            self.grism_element = 'G141'
            self.direct_filter = 'F140W'
            self.twod.im[0].header.update('FILTER', 'F140W')
        else:
            self.grism_element = self.twod.im[0].header['GRISM']
            self.direct_filter = self.twod.im[0].header['FILTER']

        
        if fix_direct_thumbnail:
            self.interpolate_direct_thumb()
            
        self.status = True
        if self.twod.im['SCI'].data.max() <= 0:
            if verbose:
                threedhst.showMessage('%s: \nNo non-zero pixels in the 2D spectrum.' %(root), warn=True)
            self.status = False
            return None
        
        if RELEASE:
            #self.oned = unicorn.reduce.Interlace1D('%s/%s/1D/FITS/%s.1D.fits' %(BASE_PATH, self.pointing, root), PNG=False)
            self.oned = unicorn.reduce.Interlace1D('%s/%s/1D/FITS/%s.1D.fits' %(BASE_PATH, self.pointing.split('-G141')[0], root), PNG=False)
            
            #print '%s/%s/1D/FITS/%s.1D.fits' %(BASE_PATH, self.pointing, root)
        else:
            #self.oned = unicorn.reduce.Interlace1D(root+'.1D.fits', PNG=False)
            self.oned = unicorn.reduce.Interlace1D(self.twod.file.replace('2D','1D'), PNG=False)
            
        if self.oned.data.flux.max() <= 0:
            print '%s: No valid pixels in 1D spectrum.' %(root)
            self.status = False
            return None
        
        ##### Fast fitting
        if fast:
            self.use_fast = True
            self.twod.init_fast_model()
            self.compute_model_function = self.twod.fast_compute_model
        else:
            self.use_fast = False
            self.compute_model_function = self.twod.compute_model
            
        #### Convert to 10**-17 ergs / s / cm**2 / A
        #self.oned.data.sensitivity *= np.diff(self.oned.data.wave)[0]
        
        self.grism_id = os.path.basename(self.twod.file.split('.2D.fits')[0])

        #### Get the photometric match.
        ra, dec = self.oned.header['RA'], self.oned.header['DEC']
        if skip_photometric:
            self.phot_zgrid = np.arange(0,6,0.1)
            self.phot_lnprob = self.phot_zgrid*0.
            self.dr = 10
            self.skip_photometric=True
            self.zout = None
            self.cat = None
            
        else:
            self.get_photometric_constraints(ra=ra, dec=dec, verbose=verbose, p_flat=p_flat)
            self.skip_photometric=False

        #### If match distance > 1", prior is flat
        if self.dr > dr_match:
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
        
        #### EAZY (K) magnitude prior, OK for H
        if use_mag_prior & (self.dr > 1):
            #### Use EAZY magnitude prior
            prior = open('%s/data/prior_K_zmax7_coeff.dat' %(os.path.dirname(threedhst.__file__))).readlines()
            prior_mag = np.cast[float](prior[0].split()[2:])
            prior_z0 = np.cast[float](prior[1].split()[1:])
            prior_gamma = np.cast[float](prior[2].split()[1:])
            z0 = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_z0, left=prior_z0[0], right=prior_z0[-1])
            gamma = np.interp(self.twod.im[0].header['MAG'], prior_mag, prior_gamma, left=prior_gamma[0], right=prior_gamma[-1])
            prior_pz = self.phot_zgrid**gamma * np.exp(-(self.phot_zgrid / z0)**gamma)
            self.phot_lnprob += np.log(prior_pz / np.trapz(prior_pz, self.phot_zgrid))
                    
        #### Initialize the continuum and emission line templates    
        self.line_free_template()
        self.linex, self.liney = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.loOIII.txt', unpack=True)
        self.linex2, self.liney2 = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.hiOIII.txt', unpack=True)
        # self.linex, self.liney = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
        # self.linex2, self.liney2 = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
        
        #### Smooth the line templates
        if True:
            threedhst.showMessage('Smooth template')
            import scipy.ndimage as nd
            
            s = 10
            xx = np.arange(100,8.e4)
            yy = nd.gaussian_filter1d(np.interp(xx, self.linex, self.liney), s)
            self.liney = np.interp(self.linex, xx, yy)
            
            yy = nd.gaussian_filter1d(np.interp(xx, self.linex2, self.liney2), s)
            self.liney2 = np.interp(self.linex2, xx, yy)
                              
        #### Sensitivity error component
        self.sens_var = 0.
        self.twod.get_sensitivity_error(scale=1./20)
        self.sens_var = self.twod.sens_var.flatten()
        
        #### Try to read the previous p(z) from the pickle file
        if RELEASE:
            pz_path = '%s/%s/ZFIT/PZ/' %(BASE_PATH, self.pointing)
        else:
            pz_path = './'
                
        if self.load_fits(path=pz_path):
            if verbose:
                print 'Read p(z) from %s.zfit.pz.fits' %(self.grism_id)
                self.get_best_fit()
        
        ### Needed for interlace_test
        self.coeffs = 0.
        
    def get_photometric_constraints(self, ra=0., dec=0., verbose=True, p_flat=0.):
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
        print self.grism_id
        cat, zout, fout = unicorn.analysis.read_catalogs(root=self.grism_id)
        #cat.kmag = 23.86-2.5*np.log10(cat.field(KTOT_COL))
        CAT_PATH = os.path.dirname(cat.filename)
        
        root = os.path.basename(zout.filename).split('.zout')[0]
        ZOUT_PATH = os.path.dirname(zout.filename)
        
        self.eazy_tempfilt, self.eazy_coeffs, self.eazy_temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE = root,                                                 OUTPUT_DIRECTORY=ZOUT_PATH, CACHE_FILE = 'Same')
        
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
            print 'Match in %s\n #%s  dr=%.2f"  z_spec=%7.3f  z_peak=%7.3f' %(cat.filename, cat.id[ix], self.dr, zout.z_spec[ix], self.z_peak)
        #                    
        self.best_fit = np.dot(self.eazy_temp_sed['temp_seds'], self.eazy_coeffs['coeffs'][:,ix])
        self.templam = self.eazy_temp_sed['templam']*1.
        
        self.best_fit_flux = self.full_eazy_template(self.eazy_temp_sed['templam'], self.eazy_temp_sed['temp_seds'], self.z_peak)
        
        # self.phot_fnu = tempfilt['fnu'][:,ix]
        # self.phot_efnu = tempfilt['efnu'][:,ix]
        # self.phot_lc = tempfilt['lc']
        # 
        # self.tempfilt = tempfilt
        # self.tempgrid = tempfilt['tempfilt']

        self.phot_zgrid = self.eazy_tempfilt['zgrid']
        self.phot_lnprob = -0.5*pz['chi2fit'][:,ix]
        self.phot_lnprob -= self.phot_lnprob.max()
        
        ### normalized probability
        self.phot_linear = np.exp(self.phot_lnprob)
        norm  = np.trapz(self.phot_linear, self.phot_zgrid)
        self.phot_linear /= norm
        self.phot_lnprob -= np.log(norm)
        
        #### Total integrated probability in the flat prior, sets a floor
        #### on the minimum photometric p(z)
        #p_flat = 1.e-4
        if p_flat > 0:
            y0 = self.phot_zgrid*0.+np.log(p_flat/np.sum(np.diff(self.phot_zgrid)))
            y0[self.phot_lnprob < -200] = -200
            self.phot_lnprob = np.maximum(self.phot_lnprob-self.phot_lnprob.max(), y0)
            norm  = np.trapz(np.exp(self.phot_lnprob), self.phot_zgrid)
            self.phot_lnprob -= np.log(norm)
            self.phot_linear = np.exp(self.phot_lnprob)
        
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
        var += (self.contam_ferror*self.twod.im['CONTAM'].data**2).flatten()
        
        var += self.sens_var
        
        # 
        # w2d = np.dot(np.ones((self.twod.im['CONTAM'].data.shape[0],1)), self.twod.im['WAVE'].data.reshape(1,-1)).flatten()
        # var[w2d < 1.08e4] *= 10
        # var[w2d > 1.68e4] *= 10

        flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        use = np.isfinite(flux)

        show = False

        #### Templates to allow blue/red tilt to the spectrum with respect to the EAZY templates.
        #### These are simple lines normalized at lam=1.4e4.  That it works appears a bit magical, but
        #### it does seem to.  Only need to compute them once since they don't change with the 
        #### redshift grid.
        if self.grism_element == 'G141':
            x0, dx = 1.4e4, 4000.
        else:
            x0, dx = 0.98e4, 2800.
                    
        tilt_red = ((self.templam_nolines-x0)/dx+1)*2
        self.compute_model_function(self.templam_nolines, tilt_red)
        tilt_red_model = self.twod.model*1.
                
        tilt_blue = (-(self.templam_nolines-x0)/dx+1)*2
        self.compute_model_function(self.templam_nolines, tilt_blue)
        tilt_blue_model = self.twod.model*1.
        #tilt_red_model = tilt_blue_model*1.
        
        
        #### Loop through redshift grid
        pmax = -1.e10
        zbest = 0
        for i in range(NZ):
            if verbose:
                print unicorn.noNewLine+'z: %.3f  %.3f' %(zgrid[i], zbest)
            
            z_test = zgrid[i]

            #### Generate the 2D spectrum model for the continuum and emission line templates
            if self.best_fit_nolines.sum() == 0.:
                self.best_fit_nolines = self.best_fit_nolines*0.+1
            
            self.compute_model_function(self.templam_nolines*(1+z_test), self.best_fit_nolines)
            continuum_model = self.twod.model*1.
                
            self.compute_model_function(self.linex*(1+z_test), self.liney)
            line_model = self.twod.model*1.
            ixl = [1]
            
            #### Initialize if on the first grid point
            if i == 0:
                templates = np.ones((5,line_model.size))
                templates[2,:] = tilt_blue_model.flatten()
                templates[3,:] = tilt_red_model.flatten()
                
                #### Line and continuum indices
                ixl = np.array([1,4])
                ixc = np.array([0,2,3])
                
                # if True:
                #     threedhst.showMessage('fit background', warn=True)
                #     templates = np.ones((6,line_model.size))
                #     templates[2,:] = tilt_blue_model.flatten()
                #     templates[3,:] = tilt_red_model.flatten()
                #     templates[5,:] = tilt_red_model.flatten()*0+0.01
                # 
                #     #### Line and continuum indices
                #     ixl = np.array([1,4])
                #     ixc = np.array([0,2,3,5])
                
            #### Put in the spectral templates
            templates[0,:] = continuum_model.flatten()
            ### Don't use slope
            # templates[2,:] = continuum_model.flatten()
            # templates[3,:] = continuum_model.flatten()
            templates[1,:] = line_model.flatten()
            
            #### Second line template
            self.compute_model_function(self.linex2*(1+z_test), self.liney2)
            templates[4,:] = (self.twod.model*1.).flatten()
                        
            ##### Probably no flux in the direct image
            if templates.max() == 0:
                self.status = False
                return False
                
            #### Compute the non-negative normalizations of the template components
            #print var.shape, use.shape, templates.shape, self.twod.model.shape
            amatrix = utils_c.prepare_nmf_amatrix(var[use], templates[:,use])
            coeffs = utils_c.run_nmf(flux[use], var[use], templates[:,use], amatrix, toler=1.e-5)
            flux_fit = np.dot(coeffs.reshape((1,-1)), templates).flatten() #.reshape(line_model.shape)

            #### If `get_model_at_z`, return the model at that redshift
            if (get_model):
                #threedhst.showMessage(coeffs.__str__(), warn=True)
                
                self.cont_model = np.dot(coeffs[ixc].reshape((1,-1)), templates[ixc,:]).reshape(line_model.shape)
                self.slope_model = np.dot(coeffs[2:4].reshape((1,-1)), templates[2:4,:]).reshape(line_model.shape)
                self.line_model = np.dot(coeffs[ixl].reshape((1,-1)), templates[ixl,:]).reshape(line_model.shape)
                #self.line_model = (coeffs[1]*templates[1,:]).reshape(line_model.shape)
                self.flux_model = flux_fit.reshape(line_model.shape)
                self.oned_wave, self.model_1D = self.twod.optimal_extract(self.flux_model)
                oned_wave_x, self.cont_1D = self.twod.optimal_extract(self.cont_model)
                oned_wave_x, self.slope_1D = self.twod.optimal_extract(self.slope_model)
                oned_wave_x, self.line_1D = self.twod.optimal_extract(self.line_model)
                
                return True
                
                #return flux_model, cont_model, line_model, oned_wave, model_1D, cont_1D, line_1D, slope_1D

            #### Log likelihood at redshift zgrid[i] (-0.5*chi2)
            #print flux.shape, flux_fit.shape, var.shape
            spec_lnprob[i] = -0.5*np.sum(((flux-flux_fit)**2/var)[use])
            if spec_lnprob[i] > pmax:
                pmax = spec_lnprob[i]
                zbest = z_test*1.
                
        #### Done!
        return zgrid, spec_lnprob
    
    def fit_zgrid_oned(self, zrange = (0.01,6), dz = 0.02, get_model_at_z=None, verbose=True):
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
        var = self.twod.oned.data['error']**2
        var[var == 0] = 1.e6
        var += self.contam_ferror*self.twod.oned.data['contam']**2
        
        sens_int = np.interp(self.twod.im['WAVE'].data, unicorn.reduce.sens_files['A']['WAVELENGTH'], unicorn.reduce.sens_files['A']['SENSITIVITY'].max()*unicorn.reduce.sens_files['A']['ERROR']/unicorn.reduce.sens_files['A']['SENSITIVITY']**2)
        
        xx, flux = self.twod.optimal_extract(np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data))
        
        var += (flux*sens_int)**2
        
        use = np.isfinite(flux) & (var < 1.e6)

        show = False

        #### Templates to allow blue/red tilt to the spectrum with respect to the EAZY templates.
        #### These are simple lines normalized at lam=1.4e4.  That it works appears a bit magical, but
        #### it does seem to.  Only need to compute them once since they don't change with the 
        #### redshift grid.
        if self.grism_element == 'G141':
            x0, dx = 1.4e4, 4000.
        else:
            x0, dx = 0.98e4, 2800.
                    
        tilt_red = ((self.templam_nolines-x0)/dx+1)*2
        self.compute_model_function(self.templam_nolines, tilt_red)
        tilt_red_model = self.twod.model*1.
        xx, tilt_red_model_1D = self.twod.optimal_extract(tilt_red_model)
        
        tilt_blue = (-(self.templam_nolines-x0)/dx+1)*2
        self.compute_model_function(self.templam_nolines, tilt_blue)
        tilt_blue_model = self.twod.model*1.
        xx, tilt_blue_model_1D = self.twod.optimal_extract(tilt_blue_model)
        
        #tilt_red_model = tilt_blue_model*1.
        #### Loop through redshift grid
        pmax = -1.e10
        zbest = 0
        for i in range(NZ):
            if verbose:
                print unicorn.noNewLine+'z: %.3f  %.3f' %(zgrid[i], zbest)
            
            z_test = zgrid[i]

            #### Generate the 2D spectrum model for the continuum and emission line templates
            if self.best_fit_nolines.sum() == 0.:
                self.best_fit_nolines = self.best_fit_nolines*0.+1
            
            self.compute_model_function(self.templam_nolines*(1+z_test), self.best_fit_nolines)
            continuum_model = self.twod.model*1.
            xx, continuum_model_1D = self.twod.optimal_extract(continuum_model)
                
            self.compute_model_function(self.linex*(1+z_test), self.liney)
            line_model = self.twod.model*1.
            xx, line_model_1D = self.twod.optimal_extract(line_model)
            ixl = [1]
            
            #### Initialize if on the first grid point
            if i == 0:
                templates = np.ones((5,line_model.size))
                templates[2,:] = tilt_blue_model.flatten()
                templates[3,:] = tilt_red_model.flatten()

                templates_1D = np.ones((5,line_model_1D.size))
                templates_1D[2,:] = tilt_blue_model_1D
                templates_1D[3,:] = tilt_red_model_1D


            #### Put in the spectral templates
            templates[0,:] = continuum_model.flatten()
            templates_1D[0,:] = continuum_model_1D
            ### Don't use slope
            # templates[2,:] = continuum_model.flatten()
            # templates[3,:] = continuum_model.flatten()
            templates[1,:] = line_model.flatten()
            templates_1D[1,:] = line_model_1D
            
            #### Second line template
            self.compute_model_function(self.linex2*(1+z_test), self.liney2)
            templates[4,:] = (self.twod.model*1.).flatten()
            xx, templates_1D[4,:] = self.twod.optimal_extract(self.twod.model)
            ixl = [1,4]
            
            ##### Probably no flux in the direct image
            if templates.max() == 0:
                self.status = False
                return False
                
            #### Compute the non-negative normalizations of the template components
            #print var.shape, use.shape, templates.shape, self.twod.model.shape
            amatrix = utils_c.prepare_nmf_amatrix(var[use], templates_1D[:,use])
            coeffs = utils_c.run_nmf(flux[use], var[use], templates_1D[:,use], amatrix, toler=1.e-5)
            flux_fit = np.dot(coeffs.reshape((1,-1)), templates) #.reshape(line_model.shape)
            flux_fit_1D = np.dot(coeffs.reshape((1,-1)), templates_1D).flatten() #.reshape(line_model.shape)
            
            #### If `get_model_at_z`, return the model at that redshift
            if (get_model):
                ixc = np.array([0,2,3])
                self.cont_model = np.dot(coeffs[ixc].reshape((1,-1)), templates[ixc,:]).reshape(line_model.shape)
                self.slope_model = np.dot(coeffs[2:4].reshape((1,-1)), templates[2:4,:]).reshape(line_model.shape)
                self.line_model = np.dot(coeffs[ixl].reshape((1,-1)), templates[ixl,:]).reshape(line_model.shape)
                #self.line_model = (coeffs[1]*templates[1,:]).reshape(line_model.shape)
                self.flux_model = flux_fit.reshape(line_model.shape)
                self.oned_wave, self.model_1D = self.twod.optimal_extract(self.flux_model)
                oned_wave_x, self.cont_1D = self.twod.optimal_extract(self.cont_model)
                oned_wave_x, self.slope_1D = self.twod.optimal_extract(self.slope_model)
                oned_wave_x, self.line_1D = self.twod.optimal_extract(self.line_model)
                
                return True
                
                #return flux_model, cont_model, line_model, oned_wave, model_1D, cont_1D, line_1D, slope_1D

            #### Log likelihood at redshift zgrid[i] (-0.5*chi2)
            #print flux.shape, flux_fit.shape, var.shape
            #spec_lnprob[i] = -0.5*np.sum(((flux-flux_fit_1D)**2/var)[use])
            spec_lnprob[i] = -0.5*np.sum(((flux-flux_fit_1D)**2 / (var))[use])
            #print use.shape, ((flux-flux_fit_1D)**2 / (var+(flux_fit_1D*sens_int)**2)).shape
            
            if spec_lnprob[i] > pmax:
                pmax = spec_lnprob[i]
                zbest = z_test*1.
                
        #### Done!
        return zgrid, spec_lnprob
        
    def fit_in_steps(self, dzfirst=0.003, zrfirst = (0.00,3.8), dzsecond=0.0002, save=True, make_plot=True, skip_second=False, oned=False):
        """
        Do two fit iterations, the first on a coarse redshift grid over the full z=(0.01,6)
        and the second refined grid around the peak found in the first iteration
        
        If `save` is set, save the results to the ascii and pickle files (.dat, .pkl).
        
        If `make_plot` is set, make the diagnostic plot.
        """
        import scipy.ndimage as nd
        zout = self.zout
        
        if oned:
            fit_function = self.fit_zgrid_oned
        else:
            fit_function = self.fit_zgrid
            
        print '\n'
        
        #################
        #### First iteration, fit full range at low dz resolution
        #################
        result = fit_function(zrange=zrfirst, dz=dzfirst)
        if result is False:
            self.status = False
            return False
        
        zgrid0, spec_lnprob0 = result
        
        #### Smooth it with a gaussian???
        wz = 0.004
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzfirst)
        yg = 1./np.sqrt(2*np.pi*wz**2)*np.exp(-xg**2/2/wz**2)
        sm = np.log(nd.convolve1d(np.exp(spec_lnprob0-spec_lnprob0.max()), yg/yg.sum(), mode='constant', cval=0.))
        #spec_lnprob0 = (sm)
        
        #### Interpolate the photometric p(z) to apply it as a prior
        phot_int0 = np.interp(zgrid0, self.phot_zgrid,  self.phot_lnprob, left=-1.e8, right=-1.e8)
        phot_int0[~np.isfinite(phot_int0)] = -1.e8
        full_prob0 = phot_int0 + spec_lnprob0 - spec_lnprob0.max()
        full_prob0 -= full_prob0.max()  ### normalize p_max = 1
        
        #print full_prob0, phot_int0, spec_lnprob0
        
        z_max = zgrid0[full_prob0 == full_prob0.max()][0]
        zsub = full_prob0 > np.log(1.e-10)

        #################
        #### Second iteration, small dz steps around the probability peak
        #################
        if not skip_second:
            #### Define the range and step size for the next iteration
            width = 0.002*(1+z_max)
            #zrange, dz = (z_max-width, z_max+width), dzsecond
            zsecond = np.arange(zrfirst[0], zrfirst[1], dzsecond)
            pzint = np.interp(zsecond, zgrid0, full_prob0)
            #
            #### Smooth it with a gaussian???
            wz = 0.008
            xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
            yg = np.exp(-xg**2/2/wz**2)
            sm = np.log(nd.convolve1d(np.exp(pzint), yg/yg.sum(), mode='constant', cval=0.))
            
            #zsub = pzint > np.log(1.e-5)
            zsub = (sm-sm.max()) > np.log(1.e-5)
            if zsub.sum() == 0:
                threedhst.showMessage('Something went wrong with the redshift grid...', warn=True)
                print pzint.max(), pzint.min(), full_prob0.max()
                return False
                
            if (zsecond[zsub].max() - zsecond[zsub].min()) < width*2:
                zrange = (z_max-width, z_max + width)
            else:
                zrange = (zsecond[zsub].min(), zsecond[zsub].max())
            
            print '\nStep 1: z_max = %.3f (%.3f,%.3f)\n' %(z_max, zrange[0],
                                                           zrange[1])
            zgrid1, spec_lnprob1 = fit_function(zrange=zrange, dz=dzsecond)
        else:
            zgrid1, spec_lnprob1 = zgrid0, spec_lnprob0
            
        #### Interpolate the photometric p(z) to apply it as a prior
        phot_int1 = np.interp(zgrid1, self.phot_zgrid,  self.phot_lnprob)
        phot_int1[~np.isfinite(phot_int1)] = -1000
        full_prob1 = phot_int1 + spec_lnprob1 - spec_lnprob1.max()
        #print full_prob1, phot_int1, spec_lnprob1
        
        full_prob1 -= full_prob1.max()
        
        #### Smooth the high resolution p(z) because otherwise get "ringing" behavior related
        #### to the discrete input pixels
        wz = 0.0005
        xg = np.arange(-5*wz, 5*wz+1.e-6, dzsecond)
        yg = np.exp(-xg**2/2./wz**2)
        sm = np.log(nd.convolve1d(np.exp(full_prob1), yg/yg.sum(), mode='constant', cval=0.))

        full_prob1 = sm-sm.max()
        
        
        self.z_max_spec = zgrid1[full_prob1 == full_prob1.max()][0]
        
        #### Save results and make the diagnostic plot
        self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1 = zgrid0, full_prob0, zgrid1, full_prob1
        
        ### Get the best 2D/1D fit
        self.get_best_fit()
        
        if save:
            ### Write the ascii file
            self.save_results()
            
            ### Save p(z) to a pickle
            self.save_fits()    
        
        if make_plot:
            ### Make diagnostic figure
            self.make_figure()
            self.twod_figure()
            self.make_spectrum_fits()
            
        return True
        
    def fit_mcmc():
        """
        Fit on a grid just fine enough to find lines.  Then run the MCMC sampler with lines
        and continuum templates free
        """
        pass
        
    def save_results(self, verbose=True):
        """
        Print the results of the fit to a file, like 
        GOODS-S-34_00622.zfit.dat.

        """
        zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
        
        z_max_spec = zgrid1[full_prob1 == full_prob1.max()][0]
        z_peak_spec = np.trapz(zgrid1*np.exp(full_prob1), zgrid1)/np.trapz(np.exp(full_prob1), zgrid1)

        fp = open(self.OUTPUT_PATH + '/' + self.grism_id+'.zfit.dat','w')
        fp.write('#  spec_id   phot_id   dr   z_spec  z_peak_phot  z_max_spec z_peak_spec\n')
        
        if self.skip_photometric:
            fp.write('#  Phot: None\n')
            file_string = '%s %d  %.3f  %.5f  %.5f  %.5f  %.5f\n' %(self.grism_id, -1, -1, -1, -1, z_max_spec, z_peak_spec)
        else:
            fp.write('#  Phot: %s\n' %(self.cat.filename))
            file_string = '%s %s  %.3f  %.5f  %.5f  %.5f  %.5f\n' %(self.grism_id, self.cat.id[self.ix], self.dr, self.zout.z_spec[self.ix], self.zout.z_peak[self.ix], z_max_spec, z_peak_spec)
            
        fp.write(file_string)
        fp.close()

        if verbose:
            print file_string
    
    def save_pickle(self):
        """

        Save the redshift/probability grid to a pickle file, like
        GOODS-S-34_00622.zfit.pkl

        """
        import cPickle as pickle

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
        import cPickle as pickle
        
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
    
    def save_fits(self):
        
        header = pyfits.Header()
        if self.skip_photometric:
            header.update('PHOTID', -1)
        else:
            header.update('PHOTID', self.cat.id[self.ix], comment='ID in %s' %(os.path.basename(self.cat.filename)))
        header.update('GRISID', self.grism_id, comment='Grism ID')
        
        hdu = [pyfits.PrimaryHDU(header=header)]

        hdu.append(pyfits.ImageHDU(data=self.phot_zgrid, name='PZGRID'))
        hdu.append(pyfits.ImageHDU(data=self.phot_lnprob, name='LN_PZ_PROB'))
                
        hdu.append(pyfits.ImageHDU(data=self.zgrid0, name='ZGRID0'))
        hdu.append(pyfits.ImageHDU(data=self.full_prob0, name='LN_PROB_0'))
        hdu.append(pyfits.ImageHDU(data=self.zgrid1, name='ZGRID1'))
        hdu.append(pyfits.ImageHDU(data=self.full_prob1, name='LN_PROB_1'))
        hduList = pyfits.HDUList(hdu)
        hduList.writeto(self.OUTPUT_PATH + '/' + self.grism_id+'.zfit.pz.fits', clobber=True, output_verify='silentfix')
    
    def load_fits(self, path='./'):
        
        file = path + self.grism_id+'.zfit.pz.fits'
        #file = self.twod.file.replace('2D','new_zfit.pz')
        #threedhst.showMessage(file, warn=True)
        
        if not os.path.exists(file):
            return False
            
        im = pyfits.open(file)
        self.zgrid0 = im['ZGRID0'].data
        self.full_prob0 = im['LN_PROB_0'].data
        self.zgrid1 = im['ZGRID1'].data
        self.full_prob1 = im['LN_PROB_1'].data
        im.close()
        
        self.z_max_spec = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        return True
        
    def get_redshift_confidence(self, limits=[2.5, 16, 84, 97.5]):
        """
        Get confidence intervals interpolated from the full redshift PDF
        """
        linear_prob = np.exp(self.full_prob1*1.)
        linear_prob /= np.trapz(linear_prob, self.zgrid1)
        linear_prob_cumsum = np.cumsum(linear_prob[1:]*np.diff(self.zgrid1))
        z0 = (self.zgrid1[:-1]+self.zgrid1[1:])/2.
        zi = np.interp(np.array(limits)/100., linear_prob_cumsum, z0)
        return zi
        
    def get_best_fit(self):
        """
        Get the 2D fit at the redshift of peak probability
        
        """
        z_max_spec = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        #self.flux_model, self.cont_model, self.line_model, self.oned_wave, self.model_1D, self.cont_1D, self.line_1D, self.slope_1D = self.fit_zgrid(get_model_at_z=z_max_spec, verbose=False)
        status = self.fit_zgrid(get_model_at_z=z_max_spec, verbose=False)
    
    def make_spectrum_fits(self, base='zfit'):
        """
        Make an output FITS file containing the continuum and line models, 1D and 2D
        """
        import astropy.io.fits as fits
        
        zout = self.zout
        #zgrid0, full_prob0, zgrid1, full_prob1 = self.zgrid0, self.full_prob0, self.zgrid1, self.full_prob1
        
        prim = fits.PrimaryHDU()
        prim.header = self.twod.im[0].header.copy()
        hdus = [prim]
        
        twod_header = self.twod.im['SCI'].header.copy()
        
        #self.model_1D
        #self.oned.data.sensitivity
        #self.cont_model
        #self.line_model
        
        hdus.append(fits.ImageHDU(data=self.cont_model, header=twod_header, name='CONT2D'))
        hdus.append(fits.ImageHDU(data=self.line_model, header=twod_header, name='LINE2D'))
        hdus.append(fits.ImageHDU(data=self.cont_1D, name='CONT1D'))
        hdus.append(fits.ImageHDU(data=self.line_1D, name='LINE1D'))
        hdus.append(fits.ImageHDU(data=self.coeffs, name='COEFFS'))
        
        fits.HDUList(hdus).writeto(self.OUTPUT_PATH + '/' + self.grism_id+'.%s.fits'%(base), clobber=True)
        
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
        
        if self.grism_element == 'G102':
            wuse = (self.oned.data.wave > 0.84e4) & (self.oned.data.wave < 1.15e4)
        #
        if self.grism_element == 'G800L':
            wuse = (self.oned.data.wave > 0.58e4) & (self.oned.data.wave < 0.92e4)
        
        yflux, ycont = self.oned.data.flux, self.oned.data.contam
        y = yflux-ycont
        #### Show fresh optimal extractions xxx
        # x0, yflux = self.twod.optimal_extract(self.twod.im['SCI'].data)
        # x0, ycont = self.twod.optimal_extract(self.twod.im['CONTAM'].data)
        # y = (yflux-ycont)[show]
        # self.oned.data.flux = yflux
        # self.oned.data.contam = ycont
        
        yerr = self.oned.data.error #[show]
        ax.fill_between(self.oned.data.wave[show]/1.e4, (y+yerr)[show], (y-yerr)[show], color='blue', alpha=0.1)
        
        ax.plot(self.oned.data.wave[show]/1.e4, yflux[show], color='black', alpha=0.1)
        ax.plot(self.oned.data.wave[show]/1.e4, (yflux-ycont)[show], color='black')
        
        ax.plot(self.oned_wave[show]/1.e4, self.model_1D[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave/1.e4, self.model_1D, color='red', alpha=0.08, linewidth=2)
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'e$^-$ / s')
        ok = np.isfinite(yflux) & np.isfinite(y)
        if wuse.sum() > 5:
            #ymax = yflux[wuse].max(); 
            ymax = y[ok & wuse].max()
        else:
            ymax = yflux[ok].max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax) 
        ax.set_xlim(1.0, 1.73)
        xint = [1.1,1.2,1.3,1.4,1.5,1.6]
        ax_int = np.array(xint)#*1.e4

        if self.grism_element == 'G102':
            ax.set_xlim(0.74, 1.17)
            xint = [0.8, 0.9, 1.0, 1.1]
            ax_int = np.array(xint)#*1.e4
        #
        if self.grism_element == 'G800L':
            ax.set_xlim(0.49, 1.08)
            xint = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            ax_int = np.array(xint)#*1.e4
        
        ax.set_xticks(ax_int)
        
        #### Spectrum in f_lambda
        self.oned.data.sensitivity /= 100
        
        ax = fig.add_subplot(142)
        ax.plot(self.oned.data.wave[show]/1.e4, yflux[show]/self.oned.data.sensitivity[show], color='black', alpha=0.1)
        
        show_flux = (yflux[show]-ycont[show])/self.oned.data.sensitivity[show]
        show_err = self.oned.data.error[show]/self.oned.data.sensitivity[show]
        ax.plot(self.oned.data.wave[show]/1.e4, show_flux, color='black')
        #ax.fill_between(self.oned.data.wave[show]/1.e4, show_flux+show_err, show_flux-show_err, color='0.5', alpha=0.2)
        
        ax.plot(self.oned_wave/1.e4, self.model_1D/self.oned.data.sensitivity, color='red', alpha=0.08, linewidth=2)
        ax.plot(self.oned_wave[show]/1.e4, self.model_1D[show]/self.oned.data.sensitivity[show], color='red', alpha=0.5, linewidth=2)
        ax.plot(self.oned_wave[show]/1.e4, self.slope_1D[show]/self.oned.data.sensitivity[show], color='orange', alpha=0.2, linewidth=1)
        
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'$f_\lambda$')
        if wuse.sum() > 5:
            ymax = ((yflux-ycont)/self.oned.data.sensitivity)[wuse & ok].max()
            #ymax = show_flux.max()
        else:
            ymax = (yflux[ok]/self.oned.data.sensitivity).max()
            
        ax.set_ylim(-0.05*ymax, 1.1*ymax)
        ax.set_xlim(1.0, 1.73)

        if self.grism_element == 'G102':
            ax.set_xlim(0.74, 1.17)
        #
        if self.grism_element == 'G800L':
            ax.set_xlim(0.49, 1.08)
        
        ax.set_xticks(ax_int)
        
        #### p(z)
        ax = fig.add_subplot(143)
        ax.plot(self.phot_zgrid, np.exp(self.phot_lnprob-self.phot_lnprob.max()), color='green')
        ax.plot(zgrid0, np.exp(full_prob0), color='blue', alpha=0.4)
        ax.plot(zgrid1, np.exp(full_prob1), color='blue')

        # ax.plot(self.phot_zgrid, (self.phot_lnprob-self.phot_lnprob.max()), color='green')
        # ax.plot(zgrid0, (full_prob0), color='blue', alpha=0.4)
        # ax.plot(zgrid1, (full_prob1), color='blue')
        
        if not self.skip_photometric:
            if zout.z_spec[self.ix] > 0:
                ax.plot(zout.z_spec[self.ix]*np.array([1,1]), [0,1], color='red')

        if self.dr < 1:
            ax.set_xlim(np.min([self.zgrid1.min(), zout.l99[self.ix]]), np.max([self.zgrid1.max(), zout.u99[self.ix]]))

        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'$p(z)$')
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune='both'))

        #### Make title text
        if not self.skip_photometric:
            deltaz = ''
            if zout.z_spec[self.ix] > 0:
                deltaz = '$\Delta z$ = %.4f' %(-(zout.z_spec[self.ix]-self.z_max_spec)/(1+zout.z_spec[self.ix]))
            #
            ax.text(-0.05, 1.1, r'%s  %s  $H_{140}=$%.2f $z_\mathrm{spec}$=%.3f  $z_\mathrm{phot}$=%.3f  $z_\mathrm{gris}$=%.3f  %s' %(self.grism_id, self.zout.id[self.ix], self.twod.im[0].header['MAG'], zout.z_spec[self.ix], zout.z_peak[self.ix], self.z_max_spec, deltaz), transform=ax.transAxes, horizontalalignment='center')
        else:
            ax.text(-0.05, 1.1, r'%s  $H_{140}=$%.2f  $z_\mathrm{gris}$=%.3f' %(self.grism_id, self.twod.im[0].header['MAG'], self.z_max_spec), transform=ax.transAxes, horizontalalignment='center')
            
            unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.zfit.%s' %(self.FIGURE_FORMAT))
            return True
            
        ####  Show the (photometric) SED with the spectrum overplotted
        ax = fig.add_subplot(144)

        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit

        temp_sed_int = np.interp(self.oned.data.wave, lambdaz, temp_sed)
        keep = (self.oned.data.wave > 1.2e4) & (self.oned.data.wave < 1.5e4)
        if self.grism_element == 'G102':
            keep = (self.oned.data.wave > 0.85e4) & (self.oned.data.wave < 1.05e4)
        #
        if self.grism_element == 'G800L':
            keep = (self.oned.data.wave > 0.58e4) & (self.oned.data.wave < 0.92e4)
        
        flux_spec = (self.oned.data.flux-self.oned.data.contam-self.slope_1D*0)/self.oned.data.sensitivity
        
        ### factor of 100 to convert from 1.e-17 to 1.e-19 flux units
        #anorm = np.sum(temp_sed_int[keep]*flux_spec[keep])/np.sum(flux_spec[keep]**2)
        #scale = 100.
        anorm = 1
        scale = 100
        
        ax.plot(lambdaz, temp_sed*scale, color='blue', alpha=0.5)
        ax.errorbar(lci, fobs*scale, efobs*scale, color='black', marker='o', ms=7, alpha=0.7, linestyle='None')
        ax.plot(self.oned.data.wave, flux_spec*anorm, color='red', alpha=0.3)
        bin = 4
        binned = unicorn.utils_c.interp_conserve(self.oned.data.wave[::4], self.oned.data.wave, flux_spec)
        ax.plot(self.oned.data.wave[::4], binned, color='red', alpha=0.7)
        
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
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.zfit.%s' %(self.FIGURE_FORMAT))
        
        self.oned.data.sensitivity *= 100
        
        #z_peaks = zgrid[1:-1][(full_prob[1:-1] > np.log(0.05)) & (np.diff(full_prob,2) < 0)]
        #zrange, dz = (z_peaks[0]-0.02*(1+z_peaks[0]), z_peaks[0]+0.02*(1+z_peaks[0])), 0.0002
    
    def twod_figure(self, vmax=None, base='zfit'):
        """
        Make a figure showing the raw twod-spectrum, plus contamination- and 
        continuum-subtracted versions.
        """
        
        sh = self.twod.im['SCI'].data.shape
        top = 0.055
        bottom = 0.07
        left = 0.06
        aspect = 3.*sh[0]/sh[1]/(1-(top+bottom-left))
        
        wave = self.twod.im['WAVE'].data
        if self.grism_element == 'G141':
            xint = [1.1,1.2,1.3,1.4,1.5,1.6]
            ax_int = np.interp(np.array(xint)*1.e4, wave, np.arange(wave.shape[0]))
        #
        if self.grism_element == 'G102':
            xint = [0.8, 0.9, 1.0, 1.1]
            ax_int = np.interp(np.array(xint)*1.e4, wave, np.arange(wave.shape[0]))
        #
        if self.grism_element == 'G800L':
            xint = [0.6, 0.7, 0.8, 0.9, 1.0]
            ax_int = np.interp(np.array(xint)*1.e4, wave, np.arange(wave.shape[0]))
        
        fig = unicorn.catalogs.plot_init(xs=5,aspect=aspect, left=left, right=0.02, bottom=bottom, top=top, NO_GUI=True)
        #plt.gray()
        fig.subplots_adjust(hspace=0.001)
        
        if vmax==None:
            #values = self.twod.im['SCI'].data.flatten()
            #vmax = values[np.argsort(values)][-10]
            vmax = self.flux_model.max()
            
        #### Raw spectrum
        ax = fig.add_subplot(311)
        ax.imshow(0-self.twod.im['SCI'].data, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto', cmap=plt.cm.gray)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks(ax_int)
        ax.text(0.98,0.96, '-Raw', ha='right', va='top', color='red', size=7,  transform=ax.transAxes)
        #ax.text(0.95,1+0.1*aspect,self.grism_id, transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')
        ax.set_title(self.grism_id)
        
        #### Contam-subtracted spectrum
        ax = fig.add_subplot(312)
        ax.imshow(0-self.twod.im['SCI'].data+self.twod.im['CONTAM'].data, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto', cmap=plt.cm.gray)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks(ax_int)
        ax.text(0.98,0.96, '-Contam.', ha='right', va='top', color='red', size=7, transform=ax.transAxes)

        #### Continuum-subtracted spectrum
        ax = fig.add_subplot(313)
        ax.imshow(0-self.twod.im['SCI'].data+self.twod.im['CONTAM'].data+self.cont_model, vmin=-vmax, vmax=0.05*vmax, interpolation='nearest', aspect='auto', cmap=plt.cm.gray)
        ax.set_yticklabels([])
        ax.set_xticklabels(xint)
        ax.set_xticks(ax_int)
        #ax.set_ylabel('-Continuum')
        ax.text(0.98,0.96, '-Continuum', ha='right', va='top', color='red', size=7,  transform=ax.transAxes)
        
        ax.set_xlabel(r'$\lambda\ (\mu\mathrm{m})$')
        
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' + self.grism_id+'.%s.2D.%s' %(base, self.FIGURE_FORMAT))
        
    def line_free_template(self):
        """
        Assuming that the v1.1 line templates were used in the fit, generate the best-fit template
        from the version *without* emission lines.
        """

        #### Read in the line-free templates and scale them with "tnorm" to match those stored
        #### in the temp_sed file
        nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed1_nolines.dat', unpack=True)
        if self.skip_photometric:
            self.templam_nolines=nlx
            self.best_fit_nolines = nly #nlx*0.+1
            self.best_fit_nolines_flux = nly #nlx*0.+1
            return True
        
        #print self.eazy_coeffs['coeffs'].shape
        
        NTEMP = self.eazy_coeffs['coeffs'].shape[0]    
        noline_temps = np.zeros((nlx.shape[0], NTEMP))
        noline_temps[:,0] = nly/self.eazy_coeffs['tnorm'][0]
        for i in range(2,7):
            nlx, nly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.0_lines/eazy_v1.0_sed%d_nolines.dat' %(i), unpack=True)
            noline_temps[:,i-1] = nly/self.eazy_coeffs['tnorm'][i-1]

        #### The last v1.1 template was a BC03 model without lines, so just use it directly    
        i = 7
        lx, ly = np.loadtxt(unicorn.GRISM_HOME+'/templates/EAZY_v1.1_lines/eazy_v1.1_sed%d.dat' %(i), unpack=True)
        noline_temps[:,6] = np.interp(nlx, lx, ly)/self.eazy_coeffs['tnorm'][6]
        
        ### dusty old template
        if NTEMP == 8:
            lx, ly = np.loadtxt(unicorn.GRISM_HOME+'/templates/Ezgal/c09_del_8.6_z_0.019_chab_age09.40_av2.0.dat', unpack=True)
            noline_temps[:,7] = np.interp(nlx, lx, ly)/self.eazy_coeffs['tnorm'][7]
            
        #### Add an additional template here, need to make NTEMP=8 above
        # lx, ly = np.loadtxt(xxx path to mattia's new template xxx, unpack=True)
        # noline_temps[:,7] = np.interp(nlx, lx, ly)/self.eazy_coeffs['tnorm'][7]

        #### noline_temps has to have same number of templates as the 
        #### eazy coeffs file associated with the "zout" file specified in
        #### unicorn.analysis.read_catalogs.
        self.best_fit_nolines = np.dot(noline_temps, self.eazy_coeffs['coeffs'][:,self.ix])
        self.templam_nolines = nlx
        
        self.best_fit_nolines_flux = self.full_eazy_template(nlx, noline_temps, self.z_peak)
        
        return True
        
    def full_eazy_template(self, wavelength, templates, zi):
        """
        Generate the full EAZY template at a given redshift
        in physical units with IGM absorption
        """
        param = eazy.EazyParam(self.zout.filename.replace('zout','param'))
        flam_factor = 10**(-0.4*(param['PRIOR_ABZP']+48.6))*3.e18/1.e-17
        zi = self.eazy_tempfilt['zgrid'][self.eazy_coeffs['izbest'][self.ix]]

        ###### Full template SED, observed frame
        lambdaz = wavelength*(1+zi)
        temp_sed = np.dot(templates, self.eazy_coeffs['coeffs'][:,self.ix])
        temp_sed /= (1+zi)**2
        temp_sed *= (1/5500.)**2*flam_factor

        ###### IGM absorption
        lim1 = np.where(wavelength < 912)
        lim2 = np.where((wavelength >= 912) & (wavelength < 1026))
        lim3 = np.where((wavelength >= 1026) & (wavelength < 1216))

        if lim1[0].size > 0: temp_sed[lim1] *= 0.
        if lim2[0].size > 0: temp_sed[lim2] *= 1.-self.eazy_temp_sed['db'][self.eazy_coeffs['izbest'][self.ix]]
        if lim3[0].size > 0: temp_sed[lim3] *= 1.-self.eazy_temp_sed['da'][self.eazy_coeffs['izbest'][self.ix]]
        
        return temp_sed
    
    def new_fit_free_emlines(self, ztry=None, verbose=True, NTHREADS=1, NWALKERS=50, NSTEP=200, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100, save_chain=False, lrange=(0.9e4,2e4), tilt_order=1):
        """
        Fit the normalization and slope directly rather than with the
        fake red/blue templates.
        """
        import emcee
        import time
        if ztry is None:
            ztry = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        if verbose:
            print 'Fit lines: z=%.4f' %(ztry)
            
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        var[var == 0] = 1.e6
        var += (self.contam_ferror*self.twod.im['CONTAM'].data**2).flatten()

        flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        use = np.isfinite(flux)
        
        sh = self.twod.im['CONTAM'].data.shape
        wave2d = np.dot(np.ones(sh[0]).reshape((-1,1)), self.twod.im['WAVE'].data.reshape((1,-1)))
        #print use.sum()
        if self.grism_element == 'G800L':
            lrange = [0.6e4, 0.9e4]
        
        lrange = unicorn.reduce.grism_wlimit[self.grism_element][0:2]
        
        wave_ok = (wave2d.flatten() >= lrange[0]) & (wave2d.flatten() <= lrange[1])
        var[~wave_ok] = 1.e6
        
        #print use.sum()
        
        (use_lines, fancy), templates = self.get_emline_templates_new(ztry=ztry, line_width=line_width0)
        
        #### Initial guess for template normalizations
        twod_templates = np.zeros((1+len(use_lines), var.size))
        for i, temp in enumerate(templates):
            self.compute_model_function(lam_spec=self.linex*(1+ztry), flux_spec=temp/self.twod.total_flux)
            twod_templates[i,:] = self.twod.model.flatten()
                
        #
        amatrix = utils_c.prepare_nmf_amatrix(var[use], twod_templates[:,use])
        coeffs = utils_c.run_nmf(flux[use], var[use], twod_templates[:,use], amatrix, toler=1.e-5)
        
        try:
            amatrix = utils_c.prepare_nmf_amatrix(var[use], twod_templates[:,use])
            coeffs = utils_c.run_nmf(flux[use], var[use], twod_templates[:,use], amatrix, toler=1.e-5)
        except:
            threedhst.showMessage('Something went wrong with setting up the line fit', warn=True)
            return False
            
        flux_fit = np.dot(coeffs.reshape((1,-1)), templates).reshape((self.linex.shape))
        
        coeffs = np.maximum(coeffs, 1.e-5)
        #coeffs[0] = 1
        #### First parameter is slope, rest are template normalizations
        
        ## Log line fluxes 
        # init = np.append(np.log(coeffs), np.zeros(tilt_order))
        # #step_sig = np.append(np.ones(len(coeffs))*np.log(1.05), 0.1**(-np.arange(tilt_order)/2.))
        # step_sig = np.append(np.ones(len(coeffs))*np.log(1.05), 0.1*np.ones(tilt_order))
        # step_sig[0] = 0.1

        ## Linear line fluxes
        init = np.append(coeffs, np.zeros(tilt_order))
        #step_sig = np.append(np.ones(len(coeffs))*np.log(1.05), 0.1**(-np.arange(tilt_order)/2.))
        step_sig = np.append(np.ones(len(coeffs))*0.5, 0.1*np.ones(tilt_order))
        step_sig[0] = 0.1
                
        #threedhst.showMessage('%d' %(len(coeffs)), warn=True)
        #threedhst.showMessage(init.__str__(), warn=True)
        #threedhst.showMessage(step_sig.__str__(), warn=True)
        
        #step_sig = np.append(np.ones(len(coeffs))*np.log(1.05), 0.0*np.ones(tilt_order))
        
        #### Fit just emission lines, z fixed to z_grism
        obj_fun = unicorn.interlace_fit._objective_lineonly_new
        sh = self.twod.im['SCI'].data.shape
        wave_flatten = np.dot(np.ones(sh[0]).reshape(sh[0],1), self.oned.data.wave.reshape(1,sh[1])).flatten()
        obj_args = [flux, var, twod_templates, np.log(wave_flatten/1.4e4), 1]
        
        self.obj_args = obj_args
        self.obj_fun = obj_fun
        self.obj_init = init
        
        model = obj_fun(init, *obj_args).reshape(sh)
        
        #### Run the Markov chain
        obj_args = [flux, var, twod_templates, np.log(wave_flatten/1.4e4), 0]
        
        ndim = len(init)
        p0 = [(init+np.random.normal(size=len(init))*step_sig) for i in xrange(NWALKERS)]

        # NTHREADS, NSTEP = 1, 100
        if verbose:
            print 'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d)' %(NWALKERS, NSTEP)
            t0 = time.time()
            
        self.sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
        result = self.sampler.run_mcmc(p0, NSTEP)
        
        param_names = ['s0']
        param_names.extend(use_lines)
        param_names.extend(['s%d' %(i+1) for i in range(tilt_order)])
        
        self.chain = unicorn.interlace_fit.emceeChain(chain=self.sampler.chain, param_names = param_names)
        
        if verbose:
            print unicorn.noNewLine+'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d) -> Done (%.1f s).  Make plot...' %(NWALKERS, NSTEP, time.time()-t0)
        
        ##### Make the plot
        fig = unicorn.catalogs.plot_init(xs=5,aspect=1./1.65, left=0.15, right=0.02, bottom=0.09, top=0.02, NO_GUI=True)
        ax = fig.add_subplot(111)
        
        if self.grism_element == 'G141':
            wok = (self.oned_wave > 1.e4)
        
        if self.grism_element == 'G102':
            wok = (self.oned_wave > 0.74e4)
        
        if self.grism_element == 'G800L':
            wok = (self.oned_wave > 0.60e4)
        
        ax.plot(self.oned_wave[wok]/1.e4, self.oned.data.flux[wok]-self.oned.data.contam[wok], color='black', alpha=0.8)
        
        #### Line equivalent widths, combined error of line and continuum
        #### normalizations
        ntemp = templates.shape[0]
        eqw = np.ones(ntemp-1)
        eqw_err = np.ones(ntemp-1)

        ok = templates[0,:] > 0
        
        for i, line in enumerate(use_lines):
            #continuum_scale = np.exp(self.chain['s0']+self.chain['s1']*lx)
            lx = np.log(fancy[line][1][0]*(1+ztry)/1.4e4)
            continuum_scale = lx*0.
            for o in range(tilt_order+1):
                continuum_scale += self.chain['s%d' %(o)]*lx**o
            #
            continuum_scale = np.exp(continuum_scale)
            #
            ## log line flux
            #eqw_post = -np.trapz(-templates[i+1,ok]/templates[0,ok], self.linex[ok]*(1+ztry))*np.exp(self.chain[line])/continuum_scale
            ## linear line flux
            eqw_post = -np.trapz(-templates[i+1,ok]/templates[0,ok], self.linex[ok]*(1+ztry))*(self.chain[line])/continuum_scale
            eqw[i] = np.median(eqw_post[:,self.chain.nburn:])
            #eqw_err[i] = threedhst.utils.biweight(eqw_post)
            eqw_err[i] = (np.percentile(eqw_post[:,self.chain.nburn:], 84) - np.percentile(eqw_post[:,self.chain.nburn:], 16))/2.
            
        #### Show model drawn from posterior
        NSHOW = 100
        rand_walk = np.cast[int](np.random.rand(NSHOW)*self.chain.nwalkers)
        rand_step = np.cast[int](np.random.rand(NSHOW)*(self.chain.nstep-self.chain.nburn))+self.chain.nburn
        
        flux_models = []
        
        for i in range(NSHOW):
            param = self.chain.chain[rand_walk[i],rand_step[i],:]
            #print param
            obj_args = [flux, var, twod_templates, np.log(wave_flatten/1.4e4), 1]
            model_step = obj_fun(param, *obj_args).reshape(sh)
            wave_model, flux_model = self.twod.optimal_extract(model_step)
            aa = ax.plot(wave_model/1.e4, flux_model, color='red', alpha=0.05, zorder=1)    
            flux_models.append(flux_model)
        
        ##### Save the EMCEE chain and the draws made for the plot
        self.chain.save_fits(self.root+'.linefit.fits', verbose=False)
        lf = pyfits.open(self.root+'.linefit.fits', mode='update')
        lf[0].header.update('pointing', self.pointing)
        lf[0].header.update('id', self.twod.id)
        lf[0].header.update('z', self.z_max_spec)
        lf.append(pyfits.ImageHDU(data=np.array(flux_models), name='DRAW1D'))
        lf.append(pyfits.ImageHDU(data=wave_model, name='WAVE1D'))
        lf.flush()
        
        #ax.plot(self.oned_wave[wok]/1.e4, self.model_1D[wok], color='green', alpha=0.5, linewidth=2, zorder=10)
        
        #### Correction factors due to offset from photometry.  The line
        #### strenghts as fit are *before* the correction, so in observed
        #### frame need to be scaled
        lam = self.oned.data.wave
        #scale = np.exp(self.chain.stats['s0']['q50'] + self.chain.stats['s1']['q50']*np.log(lam/1.4e4))
        scale = np.exp(np.sum([self.chain.stats['s%d' %(o)]['q50']*np.log(lam/1.4e4)**o for o in range(tilt_order+1)], axis=0))
                
        corr_factor = eqw*0.
        for i, line in enumerate(use_lines):
            inv_corr = np.interp(fancy[line][1][0]*(1+ztry), lam, scale)
            if inv_corr == 0:
                corr_factor[i] = 1.
            else:
                corr_factor[i] = 1./inv_corr
            
        fp = open(self.grism_id+'.linefit.dat','w')
        fp.write('# line  flux error scale_to_photom EQW_obs EQW_obs_err \n# z=%.5f\n# flux: 10**-17 ergs / s / cm**2\n' %(ztry))
        
        fp.write('# [xxx] tilt correction to photometry: s0 = %.3f [ %.3f ], s1 = %.3f [ %.3f ]\n' %(self.chain.stats['s0']['q50'], self.chain.stats['s0']['width'], self.chain.stats['s1']['q50'], self.chain.stats['s1']['width']))
        print use_lines
        for i, line in enumerate(use_lines):
            #hist = plt.hist(chain[:,i+1], bins=50)
            #stats = threedhst.utils.biweight(chain[:,i+1], both=True)
            stats = self.chain.stats[self.chain.param_names[i+1]]
            ## Log line flux
            # flux_med = np.exp(stats['q50'])*(1+ztry)#/corr_factor[i]
            # flux_err = flux_med * stats['width'] #/corr_factor[i]
            ## Linear line flux
            flux_med = (stats['q50'])*(1+ztry)#/corr_factor[i]
            flux_err = stats['width']*(1+ztry) #/corr_factor[i]
            print flux_med, flux_err
            #median = np.median(chain[:,i+1])
            #range = np.percentile(chain[:,i+1], [15.8655,100-15.8655])
            fp.write('%4s  %6.2f  %.2f  %.3f %6.2f %6.2f\n' %(line, flux_med, flux_err, corr_factor[i], eqw[i], eqw_err[i]))
            ax.fill_between([0.03,0.43],np.array([0.95,0.95])-0.07*i, np.array([0.88, 0.88])-0.07*i, color='white', alpha=0.8, transform=ax.transAxes, zorder=-21)
            #
            if flux_med/flux_err > -1000:
                ax.text(0.05, 0.95-0.07*i, '%4s  %6.1f$\pm$%.1f  %4.1f  %6.1f\n' %(fancy[line][0], flux_med, flux_err, flux_med/flux_err, eqw[i]), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, zorder=-20, fontsize=7)
            else:
                ax.text(0.05, 0.95-0.07*i, '%4s  (%0.1f,%.1f)  %4.1f  %6.1f\n' %(fancy[line][0], (stats['q05'])*(1+ztry), (stats['q95'])*(1+ztry), flux_med/flux_err, eqw[i]), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, zorder=-20, fontsize=7)
                
        ax.fill_between([0.61,0.96],np.array([0.95,0.95]), np.array([0.88, 0.88]), color='white', alpha=0.8, transform=ax.transAxes, zorder=19)
        ax.text(0.95, 0.95, self.grism_id, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, zorder=20)
        ax.fill_between([0.8,0.96],np.array([0.85,0.85]), np.array([0.78, 0.78]), color='white', alpha=0.8, transform=ax.transAxes, zorder=19)
        ax.text(0.95, 0.85, r'$z=%.4f$' %(ztry), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, zorder=20)
        
        if self.grism_element == 'G141':
            ax.set_xlim(1.0,1.73)
            
        if self.grism_element == 'G102':
            ax.set_xlim(0.74, 1.18)
        #
        if self.grism_element == 'G800L':
            ax.set_xlim(0.58, 1.05)
        
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        ax.set_ylabel(r'Flux (e - / s)')
        ymax = self.model_1D[np.isfinite(self.model_1D)].max()
        ax.set_ylim(-0.1*ymax, 1.5*ymax)
        
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' +  self.grism_id+'.linefit.'+self.FIGURE_FORMAT)
        
        fp.close()
        
        return True
        
    def fit_free_emlines(self, ztry=None, verbose=True, NTHREADS=1, NWALKERS=50, NSTEP=200, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100):
        import emcee
        
        if ztry is None:
            ztry = self.zgrid1[self.full_prob1 == self.full_prob1.max()][0]
        
        if verbose:
            print 'Fit lines: z=%.4f' %(ztry)
            
        var = np.cast[float](self.twod.im['WHT'].data**2).flatten()
        var[var == 0] = 1.e6
        var += (self.contam_ferror*self.twod.im['CONTAM'].data**2).flatten()

        flux = np.cast[float](self.twod.im['SCI'].data-self.twod.im['CONTAM'].data).flatten()
        use = np.isfinite(flux)
        
        (use_lines, fancy), templates = self.get_emline_templates(ztry=ztry, line_width=line_width0)
        
        #### Initialize fit
        amatrix = utils_c.prepare_nmf_amatrix(var[use], templates[:,use])
        coeffs = utils_c.run_nmf(flux[use], var[use], templates[:,use], amatrix, toler=1.e-7)
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
        
        unicorn.catalogs.savefig(fig, self.OUTPUT_PATH + '/' +  self.grism_id+'.linefit.'+self.FIGURE_FORMAT)
        
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
    #
    def get_emline_templates_new(self, ztry=0., use_determined_lines=True, line_width=100):
        
        line_wavelengths = {} ; line_ratios = {}
        line_wavelengths['Ha'] = [6564.61]; line_ratios['Ha'] = [1.]
        line_wavelengths['Hb'] = [4862.68]; line_ratios['Hb'] = [1.]
        line_wavelengths['Hg'] = [4341.68]; line_ratios['Hg'] = [1.]
        line_wavelengths['Hd'] = [4102.892]; line_ratios['Hd'] = [1.]
        line_wavelengths['OIIIx'] = [4364.436]; line_ratios['OIIIx'] = [1.]
        line_wavelengths['OIII'] = [5008.240, 4960.295]; line_ratios['OIII'] = [2.98, 1]
        line_wavelengths['OII'] = [3729.875]; line_ratios['OII'] = [1]
        line_wavelengths['OI'] = [6302.046]; line_ratios['OI'] = [1]

        line_wavelengths['NeIII'] = [3869]; line_ratios['NeIII'] = [1.]
        #### Fix OIII / Hb
        #line_wavelengths['OIII'] = [5008.240, 4960.295, 4862.68]; line_ratios['OIII'] = [2.98, 1, 0.332]

        #line_wavelengths['OIII'] = [5008.240, 4960.295]; line_ratios['OIII'] = [4., 1]
        line_wavelengths['SII'] = [6718.29, 6732.67]; line_ratios['SII'] = [1, 1]
        line_wavelengths['SIII'] = [9068.6, 9530.6]; line_ratios['SIII'] = [1, 2.44]
        line_wavelengths['HeII'] = [4687.5]; line_ratios['HeII'] = [1.]
        line_wavelengths['HeI'] = [5877.2]; line_ratios['HeI'] = [1.]
        #### Test line
        #line_wavelengths['HeI'] = fakeLine; line_ratios['HeI'] = [1. for line in fakeLine]
        
        line_wavelengths['MgII'] = [2799.117]; line_ratios['MgII'] = [1.]
        line_wavelengths['CIV'] = [1549.480]; line_ratios['CIV'] = [1.]
        line_wavelengths['Lya'] = [1215.4]; line_ratios['Lya'] = [1.]
        
        fancy = {}
        fancy['Ha'] = (r'H$\alpha$', line_wavelengths['Ha'] )
        fancy['Hb'] = (r'H$\beta$', line_wavelengths['Hb'] )
        fancy['Hg'] = (r'H$\gamma$', line_wavelengths['Hg'] )
        fancy['Hd'] = (r'H$\delta$', line_wavelengths['Hd'] )
        fancy['OIII'] = ( 'O III', line_wavelengths['OIII'] )
        fancy['OI'] = ( 'O I', line_wavelengths['OI'] )
        fancy['OIIIx'] = ( 'O IIIx', line_wavelengths['OIIIx'] )
        fancy['NeIII'] = ( 'Ne III', line_wavelengths['NeIII'] )
        fancy['OII']  = ( 'O II', line_wavelengths['OII'] )
        fancy['SII']  = ( 'S II', line_wavelengths['SII'] ) 
        fancy['SIII'] = ( 'S III', line_wavelengths['SIII'] )
        fancy['HeI']  = ( 'He I', line_wavelengths['HeI'] ) 
        fancy['HeII']  = ( 'He II', line_wavelengths['HeII'] ) 
        fancy['MgII'] = ( 'Mg II', line_wavelengths['MgII'] )
        fancy['CIV']  = ( 'C IV', line_wavelengths['CIV'] ) 
        fancy['Lya'] = (r'Ly$\alpha$', line_wavelengths['Lya'] )
        
        fit_wavelengths = {'G102':[8000, 1.16e4], 'G141':[1.06e4, 1.68e4], 'G800L':[5420, 1.035e4]}
        
        if use_determined_lines & (self.use_lines is not None):
            use_lines = self.use_lines
        else:
            use_lines = []
            for line in line_wavelengths.keys():
                lam = line_wavelengths[line][0]
                if (lam*(1+ztry) > np.maximum(self.oned_wave.min(), fit_wavelengths[self.twod.grism_element][0])) & (lam*(1+ztry) < np.minimum(fit_wavelengths[self.twod.grism_element][1], self.oned_wave.max())):
                    use_lines.append(line)
        
        ### Make sure have resolution at the desired wavelenghts
        newx = self.linex*1.
        for line in use_lines:
            for l0 in line_wavelengths[line]:
                #l0 = line_wavelengths[line][0]
                if np.interp(l0, self.linex[1:], np.diff(self.linex))/l0*3.e5 > 100:
                    #print line, np.interp(l0, self.linex[1:], np.diff(self.linex))/l0*3.e5
                    xg = np.arange(-5,5,0.1)*100./3.e5*l0+l0
                    newx = np.append(newx, xg)
        #
        newx = np.unique(newx)
        newx = newx[np.argsort(newx)]
        self.liney = np.interp(newx, self.linex, self.liney)
        self.linex = newx
        
        
        templates = np.ones((1+len(use_lines), self.linex.shape[0]))
        templates[0,:] = np.interp(self.linex, self.templam_nolines, self.best_fit_nolines_flux)
        
        #line_width = 100 # km/s
        xline = self.linex
        sens_int = {}
        
        has_line_flux = np.ones(1+len(use_lines)) > 0
        
        for i,line in enumerate(use_lines):
            yline = self.linex*0.
            nmult = len(line_wavelengths[line])
            #print line, nmult, line_ratios[line]
            for j in range(nmult):
                lam = line_wavelengths[line][j]
                norm = line_ratios[line][j]
                dlam = lam*line_width/3.e5
                yline += 1./np.sqrt(np.pi*2*dlam**2)*np.exp(-(xline-lam)**2/2/dlam**2)*norm
            # 
            yline /= np.sum(line_ratios[line])
            templates[i+1,:] = yline
                
        self.use_lines = use_lines
        
        return (use_lines, fancy), templates
    
    def get_emline_templates(self, ztry=0., use_determined_lines=True, line_width=100):
        
        ### Continuum at z=ztry
        status = self.fit_zgrid(get_model_at_z=ztry, verbose=False)
        
        line_wavelengths = {} ; line_ratios = {}
        line_wavelengths['Ha'] = [6564.61]; line_ratios['Ha'] = [1.]
        line_wavelengths['Hb'] = [4862.68]; line_ratios['Hb'] = [1.]
        line_wavelengths['Hg'] = [4341.68]; line_ratios['Hg'] = [1.]
        line_wavelengths['Hd'] = [4102.892]; line_ratios['Hd'] = [1.]
        line_wavelengths['OIII'] = [5008.240, 4960.295]; line_ratios['OIII'] = [2.98, 1]
        line_wavelengths['OII'] = [3729.875]; line_ratios['OII'] = [1]
        line_wavelengths['SII'] = [6718.29, 6732.67]; line_ratios['SII'] = [1, 1]
        line_wavelengths['SIII'] = [9068.6, 9530.6]; line_ratios['SIII'] = [1, 2.44]
        line_wavelengths['HeI'] = [5877.2]; line_ratios['HeI'] = [1.]
        line_wavelengths['MgII'] = [2799.117]; line_ratios['MgII'] = [1.]
        line_wavelengths['CIV'] = [1549.480]; line_ratios['CIV'] = [1.]
        line_wavelengths['CIII'] = [1908.7]; line_ratios['CIII'] = [1.]
        line_wavelengths['Lya'] = [1215.4]; line_ratios['Lya'] = [1.]
        
        fancy = {}
        fancy['Ha'] = r'H$\alpha$'
        fancy['Hb'] = r'H$\beta$'
        fancy['Hg'] = r'H$\gamma$'
        fancy['Hd'] = r'H$\delta$'
        fancy['OIII'] =  'O III'
        fancy['OII']  =  'O II' 
        fancy['SII']  =  'S II' 
        fancy['SIII'] =  'S III'
        fancy['HeI']  =  'He I' 
        fancy['MgII'] =  'Mg II'
        fancy['CIV']  =  'C IV' 
        fancy['CIII']  =  'C III' 
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
            self.compute_model_function(xline*(1+ztry), yline)
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
        self.compute_model_function()
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
    
    def interpolate_direct_thumb(self, smooth_fraction=0.3):
        """
        Interpolate bad pixels in the direct image as they will mess 
        up the interlace model.
        
        If the fraction of negative pixels within the segmentation region is 
        greater than `smooth_fraction`, smooth the thumbnail with a gaussian
        
        """
        import scipy.ndimage as nd
        #from scipy.signal import convolve2d
        #kernel = np.ones((3,3))
        # wht = self.twod.im['DWHT'].data != 0
        # bad = (wht == 0) | (self.twod.im['DSCI'].data <= 0)
        # if wht.sum() < 10:
        #     return None
            
        # npix = convolve2d(wht, kernel, boundary='fill', fillvalue=0, mode='same')
        # sum = convolve2d(self.twod.im['DSCI'].data, kernel, boundary='fill', fillvalue=0, mode='same')
        # sumwht = convolve2d(self.twod.im['DWHT'].data, kernel, boundary='fill', fillvalue=0, mode='same')
        # fill_pix = bad & (npix > 0)
        # self.twod.im['DSCI'].data[fill_pix] = (sum/npix)[fill_pix]*1.
        # self.twod.im['DWHT'].data[fill_pix] = (sumwht/npix)[fill_pix]*1.
        # self.twod.thumb = self.twod.im['DSCI'].data*1.
        
        bad = (self.twod.im['DWHT'].data == 0) & (self.twod.im['DSCI'].data <= 0)
        thumb = self.twod.im['DSCI'].data*1
        thumb[bad] = 0
        max_filter = nd.maximum_filter(thumb, (3,3))
        thumb[bad] = max_filter[bad]
        
        thumb[thumb < 0] = 0
        
        seg_mask = (self.twod.seg == self.twod.id)
        neg_pix = self.twod.im['DSCI'].data < 0
        print seg_mask.sum()*smooth_fraction, (neg_pix & seg_mask).sum()
        if (seg_mask.sum()*smooth_fraction < (neg_pix & seg_mask).sum()):
            print 'Smooth direct thumbnail'
            thumb = nd.gaussian_filter(thumb, (1.4,1.4))
                    
        self.twod.thumb = thumb*1
        self.twod.im['DSCI'].data = thumb*1
        
        # self.twod.flux = self.twod.thumb * 10**(-0.4*(26.46+48.6))* 3.e18 / 1.3923e4**2 / 1.e-17
        self.twod.flux = self.twod.thumb * 10**(-0.4*(unicorn.reduce.ZPs[self.twod.im[0].header['FILTER']]+48.6))* 3.e18 / unicorn.reduce.PLAMs[self.twod.im[0].header['FILTER']]**2 / 1.e-17

        self.twod.total_flux = np.sum(self.twod.flux*seg_mask)
        self.twod.init_model()
            
    def show_lnprob(self):
        """
        Show ln probability from photometry and spectra
        """
        fig = unicorn.plotting.plot_init(square=True, aspect=0.6, xs=5, left=0.15, right=0.02, bottom=0.08, top=0.08, NO_GUI=False, use_tex=False, fontsize=10)
        
        spec_only = self.full_prob0-np.interp(self.zgrid0, self.phot_zgrid, self.phot_lnprob)
        spec_only -= spec_only.max()
        
        ax = fig.add_subplot(111)
        ax.plot(self.phot_zgrid, self.phot_lnprob, color='blue')
        ax.plot(self.zgrid0, spec_only, color='green', alpha=0.3)
        ax.plot(self.zgrid0, self.full_prob0, color='green')
        ax.plot(self.zgrid1, self.full_prob1, color='red')
        ax.set_xlabel(r'$z$')
        ax.set_ylabel(r'ln $p(z)$')
        ax.set_ylim(-210,20)
        ax.set_title(self.root)
        
        return fig
        
        #### Test floor on phot p(z)
        xfloor = self.phot_zgrid
        y0 = xfloor*0.+np.log(1.e-3/np.sum(np.diff(self.phot_zgrid)))
        y0[self.phot_lnprob < -200] = -200
        yfloor = np.maximum(self.phot_lnprob, y0)
        plt.plot(xfloor, yfloor)
#
def _objective_lineonly_new(params, observed, var, twod_templates, wave_flatten, get_model):
    ### The "minimum" function limits the exponent to acceptable float values
    flux_fit = twod_templates[0,:]*1.    
    #s0, s1 = params[1], params[0]
    #scale = np.exp(s0 + s1*wave_flatten)
    NT = twod_templates.shape[0]
    s = np.append(params[0], params[NT:])[::-1]
    scale = np.exp(polyval(s, wave_flatten))        
    flux_fit *= scale
    
    ## log line fluxes
    # flux_fit += np.dot(np.exp(np.minimum(params[1:NT],345)).reshape((1,-1)), twod_templates[1:,:]).flatten()
    
    ## linear line fluxes
    flux_fit += np.dot(params[1:NT].reshape((1,-1)), twod_templates[1:,:]).flatten()
    
    #print scale.min(), scale.max()
    
    prior = 0
    ### prior against negative fluxes
    prior = -0.5*np.sum(np.minimum(params[1:NT], 0)**2/1.**2)
    
    lnprob = -0.5*np.sum((observed-flux_fit)**2/var) + prior
    #print params, lnprob
    
    if get_model:
        return flux_fit
    else:
        if np.isfinite(lnprob):
            return lnprob
        else:
            print params, 'Nan!'
            return -1e10
            
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

def go_MCMC_new():
    import time
    import numpy as np
    import unicorn
    import emcee
    
    model = unicorn.reduce.process_GrismModel('GOODS-S-4')
    
    id=370
    
    if not os.path.exists('GOODS-S-4_%05d.2D.fits'):
        model.twod_spectrum(id=id)
    
    self = unicorn.interlace_fit.GrismSpectrumFit(root='GOODS-S-4_%05d' %(id))
    
    t0 = time.time()
    self.fit_in_steps(dzfirst=0.004, skip_second=True)
    t1 = time.time()
    print 'Simple grid: %.1f s' %(t1-t0)
    
    observed = self.twod.im['SCI'].data
    var = self.twod.im['WHT'].data**2
    var[var == 0] = 1.e4
    
    ha = np.abs(self.linex-6563.) < 100
    haflux = np.trapz(self.liney[ha], self.linex[ha])
    self.liney /= haflux
    
    self.continuum = np.interp(self.linex, self.templam_nolines, self.best_fit_nolines_flux)
    
    init = np.array([self.z_max_spec, 0, 0, 1.])
    step_sig = np.array([self.dz, 0.1, 0.1, 1])
    obj_fun = unicorn.interlace_fit._objective_z_simple
    obj_args = [observed, var, self, 0]
        
    ndim, nwalkers = len(init), 50
    p0 = [(init+np.random.normal(size=ndim)*step_sig) 
          for i in xrange(nwalkers)]
    
    NTHREADS, NSTEP = 1, 100
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, args = obj_args, 
                                    threads=NTHREADS)
    
    t0 = time.time()
    result = sampler.run_mcmc(p0, NSTEP)
    t1 = time.time()
    print 'Sampler: %.1f s' %(t1-t0)
    
    chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names=['z','c0','c1','line'])
    
    obj_args = [observed, var, self, 1]

    h = plt.hist(chain['z'][:,50:].flatten(), bins=200, range=(chain.stats['z']['q16']-0.05, chain.stats['z']['q84']+0.05), normed=True, alpha=0.5)
    linear = np.exp(self.full_prob1)/np.trapz(np.exp(self.full_prob1), self.zgrid1)
    p = plt.plot(self.zgrid1, linear)
    if self.zout.z_spec[self.ix] > 0:
        p = plt.plot(self.zout.z_spec[self.ix]*np.array([1,1]), [0,linear.max()], color='red', alpha=0.5, linewidth=4)
    
    plt.xlim(h[1].min(), h[1].max())
    
    final = chain.median

    plt.plot(self.oned.data.wave, self.oned.data.flux/self.oned.data.sensitivity)
    
    for i in range(nwalkers):
        walk = chain.chain[i,-1,:]
        lam, model1, model2 = obj_fun(walk, *obj_args)
        fx, fy = self.twod.optimal_extract(input=model2)
        p = plt.plot(fx, fy/self.oned.data.sensitivity, color='green', alpha=0.1)

    z6 = chain['z'] > 0.6
    z5 = chain['z'] < 0.6
    
    p6 = [np.median(chain['z'][z6]), np.median(chain['c0'][z6]), np.median(chain['c1'][z6]), np.median(chain['line'][z6])]
    p5 = [np.median(chain['z'][z5]), np.median(chain['c0'][z5]), np.median(chain['c1'][z5]), np.median(chain['line'][z5])]
    

    lam, model1, model2 = obj_fun(p6, *obj_args)
    plt.plot(lam, model1, color='red')

    lam, model1, model2 = obj_fun(p5, *obj_args)
    plt.plot(lam, model1, color='blue')
    
    plt.xlim(1.e4, 1.7e4)
    
def _objective_z_simple(params, observed, var, self, show):
    """
    Do all of the calculation on the oned spectrum and then
    send it through the compute_object_model at the end.
    
    1) Take the best-fit EAZY template *in flux units*
    2) Apply linear scaling, normalization and slope
    3) Eventually add em line templates
    
    params = [z, s0, s1]
    
    scale = np.exp(s0 + s1*np.log(lam(z)/1.4e4))
    model = eazy_model * scale
    
    self is a GrismSpectrumFit object
    
    """
    
    ztry = params[0]
    s0, s1 = params[1], params[2]
    scale = np.exp(s0 + s1*np.log(self.linex*(1+ztry)/1.4e4))
    model_spec = self.continuum * scale + self.liney*np.exp(params[3])
    
    self.compute_model_function(lam_spec=self.linex*(1+ztry), flux_spec=model_spec/self.twod.total_flux)
    
    if show == 1:
        return self.linex*(1+ztry), model_spec, self.twod.model
        
    ### photo-z prior
    #prior = np.interp(ztry, self.phot_zgrid, self.phot_lnprob)
    prior = np.interp(ztry, self.zgrid0, self.full_prob0)
    lnprob = -0.5*np.sum((observed-self.twod.model)**2/var)
    # print params, lnprob, prior
    
    return lnprob+prior
    
def more_fit_testing():
    """
    Look into the fit
    """
    model = unicorn.reduce.process_GrismModel('GOODS-S-4')
    id=296
    
    model.twod_spectrum(id=id)
    model.show_2d(savePNG=True, verbose=True)
    
    gris = unicorn.interlace_fit.GrismSpectrumFit(root='GOODS-S-4_%05d' %(id))
    print '\n\n'
    gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
    
    gris.fit_free_emlines(ztry=None, verbose=False, NTHREADS=1, NWALKERS=50, NSTEP=200, FIT_REDSHIFT=True, FIT_WIDTH=False, line_width0=100)
    
def go_MCMC_fit_OLD():
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

class emceeChain():     
    def __init__(self, chain=None, file=None, param_names=[],
                       burn_fraction=0.5):
        
        self.param_names = []
        
        if chain is not None:
            self.chain = chain
        
        if file is not None:
            if 'fits' in file.lower():
                self.load_fits(file=file)            
            else:
                self.load_chain(file=file)
                
        self.process_chain(param_names = param_names,
                           burn_fraction=burn_fraction)
        
    def process_chain(self, param_names=[], burn_fraction=0.5):
        """
        Define parameter names and get parameter statistics
        """
        self.nwalkers, self.nstep, self.nparam = self.chain.shape
                
        if param_names == []:            
            if self.param_names == []:
                for i in range(self.nparam):
                    param_names.append('a%d' %(i+1))
               
                self.param_names = param_names
        
        else:
            if len(param_names) != self.nparam:
                print 'param_names must have N=%d (or zero) entries' %(self.nparam)
                return False
                        
            self.param_names = param_names
                
        self.param_dict = {}
        for i in range(self.nparam):
            self.param_dict[self.param_names[i]] = i
        
        self.nburn = int(np.round(burn_fraction*self.nstep))
        self.stats = {}
        self.median = np.zeros(self.nparam)
        for param in self.param_names:
            pid = self.param_dict[param]
            self.stats[param] = self.get_stats(pid, burn=self.nburn)
            self.median[pid] = self.stats[param]['q50']
            
    def get_stats(self, pid, burn=0, raw=False):
        """
        Get percentile statistics for a parameter in the chain
        """
        if raw:
            pchain = pid*1.
        else:
            pchain = self.chain[:,burn:,pid].flatten()
        
        stats = {}
        stats['q05'] = np.percentile(pchain, 5)
        stats['q16'] = np.percentile(pchain, 16)
        stats['q50'] = np.percentile(pchain, 50)
        stats['q84'] = np.percentile(pchain, 84)
        stats['q95'] = np.percentile(pchain, 95)
        stats['mean'] = np.mean(pchain)
        stats['std'] = np.std(pchain)
        stats['width'] = (stats['q84']-stats['q16'])/2.
        return stats
        
    def show_chain(self, param='a1', chain=None, alpha=0.15, color='blue', scale=1, diff=0, ax = None, add_labels=True, hist=False, autoscale=True, *args, **kwargs):
        """
        Make a plot of the chain for a given parameter.
        
        For plotting, multiply the parameter by `scale` and subtract `diff`.
        
        """
        if chain is None:
            pid = self.param_dict[param]
            chain = self.chain[:,:,pid]
        
        if ax is not None:
            plotter = ax
            xlabel = ax.set_xlabel
            ylabel = ax.set_ylabel
            ylim = ax.set_ylim
        else:
            plotter = plt
            xlabel = plt.xlabel
            ylabel = plt.ylabel
            ylim = plt.ylim
            
        if hist:
            h = plotter.hist(chain[:,self.nburn:].flatten(), alpha=alpha, color=color, *args, **kwargs)
            if add_labels:
                ylabel('N')
                xlabel(param)
        else:
            for i in range(self.nwalkers):
                p = plotter.plot(chain[i,:]*scale-diff, alpha=alpha, color=color, *args, **kwargs)
            if add_labels:
                xlabel('Step')
                ylabel(param)
            #
            if autoscale:
                ylim(self.stats[param]['q50'] + np.array([-8,8])*self.stats[param]['width'])
            
    def save_chain(self, file='emcee_chain.pkl', verbose=True):
        """
        Save the chain to a Pkl file
        """
        import cPickle as pickle
        
        fp = open(file,'wb')
        pickle.dump(self.nwalkers, fp)
        pickle.dump(self.nstep, fp)
        pickle.dump(self.nparam, fp)
        pickle.dump(self.param_names, fp)
        pickle.dump(self.chain, fp)
        fp.close()
        
        if verbose:
            print 'Wrote %s.' %(file)
        
    def load_chain(self, file='emcee_chain.pkl'):
        """
        Read the chain from the pickle file
        """
        import cPickle as pickle
        
        fp = open(file, 'rb')
        self.nwalkers = pickle.load(fp)
        self.nstep = pickle.load(fp)
        self.nparam = pickle.load(fp)
        self.param_names = pickle.load(fp)
        self.chain = pickle.load(fp)
        fp.close()
    
    def save_fits(self, file='emcee_chain.fits', verbose=True):
        """
        Make a FITS file of an EMCEE chain
        """
        header = pyfits.Header()
        header.update('NWALKERS', self.nwalkers)
        header.update('NSTEP', self.nstep)
        header.update('NPARAM', self.nparam)
        
        hdu = [pyfits.PrimaryHDU(header=header)]
        
        for param in self.param_names:
            header.update('PARAM', param)
            hdu.append(pyfits.ImageHDU(data=self.__getitem__(param), header=header, name=param))
        
        hduList = pyfits.HDUList(hdu)
        hduList.writeto(file, clobber=True, output_verify='silentfix')
    
        if verbose:
            print 'Wrote %s.' %(file)
    
    def load_fits(self, file='emcee_chain.fits'):
        im = pyfits.open(file)
        self.nwalkers = im[0].header['NWALKERS']
        self.nstep = im[0].header['NSTEP']
        self.nparam = im[0].header['NPARAM']
        self.param_names = []
        self.chain = np.ones((self.nwalkers, self.nstep, self.nparam))
        for i in range(self.nparam):
            self.param_names.append(im[i+1].header['PARAM'])
            self.chain[:,:,i] = im[i+1].data
        
        im.close()
    
    def parameter_correlations(self, size=8, shrink=5, show=None, file=None):
        if show is None:
            show = self.param_names
        
        NP = len(show)
        fig = unicorn.plotting.plot_init(square=True, aspect=1, xs=size, left=0.05, right=0.01, bottom=0.01, top=0.01, NO_GUI=False, use_tex=False, fontsize=7)
        fig.subplots_adjust(wspace=0.0,hspace=0.0)
        
        counter = 0
        for i in range(NP):
            for j in range(NP):
                counter = counter + 1
                ax = fig.add_subplot(NP, NP, counter)
                a = ax.plot(self[show[i]][:,self.nburn::shrink].flatten(), self[show[j]][:,self.nburn::shrink].flatten(), alpha=0.03, color='black', linestyle='None', marker=',')
                a = ax.set_xlim(self.stats[show[i]]['q50']-3*self.stats[show[i]]['std'], self.stats[show[i]]['q50']+3*self.stats[show[i]]['std'])
                a = ax.set_ylim(self.stats[show[j]]['q50']-3*self.stats[show[j]]['std'], self.stats[show[j]]['q50']+3*self.stats[show[j]]['std'])
                if i == j:
                    a = ax.text(0.5, 0.92, show[i], fontsize=8, color='red', horizontalalignment='center', verticalalignment='top', transform=ax.transAxes)

        if file is not None:
            fig.savefig(file)
    
    def draw_random(self, N=10):
        """
        Draw random sets of parameters from the chain
        """
        #ok_walk = self.sampler.acceptance_fraction > min_acceptance
        iwalk = np.cast[int](np.random.rand(N)*self.nwalkers)
        istep = self.nburn + np.cast[int](np.random.rand(N)*(self.nstep-self.nburn))
        draw = self.chain[iwalk, istep, :]
        return draw
        
    def __getitem__(self, param):
        pid = self.param_dict[param]
        return self.chain[:,:,pid]
        
def objective_twod_OLD(params, return_model=False):
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
def test_threading(nproc=4, N=5):
    """
    Testing parallel threading, nothing promising because it seems like the threading/multiprocessing bookkeeping overhead is significant compared to the "model.compute_model()" execution time.
    """
    import stsci.convolve
 
    from multiprocessing import Pool, Process
    from threading import Thread
    
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

def check_lineflux():
    import pysynphot as S
    import scipy.ndimage as nd
    
    os.chdir('%s/3DHST_VariableBackgrounds/Specz' %os.getenv('THREEDHST'))
    
    gris = unicorn.interlace_fit.GrismSpectrumFit('AEGIS-10_19259')
    gris.new_fit_free_emlines()
    
    
    clean = gris.twod.im['SCI'].data - gris.twod.im['CONTAM'].data - gris.twod.im['MODEL'].data

    kernel = np.ones((4,4))
    sm_sci = nd.convolve(clean, kernel)
    sm_var = nd.convolve(gris.twod.im['WHT'].data**2, kernel)
    
    l0 = 6563.*(1+gris.z_max_spec)
    line = S.GaussianSource(1.e-17, l0, 120./3.e5*l0)
    line = S.GaussianSource(24.5e-17, l0, 80.)
    contin = S.FlatSpectrum(22.68, fluxunits='ABMag')
    grism = S.ObsBandpass('wfc3,ir,g141')
    grism_obs = S.Observation(line+contin, grism)
    
    pix_scl = gris.twod.growx*gris.twod.growy*1.
    peak = sm_sci == sm_sci.max()
    line_flux = sm_sci[peak][0]/(grism_obs.countrate()/pix_scl)
    line_err = np.sqrt(sm_var[peak][0])/(grism_obs.countrate()/pix_scl)
    
    ensquared = 0.554
    gris.twod.compute_model(grism_obs.wave, grism_obs.flux*line_flux/ensquared/1.e-17/gris.twod.total_flux)
    
# class multiModel(Process):
#     def __init__(self, spec_2d):
#         Process.__init__(self)
#         self.spec_2d = spec_2d
#         self.status = False
#     #
#     def run(self):
#         self.spec_2d.compute_model()
#         self.model = self.spec_2d.model.copy()
#         self.status = True
#     
# class threadedModel(Thread):
#     def __init__(self, spec_2d):
#         Thread.__init__(self)
#         self.spec_2d = spec_2d
#         self.status = False
#     #
#     def run(self):
#         self.spec_2d.compute_model()
#         self.model = self.spec_2d.model.copy()
#         self.status = True

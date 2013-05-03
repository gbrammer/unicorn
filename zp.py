"""
Try monte-carloing EAZY ZP corrections
"""
import threedhst
from threedhst import catIO
from threedhst import eazyPy as eazy
import unicorn
from unicorn import utils_c

import matplotlib.pyplot as plt
import numpy as np

class IterZeropoints():
    def __init__(self, root='hudf09', PATH='OUTPUT'):
    
        if not PATH.endswith('/'):
            PATH += '/'

        ##### Read the param file
        self.param = eazy.EazyParam(PARAM_FILE=PATH+'/'+root+'.param')
        self.filter_index = {}
        for i in range(self.param.NFILT):
            self.filter_index[self.param.filters[i].fnumber] = i
            
        ##### Read template fluxes and coefficients
        self.tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH,CACHE_FILE = 'Same')
        self.lc = self.tempfilt['lc']
                
        ##### Photometry interpolator
        self.phot = eazy.TemplateInterpolator(bands=None, MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH, CACHE_FILE='Same', zout=None, f_lambda=False)
    
        ##### Template Error interpolator
        self.te = eazy.TemplateError()
        
        #### Use Z_SPEC objects
        self.zout = catIO.Readfile('%s/%s.zout' %(PATH, root))
        
        self.initialize_arrays()
        
    def initialize_arrays(self):
        
        ok = self.zout.z_spec > 0
        self.NOBJ = ok.sum()
        
        self.fnu = np.zeros((self.NOBJ, self.phot.NFILT))
        self.efnu = np.zeros((self.NOBJ, self.phot.NFILT))
        self.full_var = np.zeros((self.NOBJ, self.phot.NFILT))
        self.templates = np.zeros((self.NOBJ, self.phot.NTEMP, self.phot.NFILT))
        self.z = self.zout.z_spec[ok]
        self.coeffs = np.ones((self.NOBJ, self.phot.NTEMP))
        self.model = np.zeros((self.NOBJ, self.phot.NFILT))
        self.amatrix = np.zeros((self.NOBJ, self.phot.NTEMP, self.phot.NTEMP))
        
        for i, ix in enumerate(np.arange(len(ok))[ok]):
            self.templates[i,:,:] = self.phot.interpolate_photometry(self.z[i]).reshape(self.phot.NFILT, self.phot.NTEMP).T
            self.fnu[i,:] = self.tempfilt['fnu'][:,ix]
            self.efnu[i,:] = self.tempfilt['efnu'][:,ix]
            t_err = self.te.interpolate(self.lc, self.z[i])
            self.full_var[i,:] = self.efnu[i,:]**2+(t_err*0.5*self.fnu[i,:])**2
            #
            
        self.mask = (self.efnu > 0) & (self.fnu > -99)

        for i in range(self.NOBJ):
            m_i = self.mask[i,:]
            self.amatrix[i,:,:] = utils_c.prepare_nmf_amatrix(self.full_var[i,m_i], self.templates[i,:,m_i].T)
        
        self.zeropoints = np.ones(self.phot.NFILT)
        self.first = 1
        self.prior_sigma = 0.01
    
    def go_fit_UDS_v3(self):
        import unicorn.catalogs2 as cat2
        
        #### Make catalog with more spec-zs
        zsp = cat2.SpeczCatalog()
        cat = catIO.Readfile('uds_3dhst.v3.0.nzpcat', force_lowercase=False)
        matcher = catIO.CoordinateMatcher(cat2.cat)
        
        cat2.read_catalogs('UDS')
        dr, idx = matcher.match_list(cat.ra, cat.dec)

        zsp3d = (cat2.dq.f_cover > 0.8) & (cat2.dq.f_negative < 0.1) & (cat2.cat.m140 < 22.5) & ((cat2.lines['HALPHA_FLUX']/cat2.lines['HALPHA_FLUX_ERR'] > 3) | (cat2.lines['OIII_FLUX']/cat2.lines['OIII_FLUX_ERR'] > 3))
        ok = (dr < 0.5) & zsp3d[idx]
        self.zout.z_spec[ok] = cat2.zfit.z_max_spec[idx][ok]*1
        
        dr, idx = zsp.match_list(cat.ra, cat.dec)
        ok = (dr < 0.5) & (zsp.zspec[idx] > 0)
        #self.zout.z_spec = self.zout.z_spec*0.-99        
        self.zout.z_spec[ok] = zsp.zspec[idx][ok]*1
                
        cat.z_spec = self.zout.z_spec
        cat.write_text(filename='uds_3dhst.v3.0.nzpcat.z_3d', select=cat.z_spec > 0)
        
        #### Initial perturbed bands
        
        self.initialize_arrays()
        
        self.first = np.ones(self.phot.NFILT)
        self.prior_sigma = np.ones(self.phot.NFILT)*0.02
        
        ### u, subaru bands are a mess
        initial_adjustments = [[88, 1.3, 0.1], [122, 0.95, 0.1], [79, 0.84, 0.05], [123, 0.8, 0.1], [124, 0.81, 0.1], [125, 0.9, 0.1], [236, 1.1, 0.2], [239, 1.1, 0.2]]
        for adj in initial_adjustments:
            ix = self.filter_index[adj[0]]
            self.first[ix] = adj[1]
            self.prior_sigma[ix] = adj[2]
        #
        self.fit_emcee(NWALKERS=50, NSTEP=100, verbose=True)
        
        for i in range(self.phot.NFILT):
            FID = 'F%d' %(self.param.filters[i].fnumber)
            norm_chain = self.chain[FID]/self.chain['F265']
            print '%s  %.3f' %(FID, self.chain.get_stats(norm_chain[:, self.chain.nburn:], raw=True)['q50'])
        
        
    def go_fit(self):
        import unicorn.catalogs2 as cat2
        
        #### Get new spec-zs
        zsp = cat2.SpeczCatalog()
        cat = catIO.Readfile('hudfcat_mar3012.cat')
        kmag = 25-2.5*np.log10(cat.f_f160w_tot)
        ok = (kmag < 23) & (self.zout.z_spec < 0)
        self.zout.z_spec[ok] = self.zout.z_peak[ok]
        
        dr, idx = zsp.match_list(cat.ra, cat.dec)
        ok = (dr < 0.5) & (zsp.zspec[idx] > 0)
        #self.zout.z_spec = self.zout.z_spec*0.-99        
        self.zout.z_spec[ok] = zsp.zspec[idx][ok]*1
        
        self.first = np.ones(self.phot.NFILT)
        self.prior_sigma = np.ones(self.phot.NFILT)*0.02

        self.tempfilt['efnu'][0,:] *= 1000
        self.tempfilt['efnu'][-2,:] *= 1000
        
        self.first[0] = 1.20
        self.first[-2] = 0.8
        self.first[-1] = 0.93
        self.prior_sigma[0] = 0.1
        self.prior_sigma[-2:] = 0.1
        
        self.initialize_arrays()
        
        self.first = np.ones(self.phot.NFILT)
        self.prior_sigma = np.ones(self.phot.NFILT)*0.02
        
        self.tempfilt['efnu'][0,:] = -1 #*= 1000
        self.tempfilt['efnu'][-2,:] = -1 #*= 1000
        
        self.initialize_arrays()
        
        self.fit_emcee(NWALKERS=50, NSTEP=100, verbose=True)
        
        colors = range(self.phot.NFILT)
        for ci, i in enumerate(np.argsort(self.lc)):
            colors[i] = threedhst.utils.color_table((ci+1.)/self.tempfilt['NFILT']*250, table='rainbow.rgb')
        
        for i in range(self.phot.NFILT):
            FID = 'F%d' %(self.param.filters[i].fnumber)
            norm_chain = self.chain[FID]/self.chain['F205']
            self.chain.show_chain(chain=norm_chain, color=colors[i], label=FID)
        
        #
        for i in range(self.phot.NFILT):
            FID = 'F%d' %(self.param.filters[i].fnumber)
            norm_chain = self.chain[FID]/self.chain['F265']
            print '%s  %.3f' %(FID, self.chain.get_stats(norm_chain[:, self.chain.nburn:], raw=True)['q50'])
            
    def fit_emcee(self, NWALKERS=50, NSTEP=200, verbose=True):
        import time
        import emcee
        
        p0 = np.array([np.random.normal(size=self.phot.NFILT)*0.05+self.first for i in range(NWALKERS)])
        
        obj_fun = self._objective_loop_over_objects
        obj_args = [self] ## 
        
        NTHREADS=1
        ndim = self.phot.NFILT
        
        self.verbose=verbose
        self.emcee_counter = 0
        
        if verbose:
            print 'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d)' %(NWALKERS, NSTEP)
            t0 = time.time()
        
        self.sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
        result = self.sampler.run_mcmc(p0, NSTEP)
        
        param_names = ['F%d' %(f.fnumber) for f in self.param.filters] 
        self.chain = unicorn.interlace_fit.emceeChain(chain=self.sampler.chain, param_names = param_names)
        
        if verbose:
            print unicorn.noNewLine+'emcee MCMC fit: (NWALKERS x NSTEPS) = (%d x %d) -> Done (%.1f s).  Make plot...' %(NWALKERS, NSTEP, time.time()-t0)
        
    @staticmethod    
    def _objective_loop_over_objects(zeropoints, self):
        """
        Apply zeropoints to observed photometry and compute
        template normalizations
        """
        tweaked = self.fnu*zeropoints
        model = np.zeros((self.NOBJ, self.phot.NFILT))
        for i in range(self.NOBJ):
            m_i = self.mask[i,:]
            if m_i.sum() == 0:
                continue
            #
            t_i = self.templates[i,:,:]
            ### Warning: updating coefficients in place
            self.coeffs[i,:]  = utils_c.run_nmf(tweaked[i,m_i], self.full_var[i,m_i], t_i[:,m_i], self.amatrix[i,:,:], toler=1.e-4, MAXITER=1e6, init_coeffs=self.coeffs[i,:] )
            model[i,:] = np.dot(self.coeffs[i,:], t_i)
            
        lnprob = -0.5*np.sum(((tweaked-model)**2/self.full_var)[self.mask])
        prior = -np.sum((zeropoints-self.first)**2/2/self.prior_sigma**2)
        if self.verbose:
            self.emcee_counter += 1
            st = ' '.join(['%.3f' %(z) for z in zeropoints])
            print '%s - %f  %f [%d]' %(st, lnprob, prior, self.emcee_counter)
            
        return lnprob + prior

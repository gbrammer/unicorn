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
        self.t_err = np.zeros((self.NOBJ, self.phot.NFILT))
        self.templates = np.zeros((self.NOBJ, self.phot.NTEMP, self.phot.NFILT))
        self.z = self.zout.z_spec[ok]
        self.coeffs = np.ones((self.NOBJ, self.phot.NTEMP))
        self.model = np.zeros((self.NOBJ, self.phot.NFILT))
        self.amatrix = np.zeros((self.NOBJ, self.phot.NTEMP, self.phot.NTEMP))
        
        for i, ix in enumerate(np.arange(len(ok))[ok]):
            self.templates[i,:,:] = self.phot.interpolate_photometry(self.z[i]).reshape(self.phot.NFILT, self.phot.NTEMP).T
            self.fnu[i,:] = self.tempfilt['fnu'][:,ix]
            self.efnu[i,:] = self.tempfilt['efnu'][:,ix]
            self.t_err[i,:] = self.te.interpolate(self.lc, self.z[i])
            self.full_var[i,:] = self.efnu[i,:]**2+(self.t_err[i,:]*0.5*self.fnu[i,:])**2
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
        #initial_adjustments = [[88, 1.3, 0.1], [122, 0.95, 0.1], [79, 0.84, 0.05], [123, 0.8, 0.1], [124, 0.81, 0.1], [125, 0.9, 0.1], [236, 1.1, 0.2], [239, 1.1, 0.2]]
        initial_adjustments = [[88, 1.3, 0.1], [18, 1.0, 0.1], [19, 1.0, 0.2], [20, 1.0, 0.3], [21, 1.0, 0.3], [122, 0.94, 0.05]]
            
        for adj in initial_adjustments:
            ix = self.filter_index[adj[0]]
            self.first[ix] = adj[1]
            self.prior_sigma[ix] = adj[2]
        #
        ### Blow up errors on some filters so they don't contribute to the fit
        errors = [[88, 100], [265, 5], [18, 100], [19, 100], [20, 100], [21, 100]]
        for adj in initial_adjustments:
            ix = self.filter_index[adj[0]]
            self.efnu[:,ix] *= adj[1]
            self.full_var[:,ix] = self.efnu[:,ix]**2+(self.t_err[:,ix]*0.5*self.fnu[:,ix])**2
            
        self.fit_emcee(NWALKERS=50, NSTEP=100, verbose=True)
        
        for i in range(self.phot.NFILT):
            FID = 'F%d' %(self.param.filters[i].fnumber)
            norm_chain = self.chain[FID]/self.chain['F204']
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

def stars_zp():
    
    c = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.test.nzpcat', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.fixconv', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.HAWKI', force_lowercase=False)
    
    c = catIO.Readfile('../Catalogs/cosmos_3dhst.v4.0.test2.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/aegis_3dhst.v4.0.test2.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.test2.nzpcat', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat.fixconv', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat.HDFN', force_lowercase=False)
    c.flux_radius_F160W = c.flux_radius; c.f_F814W = c.f_F775W
    
    c = catIO.Readfile('../Catalogs/cosmos_3dhst.v4.0.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/aegis_3dhst.v4.0.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.GABODS', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.FULL', force_lowercase=False)
    c = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.fixconv', force_lowercase=False)
    
    hmag = 25-2.5*np.log10(c.f_F160W)
    plt.scatter(hmag, c.flux_radius_F160W, alpha=0.5, color='black')
    plt.xlim(11,25); plt.ylim(0,10)
    
    xstar, ystar = [16,20,24.8], [5.8,3., 2.4]
    yint = np.interp(hmag, xstar, ystar)
    is_star = (hmag > xstar[0]) & (hmag < xstar[-1]) & (c.flux_radius_F160W < yint)    
    plt.scatter(hmag[is_star], c.flux_radius_F160W[is_star], alpha=0.5, color='orange')
    
    #### Only use stars well-spanned by the Kurucz models
    if c.filename == '../Catalogs/goodss_3dhst.v4.0.nzpcat.fixconv':
        iH = -2.5*np.log10(c.f_F814Wcand/c.f_F160W)
    else:
        iH = -2.5*np.log10(c.f_F814W/c.f_F160W)
    
    #iH = -2.5*np.log10(c.f_F850LP/c.f_F160W)
    sp_type = (iH > -0.2) & (iH < 0.6)
    plt.scatter(hmag[is_star & sp_type], c.flux_radius_F160W[is_star & sp_type], alpha=0.5, color='red')
    
    c.write_text(filename=c.filename+'.star', select=is_star & sp_type)

    c.write_text(filename=c.filename+'.bright', select=~is_star & (hmag < 24) & (hmag > 17))
    
    plt.scatter(hmag[is_star], (c.f_F814W/c.e_F814W)[is_star])
    plt.scatter(hmag[is_star & sp_type], (c.f_F814W/c.e_F814W)[is_star & sp_type], color='red')
    ok = hmag < xstar[1]
    
    vj = -2.5*np.log10(c.f_F814W/c.f_F125W)
    jh = -2.5*np.log10(c.f_F125W/c.f_F160W)

    vj = -2.5*np.log10(c.f_B/c.f_z)
    jh = -2.5*np.log10(c.f_z/c.f_K)
    
    plt.scatter(vj[ok & ~is_star], jh[ok & ~is_star], color='black', alpha=0.1)
    plt.scatter(vj[is_star], jh[is_star], color='red', alpha=0.5)
    
    tempfilt, coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='pickles',                                                 OUTPUT_DIRECTORY='./OUTPUT',                                                 CACHE_FILE = 'Same')
        
    
#
def check_zeropoints(root='kurucz3', PATH='./OUTPUT/', savefig=None, adjust_zeropoints='zphot.zeropoint', fix_filter=None):
    """
    Plot the EAZY fit residuals to evaluate zeropoint updates
    """
    import threedhst
    import threedhst.eazyPy as eazy
    import os
    
    if not PATH.endswith('/'):
        PATH += '/'
    
    ##### Read the param file
    param = eazy.EazyParam('%s%s.param' %(PATH, root))
    
    ##### Read template fluxes and coefficients
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH,CACHE_FILE = 'Same')
    
    if coeffs['izbest'].max() == 0:
        STAR_FIT = True
    else:
        STAR_FIT = False
        
    param = eazy.EazyParam(PARAM_FILE=PATH+'/'+root+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
        
    zpfile = PATH+'/'+root+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])
        
    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1), np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))
    
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    obs_sed = np.zeros((tempfilt['NFILT'], tempfilt['NOBJ']), dtype=np.float)
    for i in xrange(tempfilt['NOBJ']):
        obs_sed[:,i] = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][i]], coeffs['coeffs'][:,i])
    
    zi = tempfilt['zgrid'][coeffs['izbest']]
    lc = tempfilt['lc']
    offsets = lc*0.
    
    resid = (obs_sed-tempfilt['fnu']) / obs_sed + 1
    signoise = tempfilt['fnu']/np.sqrt(tempfilt['efnu']**2+(0.01*tempfilt['fnu'])**2)
    
    #### Plot colors
    colors = range(tempfilt['NFILT'])
    for ci, i in enumerate(np.argsort(lc)):
        colors[i] = threedhst.utils.color_table((ci+1.)/tempfilt['NFILT']*250, table='rainbow.rgb')

    sigma_resid = (obs_sed-tempfilt['fnu'])/np.sqrt(tempfilt['efnu']**2 + (0.015*tempfilt['fnu'])**2)
    mag = 25-2.5*np.log10(tempfilt['fnu'])
    
    NFILT, NOBJ = mag.shape
    
    NY = np.round(np.sqrt(NFILT))
    NX = np.ceil(NFILT/NY)
    
    fig = unicorn.plotting.plot_init(square=True, xs=10, aspect=NY/NX, left=0.1, bottom=0.1)
    for i in range(NFILT):
        ax = fig.add_subplot(int(NY*100+NX*10+1+i))
        #i=0
        ax.scatter(mag[i,:], sigma_resid[i,:], alpha=0.8, color='black')
        ok = np.isfinite(mag[i,:]) 
        xm, ym, ys, nn = threedhst.utils.runmed(mag[i,ok], sigma_resid[i,ok], NBIN=ok.sum()/30)
        ax.fill_between(xm, ym+ys,ym-ys, color='0.5', alpha=0.5)
        ax.plot(xm, ym, color='0.1', alpha=0.5)
        ax.plot(xm, xm*0+1, color='red', alpha=0.8)
        ax.plot(xm, xm*0-1, color='red', alpha=0.8)
        ax.set_ylim(-5,5)
        ax.set_xlim(16,24)
        ax.text(0.5,0.95,'F%d  %s' %(param.filters[i].fnumber, os.path.basename(param.filters[i].name)), ha='center', va='top', transform=ax.transAxes)

def get_v40_zeropoints():
    import unicorn.zp
    
    ##### Fit HST first, then fix
    cat = catIO.Readfile('../Catalogs/cosmos_3dhst.v4.0.test2.nzpcat.bright')
    cat = catIO.Readfile('../Catalogs/cosmos_3dhst.v4.0.nzpcat.bright', force_lowercase=False)
    keep = cat.z_spec > 0
    cat.write_text(cat.filename+'.zspec', select=keep)
    
    os.system('eazy -p zphot.param.cosmos.step -t zphot.translate.cosmos.step')
    unicorn.zp.loop_zeropoints(root='cosmos.step', tfile='zphot.translate.cosmos.step',  zfile='zphot.zeropoint.cosmos.step', fix_zspec=True, fix_filter={205:1}, ignore_initial=['f_F814W','f_F140W'], init_filter={239: 0.98}, toler=0.01)
    
    f0, zp0 = np.loadtxt('zphot.zeropoint.cosmos.step', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
    
    #### Fix IRAC
    for id in np.arange(18,22):
        fix_hst[id] = 1.
    #
    # fix_hst = {203: 1.0264, 204: 1.0414, 205: 1.0, 236: 1.0302, 239: 1.0005}
    os.system('eazy -p zphot.param.cosmos -t zphot.translate.cosmos')
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter=fix_hst, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'], toler=0.01) #, ignore_all=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'])

    #### Initial guess to get the very discrepant bands close
    init_filter = {88:1.18, 81: 0.85, 82:1.2, 79:0.8, 190: 0.8, 192:0.7}
    
    os.system('eazy -p zphot.param.cosmos -t zphot.translate.cosmos')
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter=fix_hst, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709', 'f_R', 'f_I', 'f_UVISTA_Y', 'f_J1', 'f_J2', 'f_J3', 'f_H1', 'f_H2'], toler=0.005, init_filter=init_filter) #, ignore_all=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'])
    os.system('mv cosmos_*png float_all_cosmos/')
    os.system('mv zphot.zeropoint.cosmos float_all_cosmos/')
    #### use tweaked template at beginning
    init_filter = {88:1.18, 81: 0.85, 82:1.2, 79:0.8, 190: 0.8, 192:0.7}
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter=fix_hst, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709', 'f_R', 'f_I', 'f_UVISTA_Y', 'f_J1', 'f_J2', 'f_J3', 'f_H1', 'f_H2'], toler=0.005, init_filter=init_filter) #, ignore_all=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'])
    
    #### Use full set of zeropoints as initial guess and allow HST bands to float
    f0, zp0 = np.loadtxt('zphot.zeropoint.cosmos', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
    #
    os.system('mv zphot.zeropoint.cosmos float_all_cosmos/')
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter={205:1.}, toler=0.003, init_filter=fix_hst) #, ignore_all=['f_Zp', 'f_Ip', 
    
    ### GOODS-N
    cat = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.test2.nzpcat.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat.fixconv.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.nzpcat.HDFN.bright', force_lowercase=False)
    keep = cat.z_spec > 0
    cat.write_text(cat.filename+'.zspec', select=keep)
    
    os.system('eazy -p zphot.param.goodsn.step -t zphot.translate.goodsn.step')
    unicorn.zp.loop_zeropoints(root='goodsn.step', tfile='zphot.translate.goodsn.step',  zfile='zphot.zeropoint.goodsn.step', fix_zspec=True, fix_filter={205:1, 236:1}, ignore_initial=['f_F435W'], toler=0.01, check_uvj=False)
    
    f0, zp0 = np.loadtxt('zphot.zeropoint.goodsn.step', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
    #
    fix_hst = {203: 1.0099,
     204: 1.0129,
     205: 1.0,
     233: 1.0237,
     236: 1.0,
     238: 0.9899,
     240: 0.9877}
     
    os.system('eazy -p zphot.param.goodsn -t zphot.translate.goodsn')
    unicorn.zp.loop_zeropoints(root='goodsn', tfile='zphot.translate.goodsn',  zfile='zphot.zeropoint.goodsn', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_U','f_J'], fix_zspec=False, use_tweaked_templates=False)
    
    init_filter = {227:0.85, 231:0.97, 224:0.9, 20:1.04, 21:1.04}
    os.system('eazy -p zphot.param.goodsn -t zphot.translate.goodsn')
    unicorn.zp.loop_zeropoints(root='goodsn', tfile='zphot.translate.goodsn',  zfile='zphot.zeropoint.goodsn', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_U','f_J'], fix_zspec=False, use_tweaked_templates=False, init_filter=init_filter)

    ## HDFN
    fix=False #True
    init_filter = {227:0.85, 231:0.97, 224:0.9, 20:1.04, 21:1.04, 114:0.8, 115:0.8, 116:0.7, 117:0.75, 118:0.85}
    os.system('eazy -p zphot.param.goodsn -t zphot.translate.goodsn')
    unicorn.zp.loop_zeropoints(root='goodsn', tfile='zphot.translate.goodsn',  zfile='zphot.zeropoint.goodsn', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_U','f_J', 'f_B', 'f_V', 'f_R', 'f_I', 'f_Z'], use_tweaked_templates=False, init_filter=init_filter, fix_zspec=fix, check_uvj=(not fix))
    
    #### GOODS-S
    cat = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.GABODS.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.FULL.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.fixconv.bright', force_lowercase=False)
    keep = cat.z_spec > 0
    cat.write_text(cat.filename+'.zspec', select=keep)
    
    os.system('eazy -p zphot.param.goodss.step -t zphot.translate.goodss.step')
    unicorn.zp.loop_zeropoints(root='goodss.step', tfile='zphot.translate.goodss.step',  zfile='zphot.zeropoint.goodss.step', fix_zspec=True, fix_filter={205:1}, ignore_initial=['f_F435W'], toler=0.01, check_uvj=False)
    
    f0, zp0 = np.loadtxt('zphot.zeropoint.goodss.step', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
    
    fix_hst = {1: 1.0819,
         4: 1.0033,
         5: 0.9845,
         7: 0.9838,
         203: 1.0028,
         204: 1.0072,
         205: 1.0,
         236: 1.0038,
         239: 0.9919,
         240: 1.0022}
    #
    # os.system('eazy -p zphot.param.goodss -t zphot.translate.goodss')
    # unicorn.zp.loop_zeropoints(root='goodss', tfile='zphot.translate.goodss', fix_zspec=False,  zfile='zphot.zeropoint.goodss', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_R', 'f_H', 'f_Ks'], use_tweaked_templates=True) #, 128:0.9, 129:0.9,

    #### WIth GABODS images
    os.system('eazy -p zphot.param.goodss -t zphot.translate.goodss')
    init_filter = {46:0.72, 50:0.95, 54: 0.95, 58: 0.91, 36:1.07, 37:1.04, 103:1.07, 107:1.02}
    unicorn.zp.loop_zeropoints(root='goodss', tfile='zphot.translate.goodss', fix_zspec=False,  zfile='zphot.zeropoint.goodss', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_R', 'f_H', 'f_Ks', 'f_B', 'f_V', 'f_Rc', 'f_I', 'f_U38'], init_filter=init_filter, use_tweaked_templates=True) #, 128:0.9, 129:0.9,
    
    os.system('eazy -p zphot.param.goodss -t zphot.translate.goodss')
    unicorn.zp.loop_zeropoints(root='goodss', tfile='zphot.translate.goodss', fix_zspec=True, check_uvj=False,  zfile='zphot.zeropoint.goodss', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_R', 'f_H', 'f_Ks', 'f_B', 'f_V', 'f_Rc', 'f_I', 'f_U38'], init_filter=init_filter, use_tweaked_templates=True) #, 128:0.9, 129:0.9,
    
    ####### WIth all Subaru Medium bands!
    init_filter = {46:0.72, 50:0.95, 54: 0.95, 58: 0.91, 36:1.07, 37:1.04, 103:1.07, 107:1.02, 181:0.72, 182:0.62, 183:1.2, 184:3, 186:0.96, 187:0.9, 188:1.05, 190:0.9, 194:0.93, 195:0.92, 196:0.92, 197:1.08, 198:0.92}
    
    ##### "fixconv" catalog with some fixes
    init_filter = {46:1.0, 50:0.95, 54: 0.95, 58: 1.0, 36:1.07, 37:1.04, 103:1.07, 107:1.16, 181:0.9, 182:0.9, 183:1.2, 184:3, 185:0.94, 186:0.96, 187:0.9, 188:1.0, 189:0.9, 190:0.9, 194:0.93, 195:0.92, 196:0.88, 197:1.0, 198:0.85, 220:0.85, 222:0.75}
    
    os.system('eazy -p zphot.param.goodss -t zphot.translate.goodss')
    unicorn.zp.loop_zeropoints(root='goodss', tfile='zphot.translate.goodss', fix_zspec=False, check_uvj=True,  zfile='zphot.zeropoint.goodss', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_R', 'f_H', 'f_Ks', 'f_B', 'f_V', 'f_Rc', 'f_I', 'f_U38', 'f_IA427', 'f_IA445', 'f_IA464', 'f_IA484', 'f_IA505', 'f_IA527', 'f_IA550', 'f_IA574', 'f_IA598', 'f_IA624', 'f_IA651', 'f_IA679', 'f_IA738', 'f_IA767', 'f_IA797', 'f_IA827', 'f_IA856'], init_filter=init_filter, use_tweaked_templates=True) #, 128:0.9, 129:0.9,

    os.system('eazy -p zphot.param.goodss -t zphot.translate.goodss')
    unicorn.zp.loop_zeropoints(root='goodss', tfile='zphot.translate.goodss', fix_zspec=True, check_uvj=False,  zfile='zphot.zeropoint.goodss', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_R', 'f_H', 'f_Ks', 'f_B', 'f_V', 'f_Rc', 'f_I', 'f_U38', 'f_IA427', 'f_IA445', 'f_IA464', 'f_IA484', 'f_IA505', 'f_IA527', 'f_IA550', 'f_IA574', 'f_IA598', 'f_IA624', 'f_IA651', 'f_IA679', 'f_IA738', 'f_IA767', 'f_IA797', 'f_IA827', 'f_IA856'], init_filter=init_filter, use_tweaked_templates=True) #, 128:0.9, 129:0.9,
    
    #### AEGIS
    cat = catIO.Readfile('../Catalogs/aegis_3dhst.v4.0.test2.nzpcat.bright', force_lowercase=False)
    keep = cat.z_spec > 0
    cat.write_text(cat.filename+'.zspec', select=keep)
    
    os.system('eazy -p zphot.param.aegis.step -t zphot.translate.aegis.step')
    unicorn.zp.loop_zeropoints(root='aegis.step', tfile='zphot.translate.aegis.step',  zfile='zphot.zeropoint.aegis.step', fix_zspec=True, fix_filter={205:1}, ignore_initial=['f_F814W'], toler=0.01, check_uvj=False)
    
    f0, zp0 = np.loadtxt('zphot.zeropoint.aegis.step', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
    
    fix_hst = {203: 0.9966, 204: 0.9967, 205: 1.0, 236: 0.9999, 239: 0.9353}
    #### Fix IRAC
    #for id in np.arange(18,22):
    #    fix_hst[id] = 1.
    
    os.system('eazy -p zphot.param.aegis -t zphot.translate.aegis')
    unicorn.zp.loop_zeropoints(root='aegis', tfile='zphot.translate.aegis', fix_zspec=False,  zfile='zphot.zeropoint.aegis', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_J1', 'f_J2', 'f_J3', 'f_H1', 'f_H2', 'f_K', 'f_U', 'f_I', 'f_Z'], use_tweaked_templates=False, init_filter={88:1.2} ) #, 128:0.9, 129:0.9, 130:0.92, 131:0.9, 132:0.95, 134:0.95})
    
    ##### UDS
    cat = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.test.nzpcat.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.fixconv.bright', force_lowercase=False)
    cat = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.HAWKI.bright', force_lowercase=False)
    keep = cat.z_spec > 0
    cat.write_text(cat.filename+'.zspec', select=keep)
    
    os.system('eazy -p zphot.param.uds.step -t zphot.translate.uds.step')
    unicorn.zp.loop_zeropoints(root='uds.step', tfile='zphot.translate.uds.step',  zfile='zphot.zeropoint.uds.step', fix_zspec=True, fix_filter={205:1}, ignore_initial=['f_F814W'], toler=0.01, check_uvj=False)
    
    f0, zp0 = np.loadtxt('zphot.zeropoint.uds.step', unpack=True, dtype=np.str)
    fix_hst = {}
    for i in range(len(f0)):
        fix_hst[int(f0[i][1:])] = float(zp0[i])
        
    #### Fix HST, IRAC, U
    fix_hst = {203: 1.0051, 204: 1.025, 205: 1.0, 18:1., 19:1., 20:1., 21:1., 236: 1.0072, 239: 0.9335, 88:1.22}
    fix_hst = {203: 1.0048, 204: 1.0257, 205: 1.0, 236: 1.0079, 239: 0.9342} #, 18:1., 19:1., 20:1., 21:1.}
    fix_hst = {203: 1.0048, 204: 1.0257, 205: 1.0, 236: 1.0079, 239: 0.9342, 20:1.15, 21:1.15} #, 18:1., 19:1., 20:1., 21:1.}

    fix=False

    init_filter = {79:0.9, 123:0.9, 124:0.88, 125:0.95, 88:1.28, 264:1.04, 265:1.04}
    init_filter = {79:0.95, 123:0.86, 124:0.8,  125:0.86, 88:1.25, 264:1.04, 265:1.04, 266:2.3, 269:1.0}

    os.system('eazy -p zphot.param.uds -t zphot.translate.uds')
    unicorn.zp.loop_zeropoints(root='uds', tfile='zphot.translate.uds', fix_zspec=fix, check_uvj=not fix, zfile='zphot.zeropoint.uds', fix_filter=fix_hst, toler=0.005, ignore_initial=['f_u', 'f_V', 'f_R', 'f_i', 'f_z', 'f_J', 'f_H', 'f_IRAC3', 'f_IRAC4', 'f_Y', 'f_Ks'], use_tweaked_templates=False, init_filter=init_filter) #, 128:0.9, 129:0.9, 130:0.92, 131:0.9, 132:0.95, 134:0.95})
    
    #### Fit just U
    f0, zp0 = np.loadtxt('zphot.zeropoint.uds', unpack=True, dtype=np.str)
    fix_all = {}
    for i in range(len(f0)):
        fix_all[int(f0[i][1:])] = float(zp0[i])
    #
    fix_all = {18: 1.0,
     19: 1.0,
     20: 1.0,
     21: 1.0,
     79: 0.9271,
     122: 1.0291,
     123: 0.9243,
     124: 0.8948,
     125: 0.9715,
     203: 1.0051,
     204: 1.025,
     205: 1.0,
     236: 1.0072,
     239: 0.9335,
     263: 1.0235,
     264: 1.0543,
     265: 1.0411}

    os.system('eazy -p zphot.param.uds -t zphot.translate.uds')
    unicorn.zp.loop_zeropoints(root='uds', tfile='zphot.translate.uds', fix_zspec=False, check_uvj=True, zfile='zphot.zeropoint.uds', fix_filter=fix_all, toler=0.005, ignore_initial=[], use_tweaked_templates=False, init_filter={88:1.3}) 
    
    #### Compare Galametz and Skelton Y
    gal = pyfits.open('../Catalogs/CANDELS.UDS.F160W.v1.fits')[1].data
    ros = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.HAWKI.bright')
    mat = catIO.CoordinateMatcher(ros)
    dr, idx = mat.match_list(gal.RA, gal.Dec)
    ok = dr < 0.1
    
    hmag = 23.9-2.5*np.log10(gal.Flux_F160W_hst)
    dmag = -2.5*np.log10(gal.Flux_Y_hawki*10**(-0.4*(23.9-25))/ros.f_Y[idx])
    split = dmag < -0.933
    plt.scatter(hmag[ok & split], dmag[ok & split], alpha=0.1, color='blue', label='H')
    plt.scatter(hmag[ok & ~split], dmag[ok & ~split], alpha=0.1, color='orange', label='H')
    plt.plot([18,24], [23.9-25,23.9-25], color='black', label='23.9-25')
    dmag = -2.5*np.log10(gal.Flux_Ks_hawki*10**(-0.4*(23.9-25))/ros.f_Ks[idx])
    plt.scatter(hmag[ok], dmag[ok], alpha=0.1, color='red', label='Ks')

    dmag = -2.5*np.log10(gal.Flux_F160W_hst*10**(-0.4*(23.9-25))/ros.f_F160W[idx])
    plt.scatter(hmag[ok], dmag[ok], alpha=0.1, color='green', label='F160W', zorder=-10)

    plt.legend(loc='upper left')
    plt.ylim(-1.5,0.6)
    plt.xlim(18,24)
    plt.xlabel('F160W mag')
    plt.ylabel(r'$\Delta$mag')
    
    plt.savefig('UDS_HAWKI_compare.png')
    
    #dmag = -2.5*np.log10(gal.Flux_Y_hawki*10**(-0.4*(23.9-25))/ros.f_Y[idx])
    #plt.scatter(ros.flux_radius_F160W[idx][ok], dmag[ok], alpha=0.1, color='blue')
    #plt.scatter(ros.flux_radius_F160W[idx][ok & split], dmag[ok & split], alpha=0.1, color='red')

    plt.scatter(gal.RA[ok & split], gal.Dec[ok & split], color='blue', alpha=0.1)
    plt.scatter(gal.RA[ok & ~split], gal.Dec[ok & ~split], color='orange', alpha=0.1)
    plt.xlabel('RA'); plt.ylabel('Dec')
    plt.savefig('UDS_HAWKI_compare_spatial.png')
    
    #########
    os.system('eazy -p zphot.param.cosmos -t zphot.translate.cosmos')
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter={205:1, 239:1}, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'], toler=0.01) #, ignore_all=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'])
    
    os.system('eazy -p zphot.param.aegis -t zphot.translate.aegis')
    unicorn.zp.loop_zeropoints(root='aegis', tfile='zphot.translate.aegis',  zfile='zphot.zeropoint.aegis', fix_filter={205:1, 239:1}, toler=0.01)

    os.system('eazy -p zphot.param.goodsn -t zphot.translate.goodsn')
    unicorn.zp.loop_zeropoints(root='goodsn', tfile='zphot.translate.goodsn',  zfile='zphot.zeropoint.goodsn', fix_filter={205:1, 5:1}, toler=0.01, ignore_initial=['f_U','f_J'])

    os.system('eazy -p zphot.param.uds -t zphot.translate.uds')
    unicorn.zp.loop_zeropoints(root='uds', tfile='zphot.translate.uds',  zfile='zphot.zeropoint.uds', fix_filter={205:1, 239:1}, toler=0.01, ignore_initial=['f_u','f_i'])
    
def loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter={}, ref_filter=None, init_filter={}, ignore_initial=[], ignore_all=[], toler=0.005, PATH='./OUTPUT/', fix_zspec=False, check_uvj=True, use_tweaked_templates=True, wclip=[1200, 3.e4]):
    
    import threedhst
    import threedhst.eazyPy as eazy
    import os
    
    if not PATH.endswith('/'):
        PATH += '/'
    
    ##### Read the param file
    param = eazy.EazyParam('%s.param' %(os.path.join(PATH, root)))
    tf = eazy.TranslateFile(tfile)
    
    ##### Make an empty zeropoint file
    fp = open(zfile,'w')
    for filter in param.filters:
        if filter.fnumber in fix_filter.keys():
            fp.write('F%d  %s\n' %(filter.fnumber, fix_filter[filter.fnumber]))
        else:
            if filter.fnumber in init_filter.keys():
                fp.write('F%d  %s\n' %(filter.fnumber, init_filter[filter.fnumber]))
            else:
                fp.write('F%d  1.0\n' %(filter.fnumber))
    #
    fp.close()
    
    clean_files = [os.path.join(PATH, root)+'.zeropoint', 'tweak/tweak.dat']
    for file in clean_files:
        if os.path.exists(file):
            os.remove(file)
    #
    NT = len(param.templates)
    fp = open('tweak/spectra.param','w')
    for i in range(NT):
        wt, ft = np.loadtxt(param.templates[i], unpack=True)
        np.savetxt('tweak/%s' %(os.path.basename(param.templates[i])), np.array([wt, ft]).T, fmt='%.5e')
        fp.write('%d tweak/%s 1.0 0 1.0\n' %(i+1, os.path.basename(param.templates[i])))
    
    fp.close()
    
    param.params['GET_ZP_OFFSETS'] = True
    param.params['FIX_ZSPEC'] = fix_zspec
    
    param.write('zphot.param.iter')
    
    #### Scale huge errors for first step
    for filter in ignore_initial:
        tf.change_error(filter, 1.e6)
    
    tf.write('zphot.translate.iter')
    
    os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
    
    fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, ref_filter=ref_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, 0), wclip=wclip)
    
    #### Extract UVJ
    if check_uvj:
        for color in ['153,155', '155,161']:
            param.params['REST_FILTERS'] = color
            param.params['READ_ZBIN'] = 'y'
            param.write('zphot.param.iter')
            os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
        #
        param.params['READ_ZBIN'] = 'n'
        param.write('zphot.param.iter')
        unicorn.zp.diagnostic_UVJ(root=root, ext='UVJ_%04d' %(0))
        
    fp = open('%s_iter.log' %(root),'w')
    fp.write('\n\nIter #%d\n======\n' %(0))
    unicorn.zp.log_offsets(fp, fnumbers, lc_i, delta_i, toler)
    
    ### Now loop
    for filter in ignore_initial:
        tf.change_error(filter, 1.)
    #
    for filter in ignore_all:
        tf.change_error(filter, 1.e6)
    
    tf.write('zphot.translate.iter')
    
    param.params['GET_ZP_OFFSETS'] = True
    param.write('zphot.param.iter')
    
    for i in range(15):
        #### Use tweaked templates after two iterations
        if i == 0:
            try:
                os.remove('tweak/tweak.dat')
            except:
                pass
            #
        
        if (i == 1) & (use_tweaked_templates):
            param.params['TEMPLATES_FILE'] = 'tweak/spectra.param'
            param.write('zphot.param.iter')
        
        os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
        fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, ref_filter=ref_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, i+1), wclip=wclip)
        fp.write('\n\nIter #%d\n======\n' %(i))
        unicorn.zp.log_offsets(fp, fnumbers, lc_i, delta_i, toler)
        
        #
        #### Extract UVJ
        if check_uvj:
            for color in ['153,155', '155,161']:
                param.params['REST_FILTERS'] = color
                param.params['READ_ZBIN'] = 'y'
                param.write('zphot.param.iter')
                os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
            #
            param.params['READ_ZBIN'] = 'n'
            param.write('zphot.param.iter')
            unicorn.zp.diagnostic_UVJ(root=root, ext='UVJ_%04d' %(i+1))
        
        fixed_bands = lc_i < 0
        for fi in fix_filter.keys():
            fixed_bands = fixed_bands | (fnumbers == fi)
        
        use_bands = (lc_i > 4500) & (lc_i < 2.5e4) & ~fixed_bands
        
        if (np.abs(delta_i[use_bands]-1).max() < toler) & (i >= 1):
            break
    
    fp.close()
    
    pass
    
def log_offsets(fp, fnumbers, lc_i, delta_i, toler):
    so = np.argsort(lc_i)
    for j in so:
        if np.abs(delta_i[j]-1) > toler:
            log = '* F%d  %.4f' %(fnumbers[j], delta_i[j])
        else:
            log = '  F%d  %.4f' %(fnumbers[j], delta_i[j])
        #
        print log
        fp.write(log+'\n')

def diagnostic_UVJ(root='cosmos', ext='_UVJ_0000'):
    
    zout = catIO.Readfile('OUTPUT/%s.zout' %(root))
    uv = catIO.Readfile('OUTPUT/%s.153-155.rf' %(root))
    vj = catIO.Readfile('OUTPUT/%s.155-161.rf' %(root))
    
    UmV = -2.5*np.log10(uv.l153/uv.l155)
    VmJ = -2.5*np.log10(vj.l155/vj.l161)
    
    fig = unicorn.plotting.plot_init(aspect=0.4, left=0.10, bottom=0.11, xs=8, wspace=0.0, NO_GUI=True)
    
    zbins = [0.5,1,1.5,2]
    for i in range(len(zbins)-1):
        zbin = (zout.z_peak > zbins[i]) & (zout.z_peak <= zbins[i+1])
        ax = fig.add_subplot(131+i)
        ax.scatter(VmJ[zbin], UmV[zbin], color='black', alpha=0.1, s=5)
        ax.set_xlim(-0.2,2.1)
        ax.set_ylim(0.1,2.4)
        if i == 0:
            ax.set_ylabel('U-V')
        else:
            ax.set_yticklabels([])
        #
        ax.text(0.05,0.95, '%.1f < z < %.1f; N=%d' %(zbins[i], zbins[i+1], zbin.sum()), ha='left', va='top', transform=ax.transAxes)
        ax.set_xlabel('V-J')
    
    unicorn.plotting.savefig(fig, '%s_%s.png' %(root, ext))
    
def compare_J():
    """
    Simple test:   compare J2-F125W catalog colors as f'n of mag/z to EAZY templates
    """
    
    from threedhst import catIO
    import threedhst.eazyPy as eazy
    import matplotlib.pyplot as plt
    import numpy as np
    
    PATH = './OUTPUT/'
    
    cat = catIO.Readfile('../Catalogs/aegis_3dhst.v4.0.test2.nzpcat.bright', force_lowercase=False)
    root='aegis'
    f1, f2 = 10, 8 # J2, F125W

    root='cosmos'
    f1, f2 = 16, 14
    
    zout = catIO.Readfile('OUTPUT/%s.zout' %(root))    
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH,CACHE_FILE = 'Same')
    param = eazy.EazyParam('%s%s.param' %(PATH, root))
    
    for i in range(N):
        print '%d: %s' %(i, param.filters[i].name)
    
    
    color_obs = -2.5*np.log10(tempfilt['fnu'][f1]/tempfilt['fnu'][f2])
    color_template = -2.5*np.log10(tempfilt['tempfilt'][f1,:,:]/tempfilt['tempfilt'][f2,:,:])
    
    ok = (zout.z_peak > 0.1) & np.isfinite(color_obs) & (color_obs != 0)
    
    plt.scatter(zout.z_peak[ok], color_obs[ok], color='black', alpha=0.5)
    xm, ym, ys, nn = threedhst.utils.runmed(zout.z_peak[ok], color_obs[ok], NBIN=ok.sum()/500)
    plt.fill_between(xm, ym-ys, ym+ys, color='red', alpha=0.1)
    plt.plot(xm, ym, color='red', linewidth=2)
    
    for i in range(tempfilt['NTEMP']):
        plt.plot(tempfilt['zgrid'], color_template[i,:], color='blue', alpha=0.8, linewidth=2)
        
    plt.ylim(-1,1)
    plt.xlim(0,2.5)
    plt.ylabel('color')
    plt.xlabel('z_phot')
    plt.title('color: %s - %s' %(os.path.basename(param.filters[f1].name), os.path.basename(param.filters[f2].name)))
    
    hmag = 25-2.5*np.log10(cat.f_F125W)
    plt.scatter(hmag[ok], color_obs[ok], alpha=0.5)
    plt.ylim(-1,1)
    plt.xlim(18,26)
    plt.ylabel('J2-F125W')
    plt.xlabel('F125W')
    
def template_error(root='cosmos', PATH='OUTPUT/'):
    
    #zout = catIO.Readfile('%s/%s.zout' %(PATH, root))
    
    param = eazy.EazyParam('%s/%s.param' %(PATH, root))
    
    ##### Read template fluxes and coefficients
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH,CACHE_FILE = 'Same')
    
    if coeffs['izbest'].max() == 0:
        STAR_FIT = True
    else:
        STAR_FIT = False
        
    param = eazy.EazyParam(PARAM_FILE=os.path.join(PATH, root+'.param'))
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
        
    zpfile = PATH+'/'+root+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])
        
    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1), np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))
    
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    obs_sed = np.zeros((tempfilt['NFILT'], tempfilt['NOBJ']), dtype=np.float)
    for i in xrange(tempfilt['NOBJ']):
        obs_sed[:,i] = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][i]], coeffs['coeffs'][:,i])
    
    zi = tempfilt['zgrid'][coeffs['izbest']]
    lc = tempfilt['lc']
    offsets = lc*0.
    
    residual = (obs_sed-tempfilt['fnu'])/tempfilt['fnu']
    observed_error = tempfilt['efnu']
        
    hr = (-11,-1)
    plt.hist(np.log(residual.flatten()**2), range=hr, bins=100, alpha=0.4, color='blue', histtype='step')

    deviates = np.random.normal(size=residual.shape) #*tempfilt['efnu']/tempfilt['fnu']

    in_err = tempfilt['efnu']/tempfilt['fnu']
    plt.hist(np.log((deviates*in_err).flatten()**2), range=hr, bins=100, alpha=0.4, color='red', histtype='step')
    
    te = eazy.TemplateError()
    lc_rest = np.dot(lc.reshape((-1,1)), 1/(1+zi.reshape((1,-1))))
    terr = te.interpolate(lc_rest.flatten(), tempfilt['fnu'].flatten()*0.).reshape(lc_rest.shape)
    
    grow_err = deviates**2*(in_err**2+(terr/2)**2)
    yh_te_hist, xh = np.histogram(np.log(grow_err.flatten()), range=hr, bins=100)
    plt.plot(xh[1:], yh_te_hist, alpha=0.8, color='green')
    
    from astroML.sum_of_norms import sum_of_norms, norm
    ok = (te.te_x > 800) & (te.te_x < 8.e4)
    w_best, rms, locs, widths = sum_of_norms(np.log(te.te_x)[ok], te.te_y[ok]*0.5, 15, spacing='linear', full_output=True)
    norms = w_best * norm(np.log(te.te_x)[:, None], locs, widths)
    plt.plot(te.te_x, te.te_y*0.5, color='black')
    plt.plot(te.te_x, norms, color='red', alpha=0.1)
    plt.plot(te.te_x, norms.sum(1), color='red')
    
    signoise = tempfilt['fnu']/np.sqrt(tempfilt['efnu']**2+(0.01*tempfilt['fnu'])**2)
    
    ##### Run a chain
    init = w_best #*0.+0.01
    step_sig = w_best*0.+0.01
    init[0], step_sig[0] = 0.01, 0.001
    
    hist = ((-11,-1), 100)
    ok = (lc_rest > 500) & (lc_rest < 8.e4)
    residual_hist, xh = np.histogram(np.log(residual[ok]**2), range=hist[0], bins=hist[1])
    
    obj_fun = unicorn.zp._objective_template_error
    obj_args = [residual_hist, deviates[ok].flatten(), in_err[ok].flatten(), lc_rest[ok].flatten(), np.log(te.te_x), locs, widths, hist, None]

    obj_fun = unicorn.zp._objective_template_error_chi2
    obj_args = [residual[ok].flatten(), in_err[ok].flatten(), lc_rest[ok].flatten(), np.log(te.te_x), locs, widths, hist, None]
    
    # residual_hist_0, deviates_0, in_err_0, lc_rest_0, nx, nloc, nwid, hist, axes = obj_args
    ndim, nwalkers = len(init), 80
    p0 = [(init+np.random.normal(size=ndim)*step_sig) 
          for i in xrange(nwalkers)]
    #
    # plt.plot(xh[1:], residual_hist, color='black')
    # for i in xrange(nwalkers):
    #     obj_fun(p0[i], *obj_args)
        
    NTHREADS, NSTEP = 1, 200
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obj_fun, args = obj_args, 
                                    threads=NTHREADS)
    
    result = sampler.run_mcmc(p0, NSTEP)
    
    te0 = (w_best * norm(np.log(te.te_x)[:, None], locs, widths)).sum(1)
    plt.plot(te.te_x, te0, color='black', linewidth=2)
    te0 = (init * norm(np.log(te.te_x)[:, None], locs, widths)).sum(1)
    plt.plot(te.te_x, te0, color='red', linewidth=2)

    for i in range(nwalkers):
        te_i = (sampler.chain[i,NSTEP-1,:]* norm(np.log(te.te_x)[:, None], locs, widths)).sum(1)
        #te_i = (p0[i]* norm(np.log(te.te_x)[:, None], locs, widths)).sum(1)
        plt.plot(te.te_x, te_i, color='purple', alpha=0.05)
    
    chain = sampler.chain[:,-40:,:].reshape((-1,15))
    #chain[:,-2:] = 0.08
    mean = np.dot(chain, norm(np.log(te.te_x)[:, None], locs, widths).T)
    plt.plot(te.te_x, mean.sum(0)/2000., color='purple', alpha=0.5, linewidth=2)
    best = mean.sum(0)/2000.
    ### IR peak
    best[te.te_x > te.te_x[np.argmax(best)]] = best.max()
    ### Lyman-break
    xuv, suv = 1500, 300.
    scl = np.exp(-(te.te_x-xuv)**2/2./suv**2)
    best[te.te_x <= xuv] = scl[te.te_x <= xuv]*best[te.te_x >= xuv][0]
    
    plt.plot(te.te_x, best, color='orange', linewidth=2)
    
    np.savetxt('templates/TEMPLATE_ERROR.calibrated.vX', np.array([te.te_x, best]).T, fmt='%.6e %.4f')
    #np.savetxt('templates/TEMPLATE_ERROR.calibrated.vX', np.array([te.te_x, te0]).T, fmt='%.6e %.4f')
    
    #### Show distributions
    obj_args = [residual_hist, deviates[ok].flatten(), in_err[ok].flatten(), lc_rest[ok].flatten(), np.log(te.te_x), locs, widths, ((-11,-1), 100), plt]
    
    plt.plot(xh[1:], residual_hist, color='black', linewidth=2)
    plt.plot(xh[1:], yh_te_hist, color='green')
    for i in xrange(nwalkers):
        obj_fun(sampler.chain[i,-1,:], *obj_args)
    
    #median = np.median(chain, axis=0)*norm(np.log(te.te_x)[:, None], locs, widths))
    #plt.plot(te.te_x, median.sum(1), color='red', alpha=0.5, linewidth=2)
    
ITER = 1

def _objective_template_error_chi2(params, residuals_0, in_err_0, lc_rest_0, nx, nloc, nwid, hist, axes):
    """
    Objective function for fitting the distribution of fit residuals
    """
    from astroML.sum_of_norms import norm
    import unicorn.zp
    import scipy.stats
    
    #template_error = np.clip(params*norm(nx[:, None], nloc, nwid)).sum(1),0,1)
    template_error = (params*norm(nx[:, None], nloc, nwid)).sum(1)
    
    te_interp = np.interp(lc_rest_0, np.exp(nx), template_error)
    err = (in_err_0**2+te_interp**2)
    lnprob = scipy.stats.chi2.logpdf(np.sum(residuals_0**2/err), err.size-1)#-1e8*(template_error.min() < 0)
    #yh, xh = np.histogram(np.log(err), range=hist[0], bins=hist[1])
    #yh2, xh2 = np.histogram(np.log(residual**2), range=hist[0], bins=hist[1])
    
    # lnprob = -0.5*np.sum((residual_hist_0-yh)**2/residual_hist_0)
    # 
    # if axes is not None:
    #     axes.plot(xh[1:], yh, color='red', alpha=0.1)
    
    unicorn.zp.ITER += 1
    print '%10d %d' %(lnprob, unicorn.zp.ITER)
    
    return lnprob
        
def _objective_template_error(params, residual_hist_0, deviates_0, in_err_0, lc_rest_0, nx, nloc, nwid, hist, axes):
    """
    Objective function for fitting the distribution of fit residuals
    """
    from astroML.sum_of_norms import norm
    import unicorn.zp
    
    template_error = (params*norm(nx[:, None], nloc, nwid)).sum(1)
    te_interp = np.interp(lc_rest_0, np.exp(nx), template_error)
    err = deviates_0**2*(in_err_0**2+te_interp**2)
    yh, xh = np.histogram(np.log(err), range=hist[0], bins=hist[1])
    #yh2, xh2 = np.histogram(np.log(residual**2), range=hist[0], bins=hist[1])
    
    lnprob = -0.5*np.sum((residual_hist_0-yh)**2/residual_hist_0)
    
    if axes is not None:
        axes.plot(xh[1:], yh, color='red', alpha=0.1)
    
    unicorn.zp.ITER += 1
    print lnprob, unicorn.zp.ITER
    
    return lnprob

def compare_filters():
    """
    Compare EAZY filters with P. Capak's COSMOS compilation
    """
    import pysynphot as S
    import threedhst.eazyPy as eazy
    
    PATH='/usr/local/share/eazy-filters'
    
    ff = eazy.FilterFile(PATH+'/FILTER.RES.latest.mod')
    
    filters = {78:'Subaru/B_subaru.res', 79:'Subaru/V_subaru.res', 80:'Subaru/g_subaru.res', 81:'Subaru/r_subaru.res', 82:'Subaru/i_subaru.res', 83:'Subaru/z_subaru.res'}
    
    filters = {88:'CFHT/u_megaprime_sagem.res', 89:'CFHT/g_megaprime_sagem.res', 90:'CFHT/r_megaprime_sagem.res', 91:'CFHT/i_megaprime_sagem.res', 92:'CFHT/z_megaprime_sagem.res'}
    
    i = 78
    for i in filters.keys():
        w, t = np.loadtxt(os.path.join(PATH, filters[i]), unpack=True)
        plt.plot(w,t,color='red', alpha=0.5, linewidth=2)
        f = ff.filters[i-1]
        plt.plot(f.wavelength, f.transmission, color='blue', alpha=0.5, linewidth=2)
    
    #### Build Subaru bandpasses
    ccd = S.FileBandpass(PATH+'/Subaru/suprime_FDCCD.txt')
    optics = S.FileBandpass(PATH+'/Subaru/subaru_optics.cat')
    atm = S.FileBandpass(PATH+'/Subaru/atm_airmass1.2.cat')
    
    i, filter = 81, S.FileBandpass(PATH+'/Subaru/subaru_filt_rp.txt')
    i, filter = 78, S.FileBandpass(PATH+'/Subaru/subaru_filter_B.txt')
    i, filter = 79, S.FileBandpass(PATH+'/Subaru/subaru_filter_V.txt')
    i, filter = 82, S.FileBandpass(PATH+'/Subaru/subaru_filter_ip.txt')
    
    i, filter = 82, S.FileBandpass(PATH+'/Subaru/subaru_filt_Ic.txt')
    
    total = ccd*optics*atm*filter
    
    f = ff.filters[i-1]
    plt.plot(f.wavelength, f.transmission, color='blue', alpha=0.5, linewidth=2)
    plt.plot(total.wave, total.throughput/total.throughput.max(), color='red', alpha=0.5)
    
    #
    i, filter = 82, S.FileBandpass(PATH+'/Subaru/subaru_filt_Ic.txt')
    total = ccd*optics*atm*filter
    
    through = total.throughput/total.throughput.max()
    ok = np.arange(len(total.wave))[through > 0]
    
    np.savetxt(PATH+'/Subaru/suprime_Ic.dat', np.array([total.wave, through]).T[ok.min()-10:ok.max()+10], fmt='%.5e')

    #
    i, filter = 81, S.FileBandpass(PATH+'/Subaru/subaru_filt_Rc.txt')
    total = ccd*optics*atm*filter
    
    through = total.throughput/total.throughput.max()
    ok = np.arange(len(total.wave))[through > 0]
    
    np.savetxt(PATH+'/Subaru/suprime_Rc.dat', np.array([total.wave, through]).T[ok.min()-10:ok.max()+10], fmt='%.5e')
    
    
    
    
     
    
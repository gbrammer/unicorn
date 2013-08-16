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
    
    c = catIO.Readfile('../Catalogs/cosmos_3dhst.v4.0.test2.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/aegis_3dhst.v4.0.test2.nzpcat', force_lowercase=False)

    c = catIO.Readfile('../Catalogs/goodsn_3dhst.v4.0.test2.nzpcat', force_lowercase=False)
    
    hmag = 25-2.5*np.log10(c.f_F160W)
    plt.scatter(hmag, c.flux_radius_F160W, alpha=0.5, color='black')
    plt.xlim(11,25); plt.ylim(0,10)
    
    xstar, ystar = [16,24], [2.8,2.4]
    yint = np.interp(hmag, xstar, ystar)
    is_star = (hmag > xstar[0]) & (hmag < xstar[1]) & (c.flux_radius_F160W < yint)    
    plt.scatter(hmag[is_star], c.flux_radius_F160W[is_star], alpha=0.5, color='orange')
    
    #### Only use stars well-spanned by the Kurucz models
    iH = -2.5*np.log10(c.f_F814W/c.f_F160W)
    sp_type = (iH > -0.2) & (iH < 0.6)
    plt.scatter(hmag[is_star & sp_type], c.flux_radius_F160W[is_star & sp_type], alpha=0.5, color='red')
    
    c.write_text(filename='../Catalogs/uds_3dhst.v4.0.test.nzpcat.star', select=is_star & sp_type)

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

def loop_cosmos():
    import unicorn.zp
    
    os.system('eazy -p zphot.param.cosmos -t zphot.translate.cosmos')
    unicorn.zp.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter={205:1, 239:1}, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'], toler=0.01) #, ignore_all=['f_Zp', 'f_Ip', 'f_Rp', 'f_V', 'f_IA427', 'f_IA464', 'f_IA624', 'f_IA679', 'f_IA709'])
    
    os.system('eazy -p zphot.param.aegis -t zphot.translate.aegis')
    unicorn.zp.loop_zeropoints(root='aegis', tfile='zphot.translate.aegis',  zfile='zphot.zeropoint.aegis', fix_filter={205:1, 239:1}, toler=0.01)

    os.system('eazy -p zphot.param.goodsn -t zphot.translate.goodsn')
    unicorn.zp.loop_zeropoints(root='goodsn', tfile='zphot.translate.goodsn',  zfile='zphot.zeropoint.goodsn', fix_filter={205:1, 5:1}, toler=0.01, ignore_initial=['f_U','f_J'])

    os.system('eazy -p zphot.param.uds -t zphot.translate.uds')
    unicorn.zp.loop_zeropoints(root='uds', tfile='zphot.translate.uds',  zfile='zphot.zeropoint.uds', fix_filter={205:1, 239:1}, toler=0.01, ignore_initial=['f_u','f_i'])
    
def loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter=None, ignore_initial=[], ignore_all=[], toler=0.005, PATH='./OUTPUT/', fix_zspec=False):
    
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
    
    param.params['GET_ZP_OFFSETS'] = False
    param.params['FIX_ZSPEC'] = fix_zspec
    
    param.write('zphot.param.iter')
    
    #### Scale huge errors for first step
    for filter in ignore_initial:
        tf.change_error(filter, 1.e6)
    
    tf.write('zphot.translate.iter')
    
    os.system('eazy -p zphot.param.iter -t zphot.translate.iter')
    
    fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, 0))
    
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
        
        if i == 1:
            param.params['TEMPLATES_FILE'] = 'tweak/spectra.param'
            param.write('zphot.param.iter')
        
        os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
        fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, i+1))
        fp.write('\n\nIter #%d\n======\n' %(i))
        unicorn.zp.log_offsets(fp, fnumbers, lc_i, delta_i, toler)
        #
        use_bands = (lc_i > 4500) & (lc_i < 2.5e4)
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
    
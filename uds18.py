"""
Custom fits for the lens in UDS-18
"""
import numpy as np
import scipy.stats
import scipy.ndimage as nd

import pyfits

import emcee
import unicorn

z_lens = 1.625
z_arc = 2.263

def test():
    
    twod_model = unicorn.reduce.Interlace2D('model_31684.2D.fits')
    twod_orig = unicorn.reduce.Interlace2D('orig_31684.2D.fits')
    twod_line = unicorn.reduce.Interlace2D('acs_31684.2D.fits')
    
    tx, ty_hi = np.loadtxt('SF0_0.emline.hiOIII.txt', unpack=True)
    #tx, ty_lo = np.loadtxt('SF0_0.emline.loOIII.txt', unpack=True)

    xlens, ylens, xx, yy = np.loadtxt('31684.temp_sed', unpack=True)
    ylens = ylens / (xlens/5500.)**2
    ylens /= np.interp(1.4e4, xlens, ylens)
    twod_model.compute_model(xlens, ylens)
    
    ty = ty_hi
    
    status = twod_line.compute_model(tx*(1+z_arc), ty*3)
    
    zgrid = np.arange(2.25,2.30,0.001)
    chi2 = zgrid*0.
    
    for i in range(len(zgrid)):
        status = twod_line.compute_model(tx*(1+zgrid[i]), ty*2.5)
        mask = (twod_line.model > 1.e-5) & (twod_line.im['WHT'].data > 0)
        resid = (twod_line.im['SCI'].data-twod_line.im['CONTAM'].data - twod_line.model)
        chi2[i] = np.sum(resid[mask]**2/(twod_line.im['WHT'].data[mask]**2))
        
    z_best = zgrid[chi2 == chi2.min()][0]
    status = twod_line.compute_model(tx*(1+z_best), ty*2.5)
    
    import scipy.ndimage as nd
    kernel = twod_line.im['DSCI'].data*1
    kernel[twod_line.im['DSEG'].data != 31684] = 0
    kernel /= kernel.sum()
    resid = twod_orig.im['SCI'].data-twod_orig.im['CONTAM'].data-twod_model.model
    sm_resid = nd.correlate(resid, kernel, mode='reflect', cval=0.0)
    
    status = twod_model.compute_model()
    
def modify_model_seg():
    twod = pyfits.open('UDS-18_31684_model.2D.fits', mode='update')
    sh = twod['DSCI'].data.shape
    yc, xc = np.indices(sh)
    r = np.sqrt((xc-26)**2+(yc-26)**2)
    mask = r < 26
    twod['DSEG'].data[twod['DSEG'].data == 31684] = 0
    twod['DSEG'].data[mask & (twod['DSEG'].data == 0)] = 31684
    twod.flush()
    
def make_segmentation_for_ACS():
    """
    Modify the segmentation of the ACS thumbail to just include parts of the image where 
    the arc has flux.
    """
    import pyregion
    
    twod = pyfits.open('acs_31684.2D.fits', mode='update')
    r = pyregion.open('arc_acs.reg').as_imagecoord(header=twod['DSEG'].header)
    mask = r.get_mask(hdu=twod['DSEG'])
    new_seg = twod['DSEG'].data*0
    new_seg[mask] = 31684
    twod['DSEG'].data = new_seg
    twod.flush()
        
def mcmc_fit():
    """ 
    Run the MCMC fit for various components of the lens system
    """
    import numpy as np
    import scipy.stats
    import pyfits

    import emcee
    import unicorn
    import threedhst.dq
    
    import unicorn.uds18 as uds18
    
    z_lens = 1.625
    z_arc = 2.263
    
    twod_model = unicorn.reduce.Interlace2D('model_31684.2D.fits')
    twod_orig = unicorn.reduce.Interlace2D('orig_31684.2D.fits')
    twod_line = unicorn.reduce.Interlace2D('acs_31684.2D.fits')

    twod_flat = unicorn.reduce.Interlace2D('flat_31684.2D.fits')
    twod_nobcg = unicorn.reduce.Interlace2D('nobcg_31684.2D.fits')
    twod_nopoint = unicorn.reduce.Interlace2D('nopoint_31684.2D.fits')
    
    data = {}
    data['bcg'] = twod_flat.im['CONTAM'].data-twod_nobcg.im['CONTAM'].data
    data['point'] = twod_flat.im['CONTAM'].data-twod_nopoint.im['CONTAM'].data
    data['contam'] = twod_nobcg.im['CONTAM'].data*1-data['point']
    
    data['var'] = twod_orig.im['WHT'].data**2
    data['var'][data['var'] == 0] = 1.e10
    data['mask'] = data['var'] < 1.e10
    data['N'] = (twod_orig.im['WHT'].data > 0).sum()
    
    data['spec2d'] = twod_orig.im['SCI'].data*1
    
    #### emission line components for the arc
    tx_hi, ty_hi = np.loadtxt('SF0_0.emline.hiOIII.txt', unpack=True)
    tx_lo, ty_lo = np.loadtxt('SF0_0.emline.loOIII.txt', unpack=True)
    
    data['twod_line'] = twod_line
    data['arc_templates'] = [[tx_lo, ty_lo], [tx_hi, ty_hi]]
    
    #### Eazy SED of the lens
    xlens, ylens, xx, yy = np.loadtxt('31684.temp_sed', unpack=True)
    ylens = ylens / (xlens/5500.)**2
    ylens /= np.interp(1.4e4, xlens, ylens)
    xlens /= (1+z_lens)
    
    data['twod_lens'] = twod_model
    data['lens_template'] = (xlens, ylens)
    
    ### Helper data
    sh = twod_model.model.shape
    data['x2d'] = np.dot(np.ones(sh[0]).reshape((-1,1)), twod_model.oned.lam.reshape(1,-1))
    data['getModel'] = True
    
    ### init = [lens_a0, lens_beta, z_lens, contam_a0, contam_beta, point_a0, point_beta, bcg_a0, bcg_beta, bcg_yshift, z_arc, arc_lo, arc_hi]
    init = np.array([0, 0, 1.634, 0, 0, 0, 0, 0, 0, 0, 2.263, 0.223, 0.223])
    init = np.array([-0.1, -0.03, 1.634, 0.01, 0, 1, 0, 0.25, 0.2, 0.6, 2.261, 0.223, 0.7])
    
    ## test
    model_init = uds18._objective_full_fit(init, data)
    
    NWALKERS, NSTEP = 100, 300
    step_sigma = np.array([0.2, 0.1, 0.05, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.5, 0.02, 0.3, 0.3])
    ndim = len(init)
    
    p0 = [init + np.random.normal(size=ndim)*step_sigma for i in range(NWALKERS)]
    
    NTHREADS = 1
    
    data['getModel'] = False
    sampler = emcee.EnsembleSampler(NWALKERS, ndim, uds18._objective_full_fit, threads=NTHREADS, args=[data])
    
    result = sampler.run_mcmc(p0, NSTEP)
    
    param_names = ['lens_a0', 'lens_beta', 'z_lens', 'contam_a0', 'contam_beta', 'point_a0', 'point_beta', 'bcg_a0', 'bcg_beta', 'bcg_yshift', 'z_arc', 'arc_lo', 'arc_hi']
    chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names = param_names)
    chain.save_chain()
    
    #### Get best-fit model
    data['getModel'] = True
    model_best = uds18._objective_full_fit(chain.median, data)
    
    ### Model of contaminating objects, without the arc
    p_clean = chain.median*1
    p_clean[-2:] = -1000
    model_clean = uds18._objective_full_fit(p_clean, data)

    ### model of the arc
    p_line = chain.median*1
    p_line[0] = -1000; p_line[3] = -1000; p_line[5] = -1000; p_line[7] = -100
    model_line = uds18._objective_full_fit(p_line, data)
    
    err = np.random.normal(size=model_line.shape)*twod_line.im['WHT'].data
    
    #### Smooth cleaned spectrum by line kernel
    import scipy.ndimage as nd
    kernel = twod_line.im['DSCI'].data*1
    kernel[twod_line.im['DSEG'].data != 31684] = 0
    kernel /= kernel.sum()
    resid = twod_orig.im['SCI'].data-model_clean
    sm_resid = nd.correlate(resid, kernel, mode='reflect', cval=0.0)
    sm_line = nd.correlate(resid, kernel, mode='reflect', cval=0.0)
    sm_model_line = nd.correlate(model_line, kernel, mode='reflect', cval=0.0)
    sm_x = nd.correlate(twod_orig.im['SCI'].data, kernel, mode='reflect', cval=0.0)
    
    ##### Examples using the chain
    chain = unicorn.interlace_fit.emceeChain(file='emcee_chain.pkl')
    
    print chain.stats  # show parameter statistics
    plt.scatter(chain.chain[:,chain.nburn:,-2], chain.chain[:,chain.nburn:,-1], alpha=0.01) # scatter plots of parameters
    draw = chain.draw_random(N=100)
    
    ### Save image in ds9
    ds9.frame(1)
    ds9.view(data['spec2d'])
    ds9.frame(2)
    ds9.view(model_best)
    ds9.frame(3)
    ds9.view(data['spec2d']-model_clean)
    ds9.frame(4)
    ds9.view(model_line)
    
    
def _objective_full_fit(params, data):
    """
    Get log(chi2/2) for some set of input parameters
    
    params:  
        tilt of lens (redshift of lens?)
        tilt of contam
        redshift of arc, coeffs of two line templates spanning Hb/[OIII]
        
    data:
        2D array of contamination
        2D variance array
        Interlace2D object for the line component
        two line templates
        Interlace2D object for the lens
        lens template (wave, flux)
    """
        
    lens_a0 = np.exp(params[0])
    lens_beta = params[1]
    z_lens = params[2]
    
    contam_a0 = np.exp(params[3])
    contam_beta = params[4]
    
    point_a0 = np.exp(params[5])
    point_beta = params[6]
    
    bcg_a0 = np.exp(params[7])
    bcg_beta = params[8]
    bcg_yshift = params[9]
    
    z_arc = params[10]
    arc_lo = np.exp(params[11])
    arc_hi = np.exp(params[12])
    
    ### First component is the lens
    xlens, ylens = data['lens_template']
    xlens_i = xlens*(1+z_lens)
    ylens_i = ylens*lens_a0 * (xlens_i/1.4e4)**lens_beta
    data['twod_lens'].compute_model(xlens_i, ylens_i)
    full_model = data['twod_lens'].model*1
    
    ### Next component is (tilted) contamination
    full_model += contam_a0 * (data['x2d']/1.4e4)**(contam_beta)*data['contam']
    
    ### Tilted "point-source" contamination
    full_model += point_a0 * (data['x2d']/1.4e4)**(point_beta)*data['point']
    
    ### Add shifted, tilted BCG
    bcg = nd.shift(bcg_a0 * (data['x2d']/1.4e4)**(bcg_beta)*data['bcg'], (bcg_yshift, 0))
    full_model += bcg
    
    ### Now make the lens itself
    # low [OIII]/Hbeta
    x_lo, y_lo = data['arc_templates'][0]
    x_lo_i = x_lo * (1+z_arc)
    y_lo_i = y_lo * arc_lo
    data['twod_line'].compute_model(x_lo_i, y_lo_i)
    full_model += data['twod_line'].model
    
    x_hi, y_hi = data['arc_templates'][1]
    x_hi_i = x_hi * (1+z_arc)
    y_hi_i = y_hi * arc_hi
    data['twod_line'].compute_model(x_hi_i, y_hi_i)
    full_model += data['twod_line'].model
    
    ## return model itself?
    if data['getModel']:
        return full_model
    
    ## compute log likelihood
    chi2 = np.sum((data['spec2d']-full_model)**2/data['var'])
    lnprob = -0.5*(chi2) #+0.5*data['N']
    
    ### apply priors
    prior_z_arc = scipy.stats.norm.logpdf(z_arc, loc=2.261, scale=0.02)
    prior_z_lens = scipy.stats.norm.logpdf(z_lens, loc=1.634, scale=0.1)

    prior_lens_a0 = scipy.stats.norm.logpdf(lens_a0, loc=0, scale=2)
    prior_lens_beta = scipy.stats.norm.logpdf(lens_beta, loc=0, scale=2)

    prior_contam_a0 = scipy.stats.norm.logpdf(contam_a0, loc=0, scale=2)
    prior_contam_beta = scipy.stats.norm.logpdf(contam_beta, loc=0, scale=2)

    prior_bcg_a0 = scipy.stats.norm.logpdf(bcg_a0, loc=0, scale=2)
    prior_bcg_beta = scipy.stats.norm.logpdf(bcg_beta, loc=0, scale=2)
    prior_bcg_yshift = scipy.stats.norm.logpdf(bcg_yshift, loc=0.4, scale=1)

    prior_arc_lo = scipy.stats.norm.logpdf(np.log(arc_lo/arc_hi), loc=-0.3, scale=0.5)
        
    lnprob += prior_z_arc + prior_z_lens + prior_lens_a0 + prior_lens_beta + prior_contam_a0 + prior_contam_beta + prior_bcg_a0 + prior_bcg_beta + prior_bcg_yshift + prior_arc_lo
    
    print '%6.3f %6.3f %6.3f   %6.3f %6.3f  %6.3f %6.3f %4.1f  %6.3f   %6.3f  %6.3f  %13.2f' %(lens_a0, lens_beta, z_lens, contam_a0, contam_beta, bcg_a0, bcg_beta, bcg_yshift, z_arc, arc_lo, arc_hi, lnprob)
    
    return lnprob    
    
"""
Simulate interlaced spectra.
"""

def simspec(root='COSMOS-19'):
    """
    Root is the base image where the noise and direct images come from.
    """
    
    #### Simple model of a gaussian Ha emission line
    xflux = np.arange(1.e4,1.8e4)

    dv = 100 # km/s
    z0 = 1.0
    
    l0 = 6564.61*(1+z0)
    dlam = dv/3.e5*l0
    yline = 1./np.sqrt(2*np.pi*dlam**2)*np.exp(-(xflux-l0)**2/2/dlam**2)
    
    ### Use template emission lines rather than a single gaussian
    xflux, yline = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
    xflux *= (1+z0)

    #### Add continuum, here with level 0.1*max(line)
    ycont = yline.max()*0.1
    yflux = ycont+yline
    
    #### Normalize to F140W passband
    x_filt, y_filt = np.loadtxt(os.getenv('iref')+'/F140W.dat', unpack=True)
    y_filt_int = utils_c.interp_c(xflux, x_filt, y_filt)
    filt_norm = np.trapz(y_filt_int*yflux, xflux) / np.trapz(y_filt_int, xflux)
    #yflux /= filt_norm
    
    model = unicorn.reduce.GrismModel(root)
    ids = model.cat.id[model.cat.mag < 24]
    
    #### Generate model where every spectrum is the line template but the mag/shape of the galaxies
    #### is as observed
    for i,id in enumerate(ids):
        print unicorn.noNewLine+'%d (%d/%d)' %(id, i+1, len(ids))
        model.compute_object_model(id, lam_spec=xflux, flux_spec=yflux)
        model.model += model.object
    
    #### Get error array from the error extension    
    err = np.random.normal(size=model.model.shape)*model.gris[2].data
    mask = (err != 0) & (model.segm[0].data == 0)
    
    #### Compare background flux distributions
    #plt.hist(model.gris[1].data[mask].flatten(), range=(-0.1,0.1), bins=100, alpha=0.5)
    #plt.hist(err[mask].flatten(), range=(-0.1,0.1), bins=100, alpha=0.5)
    
    #### Store the new model in the grism image data extension so that we can fit it with the
    #### various tools (z, line strength, etc)
    old = model.gris[1].data*1.
    model.gris[1].data = model.model*(err != 0) + err
    model.get_corrected_wcs(verbose=True)
    model.init_object_spectra()
    model.model*=0
    
    ##### Try extracting a spectrum and fitting it
    id=685
    flam = np.sum(model.flux[model.segm[0].data == id])

    ### *Input* line flux, should be able to get this directly from the input spectrum and the
    ### observed magnitude, but check units.
    #plt.plot(xflux, yflux/filt_norm*flam*1.e-17)
    ha = np.abs(xflux-6564*(1+z0)) < 100
    ha_flux = np.trapz(yline[ha]/filt_norm*flam*1.e-17, xflux[ha])
    
    model.twod_spectrum(id, refine=True, verbose=True)
    model.show_2d(savePNG=True)
    spec = unicorn.reduce.Interlace1D(root+'_%05d.1D.fits' %(id), PNG=True)
    
    #### Redshift fit, set template to flat and the redshift prior to a broad gaussian centered
    #### on the input value, z0
    zgrid = np.arange(0,4,0.005)
    pz = np.exp(-(zgrid-z0)**2/2/0.5**2)
    lnprob = np.log(pz)
    
    obj='%s_%05d' %(root, id)
    gris = unicorn.interlace_fit.GrismSpectrumFit(root=obj, lowz_thresh=0.01)
    gris.zout.z_spec = gris.zout.z_spec*0.+z0
    gris.zout.l99 = gris.zout.l99*0.+z0-0.1
    gris.zout.u99 = gris.zout.l99+0.2
    gris.z_peak = 1
    gris.best_fit = gris.best_fit*0+1
    
    gris.phot_zgrid = zgrid
    gris.phot_lnprob = lnprob
    gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0005, zrfirst=(z0-0.2,z0+0.2))
    
    #### Emission line fit
    gris.fit_free_emlines(ztry=z0, verbose=True, NTHREADS=1, NWALKERS=50, NSTEP=100, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100)
    
    
    # test, FLT errors
    flt = pyfits.open('ibhm47gwq_flt.fits')
    err_flt = np.random.normal(size=flt[1].data.shape)*flt[2].data
    mask_flt = (flt[1].data < 0.1) & (err_flt != 0)
    threedhst.utils.biweight(flt[1].data[mask_flt].flatten())   
    threedhst.utils.biweight(err_flt[mask_flt].flatten())
    
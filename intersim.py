"""
Simulate interlaced spectra.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pyfits

import unicorn
import unicorn.interlace_fit
import unicorn.utils_c as utils_c
import threedhst

def sim_all():
    """
    Run in GRISM_HOME / SIMULATIONS, loop through all pointings and run the 
    spectra simulation.
    """
    import glob
    
    import unicorn
    import unicorn.intersim
    
    files = glob.glob('*G141_inter.fits')
    for file in files:
        root=file.split('-G141')[0]
        unicorn.intersim.simspec(root=root)
        
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
    yflux /= filt_norm
    yline /= filt_norm
    ycont /= filt_norm
    
    ids = [290]

    model = unicorn.reduce.GrismModel(root)
    ids = model.cat.id[model.cat.mag < 24]
    #ids = [245]
    
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
    #old = model.gris[1].data*1.
    model.gris[1].data = model.model*(err != 0) + err
    model.get_corrected_wcs(verbose=True)
    model.init_object_spectra()
    model.model*=0
    
    ##### Try extracting a spectrum and fitting it
    #id=685
    #id=343
    #id=ids[0]
    for id in ids:

        obj='%s_%05d' %(root, id)
        print '%s.linefit.png' %(obj)
        if os.path.exists('%s.linefit.png' %(obj)):
            print 'skip'
            continue
            
        flam = np.sum(model.flux[model.segm[0].data == id])
        fnu = np.sum(model.flux_fnu*(model.segm[0].data == id))

        ### *Input* line flux, should be able to get this directly from the input spectrum and the
        ### observed magnitude, but check units.
        #plt.plot(xflux, yflux/filt_norm*flam*1.e-17)
        ha = np.abs(xflux-6564*(1+z0)) < 100
        ha_flux = np.trapz(yline[ha]*flam*1.e-17, xflux[ha])
        ha_eqw = np.trapz(yline[ha]/ycont, xflux[ha])
        
        s2 = np.abs(xflux-6731*(1+z0)) < 100
        s2_flux = np.trapz(yline[s2]*flam*1.e-17, xflux[s2])
        s2_eqw = np.trapz(yline[s2]/ycont, xflux[s2])

        model.twod_spectrum(id, refine=True, verbose=True)
        if not model.twod_status:
            continue
            
        model.show_2d(savePNG=True)
        spec = unicorn.reduce.Interlace1D(root+'_%05d.1D.fits' %(id), PNG=True)

        #### Redshift fit, set template to flat and the redshift prior to a broad gaussian centered
        #### on the input value, z0
        zgrid = np.arange(0,4,0.005)
        pz = np.exp(-(zgrid-z0)**2/2/0.5**2)
        lnprob = np.log(pz)

        gris = unicorn.interlace_fit.GrismSpectrumFit(root=obj, lowz_thresh=0.01, FIGURE_FORMAT='png')
        if not gris.status:
            continue
            
        gris.zout.z_spec = gris.zout.z_spec*0.+z0
        gris.zout.l99 = gris.zout.l99*0.+z0-0.1
        gris.zout.u99 = gris.zout.l99+0.2
        gris.z_peak = 1
        gris.best_fit = gris.best_fit*0+1

        gris.phot_zgrid = zgrid
        gris.phot_lnprob = lnprob
        gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0005, zrfirst=(z0-0.2,z0+0.2))
        
        if not gris.status:
            continue
            
        #### Emission line fit
        try:
            gris.fit_free_emlines(ztry=gris.z_max_spec, verbose=True, NTHREADS=1, NWALKERS=50, NSTEP=100, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100)
        except:
            pass
            
        status = os.system('cat %s.linefit.dat' %(obj))
        print '\n -- input --\nSII  %6.2f  %6.2f' %(s2_flux/1.e-17, s2_eqw)
        print ' Ha  %6.2f  %6.2f' %(ha_flux/1.e-17, ha_eqw)
    
def results():
    import threedhst.catIO as catIO
    
    files=glob.glob('*linefit.dat')
    
    cat = None
    
    fp = open('simspec.dat','w')
    fp.write('# object mag r50 r90 z_fit continuum_sn ha_flux ha_flux_err ha_eqw ha_eq_err s2_flux s2_flux_err s2_eqw s2_eq_err\n')
    fp.close()
    
    for ii, file in enumerate(files):
        fp = open('simspec.dat','a')
        root = file.split('.linefit')[0]
        print unicorn.noNewLine+'%s (%d/%d)' %(root, ii+1, len(files))
        pointing = root.split('_')[0]
        id = int(root.split('_')[1])
        if cat is None:
            cat = threedhst.sex.mySexCat(pointing+'_inter.cat')
        else:
            if not cat.filename.startswith(pointing):
                cat = threedhst.sex.mySexCat(pointing+'_inter.cat')
        #        
        gris = unicorn.interlace_fit.GrismSpectrumFit(root, verbose=False)
        if not gris.status:
            continue
        #
        result = gris.stats()
        if result is False:
            continue
        #
        DIRECT_MAG, Q_Z, F_COVER, F_FLAGGED, MAX_CONTAM, INT_CONTAM, F_NEGATIVE = result
        #
        lwindow = (gris.oned.data.wave > 1.4e4) & (gris.oned.data.wave < 1.6e4)
        if (lwindow.sum() < 10) | (INT_CONTAM > 0.3):
            continue
        #
        continuum_sn = np.median((gris.oned.data.flux/gris.oned.data.error)[lwindow])
        #
        lfit = catIO.Readfile(root+'.linefit.dat')
        if 'Ha' in lfit.line:
            ix = np.arange(len(lfit.line))[lfit.line == 'Ha'][0]
            ha_flux, ha_flux_err, ha_eqw, ha_eqw_err = lfit.flux[ix], lfit.error[ix], lfit.eqw_obs[ix], lfit.eqw_obs_err[ix]
        else:
            ha_flux, ha_flux_err, ha_eqw, ha_eqw_err = -1,-1,-1,-1
        #
        if 'SII' in lfit.line:
            ix = np.arange(len(lfit.line))[lfit.line == 'SII'][0]
            s2_flux, s2_flux_err, s2_eqw, s2_eqw_err = lfit.flux[ix], lfit.error[ix], lfit.eqw_obs[ix], lfit.eqw_obs_err[ix]
        else:
            s2_flux, s2_flux_err, s2_eqw, s2_eqw_err = -1,-1,-1,-1
        #
        ic = np.arange(cat.nrows)[cat.id == id]
        fp.write(' %s  %6.3f  %6.2f %6.2f %6.4f %6.2f  %6.2f %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f %6.2f\n' %(root, DIRECT_MAG, float(cat.FLUX_RADIUS[ic]), float(cat.FLUX_RADIUS2[ic]), gris.z_max_spec, continuum_sn, ha_flux, ha_flux_err, ha_eqw, ha_eqw_err, s2_flux, s2_flux_err, s2_eqw, s2_eqw_err))
        #
        fp.close()
    
    stats = catIO.Readfile('simspec_full.dat')
    ha_model, s2_model = unicorn.intersim.get_line_fluxes(z0=1.0, mag=stats.mag)
    
    #### Color by r50/r90 concentration
    concentration = stats.r50/stats.r90
    msize = np.maximum((concentration/0.2)**4,4)
    mcol = np.minimum((np.maximum(concentration,0.3)-0.3)/0.3,1)
    plt.scatter(stats.mag, concentration, c=mcol, alpha=0.5)

    
    #### Continuum depth
    BINWIDTH=92
    bin_sn = np.sqrt(BINWIDTH/22)
    plt.scatter(stats.mag, stats.continuum_sn*bin_sn, alpha=0.8, c=mcol)
    plt.ylim(0.1,2000)
    plt.plot([17,24],[5,5], color='black', alpha=0.4)
    plt.xlim(17,24)
    plt.semilogy()
    
    sub = np.abs(stats.mag-23) < 0.5
    plt.scatter(stats.r50[sub], stats.continuum_sn[sub]*bin_sn, c=mcol[sub], alpha=0.8)
    
    #### Line fluxes
    plt.errorbar(ha_model, stats.ha_flux, stats.ha_flux_err, marker='o', markersize=0.1, linestyle='None', color='0.5')
    plt.scatter(ha_model, stats.ha_flux, c=mcol, zorder=100, alpha=0.5)
    #plt.scatter(stats.s2_flux, s2_model, alpha=0.8, c=mc)
    plt.plot([0.1,1000],[0.1,1000], color='black', alpha=0.5)
    plt.xlim(0.5,1000)
    plt.ylim(0.5,1000)
    plt.loglog()
    
    ha_sn = stats.ha_flux/stats.ha_flux_err
    plt.scatter(stats.ha_flux, ha_sn, c=mcol, zorder=100, alpha=0.5)
    plt.plot([0.5,1000],[5,5], color='black', alpha=0.4)
    plt.xlim(0.5,1000)
    plt.ylim(0.5,1000)
    plt.loglog()
    
    #### EQW 
    dha = stats.ha_eqw-130.
    hi = plt.hist(dha/stats.ha_eq_err, range=(-5,5), bins=50, alpha=0.7)
    
    #### redshift
    dz = (stats.z_fit-1)/2.
    
    #### surface density
    mu = stats.mag-2*np.log(stats.r90*0.06)
    plt.scatter(stats.mag, mu, c=mcol)
    
def get_line_fluxes(z0=1.0, mag=21):
    """ 
    Get emission line fluxes for a given continuum magnitude.
    """
    print z0
    
    xflux, yline = np.loadtxt(unicorn.GRISM_HOME+'/templates/dobos11/SF0_0.emline.txt', unpack=True)
    xflux *= (1+z0)

    #### Add continuum, here with level 0.1*max(line)
    ycont = yline.max()*0.1
    yflux = ycont+yline
    
    #### Normalize to F140W passband
    x_filt, y_filt = np.loadtxt(os.getenv('iref')+'/F140W.dat', unpack=True)
    y_filt_int = utils_c.interp_c(xflux, x_filt, y_filt)
    filt_norm = np.trapz(y_filt_int*yflux, xflux) / np.trapz(y_filt_int, xflux)
    yflux /= filt_norm
    yline /= filt_norm
    ycont /= filt_norm
    
    fnu = 10**(-0.4*(mag+48.6))
    flam = fnu*3.e18/(6564.*(1+z0))**2/1.e-17
    
    ha = np.abs(xflux-6564*(1+z0)) < 100
    ha_flux = np.trapz(yline[ha], xflux[ha])
    
    s2 = np.abs(xflux-6731*(1+z0)) < 100
    s2_flux = np.trapz(yline[s2], xflux[s2])
    
    return ha_flux*flam, s2_flux*flam
    
    # #### Trying to figure out units
    # plt.plot(gris.twod.im['WAVE'].data, gris.twod.im['SENS'].data)
    # plt.plot(unicorn.reduce.sens_files['A'].field('WAVELENGTH'), unicorn.reduce.sens_files['A'].field('SENSITIVITY')*1.e-17*np.median(np.diff(gris.twod.im['WAVE'].data))/2**2)
    # 
    # # test, FLT errors
    # flt = pyfits.open('ibhm47gwq_flt.fits')
    # err_flt = np.random.normal(size=flt[1].data.shape)*flt[2].data
    # mask_flt = (flt[1].data < 0.1) & (err_flt != 0)
    # threedhst.utils.biweight(flt[1].data[mask_flt].flatten())   
    # threedhst.utils.biweight(err_flt[mask_flt].flatten())
    
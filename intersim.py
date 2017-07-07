"""
Simulate interlaced spectra.
"""
import os
import glob

from pylab import cm

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np
import astropy.io.fits as pyfits

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
        try:
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0005, zrfirst=(z0-0.2,z0+0.2))
        except:
            continue
            
        if not gris.status:
            continue
            
        #### Emission line fit
        try:
            gris.fit_free_emlines(ztry=gris.z_max_spec, verbose=True, NTHREADS=1, NWALKERS=50, NSTEP=100, FIT_REDSHIFT=False, FIT_WIDTH=False, line_width0=100)
        except:
            continue
            
        status = os.system('cat %s.linefit.dat' %(obj))
        print '\n -- input --\nSII  %6.2f  %6.2f' %(s2_flux/1.e-17, s2_eqw)
        print ' Ha  %6.2f  %6.2f' %(ha_flux/1.e-17, ha_eqw)

    
def get_results(force_new=False):
    """
    Collate the results from the simulated spectra and the input catalogs into single output 
    catalogs suitable for reading and plotting.
    
    for field in ['AEGIS','COSMOS','UDS','GOODS-S']:
        os.chdir(unicorn.GRISM_HOME+'%s/PREP_FLT' %(field))
        unicorn.intersim.get_results()
    
    os.chdir(unicorn.GRISM_HOME+'SIMULATIONS')
    status = os.system('cat ../AEGIS/PREP_FLT/simspec.dat ../COSMOS/PREP_FLT/simspec.dat ../GOODS-S/PREP_FLT/simspec.dat ../UDS/PREP_FLT/simspec.dat > all_simspec.dat')
    
    """
    
    import threedhst.catIO as catIO
    
    files=glob.glob('*linefit.dat')
    
    cat = None
    
    if (not os.path.exists('simspec.dat')) | force_new:
        fp = open('simspec.dat','w')
        fp.write('# object sky_avg sky_lo sky_hi mag r50 r90 z_fit continuum_sn ha_flux ha_flux_err ha_eqw ha_eq_err s2_flux s2_flux_err s2_eqw s2_eq_err\n')
        fp.write('dummy 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n')
        fp.close()
    
    log = catIO.Readfile('simspec.dat')
    
    for ii, file in enumerate(files):
        root = file.split('.linefit')[0]
        print unicorn.noNewLine+'%s (%d/%d)' %(root, ii+1, len(files))
        if root in log.object:
            continue
        #
        fp = open('simspec.dat','a')
        pointing = root.split('_')[0]
        id = int(root.split('_')[1])
        if cat is None:
            cat = threedhst.sex.mySexCat(pointing+'_inter.cat')
            ### Get sky background
            asn = threedhst.utils.ASNFile(pointing+'-G141_asn.fits')
            bg = []
            for exp in asn.exposures:
                flt = pyfits.open(exp+'_flt.fits')
                bg.append(flt[0].header['SKYSCALE'])
            #
            bg_avg = np.mean(bg)
            bg_lo = np.min(bg)
            bg_hi = np.max(bg)
        else:
            if not cat.filename.startswith(pointing+'-'):
                cat = threedhst.sex.mySexCat(pointing+'_inter.cat')
                asn = threedhst.utils.ASNFile(pointing+'-G141_asn.fits')
                bg = []
                for exp in asn.exposures:
                    flt = pyfits.open(exp+'_flt.fits')
                    bg.append(flt[0].header['SKYSCALE'])
                #
                bg_avg = np.mean(bg)
                bg_lo = np.min(bg)
                bg_hi = np.max(bg)
        #        
        gris = unicorn.interlace_fit.GrismSpectrumFit(root, verbose=False)
        if not gris.status:
            fp.close()
            continue
        #
        result = gris.stats()
        if result is False:
            fp.close()
            continue
        #
        DIRECT_MAG, Q_Z, F_COVER, F_FLAGGED, MAX_CONTAM, INT_CONTAM, F_NEGATIVE = result
        #
        lwindow = (gris.oned.data.wave > 1.4e4) & (gris.oned.data.wave < 1.6e4)
        if (lwindow.sum() < 10) | (INT_CONTAM > 0.3):
            fp.close()
            continue
        #
        continuum_sn = np.median((gris.oned.data.flux/gris.oned.data.error)[lwindow])
        #
        lfit = catIO.Readfile(root+'.linefit.dat')
        if lfit.status is None:
            fp.close()
            continue
        #
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
        ic = np.arange(cat.nrows)[cat.id == id][0]
        fp.write(' %s  %5.2f %5.2f %5.2f  %6.3f  %6.2f %6.2f %6.4f %6.2f  %6.2f %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f %6.2f\n' %(root, bg_avg, bg_lo, bg_hi, DIRECT_MAG, float(cat.FLUX_RADIUS[ic]), float(cat.FLUX_RADIUS2[ic]), gris.z_max_spec, continuum_sn, ha_flux, ha_flux_err, ha_eqw, ha_eqw_err, s2_flux, s2_flux_err, s2_eqw, s2_eqw_err))
        #
        fp.close()
    
def show_results(use_tex=False):
    import threedhst.catIO as catIO
    
    stats = catIO.Readfile('all_simspec.dat')
    ha_model, s2_model = unicorn.intersim.get_line_fluxes(z0=1.0, mag=stats.mag)
    
    xstar = [14.5, 24.1]
    ystar = [3.00, 2.13]
    yi = np.interp(stats.mag, xstar, ystar)
    #plt.scatter(stats.mag, yi, s=0.1, color='black')
    is_star = stats.r50 < yi
    plt.scatter(stats.mag[is_star], stats.r50[is_star], alpha=0.5)
    plt.scatter(stats.mag[~is_star], stats.r50[~is_star], alpha=0.2, color='red')
    
    #### Color by r50/r90 concentration
    concentration = stats.r50/stats.r90
    msize = np.maximum((concentration/0.2)**4,4)
    mcol = np.minimum((np.maximum(concentration,0.3)-0.3)/0.2,1)
    plt.scatter(stats.mag, concentration, c=mcol, alpha=0.5)
    
    mcol = np.minimum(np.log10(stats.r50-1.1),1)
    
    stats.sky_avg += np.random.normal(size=stats.sky_avg.shape)*0.01
    sky_col = np.minimum((stats.sky_avg - 0.8)/0.8,1)
    plt.scatter(stats.mag, stats.sky_avg, c=sky_col, alpha=0.5)
    
    #### Continuum depth
    
    BINWIDTH=92
    bin_sn = np.sqrt(BINWIDTH/22)
    binned = stats.continuum_sn*bin_sn
    
    #### Get correction functions
    xm, ym, ys, nn = threedhst.utils.runmed(stats.mag, binned, NBIN=80)
    ymag = np.interp(stats.mag, xm, ym)

    sub = (stats.mag > 19) & (stats.mag < 22.5) & (stats.continuum_sn > 0) & (stats.ha_flux > 0) #& (~is_star)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.r50[sub], (binned/ymag)[sub], NBIN=20)
    ysize = np.interp(stats.r50, xm, ym)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.sky_avg[sub], (binned/ymag/ysize)[sub], NBIN=25)
    ysky = np.interp(stats.sky_avg, xm, ym)
    
    xm, ym, ys, nn = threedhst.utils.runmed(concentration[sub], (binned/ymag/ysize/ysky)[sub], NBIN=25)
    ycons = np.interp(concentration, xm, ym)
    
    fig = unicorn.catalogs.plot_init(xs=8, aspect=1./4, left=0.07, use_tex=use_tex)
    #fig.subplots_adjust(wspace=0.27, hspace=0.25, left=0.12)  # 2x2
    fig.subplots_adjust(wspace=0.38, hspace=0.25, left=0.074, bottom=0.22)
        
    si = 4
    mark = 'o'
    cmap = cm.jet
    bins = [80,80]
    
    ax = fig.add_subplot(141)
    
    #plt.scatter(stats.mag, stats.continuum_sn*bin_sn, alpha=0.5, c=mcol)
    use = np.isfinite(binned) & (binned > 0)
    #plt.scatter(stats.mag[use], (binned/ysize/ysky)[use], alpha=0.5, c=mcol[use], s=si, marker=mark)
    unicorn.intersim.show_hist_contour(stats.mag[use], (binned/ysize/ysky)[use], axrange=[[20,24],[0.5,100]], ylog=True, cmap=cmap, bins=bins)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.mag[use], (binned/ysize/ysky)[use], NBIN=80)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.ylim(0.5,100)
    plt.plot([20,24],[5,5], color='black', alpha=0.4)
    plt.xlim(20,24)
    plt.semilogy()
    if use_tex:
        plt.xlabel(r'MAG\_AUTO $m_{140}$')
    else:
        plt.xlabel(r'MAG_AUTO $m_{140}$')
        
    plt.ylabel('continuum S/N')
    ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(6, integer=True))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_yticks([1,10,100]); ax.set_yticklabels(['1','10','100'])
    
    sn5_limit = np.interp(5,ym[::-1],xm[::-1])
    print 'Continuum, S/N=5 @ %.3f' %(sn5_limit)
    print threedhst.utils.biweight(stats.r50[sub], both=True)
    
    ax = fig.add_subplot(142)
    
    #plt.scatter(stats.r50[sub], (binned/ymag/ysky)[sub], c=mcol[sub], alpha=0.5, s=si)
    unicorn.intersim.show_hist_contour(stats.r50[sub]*0.06, (binned/ymag/ysky)[sub], axrange=[[0,20*0.06],[0.3,1.7]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(stats.r50[sub]*0.06, (binned/ymag/ysky)[sub], NBIN=20)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.fill_betweenx([0,10],[1.7*0.06,1.7*0.06],[2.5*0.06,2.5*0.06], alpha=0.15, color='black')
    #plt.xlabel(r'$R_{50}$ [$0.\!\!^{\prime\prime}06$ pix]')
    plt.xlabel(r'$R_{50}$ [arcsec]')
    plt.ylabel(r'$\delta$ cont. S/N')
    plt.ylim(0.3,1.7)
    #plt.ylim(0.3,2.5)
    plt.xlim(0,15*0.06)
    majorLocator   = MultipleLocator(0.2)
    minorLocator   = MultipleLocator(0.1)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    
    # x0 = np.interp(1,ym[::-1],xm[::-1])
    # plt.plot(xm,(x0/xm), color='red') 
    # plt.plot(xm,(x0/xm)**0.5, color='red')
    
    x0 = np.interp(1,ym[::-1],xm[::-1])
    plt.plot(xm,(x0/xm)**(0.5), color='white', alpha=0.5, linewidth=2)
    plt.plot(xm,(x0/xm)**(0.5), color='red', alpha=0.8)
     
    ysize = np.interp(stats.r50*0.06, xm, ym)
    
    # plt.scatter(stats.r50[sub], (binned/ymag/ysize)[sub], c=sky_col[sub], alpha=0.5)
    # xm, ym, ys, nn = threedhst.utils.runmed(stats.r50[sub], (binned/ymag/ysize)[sub], NBIN=10)
    # plt.plot(xm, ym, linewidth=2, color='black', alpha=0.5)
    
    ax = fig.add_subplot(143)
    
    #plt.scatter(stats.sky_avg[sub], (binned/ymag/ysize)[sub], c=mcol[sub], alpha=0.5, s=si)
    unicorn.intersim.show_hist_contour(stats.sky_avg[sub], (binned/ymag/ysize)[sub], axrange=[[0.5,3.5],[0.3,1.7]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(stats.sky_avg[sub], (binned/ymag/ysize)[sub], NBIN=25)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.ylim(0.3,1.7)
    plt.xlim(0.5,3.5)
    plt.xlabel(r'Background [e$^-$ / s]')
    plt.ylabel(r'$\delta$ cont S/N')
    ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(6, integer=True))
    
    x0 = np.interp(1,ym[::-1],xm[::-1])
    plt.plot(xm,(x0/xm)**(0.5), color='white', alpha=0.5, linewidth=2)
    plt.plot(xm,(x0/xm)**(0.5), color='red', alpha=0.7)
    
    ysky = np.interp(stats.sky_avg, xm, ym)
    
    ### Very little residual trend with concentration
    ax = fig.add_subplot(144)
    
    #plt.scatter(concentration[sub], (binned/ymag/ysize/ysky)[sub], c=mcol[sub], s=si, alpha=0.5)
    unicorn.intersim.show_hist_contour(concentration[sub], (binned/ymag/ysize/ysky)[sub], axrange=[[0.25,0.60],[0.3,1.7]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(concentration[sub], (binned/ymag/ysize/ysky)[sub], NBIN=25)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.xlim(0.25,0.60)
    plt.ylim(0.3,1.7)
    #plt.ylim(0.5,1.5)
    plt.xlabel(r'$C = R_{50}/R_{90}$')
    plt.ylabel(r'$\delta$ cont S/N')
    #ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(5, prune=None))
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    
    ycons = np.interp(concentration, xm, ym)
    
    plt.savefig('grism_cont_sensitivity.pdf')
    
    # #### Test
    # plt.scatter(stats.mag, binned, alpha=0.5, c=sky_col, s=4)
    # xm, ym, ys, nn = threedhst.utils.runmed(stats.mag, binned, NBIN=80)
    # plt.errorbar(xm, ym, ys, linewidth=2, color='black', alpha=0.5)
    # plt.ylim(0.1,2000)
    # plt.plot([17,24],[5,5], color='black', alpha=0.4)
    # plt.xlim(17,24)
    # plt.semilogy()
    
    #### Line fluxes
    ha_sn = stats.ha_flux/stats.ha_flux_err
    
    show = np.isfinite(ha_sn) & (ha_sn > 0) & (stats.ha_flux > 0)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.ha_flux[~is_star & show], ha_sn[~is_star & show], NBIN=25)
    yline_flux = np.interp(stats.ha_flux, xm, ym)
    
    #sub = (stats.ha_flux > 6) & (stats.ha_flux < 100) & (stats.mag > 18) & (np.isfinite(ha_sn)) # & (~is_star)
    #sub = (stats.mag > 19) & (stats.mag < 22.5) & (stats.continuum_sn > 0) & (stats.ha_flux > 0) #& (~is_star)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.r50[sub], (ha_sn/yline_flux)[sub], NBIN=30)
    yline_r50 = np.interp(stats.r50, xm, ym)
    
    xm, ym, ys, nn = threedhst.utils.runmed(stats.sky_avg[sub], (ha_sn/yline_flux/yline_r50)[sub], NBIN=20)
    yline_sky = np.interp(stats.sky_avg, xm, ym)
    
    xm, ym, ys, nn = threedhst.utils.runmed(concentration[sub], (ha_sn/yline_flux/yline_r50/yline_sky)[sub], NBIN=10)
    yline_con = np.interp(concentration, xm, ym)
    
    plt.errorbar(ha_model, stats.ha_flux, stats.ha_flux_err, marker='o', markersize=0.1, linestyle='None', color='0.5')
    plt.scatter(ha_model, stats.ha_flux, c=mcol, zorder=100, alpha=0.5)
    #plt.scatter(stats.s2_flux, s2_model, alpha=0.8, c=mc)
    plt.plot([0.1,1000],[0.1,1000], color='black', alpha=0.5)
    plt.xlim(0.5,1000)
    plt.ylim(0.5,1000)
    plt.loglog()
     
    # 2x2   
    #fig = unicorn.catalogs.plot_init(xs=5.5, aspect=1, left=0.08)
    #fig.subplots_adjust(wspace=0.27, hspace=0.25, left=0.12)
    fig = unicorn.catalogs.plot_init(xs=8, aspect=1./4, left=0.07, use_tex=use_tex)
    fig.subplots_adjust(wspace=0.38, hspace=0.25, left=0.074, bottom=0.22)
    
    ax = fig.add_subplot(141)
    
    si = 4
    
    show = np.isfinite(ha_sn) & (ha_sn > 0) & (stats.ha_flux > 0)
    #plt.scatter(stats.ha_flux[show], ha_sn[show], c=mcol[show], s=si, zorder=100, alpha=0.3)
    unicorn.intersim.show_hist_contour(stats.ha_flux[show], (ha_sn/yline_r50/yline_sky/yline_con)[show], axrange=[[0.5,100],[0.5,100]], bins=bins, cmap=cmap, xlog=True, ylog=True)
    xm, ym, ys, nn = threedhst.utils.runmed(stats.ha_flux[~is_star & show], (ha_sn/yline_r50/yline_sky/yline_con)[~is_star & show], NBIN=25)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.plot([0.5,100],[5,5], color='black', alpha=0.4)
    plt.xlim(0.5,100)
    plt.ylim(0.5,100)
    plt.loglog()
    plt.xlabel(r'line flux [$10^{-17}$ ergs / s / cm$^2$]')
    plt.ylabel('line S/N')
    
    ax.set_yticks([1,10,100]); ax.set_yticklabels(['1','10','100'])
    ax.set_xticks([1,10,100]); ax.set_xticklabels(['1','10','100'])
    
    sn5_limit = np.interp(5,ym,xm)
    print 'Line, S/N=5 @ %.3e' %(sn5_limit)
    print threedhst.utils.biweight(stats.r50[sub], both=True)
    
    yline_flux = np.interp(stats.ha_flux, xm, ym)
    #plt.scatter(stats.ha_flux, ha_sn/yline_flux, c=mcol, alpha=0.2)

    #### Nice:  line flux with respect to concentration after taking out the overall trend with
    #### line strength
    ax = fig.add_subplot(142)
    
    #plt.scatter(stats.r50[sub], (ha_sn/yline_flux)[sub], c=mcol[sub], s=si, alpha=0.3)
    unicorn.intersim.show_hist_contour(stats.r50[sub]*0.06, (ha_sn/yline_flux/yline_sky/yline_con)[sub], axrange=[[0,15*0.06],[0.3,2.5]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(stats.r50[sub]*0.06, (ha_sn/yline_flux/yline_sky/yline_con)[sub], NBIN=30)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20*0.06],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.fill_betweenx([0,10],[1.7*0.06,1.7*0.06],[2.5*0.06,2.5*0.06], alpha=0.15, color='black')
    plt.ylim(0.3,2.5)
    plt.xlim(0,15*0.06)
    #plt.xlabel(r'$R_{50}$ [$0.\!\!^{\prime\prime}06$ pix]')
    plt.ylabel(r'$\delta$ line S/N')
    #plt.semilogy()
    # x0 = np.interp(1,ym[::-1],xm[::-1])
    # plt.plot(xm,(x0/xm), color='red') 
    # plt.plot(xm,(x0/xm)**0.5, color='red')
    
    plt.xlabel(r'$R_{50}$ [arcsec]')
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    
    x0 = np.interp(1,ym[::-1],xm[::-1])
    plt.plot(xm,(x0/xm)**(0.5), color='red', alpha=0.7)
    
    yline_r50 = np.interp(stats.r50*0.06, xm, ym)
    
    ax = fig.add_subplot(143)
    
    #plt.scatter(stats.sky_avg[sub], (ha_sn/yline_flux/yline_r50)[sub], c=mcol[sub], s=si, alpha=0.3)
    unicorn.intersim.show_hist_contour(stats.sky_avg[sub], (ha_sn/yline_flux/yline_r50/yline_con)[sub], axrange=[[0.5,3.5],[0.3,1.7]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(stats.sky_avg[sub], (ha_sn/yline_flux/yline_r50/yline_con)[sub], NBIN=20)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.ylim(0.3,1.7)
    plt.xlim(0.5,3.5)
    plt.xlabel(r'Background [e$^-$ / s]')
    plt.ylabel(r'$\delta$ line S/N')
    ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(6, integer=True))

    yline_sky = np.interp(stats.sky_avg, xm, ym)
    
    x0 = np.interp(1,ym[::-1],xm[::-1])
    plt.plot(xm,(x0/xm)**(0.5), color='red', alpha=0.7)
    
    ax = fig.add_subplot(144)

    #plt.scatter(concentration[sub], (ha_sn/yline_flux/yline_r50/yline_sky)[sub], c=mcol[sub], s=si, alpha=0.3)
    unicorn.intersim.show_hist_contour(concentration[sub], (ha_sn/yline_flux/yline_r50/yline_sky)[sub], axrange=[[0.25,0.60],[0.3,1.7]], bins=bins, cmap=cmap)
    xm, ym, ys, nn = threedhst.utils.runmed(concentration[sub], (ha_sn/yline_flux/yline_r50/yline_sky)[sub], NBIN=10)
    plt.plot(xm, ym, linewidth=2, color='white', alpha=0.5, zorder=100)
    plt.plot(xm, ym, linewidth=1, color='black', alpha=0.8, zorder=100)
    plt.plot([0,20],[1,1], linewidth=1, alpha=0.4, zorder=101, color='black')
    plt.xlim(0.25,0.60)
    plt.ylim(0.3,1.7)
    plt.xlabel(r'$C = R_{50}/R_{90}$')
    plt.ylabel(r'$\delta$ line S/N')
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    
    yline_con = np.interp(concentration, xm, ym)
    
    plt.savefig('grism_line_sensitivity.pdf')
    
    # #### Test:
    # show = (np.isfinite(ha_sn)) & (stats.ha_flux > 0)
    # plt.scatter(stats.ha_flux[show], (ha_sn/yline_sky)[show], c=mcol[show], zorder=100, alpha=0.2)
    # xm, ym, ys, nn = threedhst.utils.runmed(stats.ha_flux[show],  (ha_sn/yline_sky)[show], NBIN=25)
    # plt.plot(xm, ym, linewidth=2, color='black', alpha=0.5, zorder=100)
    # plt.plot([0.5,1000],[5,5], color='black', alpha=0.4)
    # plt.xlim(0.5,1000)
    # plt.ylim(0.5,300)
    # plt.loglog()
    
    #plt.semilogy()
    
    #
    plt.scatter(stats.mag, stats.ha_flux, c=mcol, zorder=100, alpha=0.5)
    plt.ylim(0.1,5000)
    plt.semilogy()
    
    #### EQW 
    dha = stats.ha_eqw-130.
    hy, hx, hh = plt.hist(dha/stats.ha_eq_err, range=(-5,5), bins=50, alpha=0.7)
    threedhst.utils.biweight(dha/stats.ha_eq_err, both=True)
    
    #### redshift
    dz = (stats.z_fit-1)/2.
    plt.scatter(stats.mag, dz, c=mcol, alpha=0.5)
    
    plt.scatter(stats.ha_flux, dz, c=mcol, alpha=0.5)
    plt.xlim(0.1,5000)
    plt.semilogx()
    
    #### surface density
    mu = stats.mag-2*np.log(stats.r90*0.06)
    plt.scatter(stats.mag, mu, c=mcol)

def show_hist_contour(xin, yin, axrange=None, bins=[50,50], xlog=False, ylog=False, ax=None, Vbins=[2, 4, 8, 16, 32, 64, 128, 256, 512, 4096], cmap=cm.jet, fill=True, *args, **kwargs):
    import matplotlib.colors as co
    
    if xlog:
        xdata = np.log10(xin)
    else:
        xdata = xin
    
    if ylog:
        ydata = np.log10(yin)
    else:
        ydata = yin
    
    if axrange is None:
        axrange = [[np.min(xdata),np.max(xdata)],[np.min(ydata),np.max(ydata)]]
    
    if xlog:
        for i in range(2):
            axrange[0][i] = np.log10(axrange[0][i])

    if ylog:
        for i in range(2):
            axrange[1][i] = np.log10(axrange[1][i])
                        
    hist, xedge, yedge = np.histogram2d(xdata, ydata, bins=bins, range=axrange)
    #Vbins = [2, 4, 8, 16, 32, 64, 128, 256, 512, 4096]
    values =   1.-np.arange(len(Vbins))*1./len(Vbins)
    Vcolors = []
    for i in range(len(Vbins)):
        Vcolors.append('%f' %(values[i]))
    
    if xlog:
        xx = 10**((xedge[:-1]+xedge[1:])/2.)
    else:
        xx = (xedge[:-1]+xedge[1:])/2.
    
    if ylog:
        yy = 10**((yedge[:-1]+yedge[1:])/2.)
    else:
        yy = (yedge[:-1]+yedge[1:])/2.
        
    norml = co.BoundaryNorm(Vbins, 312)
    
    if ax is None:
        if fill:
            plt.contourf(xx, yy, hist.transpose(), Vbins, linethick=2, norm=norml, cmap=cmap, *args, **kwargs)
        else:
            plt.contour(xx, yy, hist.transpose(), Vbins, linethick=2, norm=norml, cmap=cmap, *args, **kwargs)
    else:
        if fill:
            ax.contourf(xx, yy, hist.transpose(), Vbins, linethick=2, norm=norml, cmap=cmap, *args, **kwargs)
        else:
            ax.contour(xx, yy, hist.transpose(), Vbins, linethick=2, norm=norml, cmap=cmap, *args, **kwargs)
            
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
    
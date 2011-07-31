import os
import pyfits
import numpy as np
import glob
import shutil
import re

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from pyraf import iraf
from iraf import iraf

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

noNewLine = '\x1b[1A\x1b[1M'

zout = None
mcat = None
gfit = None
lines = None

def match_string_arrays(target=['b','a','d'], source=['a','b','c']):
    """
    
    Given a list of strings in `target`, return the integer indices of
    `source` where the strings match.  Fill with -1 if no match so that 
    the output array has the same size as `target`.
    
    Example:
    >>> target=np.array(['b','a','d'])
    >>> source=np.array(['a','b','c'])
    >>> found_in_source, indices = matche_string_arrays(target, source)
    >>> print target[found_in_source], source[indices]
    
    """
    
    indices = np.zeros(len(target), dtype='int')-1
    source_arr = np.array(source)
    for i,targ in enumerate(target):
        match = np.where(source_arr == targ)[0]
        if len(match) > 0:
            indices[i] = int(match[0])
    
    found_in_source = indices >= 0
    
    return found_in_source, indices
    
def test_plots():
    import unicorn.catalogs
    from unicorn.catalogs import match_string_arrays
    import cosmocalc
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/')
    
    ######################################
    #### Full photometry
    ######################################
    
    phot = catIO.Readfile('full_sextractor.cat')
    phot_goods = catIO.Readfile('goodsn_sextractor.cat')
    
    # Marijn's test of size vs mag
    ok = (phot.field == 'AEGIS') | (phot.field == 'COSMOS') | (phot.field == 'GOODS-S')
    
    xsep = np.array([14.6, 20.1, 21.7, 23.7, 24.6])
    ysep = np.array([9.6, 4.90, 3.11, 2.34, 2.14])
    yint = np.interp(phot.mag_f1392w, xsep, ysep)
    yint_goods = np.interp(phot_goods.mag_f1392w, xsep, ysep)
    
    is_star = phot.a_image < yint
    is_star_goods = phot_goods.a_image < yint_goods
    
    radius = np.sqrt(phot.a_image**2+phot.b_image**2)
    
    
    ## A_image
    plt.plot(xsep, ysep, color='red', linewidth=3, alpha=0.3)
    plt.plot(phot.mag_f1392w[ok & ~is_star], phot.a_image[ok & ~is_star], marker='.', alpha=0.08, linestyle='None', color='black', markersize=3)
    plt.plot(phot.mag_f1392w[ok & is_star], phot.a_image[ok & is_star], marker='.', alpha=0.2, linestyle='None', color='red', markersize=3)

    plt.plot(phot_goods.mag_f1392w[~is_star_goods], phot_goods.a_image[~is_star_goods], marker='.', alpha=0.08, linestyle='None', color='black', markersize=3)
    plt.plot(phot_goods.mag_f1392w[is_star_goods], phot_goods.a_image[is_star_goods], marker='.', alpha=0.2, linestyle='None', color='red', markersize=3)

    plt.xlim(14,24.8)
    plt.ylim(0,20)
    plt.xlabel(r'$m_\mathrm{140}$')
    plt.ylabel('A_IMAGE')
    
    plt.savefig('mag_vs_SEx_size.pdf')
    plt.savefig('mag_vs_SEx_size.png')
    
    ## Kron radius
    plt.plot(phot.mag_f1392w[ok & is_star], phot.kron_radius[ok & is_star], marker='o', alpha=0.05, linestyle='None', color='red', markersize=3)
    plt.plot(phot.mag_f1392w[ok & ~is_star], phot.kron_radius[ok & ~is_star], marker='o', alpha=0.05, linestyle='None', color='black', markersize=3)
    plt.xlim(14,24.8)
    plt.ylim(0,20)
    plt.xlabel(r'$m_\mathrm{140}$')
    plt.ylabel('KRON_RADIUS')
    
    ## a**2+b**2
    plt.plot(phot.mag_f1392w[ok & ~is_star], radius[ok & ~is_star], marker='o', alpha=0.05, linestyle='None', color='black', markersize=3)
    plt.plot(phot.mag_f1392w[ok & is_star], radius[ok & is_star], marker='o', alpha=0.05, linestyle='None', color='red', markersize=3)
    plt.xlim(14,24.8)
    plt.ylim(0,20)
    plt.xlabel(r'$m_\mathrm{140}$')
    plt.ylabel(r'$\mathrm{A\_IMAGE}^2+\ \mathrm{B\_IMAGE}^2$')
    
    ######################################
    #### Redshifts and masses
    ######################################
    
    zout = catIO.Readfile('full_redshift.zout')
    pointing = []
    field = []
    for obj in zout.id[0::3]:
        point = obj.split('-G141')[0]
        pointing.append(point)
        field.append(re.split('-[1-9]',point)[0])
    
    pointing = np.array(pointing)
    field = np.array(field)
    
    mcat = catIO.Readfile('full_match.cat')
        
    zsp = zout.z_spec[0::3] > 0
    dz = (zout.z_peak-zout.z_spec)/(1+zout.z_spec)
    found, idx = match_string_arrays(zout.id[0::3], mcat.id_f140w)
    idx = idx[idx >= 0]
    zsp = zsp & found

    ######################################
    #### Emission line fits
    ######################################
    
    lines = catIO.Readfile('full_emission_lines.cat')
    found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
    #lines.halpha_eqw[idx_lines] /= (1+zout.z_peak[0::3][found_lines])
    red_limit = 15
    red_sequence = -lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3]) < red_limit
    
    ##################################### 
    #### Galfit
    #####################################
    gfit = catIO.Readfile('full_galfit.cat')
    found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.object, zout.id[0::3])
    found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.object)
    
    #### Make a grid of redshifts to compute the plate scale and then interpolate it.
    zgrid = np.arange(100)/100.*4+1./100
    scale = zgrid*0.
    for i in range(100):
        cc = cosmocalc.cosmocalc(zgrid[i])
        scale[i] = cc['PS_kpc']
    
    gfit.r_e_kpc = gfit.r_e*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_err = gfit.r_e_err*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(np.abs(gfit.ba)) 
    gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(np.abs(gfit.ba)) 
    
    ##### Selection slices
    dr = mcat.rmatch[idx] < 1
    zrange = (zout.z_peak[0::3] > 0.2)
    
    keep = dr & zrange & (mcat.logm[idx] > 10.98) & (mcat.fcontam[idx] < 0.5)
    keep = dr & zrange & (mcat.mag_f140w[idx] < 24) & (mcat.fcontam[idx] < 0.5)
    
    #### refine
    # zrange = (zout.z_peak[0::3] > 2) & (zout.z_peak[0::3] < 3)
    # zrange = (zout.z_peak[0::3] > 1.5) & (zout.z_peak[0::3] < 2.)
    # zrange = (zout.z_peak[0::3] > 0.4) & (zout.z_peak[0::3] < 1.0)
    zrange = (zout.z_peak[0::3] > 1.0) & (zout.z_peak[0::3] < 1.5)
    keep = dr & zrange & (mcat.fcontam[idx] < 0.2) & (zout.q_z[0::3] < 0.1) 
    
    ### copy to macbook
    # copy = keep & (mcat.logm[idx] > 8.)
    # unicorn.catalogs.make_object_tarfiles(zout.id[0::3][copy])
    
    ##### Q_z vs dz
    plt.semilogx(zout.q_z[0::3][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='red', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(1.e-3, 10)
    plt.xlabel(r'$Q_z$')
    plt.ylabel(r'$\Delta z$')

    ##### z vs dz
    fig = plt.figure(figsize=(5.5,5), dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.15,
                        bottom=0.10,right=0.99,top=0.97)

    stats = abs(dz[0::3]) < 0.06 
    ## zphot
    ax = fig.add_subplot(211)
    
    dzmcat = (mcat.z_peak[idx]-zout.z_spec[0::3])/(1+zout.z_spec[0::3])
    plt.plot(zout.z_spec[0::3][zsp & keep], dzmcat[zsp & keep], marker='o', linestyle='None', color='orange', alpha=0.4)
    hist = np.histogram(dzmcat[zsp & keep], bins=100, range=(-0.05, 0.05))
    plt.plot(hist[0]*1./np.sum(hist[0])+0.85, hist[1][1:], linestyle='steps', alpha=0.8, color='orange')

    plt.ylim(-0.059, 0.059)
    plt.xlim(0.85, 1.7)
    plt.xlabel(r'$z_\mathrm{spec}$')
    plt.ylabel(r'$\Delta z$')
    plt.text(1.5,-0.05,r'$\sigma=%.3f$' %(threedhst.utils.biweight(dzmcat[zsp & keep & stats])), fontsize=16)
    
    ## gspec
    ax = fig.add_subplot(212)
    plt.plot(zout.z_spec[0::3][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='purple', alpha=0.4)
    hist = np.histogram(dz[0::3][zsp & keep], bins=100, range=(-0.05, 0.05))
    plt.plot(hist[0]*1./np.sum(hist[0])+0.85, hist[1][1:], linestyle='steps', alpha=0.8, color='purple')

    plt.ylim(-0.059, 0.059)
    plt.xlim(0.85, 1.7)
    plt.xlabel(r'$z_\mathrm{spec}$')
    plt.ylabel(r'$\Delta z$')
    plt.text(1.5,-0.05,r'$\sigma=%.3f$' %(threedhst.utils.biweight(dz[0::3][zsp & keep & stats])), fontsize=16)

    plt.savefig('redshift_errors.pdf')
    plt.savefig('redshift_errors.png')
    
    ##### F140W_mag vs dz
    plt.plot(mcat.mag_f140w[idx][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(15, 24)
    plt.xlabel(r'$m_\mathrm{f140w}$')
    plt.ylabel(r'$\Delta z$')
        
    ##### Contamination vs dz
    plt.plot(mcat.fcontam[idx][zsp & keep]+1.e-2, dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0, 0.5)
    plt.xlabel(r'$f_\mathrm{contam}$')
    plt.ylabel(r'$\Delta z$')
    
    ##### mass vs dz
    plt.plot(mcat.logm[idx][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(8, 12)
    plt.xlabel(r'$\log M/M_\odot$')
    plt.ylabel(r'$\Delta z$')

    ##### Redshift distribution
    unicorn.catalogs.plot_init()
    
    hist = np.histogram(zout.z_peak[0::3][keep], bins=100, range=(1,1.5))
    plt.plot((hist[1][:-1]+hist[1][1:])/2., hist[0], linestyle='steps', color='black', alpha=0.5)
    hist_wide = np.histogram(zout.z_peak[0::3][keep], bins=6, range=(0.9,1.5))
    plt.plot(hist_wide[1][1:], hist_wide[0]/20., linestyle='steps', color='black', alpha=0.5, linewidth=4)
    #plt.plot((hist_wide[1][:-1]+hist_wide[1][1:])/2., hist_wide[0]/20., linestyle='steps', color='black', alpha=0.5, linewidth=4)
    plt.xlim(1,1.5)
    plt.ylim(0,30)
    plt.xlabel(r'$z_\mathrm{grism}\ (\Delta z=0.005)$')
    plt.ylabel(r'$N\ \ (N_\mathrm{tot}=%d)$' %len(zout.z_peak[0::3][keep]))
   
    peak = np.abs(zout.z_peak[0::3][keep]-1.335) < 0.005
    zout.id[0::3][keep][peak]
    
    plt.savefig('redshift_dist.pdf')
    plt.savefig('redshift_dist.png')
    
    ##### line strength vs dz
    
    plt.semilogx(-lines.halpha_eqw[idx_lines][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0.1, 500)
    plt.xlabel(r'$\log M/M_\odot$')
    plt.ylabel(r'$\Delta z$')
    
    ##### plot mass vs equivalent width
    unicorn.catalogs.plot_init()
    
    plt.semilogy(mcat.logm[idx][~zsp & keep], -lines.halpha_eqw[idx_lines][~zsp & keep]/(1+zout.z_peak[0::3][~zsp & keep]), marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & keep], -lines.halpha_eqw[idx_lines][zsp & keep]/(1+zout.z_peak[0::3][zsp & keep]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    # oiii = keep & (zout.z_peak[0::3] > 1.2) & (-lines.oiii_eqw[idx_lines] > 25)
    # plt.semilogy(mcat.logm[idx][oiii], -lines.halpha_eqw[idx_lines][oiii], marker='s', linestyle='None', color='red', alpha=0.2, markersize=8)
    
    ## red sequence
    plt.plot([0,100],[red_limit,red_limit], linestyle='--', alpha=0.5, color='green', linewidth=4)

    plt.text(9.5,1,r'has $z_\mathrm{spec}$', color='blue', fontsize=18)
    
    ### For SDSS, read in catalogs and get selection from unicorn.paper1
    # sel, NOBJ = p1.sdss_selection(zmax=0.2)
    # hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], np.log10(-p1.sdss_line.H_ALPHA_EQW[sel]/(1+p1.sdss_info.Z[sel])), bins=100, range=[[9,11.6], [-1, 2.7]])
    # Vbins = [2, 4, 8, 16, 32, 64, 128, 256, 512, 4096]
    # values =   1.-np.arange(len(Vbins))*1./len(Vbins)
    # Vcolors = []
    # for i in range(len(Vbins)):
    #     Vcolors.append('%f' %(values[i]))
    # 
    # unicorn.catalogs.plot_init()
    # plt.contourf(xedge[1:], 10**yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    # plt.semilogy(mcat.logm[idx][keep], -lines.halpha_eqw[idx_lines][keep]/(1+zout.z_peak[0::3][keep]), marker='o', linestyle='None', color='red', alpha=0.5, markersize=5)
    f_init
    
    plt.ylim(0.1, 500)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{H}\alpha\ \mathrm{eqw}$')

    # plt.savefig('ha_mass_sdss.pdf')
    # plt.savefig('ha_mass_sdss.png')
    
    plt.savefig('ha_mass_zspec.pdf')
    plt.savefig('ha_mass_zspec.png')
    
    ####### OIII emitters
    oiii = keep & (zout.z_peak[0::3] > 1.2)
    plt.semilogy(mcat.logm[idx][oiii], -lines.halpha_eqw[idx_lines][oiii], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & oiii], -lines.halpha_eqw[idx_lines][zsp & oiii], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    plt.ylim(0.01, 500)
    plt.xlim(8.5, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{H}\alpha\ \mathrm{EQW}$')

    plt.semilogy(mcat.logm[idx][oiii], -lines.oiii_eqw[idx_lines][oiii], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & oiii], -lines.oiii_eqw[idx_lines][zsp & oiii], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    plt.ylim(0.01, 500)
    plt.xlim(8.5, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{OIII}\ \mathrm{eqw}$')
    
    #### Mass - size
    # plt.semilogy(mcat.logm[idx][keep], gfit.r_e_kpc[idx_gfit][keep], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.semilogy(mcat.logm[idx][keep & ~red_sequence], gfit.r_e_kpc[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.3, markersize=8)
    plt.semilogy(mcat.logm[idx][keep & red_sequence], gfit.r_e_kpc[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.3, markersize=8)
    
    # circularized
    unicorn.catalogs.plot_init(square=True)
    
    plt.semilogy(mcat.logm[idx][keep & ~red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.1, markersize=8)
    plt.semilogy(mcat.logm[idx][keep & red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.1, markersize=8)
    
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep], gfit.r_e_kpc_circ[idx_gfit][keep], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='black', ecolor='black', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)
    
    vd08 = catIO.Readfile('vd08.dat')
    plt.plot(vd08.logm, 10**vd08.log_re, marker='o', color='orange', alpha=0.7, markersize=12, linestyle='None')

    vd08 = catIO.Readfile('vd08_sdss.dat')
    plt.plot(vd08.logm, 10**vd08.log_re, marker='None', color='green', alpha=0.3, linewidth=14)
    
    plt.ylim(0.6, 20)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$r_\mathrm{e}\ \mathrm{(kpc)}$')
    plt.savefig('mass_size.pdf')
    plt.savefig('mass_size.png')
    
    #### Mass - n
    unicorn.catalogs.plot_init()
    
    plt.semilogy(mcat.logm[idx][keep & ~red_sequence], gfit.n[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.1, markersize=8)
    plt.semilogy(mcat.logm[idx][keep & red_sequence], gfit.n[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.1, markersize=8)

    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep], gfit.n[idx_gfit][keep], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='black', ecolor='black', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.n[idx_gfit][keep & ~red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.n[idx_gfit][keep & red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)

    plt.plot([0,100],[1,1], color='blue', linestyle='--', alpha=0.4, linewidth=4)
    plt.plot([0,100],[4,4], color='red', linestyle='--', alpha=0.4, linewidth=4)
    plt.ylim(0.2, 21)
    plt.xlim(9, 11.6)
    plt.text(9.2,1.1,r'$n=1$', fontsize=14)
    plt.text(9.2,4.4,r'$n=4$', fontsize=14)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$n$')
    
    plt.savefig('mass_n.pdf')
    plt.savefig('mass_n.png')
    
    #### Mass - b/a
    unicorn.catalogs.plot_init()
    
    plt.plot(mcat.logm[idx][keep & ~red_sequence], gfit.ba[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(mcat.logm[idx][keep & red_sequence], gfit.ba[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)

    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.ba[idx_gfit][keep & ~red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.ba[idx_gfit][keep & red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)

    plt.ylim(0., 1)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$b/a$')
    
    plt.savefig('mass_ba.pdf')
    plt.savefig('mass_ba.png')
        
    ######### Composite spectra
    eqw = -lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3])
    
    mass = keep & (mcat.logm[idx] > 10.6)
    
    eqws = eqw[mass]
    eqws.sort()
    quartiles = eqws[np.cast[int](np.round(np.array([0.25,0.5,0.75])*len(eqws)))]
    thirds = eqws[np.cast[int](np.round(np.array([1./3, 2./3])*len(eqws)))]
    
    high_eqw = (eqw > 23.7) & mass
    med_eqw = (eqw > 9.8) & (eqw <= 23.7) & mass
    low_eqw = (eqw <= 9.8) & mass
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    unicorn.catalogs.composite_spectra(lines.id[idx_lines][high_eqw], color='blue', alpha=0.02, lnorm=6.e3, NITER=3, show_lines=True)
    plt.text(7700,1.8,r'$N_\mathrm{EW>24}=%d$' %(len(lines.id[idx_lines][high_eqw])), color='blue', horizontalalignment='right')
    plt.ylim(0.7,2)
    plt.xlim(4300, 7900)

    #unicorn.catalogs.plot_init()
    unicorn.catalogs.composite_spectra(lines.id[idx_lines][med_eqw], color='green', alpha=0.02, lnorm=6.e3, NITER=3, show_lines=False)
    plt.text(7700,1.7,r'$N_\mathrm{10<EW<24}=%d$' %(len(lines.id[idx_lines][med_eqw])), color='green', horizontalalignment='right')
    plt.ylim(0.7,2)
    plt.xlim(4300, 7900)

    #unicorn.catalogs.plot_init()
    unicorn.catalogs.composite_spectra(lines.id[idx_lines][low_eqw], color='red', alpha=0.02, lnorm=6.e3, NITER=3, show_lines=False)
    plt.text(7700,1.6,r'$N_\mathrm{EW<10}=%d$' %(len(lines.id[idx_lines][low_eqw])), color='red', horizontalalignment='right')
    plt.ylim(0.7,2)
    plt.xlim(4300, 7900)
    
    plt.xlabel(r'$\lambda_\mathrm{rest}$')
    plt.ylabel(r'$f_\lambda$')
    
    plt.savefig('composite_spectra.pdf')
    plt.savefig('composite_spectra.png')
    
    ################################################
    #### Eq-W vs size
    eqw = -lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3])
    eqw[eqw < 0.2] = 0.2
    mass_limit = 10.6
    massive = mcat.logm[idx] > mass_limit
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    plt.plot(mcat.logm[idx][keep & ~red_sequence & massive], np.log10(eqw[keep & ~red_sequence & massive]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(mcat.logm[idx][keep & red_sequence & massive], np.log10(eqw[keep & red_sequence & massive]), marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)
    
    plt.plot([-10,-10],[-10,-10])
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.xlim(9, 11.6)
    x = mcat.logm[idx][keep & massive]
    plt.ylim(-1, np.log10(200))
    y = np.log10(eqw[keep & massive])
    
    plt.savefig('eqw_mass9_morph_A.pdf')
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    plt.plot(np.log10(eqw[keep & ~red_sequence & massive]), np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence & massive]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(np.log10(eqw[keep & red_sequence & massive]), np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence & massive]), marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)

    plt.xlim(-1, np.log10(200))
    plt.ylim(-0.2, 1.5)
    plt.xlabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.ylabel(r'$\log\ r_\mathrm{e}\ (\mathrm{kpc})$')
    y = np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & massive])
    x = np.log10(eqw[keep & massive])
    
    plt.savefig('eqw_size_morph_A.pdf')
        
    objects = mcat.id_f140w[idx][keep & massive]
    xrange = plt.xlim()
    yrange = plt.ylim()
    scale = 1./10**(-0.4*(mcat.mag_f140w[idx][keep & massive]-22))
    
    #scale = scale*0+1
    
    im = unicorn.catalogs.fill_image(objects, x, y, scale=scale, xrange=xrange, yrange=yrange, NX=1024*1., NY=1024*1.)
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    plt.imshow(0-im, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), aspect='auto', interpolation='nearest', vmin=-0.5, vmax=0.05)

    plt.xlabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.ylabel(r'$\log\ r_\mathrm{e}\ (\mathrm{kpc})$')
    plt.savefig('eqw_size_morph_B.pdf')

    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\log\ \mathrm{EQW\ H}\alpha$')    
    plt.savefig('eqw_mass9_morph_B.pdf')
    
def fill_image(objects, x, y, scale=None, xrange=(0,1), yrange=(0,1), NX=1024, NY=1024):
    
    im = np.zeros((NY, NX))
    NOBJ = len(objects)
    dx = xrange[1]-xrange[0]
    dy = yrange[1]-yrange[0]
    
    if scale is None:
        scale = np.ones(NOBJ)
        
    for i in range(NOBJ):
        print noNewLine+objects[i]
        thu = pyfits.open('DATA/'+objects[i]+'_thumb.fits.gz')
        shp = thu[0].data.shape
        
        xi = np.round((x[i]-xrange[0])/dx*NX)-shp[1]/2
        yi = np.round((y[i]-yrange[0])/dy*NY)-shp[0]/2
        
        x0, x1 = 0, shp[1]
        y0, y1 = 0, shp[0]
        
        if (xi+x1/2 < 0) | (xi > NX) | (yi+y1/2 < 0) | (yi > NY):
            print '\n', x[i], y[i], shp, '\n'
            continue
            
        if xi < 0:
            x0 = -xi
            xi = 0
        if (xi+shp[1]) > NX:
            x1 = NX-xi
            
        if yi < 0:
            y0 = -yi
            yi = 0
        if (yi+shp[0]) > NY:
            y1 = NY-yi
        
        print x0,x1,y0,y1,shp, xi, yi, scale[i]
        
        thu = thu[0].data[y0:y1, x0:x1]*scale[i] #/thu[0].data.sum()
        shp = thu.shape
        im[yi:yi+shp[0], xi:xi+shp[1]] += thu 
    
    return im
    
def composite_spectra(objects, color='red', alpha=0.1, lnorm=8.e3, NITER=3, show_lines=False):
    
    zout = catIO.Readfile('full_redshift.zout')
    
    #objects = ['GOODS-S-24-G141_00459','GOODS-N-27-G141_00447']
    xm = np.array([3000,8000])
    ym = np.array([1,1])
    
    ### Iterate on normalizing to the average spectrum.
    for iter in range(1,NITER+1):
      avgx = []
      avgy = []
      for object in objects:
        print noNewLine+object
        zi = zout.z_peak[0::3][zout.id[0::3] == object][0]
        #
        try:
            data = np.loadtxt('DATA/'+object+'_scaled.dat')
        except:
            continue
        #           
        lc = data[:,0]
        fnu = data[:,1]
        efnu = data[:,2]
        #
        dlam_spec = lc[-1]-lc[-2]
        is_spec = np.append(np.abs(1-np.abs(lc[1:]-lc[0:-1])/dlam_spec) < 0.05,True)
        #
        flam = fnu/(lc/5500.)**2    
        s = np.argsort(lc)
        fint = np.interp(lnorm, lc[s]/(1+zi), flam[s])
        
        mask = is_spec & (np.abs(lc/(1+zi)-6563.) > 200)
        ymint = np.interp(lc[mask]/(1+zi), xm, ym)
        fnorm = np.sum(flam[mask]*ymint*efnu[mask]**2)/np.sum(ymint**2*efnu[mask]**2)
        flam /= fnorm
        eflam = efnu/(lc/5500.)**2/fnorm
        #
        if iter == NITER:
            plt.plot(lc[is_spec]/(1+zi), flam[is_spec], color=color, alpha=alpha)
            #plt.plot(lc[~is_spec]/(1+zi), flam[~is_spec], color=color, alpha=alpha, marker='o', markersize=5, linestyle='None')
        #
        avgx.extend(lc[is_spec]/(1+zi))
        avgy.extend(flam[is_spec])
    
      avgx = np.array(avgx)
      avgy = np.array(avgy)
      xm, ym, ys, N = threedhst.utils.runmed(avgx, avgy, NBIN=100)
      if iter == NITER:
          plt.plot(xm, ym, color='white', alpha=0.5, linewidth=6)
          plt.plot(xm, ym, color=color, alpha=8, linewidth=3)
    
    if show_lines:
        for line in [4861, 4959, 5007, 5178, 5891, 6563, 6585, 6718, 6731]:
            plt.plot(line*np.array([1,1]), [0,10], color='black', linestyle='--', alpha=0.3)
            
def plot_init(square=True, xs=6, aspect=1):
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = ['Times']
    plt.rcParams['patch.edgecolor'] = 'None'
    plt.rcParams['font.size'] = 10

    if square:
        #xs=5
        lrbt = np.array([0.22,0.02,0.11,0.02])*5./xs     
        ys = (1-lrbt[1]-lrbt[0])/(1-lrbt[3]-lrbt[2])*xs*aspect
        lrbt[[2,3]] /= aspect
        fig = plt.figure(figsize=(xs,ys), dpi=100)
        #fig.subplots_adjust(left=0.13,bottom=0.10,right=0.98,top=0.98)
        fig.subplots_adjust(left=lrbt[0],bottom=lrbt[2],right=1-lrbt[1],top=1-lrbt[3])
        # plt.plot([0,2])
        # plt.xlabel('x')
        # plt.ylabel('y')        
    else:
        fig = plt.figure(figsize=(7,5), dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.10,
                        bottom=0.10,right=0.99,top=0.97)        
    return fig
    
def make_object_tarfiles(objects, thumbs=True):
    """
    Make a script to get spectra and thumbnails from unicorn.
    """
    
    ###### Thumbnails
    if thumbs:
        os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.5/images')
        line = 'tar czvf thumbs.tar.gz'
        for object in objects:
            line += ' %s_thumb.fits.gz' %(object)
        
            os.system(line)
        
        os.system('mv thumbs.tar.gz /3DHST/Spectra/Work/ANALYSIS/FIRST_PAPER/')
    
    
    #### Get photometry + scaled spectra from the tempfilt files
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FIRST_PAPER')
    
    tarline = 'tar czvf spec.tar.gz'
    for object in objects:
        print noNewLine+object
        #
        tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY='../REDSHIFT_FITS/OUTPUT', CACHE_FILE = 'Same')
        
        #### Photometry + spectrum
        fp = open(object+'_obs_sed.dat','w')
        fp.write('# lc fnu efnu obs_sed \n')

        obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][0]],\
                         coeffs['coeffs'][:,0])# /(lci/5500.)**2
        
        for i in range(tempfilt['NFILT']):
            fp.write(' %8.4e %8.4e %8.4e  %8.4e\n' %(tempfilt['lc'][i], tempfilt['fnu'][i,0], tempfilt['efnu'][i,0], obs_sed[i]))
        
        fp.close()
        
        tarline += ' %s_obs_sed.dat' %(object)
        
        #### template fit
        zi = tempfilt['zgrid'][coeffs['izbest'][0]]
        lambdaz = temp_seds['templam']*(1+zi)
        temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,0])
        temp_sed *= (lambdaz/5500.)**2/(1+zi)**2
        
        fp = open(object+'_temp_sed.dat','w')
        fp.write('# lam fnu_temp\n')
        for i in range(temp_seds['NTEMPL']):
            fp.write(' %8.4e %8.4e\n' %(lambdaz[i], temp_sed[i]))
        fp.close()
        tarline += ' %s_temp_sed.dat' %(object)
        
        #
    os.system(tarline)
    
    print 'scp $UNICORN:'+unicorn.GRISM_HOME+'ANALYSIS/FIRST_PAPER/[st]*tar.gz .'
    
def combine_sextractor_catalogs():
    os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/')

    files=glob.glob('*drz.cat')
    
    fp = open(files[0])
    lines = fp.readlines()
    fp.close()
    
    line = lines[0]
    columns = []
    i=0
    while line.startswith('#'):
        columns.append(line.split()[2])
        i+=1
        line = lines[i]
    
    header = '# ID FIELD POINTING '+' '.join(columns)+'\n'
    out_lines = []
    for file in files:
        print noNewLine+file
        if file.startswith('GN20'):
            continue
        ### GOODS-N catalogs don't have the MAG_APER columns
        if file.startswith('GOODS-N'):
            continue
        #
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        pointing = file.split('-G141')[0]
        field = re.split('-[1-9]',pointing)[0]
        if len(pointing.split(field+'-')) == 1:
            pointing_number == '1'
        else:
            pointing_number = pointing.split(field+'-')[1]
        #
        for line in lines:
            if not line.startswith('#'):
                number = int(line.split()[0])
                out_lines.append('%s-G141_%05d %s %3s ' %(pointing, number, field, pointing_number)+line)
    
    fp = open('ANALYSIS/full_sextractor.cat','w')
    fp.write(header)
    fp.writelines(out_lines)
    fp.close()
    
    ########## GOODS-N separate for now because doesn't have same columns
    files=glob.glob('GOODS-N*drz.cat')
    
    fp = open(files[0])
    lines = fp.readlines()
    fp.close()
    
    line = lines[0]
    columns = []
    i=0
    while line.startswith('#'):
        columns.append(line.split()[2])
        i+=1
        line = lines[i]
    
    header = '# ID FIELD POINTING '+' '.join(columns)+'\n'
    out_lines = []
    for file in files:
        print noNewLine+file
        if file.startswith('GN20'):
            continue
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        pointing = file.split('-G141')[0]
        field = re.split('-[1-9]',pointing)[0]
        if len(pointing.split(field+'-')) == 1:
            pointing_number = '1'
        else:
            pointing_number = pointing.split(field+'-')[1]
        #
        for line in lines:
            if not line.startswith('#'):
                number = int(line.split()[0])
                out_lines.append('%s-G141_%05d %s %3s ' %(pointing, number, field, pointing_number)+line)
    
    fp = open('ANALYSIS/goodsn_sextractor.cat','w')
    fp.write(header)
    fp.writelines(out_lines)
    fp.close()

def make_full_redshift_catalog():
    """
    Cat all individual zout files into a single file
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    files = glob.glob('OUTPUT/*G141*zout')
    if len(files) == 0:
        os.system("ls OUTPUT/ |grep G141 |grep zout |awk '{print \"x\" $1 }' |sed \"s/x/OUTPUT\//\" > files.list")
        files = np.loadtxt('files.list', dtype=np.str)
        
    fp = open(files[0])
    lines = fp.readlines()
    fp.close()
    for file in files[1:]:
        print noNewLine+file
        fp = open(file)
        lines.extend(fp.readlines()[2:])
        fp.close()
        
    fp = open('full_redshift.zout','w')
    fp.writelines(lines)
    fp.close()
    
    status = os.system('cp full_redshift.zout /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
    status = os.system('gzip /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS/full_redshift.zout')
    
def combine_matched_catalogs():
    """
    Combine all of the match catalogs in the HTML/SED directories, adding
    the correct object id with the full pointing name.
    """
    
    os.chdir('/Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
    files = glob.glob('../SED/*match.cat')
    fp = open(files[0])
    full_lines = fp.readlines()[0:2]
    fp.close()
    for file in files:
        pointing = os.path.basename(file).split('_match')[0]
        print noNewLine+pointing
        fp = open(file)
        lines = fp.readlines()[3:]
        fp.close()
        for line in lines:
            spl = line.split()
            id = int(spl[0])
            object = "%s_%05d  " %(pointing, id)
            full_lines.append(object+'  '.join(spl[1:])+'\n')
    
    fp = open('full_match.cat','w')
    fp.writelines(full_lines)
    fp.close()
    
    status = os.system('gzip full_match.cat')
    

def show_acs_spectra():
    import unicorn.catalogs
    from unicorn.catalogs import match_string_arrays
    import cosmocalc

    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/')

    ######################################
    #### Redshifts and masses
    ######################################
    
    zout = catIO.Readfile('full_redshift.zout')
    pointing = []
    field = []
    for obj in zout.id[0::3]:
        point = obj.split('-G141')[0]
        pointing.append(point)
        field.append(re.split('-[1-9]',point)[0])
    
    pointing = np.array(pointing)
    field = np.array(field)
    
    # mcat = catIO.Readfile('full_match.cat')
    #     
    # zsp = zout.z_spec[0::3] > 0
    # dz = (zout.z_peak-zout.z_spec)/(1+zout.z_spec)
    # found, idx = match_string_arrays(zout.id[0::3], mcat.id_f140w)
    # idx = idx[idx >= 0]
    # zsp = zsp & found
    # 
    # ######################################
    # #### Emission line fits
    # ######################################
    # 
    # lines = catIO.Readfile('full_emission_lines.cat')
    # found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
    # #lines.halpha_eqw[idx_lines] /= (1+zout.z_peak[0::3][found_lines])
    # 
    # ##################################### 
    # #### Galfit
    # #####################################
    # gfit = catIO.Readfile('full_galfit.cat')
    # found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.object, zout.id[0::3])
    # found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.object)
    # 
    # #### Make a grid of redshifts to compute the plate scale and then interpolate it.
    # zgrid = np.arange(100)/100.*4+1./100
    # scale = zgrid*0.
    # for i in range(100):
    #     cc = cosmocalc.cosmocalc(zgrid[i])
    #     scale[i] = cc['PS_kpc']
    # 
    # gfit.r_e_kpc = gfit.r_e*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    # gfit.r_e_kpc_err = gfit.r_e_err*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    # gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(1./np.abs(gfit.ba)) 
    # gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(1./np.abs(gfit.ba)) 
    # 
    # ##### Selection slices
    # dr = mcat.rmatch[idx] < 1
    # zrange = (zout.z_peak[0::3] > 0.2)
    # 
    # keep = dr & zrange & (mcat.logm[idx] > 10.98) & (mcat.fcontam[idx] < 0.5)
    # keep = dr & zrange & (mcat.mag_f140w[idx] < 24) & (mcat.fcontam[idx] < 0.5)
    # 
    # #### refine
    # # zrange = (zout.z_peak[0::3] > 2) & (zout.z_peak[0::3] < 3)
    # # zrange = (zout.z_peak[0::3] > 1.5) & (zout.z_peak[0::3] < 2.)
    # # zrange = (zout.z_peak[0::3] > 0.4) & (zout.z_peak[0::3] < 1.0)
    # zrange = (zout.z_peak[0::3] > .0) & (zout.z_peak[0::3] < 1.5)
    # keep = dr & zrange & (mcat.fcontam[idx] < 0.2) & (zout.q_z[0::3] < 0.1) 
    
    #### Photometry with coords
    phot = catIO.Readfile('full_sextractor.cat')
    found_phot, idx_phot = match_string_arrays(zout.id[0::3], phot.id)
    
    #### ACS
    root = 'jbhm54020'
    acs = threedhst.sex.mySexCat('ACS/'+root+'_drz.cat')
    acs.ra = np.cast[float](acs.X_WORLD)
    acs.dec = np.cast[float](acs.Y_WORLD)

    acs_match = []
    acs_dr = np.arange(len(acs.NUMBER))+1.
    
    for i in range(len(acs.NUMBER)):
        dr = np.sqrt((acs.ra[i]-phot.x_world)**2*np.cos(acs.dec[i]/360.*2*np.pi)**2+(acs.dec[i]-phot.y_world)**2)*3600.
        ma = dr == dr.min()
        acs_match.append(phot.id[ma][0])
        acs_dr[i] = dr[ma][0]
    
    acs_match = np.array(acs_match)
    acs.mag = np.cast[float](acs.MAG_F806W)
    matches = (acs_dr < 1.0) & (acs.mag < 28)
    idx = np.arange(len(acs.mag))[matches] 
    
    cat, zcat, fout = unicorn.analysis.read_catalogs(root='COSMOS')
    
    OUTPUT_DIRECTORY = os.path.dirname(zcat.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zcat.filename).split('.zout')[0]
    
    for i0, i in enumerate(idx):
        i = idx[i0]
        ra0, de0 = acs.ra[i], acs.dec[i]
        dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
        photom_idx = np.where(dr == np.min(dr))[0][0]
        drMatch = dr[photom_idx]*1.
        #if drMatch > 0.5: 
        #    continue
        
        ####
        acs_file = '%s_%05d' %(root, acs.id[i])
        wfc_file = acs_match[matches][i0]
        
        if not os.path.exists(acs_file+'.dat'):
            os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_ACS/ascii/%s.dat' %(acs_file))
        
        if not os.path.exists(wfc_file+'.dat'):    
            os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.5/ascii/%s.dat' %(wfc_file))

        if (not os.path.exists(acs_file+'.dat')) | (not os.path.exists(wfc_file+'.dat')):
            continue
            
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(photom_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = 'Same')
    
        acs_spec = catIO.Readfile(acs_file+'.dat')
        wfc_spec = catIO.Readfile(wfc_file+'.dat')

        acs_spec.flux -= acs_spec.contam
        wfc_spec.flux -= wfc_spec.contam
        
        zi = zout.z_peak[::3][zout.id[::3] == wfc_file]
        
        if len(zi) == 0:
            continue
        else:
            zi = zi[0]
        
        keep_acs = (acs_spec.lam > 6000) & (acs_spec.lam < 9000) & (acs_spec.flux > 0)
        mask_acs = (np.abs(acs_spec.lam/(1+zi)-6563.) > 500) & (np.abs(acs_spec.lam/(1+zi)-5007.) > 500)
        
        keep_wfc = (wfc_spec.lam > 1.15e4) & (wfc_spec.lam < 1.62e4) & (wfc_spec.flux > 0)
        mask_wfc = (np.abs(wfc_spec.lam/(1+zi)-6563.) > 500) & (np.abs(wfc_spec.lam/(1+zi))-5007. > 500)
        
        yint_acs = np.interp(acs_spec.lam, lambdaz, temp_sed)
        norm_acs = np.sum((acs_spec.flux*yint_acs)[keep_acs & mask_acs])/np.sum(yint_acs[keep_acs & mask_acs]**2)

        yint_wfc = np.interp(wfc_spec.lam, lambdaz, temp_sed)
        norm_wfc = np.sum((wfc_spec.flux*yint_wfc)[keep_wfc & mask_wfc])/np.sum(yint_wfc[keep_wfc & mask_wfc]**2)
        
        xg = np.arange(1000)-500
        yg = np.exp(-xg**2/2/100**2)
        yg /= np.sum(yg)
        ll = np.arange(3.8e4)+2000
        yy = np.interp(ll, lambdaz, temp_sed)
        
        tconv = np.convolve(yy, yg, mode='same')
        
        fig = plt.figure(figsize=[5,3.7],dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.1,left=0.1,
                            bottom=0.15,right=0.98,top=0.98)
        
        ax = fig.add_subplot(111)
        
        zcorr = (1+zi)
        zcorr = 1.
        
        #plt.plot(lambdaz, temp_sed, color='black', alpha=0.4)
        ax.plot(ll/zcorr, tconv, color='black', alpha=0.3, linewidth=2)
        #plt.plot(ll, tconv, color='orange', alpha=0.4)
        ax.plot(wfc_spec.lam[keep_wfc]/zcorr, wfc_spec.flux[keep_wfc]/norm_wfc, color='red', linewidth=3, alpha=0.8)
        ax.plot(acs_spec.lam[keep_acs]/zcorr, acs_spec.flux[keep_acs]/norm_acs, color='blue', linewidth=3, alpha=0.8)
        ax.errorbar(lci/zcorr, fobs, yerr=efobs, marker='o', markersize=10, color='orange', ecolor='orange', linestyle='None', alpha=0.3)
        ax.text(0.95, 0.85, r'$z=%.2f$' %(zi), ha='right', transform = ax.transAxes, fontsize=20)
        ax.text(0.95, 0.1, r'$i_{814}=%5.2f$' %(acs.mag[i]), ha='right', transform = ax.transAxes, fontsize=16)
        ax.text(0.95, 0.2, r'$H_{140}=%5.2f$' %(phot.mag_f1392w[phot.id == wfc_file][0]), ha='right', transform = ax.transAxes, fontsize=16)
        
        #ax.set_xlabel(r'$\lambda_\mathrm{rest}\ [\AA]$')
        ax.set_xlabel(r'$\lambda_\mathrm{obs}\ [\AA]$')
        
        #from pylab import *
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter

        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2500)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)

        #for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)

        ax.set_xlim(3000,1.9e4)
        ymax_acs = np.max(acs_spec.flux[keep_acs]/norm_acs)
        ymax_wfc = np.max(wfc_spec.flux[keep_wfc]/norm_wfc)
        
        ymax = np.max([ymax_acs, ymax_wfc])
        
        ax.set_ylim(-0.1*ymax, 1.5*ymax)
        
        plt.savefig('%04.2f_' %(zi) + acs_file+'_example.pdf')
        plt.close()
        
        
        
import os
import pyfits
import numpy as np
import glob
import shutil

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
    from unicorn.catalogs import match_string_arrays
    import cosmocalc
        
    ######################################
    #### Redshifts and masses
    ######################################
    
    zout = catIO.Readfile('full_redshift.zout')
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
    lines.halpha_eqw[idx_lines] /= (1+zout.z_peak[0::3][found_lines])
    found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
    
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
    gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(1./np.abs(gfit.ba)) 
    gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(1./np.abs(gfit.ba)) 
    
    ##### Selection slices
    dr = mcat.rmatch[idx] < 1
    zrange = (zout.z_peak[0::3] > 0.2)
    
    keep = dr & zrange & (mcat.logm[idx] > 10.98) & (mcat.fcontam[idx] < 0.5)
    
    keep = dr & zrange & (mcat.mag_f140w[idx] < 24) & (mcat.fcontam[idx] < 0.5)
    
    #### refine
    zrange = (zout.z_peak[0::3] > 2) & (zout.z_peak[0::3] < 3)
    zrange = (zout.z_peak[0::3] > 1.5) & (zout.z_peak[0::3] < 2.)
    zrange = (zout.z_peak[0::3] > 0.4) & (zout.z_peak[0::3] < 1.0)
    zrange = (zout.z_peak[0::3] > 1.0) & (zout.z_peak[0::3] < 1.5)
    keep = dr & zrange & (zout.q_z[0::3] < 0.1) &  (mcat.fcontam[idx] < 0.2)
    
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
    plot_init()
    
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
    plot_init()
    
    plt.semilogy(mcat.logm[idx][~zsp & keep], -lines.halpha_eqw[idx_lines][~zsp & keep], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & keep], -lines.halpha_eqw[idx_lines][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    # oiii = keep & (zout.z_peak[0::3] > 1.2) & (-lines.oiii_eqw[idx_lines] > 25)
    # plt.semilogy(mcat.logm[idx][oiii], -lines.halpha_eqw[idx_lines][oiii], marker='s', linestyle='None', color='red', alpha=0.2, markersize=8)
    
    ## red sequence
    red_limit = 15
    red_sequence = -lines.halpha_eqw[idx_lines] < red_limit
    plt.plot([0,100],[red_limit,red_limit], linestyle='--', alpha=0.5, color='green', linewidth=4)
    
    ### For SDSS, read in catalogs and get selection from unicorn.paper1
    # plot_init()
    # 
    # hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], np.log10(-p1.sdss_line.H_ALPHA_EQW[sel]), bins=100, range=[[9,11.6], [-1, 2.7]])
    # plt.contourf(xedge[1:], 10**yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    # 
    # plt.semilogy(mcat.logm[idx][keep], -lines.halpha_eqw[idx_lines][keep], marker='o', linestyle='None', color='red', alpha=0.5, markersize=5)
    
    plt.text(9.5,1,r'has $z_\mathrm{spec}$', color='blue', fontsize=18)
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
    plot_init()
    
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
    plot_init()
    
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
    plot_init()
    
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

def plot_init():
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times']
    
    fig = plt.figure(figsize=(7,5), dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.10,
                        bottom=0.10,right=0.99,top=0.97)
    return fig
    
def make_wget_script(objects):
    """
    Make a script to get spectra and thumbnails from unicorn.
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FIRST_PAPER')
    
    os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.5/images')
    line = 'tar czvf thumbs.tar.gz'
    for object in objects:
        line += ' %s_thumb.fits.gz' %(object))
    
    os.system(line)
    os.system('mv thumbs.tar.gz /3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/')
    os.chdir('/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/OUTPUT/')
    line = 'tar czvf spec.tar.gz'
    for object in objects:
        line += ' %s.tempfilt' %(object))
    
    
    
    
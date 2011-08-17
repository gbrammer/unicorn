#!/usr/bin/env python
# encoding: utf-8
"""
paper1.py

Scripts for G. Brammer's first paper with the 3D-HST data.  General utilities for redshift fitting, etc. are also in 'unicorn.analysis'.

$URL: https://subversion.assembla.com/svn/threedhst_internal/trunk/reduce.py $
$Author: gbrammer $
$Date: 2011-05-22 02:01:43 -0400 (Sun, 22 May 2011) $

"""

__version__ = " $Rev: 5 $"

import glob
import os

import pyfits
import numpy as np

import matplotlib.pyplot as plt
USE_PLOT_GUI=False
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import threedhst
import unicorn
import cosmocalc

sdss_iclass = None
sdss_info = None
sdss_line = None
sdss_totsfr = None
sdss_totspecsfr = None
sdss_logm = None
sdss_sersic = None
sdss_vagc = None

idx_mpa_vagc = None

sdss_selection = False

noNewLine = '\x1b[1A\x1b[1M'

def read_sdss(read_full_lines=True):
    """
    
    Read MPA-JHU SDSS libraries.  
    
    If 'read_full_lines' is set, then also read the full line catalog with fluxes etc.
    This is a big file so can take more time/memory.  
    """
    import unicorn.paper1
    
    sdss_path = unicorn.GRISM_HOME + 'ANALYSIS/SDSS/'
    
    print 'SDSS ICLASS...'
    if unicorn.paper1.sdss_iclass is None:
        unicorn.paper1.sdss_iclass = pyfits.open(sdss_path + 'gal_iclass_dr7_v5_2.fits.gz')[0].data
    
    print 'SDSS Info, redshifts etc....'    
    if unicorn.paper1.sdss_info is None:
        unicorn.paper1.sdss_info = pyfits.open(sdss_path + 'gal_info_dr7_v5_2.fits.gz')[1].data
    
    if read_full_lines:
        print 'SDSS Full line data...'
        if unicorn.paper1.sdss_line is None:
            unicorn.paper1.sdss_line = pyfits.open(sdss_path + 'gal_line_dr7_v5_2.fits')[1].data
    
    print 'SDSS Total SFR...'
    if unicorn.paper1.sdss_totsfr is None:
        unicorn.paper1.sdss_totsfr = pyfits.open(sdss_path + 'gal_totsfr_dr7_v5_2.fits.gz')[1].data

    print 'SDSS Total Specific SFR...'    
    if unicorn.paper1.sdss_totspecsfr is None:
        unicorn.paper1.sdss_totspecsfr = pyfits.open(sdss_path + 'gal_totspecsfr_dr7_v5_2.fits.gz')[1].data
    
    print 'SDSS Total stellar mass...'
    if unicorn.paper1.sdss_logm is None:
        unicorn.paper1.sdss_logm = pyfits.open(sdss_path + 'totlgm_dr7_v5_2.fits.gz')[1].data
    
    print 'SDSS-VAGC Sersic fits...'
    if unicorn.paper1.sdss_sersic is None:
        unicorn.paper1.sdss_sersic = pyfits.open(sdss_path + 'sersic_catalog.fits')[1].data
    #
    print 'SDSS-VAGC Coords...'
    if unicorn.paper1.sdss_vagc is None:
        unicorn.paper1.sdss_vagc = pyfits.open(sdss_path + 'object_catalog.fits')[1].data
    
    print 'Match between catalogs'
    if unicorn.paper1.idx_mpa_vagc is None:
        unicorn.paper1.idx_mpa_vagc = np.loadtxt(sdss_path + 'mpa_vagc.match', dtype=np.int)
    
    #### Keep only matches from VAGC
    unicorn.paper1.sdss_sersic = unicorn.paper1.sdss_sersic[unicorn.paper1.idx_mpa_vagc]
    
    #### Get R50 in kpc
    zgrid = np.arange(100)/100.*4+1./100
    scale = zgrid*0.
    for i in range(100):
        cc = cosmocalc.cosmocalc(zgrid[i])
        scale[i] = cc['PS_kpc']
    
    unicorn.paper1.sdss_sersic.SERSIC_R50_i_KPC = unicorn.paper1.sdss_sersic.SERSIC_R50[:,3]*np.interp(unicorn.paper1.sdss_info.Z, zgrid, scale)
    
def match_sersic():
    """
    Need to match MPA-JHU catalogs to the VAGC sersic catalog by RA/Dec.
    """
    import unicorn.paper1 as p1
    N = len(p1.sdss_info)
    match_idx = np.zeros(N, dtype=np.int)
    idx = np.arange(len(p1.sdss_vagc), dtype=np.int)
    
    in_spectro = p1.sdss_vagc.SDSS_SPECTRO_TAG >= 0
    # idx = idx[in_spectro]
    # ra_vagc= p1.sdss_vagc.RA[in_spectro]
    # de_vagc= p1.sdss_vagc.DEC[in_spectro]
    
    cosdec = np.cos(p1.sdss_info.DEC/360*2*np.pi)
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SDSS')
    fp = open('mpajhu.radec','w')
    for i in range(N):
        fp.write(' %13.6f %13.6f\n' %(p1.sdss_info.RA[i], p1.sdss_info.DEC[i]))
    fp.close()

    fp = open('vagc.radec','w')
    for i in range(len(p1.sdss_vagc)):
        fp.write(' %13.6f %13.6f\n' %(p1.sdss_vagc.RA[i], p1.sdss_vagc.DEC[i]))
    fp.close()
    
    ### gcc match.c
    ### ./a.out
    ### gzip *.radec
    
    # i=0; j=59367
    # js = np.loadtxt('mpa_vagc.match')
    # j = js[i]
    # dr = 3600.*np.sqrt(((p1.sdss_info.RA[i]-p1.sdss_vagc.RA[j])*np.cos(p1.sdss_info.DEC[i]))**2+(p1.sdss_info.DEC[i]-p1.sdss_vagc.DEC[j])**2)
    # print dr
        
        
def sdss_selection(zmin=0.05, zmax=0.2, type='GALAXY', massmin=8, massmax=12):
    import unicorn.paper1 as p1
    #
    redshift = (p1.sdss_info.Z >= zmin) & (p1.sdss_info.Z <= zmax) & (p1.sdss_info.Z_WARNING == 0)
    #
    target_type = p1.sdss_info.TARGETTYPE == type
    #
    mass_selection = (p1.sdss_logm.AVG >= massmin) & (p1.sdss_logm.AVG <= massmax)
    #
    selection = redshift & target_type & mass_selection
    #
    return selection, len(p1.sdss_logm.AVG[selection])
    
def testing():
    """
    Make some simple plots like SSFR vs M
    """
    import unicorn.paper1 as p1
    
    sel, NOBJ = p1.sdss_selection(zmax=0.2)
    
    #### Mass vs. sSFR
    plt.plot(p1.sdss_logm.AVG[sel], p1.sdss_totspecsfr.AVG[sel], marker='.', color='red', alpha=0.01, linestyle='None')
    
    plt.xlim(9,12)
    plt.ylim(-13,-8.5)
    
    #### Star-forming galaxies, S/N limit on line fluxes
    SN_line_limit = 3
    SN_lines = (p1.sdss_line.NII_6584_FLUX/p1.sdss_line.NII_6584_FLUX_ERR > SN_line_limit) & (p1.sdss_line.H_ALPHA_FLUX/p1.sdss_line.H_ALPHA_FLUX_ERR > SN_line_limit) & (p1.sdss_line.OIII_5007_FLUX/p1.sdss_line.OIII_5007_FLUX_ERR > SN_line_limit) & (p1.sdss_line.H_BETA_FLUX/p1.sdss_line.H_BETA_FLUX_ERR > SN_line_limit)
    
    gal_lines = sel & SN_lines # & (p1.sdss_totspecsfr > -10.7)
    
    ### BPT diagram
    bptx = np.log10( p1.sdss_line.NII_6584_FLUX / p1.sdss_line.H_ALPHA_FLUX )
    bpty = np.log10( p1.sdss_line.OIII_5007_FLUX / p1.sdss_line.H_BETA_FLUX )
    
    # plt.plot(bptx[gal_lines], bpty[gal_lines], marker='.', color='black', alpha=0.01, linestyle='None')
    
    hist, xedge, yedge = np.histogram2d(bptx[gal_lines], bpty[gal_lines], bins=100, range=[[-1.5,0.8], [-1.2,1.5]])
    
    # plt.imshow(hist.transpose(), interpolation='nearest')
    # plt.contour(xedge[1:], yedge[1:], hist.transpose(), [64,256,512], colors='red', alpha=1.0, linethick=2)

    Vbins = [2, 4, 8, 16, 32, 64, 128, 256, 512, 4096]
    values =   1.-np.arange(len(Vbins))*1./len(Vbins)
    Vcolors = []
    for i in range(len(Vbins)):
        Vcolors.append('%f' %(values[i]))
    
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    
    ### kauffman separation between SF galaxies and AGN
    xsf = np.arange(100)/100.*1.5-1.5
    ysf = 0.61/(xsf-0.05)+1.3
    plt.plot(xsf, ysf, color='b', alpha=0.6, linewidth=2)
    
    #### star-forming galaxies
    ysf_int = np.interp(bptx, xsf, ysf, right=-10, left=-10)
    sfg = (bpty < ysf_int) & (bptx < -0.1) & gal_lines
    qui = ((bpty > ysf_int) | (bptx > -0.1)) & sel
    agn = (bpty > ysf_int) & gal_lines
    
    plt.plot(bptx[sfg], bpty[sfg], color='blue', alpha=0.01, marker='.', linestyle='None')
    plt.plot(bptx[agn], bpty[agn], color='red', alpha=0.01, marker='.', linestyle='None')
    
    plt.xlim(-1.5, 0.8)
    plt.ylim(-1.2, 1.5)
    plt.xlabel('[NII]6584 / H'+r'$\alpha$')
    plt.ylabel('[OIII]5007 / H'+r'$\beta$')
    
    #### Mass vs. sSFR, for SFGs and AGN
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], p1.sdss_totspecsfr.AVG[sel], bins=100, range=[[9,12], [-13, -8.5]])
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    
    plt.plot(p1.sdss_logm.AVG[sfg], p1.sdss_totspecsfr.AVG[sfg], marker='.', color='blue', alpha=0.01, linestyle='None')
    plt.plot(p1.sdss_logm.AVG[agn], p1.sdss_totspecsfr.AVG[agn], marker='.', color='red', alpha=0.01, linestyle='None')
    
    plt.xlim(9,12)
    plt.ylim(-13,-8.5)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel('sSFR')
    
    #### Mass vs. H-a eqw
    #plt.plot(p1.sdss_logm.AVG[sel], p1.sdss_line.H_ALPHA_EQW[sel], marker='.', color='black', alpha=0.01, linestyle='None')
    
    ### log
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], np.log10(-p1.sdss_line.H_ALPHA_EQW[sel]), bins=100, range=[[9,11.6], [-1, 2.7]])
    plt.contourf(xedge[1:], 10**yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    plt.semilogy()
    
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], -p1.sdss_line.H_ALPHA_EQW[sel], bins=100, range=[[9,12], [-10, 80]])
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    plt.gray()
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    
    plt.plot(p1.sdss_logm.AVG[sfg], p1.sdss_line.H_ALPHA_EQW[sfg], marker='.', color='blue', alpha=0.01, linestyle='None')
    plt.plot(p1.sdss_logm.AVG[agn], p1.sdss_line.H_ALPHA_EQW[agn], marker='.', color='red', alpha=0.01, linestyle='None')
    
    plt.xlim(9,12)
    plt.ylim(-80, 10)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{EQW\ H}\alpha$')
    
    ### test colors
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], p1.sdss_line.H_ALPHA_EQW[sel], bins=100, range=[[9,12], [-80, 10]])
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=gen_rgb_colortable(Vbins=Vbins, rgb=(1,1,1), reverse=True), alpha=1.0, linethick=2)
    
    NBIN = 25
    Vbins = [2, 4, 8, 16, 32, 64, 128, 256, 512, 4096]
    Vbins = list(np.array(Vbins)*100/NBIN)
    
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[agn], p1.sdss_line.H_ALPHA_EQW[agn], bins=NBIN, range=[[9,12], [-80, 10]])
    #plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=gen_rgb_colortable(Vbins=Vbins, rgb=(1,0,0), reverse=True), alpha=0.6, linethick=2)
    plt.contour(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors='red', alpha=0.6, linethick=2)

    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sfg], p1.sdss_line.H_ALPHA_EQW[sfg], bins=NBIN, range=[[9,12], [-80, 10]])
    #plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=gen_rgb_colortable(Vbins=Vbins, rgb=(0,0,1), reverse=True), alpha=0.6, linethick=2)
    plt.contour(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors='blue', alpha=0.6, linethick=2)
    
    
    #### Eq-w vs sSFR
    plt.plot(p1.sdss_line.H_ALPHA_EQW[sel], p1.sdss_totspecsfr.AVG[sel], marker='.', color='black', alpha=0.01, linestyle='None')
    plt.xlim(-80, 10)
    plt.ylim(-13, -8.5)
    plt.ylabel(r'$\mathrm{EQW\ H}\alpha$')
    plt.xlabel(r'$\mathrm{sSFR}$')
    
    ##### Check OIII 4959/5007 line ratios for SF galaxies
    sfg_weak_oiii = sfg & (p1.sdss_line.OIII_4959_FLUX/p1.sdss_line.OIII_4959_FLUX_ERR > SN_line_limit)
    plt.plot(np.log10(p1.sdss_line.OIII_5007_FLUX[sfg_weak_oiii]), (p1.sdss_line.OIII_4959_FLUX/p1.sdss_line.OIII_5007_FLUX)[sfg_weak_oiii], marker='.', linestyle='None', alpha=0.01, color='black')
    
    avg_ratio = threedhst.utils.biweight((p1.sdss_line.OIII_4959_FLUX/p1.sdss_line.OIII_5007_FLUX)[sfg_weak_oiii & (np.log10(p1.sdss_line.OIII_5007_FLUX) > 1.5)], mean=True)
    plt.plot([0,10],[1./2.98,1./2.98], color='red', linewidth=2)
    
    plt.ylim(-0.1,1)
    plt.xlim(0,4)
    
    ##### Case B H-a/H-b
    plt.plot(np.log10(p1.sdss_line.H_ALPHA_FLUX[sfg]), (p1.sdss_line.H_BETA_FLUX/p1.sdss_line.H_ALPHA_FLUX)[sfg], marker='.', linestyle='None', alpha=0.01, color='black')
    plt.plot([0,10],[1./2.86,1./2.86], color='red', linewidth=2)
    
    plt.ylim(-0.1,1)
    plt.xlim(0,4)
    
    ##### Mass vs sersic n
    plt.plot(p1.sdss_logm.AVG[sfg], p1.sdss_sersic.SERSIC_N[sfg,3], marker='.', color='blue', alpha=0.01, linestyle='None')
    plt.plot(p1.sdss_logm.AVG[qui], p1.sdss_sersic.SERSIC_N[qui,3], marker='.', color='red', alpha=0.01, linestyle='None')
    
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], p1.sdss_sersic.SERSIC_N[sel, 3], bins=100, range=[[9,12], [0, 6]])
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{Sersic}\ n$')
    
    plt.xlim(9,12)
    plt.ylim(0,6)

    ##### Mass vs size
    plt.plot(p1.sdss_logm.AVG[sfg], (p1.sdss_sersic.SERSIC_R50_i_KPC[sfg]), marker='.', color='blue', alpha=0.01, linestyle='None')
    plt.plot(p1.sdss_logm.AVG[qui], (p1.sdss_sersic.SERSIC_R50_i_KPC[qui]), marker='.', color='red', alpha=0.01, linestyle='None')
    
    hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], np.log10(p1.sdss_sersic.SERSIC_R50[sel, 3]), bins=100, range=[[9,12], [-0.3, 1]])
    plt.contourf(xedge[1:], yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$R_{50}$')
    
    plt.xlim(9,12)
    plt.ylim(0,6)
    
def gen_rgb_colortable(Vbins = range(10), rgb = (1,0,0), reverse=False):
    values = np.arange(len(Vbins))*1./len(Vbins)
    if reverse:
        values = 1.-values
    Vcolors = []
    for i in range(len(Vbins)):
        Vcolors.append((rgb[0]*values[i], rgb[1]*values[i], rgb[2]*values[i]))
    #
    return Vcolors
    

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
import unicorn.brown_dwarf

import re

noNewLine = '\x1b[1A\x1b[1M'

root = None

left = 0.1
bottom = 0.13
dy2d = 0.67
aspect = 0.65
temp_color = (8/255.,47/255.,101/255.)
lrange = np.array([1.05e4,1.68e4])
spec_linewidth=2
pad_linewidth=2

import unicorn
unicorn.catalogs.read_catalogs()
from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp

USE_TEX = True

def run_all():
    import unicorn.object_examples
    unicorn.object_examples.agn_group()
    unicorn.object_examples.z4_quasar()
    unicorn.object_examples.big_dead_galaxy()
    unicorn.object_examples.high_signal_to_noise_galaxy()
    unicorn.object_examples.l_dwarf()
    unicorn.object_examples.t_dwarf()
    
def agn_group():
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
    
    ########   AGN/Quasars
    ### Binary quasar: GOODS-N-42-G141_00388/384
    ### z=2.2 quasar: COSMOS-1-G141_00206	
    ### very broad H-a line, z=1.22: PRIMO-1101-G141_00993
    ### Even more interesting merger/quasar, z=1.778: GOODS-N-36-G141_00991	
    ### Another mess: COSMOS-3-G141_01156, z=1.34
    ### Multiple components, z=1.27 GOODS-N-33-G141_01028/1073/1069/1055
    ### z=4.6, MgII: COSMOS-28-G141_00896
    
    ###################################################
    ####
    ####      AGN group
    ####
    ###################################################
    
    # GOODS-N-36-G141_00991 / 1005
    # for object in ['GOODS-N-36-G141_00991','GOODS-N-36-G141_01005']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/%s* DATA/' %(object))

    thumb = pyfits.open('DATA/GOODS-N-36-G141_00991_thumb.fits.gz')
    twod = pyfits.open('DATA/GOODS-N-36-G141_00991_2d.fits.gz')
    spec2d = twod[1].data
    y0, y1 = 24, 79
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

    spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.2, vmax=0.025, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect*plot_aspect/pix_aspect, 0.99-bottom-dy2d))

    ax.imshow(0-thumb[0].data[y0:y1, y0+5:y1+5], vmin=-0.8, vmax=0.1, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='GOODS-N-36-G141_00991', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-18*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.7)
    
    ## Secondary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='GOODS-N-36-G141_01005', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-18*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    #ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='orange', linewidth=1, alpha=0.7)
    
    #### 
    zspec = 1.773
    mag = phot.mag_f1392w[phot.id == 'GOODS-N-36-G141_00991'][0]
    ax.text(0.05,0.8,r'$a)\ z=%.3f,\ m_{140}=%.1f$' %(zspec, mag), transform=ax.transAxes, fontsize=11)
    
    lines = [4102, 4341, 4862, 4980*1.08]
    y0 = [0.7, 0.7, 1, 1.5]
    labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','[OIII]4959+5007']
    for i in range(len(lines)):
        ax.text(lines[i]*(1+zspec), 3*y0[i], labels[i], horizontalalignment='center')
        
    
    ax.set_ylim(-0.1,ymax*1.1)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    
    print 'Savefig'
    print os.getcwd()
    fig.savefig('agn_group.pdf')
    
    ######## Brown dwarf
    ###  AEGIS-3-G141_00195  T-type
    ###  GOODS-N-24-G141_01148  L-type
    
    ####### Massive galaxies
    ###  z=2.0, huge, old: COSMOS-26-G141_00725
    ### z=1.9, zspec, beautiful fit UDF: PRIMO-1101-G141_01022
    
    
def z4_quasar():

    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
    
    ########   AGN/Quasars
    ### Binary quasar: GOODS-N-42-G141_00388/384
    ### z=2.2 quasar: COSMOS-1-G141_00206	
    ### very broad H-a line, z=1.22: PRIMO-1101-G141_00993
    ### Even more interesting merger/quasar, z=1.778: GOODS-N-36-G141_00991	
    ### Another mess: COSMOS-3-G141_01156, z=1.34
    ### Multiple components, z=1.27 GOODS-N-33-G141_01028/1073/1069/1055
    ### z=4.6, MgII: COSMOS-28-G141_00896
    
    ###################################################
    ####
    ####      z=4.6 quasar
    ####
    ###################################################
    
    # for object in ['COSMOS-28-G141_00896']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/OUTPUT/%s* DATA/' %(object))

    thumb = pyfits.open('DATA/COSMOS-28-G141_00896_thumb.fits.gz')
    twod = pyfits.open('DATA/COSMOS-28-G141_00896_2d.fits.gz')
    spec2d = twod[1].data
    y0, y1 = 10,32
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

    spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.2, vmax=0.025, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))

    ax.imshow(0-thumb[0].data[y0-2:y1-2, y0-2:y1-2], vmin=-2.4, vmax=0.3, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='COSMOS-28-G141_00896', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-18*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert*1.15
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.7)
        
    #### 
    zspec = 4.656
    mag = phot.mag_f1392w[phot.id == 'COSMOS-28-G141_00896'][0]
    ax.text(0.05,0.8,r'$b)\ z=%.3f,\ m_{140}=%.1f$' %(zspec, mag), transform=ax.transAxes, fontsize=11)

    lines = [2799, 2326, 2439.]
    y0 = [1.5, 0.7, 0.7, 0.7]
    labels = ['Mg II', 'C II', 'Ne IV']
    
    for i in range(len(lines)):
        ax.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
        
    ax.set_ylim(-0.1,ymax*1.1)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    
    ytick = ax.set_yticks([0,1,2])
    
    print 'Savefig'
    print os.getcwd()
    fig.savefig('z4_quasar.pdf')
    
def big_dead_galaxy():
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
        
    # for object in ['COSMOS-26-G141_00725']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/OUTPUT/%s* DATA/' %(object))
    
    thumb = pyfits.open('DATA/COSMOS-26-G141_00725_thumb.fits.gz')
    twod = pyfits.open('DATA/COSMOS-26-G141_00725_2d.fits.gz')
    spec2d = twod[1].data-twod[4].data
    y0, y1 = 24, 60
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)
    
    spec2d_sub = spec2d[y0+1:y1+1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1*1.2, vmax=0.0125*1.2, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))
    
    ax.imshow(0-thumb[0].data[y0:y1, y0:y1], vmin=-2.4, vmax=0.3, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='COSMOS-26-G141_00725', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-18*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    temp_sed *= 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**-18*(lambdaz/5500.)**2
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.7)
        
    #### 
    zspec = 2.0832
    mag = phot.mag_f1392w[phot.id == 'COSMOS-26-G141_00725'][0]
    ax.text(0.05,0.8,r'$c)\ z=%.3f,\ m_{140}=%.1f$' %(zspec, mag), transform=ax.transAxes, fontsize=11)
    
    # lines = [4102, 4341, 4862, 4980]
    # y0 = [0.7, 0.7, 1, 1.5]
    # labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','O III 4959+5007']
    # for i in range(len(lines)):
    #     ax.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
    
    ax.set_ylim(-0.1,ymax*1.2)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    ytick = ax.set_yticks([0,1,2,3])
    
    #### Inset full sed
    ax = fig.add_axes((left+0.55, bottom+0.1, 0.99-left-0.6, dy2d*0.4))
    
    ax.plot(lci[is_spec], fobs[is_spec], alpha=0.9, color='black', linewidth=2)
    ax.plot(lambdaz,temp_sed, color='red', linewidth=1, alpha=0.3)
    ax.plot(lci[~is_spec], fobs[~is_spec], marker='o', linestyle='None', alpha=0.3, color='black')
    
    ax.semilogx()
    ax.set_xlim(3000,9.e4)
    ax.set_ylim(-0.1*ymax,ymax*1.2)
    
    ax.set_yticklabels([])
    ax.set_xticklabels([r'$10^4$',r'$5\times10^4$'])
    xtick = ax.set_xticks([1.e4,5.e4]); ytick = ax.set_yticks([0,1,2,3])
    
    #print os.getcwd()
    fig.savefig('big_dead_galaxy.pdf')


def high_signal_to_noise_galaxy():
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
    
    ########   AGN/Quasars
    ### Binary quasar: GOODS-N-42-G141_00388/384
    ### z=2.2 quasar: COSMOS-1-G141_00206	
    ### very broad H-a line, z=1.22: PRIMO-1101-G141_00993
    ### Even more interesting merger/quasar, z=1.778: GOODS-N-36-G141_00991	
    ### Another mess: COSMOS-3-G141_01156, z=1.34
    ### Multiple components, z=1.27 GOODS-N-33-G141_01028/1073/1069/1055
    ### z=4.6, MgII: COSMOS-28-G141_00896
    
    ###################################################
    ####
    ####      z=4.6 quasar
    ####
    ###################################################
    # 
    # for object in ['PRIMO-1101-G141_01022']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/%s* DATA/' %(object))
    
    thumb = pyfits.open('DATA/PRIMO-1101-G141_01022_thumb.fits.gz')
    twod = pyfits.open('DATA/PRIMO-1101-G141_01022_2d.fits.gz')
    spec2d = twod[1].data-twod[4].data
    y0, y1 = 31, 56
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

    spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1*0.8, vmax=0.0125*0.8, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))

    ax.imshow(0-thumb[0].data[y0:y1, y0:y1], vmin=-1.4, vmax=0.15, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='PRIMO-1101-G141_01022', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-19*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    temp_sed *= 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**-19*(lambdaz/5500.)**2
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.7)
        
    #### 
    zspec = 1.905
    mag = phot.mag_f1392w[phot.id == 'PRIMO-1101-G141_01022'][0]
    #mass = mcat.lmass
    ax.text(0.05,0.8,r'$d)\ z=%.3f,\ m_{140}=%.1f$' %(zspec, mag), transform=ax.transAxes, fontsize=11)
    
    # lines = [4102, 4341, 4862, 4980]
    # y0 = [0.7, 0.7, 1, 1.5]
    # labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','O III 4959+5007']
    # for i in range(len(lines)):
    #     ax.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
        
    ax.set_ylim(-0.1*ymax,ymax*1.3)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-19}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    ytick = ax.set_yticks([0,1,2,3,4,5])
    
    #### Inset full sed
    ax = fig.add_axes((left+0.55, bottom+0.1, 0.99-left-0.6, dy2d*0.4))
    
    ax.plot(lci[is_spec], fobs[is_spec], alpha=0.9, color='black', linewidth=2)
    ax.plot(lambdaz,temp_sed, color='red', linewidth=1, alpha=0.3)
    ax.plot(lci[~is_spec], fobs[~is_spec], marker='o', linestyle='None', alpha=0.3, color='black')

    ax.semilogx()
    ax.set_xlim(3000,9.e4)
    ax.set_ylim(-0.1*ymax,ymax*1.1)
    
    ax.set_yticklabels([])
    ax.set_xticklabels([r'$10^4$',r'$5\times10^4$'])
    xtick = ax.set_xticks([1.e4,5.e4]); ytick = ax.set_yticks([0,1,2,3,4,5])
    
    
    print 'Savefig'
    #print os.getcwd()
    fig.savefig('high_signal_to_noise_galaxy.pdf')

#
def l_dwarf():

    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
    
    ########   AGN/Quasars
    ### Binary quasar: GOODS-N-42-G141_00388/384
    ### z=2.2 quasar: COSMOS-1-G141_00206	
    ### very broad H-a line, z=1.22: PRIMO-1101-G141_00993
    ### Even more interesting merger/quasar, z=1.778: GOODS-N-36-G141_00991	
    ### Another mess: COSMOS-3-G141_01156, z=1.34
    ### Multiple components, z=1.27 GOODS-N-33-G141_01028/1073/1069/1055
    ### z=4.6, MgII: COSMOS-28-G141_00896
    
    ###################################################
    ####
    ####      L dwarf
    ####
    ###################################################
    # 
    # for object in ['GOODS-N-24-G141_01148']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/%s* DATA/' %(object))

    thumb = pyfits.open('DATA/GOODS-N-24-G141_01148_thumb.fits.gz')
    twod = pyfits.open('DATA/GOODS-N-24-G141_01148_2d.fits.gz')
    spec2d = twod[1].data-twod[4].data
    y0, y1 = 10, 30
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

    spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1*0.8, vmax=0.0125*0.8, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))

    ax.imshow(0-thumb[0].data[y0:y1, y0:y1], vmin=-1.4, vmax=0.15, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='GOODS-N-24-G141_01148', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-19*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    temp_sed *= 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**-19*(lambdaz/5500.)**2
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    
    bd = unicorn.brown_dwarf.BD_fit()
    type = ['L4']
    ii = 0
    colors = ['green','blue','red','orange']
    for temp in bd.templates:
        if temp.type[0:2] in type:
            if temp.type == 'L3+/-1':
                continue
            
            print temp.type
            yint = np.interp(lci[is_spec], temp.wave, temp.flux)
            norm = np.sum(fobs[is_spec]*yint)/np.sum(yint**2)
            ax.plot(temp.wave, temp.flux*norm, color='white', linewidth=2, alpha=0.4)
            ax.plot(temp.wave, temp.flux*norm, color=colors[ii % 4], linewidth=2, alpha=0.7)
            ax.text(0.9-ii*0.08, 0.83, temp.type, color=colors[ii % 4], transform=ax.transAxes)
            ii = ii + 1
            
    #ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    #ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.5)
        
    #### 
    zspec = 1.905
    mag = phot.mag_f1392w[phot.id == 'GOODS-N-24-G141_01148'][0]
    ax.text(0.05,0.8,r'$f)\ m_{140}=%.1f$' %(mag), transform=ax.transAxes, fontsize=11)
    
    # lines = [4102, 4341, 4862, 4980]
    # y0 = [0.7, 0.7, 1, 1.5]
    # labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','O III 4959+5007']
    # for i in range(len(lines)):
    #     ax.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
        
    ax.set_ylim(-0.1*ymax,ymax*1.3)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-19}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    ytick = ax.set_yticks([0,5,10,15])
    
    
    print 'Savefig'
    #print os.getcwd()
    fig.savefig('l_dwarf.pdf')

def t_dwarf():

    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
    
    ########   AGN/Quasars
    ### Binary quasar: GOODS-N-42-G141_00388/384
    ### z=2.2 quasar: COSMOS-1-G141_00206	
    ### very broad H-a line, z=1.22: PRIMO-1101-G141_00993
    ### Even more interesting merger/quasar, z=1.778: GOODS-N-36-G141_00991	
    ### Another mess: COSMOS-3-G141_01156, z=1.34
    ### Multiple components, z=1.27 GOODS-N-33-G141_01028/1073/1069/1055
    ### z=4.6, MgII: COSMOS-28-G141_00896
    
    ###################################################
    ####
    ####      L dwarf
    ####
    ###################################################
    # 
    # for object in ['AEGIS-3-G141_00195']:
    #     os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
    #     os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/%s* DATA/' %(object))

    thumb = pyfits.open('DATA/AEGIS-3-G141_00195_thumb.fits.gz')
    twod = pyfits.open('DATA/AEGIS-3-G141_00195_2d.fits.gz')
    spec2d = twod[1].data-twod[4].data
    y0, y1 = 10, 30
    
    if USE_TEX:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    ax = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
    ax.plot([0,1])
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75
    
    plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
    pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

    spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1*0.8, vmax=0.0125*0.8, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,y1-y0])
    
    #### Thumb
    ax = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))

    ax.imshow(0-thumb[0].data[y0-1:y1-1, y0-1:y1-1], vmin=-1.4, vmax=0.15, interpolation='nearest', zorder=2, aspect='auto')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xtick = ax.set_xticks([0,y1-y0]); ytick = ax.set_yticks([0,y1-y0])
    
    #### Spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='AEGIS-3-G141_00195', OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-18*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    temp_sed *= 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**-18*(lambdaz/5500.)**2
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    ax.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    
    bd = unicorn.brown_dwarf.BD_fit()
    type = ['T6','T5']
    ii = 0
    colors = ['green','blue','red','orange']
    for temp in bd.templates:
        if temp.type[0:2] in type:
            if temp.type == 'L3+/-1':
                continue
            
            print temp.type
            yint = np.interp(lci[is_spec], temp.wave, temp.flux)
            norm = np.sum(fobs[is_spec]*yint)/np.sum(yint**2)
            ax.plot(temp.wave, temp.flux*norm, color='white', linewidth=2, alpha=0.4)
            ax.plot(temp.wave, temp.flux*norm, color=colors[ii % 4], linewidth=2, alpha=0.7)
            ax.text(0.9-ii*0.05, 0.83, temp.type, color=colors[ii % 4], transform=ax.transAxes, fontsize=11)
            ii = ii + 1
            
    #ax.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
    #ax.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.5)
        
    #### 
    zspec = 1.905
    mag = phot.mag_f1392w[phot.id == 'AEGIS-3-G141_00195'][0]
    ax.text(0.05,0.8,r'$e)\ m_{140}=%.1f$' %(mag), transform=ax.transAxes, fontsize=11)
    
    # lines = [4102, 4341, 4862, 4980]
    # y0 = [0.7, 0.7, 1, 1.5]
    # labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','O III 4959+5007']
    # for i in range(len(lines)):
    #     ax.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
        
    ax.set_ylim(-0.1*ymax,ymax*1.3)
    ax.set_xlim(lrange[0], lrange[1])
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ [10^{-18}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    ytick = ax.set_yticks([0,1,2])
    
    
    print 'Savefig'
    #print os.getcwd()
    fig.savefig('t_dwarf.pdf')

    ######## Brown dwarf
    ###  AEGIS-3-G141_00195  T-type
    ###  GOODS-N-24-G141_01148  L-type
    
    ####### Massive galaxies
    ###  z=2.0, huge, old: COSMOS-26-G141_00725
    ### z=1.9, zspec, beautiful fit UDF: PRIMO-1101-G141_01022

def get_tdwarf_mag():
    """
    Get the broad / medium-band H magnitudes of the T dwarf
    to compare to m140
    """
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    object = 'AEGIS-3-G141_00195'
    ra = phot.x_world[phot.id == object][0]
    dec = phot.y_world[phot.id == object][0]
    m140 = phot.mag_f1392w[phot.id == object][0]
    
    nmbs_cat, nmbs_zout, nmbs_fout = unicorn.analysis.read_catalogs(root=object)
    dr = np.sqrt((nmbs_cat.ra-ra)**2*np.cos(dec/360*2*np.pi)**2+(nmbs_cat.dec-dec)**2)*3600.
    h1mag = 25-2.5*np.log10((nmbs_cat.H1*nmbs_cat.Ktot/nmbs_cat.K)[dr == dr.min()][0])
    h2mag = 25-2.5*np.log10((nmbs_cat.H2*nmbs_cat.Ktot/nmbs_cat.K)[dr == dr.min()][0])
    hmag = 25-2.5*np.log10(((nmbs_cat.H1+nmbs_cat.H2)/2.*nmbs_cat.Ktot/nmbs_cat.K)[dr == dr.min()][0])
    jmag = 25-2.5*np.log10(((nmbs_cat.J2+nmbs_cat.J3)/2.*nmbs_cat.Ktot/nmbs_cat.K)[dr == dr.min()][0])
    jmag = 25-2.5*np.log10(((nmbs_cat.J3)/1.*nmbs_cat.Ktot/nmbs_cat.K)[dr == dr.min()][0])
    
    wirds = catIO.Readfile('/Users/gbrammer/research/drg/PHOTZ/EAZY/WIRDS/WIRDS_D3-95_Ks_ugrizJHKs_141927+524056_T0002.cat.candels')
    dr = np.sqrt((wirds.ra-ra)**2*np.cos(dec/360.*2*np.pi)**2+(wirds.dec-dec)**2)*3600.
    jwirds = wirds.jtot[dr == dr.min()][0]
    hwirds = wirds.htot[dr == dr.min()][0]
    
    print '       J     H     J-H   H1    H2'
    print 'NMBS  %5.2f %5.2f %5.2f %5.2f %5.2f' %(jmag, hmag, jmag-hmag, h1mag, h2mag)
    print 'WIRDS %5.2f %5.2f %5.2f' %(jwirds, hwirds, jwirds-hwirds)
    
    #### Vrba et al. (2004)
    #absH = np.array([14.52,14.78,15.07])
    #d =
     
def misc_objects():
    
    unicorn.object_examples.general_plot('UDF-Full-G141_00624', flam_norm=-19, vscale=0.1, vthumb=(-0.08*0.3,0.01*0.3), SED_voffset=0.42, SED_hoffset=0.05, remove_contamination=False)

    
def general_plot(object='AEGIS-9-G141_00154', show_SED=True, sync=False, y0=None, y1=None, SED_voffset=0.1, SED_hoffset=0, plot_min=None, plot_max=None, yticks=None, fit_path='REDSHIFT_FITS_v1.6', dy_thumb=0, dx_thumb=0, remove_contamination=True, vscale=1, vthumb=(-1,0.1), fit_version=0, show_2D = True, show_Thumb=True, show_Fit=True, flam_norm=-18):
    
    import unicorn.object_examples
    dy2d = unicorn.object_examples.dy2d
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/OBJECT_EXAMPLES')
        
    ### F_lambda
    ## obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-18
        
    if sync:
        os.system('rsync -avz $UNICORN:/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images/%s* DATA/' %(object))
        os.system('rsync -avz $UNICORN:/3DHST/Spectra/Work/ANALYSIS/%s/OUTPUT/%s* DATA/' %(fit_path, object))
    
    zout_file = catIO.Readfile('DATA/%s.zout' %(object))
    
    thumb = pyfits.open('DATA/%s_thumb.fits.gz' %(object))
    twod = pyfits.open('DATA/%s_2d.fits.gz' %(object))
    spec2d = twod[1].data
    if remove_contamination:
        spec2d -= twod[4].data
        
    #y0, y1 = 24, 60
    if y0 is None:
        y0 = 0
    if y1 is None:
        y1 = spec2d.shape[0]
        print 'NY: %d' %(spec2d.shape[0])
        
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=aspect, left=0.12)
    
    #### Twod
    if show_2D:
        ax2D = fig.add_axes((left, bottom+dy2d, 0.99-left, 0.99-bottom-dy2d))
        ax2D.plot([0,1])

        head = twod[1].header
        lam_idx = np.arange(head['NAXIS1'])
        lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
        lam_mima = np.cast[int](np.round(np.interp(lrange, lam, lam_idx)))
        tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx) - np.interp(lrange[0], lam, lam_idx)-0.75

        plot_aspect = (bottom+dy2d)/(0.99-bottom-dy2d)/aspect
        pix_aspect = (lam_mima[1]-lam_mima[0])*1./(y1-y0)

        spec2d_sub = spec2d[y0:y1,lam_mima[0]:lam_mima[1]]
        ax2D.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1*1.2*vscale, vmax=0.0125*1.2*vscale, interpolation='nearest')
        ax2D.set_yticklabels([]); ax2D.set_xticklabels([])
        xtick = ax2D.set_xticks(tick_int); ytick = ax2D.set_yticks([0,y1-y0])
    
        #### Thumb
        if show_Thumb:
            axThumb = fig.add_axes((left, bottom+dy2d, (0.99-bottom-dy2d)*aspect, 0.99-bottom-dy2d))

            if dx_thumb is None:
                dx_thumb = dy_thumb

            axThumb.imshow(0-thumb[0].data[y0+dy_thumb:y1+dy_thumb, y0+dx_thumb:y1+dx_thumb], vmin=vthumb[0], vmax=vthumb[1], interpolation='nearest', zorder=2, aspect='auto')
            axThumb.set_yticklabels([])
            axThumb.set_xticklabels([])
            xtick = axThumb.set_xticks([0,y1-y0]); ytick = axThumb.set_yticks([0,y1-y0])
        else:
            axThumb=None
    else:
        ax2D = None
        axThumb=None
        dy2d = 0.99-bottom
        
    #### Spectrum
    axSpec = fig.add_axes((left, bottom, 0.99-left, dy2d))
    
    ## Primary
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(fit_version, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY='DATA', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**flam_norm*(lci/5500.)**2
    fobs, efobs, obs_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert
    
    temp_sed *= 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**flam_norm*(lambdaz/5500.)**2
    
    ymax = max(fobs[is_spec & (fobs > 0)])
    axSpec.plot(lci[is_spec],fobs[is_spec], color='black', linewidth=spec_linewidth)
    if show_Fit:
        axSpec.plot(lci[is_spec],obs_sed[is_spec], color='white', alpha=0.8, linewidth=pad_linewidth)
        axSpec.plot(lci[is_spec],obs_sed[is_spec], color='red', linewidth=1, alpha=0.7)
        
    #### 
    zspec = 2.0832
    #zspec = zout.z_peak[0::3][zout.id[0::3] == object]
    zspec = zout_file.z_peak[fit_version]
    
    try:
        mag = phot.mag_f1392w[phot.id == object][0]
    except:
        mag = -1
        
    axSpec.text(0.05,0.9,object, transform=axSpec.transAxes, fontsize=9)
    if mag > 0:
        axSpec.text(0.05,0.8,r'$ \ z=%.3f,\ m_{140}=%.1f$' %(zspec, mag), transform=axSpec.transAxes, fontsize=11)
    else:
        axSpec.text(0.05,0.8,r'$z=%.3f$' %(zspec), transform=axSpec.transAxes, fontsize=11)
        
    # lines = [4102, 4341, 4862, 4980]
    # y0 = [0.7, 0.7, 1, 1.5]
    # labels = [r'H$\delta$',r'H$\gamma$', r'H$\beta$','O III 4959+5007']
    # for i in range(len(lines)):
    #     axSpec.text(lines[i]*(1+zspec), 0.5*y0[i], labels[i], horizontalalignment='center')
    
    if plot_min is None:
        plot_min = -0.1*ymax
    if plot_max is None:
        plot_max = 1.2*ymax
        
    axSpec.set_ylim(plot_min,plot_max)
    axSpec.set_xlim(lrange[0], lrange[1])
    axSpec.set_xlabel(r'$\lambda\ [\AA]$')
    axSpec.set_ylabel(r'$f_\lambda\ [10^{%0d}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$' %(flam_norm))
    if yticks is not None:
        ytick = axSpec.set_yticks(yticks)
    
    #### Inset full sed
    if show_SED:
        axInset = fig.add_axes((left+0.55+SED_hoffset, bottom+SED_voffset, 0.99-left-0.6, dy2d*0.4))
        axInset.plot(lci[is_spec], fobs[is_spec], alpha=0.9, color='black', linewidth=1)
        axInset.plot(lambdaz,temp_sed, color='red', linewidth=1, alpha=0.3)
        axInset.plot(lci[~is_spec], fobs[~is_spec], marker='o', linestyle='None', alpha=0.5, color='white')
        axInset.plot(lci[~is_spec], fobs[~is_spec], marker='o', linestyle='None', alpha=0.3, color='black')

        axInset.semilogx()
        axInset.set_xlim(3000,9.e4)
        axInset.set_ylim(-0.1*ymax,ymax*1.2)

        axInset.set_xticklabels([r'$10^4$',r'$5\times10^4$'])
        xtick = axInset.set_xticks([1.e4,5.e4])
        if yticks is not None:
            axInset.set_yticklabels([])
            ytick = axInset.set_yticks(yticks)
    else:
        axInset = None
    #print os.getcwd()
    fig.savefig(object+'_display.pdf')
    return fig, ax2D, axThumb, axSpec, axInset
        
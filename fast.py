
import os
import pyfits
import numpy as np
import glob
import shutil
import time

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

from pyraf import iraf
from iraf import iraf

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

import cosmocalc

def go_all():
    import glob
    
    force=False
    files=glob.glob('FAST_OUTPUT/*fout')
    for file in files: 
        object = os.path.basename(file).split('_threedhst')[0]
        if (not os.path.exists(object+'_fast.png')) | force:
            check_fast(object=object)
        
def check_fast(object='AEGIS-1-G141_00497', wmin=2000, wmax=2.4e4, logx=False, image_type='png'):

    if object.startswith('GOODS-S') | object.startswith('WFC3') | object.startswith('GEORGE') | object.startswith('PRIMO'):
        abzp=23.86
    else:
        abzp = 25

    obs_sed = catIO.Readfile('ASCII/%s_obs_sed.dat' %(object))
    temp_sed = catIO.Readfile('ASCII/%s_temp_sed.dat' %(object))

    lc = obs_sed.lc
    dlam_spec = lc[-1]-lc[-2]
    is_spec = np.append(np.abs(1-np.abs(lc[1:]-lc[0:-1])/dlam_spec) < 0.05,True)

    obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-19
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=2./3, left=0.12)
    
    ax = fig.add_subplot(111)
    
    ax.plot(lc[~is_spec], obs_sed.fnu[~is_spec]*obs_convert[~is_spec], marker='o', color='orange', linestyle='None', markersize=15, alpha=0.7)
    ax.plot(lc[is_spec], obs_sed.fnu[is_spec]*obs_convert[is_spec], color='blue', alpha=0.5)
    ax.plot(lc[~is_spec], obs_sed.obs_sed[~is_spec]*obs_convert[~is_spec], color='red', alpha=0.7, marker='o', linestyle='None', markersize=8)

    temp_convert = 10**(-0.4*(abzp+48.6))*3.e18/temp_sed.lam**2/10.**-19
    flam_temp = temp_sed.fnu_temp*temp_convert

    ax.plot(temp_sed.lam, flam_temp, color='red', alpha=0.3)
    #fast_norm = 1./np.interp(2.2e4, wfast, tfast)*np.interp(2.2e4, temp_sed.lam, flam_temp)

    wfast, tfast = np.loadtxt('FAST_OUTPUT/BEST_FITS/%s_threedhst_1.fit' %(object), skiprows=1, unpack=True)
    ax.plot(wfast, tfast, color='green', alpha=0.5)

    wfast, tfast = np.loadtxt('FAST_OUTPUT/BEST_FITS/%s_threedhst_2.fit' %(object), skiprows=1, unpack=True)
    ax.plot(wfast, tfast, color='purple', alpha=0.5)
    
    if logx:
        ax.semilogx()
    ax.set_xlim(wmin,wmax)

    ymax = np.max(obs_sed.fnu*obs_convert)
    ax.set_ylim(-0.1*ymax,1.3*ymax)

    fout = catIO.Readfile('FAST_OUTPUT/%s_threedhst.fout' %(object))
    
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'$f_\lambda\ (10^{-19}$)')
    xtext = 2e4
    xal = 'right'
    
    ax.text(xtext,0.4*ymax,r'$z_\mathrm{gris}=%.4f$' %(fout.z[0]), horizontalalignment=xal)

    ax.text(3000,1.1*ymax, object)
    
    ax.text(xtext,0.3*ymax,r'log M: $%.1f^{\ %.1f}_{\ %.1f}$   $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.lmass[0], fout.u68_lmass[0], fout.l68_lmass[0], fout.lmass[1], fout.u68_lmass[1], fout.l68_lmass[1]), horizontalalignment=xal)

    ax.text(xtext,0.2*ymax,r'$A_V$: $%.1f^{\ %.1f}_{\ %.1f}$   $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.av[0], fout.u68_av[0], fout.l68_av[0], fout.av[1], fout.u68_av[1], fout.l68_av[1]), horizontalalignment=xal)
    
    ax.text(xtext,0.1*ymax,r'log $\tau$: $%.1f^{\ %.1f}_{\ %.1f}$    $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.ltau[0], fout.u68_ltau[0], fout.l68_ltau[0], fout.ltau[1], fout.u68_ltau[1], fout.l68_ltau[1]), horizontalalignment=xal)
    
    ax.text(xtext,0.0*ymax,r'log Age: $%.1f^{\ %.1f}_{\ %.1f}$    $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.lage[0], fout.u68_lage[0], fout.l68_lage[0], fout.lage[1], fout.u68_lage[1], fout.l68_lage[1]), horizontalalignment=xal)
    
    ### Save to PNG
    outfile = object+'_fast.'+image_type
    
    if USE_PLOT_GUI:
        fig.savefig(outfile,dpi=100,transparent=False)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)
    
    fig.close()

def make_masked_FAST_errfunc():
    """
    Make a template error function for FAST that is very large near emission lines to 
    apply an effective mask.
    """
    temperr = '/usr/local/share/FAST/FAST_v0.9b/Template_error/TEMPLATE_ERROR.fast.v0.2'

    w,t = np.loadtxt(temperr, unpack=True)
    
    lines = [4861, 5007, 6563.]
    for line in lines:
        mask = np.abs(w-line) < 250
        t[mask] = 1000
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FAST')
    fp = open(os.path.basename(temperr)+'_linemask','w')
    for i in range(len(w)):
        fp.write('%f %f\n' %(w[i], t[i]))
    
    fp.close()

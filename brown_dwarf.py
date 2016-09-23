import os
from astropy.io import fits as pyfits
import numpy as np
import glob
import shutil
import time

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.ticker as mticker

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

HAS_PHOTOMETRY = True
PHOTOMETRY_ID = None
BAD_SPECTRUM = False

SPC_FILENAME = None
SPC = None

def fit_all_brown_dwarfs():
    import glob
    import unicorn
    
    bd = unicorn.analysis.BD_fit()
    
    files = glob.glob('/Users/gbrammer/Sites_GLOBAL/P/GRISM/ascii/*-G141*.dat')
    for file in files:
        bd.fit(file)

def fit_candidates():
    import unicorn
    
    bd = unicorn.analysis.BD_fit()
    os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM/BROWN_DWARF')
    
    bd.fit('../ascii/AEGIS-3-G141_00196.dat')
    bd.fit('../ascii/GOODS-N-33-G141_00923.dat')
    bd.fit('../ascii/GOODS-N-11-G141_00668.dat')
    bd.fit('../ascii/AEGIS-11-G141_00311.dat', chi2_limit=4, trim_mtype=False)
    bd.fit('../ascii/GEORGE-G141_00825.dat', chi2_limit=100, trim_mtype=False, max_contam=0.5)
    bd.fit('../ascii/GEORGE-G141_00326.dat', chi2_limit=100, trim_mtype=False, max_contam=0.8)
    bd.fit('../ascii/COSMOS-3-G141_00799.dat', chi2_limit=100, trim_mtype=False, max_contam=0.1)
    
class BD_template():
    def __init__(self, txt):
        self.filename = txt
        self.read_template(txt)
    
    def read_template(self, txt):
        import numpy as np
        fp = open(txt)
        lines = fp.readlines()
        fp.close()
        
        wave = []
        flux = []
        err = []
        for line in lines:
            if 'type:' in line:
                self.type = line.split('type:')[1].strip()
            if not line.startswith('#'):
                spl = line.split()
                wave.append(spl[0])
                flux.append(spl[1])
                err.append(spl[2])
        
        self.wave = np.cast[float](wave)*1.e4
        self.flux = np.cast[float](flux)
        self.err = np.cast[float](err)
    
class BD_fit():    
    def __init__(self):
        self.template_path=unicorn.GRISM_HOME+'ANALYSIS/BROWN_DWARF/STANDARDS'
        self.read_bd_templates()
        
    def read_bd_templates(self):
        import glob
        temps = glob.glob(self.template_path+'/spex*txt')
        if len(temps) == 0:
            threedhst.showMessage('No BD templates found in %s' %(self.template_path), warn=True)
            
        list = []
        for temp in temps:
            list.append(BD_template(temp))
        
        self.templates = list
        self.NTEMP = len(self.templates)
        
    def fit(self, ascii_file='AEGIS-3-G141_00177.dat', chi2_limit=1.5, trim_mtype=True, max_contam=0.05, xrange=(1.1e4,1.65e4), img_type='png', flux_min=0):
        import threedhst.catIO as catIO
        import numpy as np
        
        ### handle v2.0 format spectra
        spec = catIO.Readfile(ascii_file)
        if 'trace' in spec.columns:
            spec.lam = spec.wave*1.
            spec.flux /= spec.sensitivity
            spec.error /= spec.sensitivity
            spec.contam /= spec.sensitivity
            ascii_file = ascii_file.replace('.1D','')
        else:
            spec.error /= 2.5
        
        spec.flux -= spec.contam
        
        self.spec = spec
        self.img_type = img_type
        
        chi2 = np.zeros(self.NTEMP)
        types = []
        anorm = chi2*0.

        use = (spec.lam > xrange[0]) & (spec.lam < xrange[1]) & (spec.contam/spec.flux < max_contam) & (np.isfinite(spec.flux)) & (spec.flux > flux_min)
        self.use = use
        
        if len(spec.lam[use]) < 50:
            return False
            
        for i in range(self.NTEMP):
            temp = self.templates[i]
            types.append(temp.type.strip())
            yint = np.interp(spec.lam[use], temp.wave, temp.flux)
            #
            anorm[i] = np.sum(yint*spec.flux[use]*spec.error[use]**2) / np.sum(yint**2*spec.error[use]**2)
            #
            chi2[i] = np.sum((anorm[i]*yint-spec.flux[use])**2/spec.error[use]**2)
        
        types = np.cast[str](types)
        DOF = len(yint)-1
        chi2 /= DOF
        
        min = np.where(chi2 == chi2.min())[0][0]
        
        if trim_mtype:
            trim_mtype = types[min].startswith('M')
            
        if (chi2.min() < chi2_limit) & (~trim_mtype):
            print unicorn.noNewLine + ascii_file+' * '+' %s %0.2f' %(types[min], chi2.min())
            self.spec = spec
            self.types = types
            self.chi2 = chi2
            self.anorm = anorm
            self.ascii_file = ascii_file
            self.make_plot()
        else:
            print unicorn.noNewLine + ascii_file + ' %s %0.2f' %(types[min], chi2.min())
            
    def make_plot(self):
        import matplotlib.pyplot as plt
        spec = self.spec
        types = self.types
        chi2 = self.chi2
        anorm = self.anorm
        use = self.use
        
        so = np.argsort(types)
                
        plt.rcParams['font.size'] = 10
        #fig = plt.figure(figsize=[6.5,3],dpi=100)
        fig = Figure(figsize=[6.5, 3], dpi=100)        
        fig.subplots_adjust(wspace=0.25,hspace=0.02,left=0.09,
                            bottom=0.13,right=0.98,top=0.98)
        
        #### Plot the spectrum
        ax = fig.add_subplot(121)
        ymax = spec.flux[use].max()
        ax.plot(spec.lam, spec.flux/ymax, color='black', linewidth=2, alpha=0.1)
        ax.plot(spec.lam[use], spec.flux[use]/ymax, color='black', linewidth=2)
        ax.plot(spec.lam, spec.contam/ymax, color='yellow', linewidth=1)
        
        soc = np.argsort(chi2)
        colors=['blue','green']
        for i in range(self.NTEMP):
            j = soc[i]
            ax.plot(self.templates[j].wave, self.templates[j].flux*anorm[j]/ymax, alpha=0.1)
            
        for i in range(2):
            j = soc[i]
            ax.plot(self.templates[j].wave, self.templates[j].flux*anorm[j]/ymax, color=colors[i], linewidth=1)
            ax.text(1.6e4,1.1-0.1*i,types[j], color=colors[i])
        
        ax.text(1.1e4, 1.07, os.path.splitext(os.path.basename(self.ascii_file))[0])
        
        ax.set_xticks(np.arange(1.1,1.8,0.2)*1.e4)
        ax.set_xticklabels(np.arange(1.1,1.8,0.2))
        
        ax.set_ylim(-0.1,1.2)
        ax.set_xlim(1.05e4,1.72e4)
        ax.set_ylabel(r'$f_\lambda$')
        ax.set_xlabel(r'$\lambda / \mu\mathrm{m}$')
        
        #### Plot chi2 as a function of type
        ax = fig.add_subplot(122)
        ax.plot(np.arange(self.NTEMP), chi2[so], color=(0.9, 0.3, 0.3))
        for i in range(self.NTEMP):
            j = so[i]
            tt = ax.text(i, chi2[j], types[j], fontsize=8)
        
        ax.set_xticklabels([])
        ax.set_xlim(-1,self.NTEMP+1)    
        ax.set_ylim(chi2.min()-1,chi2.min()+9)
        ax.set_xlabel('Type')
        ax.set_ylabel(r'$\chi^2_\nu$')
        
        outfile = os.path.basename(os.path.splitext(os.path.basename(self.ascii_file))[0]+'_BD.%s' %(self.img_type))

        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)

#

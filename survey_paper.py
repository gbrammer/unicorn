import os
import pyfits
import numpy as np
import glob
import shutil

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib

from pyraf import iraf
from iraf import iraf

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

import re

root = None

def throughput():
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER')
    
    xg141, yg141 = np.loadtxt('g141.dat', unpack=True)
    xf140, yf140 = np.loadtxt('f140w.dat', unpack=True)
    xf814, yf814 = np.loadtxt('f814w.dat', unpack=True)
    xg800l, yg800l = np.loadtxt('g800l.dat', unpack=True)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'

    fig = unicorn.catalogs.plot_init(square=True, xs=8, aspect=1./3, left=0.105, bottom=0.08, top=0.01, right=0.01)
    
    ax = fig.add_subplot(111)
    
    ax.plot(xg141, yg141, color='black', linewidth=2, alpha=0.5)
    ax.fill(xg141, yg141, color='red', linewidth=2, alpha=0.1)
    ax.plot(xf140, yf140, color='black', linewidth=2, alpha=0.7)
    
    ax.plot(xg800l, yg800l, color='black', linewidth=2, alpha=0.5)
    ax.fill(xg800l, yg800l, color='blue', linewidth=2, alpha=0.1)
    ax.plot(xf814, yf814, color='black', linewidth=2, alpha=0.7)
    
    em_lines = [3727, 4861, 4959, 5007, 6563.]
    offset = np.array([0,-0.05,0,0,0])
    offset = 0.25+np.array([-0.2,-0.05,0,-0.01,-0.1])
    xoffset = np.array([0,-120,200,150,0])
    line_scale = np.array([1,1,1./2.98,1,1])
    
    em_names = ['[OII]',r'H$\beta$','','[OIII]',r'H$\alpha$']
    
    dlam = 30
    zi = 1
    
    show_spectra = True
    colors=['blue','green','red']    
    if show_spectra:
      for zi in [1,2,3]:
        sedx, sedy = np.loadtxt('templates/EAZY_v1.0_lines/eazy_v1.0_sed4_nolines.dat', unpack=True)
        sedy *= 1.*sedy.max()
        dl = dlam/(1+zi)
        #dl = dlam
        for i,em_line in enumerate(em_lines):
            em_gauss = 1./np.sqrt(2*np.pi*dl**2)*np.exp(-1*(sedx-em_line)**2/2/dl**2)
            sedy += em_gauss/em_gauss.max()*0.6*line_scale[i]

        ax.plot(sedx*(1+zi), sedy*0.4+0.5, color=colors[zi-1], alpha=0.7, linewidth=2)
        ax.text(5500.,1.18-zi*0.13,r'$z=%d$' %(zi), color=colors[zi-1], fontsize=11)
        
      for i in range(len(em_lines)):
        ax.text(em_lines[i]*(1+1)+xoffset[i], 1+offset[i], em_names[i], horizontalalignment='center', fontsize=10)
    
    show_continuous = False
    if show_continuous:
        em_lines = [3727, 5007, 6563.]
        zgrid = np.arange(1000)/1000.*4
        for line in em_lines:
            ax.plot(line*(1+zgrid), zgrid/4.*0.8+0.5, linewidth=2, alpha=0.5, color='black')
        for zi in [0,1,2,3]:
            ax.plot([0.1,2.e4],np.array([zi,zi])/4.*0.8+0.5, linestyle='--', color='black', alpha=0.2)
            
    
    ax.text(5800, 0.08,'G800L',rotation=33., color='black', alpha=0.7)
    ax.text(5800, 0.08,'G800L',rotation=33., color='blue', alpha=0.4)

    ax.text(7100, 0.03,'F814W',rotation=80., color='black', alpha=0.9)

    ax.text(1.115e4, 0.17,'G141',rotation=15., color='black', alpha=0.7)
    ax.text(1.115e4, 0.17,'G141',rotation=15., color='red', alpha=0.4)

    ax.text(1.21e4, 0.03,'F140W',rotation=88., color='black', alpha=0.9)

    
    ax.set_xlim(4500, 1.79e4)
    ax.set_ylim(0,1.4)
    ax.set_xlabel(r'$\lambda$ [\AA]')
    ax.set_ylabel('throughput')
    
    #ax.set_yticklabels([]); 
    ytick = ax.set_yticks([0,0.25,0.5,0.75,1.0])
    
    fig.savefig('throughput.pdf')
    
    plt.rcParams['text.usetex'] = False
#
def throughput_v2():
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER')
    
    xg141, yg141 = np.loadtxt('g141.dat', unpack=True)
    xf140, yf140 = np.loadtxt('f140w.dat', unpack=True)
    xf814, yf814 = np.loadtxt('f814w.dat', unpack=True)
    xg800l, yg800l = np.loadtxt('g800l.dat', unpack=True)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'Serif'
    plt.rcParams['font.serif'] = 'Times'

    fig = unicorn.catalogs.plot_init(square=True, xs=8, aspect=1./3, left=0.095, bottom=0.08, top=0.065, right=0.01)
    
    #ax = fig.add_subplot(111)
    #ax = fig.add_axes(((x0+(dx+x0)*0), y0+0.5, dx, 0.5-top_panel-y0))
    ysplit = 0.65
    ax = fig.add_axes((0.06, 0.135, 0.935, (ysplit-0.135)))
    
    ax.plot(xg141, yg141, color='black', linewidth=2, alpha=0.5)
    ax.fill(xg141, yg141, color='red', linewidth=2, alpha=0.1)
    ax.plot(xf140, yf140, color='black', linewidth=2, alpha=0.7)
    
    ax.plot(xg800l, yg800l, color='black', linewidth=2, alpha=0.5)
    ax.fill(xg800l, yg800l, color='blue', linewidth=2, alpha=0.1)
    ax.plot(xf814, yf814, color='black', linewidth=2, alpha=0.7)
        
    em_names = ['[OII]',r'H$\beta$','','[OIII]',r'H$\alpha$']
    
    dlam = 30
    zi = 1
    
    yy = 0.1
    ax.text(5800, 0.12+yy,'G800L',rotation=48., color='black', alpha=0.7)
    ax.text(5800, 0.12+yy,'G800L',rotation=48., color='blue', alpha=0.4)

    ax.text(7100, 0.14+yy,'F814W',rotation=80., color='black', alpha=0.9)

    ax.text(1.115e4, 0.17+yy,'G141',rotation=15., color='black', alpha=0.7)
    ax.text(1.115e4, 0.17+yy,'G141',rotation=15., color='red', alpha=0.4)

    ax.text(1.21e4, 0.14+yy,'F140W',rotation=88., color='black', alpha=0.9)
    
    #ax.set_yticklabels([]); 
    
    #ax2 = ax.twiny()
    ax2 = fig.add_axes((0.06, ysplit+0.02, 0.935, (0.86-ysplit)))
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.set_ticks_position('top')
    ### H-alpha
    xbox = np.array([0,1,1,0,0])
    dy= 0.333333
    y0= 1.0
    ybox = np.array([0,0,1,1,0])*dy
    
    
    width_acs = np.array([6000.,9000.])
    width_wfc3 = np.array([1.1e4,1.65e4])
    
    line_names = [r'H$\alpha$',r'H$\beta$ / [OIII]','[OII]']
    
    for i, l0 in enumerate([6563.,4934.,3727]):
        zline = width_acs/l0-1
        ax2.fill(xbox*(zline[1]-zline[0])+zline[0],y0-ybox-i*dy, color='blue', alpha=0.1)
        zline = width_wfc3/l0-1
        ax2.fill(xbox*(zline[1]-zline[0])+zline[0],y0-ybox-i*dy, color='red', alpha=0.1)
        ax2.plot([0,4],np.array([0,0])+y0-(i+1)*dy, color='black')
        ax2.text(3.7,y0-(i+0.5)*dy, line_names[i], horizontalalignment='right', verticalalignment='center')
        
    ax.set_xlim(4500, 1.79e4)
    ax.set_ylim(0,0.65)
    ytick = ax.set_yticks([0,0.2,0.4,0.6])
    ax.set_xlabel(r'$\lambda$ [\AA]')
    ax.set_ylabel('Throughput')
    minorLocator   = MultipleLocator(1000)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator   = MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(minorLocator)
    
    ax2.set_xlim(0,3.8)
    ytick = ax2.set_yticks([])
    ax2.set_ylim(0,1)
    minorLocator   = MultipleLocator(0.1)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.set_xlabel(r'$z_\mathrm{line}$')
    
    fig.savefig('throughput.eps')
    plt.rcParams['text.usetex'] = False
    
def orbit_structure():
    """
    Show the POSTARG offsets in WFC3 / ACS
    """
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER')
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    wfc3_color, acs_color = 'red','blue'
    wfc3_color, acs_color = 'blue','green'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=4.4, aspect=1, left=0.09, bottom=0.08)
    ax = fig.add_subplot(111)
    
    a11 = 0.1355
    b10 = 0.1211 # arcsec / pix, from instrument HB
    
    #dxs = np.array([0,-20,-13,7]) + np.int(np.round(xsh[0]))*0
    #dys = np.array([0,-7,-20,-13]) + np.int(np.round(ysh[0]))*0
    
    x3dhst = np.array([0, 1.355, 0.881, -0.474])/a11
    y3dhst = np.array([0, 0.424, 1.212, 0.788])/b10
    
    xgoodsn = np.array([0,0.6075, 0.270, -0.3375])/a11
    ygoodsn = np.array([0,0.1815, 0.6655, 0.484])/b10
    
    #### SN fields:
    test = """
    files=`ls ibfuw1*flt.fits.gz`
    for file in $files; do result=`dfitsgz $file |fitsort FILTER APERTURE POSTARG1 POSTARG2 | grep -v POST`; echo "${file} ${result}"; done
    """
    xmarshall = np.array([-0.34, -0.540, 0.0, 0.608, 0.273])/a11
    ymarshall = np.array([-0.34, -0.243, 0.0, 0.244, 0.302])/b10
    
    xgeorge = np.array([0.0, -0.608, 0.273, -0.340, 0.540])/a11
    ygeorge = np.array([0.0, 0.244, 0.302, -0.340, -0.243])/b10
    
    x41 = np.array([0.273, -0.608, 0.540, -0.340, -0.340, 0.540])/a11
    y41 = np.array([0.302, 0.244, -0.243, -0.301, -0.301, -0.243])/b10
    
    xprimo = np.array([0.0, 0.474, 0.290, 0.764, -0.290, 0.184])/a11
    yprimo = np.array([0.0, 0.424, -0.290, 0.134, 0.290, 0.714])/b10
    
    xers = np.array([-10.012, 9.988, 9.971, -9.958])/a11
    yers = np.array([5.058, 5.050, -5.045, -5.045])/b10
    
    xcooper = np.array([0, 0.6075, 0.270 ,-0.3375])/a11
    ycooper = np.array([0, 0.1815, 0.6655, 0.484])/b10
    
    xstanford = np.array([-0.169, 0.372, 0.169, -0.372])/a11
    ystanford = np.array([-0.242, 0.06064, 0.242, 0.06064])/b10
    xstanford += 0.2
    
    xoff = x3dhst
    yoff = y3dhst
    
    print np.round(xoff*10)/10.*2
    print np.round(yoff*10)/10.*2
    
    plt.plot(np.round(xoff*10)/10. % 1, np.round(yoff*10)/10. % 1)
    plt.xlim(-0.1,1.1); plt.ylim(-0.1,1.1)
    
        
    ax.plot(xoff, yoff, marker='o', markersize=10, color=wfc3_color, alpha=0.8, zorder=10)
    
    if 1 == 1:
        for i in range(4):
            ax.text(xoff[i], yoff[i]+0.5, 'F140W + G141', horizontalalignment='center', backgroundcolor='white', zorder=20)
            
    ax.set_xlabel(r'$x$ offset [pix]')
    ax.set_ylabel(r'$y$ offset [pix]')

    scale = 4
    x0 = -2.9
    y0 = -5.5
    
    ax.fill(np.array([0,1,1,0])*scale+x0, np.array([0,0,1,1])*scale+y0, color='white', zorder=10)
    ax.fill(np.array([0,1,1,0])*scale+x0, np.array([0,0,1,1])*scale+y0, color='black', alpha=0.1, zorder=11)
    ax.plot(np.array([0.5,0.5])*scale+x0, np.array([0,1])*scale+y0, color='black', alpha=0.2, zorder=12)
    ax.plot(np.array([0,1])*scale+x0, np.array([0.5,0.5])*scale+y0, color='black', alpha=0.2, zorder=12)
    
    ax.plot(np.abs(xoff-np.cast[int](xoff))*scale+x0, np.abs(yoff-np.cast[int](yoff))*scale+y0, marker='o', markersize=10, color=wfc3_color, alpha=0.8, zorder=13)
    ax.text(x0+scale/2., y0-1, 'WFC3 Primary', horizontalalignment='center')
    
    #plt.xlim(-5,11)
    #plt.ylim(-5,11)
    
    #### ACS:
    # XPOS = x*a11 + y*a10
    # YPOS = x*b11 + y*b10
    # 
    #       a10    a11      b10     b11
    # WFC: 0.0000   0.0494  0.0494  0.0040
    a10, a11, b10, b11 =  0.0000,0.0494,0.0494,0.0040
    #xoff = np.cumsum(xpostarg)
    #yoff = np.cumsum(ypostarg)

    acsang = 45.
    acsang = 92.16- -45.123
    xpos_acs, ypos_acs = threedhst.utils.xyrot(xpostarg, ypostarg, acsang)
    xpix_acs = xpos_acs / a11
    ypix_acs = (ypos_acs-xpix_acs*b11)/b10

    x0 = 5.5
    #y0 = -4.5
    
    ax.fill(np.array([0,1,1,0])*scale+x0, np.array([0,0,1,1])*scale+y0, color='white', zorder=10)
    ax.fill(np.array([0,1,1,0])*scale+x0, np.array([0,0,1,1])*scale+y0, color='black', alpha=0.1, zorder=11)
    ax.plot(np.array([0.5,0.5])*scale+x0, np.array([0,1])*scale+y0, color='black', alpha=0.2, zorder=12)
    ax.plot(np.array([0,1])*scale+x0, np.array([0.5,0.5])*scale+y0, color='black', alpha=0.2, zorder=12)

    ax.plot(np.abs(xpix_acs-np.cast[int](xpix_acs))*scale+x0, np.abs(ypix_acs-np.cast[int](ypix_acs))*scale+y0, marker='o', markersize=10, color=acs_color, alpha=0.8, zorder=13)
    #ax.plot(np.array([0,0.5,1,0.5])*scale+x0, np.array([0,0.5,0.5,1])*scale+y0, marker='o', marker='None', color=acs_color, linestyle='--', alpha=0.6, zorder=13)
    ax.plot(np.array([0,0.5,0.5,1])*scale+x0, np.array([0,1,0.5,0.5])*scale+y0, marker='o', color=acs_color, linestyle='--', alpha=0.6, zorder=13)
    
    ax.text(x0+scale/2., y0-1, 'ACS Parallel', horizontalalignment='center')
    
    #plt.grid(alpha=0.5, zorder=1, markevery=5)
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    ax.set_xlim(-5.9,12.5)
    #ax.set_ylim(-5.9,12.5)
    ax.set_ylim(-6.9,11.5)

    majorLocator   = MultipleLocator(5)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(1)

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.grid(alpha=0.5, zorder=1, which='major')
    ax.xaxis.grid(alpha=0.2, zorder=1, which='minor')

    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.grid(alpha=0.5, zorder=1, which='major')
    ax.yaxis.grid(alpha=0.2, zorder=1, which='minor')
        
    fig.savefig('dither_box.pdf')
    
    plt.rcParams['text.usetex'] = False
    
def exptimes():
    """
    Extract the range of exposure times from the formatted Phase-II files
    """
    os.system('grep F814W 12???.pro   |grep " S " > f814w.exptime')
    os.system('grep G800L 12???.pro   |grep " S " > g800l.exptime')
    os.system('grep F140W *.pro |grep MULTIAC |grep NSAMP > f140w.exptime')
    os.system('grep G141 *.pro |grep MULTIAC |grep NSAMP > g141.exptime')
    
    ### G141
    fp = open('g141.exptime')
    lines = fp.readlines()
    fp.close()
    nsamp = []
    object = []
    for line in lines:
        nsamp.append(line.split('NSAMP=')[1][0:2])
        object.append(line.split()[1])
        
    expsamp = np.zeros(17)
    expsamp[12] = 1002.94
    expsamp[13] = 1102.94
    expsamp[14] = 1202.94
    expsamp[15] = 1302.94
    expsamp[16] = 1402.94
    
    nsamp = np.cast[int](nsamp)
    nsamp+=1 ### this seems to be the case for all actual observations !!
    objects = object[::4]
    NOBJ = len(nsamp)/4
    texp = np.zeros(NOBJ)
    for i in range(NOBJ):
        texp[i] = np.sum(expsamp[nsamp[i*4:(i+1)*4]])
        print objects[i], texp[i]
        
    print 'G141:  %.1f - %.1f' %(texp.min(), texp.max())
    
def spectral_features():
    wmin, wmax = 1.1e4, 1.65e4
    
    print '%-8s  %.1f -- %.1f' %('Halpha', wmin/6563.-1, wmax/6563.-1)
    print '%-8s  %.1f -- %.1f' %('OIII', wmin/5007.-1, wmax/5007.-1)
    print '%-8s  %.1f -- %.1f' %('OII', wmin/3727.-1, wmax/3727.-1)
    print '%-8s  %.1f -- %.1f' %('4000', wmin/4000.-1, wmax/4000.-1)

    wmin, wmax = 0.55e4, 1.0e4
    print '\n\nACS\n\n'
    print '%-8s  %.1f -- %.1f' %('Halpha', wmin/6563.-1, wmax/6563.-1)
    print '%-8s  %.1f -- %.1f' %('OIII', wmin/5007.-1, wmax/5007.-1)
    print '%-8s  %.1f -- %.1f' %('OII', wmin/3727.-1, wmax/3727.-1)
    print '%-8s  %.1f -- %.1f' %('4000', wmin/4000.-1, wmax/4000.-1)

def aXe_model():
    import copy
    import scipy.ndimage as nd
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER')
    
    dir = pyfits.open('/research/HST/GRISM/3DHST/GOODS-S/PREP_FLT/UDF-F140W_drz.fits')
    gri = pyfits.open('/research/HST/GRISM/3DHST/GOODS-S/PREP_FLT/UDF-FC-G141_drz.fits')
    mod = pyfits.open('/research/HST/GRISM/3DHST/GOODS-S/PREP_FLT/UDF-FC-G141CONT_drz.fits')

    #### rotate all the images so that dispersion axis is along X
    angle = gri[1].header['PA_APER']#+180
    direct = nd.rotate(dir[1].data, angle, reshape=False)
    grism = nd.rotate(gri[1].data, angle, reshape=False)
    model = nd.rotate(mod[1].data, angle, reshape=False)
    
    xc, yc = 1877, 2175
    NX, NY = 1270, 365
    aspect = 3.*NY/NX
    
    xc, yc = 1731, 977
    NX, NY = 882, 467
    
    aspect = 1.*NY/NX
    
    plt.gray()
    
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'

    fig = unicorn.catalogs.plot_init(square=True, xs=8, aspect=aspect, left=0.12)
    
    fig.subplots_adjust(wspace=0.0,hspace=0.0,left=0.01,
                        bottom=0.005,right=0.99,top=0.995)
    
    fs1 = 12 ### Font size of label
    xlab, ylab = 0.04*NX/3., NY-0.02*NY-0.02*NY
    
    xd, yd = 65, 412 ### my direct images don't line up any more (deleted old ones?)
    xd, yd = 412, 65*2
    
    ax = fig.add_subplot(221)
    ax.imshow(0-direct[yc-NY/2+yd:yc+NY/2+yd, xc-NX/2+xd:xc+NX/2+xd], vmin=-0.2, vmax=0.02, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks([0,NX]); ytick = ax.set_yticks([0,NY])
    ax.text(xlab, ylab, 'a) Direct F140W', fontsize=fs1, verticalalignment='top')

    #ax.text(xlab, ylab, r'$%d\times\ $$%d^{\prime\prime}$' %(NX*0.06, NY*0.06), fontsize=18, backgroundcolor='white', verticalalignment='top')
    
    ax = fig.add_subplot(222)
    ax.imshow(0-grism[yc-NY/2:yc+NY/2, xc-NX/2:xc+NX/2], vmin=-0.04, vmax=0.004, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks([0,NX]); ytick = ax.set_yticks([0,NY])
    ax.text(xlab, ylab, 'b) Grism G141', fontsize=fs1, verticalalignment='top')

    ax = fig.add_subplot(224)
    diff = grism-model
        
    ### Flag em lines and 0th order
    dy0 = 20
    
    emx, emy = [223, 272, 487, 754, 520, 850, 565, 558, 51, 834, 345, 495], [122, 189, 83, 240, 148, 124, 336, 418, 338, 225, 197, 268]
    ax.plot(np.array(emx)+(xc-1731), np.array(emy)-dy0+(yc-977), marker='^', markersize=6, linestyle='None', color='blue', alpha=0.9)
    
    zx, zy = [301, 393, 648], [183, 321, 446]
    ax.plot(np.array(zx)+(xc-1731), np.array(zy)-dy0+(yc-977), marker='^', markersize=6, linestyle='None', markeredgecolor='black', markerfacecolor='None', alpha=0.9, markeredgewidth=1.2)
    
    fonts = matplotlib.font_manager.FontProperties()
    fonts.set_size(9)
    ax.legend(['Emission',r'0th order'], numpoints=1, prop=fonts, handletextpad=0.001, borderaxespad=0.001)
    # ax.text(0.04*NX/3., 0.02*NY, 'Em.', fontsize=fs1*0.8, backgroundcolor='white', verticalalignment='bottom', color='green')
    # ax.text(0.3*NX/3., 0.02*NY, r'0th', fontsize=fs1*0.8, backgroundcolor='white', verticalalignment='bottom', color='red')
    
    ax.imshow(0-diff[yc-NY/2:yc+NY/2, xc-NX/2:xc+NX/2], vmin=-0.02, vmax=0.002, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks([0,NX]); ytick = ax.set_yticks([0,NY])
    ax.text(xlab, ylab, r'd) Model-subtracted grism', fontsize=fs1, verticalalignment='top')
    
    ax = fig.add_subplot(223)
    ax.imshow(0-model[yc-NY/2:yc+NY/2, xc-NX/2:xc+NX/2], vmin=-0.04, vmax=0.004, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks([0,NX]); ytick = ax.set_yticks([0,NY])
    ax.text(xlab, ylab, r'c) aXe Model', fontsize=fs1, verticalalignment='top')
    
    fig.savefig('grism_model.pdf')
    plt.rcParams['text.usetex'] = False
    
def sync():
    pass
    """
paths="RAW HTML/scripts PREP_FLT"
dirs="AEGIS COSMOS ERS GOODS-N GOODS-S SN-GEORGE SN-MARSHALL SN-PRIMO UDS"
for dir in $dirs; do 
    mkdir ${dir}
    mkdir ${dir}/HTML
    for path in $paths; do 
        # du -sh ${dir}/${path}
        mkdir ${dir}/${path}
        rsync -avz --progress $UNICORN:/3DHST/Spectra/Work/${dir}/${path} ${dir}/${path} 
    done
done


    """
def all_pointings():
    pointings(ROOT='GOODS-SOUTH')
    pointings(ROOT='COSMOS')
    pointings(ROOT='AEGIS')
    pointings(ROOT='UDS')
    pointings(ROOT='GOODS-N')

def all_pointings_width():
    """
    This is the mosaic figure from the paper.  The individual fields are combined 
    manually with Adobe Illustrator.
    """
    from unicorn.survey_paper import pointings
    
    fs, left = 10, 0.22
    fs, left = 12, 0.25
    fs, left = 14, 0.28
    fs, left = 16, 0.32
    
    pointings(ROOT='GOODS-SOUTH', width=7, corner='ll', fontsize=fs, left=left)
    pointings(ROOT='COSMOS', width=6, corner='lr', fontsize=fs, left=left, right=0.03, bottom=0.115)
    pointings(ROOT='AEGIS', width=7, corner='ll', fontsize=fs, left=left, right=0.045)
    pointings(ROOT='UDS', width=9, corner='lr', fontsize=fs, left=left-0.02, right=0.04, top=0.02)
    pointings(ROOT='GOODS-N', width=6, corner='ur', fontsize=fs, left=left, bottom=0.115)
    
def pointings_with_status():
    """
    Highlight pointings that have been observed (status == 'Archived') 
    """
    from unicorn.survey_paper import pointings
    
    pointings(ROOT='GOODS-SOUTH', width=7, corner='ll', use_status=True)
    pointings(ROOT='GOODS-SOUTH', width=7, corner='ll', use_status=True, show_acs=False)
    pointings(ROOT='GOODS-SOUTH', width=7, corner='ll', use_status=True, show_wfc3=False)
    pointings(ROOT='COSMOS', width=6, corner='lr', use_status=True)
    pointings(ROOT='COSMOS', width=6, corner='lr', use_status=True, show_acs=False)
    pointings(ROOT='COSMOS', width=6, corner='lr', use_status=True, show_wfc3=False)
    pointings(ROOT='AEGIS', width=7, corner='ll', use_status=True)
    pointings(ROOT='AEGIS', width=7, corner='ll', use_status=True, show_acs=False)
    pointings(ROOT='AEGIS', width=7, corner='ll', use_status=True, show_wfc3=False)
    pointings(ROOT='UDS', width=9, corner='lr', use_status=True)
    pointings(ROOT='UDS', width=9, corner='lr', use_status=True, show_acs=False)
    pointings(ROOT='UDS', width=9, corner='lr', use_status=True, show_wfc3=False)

    pointings(ROOT='UDS', width=9, corner='lr', show_sn_fields=True, use_status=True)
    pointings(ROOT='GOODS-SOUTH', width=7, corner='ll', show_sn_fields=True, use_status=True)
           
def pointings(ROOT='GOODS-SOUTH', width=None, corner='lr', use_status=False, show_acs=True, show_wfc3=True, show_sn_fields=False, fontsize=10, left=22, right=0.02, top=0.01, bottom=0.11):
    """ 
    Make a figure showing the 3D-HST pointing poisitions, read from region files
    """
    import unicorn.survey_paper as sup
    
    plt.rcParams['lines.linewidth'] = 0.3
    wfc3_color = 'blue'
    acs_color = 'green'
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/SURVEY_PAPER')
    yticklab = None
    
    dx_ref = (53.314005633802822-52.886197183098595)*np.cos(-27.983151830808076/360*2*np.pi)
    CANDELS_PATH = '/research/HST/GRISM/3DHST/REGIONS/CANDELS/'
    candels_files = []
    candels_alpha = 0.3
    candels_color = '0.1'
    
    if use_status:
        pointing_list, pointing_status = np.loadtxt('/research/HST/GRISM/3DHST/REGIONS/pointing_status.dat', dtype=np.str, unpack=True)
        pointing_list, pointing_status = np.array(pointing_list), np.array(pointing_status)
    else:
        pointing_list, pointing_status = np.array([]), np.array([])
        
    #### GOODS-S
    if ROOT=='GOODS-SOUTH':
        x0, x1 = 53.314005633802822,  52.886197183098595
        y0, y1 = -27.983151830808076, -27.654474431818176
        xticklab = [r'$3^\mathrm{h}33^\mathrm{m}00^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}30^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(3,33,00, hours=True), degrees(3,32,30, hours=True), degrees(3,32,00, hours=True)]
        yticklab = [r'$-27^\circ40^\prime00^{\prime\prime}$', r'$45^\prime00^{\prime\prime}$', r'$-27^\circ50^\prime00^{\prime\prime}$', r'$55^\prime00^{\prime\prime}$']
        ytickv = [degrees(-27, 40, 00, hours=False), degrees(-27, 45, 00, hours=False), degrees(-27, 50, 00, hours=False), degrees(-27, 55, 00, hours=False)]
        candels_files = glob.glob(CANDELS_PATH+'/GOODS-S*reg')
        candels_files.extend(glob.glob(CANDELS_PATH+'/GOODS-W*reg'))
    
    #### COSMOS
    if ROOT=='COSMOS':
        x1, x0 = 149.99120563380279, 150.23823661971829
        y0, y1 = 2.1678109478476815, 2.5996973302980129
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$00^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(10,00,30, hours=True), degrees(10,00,00, hours=True)]
        yticklab = [r'$+2^\circ15^\prime00^{\prime\prime}$', r'$25^\prime00^{\prime\prime}$', r'$35^\prime00^{\prime\prime}$']
        ytickv = [degrees(02, 15, 00, hours=False), degrees(02, 25, 00, hours=False), degrees(02, 35, 00, hours=False)]
        candels_files = glob.glob(CANDELS_PATH+'/COSMOS*reg')
    
    #### AEGIS
    if ROOT=='AEGIS':
        x1, x0 = 214.49707154104345, 215.12704734584406
        y0, y1 = 52.680946433013482, 53.01597137966467
        xticklab = [r'$18^\mathrm{m}00^\mathrm{s}$', r'$14^\mathrm{h}19^\mathrm{m}00^\mathrm{s}$', r'$20^\mathrm{m}00^\mathrm{s}$']
        xticklab = [r'$18^\mathrm{m}$', r'$14^\mathrm{h}19^\mathrm{m}$', r'$20^\mathrm{m}$']
        xtickv = [degrees(14,18,00, hours=True), degrees(14,19,00, hours=True), degrees(14,20,00, hours=True)]
        yticklab = [r'$+52^\circ45^\prime00^{\prime\prime}$', r'$50^\prime00^{\prime\prime}$', r'$55^\prime00^{\prime\prime}$']
        ytickv = [degrees(52, 45, 00, hours=False), degrees(52, 50, 00, hours=False), degrees(52, 55, 00, hours=False)]
        candels_files = glob.glob(CANDELS_PATH+'/EGS*reg')
    
    #### UDS
    if ROOT=='UDS':
        x1, x0 = 34.116935194128146, 34.51871547581829
        y0, y1 = -5.2957542206957582, -5.0834327182123147+2./3600
        xticklab = [r'$18^\mathrm{m}00^\mathrm{s}$', r'$2^\mathrm{h}17^\mathrm{m}30^\mathrm{s}$', r'$17^\mathrm{m}00^\mathrm{s}$', r'$16^\mathrm{m}30^\mathrm{s}$']
        xtickv = [degrees(2,18,00, hours=True), degrees(2,17,30, hours=True), degrees(2,17,00, hours=True), degrees(2,16,30, hours=True)]
        yticklab = [r'$05^\prime00^{\prime\prime}$', r'$-5^\circ10^\prime00^{\prime\prime}$', r'$15^\prime00^{\prime\prime}$']
        ytickv = [degrees(-5, 05, 00, hours=False), degrees(-5, 10, 00, hours=False), degrees(-5, 15, 00, hours=False)]
        candels_files = glob.glob(CANDELS_PATH+'/UDS*reg')
    
    #### GOODS-N
    if ROOT=='GOODS-N':
        wfc3_color = 'orange'
        acs_color = None
        x1, x0 = 188.9139017749491, 189.44688055895648
        y0, y1 = 62.093791549511998, 62.384068625281309
        xticklab = [r'$12^\mathrm{h}37^\mathrm{m}30^\mathrm{s}$', r'$37^\mathrm{m}00^\mathrm{s}$', r'$36^\mathrm{m}30^\mathrm{s}$', r'$36^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(12,37,30, hours=True), degrees(12,37,00, hours=True), degrees(12,36,30, hours=True), degrees(12,36,00, hours=True)]
        yticklab = [r'$+62^\circ10^\prime00^{\prime\prime}$', r'$15^\prime00^{\prime\prime}$', r'$20^\prime00^{\prime\prime}$']
        ytickv = [degrees(62, 10, 00, hours=False), degrees(62, 15, 00, hours=False), degrees(62, 20, 00, hours=False)]
        candels_files = glob.glob(CANDELS_PATH+'/GOODSN-OR*reg')
        candels_files.extend(glob.glob(CANDELS_PATH+'/GOODSN-SK*reg'))
        
    #### Make square for given plot dimensions
    dx = np.abs(x1-x0)*np.cos(y0/360*2*np.pi)
    dy = (y1-y0)
    
    if width is None:
        width = 7*dx/dx_ref
    
    print '%s: plot width = %.2f\n' %(ROOT, width)
        
    fig = unicorn.catalogs.plot_init(square=True, xs=width, aspect=dy/dx, fontsize=fontsize, left=left, right=right, top=top, bottom=bottom)
    #fig = unicorn.catalogs.plot_init(square=True)
    
    ax = fig.add_subplot(111)
    
    polys = []
    for file in candels_files:
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        #
        polys.append(sup.polysplit(lines[1], get_shapely=True))
        #fi = ax.fill(wfcx, wfcy, alpha=candels_alpha, color=candels_color)
    
    sup.polys = polys
    
    un = polys[0]
    for pp in polys[1:]:
        un = un.union(pp)
    
    if un.geometryType() is 'MultiPolygon':
        for sub_poly in un.geoms:
            x,y = sub_poly.exterior.xy
            ax.plot(x,y, alpha=candels_alpha, color=candels_color, linewidth=1)
            ax.fill(x,y, alpha=0.1, color='0.7')
    else:        
        x,y = un.exterior.xy
        ax.plot(x,y, alpha=candels_alpha, color=candels_color, linewidth=1)
        ax.fill(x,y, alpha=0.1, color='0.7')
    
    files=glob.glob(unicorn.GRISM_HOME+'REGIONS/'+ROOT+'-[0-9]*reg')
    
    
    if ROOT == 'UDS':
        p18 = files.pop(9)
        print '\n\nPOP %s\n\n' %(p18)
        
    if show_sn_fields:
        files.extend(glob.glob(unicorn.GRISM_HOME+'REGIONS/SN*reg'))
        files.extend(glob.glob(unicorn.GRISM_HOME+'REGIONS/ERS*reg'))
        
    wfc3_polys = []
    acs_polys = []

    for file in files:
        #            
        base = os.path.basename(file.split('.reg')[0])
        #print base, base in pointing_list
        
        if base in pointing_list:
            status = pointing_status[pointing_list == base][0] == 'Archived'
        else:
            status = False
        
        if not use_status:
            status = True
            
        field = re.split('-[0-9]', file)[0]
        pointing = file.split(field+'-')[1].split('.reg')[0]
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        
        if base.startswith('SN') | base.startswith('ERS'):
            wfc3_color = 'purple'
            acs_color = None
            pointing = os.path.basename(field)
            status = True
            
        #
        wfc3_polys.append(sup.polysplit(lines[1], get_shapely=True)) 
        acs_polys.append(sup.polysplit(lines[2], get_shapely=True)) 
        acs_polys.append(sup.polysplit(lines[3], get_shapely=True)) 
        #
        wfcx, wfcy = sup.polysplit(lines[1])
        if show_wfc3:
            if status:
                fi = ax.fill(wfcx, wfcy, alpha=0.2, color=wfc3_color)
                fi = ax.plot(wfcx, wfcy, alpha=0.8, color=wfc3_color)
            else:
                fi = ax.fill(wfcx, wfcy, alpha=0.05, color=wfc3_color)
                fi = ax.plot(wfcx, wfcy, alpha=0.8, color=wfc3_color)
        #
        if acs_color is not None:
            acsx1, acsy1 = sup.polysplit(lines[2])
            acsx2, acsy2 = sup.polysplit(lines[3])
            #
            if show_acs:
                if show_wfc3:
                    afact = 3
                else:
                    afact = 3
                    
                if status:
                    fi = ax.fill(acsx1, acsy1, alpha=0.05*afact, color=acs_color)
                    fi = ax.fill(acsx2, acsy2, alpha=0.05*afact, color=acs_color)
                    #
                    pl = ax.plot(acsx1, acsy1, alpha=0.1*afact, color=acs_color)
                    pl = ax.plot(acsx2, acsy2, alpha=0.1*afact, color=acs_color)
                else:
                    pl = ax.plot(acsx1, acsy1, alpha=0.3*afact, color=acs_color)
                    pl = ax.plot(acsx2, acsy2, alpha=0.3*afact, color=acs_color)
                
        #
        xoff, yoff = 0.0, 0.0
        if ROOT=='GOODS-SOUTH':
            #print pointing
            if pointing == '36':
                xoff, yoff = 0.002,0.0075
            if pointing == '37':
                xoff, yoff = -0.005,-0.007
            if pointing == '38':
                xoff, yoff = 0.007,-0.007
        #    
        if show_wfc3:
            te = ax.text(np.mean(wfcx[:-1])+xoff, np.mean(wfcy[:-1])+yoff, pointing, va='center', ha='center', fontsize=13)
    
    #### Get field area from full WFC3 polygons
    un_wfc3 = wfc3_polys[0]
    for pp in wfc3_polys[1:]:
        un_wfc3 = un_wfc3.union(pp)
    #
    wfc3_union= []
    
    if un_wfc3.geometryType() is 'MultiPolygon':
        total_area = 0
        xavg, yavg, wht = 0, 0, 0
        for sub_poly in un_wfc3.geoms:
            area_i = sub_poly.area*np.cos(y0/360.*2*np.pi)
            total_area += area_i
            x,y = sub_poly.exterior.xy
            wfc3_union.append(sub_poly)
            xavg += np.mean(x)*area_i**2
            yavg += np.mean(y)*area_i**2
            wht += area_i**2
            #ax.plot(x,y, alpha=0.8, color='orange', linewidth=1)
        xavg, yavg = xavg/wht, yavg/wht
    else:        
        total_area = un_wfc3.area*np.cos(y0/360.*2*np.pi)
        x,y = un_wfc3.exterior.xy
        wfc3_union.append(un_wfc3)
        #ax.plot(x,y, alpha=0.8, color='orange', linewidth=1)
        xavg, yavg = np.mean(x), np.mean(y)
    
    #plt.plot([xavg,xavg],[yavg,yavg], marker='x', markersize=20)
    
    #### Get ACS overlap fraction
    if ROOT != 'GOODS-N':
        un_acs = acs_polys[0]
        for pp in acs_polys[1:]:
            un_acs = un_acs.union(pp)
        
        acs_union = []
        if un_acs.geometryType() is 'MultiPolygon':
            for sub_poly in un_acs.geoms:
                x,y = sub_poly.exterior.xy
                acs_union.append(sub_poly)
        else:        
            x,y = un_acs.exterior.xy
            acs_union.append(un_acs)
        wfc3_area = 0.
        acs_overlap_area = 0.
        for wun in wfc3_union:
            wfc3_area += wun.area*np.cos(y0/360.*2*np.pi)*3600.
            for aun in acs_union:
                overlap = wun.intersection(aun)
                acs_overlap_area += overlap.area*np.cos(y0/360.*2*np.pi)*3600.

        print '== Combined areas ==\nWFC3, ACS, frac: %.1f %.1f %.1f' %(wfc3_area, acs_overlap_area ,acs_overlap_area/wfc3_area*100)
        
        dummy = """
        wf = 147.3 + 122.2 + 121.9 + 114.0
        ac = 134.6 + 112.7 + 102.4 + 102.8
        print ac/wf*100.
        """
    #    
    if yticklab is not None:
        ax.set_xticklabels(xticklab)
        xtick = ax.set_xticks(xtickv)
        ax.set_yticklabels(yticklab)
        ytick = ax.set_yticks(ytickv)
    
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    
    fsi = '20'
    
    if ROOT == 'GOODS-SOUTH':
        field_label = 'GOODS-S'
    else:
        field_label = ROOT
        
    if corner=='lr':
        ax.text(0.95, 0.05,r'$\mathit{%s}$' %(field_label),
            horizontalalignment='right',
            verticalalignment='bottom',
            transform = ax.transAxes, fontsize=fsi)
    #
    if corner=='ll':
        ax.text(0.05, 0.05,r'$\mathit{%s}$' %(field_label),
            horizontalalignment='left',
            verticalalignment='bottom',
            transform = ax.transAxes, fontsize=fsi)
    #
    if corner=='ur':
        ax.text(0.95, 0.95,r'$\mathit{%s}$' %(field_label),
            horizontalalignment='right',
            verticalalignment='top',
            transform = ax.transAxes, fontsize=fsi)
    #
    if corner=='ul':
        ax.text(0.05, 0.95,r'$\mathit{%s}$' %(field_label),
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes, fontsize=fsi)
    
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    # print 'RA  - ', hexagesimal(x0), hexagesimal(x1)
    # print 'Dec - ', hexagesimal(y0, hours=False), hexagesimal(y1, hours=False)
    print 'RA  - ', hexagesimal(xavg)
    print 'Dec - ', hexagesimal(yavg, hours=False)
    print 'Area: %.1f\n' %(total_area*3600.)
    
    tag = ''
    if not show_acs:
        tag += '_noacs'
    if not show_wfc3:
        tag += '_nowfc3'
    if show_sn_fields:
        tag += '_sn'
        
    if use_status:
        plt.savefig('%s_pointings_status%s.pdf' %(ROOT, tag))
    else:
        plt.savefig('%s_pointings%s.pdf' %(ROOT, tag))
        
def get_UDF_center():
    file='/research/HST/GRISM/3DHST/REGIONS/GOODS-SOUTH-38.reg'
    field = re.split('-[0-9]', file)[0]
    pointing = file.split(field+'-')[1].split('.reg')[0]
    fp = open(file)
    lines = fp.readlines()
    fp.close()
    #
    px, py = sup.polysplit(lines[1], get_shapely=False)
    print 'UDF: ', hexagesimal(np.mean(px[:-1]), hours=True), hexagesimal(np.mean(py[:-1]), hours=False)
    
def degrees(deg, min, sec, hours=True):
    
    adeg = np.abs(deg)
    degrees = adeg + min/60. + sec/3600.
    if deg < 0:
        degrees *= -1
    
    if hours:
        degrees *= 360./24
        
    return degrees
    
def hexagesimal(degrees, hours=True, string=True):
    if hours:
        degrees *= 24/360.
    
    if degrees < 0:
        sign = -1
        si = '-'
    else:
        sign = 1
        si = ''
        
    degrees = np.abs(degrees)
    
    deg = np.int(degrees)
    min = np.int((degrees-deg)*60)
    sec = (degrees-deg-min/60.)*3600
    
    if string:
        return '%s%02d:%02d:%05.2f' %(si, deg, min, sec)
    else:
        return sign*deg, min, sec

def polysplit(region='polygon(150.099223,2.391097,150.086084,2.422515,150.050573,2.407277,150.064586,2.376241)', get_shapely=False):
    
    spl = region[region.find('(')+1:region.find(')')].split(',')
    px = spl[0::2]
    py = spl[1::2]
    px.append(px[0])
    py.append(py[0])
    px, py = np.cast[float](px), np.cast[float](py)
    
    if get_shapely:
        from shapely.geometry import Polygon
        list = []
        for i in range(len(px)):
            list.append((px[i], py[i]))
        
        poly = Polygon(tuple(list))
        return poly
    else:
        return px, py
#
def demo_background_subtract(root='COSMOS-13'):
    """
    Make a figure demonstrating the background subtraction of the grism images
    """
    import threedhst
    import threedhst.prep_flt_files
    import threedhst.grism_sky as bg
    
    import unicorn
    import unicorn.survey_paper as sup
    import pyfits
    import iraf
    
    path = unicorn.analysis.get_grism_path(root)
    os.chdir(path)
    if not os.path.exists('EXAMPLE'):
        os.system('mkdir EXAMPLE')
    
    os.chdir('EXAMPLE')
    files = glob.glob('../PREP_FLT/%s-G141*' %(root))
    files.append('../PREP_FLT/%s-F140W_tweak.fits' %(root))
    
    for file in files: 
        if 'drz.fits' not in file:
            os.system('cp %s .' %(file))
            print file
    
    threedhst.process_grism.fresh_flt_files(root+'-G141_asn.fits')
    asn = threedhst.utils.ASNFile(root+'-G141_asn.fits')
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    angle = flt[1].header['PA_APER']
    
    #### First run on uncorrected images
    threedhst.prep_flt_files.startMultidrizzle(root+'-G141_asn.fits', 
            use_shiftfile=True, skysub=False,
            final_scale=0.128254, pixfrac=0.8, driz_cr=True,
            updatewcs=True, median=True, clean=True, final_rot=angle)
        
    os.system('mv %s-G141_drz.fits %s-G141_drz_first.fits' %(root, root))
    
    sup.root = root
    sup.first_prof = []
    for exp in asn.exposures:
        xp, yp = threedhst.grism_sky.profile(exp+'_flt.fits', flatcorr=False, biweight=True)
        sup.first_prof.append(yp)
    
    # for i,exp in enumerate(asn.exposures):
    #     plt.plot(xp, prof[i])
        

    #### Now divide by the flat
    threedhst.process_grism.fresh_flt_files(root+'-G141_asn.fits')
    sup.flat_prof = []
    for exp in asn.exposures:
        bg.remove_grism_sky(flt=exp+'_flt.fits', list=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits'],  path_to_sky = '../CONF/', out_path='./', verbose=False, plot=False, flat_correct=True, sky_subtract=False, second_pass=False, overall=False)
        xp, yp = threedhst.grism_sky.profile(exp+'_flt.fits', flatcorr=False, biweight=True)
        sup.flat_prof.append(yp)
    
    threedhst.prep_flt_files.startMultidrizzle(root+'-G141_asn.fits', 
            use_shiftfile=True, skysub=False,
            final_scale=0.128254, pixfrac=0.8, driz_cr=True,
            updatewcs=True, median=True, clean=True, final_rot=angle)

    os.system('mv %s-G141_drz.fits %s-G141_drz_flat.fits' %(root, root))
    
    #### Divide by the sky
    threedhst.process_grism.fresh_flt_files(root+'-G141_asn.fits')
    sup.sky_prof = []
    for exp in asn.exposures:
        print exp
        bg.remove_grism_sky(flt=exp+'_flt.fits', list=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits'],  path_to_sky = '../CONF/', out_path='./', verbose=False, plot=False, flat_correct=True, sky_subtract=True, second_pass=False, overall=False)
        xp, yp = threedhst.grism_sky.profile(exp+'_flt.fits', flatcorr=False, biweight=True)
        sup.sky_prof.append(yp)
    
    threedhst.prep_flt_files.startMultidrizzle(root+'-G141_asn.fits', 
            use_shiftfile=True, skysub=False,
            final_scale=0.128254, pixfrac=0.8, driz_cr=True,
            updatewcs=True, median=True, clean=True, final_rot=angle)

    os.system('mv %s-G141_drz.fits %s-G141_drz_sky.fits' %(root, root))

    #### Last pass along columns
    threedhst.process_grism.fresh_flt_files(root+'-G141_asn.fits')
    sup.final_prof = []
    for exp in asn.exposures:
        print exp
        bg.remove_grism_sky(flt=exp+'_flt.fits', list=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits'],  path_to_sky = '../CONF/', out_path='./', verbose=False, plot=False, flat_correct=True, sky_subtract=True, second_pass=True, overall=True)
        xp, yp = threedhst.grism_sky.profile(exp+'_flt.fits', flatcorr=False, biweight=True)
        sup.final_prof.append(yp)
    
    threedhst.prep_flt_files.startMultidrizzle(root+'-G141_asn.fits', 
            use_shiftfile=True, skysub=False,
            final_scale=0.128254, pixfrac=0.8, driz_cr=True,
            updatewcs=True, median=True, clean=True, final_rot=angle)
    
    # ##### Make segmentation images
    # run = threedhst.prep_flt_files.MultidrizzleRun(root+'-G141')
    # for i,exp in enumerate(asn.exposures):
    #     run.blot_back(ii=i, copy_new=(i is 0))
    #     threedhst.prep_flt_files.make_segmap(run.flt[i])
    
    os.system('mv %s-G141_drz.fits %s-G141_drz_final.fits' %(root, root))
    
    make_background_demo(root=root)
    
def make_background_demo(root='AEGIS-11', range1=(0.90,1.08), range2=(-0.02, 0.02)):
    import unicorn.survey_paper as sup
    
    if sup.root != root:
        print "Need to run sup.demo_background_subtract(root='%s')." %(root)
        
    path = unicorn.analysis.get_grism_path(root)
    os.chdir(path+'EXAMPLE')
    
    first = pyfits.open('%s-G141_drz_first.fits' %(root))
    flat = pyfits.open('%s-G141_drz_flat.fits' %(root))
    sky = pyfits.open('%s-G141_drz_sky.fits' %(root))
    final = pyfits.open('%s-G141_drz_final.fits' %(root))
    
    im_shape = first[1].data.shape
    sup.im_shape = im_shape
    
    med = threedhst.utils.biweight(flat[1].data, mean=True)
    
    ysize = 3.
    #fig = plt.figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*4,ysize], dpi=100)
    
    top_panel = 0.2
    NPANEL = 4
    
    #plt.hot()
    plt.gray()
    plt.close()
    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['image.interpolation'] = 'nearest'
    
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*NPANEL*(1-top_panel)/2,ysize*2],dpi=100)
    else:
        fig = Figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*NPANEL*(1-top_panel)/2.,ysize*2], dpi=100)
    
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.02,
                        bottom=0.07,right=0.99,top=0.97)
    #
    
    plt.rcParams['lines.linewidth'] = 1
    
    vmin, vmax = -0.15, 0.075
    vmin, vmax= -0.08, 0.08
    
    x0 = 0.005*2
    y0 = x0/2.
    dx = (1.-(NPANEL+1)*x0)/NPANEL*2
    
    top_panel/=2.
    
    #ax = fig.add_subplot(141)
    ax = fig.add_axes(((x0+(dx+x0)*0), y0+0.5, dx, 0.5-top_panel-y0))
    
    ax.imshow((first[1].data-threedhst.utils.biweight(first[1].data, mean=True)), interpolation='nearest',aspect='auto',vmin=vmin-0.1*0,vmax=vmax+0.15*0)    
    sup.axis_imshow(ax, text='a)\ Raw')
    ax.text(0.12, 0.85, r'$\mathrm{%s}$' %(root), horizontalalignment='left', verticalalignment='center',
                 transform = ax.transAxes, color='black', fontsize=14)
    #
    #show_limits(ax, -(vmax+0.15)+med, -(vmin-0.1)+med)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*0), (0.5-top_panel)+0.5, dx, top_panel-2*y0))
    pp = sup.first_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.first_prof[i])
        pp += sup.first_prof[i]
    ax.plot(pp/4., color='black')
    sup.axis_profile(ax, yrange=range1, text='a)\ Raw')
    
    #ax = fig.add_subplot(142)
    ax = fig.add_axes(((x0+(dx+x0)*1)+x0, y0+0.5, dx, 0.5-top_panel-y0))
    ax.imshow((flat[1].data-med), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    sup.axis_imshow(ax, text='b)\ Flat')
    #show_limits(ax, -vmax+med, -vmin+med)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*1)+x0, (0.5-top_panel)+0.5, dx, top_panel-2*y0))
    pp = sup.flat_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.flat_prof[i])
        pp += sup.flat_prof[i]
    ax.plot(pp/4.+1, color='black')
    sup.axis_profile(ax, yrange=range1, text='b)\ Flat')
    
    ###########
    #ax = fig.add_subplot(143)
    ax = fig.add_axes(((x0+(dx+x0)*0), y0, dx, 0.5-top_panel-y0))
    ax.imshow((sky[1].data-med), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    sup.axis_imshow(ax, text='c)\ Background')
    #show_limits(ax, -vmax+med, -vmin+med)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*0), (0.5-top_panel), dx, top_panel-2*y0))
    pp = sup.sky_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.sky_prof[i])
        pp += sup.sky_prof[i]
    ax.plot(pp/4., color='black')
    sup.axis_profile(ax, yrange=range2, text='c)\ Background')
    
    #ax = fig.add_subplot(144)
    ax = fig.add_axes(((x0+(dx+x0)*1)+x0, y0, dx, 0.5-top_panel-y0))
    ax.imshow((final[1].data), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    sup.axis_imshow(ax, text='d)\ Final')
    #show_limits(ax, -vmax, -vmin)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*1)+x0, (0.5-top_panel), dx, top_panel-2*y0))
    pp = sup.final_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.final_prof[i])
        pp += sup.final_prof[i]
    ax.plot(pp/4., color='black')
    sup.axis_profile(ax, yrange=range2, text='d)\ Final')
    
    outfile = '%s-G141_demo.pdf' %(root)
    
    if USE_PLOT_GUI:
        fig.savefig(outfile,dpi=100,transparent=False)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)

def show_limits(ax, vmin, vmax):
    ax.text(0.98, -0.02,r'$\mathrm{[%.1f,\ %.1f]}$' %(vmin, vmax),
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform = ax.transAxes, color='black', fontsize=12)

def axis_profile(ax, yrange=None, prof=None, text=''):
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xlim(0,1000)
    if yrange is not None:
        ax.set_ylim(yrange[0], yrange[1])
        ylimits = yrange
    else:
        ylimits = ax.get_ylim()
    
    #
    if text is  not '':
        ax.text(0.5, 0.06,r'$\mathrm{'+text+'}$',
                     horizontalalignment='center',
                     verticalalignment='bottom',
                     transform = ax.transAxes, color='black', fontsize=10)
    
    ax.text(0.98, 0.08,r'$\mathrm{[%.2f,\ %.2f]}$' %(ylimits[0], ylimits[1]),
                 horizontalalignment='right',
                 verticalalignment='bottom',
                 transform = ax.transAxes, color='black', fontsize=8)
    
        
def axis_imshow(ax, text='', shape=None):
    import numpy as np
    import unicorn.survey_paper as sup
    
    if shape is None:
        shape = sup.im_shape
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks([0,shape[1]])
    ax.set_xticklabels([])
    ytick = ax.set_yticks([0,shape[0]])
    #
    # if text is  not '':
    #     ax.text(0.5, 1.02,r'$\mathrm{'+text+'}$',
    #                  horizontalalignment='center',
    #                  verticalalignment='bottom',
    #                  transform = ax.transAxes, color='black', fontsize=12)

def compare_sky():
    """
    Make a figure showing the aXe default and derived sky images.
    """
    import unicorn.survey_paper as sup
    
    sco = pyfits.open('../CONF/sky_cosmos.fits')
    shi = pyfits.open('../CONF/sky_goodsn_hi.fits')
    slo = pyfits.open('../CONF/sky_goodsn_lo.fits')
    svh = pyfits.open('../CONF/sky_goodsn_vhi.fits')
    
    shape = sco[0].data.shape
    
    figsize = 6
    fig = Figure(figsize=[figsize,figsize], dpi=200)
    #fig = plt.figure(figsize=[figsize,figsize], dpi=100)
    
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.02,
                        bottom=0.02,right=0.98,top=0.98)

    ####### COSMOS
    ax = fig.add_subplot(221)
    ax.imshow(sco[0].data, interpolation='nearest',aspect='auto',vmin=0.95,vmax=1.05)    
    sup.axis_imshow(ax, shape=shape)
    
    ax.fill_between([850,950],[50,50],[150,150], color='white', alpha=0.8)
    ax.text(900/1014., 100./1014,r'$\mathrm{a)}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, color='black', fontsize=12)
    
    ####### GOODS-N Lo
    ax = fig.add_subplot(222)
    ax.imshow(slo[0].data, interpolation='nearest',aspect='auto',vmin=0.95,vmax=1.05)    
    sup.axis_imshow(ax, shape=shape)

    ax.fill_between([850,950],[50,50],[150,150], color='white', alpha=0.8)
    ax.text(900/1014., 100./1014,r'$\mathrm{b)}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, color='black', fontsize=12)

    ####### GOODS-N Hi
    ax = fig.add_subplot(223)
    ax.imshow(shi[0].data, interpolation='nearest',aspect='auto',vmin=0.95,vmax=1.05)    
    sup.axis_imshow(ax, shape=shape)
    
    ax.fill_between([850,950],[50,50],[150,150], color='white', alpha=0.8)
    ax.text(900/1014., 100./1014,r'$\mathrm{c)}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, color='black', fontsize=12)

    ####### GOODS-N Very hi
    ax = fig.add_subplot(224)
    ax.imshow(svh[0].data, interpolation='nearest',aspect='auto',vmin=0.95,vmax=1.05)    
    sup.axis_imshow(ax, shape=shape)
    
    ax.fill_between([850,950],[50,50],[150,150], color='white', alpha=0.8)
    ax.text(900/1014., 100./1014,r'$\mathrm{d)}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, color='black', fontsize=12)
    
    #### Done
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure('sky_backgrounds.pdf', dpi=100, transparent=False)

def axeFlat(flat_file='/research/HST/GRISM/3DHST/CONF/WFC3.IR.G141.flat.2.fits', wave=1.4e4):
    """
    Compute the aXe flat-field image at a specified wavelength.
    """
    flat = pyfits.open(flat_file)
    
    wmin = flat[0].header['WMIN']
    wmax = flat[0].header['WMAX']

    x = (wave-wmin)/(wmax-wmin)
    img = np.zeros((1014,1014), dtype='float')
    for i in range(len(flat)):
        img += flat[i].data*x**i
    
    return img
    
def get_flat_function(x=507, y=507, wave=np.arange(1.1e4,1.6e4,500), flat_file='/research/HST/GRISM/3DHST/CONF/WFC3.IR.G141.flat.2.fits'):
    
    #wave = np.arange(1.1e4, 1.6e4, .5e3)
    
    flat = pyfits.open(flat_file)
    
    wmin = flat[0].header['WMIN']
    wmax = flat[0].header['WMAX']

    xx = (wave-wmin)/(wmax-wmin)
    flat_func = xx*0.
    
    for i in range(len(flat)):
        flat_func += flat[i].data[y,x]*xx**i
    
    return flat_func    
    
def show_flat_function():
    
    wave = np.arange(1.05e4, 1.7e4, 250.)
    color='blue'
    for xi in range(50,951,50):
        print unicorn.noNewLine+'%d' %(xi)
        for yi in range(50,951,50):
            ffunc = unicorn.survey_paper.get_flat_function(x=xi, y=yi, wave=wave)
            ffunc /= np.interp(1.4e4, wave, ffunc)
            p = plt.plot(wave, ffunc , alpha=0.05, color=color)
        
        
def grism_flat_dependence():
    """
    Compute the higher order terms for the grism flat-field
    """
    
    import unicorn
    import threedhst
    
    # f140 = threedhst.grism_sky.flat_f140[1].data[5:-5, 5:-5]
    # 
    # flat = pyfits.open(unicorn.GRISM_HOME+'CONF/WFC3.IR.G141.flat.2.fits')
    # wmin, wmax = flat[0].header['WMIN'], flat[0].header['WMAX']
    # 
    # a0 = flat[0].data
    # 
    # lam = 1.1e4
    # x = (lam-wmin)/(wmax-wmin)
    # 
    # aX = a0*0.
    # for i,ext in enumerate(flat[1:]):
    #     print i
    #     aX += ext.data*x**(i+1)
        
    f105 = pyfits.open(os.getenv('iref')+'/uc72113oi_pfl.fits')[1].data[5:-5,5:-5]
    #f105 = pyfits.open(os.getenv('iref')+'/uc72113ni_pfl.fits')[1].data[5:-5,5:-5] # F098M
    f140 = pyfits.open(os.getenv('iref')+'/uc721143i_pfl.fits')[1].data[5:-5,5:-5]
    yi, xi= np.indices(f140.shape)
    death_star = (f140 < 0.65) & (xi < 390) & (yi < 80) & (xi > 330) & (yi > 30)
    REF = 'F140W'
    f160 = pyfits.open(os.getenv('iref')+'/uc721145i_pfl.fits')[1].data[5:-5,5:-5]
    
    #### Narrow bands
    #f140 = pyfits.open(os.getenv('iref')+'/PFL/uc72113si_pfl.fits')[1].data[5:-5,5:-5]
    #REF = 'F127M'
    
    #f140 = pyfits.open(os.getenv('iref')+'/PFL/uc721140i_pfl.fits')[1].data[5:-5,5:-5]
    #f160 = pyfits.open(os.getenv('iref')+'/PFL/uc721146i_pfl.fits')[1].data[5:-5,5:-5]
    
    plt.rcParams['patch.edgecolor'] = 'None'
    #plt.rcParams['font.size'] = 12

    plt.rcParams['image.origin'] = 'lower'
    plt.rcParams['image.interpolation'] = 'nearest'
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    xs = 8
    fig = plt.figure(figsize=(xs,xs/3.), dpi=100)
    fig.subplots_adjust(wspace=0.01,hspace=0.01,left=0.01, bottom=0.01,right=0.99,top=0.99)        
    
    vmin, vmax = 0.95, 1.05
    #vmin, vmax = 0.9, 1.1
    ### put scale within the box
    NX = 100
    y0, y1 = 1014-1.5*NX, 1014-1*NX
    y0 -= NX; y1 -= NX
    
    ### F140W flat
    
    textbbox = dict(facecolor='white', alpha=0.6, edgecolor='white')
    
    ## correct for pixel area map
    PIXEL_AREA = True
    if PIXEL_AREA:
        pam = pyfits.open(os.getenv('iref')+'/ir_wfc3_map.fits')[1].data
    else:
        pam = np.ones((1014,1014),dtype='float')
    
    ax = fig.add_subplot(131)
    ax.imshow(f140/pam, vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.text(50,950, REF, verticalalignment='top', fontsize=14, bbox=textbbox)
    #ax.text(50,950, REF, verticalalignment='top', fontsize=14)
    ax.set_xlim(0,1014)
    ax.set_ylim(0,1014)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    ### F105W/F140W, with label
    ratio = f105/f140
    label = 'F105W / ' + REF
    
    ratio = unicorn.survey_paper.axeFlat(wave=1.1e4)/unicorn.survey_paper.axeFlat(wave=1.6e4)
    label = r'aXe 1.1 $\mu$m / 1.6 $\mu$m'
    
    ratio[death_star] = 0.
    #### Color bar for label
    x0 = 300
    # ratio[y0:y1,x0-2.5*NX:x0-1.5*NX] = vmin
    # ratio[y0:y1,x0-1.5*NX:x0-0.5*NX] = (vmin+1)/2.
    # ratio[y0:y1,x0-0.5*NX:x0+0.5*NX] = 1.
    # ratio[y0:y1,x0+0.5*NX:x0+1.5*NX] = (vmax+1)/2.
    # ratio[y0:y1,x0+1.5*NX:x0+2.5*NX] = vmax
    
    NSPLIT = 5
    NXi = NX*2./NSPLIT
    for i in range(1,NSPLIT+1):
        #print i,NXi, 1+(vmin-1)*i/NSPLIT, x0-(i-0.5)*NXi
        ratio[y0:y1,x0-(i+0.5)*NXi:x0-(i-0.5)*NXi] = 1+(vmin-1)*i/NSPLIT
    #
    ratio[y0:y1,x0-0.5*NXi:x0+0.5*NXi] = 1
    
    for i in range(1,NSPLIT+1):
        #print i,NXi, 1+(vmin-1)*i/NSPLIT, x0-(i-0.5)*NXi
        ratio[y0:y1,x0+(i-0.5)*NXi:x0+(i+0.5)*NXi] = 1+(vmax-1)*i/NSPLIT
    
    xbox = np.array([0,1,1,0,0])*NXi
    ybox = np.array([0,0,1,1,0])*NX/2
    
    ax = fig.add_subplot(132)
    ax.imshow(ratio, vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.plot(xbox+x0-0.5*NXi, ybox+y1-0.5*NX, color='0.6', alpha=0.1)
    
    fs = 9
    ax.text(x0-2*NX, y0-0.5*NX, '%.2f' %(vmin), horizontalalignment='center', verticalalignment='center', fontsize=fs)
    ax.text(x0-0*NX, y0-0.5*NX, '%.2f' %(1), horizontalalignment='center', verticalalignment='center', fontsize=fs)
    ax.text(x0+2*NX, y0-0.5*NX, '%.2f' %(vmax), horizontalalignment='center', verticalalignment='center', fontsize=fs)

    ax.text(50,950, label, verticalalignment='top', fontsize=14, bbox=textbbox)
    ax.set_xlim(0,1014)
    ax.set_ylim(0,1014)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
        
    ### F160W/F140W
    ax = fig.add_subplot(133)
    
    ratio = f160/f140
    label = 'F160W / '+REF

    ratio = f105/f160
    label = 'F105W / F160W'

    # ratio = unicorn.survey_paper.axeFlat(wave=1.0552e4)/unicorn.survey_paper.axeFlat(wave=1.392e4)
    
    ratio[death_star] = 0.
    
    ax.imshow(ratio, vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.text(50,950, label, verticalalignment='top', fontsize=14,bbox=textbbox)
    ax.set_xlim(0,1014)
    ax.set_ylim(0,1014)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    fig.savefig('compare_flats_v2.pdf')

def process_sky_background():
    import threedhst.catIO as catIO
    import re
    
    info = catIO.Readfile('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/sky_background.dat')
    
    field = []
    for targ in info.targname:
        targ = targ.replace('GNGRISM','GOODS-NORTH-')
        field.append(re.split('-[1-9]',targ)[0].upper())
    
    field = np.array(field)   
    is_UDF = field == 'xxx'
    for i,targ in enumerate(info.targname):
        m = re.match('GOODS-SOUTH-3[4678]',targ)
        if m is not None:
            is_UDF[i] = True
            
    fields = ['AEGIS','COSMOS','GOODS-SOUTH','GOODS-NORTH','UDS']
    colors = ['red','blue','green','purple','orange']
    
    for i in range(len(fields)):
        match = (field == fields[i]) & (info.filter == 'G141')
        h = plt.hist(info.bg_mean[match], range=(0,5), bins=50, color=colors[i], alpha=0.5)
        bg = threedhst.utils.biweight(info.bg_mean[match], both=True)
        print '%-14s %.3f %.3f ' %(fields[i], bg[0], bg[1])
        
    match = is_UDF & (info.filter == 'G141')
    bg = threedhst.utils.biweight(info.bg_mean[match], both=True)
    print '%-14s %.3f %.3f ' %('HUDF09', bg[0], bg[1])
    
    plt.xlim(0,3.5)
    
def get_background_level():
    """
    Get the sky background levels from the raw FLT images, with an object and dq mask
    """
    
    xi, yi = np.indices((1014, 1014))
    DPIX = 300
    
    flat = pyfits.open(os.getenv('iref')+'uc721143i_pfl.fits')[1].data[5:-5,5:-5]
    
    fp = open(unicorn.GRISM_HOME+'ANALYSIS/sky_background.dat','w')
    fp.write('# file targname filter date_obs bg_mean bg_sigma\n')
    
    for path in ['AEGIS','COSMOS','GOODS-S','GOODS-N','UDS']:
        os.chdir(unicorn.GRISM_HOME+path+'/PREP_FLT')
        info = catIO.Readfile('files.info')
        print '\n\n%s\n\n' %(path)
        for ii, file in enumerate(info.file):
            file = file.replace('.gz','')
            print file
            #
            try:
                flt_raw = pyfits.open(threedhst.utils.find_fits_gz('../RAW/'+file))
                flt_fix = pyfits.open(threedhst.utils.find_fits_gz(file))
                seg = pyfits.open(threedhst.utils.find_fits_gz(file.replace('fits','seg.fits')))
            except:
                continue
            #    
            mask = (seg[0].data == 0) & ((flt_fix[3].data & (4+32+16+512+2048+4096)) == 0) & (np.abs(xi-507) < DPIX) & (np.abs(yi-507) < DPIX)
            #
            if info.filter[ii].startswith('G'):
                flt_raw[1].data /= flat
            #    
            background = threedhst.utils.biweight(flt_raw[1].data[mask], both=True)
            fp.write('%s   %16s   %5s   %s   %.3f %.3f\n' %(file, info.targname[ii], info.filter[ii], info.date_obs[ii], background[0], background[1]))
        
    fp.close()
    
def get_spec_signal_to_noise():
    """
    Measure the S/N directly from the spectrum of all objects.
    """
    
    fp = open(unicorn.GRISM_HOME+'ANALYSIS/spec_signal_to_noise.dat','w')
    fp.write('# object mag_auto flux_radius sig_noise n_bins\n')

    for path in ['AEGIS','COSMOS','GOODS-S','GOODS-N','UDS']:
        os.chdir(unicorn.GRISM_HOME+path)
        SPC_files = glob.glob('DRIZZLE_G141/*opt.SPC.fits')
        for file in SPC_files:
            pointing = os.path.basename(file).split('_2_opt')[0]
            #
            SPC = threedhst.plotting.SPCFile(file,axe_drizzle_dir='./')
            cat = threedhst.sex.mySexCat('DATA/%s_drz.cat' %(pointing))
            try:
                mag_auto = np.cast[float](cat.MAG_F1392W)
                flux_radius = np.cast[float](cat.FLUX_RADIUS)
            except:
                continue
            #
            for id in SPC._ext_map:
                print unicorn.noNewLine+'%s_%05d' %(pointing, id)
                #
                spec = SPC.getSpec(id)
                mask = (spec.LAMBDA > 1.15e4) & (spec.LAMBDA < 1.6e4) & (spec.CONTAM/spec.FLUX < 0.1) & np.isfinite(spec.FLUX)
                if len(mask[mask]) > 1:
                    signal_noise = threedhst.utils.biweight(spec.FLUX[mask]/spec.FERROR[mask], mean=True)
                    mat = cat.id == id
                    fp.write('%s_%05d %8.3f %8.3f   %13.3e %-0d\n' %(pointing, id, mag_auto[mat][0], flux_radius[mat][0], signal_noise, len(mask[mask])))
                #   
                else:
                    continue
    
    fp.close()
          
def process_signal_to_noise():
    import threedhst.catIO as catIO
    import re

    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/')
    
    info = catIO.Readfile('spec_signal_to_noise.dat')

    field = []
    for targ in info.object:
        targ = targ.replace('GNGRISM','GOODS-NORTH-')
        field.append(re.split('-[1-9]',targ)[0].upper())

    field = np.array(field)   
    
    snscale = 2.5 ### aXe errors too large
    
    ma = 'o'
    ms = 5
    
    ##### S/N vs mag.
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=1, left=0.10)
    fig.subplots_adjust(wspace=0.2,hspace=0.24,left=0.09, bottom=0.08,right=0.975,top=0.99)
    
    ax = fig.add_subplot(211)
    
    ff = (field == 'COSMOS') & (info.sig_noise > 0)
    ax.plot(info.mag_auto[ff], info.sig_noise[ff]*snscale, marker=ma, markersize=ms, color='red', alpha=0.1, linestyle='None')
    xm, ym, ys, ns = threedhst.utils.runmed(info.mag_auto[ff], info.sig_noise[ff]*snscale, NBIN=30)

    ff = (field == 'AEGIS') & (info.sig_noise > 0)
    ax.plot(info.mag_auto[ff], info.sig_noise[ff]*snscale, marker=ma, markersize=ms, color='blue', alpha=0.1, linestyle='None')

    ax.plot(xm,ym,color='white',linewidth=6,alpha=0.5)
    ax.plot(xm,ym,color='red',linewidth=3,alpha=0.9)

    xm, ym, ys, ns = threedhst.utils.runmed(info.mag_auto[ff], info.sig_noise[ff]*snscale, NBIN=30)
    ax.plot(xm,ym,color='white',linewidth=6,alpha=0.5)
    ax.plot(xm,ym,color='blue',linewidth=3,alpha=0.9)

    ax.semilogy()
    ax.plot([12,30],[3,3], linewidth=3, alpha=0.4, color='black', linestyle='--')
    ax.set_ylim(0.2,300)
    ax.set_xlim(16,25)
    ax.set_yticklabels(['1','3','10','100']) #; ax.set_xticklabels([])
    ytick = ax.set_yticks([1,3,10,100]) #; xtick = ax.set_xticks([0,NX]); 
    if plt.rcParams['text.usetex']:
        ax.set_xlabel('MAG\_AUTO (F140W)')
    else:
        ax.set_xlabel('MAG_AUTO (F140W)')
                
    ax.set_ylabel('S / N')
    
    ax.text(16.5,1,'AEGIS',color='blue', fontsize=12)
    ax.text(16.5,0.5,'COSMOS',color='red', fontsize=12)
    
    ##### S/N vs size for a mag bin
    ax = fig.add_subplot(212)
    
    m0, m1 = 22., 22.5
    
    mm = (info.mag_auto > m0) & (info.mag_auto < m1) & (info.sig_noise > 0)
    
    ax.plot(info.flux_radius[mm & (field == 'COSMOS')], info.sig_noise[mm & (field == 'COSMOS')]*snscale, marker=ma, markersize=ms, color='red', alpha=0.3, linestyle='None')
    xm, ym, ys, ns = threedhst.utils.runmed(info.flux_radius[mm & (field == 'COSMOS')], info.sig_noise[mm & (field == 'COSMOS')]*snscale, NBIN=5)
    
    ax.plot(info.flux_radius[mm & (field == 'AEGIS')], info.sig_noise[mm & (field == 'AEGIS')]*snscale, marker=ma, markersize=ms, color='blue', alpha=0.3, linestyle='None')

    ax.plot(xm,ym,color='white',linewidth=6,alpha=0.5)
    ax.plot(xm,ym,color='red',linewidth=3,alpha=0.9)
    
    xm, ym, ys, ns = threedhst.utils.runmed(info.flux_radius[mm & (field == 'AEGIS')], info.sig_noise[mm & (field == 'AEGIS')]*snscale, NBIN=5)
    ax.plot(xm,ym,color='white',linewidth=6,alpha=0.5)
    ax.plot(xm,ym,color='blue',linewidth=3,alpha=0.9)
    #plt.plot(xm,ym/np.sqrt(2),color='blue',linewidth=3,alpha=0.5)
    
    ax.semilogy()
    #ax.plot([0,30],[3,3], linewidth=3, alpha=0.5, color='black')
    ax.set_ylim(2,15)
    ax.set_xlim(1.5,12)
    ax.set_yticklabels(['3','5','10']) #; ax.set_xticklabels([])
    ytick = ax.set_yticks([3,5,10]) #; xtick = ax.set_xticks([0,NX]); 
    ax.set_xlabel(r'R$_{50}$ [$0.06^{\prime\prime}$ pix]')
    
    ax.set_ylabel('S / N')
    if plt.rcParams['text.usetex']:
        ax.text(11.8,13, r'$%.1f < H_{140} < %.1f$' %(m0, m1), horizontalalignment='right', verticalalignment='top')
    else:
        ax.text(11.8,13, r'%.1f < $m_{140}$ < %.1f' %(m0, m1), horizontalalignment='right', verticalalignment='top')
        
    fig.savefig('spec_signal_to_noise.pdf')
    plt.rcParams['text.usetex'] = False
    
def clash_empty_apertures():
    
    for cluster in ['a2261','a383','macs1149','macs1206','macs2129']:
        os.chdir('/Users/gbrammer/CLASH/%s' %(cluster))
        images = glob.glob('*drz.fits')
        for image in images:
            wht = image.replace('drz','wht')
            head = pyfits.getheader(image)
            zp=-2.5*np.log10(head['PHOTFLAM']) - 21.10 - 5 *np.log10(head['PHOTPLAM']) + 18.6921 
            unicorn.candels.clash_make_rms_map(image=wht, include_poisson=False)
            unicorn.survey_paper.empty_apertures(SCI_IMAGE=image, SCI_EXT=0, WHT_IMAGE=wht.replace('wht','rms'), WHT_EXT=0, aper_params=(0.4/0.065/2.,0.4/0.065/2.+1,2), ZP=zp, make_plot=False, NSIM=1000, MAP_TYPE='MAP_RMS')
            
    #### Sequence of apertures for measuring Beta
    for cluster in ['a2261','a383','macs1149','macs1206','macs2129']:
        os.chdir('/Users/gbrammer/CLASH/%s' %(cluster))
        images = glob.glob('*_f160w*_drz.fits')
        for image in images:
            wht = image.replace('drz','wht')
            head = pyfits.getheader(image)
            zp=-2.5*np.log10(head['PHOTFLAM']) - 21.10 - 5 *np.log10(head['PHOTPLAM']) + 18.6921 
            #unicorn.candels.clash_make_rms_map(image=wht, include_poisson=False)
            unicorn.survey_paper.empty_apertures(SCI_IMAGE=image, SCI_EXT=0, WHT_IMAGE=wht.replace('wht','rms'), WHT_EXT=0, aper_params=(0.2/0.065,3.05/0.065,0.2/0.065), ZP=zp, make_plot=False, NSIM=500, MAP_TYPE='MAP_RMS')
    
    
    for cluster in ['a2261','a383','macs1149','macs1206','macs2129'][:-1]:
        os.chdir('/Users/gbrammer/CLASH/%s' %(cluster))
        files = glob.glob('*empty.fits')
        print '\n--------------\n%s\n--------------\n' %(cluster.center(14))
        for file in files:
            head = pyfits.getheader(file.replace('_empty',''))
            zp=-2.5*np.log10(head['PHOTFLAM']) - 21.10 - 5 *np.log10(head['PHOTPLAM']) + 18.6921 
            em = pyfits.open(file)
            print '%-7s %.2f' %(file.split('_')[5], zp-2.5*np.log10(5*np.std(em[2].data)))
        
def run_empty_apertures_fields():
    import glob
    import os
    import unicorn
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/EMPTY_APERTURES/')
    
    files = glob.glob('/3DHST/Spectra/Work/COSMOS/PREP_FLT/COSMOS-*-F140W_drz.fits')
    files = glob.glob('/3DHST/Spectra/Work/GOODS-N/PREP_FLT/GOODS-N-*-F140W_drz.fits')
    
    for file in files[1:]:
        unicorn.survey_paper.empty_apertures(SCI_IMAGE=file, SCI_EXT=1, WHT_IMAGE=file, WHT_EXT=2, aper_params=(1,17,1), NSIM=1000, ZP=26.46, make_plot=True)
    
    
    
def grism_apertures_plot(SHOW_FLUX=False):
    
    test = """
    The following derives the correction, R, needed to scale the pixel standard
    deviations in the drizzled images with small pixels, to the noise in 
    "nominal" pixels.
    
    files = glob.glob('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/EMPTY_APERTURES/*pix*empty.fits')
    
    ratio = []
    for file in files:
        impix = pyfits.open(file)
        im = pyfits.open(file.replace('pix',''))
        print threedhst.utils.biweight(im[2].data[:,2]) /  threedhst.utils.biweight(impix[2].data[:,2])
        ratio.append(threedhst.utils.biweight(im[2].data[:,2]) /  threedhst.utils.biweight(impix[2].data[:,2]))
    
    print 'R ~ %.2f' %(1./np.mean(ratio))
    
    # Emission line flux
    
    texp, Nexp, R, dark, sky, ee_fraction = 1277, 4, 2, 0.05, 1.5, 0.75
    area = np.pi*R**2
    
    total_counts_resel = (xarr+dark)*texp*area*Nexp # 
    eq1_cts = np.sqrt(total_counts_resel+rn**2*area*Nexp)
    eq1_cps = eq1_cts/texp/Nexp
    #eq1_flam = eq1_cps/(sens_14um*46.5)
    print eq1_cps
    
    """
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/SURVEY_PAPER')
    
    # files = glob.glob('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/EMPTY_APERTURES_CIRCULAR/*G141_drz_empty.fits')
    # aper_use = 0
    
    files = glob.glob('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/EMPTY_APERTURES/*G141_drz_empty.fits')
    aper_use = 1

    #### Parameters of the plot
    lam_int = 1.4e4
    SN_show = 5
    
    if SHOW_FLUX:
        aper_use = 0
        
    bg = threedhst.catIO.Readfile('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/sky_background.dat')
    
    sens = pyfits.open(unicorn.GRISM_HOME+'CONF/WFC3.IR.G141.1st.sens.2.fits')[1].data
    sens_14um = np.interp(lam_int,sens.WAVELENGTH,sens.SENSITIVITY)
    
    colors = {'AEGIS':'red','COSMOS':'blue','UDS':'purple', 'GNGRISM':'orange','GOODSS':'green'}
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(xs=4, left=0.125, bottom=0.085, right=0.01, top=0.01, square=True, fontsize=12)
    
    ax = fig.add_subplot(111)
    
    for file in files:
        root = os.path.basename(file).split('-G141')[0].replace('GOODS-N-','GNGRISM').replace('S-S','S-SOUTH')
        print root
        mat = (bg.targname == root) & (bg.filter == 'G141')
        if len(mat[mat]) == 0:
            continue
        mean_bg = np.mean(bg.bg_mean[mat])
        err_bg = np.std(bg.bg_mean[mat])
        aps = pyfits.open(file)
        #
        fluxes = aps[2].data
        stats = threedhst.utils.biweight(fluxes, both=True)
        sigmas = fluxes[0,:].flatten()*0
        for i in range(len(sigmas)):
            sigmas[i] = threedhst.utils.biweight(fluxes[:,i])
            #sigmas[i] = threedhst.utils.nmad(fluxes[:,i])
            #sigmas[i] = np.std(fluxes[:,i])

        #sigmas *= 1.34  ### scale factor for correlated pixels, determined empirically
        inv_sens_flam = sigmas/(sens_14um*22*4)
        
        field = os.path.basename(file).replace('GOODS-N','GNGRISM').replace('GOODS-S','GOODSS').split('-')[0]
        print field, aps[1].data[aper_use]
        
        if SHOW_FLUX:
            sig3_flux = inv_sens_flam*SN_show
            sig3_flux /= 0.75**2  ### Include for encircled energy within 3 nominal spatial pixles
            sig3_ab = sig3_flux*2*46.5/1.e-17
        else:
            inv_sens_fnu = inv_sens_flam*lam_int**2/3.e18
            sig3_flux = inv_sens_fnu*SN_show
            sig3_flux /= 0.75  ### Include for encircled energy within 3 nominal spatial pixles
            sig3_ab = -2.5*np.log10(sig3_flux)-48.6#-2.5*np.log10(np.sqrt(2))
        
        p = ax.errorbar(mean_bg,sig3_ab[aper_use], xerr=err_bg, marker='o', ms=8, alpha=0.4, color=colors[field], ecolor=colors[field])
                
        
    xarr = np.arange(0,5,0.02) ###  background rate, cts / s
    yarr = 22.9-2.5*np.log10(np.sqrt(xarr/2.))
    
    #### Eq. 1 from paper
    scale = 0.06/0.128254
    Nexp, texp, rn, R, dark = 4, 1277, 20, 2*scale, 0.05
    
    #dark += 0.181/3.
    #area = np.pi*R**2
    area = 3
    ee_fraction = 0.75 # for 3 pix aperture in spatial direction, nominal pixels
    
    total_counts_resel = (xarr+dark)*texp*area*2*Nexp # 
    eq1_cts = np.sqrt(total_counts_resel+rn**2*area*Nexp)
    eq1_cps = eq1_cts/texp/Nexp/ee_fraction/2
    eq1_flam = eq1_cps/(sens_14um*46.5)
    eq1_fnu = eq1_flam*lam_int**2/3.e18
    eq1_ab = -2.5*np.log10(eq1_fnu*SN_show)-48.6#-2.5*np.log10(np.sqrt(0.5))
    
    plt.plot(xarr, eq1_ab+2.5*np.log10(1.35), color='black', alpha=0.4, linewidth=2, linestyle='--')
    plt.plot(xarr, eq1_ab, color='black', alpha=0.4, linewidth=2)
    plt.text(1.,24.21-2.5*np.log10(SN_show/3.),'ETC', rotation=-44, color='black', alpha=0.4,horizontalalignment='center')
    plt.arrow(1.5, 24.1-2.5*np.log10(SN_show/3.), 0, 0.1, color='0.6', alpha=1, fill=True, width=0.02, head_width=0.06, head_length=0.02, overhang=0.05)
    plt.text(1.5+0.05, 24.1-2.5*np.log10(SN_show/3.)+0.05,r'$R_\mathrm{driz}$', verticalalignment='center', alpha=0.4)
    
    #plt.plot(xarr, eq1_ab-2.5*np.log10(np.sqrt(2*46.5/22.5)), color='purple', linewidth=2)
    
    ### Predict source counts
    mag = 23
    source_fnu = 10**(-0.4*(mag+48.6))
    source_flam = source_fnu/lam_int**2*3.e18
    source_cps = source_flam*sens_14um*46.5*ee_fraction
    source_counts = source_cps*Nexp*texp*2
    
    print np.interp(1.88,xarr,total_counts_resel), np.interp(1.88,xarr,np.sqrt(total_counts_resel+rn**2*area*Nexp)), np.interp(1.88,xarr,eq1_ab), source_counts, source_counts / np.interp(1.88,xarr,np.sqrt(total_counts_resel+rn**2*area*Nexp))
    
    #ax.plot(xarr, yarr, color='black', alpha=0.2, linewidth=3)
    ax.set_xlim(0.6,3.2)
    if SHOW_FLUX:
        ax.set_ylim(1,5.1)
        ax.set_ylabel(r'$%0d\sigma$ emission line sensitivity ($10^{-17}\,\mathrm{erg\,s^{-1}\,cm^{-2}}$)' %(SN_show))
        #ax.semilogy()
    else:
        ax.set_ylim(23.8-2.5*np.log10(SN_show/3.),25.0-2.5*np.log10(SN_show/3.))
        ax.set_ylabel(r'$%0d\sigma$ continuum depth (1.4$\mu$m, $\Delta=92\,$\AA)' %(SN_show))
        
    ax.set_xlabel(r'Background level [electrons / s]')
    #ax.set_ylabel(r'$3\sigma$ continuum depth @ 1.4$\mu$m, $D_\mathrm{ap}=0.24^{\prime\prime}/\,90\,$\AA')
    
    x0, y0, dy = 3, 24.8-2.5*np.log10(SN_show/3.), 0.1
    for i,field in enumerate(colors.keys()):
        field_txt = field.replace('GNGRISM','GOODS-N').replace('GOODSS','GOODS-S')
        ax.text(x0, y0-i*dy, field_txt, color=colors[field], horizontalalignment='right')
     
    plt.savefig('grism_empty_apertures.pdf')   
    plt.close()
    
def grism_empty_apertures():
    """
    Try simple empty apertures routine to measure depth of grism exposures, 
    to compare with the values measured directly from the spectra.
    """
    
    unicorn.survey_paper.empty_apertures(SCI_IMAGE='GOODS-S-34-G141_drz.fits', WHT_IMAGE='GOODS-S-34-G141_drz.fits', aper_params=(2,8.1,2), NSIM=500, ZP=25, make_plot=False, verbose=True, threshold=0.8, is_grism=True)
    
    for field in ['AEGIS','GOODS-S','UDS']: #,'COSMOS','GOODS-N']:
        os.chdir(unicorn.GRISM_HOME+field+'/PREP_FLT/')
        images = glob.glob(field+'*[0-9]-G141_drz.fits')
        print images
        for image in images[1:]:
            unicorn.survey_paper.empty_apertures(SCI_IMAGE=image, WHT_IMAGE=image, aper_params=(2,8.1,2), NSIM=1000, ZP=25, make_plot=False, verbose=True, threshold=0.8, is_grism=True, rectangle_apertures = [(4,4),(2,6),(4,6)])
            #
            #### Make a version with nominal pixels for testing
            threedhst.shifts.make_grism_shiftfile(image.replace('drz','asn').replace('G141','F140W'), image.replace('drz','asn'))
            threedhst.utils.combine_asn_shifts([image.replace('drz','asn')], out_root=image.split('_drz')[0]+'pix')
            threedhst.prep_flt_files.startMultidrizzle(image.split('_drz')[0] +'pix_asn.fits',
                     use_shiftfile=True, skysub=False,
                     final_scale=0.128254, pixfrac=0.01, driz_cr=False,
                     updatewcs=False, clean=True, median=False)
            new = image.replace('G141','G141pix')
            unicorn.survey_paper.empty_apertures(SCI_IMAGE=new, WHT_IMAGE=new, aper_params=(2,8.1,2), NSIM=500, ZP=25, make_plot=False, verbose=True, threshold=1.0, is_grism=True, rectangle_apertures = [(2,2),(1,3),(2,3),(4,4)])
            
    aps = pyfits.open('GOODS-S-34-G141_drz_empty.fits')
    fluxes = aps[2].data.flatten()
    stats = threedhst.utils.biweight(fluxes, both=True)
    
    sens = pyfits.open('../../CONF/WFC3.IR.G141.1st.sens.2.fits')[1].data
    wave = sens.WAVELENGTH
    inv_sens_flam = 1./(sens.SENSITIVITY*0.06/0.128254*46.5)
    inv_sens_fnu = inv_sens_flam*wave**2/3.e18
    
    sig3 = inv_sens_fnu*3*stats[1]
    sig3_ab = -2.5*np.log10(sig5)-48.6
    plt.plot(wave, sig5_ab)
    
    mag = sig5_ab*0+23
    input_fnu = 10**(-0.4*(mag+48.6))
    input_flam = input_fnu*3.e18/wave**2
    input_counts = input_flam * sens.SENSITIVITY * 46.5
    
    #plt.plot(sens.WAVELENGTH, sens.SENSITIVITY)
    
def empty_apertures(SCI_IMAGE='PRIMO_F125W_drz.fits', SCI_EXT=1, WHT_IMAGE='PRIMO_F125W_drz.fits', WHT_EXT=2, aper_params=(1,17,0.5), NSIM=1000, ZP=26.25, make_plot=True, verbose=True, MAP_TYPE='MAP_WEIGHT', threshold=1.5, is_grism=False, rectangle_apertures = None):
    """
    1) Run SExtractor on the input image to generate a segmentation map.
    
    2) Place `NSIM` empty apertures on the image, avoiding objects and areas with 
       zero weight as defined in the `WHT_IMAGE`.  
       
       The list of aperture radii used are np.arange(`aper_params`).
    
    3) Store the results in a FITS file `SCI_IMAGE`_empty.fits.
    
    Circular apertures are placed on the science image with the fractional pixel 
    coverage determined with a polygon approximation, and this can take a while.
    
    """
    
    from shapely.geometry import Point, Polygon
    
    if SCI_EXT == 0:
        SCI_EXT_SEX = 1
    else:
        SCI_EXT_SEX = SCI_EXT
        
    if WHT_EXT == 0:
        WHT_EXT_SEX = 1
    else:
        WHT_EXT_SEX = WHT_EXT
        
    ROOT = os.path.basename(SCI_IMAGE).split('.fits')[0]
    
    #### Open the science image
    if verbose:
        print 'Read images...'
        
    img = pyfits.open(SCI_IMAGE)
    img_data = img[SCI_EXT].data
    img_head = img[SCI_EXT].header
    img_shape = img_data.shape
    
    #### Setup SExtractor and run to generate a segmentation image
    if is_grism:
        threedhst.sex.USE_CONVFILE = 'grism.conv'
    else:
        threedhst.sex.USE_CONVFILE = 'gauss_4.0_7x7.conv'
    
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = '%s_empty.cat' %(ROOT)
    se.options['CHECKIMAGE_NAME'] = '%s_empty_seg.fits' %(ROOT)
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    
    if WHT_IMAGE is None:
        se.options['WEIGHT_TYPE']     = 'NONE'
        img_wht = img_data*0.
        img_wht[img_data != 0] = 1
    else:
        se.options['WEIGHT_TYPE']     = MAP_TYPE
        se.options['WEIGHT_IMAGE']    = '%s[%d]' %(WHT_IMAGE, WHT_EXT_SEX-1)
        wht = pyfits.open(WHT_IMAGE)
        img_wht = wht[WHT_EXT].data
    
    ##### Needed for very faint limits
    se.options['MEMORY_OBJSTACK'] = '8000'
    se.options['MEMORY_PIXSTACK'] = '800000'
    
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '%.2f' %(threshold)
    se.options['ANALYSIS_THRESH']  = '%.2f' %(threshold)
    se.options['MAG_ZEROPOINT'] = '%.2f' %(ZP)  ### arbitrary, actual mags don't matter
    status = se.sextractImage('%s[%d]' %(SCI_IMAGE, SCI_EXT_SEX-1))
    
    #### Read the Segmentation image
    segim = pyfits.open('%s_empty_seg.fits' %(ROOT))
    seg = segim[0].data
    segim.close()
    
    #### Set up the apertures
    
    #NSIM = 1000
    #apertures = np.arange(1,17,0.5)
    apertures = np.arange(aper_params[0], aper_params[1], aper_params[2])
    if rectangle_apertures is not None:
        IS_RECTANGLE = True
        apertures = rectangle_apertures   ### list of tuples with (dx,dy) sizes
    else:
        IS_RECTANGLE = False
    
    fluxes = np.zeros((NSIM, len(apertures)))
    centers = np.zeros((NSIM, len(apertures), 2))
    
    #### Loop throuth the desired apertures and randomly place NSIM of them
    aper = np.zeros(img_shape, dtype=np.float)
    for iap, ap in enumerate(apertures):
        #aper_image = np.zeros(img_shape)
        icount = 0
        if IS_RECTANGLE:
            print 'Aperture %.2f x %.2f pix\n' %(ap[0], ap[1])
            rap = (ap[0]/2.,ap[1]/2.)
        else:
            print 'Aperture radius: %.2f pix\n' %(ap)
            rap = (ap, ap)
            
        while icount < NSIM:
            #### Random coordinate
            xc = np.random.rand()*(img_shape[1]-4*rap[0])+2*rap[0]
            yc = np.random.rand()*(img_shape[0]-4*rap[1])+2*rap[1]
            
            #### Quick test to see if the coordinate is within an object or 
            #### where weight is zero
            if (seg[int(yc), int(xc)] != 0) | (img_wht[int(yc), int(xc)] <= 0) | (img_data[int(yc), int(xc)] == 0):
                continue
            
            #### Shapely point + buffer to define the circular aperture
            if IS_RECTANGLE:
                aperture_polygon = Polygon(((xc-ap[0]/2.,yc-ap[1]/2.), (xc+ap[0]/2.,yc-ap[1]/2.), (xc+ap[0]/2.,yc+ap[1]/2.), (xc-ap[0]/2.,yc+ap[1]/2.)))
            else:
                point = Point(xc, yc)
                aperture_polygon = point.buffer(ap, resolution=16)
            
            #### initialize the aperture
            aper*=0
            
            #### Loop through pixels to compute fractional pixel coverage within
            #### the circular aperture using the intersection of Shapely polygons
            smax = 0
            wmin = 1.e10
            for i in range(int(np.floor(xc-rap[0])),int(np.ceil(xc+rap[0]))):
                for j in range(int(np.floor(yc-rap[1])),int(np.ceil(yc+rap[1]))):
                    pix = Polygon(((i+0.5,j+0.5), (i+1.5,j+0.5), (i+1.5,j+1.5), (i+0.5,j+1.5)))
                    isect = pix.intersection(aperture_polygon)
                    aper[j,i] = isect.area
                    if isect.area > 0:
                        smax = np.array([smax, seg[j,i]]).max()
                        wmin = np.array([wmin, img_wht[j,i]]).min()
                        
            #### Only keep the result if the aperture doesn't intersect with an object
            #### as defined in the segmention image and if all weights within the 
            #### aperture are greater than zero
            if (smax == 0) & (wmin > 0):
                fluxes[icount, iap] = (aper*img_data).sum()
                centers[icount, iap, : ]  = np.array([xc, yc])
                #aper_image += aper
                print unicorn.noNewLine+'%d' %(icount)
                icount += 1
            else:
                print unicorn.noNewLine+'Skip: %f %f' %((seg*aper).max(), (img_wht*aper).min())
                continue
    
    #### Make the output FITS file.  List of aperture radii in extension 1, aperture
    #### fluxes in extension 2.
    
    ap_head = pyfits.Header()
    ap_head.update('NSIM',NSIM, comment='Number of apertures')
    ap_head.update('SCI_IMG',SCI_IMAGE, comment='Science image')
    ap_head.update('SCI_EXT',SCI_EXT, comment='Science extension')
    ap_head.update('WHT_IMG',WHT_IMAGE, comment='Weight image')
    ap_head.update('WHT_EXT',WHT_EXT, comment='Weight extension')
    
    prim = pyfits.PrimaryHDU(header=ap_head)
    ap_hdu = pyfits.ImageHDU(data=np.array(apertures))
    fl_hdu = pyfits.ImageHDU(data=fluxes)
    ce_hdu = pyfits.ImageHDU(data=centers)
    pyfits.HDUList([prim, ap_hdu, fl_hdu, ce_hdu]).writeto('%s_empty.fits' %(ROOT), clobber='True')
    
    if make_plot is True:
        make_empty_apertures_plot(empty_file='%s_empty.fits' %(ROOT), ZP=ZP)
        
def make_empty_apertures_plot(empty_file='PRIMO_F125W_drz_empty.fits', ZP=26.25, NSIG=5):
    """
    Plot the results from the `empty_apertures` routine.
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import polyfit, polyval

    import threedhst
    import unicorn
    
    im = pyfits.open(empty_file)
    apertures = im[1].data
    fluxes = im[2].data
    
    ROOT = empty_file.split('_empty')[0]
    
    sigma = apertures*0.
    means = apertures*0.
    for iap, ap in enumerate(apertures):
        sigma[iap] = threedhst.utils.biweight(fluxes[:,iap])
        means[iap] = threedhst.utils.biweight(fluxes[:,iap], mean=True)
    
    #plt.plot(apertures, means/np.pi/apertures**2)
    #plt.ylim(0,1.5*np.max(means/np.pi/apertures**2))
    #threedhst.utils.biweight(img_data[(seg == 0) & (img_wht > 0)], both=True)

    fig = unicorn.catalogs.plot_init(xs=6, aspect=0.5, left=0.12)
    
    fig.subplots_adjust(wspace=0.24,hspace=0.0,left=0.095,
                        bottom=0.17,right=0.97,top=0.97)
    
    ax = fig.add_subplot(121)
    
    ################################## plot sigma vs radius  
    ax.plot(apertures, sigma, marker='o', linestyle='None', color='black', alpha=0.8, markersize=5)
    
    coeffs = polyfit(np.log10(apertures), np.log10(sigma), 1)
    ax.plot(apertures, 10**polyval(coeffs, np.log10(apertures)), color='red')
    
    xint = 1
    y2 = apertures**2
    y2 = y2*np.interp(xint,apertures,sigma) / np.interp(xint,apertures,y2)
    ax.plot(apertures, y2, linestyle='--', color='black', alpha=0.3)         
    
    y1 = apertures**1
    y1 = y1*np.interp(xint,apertures,sigma) / np.interp(xint,apertures,y1)
    ax.plot(apertures, y1, linestyle='--', color='black', alpha=0.3)         
    
    ax.set_xlabel(r'$R_\mathrm{aper}$ [pix]')
    ax.set_ylabel(r'$\sigma_\mathrm{biw}$')
    #ax.text(apertures.max(), 0.1*sigma.max(), r'$\beta=%.2f$' %(coeffs[0]), horizontalalignment='right')
    ax.text(0.08, 0.85, r'$N_\mathrm{ap}=%d$' %(im[0].header['NSIM']), transform=ax.transAxes)
    ax.set_ylim(0,2)
    
    ################################# Plot AB depth
    ax = fig.add_subplot(122)
    
    ax.plot(apertures, ZP-2.5*np.log10(sigma*NSIG), marker='o', linestyle='None', color='black', alpha=0.8, markersize=5)
    
    ax.plot(apertures, ZP-2.5*np.log10(10**polyval(coeffs, np.log10(apertures))*NSIG), color='red')
    
    ax.plot(apertures, ZP-2.5*np.log10(y2*NSIG), linestyle='--', color='black', alpha=0.3)         
    
    ax.plot(apertures, ZP-2.5*np.log10(y1*NSIG), linestyle='--', color='black', alpha=0.3)         
    
    ax.set_xlabel(r'$R_\mathrm{aper}$ [pix]')
    ax.set_ylabel(r'Depth AB mag (%d$\sigma$)' %(NSIG))
    
    ax.text(0.95, 0.9, ROOT, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    
    ax.text(0.95, 0.8, r'$\beta=%.2f$' %(coeffs[0]), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    
    ax.set_ylim(23, 30)
    
    ################################## Save the result
    outfile = ROOT+'_empty.pdf'
    if USE_PLOT_GUI:
        fig.savefig(outfile,dpi=100,transparent=False)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)
    
    print ROOT+'_empty.pdf'

def show_background_flux_distribution():
    """
    Extract the SKYSCALE parameter from the G141 FLT images and plot their distribution
    by field.
    """
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    # for field in ['AEGIS','COSMOS','GOODS-S','UDS']:
    #     os.chdir(unicorn.GRISM_HOME+'%s/PREP_FLT' %(field))
    #     print field
    #     status = os.system('dfits *flt.fits |fitsort TARGNAME FILTER SKYSCALE |grep G141 > %s_skyscale.dat' %(field))
    # 
    # os.chdir(unicorn.GRISM_HOME)
    # status = os.system('cat AEGIS/PREP_FLT/AEGIS_skyscale.dat COSMOS/PREP_FLT/COSMOS_skyscale.dat GOODS-S/PREP_FLT/GOODS-S_skyscale.dat UDS/PREP_FLT/UDS_skyscale.dat > 3DHST_skyscale.dat') 
    
    ### 
    os.chdir("/research/HST/GRISM/3DHST/SIMULATIONS")
    sky = catIO.Readfile('3DHST_skyscale.dat')
    pointings = np.unique(sky.pointing)
    
    colors = {}
    colors['UDS'] = 'purple'; colors['COSMOS'] = 'blue'; colors['AEGIS'] = 'red'; colors['GOODS-SOUTH'] = 'green'
    
    off = {}
    off['UDS'] = 1
    nexp = np.arange(1,5)
    fields = []
    for pointing in pointings:
        this = sky.pointing == pointing
        field = '-'.join(pointing.split('-')[:-1])
        fields.extend([field]*4)
        #print pointing, field
        #p = plt.plot(sky.skyscale[this], nexp, marker='o', linestyle='-', alpha=0.5, color=colors[field])
    
    fields = np.array(fields)
    
    plt.ylim(0,5)
    
    fig = unicorn.catalogs.plot_init(xs=3.8, left=0.06, bottom=0.08, use_tex=True, fontsize=10)
    ax = fig.add_subplot(111)
    for i, field in enumerate(np.unique(fields)[::-1]):
        this = fields == field
        yh, xh = np.histogram(sky.skyscale[this], range=(0,4), bins=40)
        p = ax.plot(xh[1:], yh*1./yh.max()+i*3.5, marker='None', linestyle='steps-', color=colors[field], alpha=0.5, linewidth=3)
        #
        pointings = np.unique(sky.pointing[fields == field])
        for pointing in pointings:
            this = sky.pointing == pointing
            p = ax.plot(sky.skyscale[this], np.arange(4)/4.*2+1.5+i*3.5, marker='o', linestyle='-', alpha=0.3, color=colors[field], ms=4)
        #
        ax.text(0.1, i*3.5+1.5, field.replace('SOUTH','S'), horizontalalignment='left')
        
    ax.set_xlim(0,3.4)
    ax.yaxis.set_major_locator(MultipleLocator(3.5))
    ax.set_yticklabels([])
    ax.set_ylabel('Field')
    ax.set_xlabel(r'Background per exposure [e$^-$/ s]')
    
    unicorn.catalogs.savefig(fig, 'pointing_backgrounds.pdf')

def eqw_as_fn_mag():
    """
    For a given line flux, plot the equivalent width 
    as a function of magnitude to show the equivalent width sensitivity.
    """
    
    lam = np.arange(1.e4,1.8e4,0.1)
    
    l0 = 1.4e4
    dv = 120 # km/s
    line = 1./np.sqrt(2*np.pi*(dv/3.e5*l0)**2)*np.exp(-(lam-l0)**2/2/(dv/3.e5*l0)**2)    
    
    continuum = lam*0.+line.max()*0.1
    
    line_fnu = line*lam**2/3.e18
    continuum_fnu = continuum*lam**2/3.e18
    
    xfilt, yfilt = np.loadtxt(os.getenv('iref')+'/F140W.dat', unpack=True)
    yfilt_int = np.interp(lam, xfilt, yfilt) #/np.trapz(yfilt, xfilt)
    
    ## filter width
    piv = np.sqrt(np.trapz(yfilt*xfilt, xfilt)/np.trapz(yfilt/xfilt, xfilt))
    
    INT, SQRT, LAM, THRU, LN = np.trapz, np.sqrt, xfilt, yfilt, np.log
    BARLAM = INT(THRU * LN(LAM) / LAM, LAM) / INT(THRU / LAM, LAM)
    BANDW = BARLAM * SQRT(INT(THRU * LN(LAM / BARLAM)**2 / LAM, LAM)) / INT(THRU / LAM, LAM)
    barlam = np.trapz(yfilt*np.log(xfilt)/xfilt, xfilt) / np.trapz(yfilt/xfilt, xfilt)
    bandw = barlam*np.sqrt(np.trapz(yfilt*np.log(xfilt/barlam)**2/xfilt, xfilt)) / np.trapz(yfilt/xfilt, xfilt)
    nu = 3.e8/(lam*1.e-10)
    
    bigL = -np.trapz(line_fnu*yfilt_int, nu)
    bigC = -np.trapz(continuum_fnu*yfilt_int, nu)
    bigF = -np.trapz(yfilt_int, nu)
    
    bigW = np.trapz(line/continuum, lam)

    xfilt_125, yfilt_125 = np.loadtxt(os.getenv('iref')+'/F125W.dat', unpack=True)
    yfilt_int_125 = np.interp(lam, xfilt_125, yfilt_125) #/np.trapz(yfilt, xfilt)
        
    bigL_125 = -np.trapz(line_fnu*yfilt_int_125, nu)
    bigC_125 = -np.trapz(continuum_fnu*yfilt_int_125, nu)
    bigF_125 = -np.trapz(yfilt_int_125, nu)
        
    # integrated line flux
    alpha = 5.e-17
    
    ### plot trend
    mag = np.arange(20, 25.5, 0.05)
    fnu = 10**(-0.4*(mag+48.6))
    
    EQW = bigW*bigC / (bigF/alpha*fnu-bigL)
    #plt.plot(mag, EQW)
    # 
    # EQW2 = bigW*bigC*alpha/bigF / (fnu-alpha/bigF*bigL)
    # plt.plot(mag, EQW2)
    # 
    # EQW3 = 3.19e-27 / (fnu - 8.04e-31)
    # plt.plot(mag, EQW3)
    # 
    # EQW4 = 3.19e-27 / (10**(-0.4*mag)*3.63e-20 - 8.04e-31)
    # plt.plot(mag, EQW4)
    # 
    # EQW5 = 8.78e-8*(alpha/5.e-17) / (10**(-0.4*mag)-2.21e-11*(alpha/5.e-17))
    # plt.plot(mag, EQW5)
    # 
    # m0 = 23
    # EQW6 = 8.78e-8*(alpha/5.e-17) / (10**(-0.4*m0)-2.21e-11*(alpha/5.e-17))

    #### above equation reaches a limit when continuum = 0, mag can't be less than 
    #### that of just the integrated line
    line_only = (line+continuum*0)
    line_only_fnu = line_only*lam**2/3.e18
    fnu_filt = np.trapz(line_only_fnu*yfilt_int, nu) / np.trapz(yfilt_int, nu)
    mag_limit = -2.5*np.log10(alpha*fnu_filt)-48.6
    mag_limit2 = -2.5*np.log10(alpha/5.e-17)-2.5*np.log10(5.e-17)-2.5*np.log10(fnu_filt)-48.6
    mag_limit3 = -2.5*np.log10(alpha/5.e-17)+26.64
    
    ### test if test case comes out right
    spec_obs = alpha*(line+continuum)
    spec_obs_fnu = spec_obs*lam**2/3.e18
    fnu_filt = np.trapz(spec_obs_fnu*yfilt_int, nu) / np.trapz(yfilt_int, nu)
    mag_obs = -2.5*np.log10(fnu_filt)-48.6
    eqw_obs = np.trapz(line/continuum, lam)
    #plt.plot([mag_obs,mag_obs], [eqw_obs,eqw_obs], marker='o', ms=15)
    
    #### Make figure
    mag = np.arange(20, 26, 0.1)
    fnu = 10**(-0.4*(mag+48.6))
    
    #fig = unicorn.catalogs.plot_init(xs=3.8, left=0.10, bottom=0.08, use_tex=True, square=True, fontsize=11)
    fig = unicorn.catalogs.plot_init(left=0.11, bottom=0.08, xs=3.8, right=0.09, top=0.01, use_tex=True)
    
    ax = fig.add_subplot(111)
    lst = ['-.','-','--',':']
    
    #### Show Arjen's sample
    elg = catIO.Readfile('elg_mag.txt')
    elg.oiii = (5+15)*1./(5+15+23)*elg.ew_obs_tot
    elg.avmag = -2.5*np.log10(0.5*10**(-0.4*elg.j)+0.5*10**(-0.4*elg.h))
    
    ax.scatter(elg.j, elg.oiii, alpha=0.4, s=10, color='red', label=r'van der Wel et al. 2011 ($J_{125}$)')
    EQW_125 = bigW*bigC_125 / (bigF_125/5.e-17*fnu-bigL_125)
    ax.plot(mag, EQW_125, color='red', alpha=0.5)
    
    for ii, limit in enumerate([1.e-17, 3.e-17, 5.e-17, 1.e-16][::-1]):
        EQW = bigW*bigC / (bigF/limit*fnu-bigL)
        ll = np.log10(limit)
        l0 = np.round(10**(ll-np.floor(ll)))
        l10 = np.floor(ll)
        ax.plot(mag, EQW, label=r'$f_{\lambda,\mathrm{line}} = %0d\times10^{%0d}$' %(l0,l10), linestyle=lst[ii], color='black')
    
    ax.semilogy()
    ax.set_xlabel(r'$m_{140}$')
    ax.set_ylabel(r'Equivalent width (\AA)')
    ax.set_ylim(1,1.e4)
    ax.set_xlim(20,26)
    
    #8.78e-8*(alpha/5.e-17) / (10**(-0.4*mag)-2.21e-11*(alpha/5.e-17))
    ax.text(23,1.7,r'$\mathrm{EQW} = \frac{8.78\times10^{-8}\left(f_\mathrm{line}/5\times10^{-17}\right)}{10^{-0.4\,m_{140}}-2.21\times10^{-11}\left(f_\mathrm{line}/5\times10^{-17}\right)}$', horizontalalignment='center')
    
    ax.legend(prop=matplotlib.font_manager.FontProperties(size=8), loc=2, bbox_to_anchor=(0.02,0.98), frameon=False)
    
    unicorn.catalogs.savefig(fig, 'eqw_as_fn_mag.pdf')
    
    #### compute where poisson error of source counts (1.4 um) is similar to sky background error
    sky, texp = 1.6, 1200.
    var = sky*texp + 20**2
    bg_error = np.sqrt(var)/1200.
    
    sens = np.interp(1.4e4, unicorn.reduce.sens_files['A'].field('WAVELENGTH'), unicorn.reduce.sens_files['A'].field('SENSITIVITY'))*46.5
    mag = np.arange(14,25,0.1)
    fnu = 10**(-0.4*(mag+48.6))
    ctrate = fnu*3.e18/1.4e4**2*sens
    
    re = 3 # pix
    peak = 1./np.sqrt(2*np.pi*re**2)
    
    poisson = np.sqrt(ctrate*peak*texp)/texp
    m_crit = np.interp(bg_error, poisson[::-1], mag[::-1])
    
    #### SIMULATIONS
    stats = catIO.Readfile('all_simspec.dat')
    plt.scatter(stats.mag, stats.ha_eqw, marker='s', alpha=0.1, s=4)
    
def make_star_thumbnails():
    """ 
    Extract thumbnails for isolated stars in COSMOS
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/SURVEY_PAPER')
    
    ######### Make full COSMOS catalog    
    file=unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-F140W_drz.fits'
    ROOT_GRISM = os.path.basename(file).split('_drz.fits')[0]
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = file+'[1]'
    se.options['FILTER']    = 'Y'
    se.options['DETECT_THRESH']    = '1.4'
    se.options['ANALYSIS_THRESH']  = '1.4'
    se.options['MAG_ZEROPOINT'] = '26.46'
    status = se.sextractImage(file+'[0]', mode='direct')
    
    cat = threedhst.sex.mySexCat('COSMOS-F140W_drz.cat')
        
    mag, radius = np.cast[float](cat.MAG_AUTO), np.cast[float](cat.FLUX_RADIUS)
    xpix, ypix = np.cast[float](cat.X_IMAGE), np.cast[float](cat.Y_IMAGE)
    ra, dec = np.cast[float](cat.X_WORLD), np.cast[float](cat.Y_WORLD)
    
    #### Find isolated point sources
    points = (mag > 17) & (mag < 22) & (radius < 2.7)
    # plt.plot(mag, radius, marker='o', linestyle='None', alpha=0.5, color='blue')
    # plt.plot(mag[points], radius[points], marker='o', linestyle='None', alpha=0.8, color='red')
    # plt.ylim(0,20)
    # plt.xlim(14,26)
    
    idx = np.arange(len(points))
    isolated = mag > 1.e10
    
    buff = 3 ## buffer, in arcsec
    dmag = 2.5
    scale = 0.06
    
    for i in idx[points]:
        dr = np.sqrt((xpix[i]-xpix)**2+(ypix[i]-ypix)**2)*scale
        near = (dr > 0) & (dr < buff) & (mag < (mag[i]+dmag))
        if len(near[near]) == 0:
            isolated[i] = True
        else:
            isolated[i] = False
    
    #### Make thumbnails
    img = pyfits.open(unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-F140W_drz.fits')
    img_data = img[1].data
    img_wht = img[2].data
    
    NPIX = int(np.ceil(buff/scale))

    prim = pyfits.PrimaryHDU()
    list_d = [prim]
    list_w = [prim]
    
    head = img[1].header
    head['CRPIX1'], head['CRPIX2'] = NPIX, NPIX
    
    for i in idx[points & isolated]:
        print unicorn.noNewLine+'%d' %(i)
        id = np.int(cat.NUMBER[i])
        xi, yi = int(np.round(xpix[i])), int(np.round(ypix[i]))
        sub_data = img_data[yi-NPIX:yi+NPIX, xi-NPIX: xi+NPIX]
        sub_wht = img_wht[yi-NPIX:yi+NPIX, xi-NPIX: xi+NPIX]
        #
        head['CRVAL1'], head['CRVAL2'] = ra[i], dec[i]
        head.update('MAG',mag[i])
        head.update('RADIUS',radius[i])
        head.update('XCENTER',xpix[i]-xi+NPIX)
        head.update('YCENTER',ypix[i]-yi+NPIX)
        #
        list_d.append(pyfits.ImageHDU(sub_data, header=head))
        list_w.append(pyfits.ImageHDU(sub_wht, header=head))
    
    pyfits.HDUList(list_d).writeto('stars_sci.fits', clobber=True)
    pyfits.HDUList(list_w).writeto('stars_wht.fits', clobber=True)
    
def curve_of_growth():
    """
    Some code to evaluate the curve of growth of the F140W PSF, the 
    optimal aperture taken to be the ratio of the CoG divided by the empty aperture
    sigmas, and the overall depth within some aperture.
    """
    import threedhst
    import unicorn
    
    sci = pyfits.open('stars_sci.fits')
    wht = pyfits.open('stars_wht.fits')
    
    apers = np.arange(1,25,0.5)
    lstep = np.arange(0, np.log10(25), 0.05)
    apers = 10**lstep
        
    NOBJ = len(sci)-1
    count = 0
    average_fluxes = apers*0.
    stack = sci[1].data*0.
    
    for i in range(NOBJ):
        print unicorn.noNewLine+'%d' %(i)
        star = sci[i+1].data
        yy, xx = np.indices(star.shape)
        center = (np.abs(xx-50) < 5) & (np.abs(yy-50) < 5)
        xc = np.sum((star*xx)[center])/np.sum(star[center])
        yc = np.sum((star*yy)[center])/np.sum(star[center])
        #xc, yc = sci[i+1].header['XCENTER'], sci[i+1].header['YCENTER']
        #
        bg = threedhst.utils.biweight(star, both=True)
        bg = threedhst.utils.biweight(star[star < (bg[0]+4*bg[1])], both=True)
        star = star-bg[0]
        stack = stack + star
        #        
        NAP = len(apers)
        fluxes = np.zeros(apers.shape)
        for i in range(NAP):
            #print unicorn.noNewLine+'%.2f' %(apers[i])
            fluxes[i] = unicorn.survey_paper.aper_phot(star, xc, yc, apers[i])
        #
        pp = plt.plot(apers, fluxes/fluxes[20], alpha=0.2, color='blue')
        average_fluxes += fluxes/fluxes[20]
        count = count + 1
    
    stack = stack / count
    stack_fluxes = np.zeros(apers.shape)
    for i in range(NAP):
        print unicorn.noNewLine+'%.2f' %(apers[i])
        stack_fluxes[i] = unicorn.survey_paper.aper_phot(star, xc, yc, apers[i])
    
    plt.xlabel(r'$R_\mathrm{aper}$')
    plt.ylabel(r'$f/f_{10}$')
    plt.text(15,0.4,'$N=%d$' %(count))
    
    plt.savefig('curve_of_growth.pdf')
    plt.close()
    
    #plt.plot(apers, average_fluxes/count, color='black', linewidth=2)
    #plt.plot(apers, stack_fluxes/stack_fluxes[20], color='red', linewidth=2)
    # plt.plot(apers, average_fluxes/count / (stack_fluxes/stack_fluxes[20]))
    
    fp = open('curve_of_growth.dat','w')
    for i in range(len(apers)):
        fp.write('%.2f %.3f\n' %(apers[i], (average_fluxes/count)[i]))
    
    fp.close()
    
    #### optimal color aperture:
    empty = pyfits.open('../../EMPTY_APERTURES/COSMOS-1-F140W_drz_empty.fits')
    files = glob.glob('../../EMPTY_APERTURES/*empty.fits')
    for file in files:
        empty = pyfits.open(file)
        apertures = empty[1].data
        fluxes = empty[2].data
        #
        sigma = apertures*0.
        means = apertures*0.
        for iap, ap in enumerate(apertures):
            sigma[iap] = threedhst.utils.biweight(fluxes[:,iap])
            means[iap] = threedhst.utils.biweight(fluxes[:,iap], mean=True)
        #
        ycog = np.interp(apertures, apers, (average_fluxes/count))
        pp = plt.plot(apertures, ycog/(sigma/np.interp(6, apertures, sigma)), alpha=0.5)
    
    plt.xlabel(r'$R_\mathrm{aper}$')
    plt.ylabel(r'CoG / $\sigma$')
    
    plt.plot(apertures, ycog/ycog.max(), color='black', linewidth=2)
    plt.plot(apertures, (sigma/np.interp(6, apertures, sigma))*0.3, color='black', alpha=0.4, linewidth=2)
    plt.savefig('optimal_aperture.pdf', dpi=100)
    plt.close()
    
    #### Actual calculation of the depth    
    from scipy import polyfit, polyval
    
    APER = 0.5 # arcsec, diameter
    
    ycog = average_fluxes/count
    ycog = ycog / ycog.max()
    
    print 'Aperture, D=%.2f"' %(APER)
    
    files = glob.glob('../../EMPTY_APERTURES/[CG]*empty.fits')
    for file in files:
        empty = pyfits.open(file)
        apertures = empty[1].data
        fluxes = empty[2].data
        #
        sigma = apertures*0.
        means = apertures*0.
        for iap, ap in enumerate(apertures):
            sigma[iap] = threedhst.utils.biweight(fluxes[:,iap])
            means[iap] = threedhst.utils.biweight(fluxes[:,iap], mean=True)
        #
        #plt.plot(apertures, 26.46-2.5*np.log10(5*sigma), marker='o', color='black', linestyle='None')
        coeffs = polyfit(np.log10(apertures), np.log10(sigma), 1)
        yfit = 10**polyval(coeffs, np.log10(apertures))
        pp = plt.plot(apertures, 26.46-2.5*np.log10(5*yfit), color='red')
        #
        apcorr = -2.5*np.log10(np.interp(APER/0.06/2, apers, ycog))
        #
        sig_at_aper = 10**polyval(coeffs, np.log10(APER/0.06/2))
        depth = 26.46-2.5*np.log10(5*sig_at_aper)-apcorr
        print '%s - %.2f' %(os.path.basename(file).split('-F14')[0], depth)
        
    plt.ylim(23,30)
    
def aper_phot(array, xc, yc, aper_radius):
    """
    Aperture photometry on an array
    """
    from shapely.geometry import Point, Polygon
    
    point = Point(xc, yc)
    buff = point.buffer(aper_radius, resolution=16)
        
    #### Make the aperture
    im_aper = array*0.
    
    yy, xx = np.indices(array.shape)
    dr = np.sqrt((xx-xc)**2+(yy-yc)**2)
    
    #### these are obviously in the aperture
    solid = dr < (aper_radius-1.5)
    im_aper[solid] = 1.
    
    #### This is the edge
    edge = (dr <= (aper_radius+1.5)) & (dr >= (aper_radius-1.5))
    # for i in range(int(np.floor(xc-aper_radius)),int(np.ceil(xc+aper_radius))):
    #     for j in range(int(np.floor(yc-aper_radius)),int(np.ceil(yc+aper_radius))):
    for i, j in zip(xx[edge], yy[edge]):
        pix = Polygon(((i+0.5,j+0.5), (i+1.5,j+0.5), (i+1.5,j+1.5), (i+0.5,j+1.5)))
        isect = pix.intersection(buff)
        im_aper[j,i] = isect.area
    
    return np.sum(array*im_aper)
    
def make_examples():
    import unicorn
    unicorn.survey_paper.redshift_fit_example(id='GOODS-N-33-G141_00946')
    unicorn.survey_paper.redshift_fit_example(id='GOODS-N-17-G141_00573')
    unicorn.survey_paper.redshift_fit_example(id='GOODS-N-33-G141_01028')
    unicorn.survey_paper.redshift_fit_example(id='COSMOS-1-G141_00252')
    unicorn.survey_paper.redshift_fit_example(id='AEGIS-4-G141_00266')
    unicorn.survey_paper.redshift_fit_example(id='COSMOS-5-G141_00751')

    unicorn.survey_paper.redshift_fit_example(id='PRIMO-1101-G141_01022')
    unicorn.survey_paper.redshift_fit_example(id='GOODS-S-24-G141_00029')
    
    #### Examples
    unicorn.survey_paper.redshift_fit_example(id='COSMOS-14-G141_00100')
    unicorn.survey_paper.redshift_fit_example(id='COSMOS-18-G141_00485')
    
    
    import unicorn
    import unicorn.catalogs
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    mat = lines.id == 'COSMOS-14-G141_00100'
    print lines.id[mat][0], lines.halpha_eqw[mat][0], lines.halpha_eqw_err[mat][0], lines.halpha_flux[mat][0]

    mat = lines.id == 'COSMOS-18-G141_00485'
    print lines.id[mat][0], lines.halpha_eqw[mat][0], lines.halpha_eqw_err[mat][0], lines.halpha_flux[mat][0]
    
def redshift_fit_example(id='COSMOS-18-G141_00485', force=False):
    """ 
    Make a big plot showing how the whole redshift fitting works
    """    
    #id = 'COSMOS-14-G141_00100'  ### aligned along dispersion axis, weak line
    #id = 'GOODS-N-33-G141_00946' ### classic spiral
    #id = 'GOODS-N-17-G141_00573'
    #id = 'COSMOS-18-G141_00485' ### asymmetric line, spiral
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER/EXAMPLE_FITS')
    
    #### Get the necessary files from unicorn
    if (not os.path.exists('%s_thumb.fits.gz' %(id))) | force:
        os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.6/images/%s_thumb.fits.gz' %(id))
        os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.6/images/%s_2D.fits.gz' %(id))
        os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.6/ascii/%s.dat' %(id))
        os.system('rsync -avz --progress $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/OUTPUT/%s* OUTPUT/ ' %(id))
        os.system('rsync -avz --progress $UNICORN:/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS/%s* ./ ' %(id))
    
    zo = threedhst.catIO.Readfile('OUTPUT/%s.zout' %(id))
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    fig = plt.figure(figsize=(6,6))
    dsep = 0.05
    xsep = 0.6
    left = 0.085
    bottom = 0.07
    spec_color = 'purple'
    dy2d = 0.13
    #spec_color = 'blue'
    spec_color = (8/255.,47/255.,101/255.)
    spec_color = 'red'
    phot_color = 'orange'
    #phot_color = (78/255.,97/255.,131/255.)
    #phot_color = '0.7'
    
    spec_color = 'black'
    phot_color = '0.7'
    temp_color = (8/255.,47/255.,101/255.)
    temp_color = 'red'
    
    ########### Full spectrum
    
    ax = fig.add_axes((left, 0.5+bottom+dy2d, 0.99-left-(1-xsep), 0.49-bottom-dy2d))

    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s' %(id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=id, OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    obs_convert = 10**(-0.4*(25+48.6))*3.e18/lci**2/10.**-19*(lci/5500.)**2
    temp_convert = 10**(-0.4*(25+48.6))*3.e18/lambdaz**2/10.**-19*(lambdaz/5500.)**2
    fobs, efobs, obs_sed, temp_sed = fobs*obs_convert, efobs*obs_convert, obs_sed*obs_convert, temp_sed*temp_convert
    
    ymax = max(fobs[is_spec & (fobs > 0)])
        
    ax.semilogx([1],[1])
    
    ## photometry
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='black', linestyle='None', markersize=6, alpha=0.2)
    ## best-fit SED
    ## Spectrum + convolved fit
    #ax.plot(lci[is_spec], obs_sed[is_spec], color='black', markersize=6, alpha=0.7, linewidth=1)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.8, color=spec_color, linewidth=2)
    ax.plot(lambdaz, temp_sed, color='white', linewidth=3, alpha=0.6)
    ax.plot(lambdaz, temp_sed, color=temp_color, alpha=0.6)
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.6, color=phot_color, markersize=10)
        
    ax.set_yticklabels([])
    ax.set_ylabel(r'$f_\lambda$')
    ax.set_xlabel(r'$\lambda$')
    xtick = ax.set_xticks(np.array([0.5, 1., 2, 4])*1.e4)
    ax.set_xticklabels(np.array([0.5, 1., 2, 4]))
    #ax.set_xlim(3000,9.e4)
    ax.set_xlim(3290,2.5e4)
    ax.set_ylim(-0.1*ymax, 1.2*ymax)
    
    ############# Sub spectrum
    ax = fig.add_axes((left, bottom, 0.99-left, 0.49-bottom))
    
    obs_sed_continuum = np.dot(tempfilt['tempfilt'][:,0:7,coeffs['izbest'][0]],coeffs['coeffs'][0:7,0])/(lci/5500.)**2*obs_convert
    
    temp_sed_continuum = np.dot(temp_seds['temp_seds'][:,0:7],coeffs['coeffs'][0:7,0])/(1+zo.z_peak[0])**2*temp_convert
    
    ymax = max(fobs[is_spec & (fobs > 0)]-obs_sed_continuum[is_spec & (fobs > 0)])
    #ymin = min(fobs[is_spec & (fobs > 0)])
    
    # ax.semilogx([1],[1])
    
    ## photometry
    ax.plot(lci[~is_spec], obs_sed[~is_spec]-obs_sed_continuum[~is_spec], marker='o', color='black', linestyle='None', markersize=6, alpha=0.2, zorder=10)
    ## best-fit SED
    ax.plot(lci[is_spec], fobs[is_spec]-obs_sed_continuum[is_spec], marker='None', alpha=0.8, color=spec_color, linewidth=2, zorder=10)
    ax.plot(lambdaz, temp_sed-temp_sed_continuum, color=temp_color, alpha=0.3, zorder=10)
    ## Spectrum + convolved fit
    ax.plot(lci[is_spec], obs_sed[is_spec]-obs_sed_continuum[is_spec], color='white', markersize=6, alpha=0.7, linewidth=4, zorder=10)
    ax.plot(lci[is_spec], obs_sed[is_spec]-obs_sed_continuum[is_spec], color=temp_color, markersize=6, alpha=0.7, linewidth=1, zorder=10)
    #ax.plot(lci[is_spec], obs_sed_continuum[is_spec]-obs_sed_continuum[is_spec], color='black', markersize=6, alpha=0.3, linewidth=2)
    
    ax.errorbar(lci[~is_spec], fobs[~is_spec]-obs_sed_continuum[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.6, color=phot_color, markersize=10)
    
    #ax.set_yticklabels([])
    #ax.set_ylabel(r'$f_\lambda-\ \mathrm{continuum}$')
    ax.set_ylabel(r'$f_\lambda - f_{\lambda,\ \mathrm{cont.}}\ [10^{-19}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}]$')
    ax.set_xlabel(r'$\lambda\ [\mu\mathrm{m}]$')
    xtick = ax.set_xticks(np.array([1.2, 1.4,1.6])*1.e4)
    ax.set_xticklabels(np.array([1.2, 1.4,1.6]))
    #ax.set_xlim(3000,9.e4)
    ax.set_xlim(1.05e4,1.7e4)
    ax.set_ylim(-0.2*ymax, 1.2*ymax)
        
    ########### p(z)
    ax = fig.add_axes((xsep+left, 0.5+bottom+dy2d, 0.99-left-xsep, 0.49-bottom-dy2d))
    
    colors = [spec_color,phot_color,'blue']
    alpha = [0.5, 0.5, 0.2]
    zmin = 4
    zmax = 0
    ymax = 0
    for i in range(2):
        zgrid, pz = eazy.getEazyPz(i, MAIN_OUTPUT_FILE='%s' %(id), 
                          OUTPUT_DIRECTORY='./OUTPUT', 
                          CACHE_FILE='Same')
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        #
        if pz.max() > ymax:
            ymax = pz.max()
        #
        if zgrid[pz > 1.e-3].min() < zmin:
            zmin = zgrid[pz > 1.e-2].min()
        #
        if zgrid[pz > 1.e-6].max() > zmax:
            zmax = zgrid[pz > 1.e-2].max()
    
    ax.plot(zo.z_spec[0]*np.array([1,1]),[0,1.e4], color='green', linewidth=1)
    
    ax.set_yticklabels([])
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$p(z)$')
    ax.xaxis.set_major_locator(unicorn.analysis.MyLocator(4, prune='both'))
    
    ### Plot labels
    #ax.text(0.5, 0.9, '%s' %(id), transform = ax.transAxes, horizontalalignment='center')
    xtxt, align = 0.95,'right'
    xtxt, align = 0.5,'right'
    
    fs, dyt = 9, 0.1
    fs, dyt = 10,0.13
    
    ax.text(xtxt, 0.8, r'$z_\mathrm{phot}=$'+'%5.3f' %(zo.z_peak[1]), transform = ax.transAxes, horizontalalignment=align, fontsize=fs)
    ax.text(xtxt, 0.8-dyt, r'$z_\mathrm{gris}=$'+'%5.3f' %(zo.z_peak[0]), transform = ax.transAxes, horizontalalignment=align, fontsize=fs)
    if zo.z_spec[0] > 0:
        ax.text(xtxt, 0.8-2*dyt, r'$z_\mathrm{spec}=$'+'%5.3f' %(zo.z_spec[0]), transform = ax.transAxes, horizontalalignment=align, fontsize=fs)
        
    ax.set_xlim(zmin, zmax)
    #ax.set_xlim(zgrid.min(), zgrid.max())
    ax.set_ylim(0,1.1*ymax)
    
    #################### 2D spectrum
    thumb = pyfits.open('%s_thumb.fits.gz' %(id))
    thumb_data = thumb[0].data
    #thumb_data[10,:] = 1000
    profile = np.sum(thumb_data, axis=1)
        
    NSUB = int(np.round(0.5*thumb_data.shape[0]))/2
    yc = thumb_data.shape[0]/2
     
    dx = NSUB*2*22/(ax.get_xlim()[1]-ax.get_xlim()[0])*(0.98-left)
    dx = dy2d

    ax = fig.add_axes((left, 0.49, 0.99-left, dy2d))
    #ax.errorbar(lci[~is_spec], fobs[~is_spec]-obs_sed_continuum[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.6, color=phot_color, markersize=10)
    twod_file = '%s_2D.fits.gz' %(id)
    twod = pyfits.open(twod_file)
    spec2d = twod[1].data-twod[4].data
    
    head = twod[1].header
    lam_idx = np.arange(head['NAXIS1'])
    lam = (lam_idx+1-head['CRPIX1'])*head['CDELT1']+head['CRVAL1']
    lam_mima = np.cast[int](np.round(np.interp(np.array([1.05e4,1.7e4]), lam, lam_idx)))
    tick_int = np.interp(np.array([1.2,1.4,1.6])*1.e4, lam, lam_idx)-np.interp(1.05e4, lam, lam_idx)
    
    spec2d_sub = spec2d[yc-NSUB:yc+NSUB,lam_mima[0]:lam_mima[1]]
    ax.imshow(0-spec2d_sub, aspect='auto', vmin=-0.1, vmax=0.01, interpolation='nearest')
    ax.set_yticklabels([]); ax.set_xticklabels([])
    xtick = ax.set_xticks(tick_int); ytick = ax.set_yticks([0,2*NSUB])
    
    ########### Thumbnail
    
    #ax = fig.add_axes((left+left*0.3, 0.49-dx-left*0.3, dx, dx))
    ax = fig.add_axes((left, 0.49, dx, dx))
    #ax.imshow(thumb_data[yc-NSUB:yc+NSUB, yc-NSUB:yc+NSUB], vmin=-0.05, vmax=0.5, interpolation='nearest')
    ax.imshow(0-thumb_data[yc-NSUB:yc+NSUB, yc-NSUB:yc+NSUB], vmin=-0.7, vmax=0.05, interpolation='nearest', zorder=2)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    #ax = fig.add_axes((left+left*0.3*2+dx, 0.49-dx-left*0.3, dx, dx))

    #profile = np.sum(thumb_data[yc-NSUB:yc+NSUB, yc-NSUB:yc+NSUB], axis=0)
    #ax.plot(profile/profile.max(), color='black', alpha=0.4)

    size = thumb[0].data.shape
    twod_file = '%s_2D.fits.gz' %(id)
    twod = pyfits.open(twod_file)
    model1D = np.matrix(twod[5].data.sum(axis=1))
    model1D /= np.max(model1D)
    model2D = np.array(np.dot(np.transpose(model1D),np.ones((1,size[0]))))
    thumb_data *= model2D
    
    profile = np.sum(thumb_data[yc-NSUB:yc+NSUB, yc-NSUB:yc+NSUB], axis=0)
    ax.plot(profile/profile.max()*2*NSUB*0.8, color='black', alpha=0.3, zorder=2)
    ax.set_xlim(0,2*NSUB); ax.set_ylim(0,2*NSUB)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    xtick = ax.set_xticks([0,2*NSUB]); ytick = ax.set_yticks([0,2*NSUB])
    
    fig.savefig('%s_example.pdf' %(id))
    
def equivalent_width_errors():
    """
    Compute the limiting equivalent widths as a function of magnitude or mass
    or something
    """
    import unicorn
    import unicorn.catalogs
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/SURVEY_PAPER')
    
    keep = unicorn.catalogs.run_selection(zmin=0.8, zmax=5.5, fcontam=0.2, qzmin=0., qzmax=0.1, dr=1.0, has_zspec=False, fcovermin=0.9, fcovermax=1.0, massmin=8.5, massmax=15, magmin=17, magmax=23.5)
    keep_22 = unicorn.catalogs.run_selection(zmin=0.8, zmax=5.5, fcontam=0.2, qzmin=0., qzmax=0.1, dr=1.0, has_zspec=False, fcovermin=0.9, fcovermax=1.0, massmin=8.5, massmax=15, magmin=21.7, magmax=22.3)
    
    halpha_sn = lines.halpha_eqw / lines.halpha_eqw_err
    #halpha_sn[(halpha_sn > 0) & (halpha_sn < 1)] = 2
    keep_ha = keep & (halpha_sn[lines.idx] > 0)
    
    oiii_sn = lines.oiii_eqw / lines.oiii_eqw_err
    keep_oiii = keep & (oiii_sn[lines.idx] > 0)
    
    #plt.scatter(phot.mag_f1392w[phot.idx][keep], phot.flux_radius[phot.idx][keep], marker='o')
    marker_size = phot.flux_radius[phot.idx]**1.5
    colors = 'purple'
    # colors = (phot.mag_f1392w[phot.idx]-17)
    # colors[colors < 0] = 0
    # colors[colors > 5] = 5
    
    ##### FLux
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=5/4., left=0.12)
    fig.subplots_adjust(wspace=0.20,hspace=0.24,left=0.12,
                        bottom=0.08,right=0.98,top=0.98)
    
    plt.rcParams['patch.edgecolor'] = 'k'
    ax = fig.add_subplot(211)
        
    ax.scatter(lines.halpha_flux[lines.idx][keep_ha], halpha_sn[lines.idx][keep_ha], marker='o', c='purple', alpha=0.1, s=marker_size[keep_ha])
    ax.scatter(lines.oiii_flux[lines.idx][keep_oiii], oiii_sn[lines.idx][keep_oiii], marker='o', c='orange', alpha=0.1, s=marker_size[keep_oiii])
    
    xm, ym, ys, ns = threedhst.utils.runmed(lines.halpha_flux[lines.idx][keep_ha], halpha_sn[lines.idx][keep_ha], NBIN=20, median=True)
    ax.plot(xm, ym, color='white', alpha=0.6, linewidth=4)
    ax.plot(xm, ym, color='purple', alpha=0.8, linewidth=3)
    xm, ym, ys, ns = threedhst.utils.runmed(lines.oiii_flux[lines.idx][keep_oiii], oiii_sn[lines.idx][keep_oiii], NBIN=20, median=True)
    ax.plot(xm[:-1], ym[:-1], color='white', alpha=0.6, linewidth=4)
    ax.plot(xm[:-1], ym[:-1], color='orange', alpha=0.8, linewidth=3)
    
    ## label
    for si in [2,4,8,16]:
        ax.scatter(np.array([1,1])*2.e-17, np.array([1,1])*25*si**0.4, s=si**1.5, color='black', alpha=0.2)
        ax.text(2.e-17*1.3, 25*si**0.4, '%.1f' %(si*0.06), verticalalignment='center')
    
    nha = len(lines.halpha_flux[lines.idx][keep_ha])
    noiii = len(lines.halpha_flux[lines.idx][keep_oiii])
    ax.text(2.e-17*1.15, 25*(0.5)**0.4, r'$N_\mathrm{H\alpha}=%d$' %(nha), color='purple', horizontalalignment='center')
    ax.text(2.e-17*1.15, 25*(0.5/3)**0.4, r'$N_\mathrm{O III}=%d$' %(noiii), color='orange', horizontalalignment='center')
    
    
    ax.semilogy()
    ax.semilogx()
    ax.set_ylim(1,100)
    ax.set_xlim(1.e-17,2.e-15)
    ax.set_yticklabels([1,3,10,30,100])
    ytick = ax.set_yticks([1,3,10,30,100])
    ax.set_ylabel('line S / N')
    ax.set_xlabel(r'line flux $\mathrm{[ergs\ s^{-1}\ cm^{-2}]}$')
        
    #### EqW
    #fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=1, left=0.12)
    #plt.rcParams['patch.edgecolor'] = 'k'
    ax = fig.add_subplot(212)
    
    marker_size = 10**(-0.4*(18-phot.mag_f1392w[phot.idx]))**0.8
    
    zz = lines.z_grism[lines.idx]*0
    zz = lines.z_grism[lines.idx]
    
    ax.scatter(lines.halpha_eqw[lines.idx][keep_ha]*(1+zz[keep_ha]), halpha_sn[lines.idx][keep_ha], marker='o', c='purple', alpha=0.1, s=marker_size[keep_ha])
    ax.scatter(lines.oiii_eqw[lines.idx][keep_oiii]*(1+zz[keep_oiii]), oiii_sn[lines.idx][keep_oiii], marker='o', c='orange', alpha=0.1, s=marker_size[keep_oiii])
    
    xm, ym, ys, ns = threedhst.utils.runmed(lines.halpha_eqw[lines.idx][keep_ha]*(1+zz[keep_ha]), halpha_sn[lines.idx][keep_ha], NBIN=20, median=False)
    ax.plot(xm, ym, color='white', alpha=0.6, linewidth=4)
    ax.plot(xm, ym, color='purple', alpha=0.8, linewidth=3)    
    xm, ym, ys, ns = threedhst.utils.runmed(lines.oiii_eqw[lines.idx][keep_oiii]*(1+zz[keep_oiii]), oiii_sn[lines.idx][keep_oiii], NBIN=20, median=True)
    ax.plot(xm, ym, color='white', alpha=0.6, linewidth=4)
    ax.plot(xm, ym, color='orange', alpha=0.8, linewidth=3)
    
    for si, mag in enumerate([19, 21, 23]):
        ax.scatter(np.array([1,1])*10, np.array([1,1])*25*(2**(si+1))**0.4, s=10**(-0.4*(18-mag))**0.8, color='black', alpha=0.2)
        ax.text(10*1.3, 25*(2**(si+1))**0.4, '%d' %(mag), verticalalignment='center')
    
    ax.semilogy()
    ax.semilogx()
    ax.set_ylim(1,100)
    ax.set_xlim(5,1000)
    ax.set_yticklabels([1,3,10,30,100])
    ytick = ax.set_yticks([1,3,10,30,100])
    ax.set_xticklabels([5,10,100,500])
    xtick = ax.set_xticks([5,10,100, 500])
    ax.set_ylabel('line S / N')
    if plt.rcParams['text.usetex']:
        ax.set_xlabel(r'Equivalent width [\AA]')
    else:
        ax.set_xlabel(r'Equivalent width [$\AA$]')
        
    fig.savefig('equivalent_width_errors.pdf')
    plt.rcParams['text.usetex'] = False

def zphot_zspec_plot():
    
    import unicorn
    import unicorn.catalogs
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    if unicorn.catalogs.zsp is None:
        unicorn.catalogs.make_specz_catalog()
    
    zsp = unicorn.catalogs.zsp
    
    
    USE_NEW_FITS=True
    if USE_NEW_FITS:
        ##### Refit redshifts gets rid of the offset
        zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_fixed_noTilt.cat')
        #zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_origTemp_noTilt.cat')
        zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_scaleSpecErr3_noTilt.cat')
        zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_scaleSpecErr2_noTilt.cat')
        zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_scaleSpecErr2_yesTilt.cat')
        refit = zout.id[0::3] == 'x'
        refit_idx = zout.z_peak[0::3]*0.
        for i in range(len(zout.id[0::3])):
            print unicorn.noNewLine+'%d' %(i)
            if zout.id[i*3] in zout_new.id:
                refit[i] = True
                refit_idx[i] = np.where(zout_new.id[0::3] == zout.id[i*3])[0][0]
        refit_idx = np.cast[int](refit_idx)

        zphot = zout_new.z_peak[0::3][refit_idx]
        qz = zout_new.q_z[0::3][refit_idx]
        qz2 = zout_new.q_z[2::3][refit_idx]
    else:
        zphot = zout.z_peak[0::3]
        qz = zout.q_z[0::3]
        qz2 = zout.q_z[2::3]
        
    maglim = 24
    qzmax = 0.2
    contam_max = 0.05
    stats_zmin = 0.7
    
    keep = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < contam_max) & (qz < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)
    
    #### Same selection but nothing on specz
    keep_nospec = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < 0.05) & (qz < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zout.z_peak[0::3] > stats_zmin)
    
    keep_nospec_goods = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < 0.05) & (qz < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zout.z_peak[0::3] > stats_zmin) & ((phot.field[phot.idx] == 'GOODS-N') | (phot.field[phot.idx] == 'GOODS-X'))
    
    keep_hasspec = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < 0.05) & (qz < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zout.z_peak[0::3] > stats_zmin) & (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)

    keep_hasspec_goods = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < 0.05) & (qz < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zout.z_peak[0::3] > stats_zmin) & (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1) & ((phot.field[phot.idx] == 'GOODS-N') | (phot.field[phot.idx] == 'GOODS-X'))
    
    #### Spectroscopic redshift ratio by field
    for field in ['GOODS-N', 'GOODS-S', 'COSMOS', 'AEGIS']:
        print '%s %.2f' %(field, len(keep[keep_hasspec & (phot.field[phot.idx] == field)])*1. / len(keep[keep_nospec & (phot.field[phot.idx] == field)]))
        
    
    print len(keep[keep_hasspec])*1./len(keep[keep_nospec]), len(keep[keep_nospec])
    
    #### Only way to get out a few objects where the photometry wasn't found for the fit
    keep = keep & (qz != qz2)
    if USE_NEW_FITS:
        keep = keep & (refit_idx > 0)
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
            
    fig = unicorn.catalogs.plot_init(left=0.07, xs=3, bottom=0.07)
    ax = fig.add_subplot(111)
    
    zsplit = stats_zmin
    ms=2
    
    ax.plot(np.log10(1+zsp.zspec[zsp.mat_idx][keep & (zsp.zspec[zsp.mat_idx] > zsplit)]), np.log10(1+zphot[keep & (zsp.zspec[zsp.mat_idx] > zsplit)]), marker='o', linestyle='None', alpha=0.2, color='black', markersize=ms)
    ax.plot(np.log10(1+zsp.zspec[zsp.mat_idx][keep & (zsp.zspec[zsp.mat_idx] < zsplit)]), np.log10(1+zphot[keep & (zsp.zspec[zsp.mat_idx] < zsplit)]), marker='o', linestyle='None', alpha=0.2, color='0.9', markersize=ms)
    
    ax.plot([0,5],[0,5], color='white', alpha=0.2, linewidth=3)
    #ax.plot([0,5],[0,5], color='black', alpha=0.3, linewidth=1)
    
    zz = np.array([0,4])
    ax.plot(np.log10(1+zz), np.log10(1+zz+0.1*(1+zz)), linestyle='--', color='0.8', alpha=0.5)
    ax.plot(np.log10(1+zz), np.log10(1+zz-0.1*(1+zz)), linestyle='--', color='0.8', alpha=0.5)
    
    ax.set_xticklabels(['0','1','2','3','4'])
    ax.set_xticks(np.log10(1+np.array([0,1,2,3,4])))
    ax.set_yticklabels(['0','1','2','3','4'])
    ax.set_yticks(np.log10(1+np.array([0,1,2,3,4])))
    ax.set_xlim(np.log10(1+0),np.log10(4+1))
    ax.set_ylim(np.log10(1+0),np.log10(4+1))
    
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{G141+phot}$')
    
    dz = (zphot - zsp.zspec[zsp.mat_idx])/(1+zsp.zspec[zsp.mat_idx])
    clip = np.abs(dz) < 0.1
    sigma_gt1 = threedhst.utils.nmad(dz[keep & (zout.z_spec[0::3] > 1)])
    sigma_gt1_clip = threedhst.utils.nmad(dz[keep & (zout.z_spec[0::3] > 1) & clip])
    sigma_gt0_biw = threedhst.utils.biweight(dz[keep & (zout.z_spec[0::3] > stats_zmin)])
    sigma_gt0 = threedhst.utils.nmad(dz[keep & (zout.z_spec[0::3] > stats_zmin)])
    sigma_gt0_clip = threedhst.utils.nmad(dz[keep & (zout.z_spec[0::3] > stats_zmin) & clip])
    NOUT = len(dz[keep & (zout.z_spec[0::3] > stats_zmin) & ~clip])*1./len(dz[keep & (zout.z_spec[0::3] > stats_zmin)])
    
    fs = 9
    
    print sigma_gt0, sigma_gt0_clip, sigma_gt1, sigma_gt1_clip, NOUT
    ax.text(0.1,0.9,r'$H_{140} <\ %.1f,\ z_\mathrm{spec} >\ %.1f,\ Q_z <\ %.2f$' %(maglim, stats_zmin, qzmax), transform=ax.transAxes, fontsize=fs)
    ax.text(0.1,0.81,r'$N=%d$' %(len(dz[keep & (zout.z_spec[0::3] > stats_zmin)])), transform=ax.transAxes, fontsize=fs)
    ax.text(0.1,0.72,r'$\sigma_\mathrm{NMAD}=%.4f$' %(sigma_gt0), transform=ax.transAxes, fontsize=fs)
    pct = '\%'
    ax.text(0.1,0.63,r'$f_\mathrm{>0.1}=%.1f%s$' %(NOUT*100,pct), transform=ax.transAxes, fontsize=fs)
    
    # zbox = np.log10(1+stats_zmin)
    # ax.fill_between([0,zbox],[0,0],[zbox,zbox], color='red', alpha=0.1)
    ax.set_xlim(np.log10(0.0+1),np.log10(3.5+1))
    ax.set_ylim(np.log10(0.0+1),np.log10(3.5+1))
    
    fig.savefig('zphot_zspec.pdf')
    plt.rcParams['text.usetex'] = False
    
    
    ##### Show line misidentifications
    # zha = np.log10(np.array([1.05e4,1.68e4])/6563.)
    # ax.plot(zha, zha+np.log10(6563./5007), color='green', alpha=0.5)
    # ax.plot(zha, zha+np.log10(6563./3727), color='purple', alpha=0.5)
    # ax.plot(zha, zha+np.log10(6563./4863), color='orange', alpha=0.5)
    # zhb = np.log10(np.array([1.05e4,1.68e4])/4861.)
    # ax.plot(zhb, zhb+np.log10(4861./3727), color='blue', alpha=0.5)
    # zoiii = np.log10(np.array([1.05e4,1.68e4])/3727.)
    # ax.plot(zoiii, zoiii+np.log10(3727./4861), color='blue', alpha=0.5)
    
    # plt.xlim(0,np.log10(1+5))
    # plt.ylim(0,np.log10(1+5))
    
    #### Show dz as a function of parameters
    if 1 == 0:
        """
        Make plots to see how the redshift residuals depend on things like mag, 
        contamination fraction, Qz.
        """
        keep = (phot.mag_f1392w[phot.idx] < 25) & (phot.fcontam[phot.idx] < 1) & (zout.q_z[0::3] < 1000) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zsp.zspec[zsp.mat_idx] > stats_zmin) & (zsp.dr < 1)
        
        keep = keep & (zout.q_z[0::3] != zout.q_z[2::3])
        
        dz = (zphot - zsp.zspec[zsp.mat_idx])/(1+zsp.zspec[zsp.mat_idx])
        
        yr = (-0.5,0.5)
        alpha, ms, color = 0.5, 2,'black'
        
        fig = unicorn.catalogs.plot_init(xs=8,aspect=0.7,left=0.12)
        
        #### Mag
        ax = fig.add_subplot(221)
        ax.plot(phot.mag_f1392w[phot.idx][keep], dz[keep], marker='o', alpha=alpha, linestyle='None', ms=ms, color=color)
        xm, ym, ys, ns = threedhst.utils.runmed(phot.mag_f1392w[phot.idx][keep], dz[keep], NBIN=20)
        ax.plot(xm, ys*10, color='red', linewidth=2)
        
        ax.set_ylim(yr[0], yr[1])
        ax.set_xlim(19,25)
        ax.set_xlabel(r'$H_{140}$')
        
        #### Contam
        ax = fig.add_subplot(222)
        ax.plot(phot.fcontam[phot.idx][keep], dz[keep], marker='o', alpha=alpha, linestyle='None', ms=ms, color=color)
        xm, ym, ys, ns = threedhst.utils.runmed(phot.fcontam[phot.idx][keep], dz[keep], NBIN=20)
        ax.plot(xm, ys*10, color='red', linewidth=2)
        #ax.semilogx()
        ax.set_ylim(yr[0], yr[1])
        ax.set_xlim(0.01,1)
        ax.set_xlabel(r'$f_\mathrm{cont}\ \mathrm{at}\ 1.4\ \mu m$')

        #### Q_z
        ax = fig.add_subplot(223)
        ax.plot(zout.q_z[0::3][keep], dz[keep], marker='o', alpha=alpha, linestyle='None', ms=ms, color=color)
        xm, ym, ys, ns = threedhst.utils.runmed(zout.q_z[0::3][keep], dz[keep], NBIN=20)
        ax.plot(xm, ys*10, color='red', linewidth=2)
        ax.semilogx()
        ax.set_ylim(yr[0], yr[1])
        ax.set_xlim(0.001,10)
        ax.set_xlabel(r'$Q_z$')
        
        
    #### Offset near z=1, appears to be due to the tilted slope of the spectrum being a bit too steep.  If you just use a 0th order offset correction to the spectrum (TILT_ORDER=0), the redshifts for many of these objects become correct
    #keep = keep & (np.array(zsp.source)[zsp.mat_idx] == 'Barger')
    
    if (1==0):
        dzlog = np.log10(1+zout.z_peak[0::3]) - np.log10(1+zsp.zspec[zsp.mat_idx])
    
        bad = (dzlog > 0.027) & (dzlog < 0.047 ) & (np.log10(1+zsp.zspec[zsp.mat_idx]) > 0.2) & keep

        bad = (dzlog > 0.18) & keep
        bad = (np.abs(dz) > 0.1) & (zsp.zspec[zsp.mat_idx] > stats_zmin) & keep

        print np.array(zsp.source)[zsp.mat_idx][bad]
        print phot.id[phot.idx][bad]

        for id in phot.id[phot.idx][bad]:
            os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.6/EAZY/%s_eazy.png' %(id))
        
        fig = unicorn.catalogs.plot_init(left=0.12)
        ax = fig.add_subplot(111)
        ax.plot(np.log10(1+zsp.zspec[zsp.mat_idx][keep]), dzlog[keep], marker='o', linestyle='None', alpha=0.5, color='black', markersize=5)
        ax.set_xlim(0,np.log10(4+1))
        ax.set_ylim(-1,1)

        #### nearby offset at z~1 ~ 0.035 in log(1+z)
        offset = 0.036
        print 6563*10**(-offset), 5007*10**(-offset), 4861*10**(-offset), 3727*10**(-offset)

def zphot_zspec_lines():
    """ 
    Investigate how the redshfit errors depend on the emission line signal to noise.
    """
    
    import unicorn
    import unicorn.catalogs
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    if unicorn.catalogs.zsp is None:
        unicorn.catalogs.make_specz_catalog()
    
    zsp = unicorn.catalogs.zsp
    
    zphot = zout.z_peak[0::3]
    
    ##### Refit redshifts gets rid of the offset
    zout_new = catIO.Readfile('/research/HST/GRISM/3DHST/UDF/CATALOGS/LINE_TEMPLATES/full_redshift_fixed_centering.cat')
    refit = zout.id[0::3] == 'x'
    refit_idx = zout.z_peak[0::3]*0.
    for i in range(len(zout.id[0::3])):
        print unicorn.noNewLine+'%d' %(i)
        if zout.id[i*3] in zout_new.id:
            refit[i] = True
            refit_idx[i] = np.where(zout_new.id[0::3] == zout.id[i*3])[0][0]
    refit_idx = np.cast[int](refit_idx)
    
    zphot = zout_new.z_peak[0::3][refit_idx]
    
    maglim = 24
    qzmax = 0.2
    contam_max = 0.05
    
    keep = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < contam_max) & (zout.q_z[0::3] < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)
    
    keep = keep & (zout.q_z[0::3] != zout.q_z[2::3])
    
    lmin, lmax = 1.2e4, 1.6e4
    z_ha = (zsp.zspec[zsp.mat_idx] > (lmin/6563.-1)) & (zsp.zspec[zsp.mat_idx] < (lmax/6563.-1))
    z_oiii = (zsp.zspec[zsp.mat_idx] > (lmin/5007.-1)) & (zsp.zspec[zsp.mat_idx] < (lmax/5007.-1))
    
    dz = (zphot - zsp.zspec[zsp.mat_idx])/(1+zsp.zspec[zsp.mat_idx])
    
    halpha_eqw = lines.halpha_eqw*1.
    oiii_eqw = lines.oiii_eqw*1.

    eqw_min = 0.5
    rnd_halpha = np.random.rand(len(halpha_eqw))*2+1
    rnd_oiii = np.random.rand(len(oiii_eqw))*2+1
    
    halpha_eqw[halpha_eqw < eqw_min] = eqw_min*rnd_halpha[halpha_eqw < eqw_min]
    oiii_eqw[oiii_eqw < eqw_min] = eqw_min*rnd_oiii[oiii_eqw < eqw_min]
    
    ha_color, oiii_color = 'black', 'orange'
    
    
    fig = unicorn.catalogs.plot_init(left=0.15, bottom=0.075, xs=4, right=0.01, top=0.01, square=True)
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    
    fig.subplots_adjust(wspace=0.0)
    
    ################   Ha eqw
    # ax = fig.add_subplot(122)
    # unicorn.survey_paper.dz_trend(halpha_eqw[lines.idx][keep & z_ha], dz[keep & z_ha], xrange=[0.8*eqw_min,1000], yrange=[-0.015, 0.02], xlog=True, ax=ax, xlabel=r'EQW H$\alpha$')
    # yticks = [r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$10^{2}$']
    # ax.set_yticklabels([])
    # xtick = ax.set_xticks([1,10,100,1000])
    # ax.set_xticklabels([1,10,100,1000])
    halpha_sn = lines.halpha_eqw / lines.halpha_eqw_err
    sn_min = 0.2
    rnd_halpha = np.random.rand(len(halpha_sn))*2+1
    
    halpha_sn[halpha_sn < sn_min] = sn_min*rnd_halpha[halpha_sn < sn_min]
    
    ax = fig.add_subplot(122)
    unicorn.survey_paper.dz_trend(halpha_sn[lines.idx][keep & z_ha], dz[keep & z_ha], xrange=[0.1,300], yrange=[-0.015, 0.015], xlog=True, ax=ax, xlabel=r'H$\alpha$ S/N')
    yticks = [r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$10^{2}$']
    ax.set_yticklabels([])
    xtick = ax.set_xticks([1,10,100])
    ax.set_xticklabels([1,10,100])
    
    ################   Mag F140W
    ax = fig.add_subplot(121)
    unicorn.survey_paper.dz_trend(phot.mag_f1392w[phot.idx][keep & z_ha], dz[keep & z_ha], xrange=[19,24], yrange=[-0.015, 0.015], ax=ax, xlabel=r'$m_{140}$')

    ax.text(0.08,0.9,r'H$\alpha$, $%.1f < z < %.1f$' %((lmin/6563.-1), (lmax/6563.-1)), color='black', transform=ax.transAxes, fontsize=12)
    ax.text(0.08,0.83,r'$N=%d$' %(len(z_ha[keep & z_ha])), color='black', transform=ax.transAxes, fontsize=12)

    ax.text(0.08,0.12,r'$\sigma_\mathrm{NMAD}=0.0025$', color='black', transform=ax.transAxes, alpha=0.8)
    ax.text(0.08,0.12,r'$\sigma_\mathrm{NMAD}=0.0025$', color='orange', transform=ax.transAxes, alpha=0.8)
    ax.text(0.08,0.07,r'$\sigma_\mathrm{NMAD}=0.0050$', color='red', transform=ax.transAxes)

    # ################   z 
    # ax = fig.add_subplot(131)
    # unicorn.survey_paper.dz_trend(zsp.zspec[zsp.mat_idx][keep & z_ha], dz[keep & z_ha], xrange=[0.7,1.5], yrange=[-0.015, 0.02], ax=ax)
    
    fig.savefig('zphot_zspec_lines.pdf')
    
def dz_trend(xin, yin, xrange=[0.7,1.5], yrange=[-0.015, 0.015], xlabel=r'$z_\mathrm{spec}$', xlog=False, ax=None, ms=3):
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    ax.plot(xin, yin, linestyle='None', marker='o', alpha=0.2, color='white', zorder=200, ms=ms)
    ax.plot(xin, yin, linestyle='None', marker='o', alpha=0.2, color='black', zorder=201, ms=ms)
    
    xm, ym, ys, ns = threedhst.utils.runmed(xin, yin, NBIN=12, use_nmad=True)
    xm_, ym_, ys_, ns_ = threedhst.utils.runmed(xin, yin, NBIN=10, use_nmad=True)
    xm[0], ym[0], ys[0], ns[0] = xrange[0]*0.8, ym_[0], ys_[0], ns_[0]
    xm[1:11], ym[1:11], ys[1:11], ns[1:11] = xm_, ym_, ys_, ns_
    xm[-1], ym[-1], ys[-1], ns[-1] = xrange[1]*1.2, ym_[-1], ys_[-1], ns_[-1]

    ax.plot(xm, ym, color='black', alpha=0.9, zorder=101, marker='o', linewidth=3)
    ax.fill_between(xm, ym+ys, ym-ys, color='black', alpha=0.4, zorder=100)
    yx = ys*0+0.0025
    ax.plot(xm, ym+yx, color='orange', alpha=0.9, zorder=101, linewidth=3)
    ax.plot(xm, ym-yx, color='orange', alpha=0.9, zorder=101, linewidth=3)
    yx = ys*0+0.005
    ax.plot(xm, ym+yx, color='red', alpha=0.9, zorder=101, linewidth=3)
    ax.plot(xm, ym-yx, color='red', alpha=0.9, zorder=101, linewidth=3)

    ax.plot(xm, ym*0, color='white', alpha=0.8, zorder=301, linewidth=3, linestyle='--')
    ax.plot(xm, ym*0, color='black', alpha=0.8, zorder=302, linewidth=3, linestyle='--')
    
    if xlog:
        ax.semilogx()
        
    ax.set_xlim(xrange[0],xrange[1])
    ax.set_ylim(yrange[0],yrange[1])
    ax.set_ylabel(r'$\Delta z / (1+z)$')
    ax.set_xlabel(xlabel)
    
def compute_SFR_limits():
    import cosmocalc as cc
    
    limiting_flux = 4.e-17
    ### z=1, H-alpha
    cosmo = cc.cosmocalc(z=1.0)
    SFR_ha = 7.9e-42 * limiting_flux * 4 * np.pi * cosmo['DL_cm']**2
    
    ### z=2, OII
    cosmo = cc.cosmocalc(z=2.0)
    SFR_oii = 1.4e-41 * limiting_flux * 4 * np.pi * cosmo['DL_cm']**2
    
    print SFR_ha, SFR_oii
    
def number_counts():
    import unicorn
    import unicorn.catalogs
    import copy
        
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
        
    keep = unicorn.catalogs.run_selection(zmin=0, zmax=8, fcontam=1, qzmin=0., qzmax=10, dr=1.0, has_zspec=False, fcovermin=0.5, fcovermax=1.0, massmin=0, massmax=15, magmin=12, magmax=27)
    
    fields = (phot.field[phot.idx] == 'AEGIS') | (phot.field[phot.idx] == 'COSMOS') | (phot.field[phot.idx] == 'GOODS-N') | (phot.field[phot.idx] == 'GOODS-S')
    
    #fields = (phot.field[phot.idx] == 'COSMOS') | (phot.field[phot.idx] == 'AEGIS') | (phot.field[phot.idx] == 'GOODS-S')

    #fields = (phot.field[phot.idx] == 'COSMOS') | (phot.field[phot.idx] == 'AEGIS')
    
    #fields = (phot.field[phot.idx] == 'GOODS-N')
    
    fields = fields & (phot.fcover[phot.idx] > 0.5)
    
    pointings = []
    for field, pointing in zip(phot.field, phot.pointing):
        pointings.append('%s-%d' %(field, pointing))
    pointings = np.array(pointings)
    NPOINT = len(np.unique(pointings[phot.idx][fields]))
    
    xrange = (12,25)
    nbin = np.int(np.round((xrange[1]-xrange[0])*10.))
    binwidth = (xrange[1]-xrange[0])*1./nbin
    
    normal = 1./binwidth/NPOINT
    
    cumul = True
    normal = 1./NPOINT
    
    #normal = 1./NPOINT*148
    
    ##### OFFSET TO TOTAL!
    m140 = phot.mag_f1392w - 0.22
    #m140 = phot.mag_f1392w
    
    #### Full histogram
    y_full, x_full = np.histogram(m140[phot.idx][fields], bins=nbin, range=xrange)
    x_full = (x_full[1:]+x_full[:-1])/2.
    if cumul:
        y_full = np.cumsum(y_full)
    
    y_full, x_full = np.histogram(m140[phot.idx][fields], bins=nbin, range=xrange)
    x_full = (x_full[1:]+x_full[:-1])/2.
    if cumul:
        y_full = np.cumsum(y_full)
    
    lo_full, hi_full = threedhst.utils.gehrels(y_full)
    
    #### Matched in photometric catalogs
    matched = mcat.rmatch[mcat.idx] < 1.
    matched = zout.z_peak[0::3] != zout.z_peak[2::3]
    
    y_matched, x_matched = np.histogram(m140[phot.idx][fields & matched], bins=nbin, range=xrange)
    x_matched = (x_matched[1:]+x_matched[:-1])/2.
    if cumul:
        y_matched = np.cumsum(y_matched)
    
    #### point sources
    xpoint, ypoint = np.array([14,18,23]), np.array([6,3.18, 2.8])
    ypoint_int = np.interp(m140, xpoint, ypoint)
    points = (phot.flux_radius[phot.idx] < ypoint_int[phot.idx]) #& (m140[phot.idx] < 23)
    
    y_points, x_points = np.histogram(m140[phot.idx][fields & matched & points], bins=nbin, range=xrange)
    x_points = (x_points[1:]+x_points[:-1])/2.
    if cumul:
        y_points = np.cumsum(y_points)
    
    #### Low contamination
    contam = phot.fcontam[phot.idx] < 0.1
    y_contam, x_contam = np.histogram(m140[phot.idx][fields & contam & matched], bins=nbin, range=xrange)
    x_contam = (x_contam[1:]+x_contam[:-1])/2.
    if cumul:
        y_contam = np.cumsum(y_contam)
    
    #### z > 1
    z1 = (zout.z_peak[0::3] > 1) & (zout.q_z[0::3] < 50.5) & ~points
    y_z1, x_z1 = np.histogram(m140[phot.idx][fields & matched & z1], bins=nbin, range=xrange)
    x_z1 = (x_z1[1:]+x_z1[:-1])/2.
    if cumul:
        y_z1 = np.cumsum(y_z1)
    
    lo_z1, hi_z1 = threedhst.utils.gehrels(y_z1)
    #wx, wy = np.loadtxt('whitaker_completeness.dat', unpack=True)
    #wscale = np.interp(x_z1, wx, wy)
    wscale = y_matched*1. / y_full
    wscale[~np.isfinite(wscale)] = 1
    wscale[wscale > 1] = 1
    wscale[wscale == 0] = 1
    #hi_z1 /= wscale
    # lo_z1 /= wscale
    # y_z1 /= wscale
    
    #### No cut on Q_z
    z1q = (zout.z_peak[0::3] > 1) & (zout.q_z[0::3] < 100)  & ~points
    y_z1q, x_z1q = np.histogram(m140[phot.idx][fields & matched & z1q], bins=nbin, range=xrange)
    x_z1q = (x_z1q[1:]+x_z1q[:-1])/2.
    if cumul:
        y_z1q = np.cumsum(y_z1q)
    
    #### Total number at z>1
    print 'NPOINT: %d' %(NPOINT)
    
    #z1q_mag = unicorn.catalogs.run_selection(zmin=1, zmax=5.5, fcontam=1, qzmin=0., qzmax=100, dr=1.0, has_zspec=False, fcovermin=0.5, fcovermax=1.0, massmin=0, massmax=15, magmin=0, magmax=23)
    
    z1q_mag = z1q & fields & (m140[phot.idx] <= 23.8)  & ~points
    N_z1_total = len(z1q_mag[z1q_mag])*1./NPOINT*149.
    N_total = len(z1q_mag[matched & fields & (m140[phot.idx] <= 23.8)])*1./NPOINT*149.
    print 'N (z>1, m<23) = %d, N_total = %d' %(N_z1_total, N_total)
    print 'N (z>1, m<23) = %d' %(np.interp(23.8, x_z1, y_z1*149./NPOINT))

    #### z > 2
    z2 = (zout.z_peak[0::3] > 2) & (zout.q_z[0::3] < 50.5)  & ~points
    y_z2, x_z2 = np.histogram(m140[phot.idx][fields & matched & z2], bins=nbin, range=xrange)
    x_z2 = (x_z2[1:]+x_z2[:-1])/2.
    if cumul:
        y_z2 = np.cumsum(y_z2)
    
    lo_z2, hi_z2 = threedhst.utils.gehrels(y_z2)
    #hi_z2 /= wscale
    
    #### Tail of bright objects in the z>2 set
    tail = (zout.z_peak[0::3] > 2) & (zout.q_z[0::3] < 50.5)  & ~points & fields & matched & (m140[phot.idx] < 21)
    print 'z2 tail:', zout.id[0::3][tail], mcat.rmatch[mcat.idx][tail], phot.flux_radius[phot.idx][tail], np.interp(m140[phot.idx][tail], xpoint, ypoint)
    
    #### No cut on Q_z
    z2q = (zout.z_peak[0::3] > 2) & (zout.q_z[0::3] < 100) & ~points
    y_z2q, x_z2q = np.histogram(m140[phot.idx][fields & matched & z2q], bins=nbin, range=xrange)
    x_z2q = (x_z2q[1:]+x_z2q[:-1])/2.
    if cumul:
        y_z2q = np.cumsum(y_z2q)
        
    #### NMBS comparison
    cat_nmbs, zout_nmbs, fout_nmbs = unicorn.analysis.read_catalogs(root='COSMOS-1')
    #nmbs_hmag = 25-2.5*np.log10(cat_nmbs.H1*cat_nmbs.Ktot/cat_nmbs.K)
    #nmbs_hmag = 25-2.5*np.log10((cat_nmbs.H1+cat_nmbs.J3+cat_nmbs.J2+cat_nmbs.H2)/4.*cat_nmbs.Ktot/cat_nmbs.K)
    nmbs_hmag = 25-2.5*np.log10((cat_nmbs.H1+cat_nmbs.J3)/2.*cat_nmbs.Ktot/cat_nmbs.K)
    keep_nmbs = (cat_nmbs.wmin > 0.3)
    
    y_nmbs, x_nmbs = np.histogram(nmbs_hmag[keep_nmbs], bins=nbin, range=xrange)
    x_nmbs = (x_nmbs[1:]+x_nmbs[:-1])/2.
    if cumul:
        y_nmbs = np.cumsum(y_nmbs)
    
    y_nmbs *= 1./(0.21*2*3600.)*4*NPOINT
    
    z1_nmbs = (zout_nmbs.z_peak > 1) & (cat_nmbs.star_flag == 0)
    
    y_nmbs_z1, x_nmbs_z1 = np.histogram(nmbs_hmag[keep_nmbs & z1_nmbs], bins=nbin, range=xrange)
    x_nmbs_z1 = (x_nmbs_z1[1:]+x_nmbs_z1[:-1])/2.
    if cumul:
        y_nmbs_z1 = np.cumsum(y_nmbs_z1)
    
    y_nmbs_z1 *= 1./(0.21*2*3600)*4*NPOINT
    
    z2_nmbs = (zout_nmbs.z_peak > 2) & (cat_nmbs.star_flag == 0)
    
    y_nmbs_z2, x_nmbs_z2 = np.histogram(nmbs_hmag[keep_nmbs & z2_nmbs], bins=nbin, range=xrange)
    x_nmbs_z2 = (x_nmbs_z2[1:]+x_nmbs_z2[:-1])/2.
    if cumul:
        y_nmbs_z2 = np.cumsum(y_nmbs_z2)
    
    y_nmbs_z2 *= 1./(0.21*2*3600)*4*NPOINT
    
    #
    nmbs_stars = (cat_nmbs.star_flag == 1)
    
    y_nmbs_stars, x_nmbs_stars = np.histogram(nmbs_hmag[keep_nmbs & nmbs_stars], bins=nbin, range=xrange)
    x_nmbs_stars = (x_nmbs_stars[1:]+x_nmbs_stars[:-1])/2.
    if cumul:
        y_nmbs_stars = np.cumsum(y_nmbs_stars)
    
    y_nmbs_stars *= 1./(0.21*2*3600)*4*NPOINT
    
    #### Make the plot
    fig = unicorn.catalogs.plot_init(left=0.11, bottom=0.08, xs=3.8, right=0.09, top=0.01)
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times'
    ax = fig.add_subplot(111)
    
    ax.plot(x_full, y_full*normal, color='black')
    ax.fill_between(x_full,lo_full*normal, hi_full*normal, color='black', alpha=0.4)

    ax.plot(x_matched, y_matched*normal, color='blue',alpha=0.8)
    ax.plot(x_contam, y_contam*normal, color='green',alpha=0.8)
    ax.plot(x_points[x_points <= 23], y_points[x_points < 23]*normal, color='purple',alpha=0.8)
    ax.plot(x_points[x_points >= 23], y_points[x_points >= 23]*normal, color='purple',alpha=0.8, linestyle=':')
    
    ax.plot(x_z1, y_z1*normal, color='orange',alpha=0.7)
    ax.fill_between(x_z1,lo_z1*normal, hi_z1*normal, color='orange', alpha=0.4)
    ax.plot(x_z1q, y_z1q*normal, color='orange',alpha=0.7, linestyle='--')
    
    ax.plot(x_z2, y_z2*normal, color='red',alpha=0.7)
    ax.fill_between(x_z2,lo_z2*normal, hi_z2*normal, color='red', alpha=0.4)
    ax.plot(x_z2q, y_z2q*normal, color='red',alpha=0.7, linestyle='--')

    # ax.plot(x_nmbs, y_nmbs*normal, color='black',alpha=0.8, linewidth=3, alpha=0.2)
    # ax.plot(x_nmbs_z1, y_nmbs_z1*normal, color='orange',alpha=0.8, linewidth=3, alpha=0.2)
    # ax.plot(x_nmbs_z2, y_nmbs_z2*normal, color='red',alpha=0.8, linewidth=3, alpha=0.2)
    # ax.plot(x_nmbs_stars, y_nmbs_stars*normal, color='purple',alpha=0.8, linewidth=3, alpha=0.2)
    
    #ax.text(0.05,0.92,r'%s ($N=%d$)' %(', '.join(np.unique(phot.field[phot.idx][fields])), NPOINT), color='black', transform=ax.transAxes)
    ax.text(0.05,0.92,r'Total, from $N=%d$ pointings' %(NPOINT), color='black', transform=ax.transAxes)
    ax.text(0.05,0.87,r'matched', color='blue', transform=ax.transAxes)
    ax.text(0.05,0.82,r'$f_\mathrm{contam} < 10\%$', color='green', transform=ax.transAxes)
    ax.text(0.05,0.77,r'point sources', color='purple', transform=ax.transAxes)
    ax.text(0.05,0.72,r'$z > 1$', color='orange', transform=ax.transAxes)
    ax.text(0.05,0.67,r'$z > 2$', color='red', transform=ax.transAxes)

    ax.set_xlabel('MAG\_AUTO (F140W $\sim$ $H$)')
    if cumul:
        ax.set_ylabel('N($<m$) per WFC3 pointing')
    else:
        ax.set_ylabel('N / pointing / mag')
    ax.semilogy()
    
    yticks = [r'$0.01$',r'$0.1$',r'$1$',r'$10$',r'$10^{2}$']
    ax.set_yticklabels(yticks)
    ytick = ax.set_yticks([0.01,0.1,1,10,100])
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    minorLocator   = MultipleLocator(1)
    ax.xaxis.set_minor_locator(minorLocator)
        
    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(0.01, 500)
    
    ax2 = ax.twinx()
    ax2.semilogy()
    yticks = [r'$10$',r'$10^{2}$',r'$10^{3}$',r'$10^{4}$']
    ax2.set_yticklabels(yticks)
    ytick = ax2.set_yticks([10,100,1000,1.e4])

    ax2.set_ylim(0.01*149, 500*149)

    ax2.set_ylabel('N($<m$), full survey')
    ax2.set_xlim(xrange[0], xrange[1])
    ax2.xaxis.set_minor_locator(minorLocator)
    
    ### Grid
    ax.xaxis.grid(alpha=0.35, zorder=1, which='major')
    ax.xaxis.grid(alpha=0.2, zorder=1, which='minor')
    ax2.yaxis.grid(alpha=0.35, zorder=1, which='major')
    
    fig.savefig('number_counts.pdf')
    plt.rcParams['text.usetex'] = False

def ancillary_matches():
    """
    Get an idea of how the matching to the ancillary catalogs depends on mag:
        matched fraction
        multiple matches
    """
    import unicorn
    import unicorn.catalogs
    import copy
        
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
        
    keep = unicorn.catalogs.run_selection(zmin=0, zmax=8, fcontam=1, qzmin=0., qzmax=10, dr=1.0, has_zspec=False, fcovermin=0.5, fcovermax=1.0, massmin=0, massmax=15, magmin=12, magmax=27)
    
    fields = (phot.field[phot.idx] == 'AEGIS') | (phot.field[phot.idx] == 'COSMOS') | (phot.field[phot.idx] == 'GOODS-N') | (phot.field[phot.idx] == 'GOODS-S')
    
    phot_dr = np.zeros(phot.field.shape)+100
    phot_id = np.zeros(phot.field.shape)
    phot_kmag = np.zeros(phot.field.shape)
    idx = np.arange(phot.field.shape[0])
    
    #### Do separate matching again on every object in photometric catalog
    for field in ['COSMOS','AEGIS','GOODS-S','GOODS-N']:
        this = phot.field == field
        cat, zout, fout = unicorn.analysis.read_catalogs(field+'-1')
        cos_dec = np.cos(np.median(cat.dec)/360.*2*np.pi)**2
        for i in idx[this]:
            print unicorn.noNewLine+'%d / %d' %(i, idx[this][-1])
            dr = np.sqrt((cat.ra-phot.x_world[i])**2*cos_dec+(cat.dec-phot.y_world[i])**2)*3600.
            ma = dr == dr.min()
            phot_dr[i] = dr.min()
            phot_id[i] = cat.id[ma][0]
            phot_kmag[i] = cat.kmag[ma][0]
    
    #### Ask, "what fraction of F140W objects have multiple matches to the same ancillary object"
    n_match = phot_dr*0
    n_brighter = phot_dr*0.
    i = 0
    base_selection = (phot.fcover > 0.5) & (phot.has_spec == 1) & (phot_dr < 1.0) 
    for f,id,m,p in zip(phot.field, phot_id, phot.mag_f1392w,phot.pointing):
        print unicorn.noNewLine+'%d / %d' %(i, len(phot_id))
        mat = (phot.field == f) & (phot_id == id) & (phot.pointing == p) & base_selection
        n_match[i] = mat.sum()
        brighter = mat & (phot.mag_f1392w-m < 0.75)
        n_brighter[i] = brighter.sum()-1
        i = i+1
    
    use = n_match > 0
    yh_full, xh_full = np.histogram(phot.mag_f1392w[use], range=(12,26), bins=14*4)
    
    fig = unicorn.plotting.plot_init(square=True, use_tex=True, left=0.09, bottom=0.07, xs=3.5)
    ax = fig.add_subplot(111)
    
    yh_n, xh_n = np.histogram(phot.mag_f1392w[use & (n_match > 1)], range=(12,26), bins=14*4)
    ax.plot(xh_n[1:]-0.22, yh_n*1./yh_full, color='blue', linestyle='steps', label=r'$N_{\rm match} > 1$')

    yh_n, xh_n = np.histogram(phot.mag_f1392w[use & (n_match > 1) & (n_brighter == 1)], range=(12,26), bins=14*4)
    ax.plot(xh_n[1:]-0.22, yh_n*1./yh_full, color='red', linestyle='steps', label=r'$N_{\rm brighter} = 1$')
    
    yh_n, xh_n = np.histogram(phot.mag_f1392w[use & (n_match > 1) & (n_brighter > 1)], range=(12,26), bins=14*4)
    ax.plot(xh_n[1:]-0.22, yh_n*1./yh_full, color='orange', linestyle='steps', label=r'$N_{\rm brighter} > 1$')
    
    ax.set_xlabel(r'$m_{140}$')
    ax.set_ylabel('fraction')
    
    ax.set_xlim(19,24.5)
    ax.set_ylim(0,0.21)
    ax.legend(loc='upper left', frameon=False)
    unicorn.plotting.savefig(fig, 'ancillary_matched_from_f140w.pdf')
    
    #### Check these cases of n_brighter == 1
    test = (n_match > 0) & (n_brighter == 1)
    idx = np.arange(len(n_match))[test]
    i = 0
    id = phot_id[idx][i]
    mat = base_selection & (phot.field == phot.field[idx][i]) & (phot_id == phot_id[idx][i])
    
    ### Some pointings, such as GOODS-S flanking fields don't overlap with photom. catalog
    test_field_goodsn = (phot.field == 'GOODS-N') 
    test_field_goodss = (phot.field == 'GOODS-S') & (phot.pointing != 1) & (phot.pointing != 28)
    test_field_cosmos = phot.field == 'COSMOS'
    
    test_field_aegis = phot.field == 'AEGIS'  ### out of NMBS
    for i in [11,2,1,6]:
        test_field_aegis = test_field_aegis & (phot.pointing != i)
    
    fig = unicorn.plotting.plot_init(square=True, use_tex=True, left=0.09, bottom=0.07, xs=3.5)
    ax = fig.add_subplot(111)
    
    #### Make a plot showing the fraction of matched galaxies
    for test_field, c in zip([test_field_goodsn, test_field_goodss, test_field_cosmos, test_field_aegis], ['orange','red','blue','green']):
        base_selection = (phot.fcover > 0.5) & test_field & (phot.has_spec == 1)
        has_match = phot_dr < 1.0
        yh_full, xh_full = np.histogram(phot.mag_f1392w[base_selection], range=(12,26), bins=14*4)
        yh_mat, xh_mat = np.histogram(phot.mag_f1392w[base_selection & has_match], range=(12,26), bins=14*4)
        yh_full, yh_mat = np.maximum(yh_full, 0.01), np.maximum(yh_mat, 0.01)
        # plt.plot(xh_full[1:], yh_full, linestyle='steps', color='blue', alpha=0.5)
        # plt.plot(xh_mat[1:], yh_mat, linestyle='steps', color='red', alpha=0.5)
        # plt.semilogy()
        # plt.ylim(0.5,500)
        #
        ax.plot(xh_full[1:]-0.22, yh_mat/yh_full, linestyle='-', linewidth=3, color=c, alpha=0.5, label=np.unique(phot.field[test_field])[0])
    
    ax.legend(loc='lower left')
    ax.plot([0,100],[1,1], color='black', linestyle='-', alpha=0.8, linewidth=2)
    ax.plot([0,100],[0.9,0.9], color='black', linestyle=':', alpha=0.8, linewidth=2)
    ax.set_ylim(0,1.1)
    ax.set_xlim(21,25.)
    ax.set_xlabel(r'$m_{140}$')
    ax.set_ylabel(r'Matched fraction')
    
    unicorn.plotting.savefig(fig, 'ancillary_matched_fraction.pdf')
    
    #### Look at multiple matches
    base_selection = (phot.fcover > 0.5) & (phot.has_spec == 1)
    
    use = base_selection & test_field_cosmos & (phot_dr < 1.0)
    plt.scatter(phot.mag_f1392w[use], phot_kmag[use], color='blue', alpha=0.1, s=10)
    
    matched_id = np.unique(phot_id[use])
    kmag = matched_id*0.
    dmag1 = matched_id*0.+100
    dmag2 = matched_id*0.+100
    N = matched_id*0
    for ii, id in enumerate(matched_id):
        print unicorn.noNewLine+'%d / %d' %(ii, len(matched_id))
        this = (phot_id == id) & use
        dmag = phot.mag_f1392w[this]-phot_kmag[this]
        kmag[ii] = phot_kmag[this][0]
        dmag1[ii] = dmag[0]
        N[ii] = this.sum()
        if this.sum() > 1:
            so = np.argsort(dmag)
            dmag2[ii] = dmag[so][1] 
    #
    fig = unicorn.plotting.plot_init(square=True, use_tex=True, left=0.09, bottom=0.07, xs=3.5)
    ax = fig.add_subplot(111)
    ax.scatter(kmag, dmag1-0.22, color='blue', alpha=0.2, s=10, label='1st matched')
    ax.scatter(kmag, dmag2-0.22, color='red', alpha=0.2, s=10, label='2nd matched')
    ax.set_xlim(17,24)
    ax.set_ylim(-1,5)
    ax.legend(loc='upper left')
    ax.set_xlabel(r'$K_\mathrm{matched}$')
    ax.set_ylabel(r'$m_{140} - K_\mathrm{matched}$')
    unicorn.plotting.savefig(fig,'ancillary_delta_mag.pdf')
    
    ### Show fraction of ancillary objects that have multiple matches a function of magnitude
    fig = unicorn.plotting.plot_init(square=True, use_tex=True, left=0.09, bottom=0.07, xs=3.5)
    ax = fig.add_subplot(111)
    
    yh_full, xh_full = np.histogram(kmag, range=(17,24), bins=7*4)
    
    yh, xh = np.histogram(kmag[dmag2 < 1.2], range=(17,24), bins=7*4)
    ax.plot(xh[1:], yh*1./yh_full, linestyle='steps', color='red', linewidth=3, alpha=0.5, label=r'$\Delta 2^{\rm nd} < 1.2$')
    
    yh, xh = np.histogram(kmag[N > 1], range=(17,24), bins=7*4)
    ax.plot(xh[1:], yh*1./yh_full, linestyle='steps', color='red', linewidth=3, alpha=0.3, label=r'$N_\mathrm{match} > 1$')

    yh, xh = np.histogram(kmag[N > 3], range=(17,24), bins=7*4)
    ax.plot(xh[1:], yh*1./yh_full, linestyle='steps', color='blue', linewidth=3, alpha=0.5, label=r'$N_\mathrm{match} > 3$')
    
    ax.set_xlabel(r'$K_\mathrm{matched}$')
    ax.set_ylabel(r'fraction')
    ax.legend(loc='upper left', prop=matplotlib.font_manager.FontProperties(size=9))
    
    unicorn.plotting.savefig(fig,'ancillary_multiple_fraction.pdf')
    
    
def get_iband_mags():
    """ 
    On Unicorn, loop through the ascii spectra to retrieve the iband mags, should all be ZP=25.
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/')
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    ids = zout.id[0::3]
    fields = phot.field[phot.idx]
    
    iflux = zout.z_peak[0::3]*0.-1
    imod = iflux*1.
    lc_i = iflux*1.
    hflux = iflux*1
    hmod = iflux*1.
    lc_h = iflux*1
    
    count = 0
    for id, field in zip(ids, fields):
        path = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS_v1.6/ASCII/%s/%s_obs_sed.dat' %(field, id)
        if os.path.exists(path):
            print unicorn.noNewLine+id
            obs = catIO.Readfile(path)
            dlam_spec = obs.lc[-1]-obs.lc[-2]
            is_spec = np.append(np.abs(1-np.abs(obs.lc[1:]-obs.lc[0:-1])/dlam_spec) < 0.05,True)
            dl_i = np.abs(obs.lc-7688.1)
            dl_h = np.abs(obs.lc[~is_spec]-1.6315e4)
            ix_i = np.where(dl_i == dl_i.min())[0][0]
            ix_h = np.where(dl_h == dl_h.min())[0][0]
            iflux[count] = obs.fnu[ix_i]
            imod[count] = obs.obs_sed[ix_i]
            lc_i[count] = obs.lc[ix_i]
            hflux[count] = obs.fnu[ix_h]
            hmod[count] = obs.obs_sed[ix_h]
            lc_h[count] = obs.lc[ix_h]
        #    
        count = count+1
    
    fp = open('full_imag_hmag.dat','w')
    fp.write('# id iflux imodel lc_i hflux hmodel lc_h\n')
    for i in range(len(ids)):
        fp.write('%s %.5e %.5e %.1f %.5e %.5e %.1f\n' %(ids[i], iflux[i], imod[i], lc_i[i], hflux[i], hmod[i], lc_h[i]))
    fp.close()
    
def zspec_colors():
    """ 
    Show as a function if H / (i-H) where the galaxies with zspec fall
    """
    import unicorn
    import unicorn.catalogs
    import copy
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    if unicorn.catalogs.zsp is None:
        unicorn.catalogs.make_specz_catalog()
    
    zsp = unicorn.catalogs.zsp
    
    maglim = 25
    qzmax = 200
    contam_max = 0.5
    
    
    ###### Selection criteria
    
    keep = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < contam_max) & (zout.q_z[0::3] < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) #& (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)
    keep = keep & (zout.q_z[0::3] != zout.q_z[2::3])
    
    has_specz = (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)
    
    #mag, radius = np.cast[float](cat.MAG_AUTO), np.cast[float](cat.FLUX_RADIUS)
    #### Find isolated point sources
    points = (phot.flux_radius[phot.idx] < 2.7)
    keep = keep & (~points)
    
    zphot = zout.z_peak[0::3]
    
    ##### H mags from F140W and matches
    icat = catIO.Readfile('full_imag_hmag.dat')
    IH = -2.5*np.log10(icat.iflux / icat.hflux)
    
    phot_zp = zphot*0.+25
    phot_zp[(phot.field[phot.idx] == 'GOODS-S') | (phot.field[phot.idx] == 'PRIMO') | (phot.field[phot.idx] == 'WFC3-ERSII-G01') | (phot.field[phot.idx] == 'GEORGE')] = 23.86
    m140 = phot.mag_f1392w[phot.idx]-0.22  #### Offset to total in catalogs!
    hmag = phot_zp-2.5*np.log10(icat.hflux)
    
    fin = np.isfinite(hmag) & (icat.iflux > 0) & (mcat.rmatch[mcat.idx] < 1)
    #### Few wierd objects with very discrepant H mags in GOODS-N
    bad = (zout.z_peak[0::3] < 1) & (IH > 3.5)
    fin  = fin & (~bad)
    
    ######### Compare h mags    
    # use = fin
    # 
    # use = (phot.field[phot.idx] == 'GOODS-S') | (phot.field[phot.idx] == 'PRIMO') | (phot.field[phot.idx] == 'WFC3-ERSII-G01') | (phot.field[phot.idx] == 'GEORGE')
    # use = phot.field[phot.idx] == 'GOODS-N'
    # 
    # dmag = m140-hmag
    # plt.plot(m140[use & fin], dmag[use & fin], marker='o', linestyle='None', alpha=0.5)
    # plt.plot([0,30],[0,0], color='black', alpha=0.5)
    # plt.xlim(15,25)
    # plt.ylim(-2,2)
    # 
    # plt.plot(phot.kron_radius[phot.idx][use & fin], dmag[use & fin], marker='o', linestyle='None', alpha=0.2)
    # xm, ym, ys, ns = threedhst.utils.runmed(phot.kron_radius[phot.idx][use & fin & (m140 < 23)], dmag[use & fin & (m140 < 23)], NBIN=30)
    # plt.plot(xm, ym, color='orange', linewidth=2)
    # 
    # plt.xlim(3,8)
    # plt.ylim(-2,2)
    # 
    # 
    # plt.plot(phot.kron_radius[phot.idx][use & fin], m140[use & fin], marker='o', linestyle='None', alpha=0.2)
    # plt.xlim(3,8)
    # plt.ylim(16,25)
    
    ########## H vs I-H
    field = phot.field[phot.idx] != 'xx'
    
    fields = {'COSMOS': ['COSMOS'], 'AEGIS': ['AEGIS'], 'GOODS-N':['GOODS-N'], 'GOODS-S':['GOODS-S','PRIMO','WFC3-ERSII-G01','GEORGE']}
    
    field_use = 'COSMOS'
    
    ix = 220
    
    fig = unicorn.catalogs.plot_init(square=True, xs=8, aspect=1./2, left=0.1, right=0.12, bottom=0.10, top=0.01, fontsize=10)
    fig.subplots_adjust(wspace=0.01,hspace=0.02,  left=0.05, right=0.94, bottom=0.10, top=0.98)
    
    #ax = fig.add_subplot(111)
    
    for field_use in ['AEGIS','COSMOS','GOODS-N','GOODS-S']:
        
        ix += 1
        ax = fig.add_subplot(ix)
        
        field = phot.field[phot.idx] == 'xx'
        print field_use
        for mat in fields[field_use]:
            field = field | (phot.field[phot.idx] == mat)

        ms = 6
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Times'

        ax2 = ax.twinx()

        ax.plot(m140[fin & keep & ~has_specz & field], IH[fin & keep & ~has_specz & field], marker='.', linestyle='None', color='black', alpha=0.1, ms=ms)
        ax.plot(m140[fin & keep & has_specz & field], IH[fin & keep & has_specz & field], marker='.', linestyle='None', color='green', alpha=0.5, ms=ms)
        #ax.plot(m140[fin & keep & field & (zout.z_peak[0::3] > 1)], IH[fin & keep & field & (zout.z_peak[0::3] > 1)], marker='o', linestyle='None', color='orange', alpha=0.5, ms=ms)
        ax.plot(np.array([10,30]), 22.5-np.array([10,30]), color='black', alpha=0.5, linewidth=3, linestyle='--')
        #ax.plot(np.array([10,30]), 24-np.array([10,30]), color='orange', alpha=0.5, linewidth=3)
        ax.plot(np.array([10,30]), [2.25, 2.25], color='purple', alpha=0.8, linewidth=3)

        #### Fraction histograms
        z1_red = IH > 2.25
        yh_a, xh_a = np.histogram(m140[fin & keep & field & ~z1_red], range=(16,25), bins=18)
        yh_z, xh_z = np.histogram(m140[fin & keep & field & has_specz & ~z1_red], range=(16,25), bins=18)
        show = yh_a > 0
        ax2.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='white', linewidth=4, alpha=0.7, linestyle='steps-mid')
        ax2.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='blue', linewidth=3, alpha=0.8, linestyle='steps-mid')

        yh_a, xh_a = np.histogram(m140[fin & keep & field & z1_red], range=(16,25), bins=18)
        yh_z, xh_z = np.histogram(m140[fin & keep & field & has_specz & z1_red], range=(16,25), bins=18)
        show = yh_a > 0
        ax2.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='white', linewidth=4, alpha=0.7, linestyle='steps-mid')
        ax2.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='red', linewidth=3, alpha=0.8, linestyle='steps-mid')

        ax.text(0.95, 0.88, field_use, transform=ax.transAxes, fontsize=12, backgroundcolor='white', horizontalalignment='right')
        if field_use == 'AEGIS':
            ax.text(19,3.5, r'$i=22.5$', horizontalalignment='center', verticalalignment='bottom', rotation=-30)
        
        if field_use == 'GOODS-N':
            ax.text(17.5,2.4, r'$\uparrow\ z > 1$, red $\uparrow$', horizontalalignment='left', verticalalignment='bottom')
            
        minorLocator   = MultipleLocator(1)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_xlim(17,24.5)
        ax.set_ylim(-0.5,5.2)

        minorLocator   = MultipleLocator(0.1)
        ax2.yaxis.set_minor_locator(minorLocator)
        ax2.set_xlim(17.1,24.5)
        ax2.set_ylim(-0.1,1.1)
        
        if field_use in ['AEGIS','COSMOS']:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(r'$m_{140}\ \sim\ H$')
            
        if field_use in ['COSMOS','GOODS-S']:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel(r'$(i-H)$')
            
        #
        if field_use in ['AEGIS','GOODS-N']:
            ax2.set_yticklabels([])
        else:
            ax2.set_ylabel(r'$f(z_\mathrm{spec})$')
        
        #fig.savefig('zspec_fraction_%s.pdf' %(field_use))
    
    fig.savefig('zspec_fraction_all.pdf')
    
    #plt.plot(zout.z_peak[0::3][fin], IH[fin], marker='o', linestyle='None', color='black', alpha=0.1, ms=4)
    # plt.plot(zout.z_peak[0::3][fin & keep & field], IH[fin & keep & field], marker='o', linestyle='None', color='black', alpha=0.1, ms=ms)
    # plt.plot(zout.z_peak[0::3][fin & keep & has_specz & field], IH[fin & keep & has_specz & field], marker='o', linestyle='None', color='blue', alpha=0.5, ms=ms)
    # plt.plot(np.array([0,30]), [2.25, 2.25], color='red', alpha=0.8, linewidth=3)
    # 
    # 
    # z1_red = IH > 2.25
    # yh_a, xh_a = np.histogram(zout.z_peak[0::3][fin & keep & field & ~z1_red], range=(0,4), bins=8)
    # yh_z, xh_z = np.histogram(zout.z_peak[0::3][fin & keep & field & has_specz & ~z1_red], range=(0,4), bins=8)
    # show = yh_a > 0
    # plt.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='blue', linewidth=3, alpha=0.8)
    # 
    # yh_a, xh_a = np.histogram(zout.z_peak[0::3][fin & keep & field & z1_red], range=(0,4), bins=8)
    # yh_z, xh_z = np.histogram(zout.z_peak[0::3][fin & keep & field & has_specz & z1_red], range=(0,4), bins=8)
    # show = yh_a > 0
    # plt.plot((xh_a[1:]+xh_a[:-1])[show]/2., (yh_z*1./yh_a)[show], color='red', linewidth=3, alpha=0.8)
    # 
    # plt.xlim(0,4)
    # plt.ylim(-0.5,5)
    
def find_brown_dwarf():
    import unicorn
    import unicorn.catalogs
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/SURVEY_PAPER')
        
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    int_lam = np.array([0.77e4, 1.25e4, 2.1e4])
    
    
    fp = open('stars_ijk.dat','w')
    fp.write('# id ii jj kk\n')
    
    ##### Known Brown dwarf
    object = 'AEGIS-3-G141_00195'
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s' %(object), OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/', CACHE_FILE = 'Same')
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    so = np.argsort(lci[~is_spec])
    yint = np.interp(int_lam, lci[~is_spec][so], fobs[~is_spec][so])/(int_lam/5500.)**2
    fp.write('%s %.3e %.3e %.3e\n' %(object, yint[0], yint[1], yint[2]))
    
    ###### Loop through all point sources
    stars = (phot.flux_radius[phot.idx] < 3) & (phot.mag_f1392w[phot.idx] < 24) & (mcat.rmatch[mcat.idx] < 0.5)
    for object in phot.id[phot.idx][stars]:
        print unicorn.noNewLine+'%s' %(object)
        try:
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s' %(object), OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'/ANALYSIS/REDSHIFT_FITS_v1.6/OUTPUT/', CACHE_FILE = 'Same')
        except:
            pass
        #
        dlam_spec = lci[-1]-lci[-2]
        is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        so = np.argsort(lci[~is_spec])
        yint = np.interp(int_lam, lci[~is_spec][so], fobs[~is_spec][so])/(int_lam/5500.)**2
        fp.write('%s %.3e %.3e %.3e\n' %(object, yint[0], yint[1], yint[2]))
        
    fp.close()
    
    if (1 == 0):
        ijk = catIO.Readfile('stars_ijk.dat')
        ij = -2.5*np.log10(ijk.ii/ijk.jj)
        jk = -2.5*np.log10(ijk.jj/ijk.kk)
        plt.plot(ij, jk, marker='o', markersize=3, color='black', alpha=0.8, linestyle='None')
        plt.plot(ij[0], jk[0], marker='o', markersize=8, color='red', alpha=0.5, linestyle='None')
        
        mat = ijk.id == 'AEGIS-11-G141_00314'
        plt.plot(ij[mat], jk[mat], marker='o', markersize=8, color='orange', alpha=0.5, linestyle='None')
        
        bd = phot.id[phot.idx] == 'x'
        for i, obj in enumerate(ijk.id):
            bd[phot.id[phot.idx] == obj] = (ij[i] > -0.) & (jk[i] < -1.7)
        
        #
        bdf = unicorn.analysis.BD_fit()
        
        for obj in phot.id[phot.idx][bd]:
            bdf.fit('/Users/gbrammer/Sites_GLOBAL/P/GRISM/ascii/%s.dat' %(obj), chi2_limit=100, trim_mtype=False, max_contam=0.8)
            
        
        unicorn.catalogs.make_selection_catalog(bd, filename='massive_lines.cat', make_html=True)
        os.system('rsync -avz massive_lines* ~/Sites_GLOBAL/P/GRISM_v1.6/ANALYSIS/')
            
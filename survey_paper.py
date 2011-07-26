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

import re

noNewLine = '\x1b[1A\x1b[1M'

root = None

def all_pointings():
    pointings(ROOT='GOODS-SOUTH')
    pointings(ROOT='COSMOS')
    pointings(ROOT='AEGIS')
    pointings(ROOT='UDS')
    pointings(ROOT='GOODS-N')
    
def pointings(ROOT='GOODS-SOUTH'):
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
    
    #### GOODS-S
    if ROOT=='GOODS-SOUTH':
        x0, x1 = 53.314005633802822,  52.886197183098595
        y0, y1 = -27.983151830808076, -27.654474431818176
        xticklab = [r'$3^\mathrm{h}33^\mathrm{m}00^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}30^\mathrm{s}$', r'$3^\mathrm{h}32^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(3,33,00, hours=True), degrees(3,32,30, hours=True), degrees(3,32,00, hours=True)]
        yticklab = [r'$-27^\circ40^\prime00^{\prime\prime}$', r'$45^\prime00^{\prime\prime}$', r'$-27^\circ50^\prime00^{\prime\prime}$', r'$55^\prime00^{\prime\prime}$']
        ytickv = [degrees(-27, 40, 00, hours=False), degrees(-27, 45, 00, hours=False), degrees(-27, 50, 00, hours=False), degrees(-27, 55, 00, hours=False)]
        
    #### COSMOS
    if ROOT=='COSMOS':
        x1, x0 = 149.99120563380279, 150.23823661971829
        y0, y1 = 2.1678109478476815, 2.5996973302980129
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$00^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(10,00,30, hours=True), degrees(10,00,00, hours=True)]
        yticklab = [r'$+2^\circ15^\prime00^{\prime\prime}$', r'$25^\prime00^{\prime\prime}$', r'$35^\prime00^{\prime\prime}$']
        ytickv = [degrees(02, 15, 00, hours=False), degrees(02, 25, 00, hours=False), degrees(02, 35, 00, hours=False)]
    
    #### AEGIS
    if ROOT=='AEGIS':
        x1, x0 = 214.49707154104345, 215.12704734584406
        y0, y1 = 52.680946433013482, 53.01597137966467
        xticklab = [r'$18^\mathrm{m}00^\mathrm{s}$', r'$14^\mathrm{h}19^\mathrm{m}00^\mathrm{s}$', r'$20^\mathrm{m}00^\mathrm{s}$']
        xtickv = [degrees(14,18,00, hours=True), degrees(14,19,00, hours=True), degrees(14,20,00, hours=True)]
        yticklab = [r'$+52^\circ45^\prime00^{\prime\prime}$', r'$50^\prime00^{\prime\prime}$', r'$55^\prime00^{\prime\prime}$']
        ytickv = [degrees(52, 45, 00, hours=False), degrees(52, 50, 00, hours=False), degrees(52, 55, 00, hours=False)]
    
    #### UDS
    if ROOT=='UDS':
        x1, x0 = 34.116935194128146, 34.51871547581829
        y0, y1 = -5.2957542206957582, -5.0834327182123147+2./3600
        xticklab = [r'$18^\mathrm{m}00^\mathrm{s}$', r'$2^\mathrm{h}17^\mathrm{m}30^\mathrm{s}$', r'$17^\mathrm{m}00^\mathrm{s}$', r'$16^\mathrm{m}30^\mathrm{s}$']
        xtickv = [degrees(2,18,00, hours=True), degrees(2,17,30, hours=True), degrees(2,17,00, hours=True), degrees(2,16,30, hours=True)]
        yticklab = [r'$05^\prime00^{\prime\prime}$', r'$-5^\circ10^\prime00^{\prime\prime}$', r'$15^\prime00^{\prime\prime}$']
        ytickv = [degrees(-5, 05, 00, hours=False), degrees(-5, 10, 00, hours=False), degrees(-5, 15, 00, hours=False)]
    
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
        
    #### Make square for given plot dimensions
    dx = np.abs(x1-x0)*np.cos(y0/360*2*np.pi)
    dy = (y1-y0)
    
    fig = unicorn.catalogs.plot_init(square=True, xs=7*dx/dx_ref, aspect=dy/dx)
    #fig = unicorn.catalogs.plot_init(square=True)
    
    ax = fig.add_subplot(111)
      
    files=glob.glob(unicorn.GRISM_HOME+'REGIONS/'+ROOT+'-[0-9]*reg')
    for file in files[15:16]:
        #
        field = re.split('-[0-9]', file)[0]
        pointing = file.split(field+'-')[1].split('.reg')[0]
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        #
        wfcx, wfcy = sup.polysplit(lines[1])
        fi = ax.fill(wfcx, wfcy, alpha=0.2, color=wfc3_color)
        #
        if acs_color is not None:
            acsx1, acsy1 = sup.polysplit(lines[2])
            acsx2, acsy2 = sup.polysplit(lines[3])
            #
            fi = ax.fill(acsx1, acsy1, alpha=0.05, color=acs_color)
            fi = ax.fill(acsx2, acsy2, alpha=0.05, color=acs_color)
            pl = ax.plot(acsx1, acsy1, alpha=0.1, color=acs_color)
            pl = ax.plot(acsx2, acsy2, alpha=0.1, color=acs_color)
        #
        te = ax.text(np.mean(wfcx[:-1]), np.mean(wfcy[:-1]), pointing, va='center', ha='center', fontsize=8)
    
    #    
    if yticklab is not None:
        ax.set_xticklabels(xticklab)
        xtick = ax.set_xticks(xtickv)
        ax.set_yticklabels(yticklab)
        ytick = ax.set_yticks(ytickv)
    
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\delta$')
    
    ax.text(0.95, 0.05,ROOT,
         horizontalalignment='right',
         verticalalignment='bottom',
         transform = ax.transAxes, fontsize='14')
    
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    
    print 'RA  - ', hexagesimal(x0), hexagesimal(x1)
    print 'Dec - ', hexagesimal(y0, hours=False), hexagesimal(y1, hours=False)
    
    plt.savefig('%s_pointings.pdf' %(ROOT))

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

def polysplit(region='polygon(150.099223,2.391097,150.086084,2.422515,150.050573,2.407277,150.064586,2.376241)'):
    spl = region[region.find('(')+1:region.find(')')].split(',')
    px = spl[0::2]
    py = spl[1::2]
    px.append(px[0])
    py.append(py[0])
    return np.cast[float](px), np.cast[float](py)
    
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
    
    
    make_plot(root=root)
    
def make_plot(root='AEGIS-11', range1=(0.90,1.08), range2=(-0.02, 0.02)):
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
    
    med = threedhst.utils.biweight(sky[1].data, mean=True)

    ysize = 3.
    #fig = plt.figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*4,ysize], dpi=100)
    
    top_panel = 0.2
    NPANEL = 4
    
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*NPANEL*(1-top_panel),ysize],dpi=100)
    else:
        fig = Figure(figsize=[ysize*im_shape[1]*1./im_shape[0]*NPANEL*(1-top_panel),ysize], dpi=100)
    
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.01,
                        bottom=0.07,right=0.99,top=0.97)
    #
    
    vmin, vmax = -0.2, 0.1
    
    x0 = 0.005
    y0 = x0
    dx = (1.-(NPANEL+1)*x0)/NPANEL
    
    #ax = fig.add_subplot(141)
    ax = fig.add_axes(((x0+(dx+x0)*0), y0, dx, 1-top_panel-y0))
    
    ax.imshow(0-(first[1].data-med), interpolation='nearest',aspect='auto',vmin=vmin-0.1,vmax=vmax+0.15)    
    axis_imshow(ax, text='a)\ Raw')
    ax.text(0.12, 0.85, r'$\mathrm{%s}$' %(root), horizontalalignment='left', verticalalignment='center',
                 transform = ax.transAxes, color='black', fontsize=14)
    #
    #show_limits(ax, -(vmax+0.15)+med, -(vmin-0.1)+med)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*0), (1-top_panel), dx, top_panel-2*y0))
    pp = sup.first_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.first_prof[i])
        pp += sup.first_prof[i]
    ax.plot(pp/4., color='black')
    axis_profile(ax, yrange=range1, text='a)\ Raw')
    
    #ax = fig.add_subplot(142)
    ax = fig.add_axes(((x0+(dx+x0)*1), y0, dx, 1-top_panel-y0))
    ax.imshow(0-(flat[1].data-med), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    axis_imshow(ax, text='b)\ Flat')
    #show_limits(ax, -vmax+med, -vmin+med)

    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*1), (1-top_panel), dx, top_panel-2*y0))
    pp = sup.flat_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.flat_prof[i])
        pp += sup.flat_prof[i]
    ax.plot(pp/4., color='black')
    axis_profile(ax, yrange=range1, text='b)\ Flat')
    
    #ax = fig.add_subplot(143)
    ax = fig.add_axes(((x0+(dx+x0)*2), y0, dx, 1-top_panel-y0))
    ax.imshow(0-(sky[1].data-med), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    axis_imshow(ax, text='c)\ Sky')
    #show_limits(ax, -vmax+med, -vmin+med)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*2), (1-top_panel), dx, top_panel-2*y0))
    pp = sup.sky_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.sky_prof[i])
        pp += sup.sky_prof[i]
    ax.plot(pp/4., color='black')
    axis_profile(ax, yrange=range1, text='c)\ Sky')

    #ax = fig.add_subplot(144)
    ax = fig.add_axes(((x0+(dx+x0)*3), y0, dx, 1-top_panel-y0))
    ax.imshow(0-(final[1].data), interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax)    
    axis_imshow(ax, text='d)\ Final')
    #show_limits(ax, -vmax, -vmin)
    
    #### Show profiles
    ax = fig.add_axes(((x0+(dx+x0)*3), (1-top_panel), dx, top_panel-2*y0))
    pp = sup.final_prof[0]*0.
    for i in range(4):
        #ax.plot(sup.final_prof[i])
        pp += sup.final_prof[i]
    ax.plot(pp/4., color='black')
    axis_profile(ax, yrange=range2, text='d)\ Final')
    
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
    
    ax.text(0.98, 0.08,r'$\mathrm{[%.3f,\ %.3f]}$' %(ylimits[0], ylimits[1]),
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

    
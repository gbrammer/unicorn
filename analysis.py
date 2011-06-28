import os
import pyfits
import numpy as np
import glob

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
 
def get_grism_path(root):
    """ 
    Given a rootname for a grism pointing, get the path to the 
    grism files
    """
    PATH = './'
    if root.startswith('COSMOS'):
        PATH = unicorn.GRISM_HOME+'COSMOS/'
    if root.startswith('AEGIS'):
        PATH = unicorn.GRISM_HOME+'AEGIS/'
    if root.startswith('GOODS-N'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('GOODS-S'):
        PATH = unicorn.GRISM_HOME+'GOODS-S/'
    if root.startswith('MARSHALL'):
        PATH = unicorn.GRISM_HOME+'SN-MARSHALL/'
    if root.startswith('PRIMO'):
        PATH = unicorn.GRISM_HOME+'SN-PRIMO/'
    if root.startswith('GEORGE'):
        PATH = unicorn.GRISM_HOME+'SN-GEORGE/'
    
    return PATH
    
def read_catalogs(root='', cosmos=False, aegis=False, goodsn=False, cdfs=False, ecdfs=False, uds=False):
    """
    
    Read photometry, redshift, SPS catalogs for a given field.
    
    Usage: cat, zout, fout = read_catalogs(cosmos=True)
    
    """        
    import threedhst.catIO as catIO
    import numpy as np
    
    KTOT_COL = None
    ### Optionally supply a grism basename and determine the field from that
    if root.startswith('COSMOS'):
        cosmos=True
    if root.startswith('AEGIS'):
        aegis=True
    if root.startswith('GOODS-N'):
        goodsn=True
    if root.startswith('GOODS-S'):
        cdfs=True
    if root.startswith('MARSHALL'):
        uds=True
    if root.startswith('PRIMO'):
        cdfs=True
    if root.startswith('GEORGE'):
        cdfs=True
    
    aegis_wirds=False
    if root.startswith('AEGIS-11') | root.startswith('AEGIS-2-') | root.startswith('AEGIS-1-'):
        aegis=False
        aegis_wirds=True
    
    CAT_FILE = None
    ZOUT_FILE = None
    FOUT_FILE = None
    
    MAGS = False
    
    #### Define paths
    if cosmos:
        GRISM_PATH=unicorn.GRISM_HOME+'COSMOS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'        
        CAT_FILE = CAT_PATH + '../cosmos-1.v4.6.cat'
        ZOUT_FILE = CAT_PATH + 'cosmos-1.v4.6.zout'
        FOUT_FILE = CAT_PATH+'../cosmos-1.bc03.v4.6.fout'
        KTOT_COL = 'ktot'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/COSMOS/NMBS/Photometry/'
            CAT_FILE = CAT_PATH + 'cosmos-1.deblend.v5.1.cat'
            ZOUT_FILE = CAT_PATH + '/cosmos-1.deblend.redshifts/cosmos-1.deblend.v5.1.zout'
            FOUT_FILE = CAT_PATH+'cosmos-1.deblend.sps/cosmos-1.bc03.del.deblend.v5.1.fout'
            KTOT_COL = 'K'
            
    if aegis:
        GRISM_PATH=unicorn.GRISM_HOME+'AEGIS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
        CAT_FILE = CAT_PATH + '../aegis-n2.v4.6.cat'
        ZOUT_FILE = CAT_PATH + 'aegis-n2.v4.6.zout'
        FOUT_FILE = CAT_PATH+'/../aegis-n2.bc03.v4.6.fout'
        KTOT_COL = 'ktot'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/AEGIS/NMBS/Photometry/'
            CAT_FILE = CAT_PATH + 'aegis-n2.deblend.v5.1.cat'
            ZOUT_FILE = CAT_PATH + '/aegis-n2.deblend.redshifts/aegis-n2.deblend.v5.1.zout'
            FOUT_FILE = CAT_PATH+'aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout'
            KTOT_COL = 'K'
    
    if aegis_wirds:
        GRISM_PATH=unicorn.GRISM_HOME+'AEGIS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/WIRDS/FAST/'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/AEGIS/WIRDS/Photometry/FAST/'
        #
        ZOUT_FILE = CAT_PATH + '../EAZY/OUTPUT/egs_candels.zout'
        FOUT_FILE = CAT_PATH+'egs_candels.fout'
        CAT_FILE = CAT_PATH + 'egs_candels.cat'
        KTOT_COL = 'ks'
        CAT_FILE = CAT_PATH + '../WIRDS_D3-95_Ks_ugrizJHKs_141927+524056_T0002.cat.candels'
        KTOT_COL = 'kstot'
        MAGS = True
    
    if goodsn:
        GRISM_PATH=unicorn.GRISM_HOME+'GOODS-N/'
        CAT_PATH = '/research/HST/GRISM/3DHST/GOODS-N/MODS/FAST/'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/GOODS-N/MODS/Photometry/FAST/'
        #
        CAT_FILE = CAT_PATH+'mods.cat'
        ZOUT_FILE = CAT_PATH+'../EAZY/OUTPUT/photz.zout'
        FOUT_FILE = CAT_PATH+'mods.bc03.fout'
        KTOT_COL = 'ks_ap'
        
    if cdfs:
        GRISM_PATH=unicorn.GRISM_HOME+'GOODS-S/'
        CAT_PATH = GRISM_PATH+'FIREWORKS/'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/GOODS-S/FIREWORKS/FAST/'
        #
        CAT_FILE = CAT_PATH+'fireworks.cat'
        ZOUT_FILE = CAT_PATH+'fireworks.zout'
        FOUT_FILE = CAT_PATH + 'fireworks.fout'
        KTOT_COL = 'Ks_totf'
        
    if uds:
        GRISM_PATH=unicorn.GRISM_HOME+'SN-MARSHALL/'
        CAT_PATH = GRISM_PATH+'UDSPHOT/'
        if unicorn.hostname().startswith('uni'):
            CAT_PATH = '/3DHST/Ancillary/UDS/UKIDSS/FAST/'
        #
        CAT_FILE = CAT_PATH+'uds.cat'
        ZOUT_FILE = CAT_PATH+'uds.zout'
        FOUT_FILE = CAT_PATH + 'uds.fout'
        KTOT_COL = 'K_totf'
    
    if KTOT_COL is None:
        """
        All catalog flags are False
        """
        return None, None, None
    
    #### Read the catalogs
    cat = catIO.ReadASCIICat(CAT_FILE)
    if MAGS:
        cat.kmag = cat.field(KTOT_COL)
    else:
        cat.kmag = 25.-2.5*np.log10(cat.field(KTOT_COL))
    
    cat.MAGS_UNIT = MAGS
    
    if 'star_flag' in cat.names:
        cat.star_flag = cat.field('star_flag')
    else:
        cat.star_flag = cat.kmag*0
    
    #### Different column names
    if cdfs:
        cat.ra = cat.field('RA')
        cat.dec = cat.field('DEC')
        cat.id = np.cast[int](cat.field('ID'))
    
    cat.filename=CAT_FILE
    
    if ZOUT_FILE:
        zout = catIO.ReadASCIICat(ZOUT_FILE)
        zout.filename=ZOUT_FILE
    else:
        zout = None

    if FOUT_FILE:
        fout = catIO.ReadASCIICat(FOUT_FILE)
        fout.filename=FOUT_FILE
    else:
        fout = None
    
    return cat, zout, fout

def read_grism_files(root='COSMOS-3-G141', BASE_PATH='', GRISM_NAME='G141'):
    """
    Read root+'_drz.cat' and the associated SPC file.
    """
    import threedhst
    
    grismCat, SPC = None, None
    
    if not BASE_PATH:
        BASE_PATH = get_grism_path(root)
        
    ##### catalog
    grismCat = threedhst.sex.mySexCat(BASE_PATH+'DATA/'+root+'_drz.cat')
    for col in grismCat.column_names:
        if col.startswith('MAG_F'):
            grismCat.MAG = grismCat[col]
            grismCat.DETECT_FILTER = col
            break
            
    ##### SPC file        
    try:
        SPC = threedhst.plotting.SPCFile(root+'_2_opt.SPC.fits',
                axe_drizzle_dir=BASE_PATH+'DRIZZLE_'+GRISM_NAME)
    except:
        SPC = None
        
    return grismCat, SPC
    
def make_SED_plots(grism_root='COSMOS-3-G141'):
    import unicorn.analysis
    
    PATH = unicorn.analysis.get_grism_path(grism_root)
    print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root)
    
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root)
        
    print 'Matched catalog'
    unicorn.analysis.match_grism_to_phot(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat',
                  OUTPUT_DIRECTORY=OUTPUT_DIRECTORY,
                  MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE)
    
    ## make figures
    for id in grismCat.id:
        unicorn.analysis.specphot(id=id, grism_root=grism_root, SPC = SPC, 
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')
    
def convolveWithThumb(id, lambdaz, temp_sed, SPC, oned=True, xint=None, verbose=False):
    """ 
    
    Convolve the best-fit eazy spectrum with the object shape
    
    """
    from scipy import convolve as conv
    
    thumb_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_thumb.fits.gz'
    try:
        thumb = pyfits.open(thumb_file)
    except:
        return lambdaz, temp_sed
    
    size = thumb[0].data.shape
    DLAM = np.sqrt(thumb[0].header['CD1_1']**2+thumb[0].header['CD1_2']**2)*3600./0.128254*46.5
    
    twod_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_2D.fits.gz'
    twod = pyfits.open(twod_file)
    model1D = np.matrix(twod[5].data.sum(axis=1))
    model2D = np.array(np.dot(np.transpose(model1D),np.ones((1,size[0]))))
    
    if oned:
        profile = np.sum(thumb[0].data*model2D,axis=0)
        profile /= profile.sum()
        if verbose: 
            xprof = np.arange(len(profile))
            print 'Profile center:', np.trapz(xprof, xprof*profile)/np.trapz(xprof, profile), np.mean(xprof)
    else:
        profile = thumb[0].data*model2D
        #for i in range(size[0]):
        #    profile[i,:] /= profile[i,:].sum()
        profile /= profile.sum()
        
    LSM = size[0]/2
    xmin = 3000
    xmax = 2.4e4
    
    q = np.where((lambdaz > (xmin-LSM)) & (lambdaz < (xmax+LSM)))[0]
    lambdaz = lambdaz[q]
    temp_sed = temp_sed[q]
    
    ### convolve with some factor of the object size
    #LSM = np.sqrt(LSM**2+(0.1*R[grism_idx]/0.128*46.5)**2)
    
    temp_sed_sm = temp_sed*1.
    # for i in range(len(lambdaz)):
    #     temp_sed_sm[i] = np.trapz(1./np.sqrt(2*np.pi*50**2)*np.exp(-0.5*(lambdaz-lambdaz[i])**2/50**2)*temp_sed, lambdaz)
    
    if verbose:
        print 'DLAM: ', DLAM
        
    #
    if xint is None:
        xint = np.arange(xmin-LSM*DLAM,xmax+LSM*DLAM,DLAM)
    
    yint = np.interp(xint, lambdaz, temp_sed_sm)
    
    #### Convolve with a gaussian
    xgauss = np.arange(20)*DLAM-10*DLAM
    ygauss = np.exp(-1*xgauss**2/2/50**2)
    ygauss /= np.sum(ygauss)
    yint = conv(yint, ygauss, mode='same')
    
    if oned:
        temp_sed_conv = conv(yint, profile, mode='same')
        # temp_sed_conv = yint*0.
        # for i in range(LSM,len(xint)-LSM):
        #     temp_sed_conv[i] = np.sum(yint[i-LSM:i+LSM]*profile)
    else:
        NX, NY = len(yint), size[0]
        temp_sed_conv = np.zeros((NY,NX))
        for i in range(size[0]):
            temp_sed_conv[i,:] = conv(yint, profile[i,:].flatten(), mode='same') #np.dot(yint[i-LSM:i+LSM],profile).reshape((NY,))
            
    
    thumb.close()
    
    return xint, temp_sed_conv
    
def specphot(id=69, grism_root='ibhm45030',
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
    CACHE_FILE = 'Same', Verbose=False,
    SPC = None, cat=None, grismCat = None,
    zout = None, fout = None, OUT_PATH='/tmp/', OUT_FILE_FORMAT=True,
    OUT_FILE='junk.png', GET_SPEC_ONLY=False):
    """
specphot(id)
    
    Get photometry/SED fit and compare G141 spectrum
    """
    #import scipy.interpolate as interpol
       
    ### 69, 54!
    
    xxx = """
    id=199
    grism_root='ibhm48'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    CACHE_FILE = 'Same'
    """
    
    #### Get G141 spectrum
    if Verbose:
        print 'Read SPC'
    
    if SPC is None:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                    axe_drizzle_dir='DRIZZLE_G141')
                    
    spec = SPC.getSpec(id)
    if spec is False:
        return False
        
    xmin = 3000
    xmax = 2.4e4
    
    lam = spec.field('LAMBDA')
    flux = spec.field('FLUX')
    ffix = flux-spec.field('CONTAM')
    ferr = spec.field('FERROR') #*0.06/0.128254
        
    if Verbose:
        print 'Read grism catalog'
        
    #### Read the grism catalog and get coords of desired object
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    #### Source size
    R = np.sqrt(np.cast[float](grismCat.A_IMAGE)*np.cast[float](grismCat.B_IMAGE))
    grism_idx = np.where(grismCat.id == id)[0][0]
    
    Rmatch = R[grism_idx]*1.
    
    #print 'R=%f"' %(Rmatch)
    ra0 = grismCat.ra[grismCat.id == id][0]
    de0 = grismCat.dec[grismCat.id == id][0]
    
    #### Read EAZY outputs and get info for desired object
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
    
    
    photom_idx = np.where(dr == np.min(dr))[0][0]
    
    drMatch = dr[photom_idx]*1.
    #print 'dr = %7.2f\n' %(drMatch)
    #print drMatch, np.min(dr)
    
    if drMatch > 2:
        return False
        
    if Verbose:
        print 'Read zout'
    if zout is None:    
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
        
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    if Verbose:
        print 'Read binaries'
        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(photom_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE)
         
    try:
        lambdaz, temp_sed_sm = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
    except:
        temp_sed_sm = temp_sed*1.
        
    if Verbose: 
        print 'Normalize spectrum'
        
    #### Normalize G141 spectrum
    #interp = interpol.interp1d(lambdaz, temp_sed_sm, kind='linear')

    q = np.where((lam > 1.08e4) & (lam < 1.68e4) & (flux > 0))[0]
    #### G102
    if lam.min() < 9000:
        q = np.where((lam > 0.8e4) & (lam < 1.13e4) & (flux > 0))[0]
    
    #### ACS G800L
    if lam.min() < 5000:
        q = np.where((lam > 0.55e4) & (lam < 1.0e4) & (flux > 0))[0]
        
    if len(q) == 0:
        return False

    yint = np.interp(lam[q], lambdaz, temp_sed_sm)
        
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    
    if GET_SPEC_ONLY:
        return lam, ffix*anorm, np.sqrt(ferr**2+spec.field('CONTAM')**2)*anorm, lci, fobs, efobs, photom_idx
         
    if Verbose:
        print 'Start plot'
        
    #### Make the plot
    threedhst.plotting.defaultPlotParameters()
    
    xs=5.8
    ys = xs/4.8*3.2
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[xs,ys],dpi=100)
    else:
        fig = Figure(figsize=[xs,ys], dpi=100)
    
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.13*4.8/xs,
                        bottom=0.15*4.8/xs,right=1.-0.02*4.8/xs,top=1-0.10*4.8/xs)
    
    ax = fig.add_subplot(111)
    
    ymax = np.max((ffix[q])*anorm)
    
    if Verbose:
        print 'Make the plot'
        
    ax.plot(lambdaz, temp_sed_sm, color='red')
    # plt.errorbar(lam[q], ffix[q]*anorm, yerr=ferr[q]*anorm, color='blue', alpha=0.8)
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.2, linewidth=1)
    
    #### Show own extraction
    sp1d = threedhst.spec1d.extract1D(id, root=grism_root, path='./HTML', show=False, out2d=False)
    lam = sp1d['lam']
    flux = sp1d['flux']
    ffix = sp1d['flux']-sp1d['contam'] #-sp1d['background']
    ferr = sp1d['error']
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.6, linewidth=1)
    
    #### Show photometry + eazy template
    ax.errorbar(lci, fobs, yerr=efobs, color='orange', marker='o', markersize=10, linestyle='None', alpha=0.4)
    ax.plot(lambdaz, temp_sed_sm, color='red', alpha=0.4)

    ax.set_ylabel(r'$f_{\lambda}$')
    
    if plt.rcParams['text.usetex']:
        ax.set_xlabel(r'$\lambda$ [\AA]')
        ax.set_title('%s: \#%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
    else:
        ax.set_xlabel(r'$\lambda$ [$\AA$]')
        ax.set_title('%s: #%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
        
    #kmag = 25-2.5*np.log10(cat.ktot[photom_idx])
    kmag = cat.kmag[photom_idx]
    
    ##### Labels
    label = 'ID='+r'%s   K=%4.1f  $\log M$=%4.1f' %(np.int(cat.id[photom_idx]),
        kmag, fout.field('lmass')[photom_idx])
        
    ax.text(5e3,1.08*ymax, label, horizontalalignment='left',
      verticalalignment='bottom')
    
    
    label = 'R=%4.1f"' %(drMatch)
    if drMatch > 1.1:
        label_color = 'red'
    else:
        label_color = 'black'
    ax.text(2.2e4,1.08*ymax, label, horizontalalignment='right',
      color=label_color, verticalalignment='bottom')
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.1*ymax,1.2*ymax)
    
    if Verbose:
        print 'Save the plot'
    
    if OUT_FILE_FORMAT:
        out_file = '%s_%05d_SED.png' %(grism_root, id)
    else:
        out_file = OUT_FILE
        
    if USE_PLOT_GUI:
        fig.savefig(OUT_PATH+'/'+out_file,dpi=100,transparent=False)
        plt.close()
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(OUT_PATH+'/'+out_file, dpi=100, transparent=False)
    
    noNewLine = '\x1b[1A\x1b[1M'
    print noNewLine+OUT_PATH+'/'+out_file
    
    if Verbose:
        print 'Close the plot window'
        
    
def match_grism_to_phot(grism_root='ibhm45',
        MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
        OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
        CACHE_FILE = 'Same', Verbose=False,
        SPC = None, cat=None, grismCat = None,
        zout = None, fout = None, OUTPUT='/tmp/match.cat'):

    
    import unicorn.analysis
    # 
    # MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    # OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    # CACHE_FILE = 'Same'
    
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    if zout is None:
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    #rfUV = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.153-155.rf')
    
    z_peak = zout.field('z_peak')
    # kmag = 25-2.5*np.log10(cat.field('ktot'))
    kmag = cat.kmag
    
    #uv = -2.5*np.log10(rfUV.field('L153')/rfUV.field('L155'))
    
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    drs = grismCat.ra*1.
    ids = np.cast[int](grismCat.ra*1.)
    
    cosfact = np.cos(np.median(grismCat.dec)/360*2*np.pi)
    
    for i in range(len(drs)):
        dr = np.sqrt((cat.ra-grismCat.ra[i])**2*cosfact**2+
                     (cat.dec-grismCat.dec[i])**2)*3600
        mat = np.where(dr == np.min(dr))[0][0]
        drs[i] = dr[mat]
        ids[i] = mat
    
    fp = open(OUTPUT,'w')
    fp.write("# id_f140w  mag_f140w  id_phot  mag_Ktot  Rmatch  z_peak  logM star_flag fcontam\n# Rmatch = match distance in arcsec (should be < 1)\n# %s\n" %(MAIN_OUTPUT_FILE))
    for i in range(len(drs)):
        j = ids[i]
        line = "%5d %6.2f %8d %6.2f %6.2f %6.2f %6.2f %d %5s\n" %(grismCat.id[i],
                  np.float(grismCat.MAG[i]), np.int(cat.id[j]),
                  kmag[j], drs[i], 
                  z_peak[j], fout.field('lmass')[j], 
                  np.int(cat.star_flag[j]), grismCat['FCONTAM'][i])
        if grismCat.id[i] in SPC._ext_map:
            fp.write(line)
    fp.close()

def make_multiple_fluximage(grism_root='COSMOS-3-G141'):
    import unicorn.analysis
    waves = [1.1e4,1.25e4,1.6e4]
    for wave in waves:
        unicorn.analysis.make_fluximage(grism_root=grism_root, wavelength=wave)

def make_fluximage(grism_root='COSMOS-3-G141', wavelength=1.1e4, direct_image=None, match_toler=1, verbose=True):
    """
    1) Read a SExtractor 3DHST catalog
    2) Open the corresponding direct image
    3) Match each object in the 3DHST catalog to an external catalog, which
       has redshifts / SED fits
    4) Get the ('wavelength' - detection band) color for each object and 
       scale its segmentation image accordingly
    5) If not match within 'match_toler', just use the detected flux
    6) Write out the new scaled image
    """
    
    out_image = 'DATA/'+grism_root.replace('G141','f%03d' %(wavelength/100))+'.fits'
    
    ##### Get the path and read the catalogs
    PATH = unicorn.analysis.get_grism_path(grism_root)
    PWD=os.getcwd()
    print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root)
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root)
    
    detect_wlen = np.float(grismCat.DETECT_FILTER.strip('MAG_FW'))*10
    
    #### Detection and segmentation images
    if direct_image is None:
        direct_image = glob.glob('PREP_FLT/'+grism_root.replace('G141','*')+'_drz.fits')[0]
    
    seg_file = grismCat.filename.replace('drz.cat','seg.fits')
    seg_file = threedhst.utils.find_fits_gz(seg_file)
    
    direct = pyfits.open(direct_image)
    seg = pyfits.open(seg_file)
    
    #### Loop through objects in the catalog
    cosfact = np.cos(np.median(grismCat.dec)/360*2*np.pi)
    xint = np.array([wavelength, detect_wlen])

    #### If a Grism SPC file exists, only use IDs defined there
    #### Otherwise use all objects in the SEx. catalog
    if SPC is not None:
        ids = SPC._ext_map
    else:
        ids = grismCat.id
        
    noNewLine = '\x1b[1A\x1b[1M'
    for j, id in enumerate(ids):
        progress = '%2d' %(np.int(j*100./len(ids))) + '%'
        print noNewLine+out_image+':  '+progress
            
        i = np.arange(grismCat.nrows)[grismCat.id == id][0]
           
        dr = np.sqrt((cat.ra-grismCat.ra[i])**2*cosfact**2+
                         (cat.dec-grismCat.dec[i])**2)*3600
        mat = np.where(dr == np.min(dr))[0][0]
        
        scale = 1.
        
        if dr[mat] < match_toler:
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
                eazy.getEazySED(mat, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                     OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                     CACHE_FILE = 'Same')
                        
            #### Smooth to grism resolution
            try:
                lambdaz, temp_sed = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
            except:
                pass
                # lambdaz2, temp_sed2 = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC, oned=False)
            
            yint = np.interp(xint, lambdaz, temp_sed)
            scale = yint[0]/yint[1]    
            direct[1].data[seg[0].data == id] *= scale
    
    ### f_nu
    direct[1].data *= wavelength**2/detect_wlen**2
    
    ### Write the image, but keep just the SCI extension
    direct[1].writeto('/tmp/fluximage.fits', clobber=True)
    try:
        os.remove(out_image)
    except:
        pass
    
    os.chdir(PWD)
    iraf.imcopy('/tmp/fluximage.fits[1]',out_image)
    
    return out_image
        
        
def show_massive_galaxies(masslim=10.5, maglim=23.5, zrange=(0,5), 
    use_kmag=False, contam=0.5, coverage=0.9):        
    
    if unicorn.hostname().startswith('unicorn'):
        os.chdir('/Library/WebServer/Documents/P/GRISM_v1.5/ANALYSIS')
        scripts="../scripts"
        matches = glob.glob('../SED/*match.cat')
    
    else:
        os.chdir(unicorn.GRISM_HOME+'ANALYSIS')
        scripts="http://localhost/~gbrammer/COSMOS/scripts"

        matches = []
        matches.extend(glob.glob('../AEGIS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../COSMOS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../GOODS-S/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../GOODS-N/HTML/SED/*match.cat'))
        # matches.extend(glob.glob('../SN-GEORGE/HTML/SED/*match.cat'))
        # matches.extend(glob.glob('../SN-PRIMO/HTML/SED/*match.cat'))    
        # matches.extend(glob.glob('../SN-MARSHALL/HTML/SED/*match.cat'))
        # matches.extend(glob.glob('../ERS/HTML/SED/*match.cat'))
    
    print matches
    
    fplist = open('massive.dat','w')
    fplist.write('# ID   ra   dec  z lmass\n')
    fplist.write('# mag > %.1f, mass > %.1f, z=(%.1f,%.1f), contam=%.2f\n' %(maglim, masslim, zrange[0], zrange[1], contam))
    
    fpreg = open('massive.reg', 'w')
    fpreg.write('fk5\n')
    
    fp = open('massive.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="%s/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="%s/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="%s/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[2,2]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                        6: {
                                sorter: false
                        },
                        7: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> Mag_WFC3 </th>
        <th> z </th>
        <th> logM </th>
        <th> Thumb </th>
        <th> 2D </th>
        <th> 1D </th>
        <th> SED </th>
    </thead>
    <tbody>
    """ %(scripts, scripts, scripts))
    
    NUSE = 0
    for match in matches:
        root_path = os.path.dirname(match).split('SED')[0]
        root = os.path.basename(match).split('_match.cat')[0]
        
        c = catIO.Readfile(match)
        xml_file = root_path+root+'.xml'
        xml = catIO.readMarkerXML(xml_file)
        #print c.keys()
        
        #### SExtractor + contamination catalog
        sex = threedhst.sex.mySexCat(root_path+root+'_drz.cat')
        mat = c.id_f140w*0
        for i in range(len(c.id_f140w)):
            mat[i] = np.where(np.cast[int](sex.NUMBER) == c.id_f140w[i])[0][0]
        
        try:
            fcover = np.cast[float](sex.FCOVER)[mat]
        except:
            fcover = c.id_f140w*0.+1
            
        select_mag = c.mag_f140w
        if use_kmag:
            #print c.keys()
            select_mag = c.mag_ktot
        
        #### Had used wrong zeropoint in FAST    
        # if ('ERS' in root_path) | ('GOODS-S' in root_path):
        #     c.logm -= 0*-0.4*(23.86-25)
            
        use = (c.star_flag < 1) & (c.logm > masslim) & (c.rmatch < 1) & (select_mag < maglim) & (c.z_peak > zrange[0]) & (c.z_peak < zrange[1]) & (fcover > coverage)
        
        if 'fcontam' in c.keys():
            if contam < 0:
                use = use & (c.fcontam > -contam)
            else:
                use = use & (c.fcontam < contam)

        use = np.where(use == True)[0]
        use = use[np.argsort(c.logm[use])]
        NUSE += len(use)
        
        for i in use:
            ID = c.id_f140w[i]
            
            file="%s_%05d" %(root, ID)
            
            print "%4d %4d %s" %(ID, len(xml), file)
            
            fplist.write("%15s %14.6f %14.6f %8.3f %5.2f\n" %(file, xml[ID].ra, xml[ID].dec, c.z_peak[i], c.logm[i]))
            
            fpreg.write('circle(%14.6f,%14.6f,1") # text={%s} color=magenta  \n' %(xml[ID].ra, xml[ID].dec, file))
            
            fp.write("""
        <tr>
            <td> %s<br> %13.6f %13.6f <br> %4.1f </td>
            <td> %5.1f </td>
            <td> %5.2f </td>
            <td> %5.2f </td>
            <td> <img src=%s/images/%s_thumb.png height=180px> </td>
            <td> <img src=%s/images/%s_2D.png height=180px> </td>
            <td> <img src=%s/images/%s_1D.png height=180px> </td>
            <td> <img src=%s/SED/%s_SED.png height=180px> </td>
        </tr>
        """ %(file, xml[ID].ra, xml[ID].dec, select_mag[i],
              c.mag_f140w[i], c.z_peak[i], c.logm[i],
              root_path, file,
              root_path, file,
              root_path, file,
              root_path, file))
            
    fp.write("</tbody></table></body></html>")
    fp.close()
    fplist.close()
    fpreg.close()
    
    print 'N = %d' %NUSE

def dn4000():
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    import numpy as np
    import shutil
    
    cat = catIO.Readfile('massive.dat')
    use = (cat.lmass < 12) & (cat.lmass > 10.2) & (cat.z > 1.85) & (cat.z < 3.6)
    idx = np.arange(cat.N)[use]
    cat.dn = cat.z*0.-1
    
    for i in idx:
        PATH = unicorn.analysis.get_grism_path(root=cat.id[i].split('G141')[0]) 
        spec_path = PATH+'/HTML/ascii/%s.dat' %(cat.id[i])
        if os.path.exists(spec_path):
            #print spec_path
            sp = catIO.Readfile(spec_path)
            # yint = np.interp(lam, sp.lam/(1+cat.z[i]), sp.flux2)
            # plt.plot(sp.lam/(1+cat.z[i]), sp.flux2/yint, color='blue', alpha=0.01, marker='.', markersize=10, linestyle='none')
            lrest = sp.lam/(1+cat.z[i])
            blue = (lrest >= 3850) & (lrest <= 3950)
            red = (lrest >= 4000) & (lrest <= 4100)
            fnu = sp.flux2*lrest**2
            cat.dn[i] = np.trapz(fnu[blue],lrest[blue])/np.trapz(fnu[red],lrest[red])
    #
    gf = catIO.Readfile('massive.galfit')
    rcirc = gf.re*0.06*np.sqrt(gf.ba)    
    zgrid = np.array([0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
    scale = np.array([3.268, 4.421, 5.343, 5.733,  6.082, 6.394, 6.673, 6.922, 7.144, 7.518, 7.812, 8.041, 8.216, 8.346, 8.439, 8.502, 8.539, 8.556, 8.555, 8.540, 8.512, 8.475, 8.430, 8.377, 8.320, 8.257, 8.192])
    sint = np.interp(cat.z, zgrid, scale)
    rkpc = rcirc*sint
    
    mass = cat.lmass + np.random.randn(cat.lmass.size)*0.05
    plt.plot(mass[use], cat.dn[use], marker='.', color='red', linestyle='none', alpha=0.5, markersize=15)
    plt.xlim(10,12)
    plt.ylim(0,2)

    plt.plot(rkpc[use], cat.dn[use], marker='.', color='red', linestyle='none', alpha=0.5, markersize=15)
    plt.xlim(0,5)
    plt.ylim(0,2)
         
    return cat.dn
    
def open_spec(use):
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    import numpy as np
    import shutil
    
    cat = catIO.Readfile('massive.dat')
    gf = catIO.Readfile('massive.galfit')
    idx = np.arange(cat.N)[use]
    list = []
    for i in idx:
        PATH = unicorn.analysis.get_grism_path(root=cat.id[i].split('G141')[0]) 
        spec_path = PATH+'/HTML/SED/%s_SED.png' %(cat.id[i])
        galfit_file = cat.id[i]+'_galfit.png'
        if os.path.exists(spec_path) & os.path.exists(galfit_file):
            list.append(spec_path)
            list.append(galfit_file)
        
    os.system('open '+' '.join(list))
    
    
def chapman_smgs():
    import os
    import numpy as np
    import glob
    import threedhst
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/SCUBA/')
    lines = open('Chapman_table2.dat').readlines()
    N = len(lines)
    ra = np.arange(N)*1.; dec = ra*0.; zsmg=ra*0.; idsmg = []
    for i in range(N):
        rad = lines[i][5:23]
        ra[i] = (float(rad[0:2])+float(rad[2:4])/60.+float(rad[4:9])/3600.)/24.*360.
        dec[i] = float(rad[9:12])+float(rad[12:14])/60.+float(rad[14:])/3600.
        zsmg[i] = float(lines[i][61:66])
        idsmg.append(rad)
    #
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/SCUBA/')
    lines = open('Chapman10_tables1.dat').readlines()
    N = len(lines)
    ra = np.arange(N)*1.; dec = ra*0.; zsmg=ra*0.; idsmg = []
    for i in range(N):
        rad = lines[i]
        ra[i] = (float(rad[2:4])+float(rad[5:7])/60.+float(rad[8:14])/3600.)/24.*360.
        dec[i] = float(rad[15:18])+float(rad[19:21])/60.+float(rad[22:27])/3600.
        zsmg[i] = float(lines[i][37:43])
        idsmg.append(lines[i][28:35])
    
    #
    lines = open('../VLA/Morrison_table2.dat').readlines()
    N = len(lines)
    ra = np.arange(N)*1.; dec = ra*0.; zsmg=ra*0.; idsmg = []
    for i in range(N):
        rad = lines[i]
        ra[i] = (float(rad[6:8])+float(rad[9:11])/60.+float(rad[12:17])/3600.)/24.*360.
        dec[i] = float(rad[24:26])+float(rad[27:29])/60.+float(rad[30:34])/3600.
        zsmg[i] = 0.
        idsmg.append(lines[i][0:6])
        
    cats = glob.glob('../DATA/*.cat')
    list = []
    smg = []
    for cat in cats:
        ss = threedhst.sex.mySexCat(cat)
        cosf = np.cos(np.mean(ss.dec)/360*2*3.14159)
        root = os.path.basename(ss.filename).split('_drz')[0]
        #
        for i in range(N):
            dr = np.sqrt((ra[i]-ss.ra)**2*cosf**2+(dec[i]-ss.dec)**2)*3600
            mat = dr < 2
            if len(dr[mat]) > 0:
                for id in ss.id[mat]:
                    sed_file = '../HTML/SED/'+root+'_%05d' %(id)+'_SED.png'
                    if os.path.exists(sed_file):
                        list.append(sed_file)
                        smg.append((idsmg[i], zsmg[i]))
    
    os.system('open '+' '.join(list))
    for i in range(len(list)):
        print list[i], smg[i][0], smg[i][1]
        
def galfit_massive_galaxies():
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    import numpy as np
    import shutil
    import unicorn.analysis
    
    unicorn.analysis.show_massive_galaxies(contam=0.1, masslim=10.0, maglim=24, zrange=(0.4,5))
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/GALFIT/')
    shutil.copy('../massive.dat','.')
    
    cat = catIO.Readfile('massive.dat')
    use = cat.z > -100
    fit_sky=True
    force=False
    
    #### Refit
    gf = catIO.Readfile('massive.galfit')
    use = gf.n > 8
    fit_sky=False
    force=True
    
    for obj in cat.id[use]:
        if (not os.path.exists(obj+'_galfit.png')) | (force):
            threedhst.galfit.fit_3dhst_object(object=obj, fit_sky=fit_sky)
    
    fp = open('massive.galfit','w')
    fp.write('# id re ba n chi2\n')
        
    for obj in cat.id:
        if os.path.exists(obj+'_galfit.log'):
            log = threedhst.galfit.GalfitLogfile(obj+'_galfit.log')
            if 'psf' in log.components:
                fp.write('%s %6.2f %4.2f %6.2f %4.2f\n' %(obj, -1,-1,-1, np.log10(log.chi2)))
            else:
                c = log.list[0]
                fp.write('%s %6.2f %4.2f %6.2f %4.2f\n' %(obj, c.re.value,c.ba.value,c.n.value, np.log10(log.chi2)))
        else:
            fp.write('%s %6.2f %4.2f %6.2f %4.2f\n' %(obj, -2,-2,-2, -2))
    
    fp.close()
    
    #### Make a plot
    gf = catIO.Readfile('massive.galfit')
    rcirc = gf.re*0.06*np.sqrt(gf.ba)
    
    zgrid = np.array([0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])
    scale = np.array([3.268, 4.421, 5.343, 5.733,  6.082, 6.394, 6.673, 6.922, 7.144, 7.518, 7.812, 8.041, 8.216, 8.346, 8.439, 8.502, 8.539, 8.556, 8.555, 8.540, 8.512, 8.475, 8.430, 8.377, 8.320, 8.257, 8.192])
    sint = np.interp(cat.z, zgrid, scale)
    
    rkpc = rcirc*sint

    xvd = np.array([0, 0.6, 1.1, 1.6, 2.0])
    revd = np.array([12.4, 8.0, 5.3, 4.1, 3.0])
    revde = np.array([1.5, 0.8, 0.2, 0.3, 0.3])
    nvd = np.array([5.9, 4.0, 2.9, 2.5, 2.1])
    nvde = np.array([0.7, 0.4, 0.2, 0.2, 0.5])
    
    use = (cat.lmass >= 11) & (gf.chi2 > -1) & (gf.chi2 < 0.5)
    use = (cat.lmass >= 10.6) & (cat.lmass < 11) & (gf.chi2 > -1) & (gf.chi2 < 0.5)
    use = (cat.lmass >= 10.2) & (cat.lmass < 10.6) & (gf.chi2 > -1) & (gf.chi2 < 0.5)

    fig = plt.figure(figsize=[7,3],dpi=100)
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.1,
                        bottom=0.17,right=0.98,top=0.98)
    
    plt.plot(cat.z[use], rkpc[use], color='red', alpha=0.5, marker='.', linestyle='none', markersize=15)
    xm, ym, ys, N = threedhst.utils.runmed(cat.z[use], rkpc[use], NBIN=8)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(N), marker='.', color='orange', ecolor='orange', markersize=20, alpha=0.8)

    plt.errorbar(xvd, revd, yerr=revde, color='black', alpha=0.2, linewidth=4)
    
    plt.xlim(0,3)
    plt.ylim(0,15)
    
    plt.semilogy(np.arange(10)-100)
    plt.xlim(0,3)
    plt.ylim(0.8,15)
    
    plt.xlabel(r'$z$')
    plt.ylabel(r'$r_e^\mathrm{circ}$ [kpc]')

    plt.savefig('z_vs_re.pdf')
    plt.close()
        
    tiny = use & (rkpc < 2) & (cat.z < 1.5)
    plt.plot(cat.z[tiny], rkpc[tiny], color='blue', alpha=0.5, marker='.', linestyle='none', markersize=15)
    
    ##### Open subselections

    small = use & (cat.z > 1.5) & (cat.z < 2.4) & (rkpc < 2)
    plt.plot(cat.z[small], rkpc[small], color='blue', alpha=0.5, marker='.', linestyle='none', markersize=15)
    unicorn.analysis.open_spec(small)

    big = use & (cat.z > 1.5) & (cat.z < 2.4) & (rkpc > 4) & (rkpc < 20)
    plt.plot(cat.z[big], rkpc[big], color='blue', alpha=0.5, marker='.', linestyle='none', markersize=15)
    unicorn.analysis.open_spec(big)
    
    keep = (cat.z > 1.4) & (cat.z < 2.5) & (cat.lmass > 10.99)
    disk = (cat.z < 1.6) & (cat.z > 1.0) & (cat.lmass > 10.6) & (gf.n > 3)
    
    unicorn.analysis.open_spec(disk)
        
    ### z vs n
    fig = plt.figure(figsize=[7,3],dpi=100)
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.1,
                        bottom=0.17,right=0.98,top=0.98)
    
    plt.plot(cat.z[use], gf.n[use], color='red', alpha=0.5, marker='.', linestyle='none', markersize=15)
    xm, ym, ys, N = threedhst.utils.runmed(cat.z[use], gf.n[use], NBIN=4)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(N), marker='.', color='red', ecolor='red', markersize=30)
    plt.errorbar(xvd, nvd, nvde, color='black', alpha=0.2, linewidth=4)
    plt.ylim(0,10)
    plt.xlim(0,3)
    plt.ylabel(r'$n$')
    plt.xlabel(r'$z$')

    plt.savefig('z_vs_n.pdf')
    plt.close()
    
    #### z vs b/a
    fig = plt.figure(figsize=[7,3],dpi=100)
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.1,
                        bottom=0.17,right=0.98,top=0.98)
    
    plt.plot(cat.z[use], gf.ba[use], color='red', alpha=0.5, marker='.', linestyle='none', markersize=15)
    xm, ym, ys, N = threedhst.utils.runmed(cat.z[use], gf.ba[use], NBIN=4)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(N), marker='.', color='red', ecolor='red', markersize=30)
    plt.ylim(0,1)
    plt.xlim(0,3)
    plt.ylabel(r'$b/a$')
    plt.xlabel(r'$z$')
    
    plt.savefig('z_vs_ba.pdf')
    plt.close()
    
    #### mass-size
    mass = cat.lmass + np.random.randn(cat.lmass.size)*0.05
    use = (cat.z >= 0.9) & (cat.z < 1.5) & (gf.chi2 > -1) & (gf.chi2 < 0.5)
    plt.plot(mass[use], rkpc[use], color='red', alpha=0.5, marker='.', linestyle='none', markersize=10)

    use = (cat.z >= 1.5) & (cat.z < 2.0) & (gf.chi2 > -1) & (gf.chi2 < 0.5)
    plt.plot(mass[use], rkpc[use], color='orange', alpha=0.5, marker='.', linestyle='none', markersize=10)

    use = (cat.z >= 2.0) & (cat.z < 2.5) & (gf.chi2 > -1) & (gf.chi2 < 0.5)
    plt.plot(mass[use], rkpc[use], color='blue', alpha=0.5, marker='.', linestyle='none', markersize=10)

    plt.xlim(10.1,11.8)
    plt.ylim(0,15)
    
def testMatch():
    """
    Match the z=2.2 GOODS-N H-a emitters from Tadaki et al. 
    http://arxiv.org/abs/1012.4860v1
    """
    import unicorn.analysis
    
    ra_haem = [189.0141,188.9930,189.0906,189.2159,189.2582,189.2364,189.2204,189.2235,189.2152,189.2709,189.3259,189.3564]
    dec_haem = [62.1760,62.1988,62.2480,62.2513,62.2639,62.2788,62.2907,62.2901,62.3071,62.2908,62.2815,62.3194]
    
    #### LAE - HU et al.
    # ra_haem = [189.215240,189.216751,189.056106,189.254120,189.399750,189.324677,189.033004,189.456543,189.342285,189.045471,189.366013,189.320419]
    # dec_haem = [62.32683,62.36460,62.12994,62.35397,62.23944,62.29974,62.14394,62.22942,62.26277,62.17144,62.19613,62.23344]
    
    matches = []
    for i in range(len(ra_haem)):
        matches.extend(unicorn.analysis.matchObject(ra_haem[i], dec_haem[i]))
        print i, len(matches)
    
    fp = open('ha_emitters.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="http://localhost/~gbrammer/COSMOS/scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[1,1]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        3: {
                                sorter: false
                        },
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                        6: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> RA </th>
        <th> Dec </th>
        <th> Thumb </th>
        <th> 2D </th>
        <th> 1D </th>
        <th> SED </th>
    </thead>
    <tbody>
    """)
    
    for match in matches:
            fp.write("""
        <tr>
            <td> %s </td>
            <td> %f </td>
            <td> %f </td>
            <td> <img src=%s/images/%s_thumb.png height=180px> </td>
            <td> <img src=%s/images/%s_2D.png height=180px> </td>
            <td> <img src=%s/images/%s_1D.png height=180px> </td>
            <td> <img src=%s/SED/%s_SED.png height=180px> </td>
        </tr>
        """ %(match[1],
              match[3], match[4],
              match[0], match[1],
              match[0], match[1],
              match[0], match[1],
              match[0], match[1]))
    
    fp.write("</tbody></table></body></html>")
    fp.close()
        
    
def matchObject(ra_in, dec_in):
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS')
    
    catalogs = glob.glob('../GOODS-N/HTML_v1.0/*.cat')
    catalogs.extend(glob.glob('../COSMOS/HTML_v1.0/*.cat'))
    catalogs.extend(glob.glob('../ERS/HTML/*.cat'))
    catalogs.extend(glob.glob('../GOODS-S-SN/HTML/*.cat'))
    
    out = []
    for catalog in catalogs:
        
        cat = threedhst.sex.mySexCat(catalog)
        dir_root = os.path.dirname(catalog) #.split('/')[1]
        cat_root = os.path.basename(catalog).split('_drz')[0]
        
        dr = np.sqrt((cat.ra-ra_in)**2*np.cos(dec_in/360.*2*np.pi)**2+(cat.dec-dec_in)**2)*3600.
        idx = np.where(dr == np.min(dr))[0][0]
        drMatch = dr[idx]*1.
        if drMatch < 1:
            out.append([dir_root, '%s_%05d' %(cat_root, np.int(cat.NUMBER[idx])),
                        drMatch, ra_in, dec_in])
    
    return out

def test_line_histogram():
    
    import threedhst
    import matplotlib.pyplot as plt
    
    f = threedhst.spec1d.readLinesDat('HTML_v1.0/images/orient1_1D_lines.dat')
    q = (f.sn > 1.3) & (f.sn < 20) & (f.ndet < 5) & (f.sigma > 20) & (f.sigma < 100)
    
    q2 = q & (f.wave > 1.45e4) & (f.wave < 1.465e4)
    

def check_masses():
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    
    #### COSMOS
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    cat_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    cat_cosmos.kmag = 25.-2.5*np.log10(cat_cosmos.field('ktot'))
    zout_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    fout_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.bc03.v4.6.fout')
    
    #### GOODS-N
    MAIN_OUTPUT_FILE = 'photz'
    OUTPUT_DIRECTORY = unicorn.GRISM_HOME+'GOODS-N/MODS/EAZY/OUTPUT/'
    
    MODS = unicorn.GRISM_HOME+'GOODS-N/MODS/'
    cat_goodsn = catIO.Readfile(MODS+'/FAST/mods.cat')
    cat_goodsn.addColumn('kmag',25-2.5*np.log10(cat_goodsn.ks_ap))
    zout_goodsn = catIO.Readfile(MODS+'/EAZY/OUTPUT/mods.zout')
    fout_goodsn = catIO.ReadASCIICat(MODS+'/FAST/mods.bc03.fout')
    
    #### GOODS-S
    FIREWORKS = unicorn.GRISM_HOME+'GOODS-S-SN/FIREWORKS/'
    cat_fw = unicorn.analysis.catIO.ReadASCIICat(FIREWORKS+'/FIREWORKS_phot.cat')        
    cat_fw.kmag = 23.86-2.5*np.log10(cat_fw.field('Ks_totf'))
    zout_fw = unicorn.analysis.catIO.ReadASCIICat(FIREWORKS+'/fireworks.zout')
    fout_fw = unicorn.analysis.catIO.ReadASCIICat(FIREWORKS+'/fireworks.fout')
    
    #### Some K limit
    K_limit = 22.8
    q_cosmos = (cat_cosmos.kmag < K_limit) & (cat_cosmos.star_flag == 0) & (cat_cosmos.nchild <= 0)
    
    q_goodsn = cat_goodsn.kmag < K_limit
    q_fw = cat_fw.kmag < K_limit
    
    plt.plot(fout_cosmos.z[q_cosmos], fout_cosmos.lmass[q_cosmos], color='red', alpha=0.1, linestyle='None', marker='o')

    plt.plot(fout_goodsn.z[q_goodsn], fout_goodsn.lmass[q_goodsn], color='blue', alpha=0.3, linestyle='None', marker='o')

    plt.plot(fout_fw.z[q_fw], fout_fw.lmass[q_fw], color='green', alpha=0.3, linestyle='None', marker='o')
    
    plt.ylim(6,13)
    plt.xlim(0,4)

def go_field_area():
    import unicorn.analysis
    import threedhst.catIO as catIO
    
    #### COSMOS
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    MAIN_OUTPUT_FILE = 'aegis-n2.v4.6'
    cat_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    q = cat_cosmos.wmin > 0.3
    area = unicorn.analysis.field_area(cat_cosmos.ra[q], cat_cosmos.dec[q])
    print area
    
def field_area(ra_in, dec_in):
    from shapely.geometry import MultiPoint, Point, MultiPolygon
    import numpy as np
    import matplotlib.pyplot as plt
    from descartes import PolygonPatch
    from shapely.ops import cascaded_union
    
    ra = ra_in*np.cos(np.mean(dec_in)/360.*2*np.pi)
    
    coords = []
    for i in xrange(len(ra)):
        coords.append((ra[i], dec_in[i]))
    
    env = MultiPoint(coords).envelope 
    env_area = env.area
    
    hull = MultiPoint(coords).convex_hull 
    hull_area = hull.area
    
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    
    ax.plot(ra, dec_in, marker='.', linestyle='None', alpha=0.2, color='blue')
    x_e, y_e = env.boundary.xy
    ax.plot(x_e, y_e, color='red')

    x_h, y_h = hull.boundary.xy
    ax.plot(x_h, y_h, color='yellow')
        
    # plt.show()
    
    return hull_area

def LAE_morphologies():
    """
    Have a catalog of GOODS-S LAEs from Lucia via Harold.  Cut out thumbnails from my CANDELS mosaic.
    """
    import pywcs
    import matplotlib.pyplot as plt
    
    os.chdir("/research/HST/GRISM/3DHST/ANALYSIS/LAE")
    print 'Reading images...'
    f125w = pyfits.open('/research/HST/CANDELS/GOODS-S/PREP_FLT/goods-s_f125w_drz.fits.gz')
    f160w = pyfits.open('/research/HST/CANDELS/GOODS-S/PREP_FLT/goods-s_f160w_drz.fits')
    print 'Done.'
    
    lae = catIO.Readfile('GOODS-S/LuciaALLsample_34red_NB.txt')
    
    wcs = pywcs.WCS(header=f125w[1].header)
    
    xy = wcs.wcs_sky2pix(lae.ra, lae.dec, 1)
    NSUB = 16
    i = np.arange(lae.N)[lae.n20 == 135][0]
    
    plt.gray()
    
    for i in range(lae.N):
        if (xy[0][i] > NSUB) & (xy[0][i] < (wcs.naxis1-NSUB)) & (xy[1][i] > NSUB) & (xy[1][i] < (wcs.naxis2-NSUB)):
             
            print 'LAE_N20-%d.png' %lae.n20[i]
                       
            xi = xy[0][i]
            yi = xy[1][i]
            sub125 = f125w[1].data[yi-NSUB:yi+NSUB,xi-NSUB:xi+NSUB]
            sub160 = f160w[1].data[yi-NSUB:yi+NSUB,xi-NSUB:xi+NSUB]
            
            if np.std(sub125) == 0:
                continue
                
            fig = plt.figure(figsize=[6,3.4],dpi=100)
            fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.02,
                                bottom=0.02,right=0.98,top=0.86)
            
            scale = 0.128254
            size=NSUB
            asec_pix = 1./scale
            nasec = int(size/asec_pix)

            ax = fig.add_subplot(121)
            ax.imshow(sub125, vmin=-0.05,vmax=0.1, interpolation='nearest')
            
            ax.set_yticklabels([])
            xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            ax.set_xticklabels([])
            ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            
            ax.set_title('F125W')
            
            ax = fig.add_subplot(122)
            ax.imshow(sub160, vmin=-0.05,vmax=0.1, interpolation='nearest')

            ax.set_yticklabels([])
            xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            ax.set_xticklabels([])
            ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            
            ax.set_title('F160W')
            
            ax.set_text(-0.01, 1.12,'LAE N20-%d' %lae.n20[i],
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = ax.transAxes, fontsize=16)
            
            plt.savefig('LAE_N20-%d.png' %lae.n20[i])
            plt.close()

def fit_all_brown_dwarfs():
    import glob
    import unicorn
    
    bd = unicorn.analysis.BD_fit()
    
    files = glob.glob('/Users/gbrammer/Sites_GLOBAL/P/GRISM/ascii/*-G141*.dat')
    for file in files:
        bd.fit(file)
        
class BD_template():
    def __init__(self, txt):
        self.filename = txt
        self.read_template(txt)
    
    def read_template(self, txt):
        import numpy as np
        lines = open(txt).readlines()
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
    def __init__(self,template_path='./STANDARDS'):
        self.template_path=template_path
        self.read_bd_templates()
        
    def read_bd_templates(self):
        import glob
        temps = glob.glob(self.template_path+'/spex*txt')
        list = []
        for temp in temps:
            list.append(BD_template(temp))
        
        self.templates = list
        self.NTEMP = len(self.templates)
        
    def fit(self, ascii_file='AEGIS-3-G141_00177.dat'):
        import threedhst.catIO as catIO
        import numpy as np
        
        spec = catIO.Readfile(ascii_file)
        spec.error /= 2.5
        
        self.spec = spec
        
        chi2 = np.zeros(self.NTEMP)
        types = []
        anorm = chi2*0.

        use = (spec.lam > 1.1e4) & (spec.lam < 1.65e4) & (spec.contam/spec.flux < 0.05) & (np.isfinite(spec.flux)) & (spec.flux > 0)
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
        
        noNewLine = '\x1b[1A\x1b[1M'
        
        min = np.where(chi2 == chi2.min())[0][0]
        if (chi2.min() < 1.5) & (~types[min].startswith('M')):
            print noNewLine + ascii_file+' * '
            self.spec = spec
            self.types = types
            self.chi2 = chi2
            self.anorm = anorm
            self.ascii_file = ascii_file
            self.make_plot()
        else:
            print noNewLine + ascii_file + ' %0.2f' %(chi2.min())
            
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
        
        ax.text(1.1e4, 1.07, os.path.basename(self.ascii_file.split('.dat')[0]))
        
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
        
        outfile = os.path.basename(self.ascii_file.replace('.dat','_BD.png'))
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)

#

def eazy_unicorn():
    import unicorn
    unicorn.analysis.run_eazy_fit(root='COSMOS-3-G141', id=874, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='AEGIS-3-G141', id=683, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-25-G141', id=597, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')
    
    unicorn.analysis.run_eazy_fit(root='GOODS-N-33-G141', id=966, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-33-G141', id=1033, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-18-G141', id=1043, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')
    
    unicorn.analysis.run_eazy_fit(root='GOODS-N-23-G141', id=949, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-13-G141', id=535, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-15-G141', id=1139, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-N-14-G141', id=848, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='MARSHALL-225-G141', id=661, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')
    
    unicorn.analysis.run_eazy_fit(root='GOODS-S-24-G141', id=39, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')

    unicorn.analysis.run_eazy_fit(root='GOODS-S-24-G141', id=222, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')
    
    unicorn.analysis.run_eazy_fit(root='AEGIS-2-G141', id=100, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest')
    
    
    
def make_eazy_inputs(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', check=False, bin_spec=1, spec_norm=1.):
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    ORIG_PATH = os.getcwd()
    
    PATH = unicorn.analysis.get_grism_path(root)
    #print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=root)
        
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    grismCat, SPC = unicorn.analysis.read_grism_files(root=root)
    
    match = threedhst.catIO.Readfile('HTML/SED/'+root+'_match.cat')
    
    lam, spflux, sperr, lci, fobs, efobs, photom_idx = unicorn.analysis.specphot(id=id,
        grism_root=root, SPC = SPC, 
        cat = cat,
        grismCat = grismCat, zout = zout, fout = fout, 
        OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
        MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
        OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
        CACHE_FILE = 'Same',
        GET_SPEC_ONLY = True)
        
    use = (lam > 1.1e4) & (lam < 1.65e4) & (spflux != 0.0) & np.isfinite(spflux) & np.isfinite(sperr)
    
    #### allow additional normalization term
    spflux *= spec_norm
    sperr *= spec_norm
    
    ## check
    if check:
        plt.semilogx([1],[1])
        plt.errorbar(lci, fobs, yerr=efobs, marker='o', linestyle='None', color='blue', ecolor='blue')
        ma = fobs.max()
        plt.plot(lam[use], spflux[use], color='red')
        plt.ylim(-0.1*ma, 1.5*ma)
        plt.xlim(3000, 9.e4)
    ##
    
    #### Rebin if bin_spec is specified
    bin_spec = int(bin_spec)
    if bin_spec > 1:
        NUSE = len(lam[use])
        lbin = np.zeros(NUSE/bin_spec)
        fbin = np.zeros(NUSE/bin_spec)
        ebin = np.zeros(NUSE/bin_spec)
        for i in range(0,NUSE,bin_spec):
            try:
                lbin[i/bin_spec] = np.mean(lam[use][i:i+bin_spec])
                ebin[i/bin_spec] = 1./np.sum(1./sperr[use][i:i+bin_spec]**2)
                fbin[i/bin_spec] = np.sum(spflux[use][i:i+bin_spec]/sperr[use][i:i+bin_spec]**2)*ebin[i/bin_spec]
            except:
                pass

        use = lbin > 0
        lam = lbin; spflux = fbin; sperr = ebin
    
    #### Take a wavelength bin and convolve it with the thumbnail
    dlam = (lam[use][1]-lam[use][0])*bin_spec
    NBIN = 100
    DX = 100
    xarr = (np.arange(NBIN*2)-NBIN)*1./NBIN*3*dlam
    yarr = xarr*0
    yarr[np.abs(xarr) < dlam/2.] = 1
    
    os.chdir(PATH)
    lsm, ysm = unicorn.analysis.convolveWithThumb(id, xarr+lam[use][0], yarr, SPC, verbose=False)
    #### Recenter the convolved response
    xoff = lam[use][0]-np.trapz(lsm, ysm*lsm)/np.trapz(lsm, ysm)
    lsm, ysm = unicorn.analysis.convolveWithThumb(id, xarr+lam[use][0]+xoff, yarr, SPC, verbose=False)
    keep = ysm > (1.e-4*np.max(ysm))

    ###################################
    #### Make a RES file for the convolved resolution and make a catalog file with
    #### photometry and spectroscopy
    
    os.chdir(ORIG_PATH)
    
    res_lines = open(OLD_RES).readlines()
    nfilt=0
    res_name = []
    res_lc = []
    for line in res_lines:
        if line.find('lambda_c') > 0:
            res_lc.append(line.split('lambda_c=')[1].split()[0])
            res_name.append(line.split()[1])
            nfilt += 1
    #
    res_lc = np.cast[float](res_lc)
    
    #### Read the parameter file associated with the zout file to get the right
    #### filters, which should be in the same order as the vectors read from the 
    #### binary files
    eazy_param = threedhst.eazyPy.EazyParam(zout.filename.replace('zout','param'))
    
    cat_head = '# phot_id id ra dec z_spec'
    cat_line = ' %d %s_%05d %f %f %8.3f' %(cat.id[photom_idx], root, id, cat.ra[photom_idx], cat.dec[photom_idx], zout.z_spec[photom_idx])
    no_spec = cat_line+''
    no_phot = cat_line+''
    for i in range(len(lci)):
        fnumber = eazy_param.filters[i].fnumber
        #print res_lc[fnumber-1], lci[i]
        cat_head += ' F%0d E%0d' %(fnumber, fnumber)
        fnui = (lci[i]/5500.)**2
        cat_line += ' %9.3e %9.3e' %(fobs[i]*fnui, efobs[i]*fnui)
        no_spec += ' %9.3e %9.3e' %(fobs[i]*fnui, efobs[i]*fnui)
        no_phot += ' -99 -99' 
    
    NKEEP = len(lsm[keep])    
    fnu_spec = spflux * (lam/5500.)**2
    efnu_spec = sperr * (lam/5500.)**2
    
    for i, l0 in enumerate(lam[use]):
        ### CAT
        cat_head += ' F%0d E%0d' %(nfilt+i+1, nfilt+i+1)
        cat_line += ' %9.3e %9.3e' %(fnu_spec[use][i], efnu_spec[use][i]/3)
        no_spec += ' -99  -99'
        no_phot += ' %9.3e %9.3e' %(fnu_spec[use][i], efnu_spec[use][i]/3)
        #
        ### RES
        lsm0 = lsm[keep]-lam[use][0]+l0
        ysmk = ysm[keep]
        res_lines.append('%d spec lambda_c=%f\n' %(NKEEP, l0))
        for j in range(NKEEP):
            res_lines.append(' %5i %13.5e %13.5e\n' %(j+1, lsm0[j], ysmk[j]))
    
    fp = open(OUT_RES,'w')
    fp.writelines(res_lines)
    fp.close()
    
    fp = open('threedhst.cat','w')
    fp.write(cat_head+'\n')
    fp.write(cat_line+'\n')
    fp.write(no_spec+'\n')
    fp.write(no_phot+'\n')
    fp.close()
    
    #### Set some EAZY parameters
    eazy_param.params['READ_ZBIN'] = 'n'
    eazy_param.params['TEMPLATES_FILE'] = 'templates/eazy_v1.1_lines.spectra.param'
    eazy_param.params['FILTERS_RES'] = 'THREEDHST.RES'
    eazy_param.params['OUTPUT_DIRECTORY'] = 'OUTPUT'
    eazy_param.params['WAVELENGTH_FILE'] = 'templates/EAZY_v1.1_lines/lambda_v1.1.def'
    eazy_param.params['CATALOG_FILE'] = 'threedhst.cat'
    eazy_param.params['MAIN_OUTPUT_FILE'] = 'threedhst'
    eazy_param.params['CACHE_FILE'] = 'threedhst.tempfilt'
    eazy_param.params['GET_ZP_OFFSETS'] = 0
    eazy_param.params['Z_STEP'] /= 4
    eazy_param.params['Z_MIN'] = zout.l99[photom_idx]*0.8
    eazy_param.params['Z_MAX'] = zout.u99[photom_idx]/0.8
    eazy_param.params['MAGNITUDES'] = 0.0
    eazy_param.write(file='threedhst.eazy.param')
    
    
def run_eazy_fit(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy'):
    
    import matplotlib.pyplot as plt
    import threedhst.eazyPy as eazy
    
    if run:
        unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm)
        os.system(eazy_binary + ' -p threedhst.eazy.param '+pipe)
    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(0, MAIN_OUTPUT_FILE='threedhst', \
                          OUTPUT_DIRECTORY='OUTPUT', \
                          CACHE_FILE = 'Same')
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    ################
    ### Make plot
    ################
    
    ### setup
    plt.rcParams['font.size'] = 10
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[7,3],dpi=100)
    else:
        fig = Figure(figsize=[7,3],dpi=100)
        
    fig.subplots_adjust(wspace=0.04,hspace=0.12,left=0.05,
                        bottom=0.15,right=0.98,top=0.98)
    
    #################################### Broad-band + spectrum
    ax = fig.add_subplot(131)
    
    ax.semilogx([1],[1])
    ymax = max(fobs)
    ymax = max(fobs[is_spec])
    
    ## photometry
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.9, color='orange', markersize=6)
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='red', linestyle='None', markersize=6, alpha=0.2)
    
    ## best-fit SED
    ax.plot(lambdaz, temp_sed, color='red', alpha=0.6)
    
    ## Spectrum + convolved fit
    ax.plot(lci[is_spec], obs_sed[is_spec], color='red', markersize=6, alpha=0.7, linewidth=1)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.8, color='blue', linewidth=1)
    
    ## check normalization
    spec_norm = np.sum(obs_sed[is_spec]*fobs[is_spec]/efobs[is_spec]**2)/np.sum(obs_sed[is_spec]**2/efobs[is_spec]**2)
    print spec_norm
    
    ax.set_yticklabels([])
    ax.set_ylabel(r'$f_\lambda$')
    ax.set_xlabel(r'$\lambda$')
    xtick = ax.set_xticks(np.array([0.5, 1., 2, 4])*1.e4)
    ax.set_xticklabels(np.array([0.5, 1., 2, 4]))
    
    ax.set_xlim(3000,9.e4)
    ax.set_ylim(-0.1*ymax, 1.2*ymax)
    
    ymax = max(fobs[is_spec])
    ymin = min(fobs[is_spec])
    
    ################################## Spectrum only
    ax = fig.add_subplot(132)
    ## photometry
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.9, color='orange', markersize=10)
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='red', linestyle='None', markersize=10, alpha=0.2)
    
    ## best-fit template
    ax.plot(lambdaz, temp_sed, color='red', alpha=0.2)
    
    ## spectrum
    ax.plot(lci[is_spec], obs_sed[is_spec], color='red', markersize=6, alpha=0.7, linewidth=2)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.6, color='blue', linewidth=2)
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.array([1.,1.2, 1.4, 1.6])*1.e4)
    ax.set_xticklabels(np.array([1.,1.2, 1.4, 1.6]))
    ax.set_xlabel(r'$\lambda$')
    
    ax.set_xlim(0.95e4,1.78e4)
    ax.set_ylim(0.8*ymin, ymax*1.1)
    
    #################################### p(z) for combined, photometry, and spec 
    ax = fig.add_subplot(133)
    
    colors = ['purple','orange','blue']
    alpha = [0.5, 0.5, 0.2]
    zo = threedhst.catIO.Readfile('OUTPUT/threedhst.zout')
    zmin = 4
    zmax = 0
    ymax = 0
    for i in range(3):
        zgrid, pz = eazy.getEazyPz(i, MAIN_OUTPUT_FILE='threedhst', 
                          OUTPUT_DIRECTORY='./OUTPUT', 
                          CACHE_FILE='Same')
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        #
        if pz.max() > ymax:
            ymax = pz.max()
        #
        if zgrid[pz > 1.e-6].min() < zmin:
            zmin = zgrid[pz > 1.e-6].min()
        #
        if zgrid[pz > 1.e-6].max() > zmax:
            zmax = zgrid[pz > 1.e-6].max()
    
    ax.plot(zo.z_spec[0]*np.array([1,1]),[0,1.e4], color='green', linewidth=1)
    
    
    ax.set_yticklabels([])
    ax.set_xlabel(r'$z$')
    
    ### Plot labels
    ax.text(0.5, 0.9, root+'_%05d' %(id), transform = ax.transAxes, horizontalalignment='center')
    ax.text(0.95, 0.8, r'$z_\mathrm{phot}=$'+'%5.3f' %(zo.z_peak[1]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
    ax.text(0.95, 0.7, r'$z_\mathrm{gris}=$'+'%5.3f' %(zo.z_peak[0]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
    if zo.z_spec[0] > 0:
        ax.text(0.95, 0.6, r'$z_\mathrm{spec}=$'+'%5.3f' %(zo.z_spec[0]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
        
    ax.set_xlim(zmin-0.2, zmax+0.2)
    ax.set_xlim(zgrid.min(), zgrid.max())
    ax.set_ylim(0,1.1*ymax)
    
    #### Save an image
    outfile = '%s_%05d_eazy.png' %(root, id)
    if USE_PLOT_GUI:
        plt.savefig(outfile)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)
    
    print outfile
    
def find_brown_dwarfs():
    import glob
    import threedhst.catIO as catIO
    
    obj = catIO.Readfile('AEGIS-3-G141_00177.dat')
    
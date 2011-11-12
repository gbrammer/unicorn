import os
import pyfits
import numpy as np
import glob
import shutil
import time

import matplotlib.pyplot as plt

# testing by Britt
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

HAS_PHOTOMETRY = True
PHOTOMETRY_ID = None
BAD_SPECTRUM = False

SPC_FILENAME = None
SPC = None

noNewLine = '\x1b[1A\x1b[1M'
 
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
    if root.startswith('UDS'):
        PATH = unicorn.GRISM_HOME+'UDS/'
    if root.startswith('GN20'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('G850.1'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('GOODS-S'):
        PATH = unicorn.GRISM_HOME+'GOODS-S/'
    if root.startswith('UDF'):
        PATH = unicorn.GRISM_HOME+'UDF/'
    if root.startswith('WFC3-ERS'):
        PATH = unicorn.GRISM_HOME+'ERS/'
    if root.startswith('MARSHALL'):
        PATH = unicorn.GRISM_HOME+'SN-MARSHALL/'
    if root.startswith('PRIMO'):
        PATH = unicorn.GRISM_HOME+'SN-PRIMO/'
    if root.startswith('GEORGE'):
        PATH = unicorn.GRISM_HOME+'SN-GEORGE/'
    #
    return PATH
    
def read_catalogs(root='', cosmos=False, aegis=False, goodsn=False, cdfs=False, ecdfs=False, uds=False, udf=False):
    """
    
    Read photometry, redshift, SPS catalogs for a given field.
    
    Usage: cat, zout, fout = read_catalogs(cosmos=True)
    
    """        
    import threedhst.catIO as catIO
    import numpy as np
    
    KTOT_COL = None
    ### Optionally supply a grism basename and determine the field from that
    if root.startswith('UDS'):
        uds=True
    if root.startswith('COSMOS'):
        cosmos=True
    if root.startswith('AEGIS'):
        aegis=True
    if root.startswith('GOODS-N'):
        goodsn=True
    if root.startswith('GN20') | root.startswith('G850.1'):
        goodsn=True
    if root.startswith('GOODS-S'):
        cdfs=True

    if root.startswith('GOODS-S-34') | root.startswith('GOODS-S-36') | root.startswith('GOODS-S-37') | root.startswith('GOODS-S-38') | root.startswith('UDF'):
        udf=True
        
    if root.startswith('MARSHALL'):
        uds=True
    if root.startswith('PRIMO'):
        cdfs=True
    if root.startswith('GEORGE'):
        cdfs=True
    if root.startswith('WFC3-ERS'):
        cdfs=True
    
    aegis_wirds=False
    if root.startswith('AEGIS-11') | root.startswith('AEGIS-2-') | root.startswith('AEGIS-1-') | root.startswith('AEGIS-6-'):
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
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
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
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            CAT_PATH = '/3DHST/Ancillary/AEGIS/NMBS/Photometry/'
            CAT_FILE = CAT_PATH + 'aegis-n2.deblend.v5.1.cat'
            ZOUT_FILE = CAT_PATH + '/aegis-n2.deblend.redshifts/aegis-n2.deblend.v5.1.zout'
            FOUT_FILE = CAT_PATH+'aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout'
            KTOT_COL = 'K'
    
    if aegis_wirds:
        GRISM_PATH=unicorn.GRISM_HOME+'AEGIS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/WIRDS/FAST/'
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
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
    
    # ##### Use Ros' new catalog (v1.7)
    # if goodsn:
    #     GRISM_PATH=unicorn.GRISM_HOME+'GOODS-N/'
    #     CAT_PATH = '/research/HST/GRISM/3DHST/GOODS-N/MODS/FAST/'
    #     if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
    #         CAT_PATH = '/3DHST/Photometry/Release/GOODS-N/v1.7/'
    #     #
    #     CAT_FILE = CAT_PATH+'goodsn_v1.7.fullz.cat'
    #     ZOUT_FILE = CAT_PATH+'goodsn_v1.7_eazy/photz_v1.7.fullz.zout'
    #     FOUT_FILE = CAT_PATH+'FAST/v1.7/goodsn_v1.7.fullz_ed.fout'
    #     KTOT_COL = 'f_Ks'
    
    # #### Force UDF to use FIREWORKS
    # if udf:
    #     cdfs = True
    #     udf = False
          
    if cdfs:
        GRISM_PATH=unicorn.GRISM_HOME+'GOODS-S/'
        CAT_PATH = GRISM_PATH+'FIREWORKS/'
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            CAT_PATH = '/3DHST/Ancillary/GOODS-S/FIREWORKS/FAST/'
        #
        CAT_FILE = CAT_PATH+'fireworks.cat'
        ZOUT_FILE = CAT_PATH+'fireworks.zout'
        FOUT_FILE = CAT_PATH + 'fireworks.fout'
        KTOT_COL = 'Ks_totf'
    #
    if udf and (unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp')):
        GRISM_PATH=unicorn.GRISM_HOME+'UDF/'
        CAT_PATH = '/3DHST/Ancillary/GOODS-S/HUDF09/CAT/'

        CAT_FILE = CAT_PATH+'hudf09.cat'
        ZOUT_FILE = CAT_PATH+'hudf09_noU.zout'
        FOUT_FILE = CAT_PATH + 'hudf09.fout'
        KTOT_COL = 'f_kv2_tot'
       
    if uds:
        GRISM_PATH=unicorn.GRISM_HOME+'SN-MARSHALL/'
        CAT_PATH = GRISM_PATH+'UDSPHOT/'
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
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
    import unicorn.analysis
    
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
    if root+'_2_opt.SPC.fits' == unicorn.analysis.SPC_FILENAME:
        SPC = unicorn.analysis.SPC
    else:
        try:
            try:
                unicorn.analysis.SPC.fits.close()
            except:
                pass
            SPC = threedhst.plotting.SPCFile(root+'_2_opt.SPC.fits',
                axe_drizzle_dir=BASE_PATH+'DRIZZLE_'+GRISM_NAME)
            unicorn.analysis.SPC = SPC
            unicorn.analysis.SPC_FILENAME = SPC.filename
        except:
            SPC = None
            unicorn.analysis.SPC_FILENAME = None
            unicorn.analysis.SPC = None
    
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
    
    profile = profile[::-1]
    
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
        if drMatch > 1:
            return False
        else:
            return lam, ffix*anorm, np.sqrt(ferr**2+(1.0*spec.field('CONTAM'))**2)*anorm, lci, fobs, efobs, photom_idx
         
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
    
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.13*4.8/xs, bottom=0.15*4.8/xs,right=1.-0.02*4.8/xs,top=1-0.10*4.8/xs)
    
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
    
    print noNewLine+OUT_PATH+'/'+out_file
    
    if Verbose:
        print 'Close the plot window'
        
    
def match_grism_to_phot(grism_root='ibhm45', MAIN_OUTPUT_FILE = 'cosmos-1.v4.6', OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/', CACHE_FILE = 'Same', Verbose=False, SPC = None, cat=None, grismCat = None, zout = None, fout = None, OUTPUT='/tmp/match.cat'):

    """
    Generate the matched catalogs of 3D-HST objects and objects in ancillary catalogs
    """
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
    
    (2011-10-11, this is old code where I tried to make a fake fluxcube with observed
    colors from photometric catalogs.  Never really used and deprecated with the
    availability of CANDELS imaging.)
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
    use_kmag=False, contam=0.5, coverage=0.9, skip_goodsn=False):        
    """
    Make a webpage showing objects selected on mass, redshift, contamination, mag...
    """
    
    if unicorn.hostname().startswith('unicorn'):
        os.chdir('/Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
        scripts="../scripts"
        files = glob.glob('../SED/*match.cat')
        matches = []
        for file in files:
            if skip_goodsn:
                if not 'GOODS-N' in file:
                    matches.append(file)
            else:
                matches.append(file)
    else:
        os.chdir(unicorn.GRISM_HOME+'ANALYSIS')
        scripts="http://localhost/~gbrammer/COSMOS/scripts"

        matches = []
        matches.extend(glob.glob('../AEGIS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../COSMOS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../GOODS-S/HTML/SED/*match.cat'))
        if not skip_goodsn:
            matches.extend(glob.glob('../GOODS-N/HTML/SED/*match.cat'))
        
        # matches.extend(glob.glob('../SN-GEORGE/HTML/SED/*match.cat'))
        # matches.extend(glob.glob('../SN-PRIMO/HTML/SED/*match.cat'))    
        # matches.extend(glob.glob('../SN-MARSHALL/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../ERS/HTML/SED/*match.cat'))
    
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
    <p>
        mag_F140W < %.2f
        <br>log M/Msun > %.2f
        <br>%.2f < z < %.2f
        <br> contam(1.4mum) < %.2f
        <br> coverage > %.2f
    </p>
    
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
    """ %(scripts, scripts, scripts,
          maglim, masslim, zrange[0], zrange[1], contam, coverage))
    
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
        
def field_area(ra_in, dec_in):
    """
    Given a list of ra / dec coordinates, compute the field area covered by those
    objects.
    """
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


def make_full_selection(zmin=None, zmax=None):
    import unicorn
    
    ########## Full selections to get everything
    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=23.5, zrange=(0.2,3.5),  use_kmag=False, contam=0.1, coverage=0.95)
    out='full_faint.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=22, zrange=(0.2,3.5),  use_kmag=False, contam=0.5, coverage=0.8)
    out='full_bright.html'
    unicorn.analysis.run_eazy_products_on_html(out)

    unicorn.analysis.show_massive_galaxies(masslim=10.49, maglim=24, zrange=(0.2,3.5),  use_kmag=False, contam=0.5, coverage=0.8)
    out='full_massive.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    ########## Bright galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=21., zrange=(0.7,2.8),  use_kmag=False, contam=0.05, coverage=0.9)
    out='full_z_F140W_lt_21.html'

    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=22., zrange=(1.0,1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='full_z_F140W_lt_22.html'
    
    ########## Massive galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=22., zrange=(0.7,2.8),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_11_F140W_22.html'

    ########## Massive z>1 galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=24., zrange=(1, 1.5),  use_kmag=False, contam=0.2, coverage=0.9)
    out='mass_11_z_1_1.5.html'
    unicorn.analysis.run_eazy_products_on_html(out)

    unicorn.analysis.show_massive_galaxies(masslim=10.49, maglim=23., zrange=(1, 1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_10.5_z_1_1.5.html'

    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=23., zrange=(1.5, 2.5),  use_kmag=False, contam=0.2, coverage=0.9)
    out='mass_11_z_1.5_2.5.html'
    
    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=22., zrange=(1.0,1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_10_mag_22_z_1.0_1.5.html'

    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=21.5, zrange=(0.8,1.4),  use_kmag=False, contam=0.05, coverage=0.95, skip_goodsn=False)
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=22.0, zrange=(0.8,1.5),  use_kmag=False, contam=0.05, coverage=0.95, skip_goodsn=False)
    out='test.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    unicorn.analysis.show_massive_galaxies(masslim=9, maglim=22, zrange=(1.1,4),  use_kmag=False, contam=0.1, coverage=0.8)
    out='test.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
def run_eazy_products_on_html(out):
    import unicorn
    import os
    import shutil
    
    if out is not 'massive.html':
        shutil.move('massive.html', out)
        
    shutil.copy('/Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS/'+out, unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    #out='full_faint.html'
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    
    unicorn.analysis.process_eazy_redshifts(html=out, zmin=0.2, zmax=3.8, compress=0.8)
    unicorn.analysis.show_triple_zphot_zspec(zout=out.replace('html','zout'), zmin=0, zmax=2.5)
    
    #status = os.system('rsync -avz *.png *.html *.pdf ~/Sites_GLOBAL/P/GRISM_v1.6/EAZY/')
    status = os.system('rsync -avz *zz.png *.html *.pdf ~/Sites_GLOBAL/P/GRISM_v1.6/EAZY/')
    
def process_eazy_redshifts(html='massive.html', zmin=None, zmax=None, compress=1.0):
    import unicorn
    
    fp = open(html)
    lines = fp.readlines()
    fp.close()
    
    # example_zout = glob.glob('OUTPUT/*.zout')[0]
    # os.system('head -1 ' + example_zout + ' > '+html.replace('html','zout'))
    fp = open(html.replace('html','zout'),'w')
    fp.write('# id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68  l95 u95  l99 u99  nfilt q_z z_peak peak_prob z_mc\n')
    fp.close()
    
    for i,line in enumerate(lines):
        failed = False
        
        if 'G141' in line and 'png' not in line:
            object = line.split('<br>')[0].split('<td>')[1].strip()
            root=object.split('G141_')[0]+'G141'
            id = int(object.split('G141_')[1])
            
            if not os.path.exists('OUTPUT/%s.zout' %(object)):
                try:
                    print 'MISSING! '+object
                    #unicorn.analysis.run_eazy_fit(root=root, id=id, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest', zmin=zmin, zmax=zmax, compress=compress, COMPUTE_TILT=True, TILT_ORDER=1)
                    #os.system('cat OUTPUT/threedhst.zout > OUTPUT/%s.zout' %(object))
                except:
                    status = os.system('echo "" > OUTPUT/%s.zout' %(object))
                    failed = True
                    pass
            
            if not failed:
                status = os.system('grep -v "#" OUTPUT/%s.zout >> %s' %(object, html.replace('html','zout')))
        
        if ('/SED/' in line) & (not failed):
            lines[i] = line.replace('..//SED/','./').replace('SED.png height=180','eazy.png height=250')
        
        if "<table id=" in line:
            insert_line = i
    
    #### photz-specz plot
    unicorn.analysis.show_triple_zphot_zspec(zout=html.replace('html','zout'), zmin=0, zmax=2.5)
    lines.insert(insert_line, '<a href=%s_zz.pdf><img src=%s_zz.png></a>\n' %(html.replace('.html',''), html.replace('.html','')))
    
    fp = open(html.replace('.html','_eazy.html'),'w')
    fp.writelines(lines)
    fp.close()
    
def show_triple_zphot_zspec(zout='massive.zout', zmin=0, zmax=4):
    import threedhst.catIO as catIO
    zo = catIO.Readfile(zout)
    dz = (zo.z_peak-zo.z_spec)/(1+zo.z_spec)
    dz = dz[zo.z_spec > 0]
    
    threedhst.utils.biweight(dz[0::3])
    
    ### setup
    plt.rcParams['font.size'] = 10
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[7,2.6],dpi=100)
    else:
        fig = Figure(figsize=[7,2.6],dpi=100)
    
    #
    #fig = plt.figure(figsize=[8, 2.6],dpi=100)
    
    fig.subplots_adjust(wspace=0.27,hspace=0.12,left=0.07,
                        bottom=0.15,right=0.98,top=0.98)
    
    #################################### Broad-band + spectrum
    ax = fig.add_subplot(131)
    
    point_size = 4
    
    ax.plot(zo.z_spec[1::3], zo.z_peak[1::3], marker='o', linestyle='None', alpha=0.5, color='orange', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[1::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.text(0.5, 0.9, r'$N=$'+'%0d' %(len(dz[1::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{phot}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)

    ax = fig.add_subplot(132)
    
    ax.plot(zo.z_spec[0::3], zo.z_peak[0::3], marker='o', linestyle='None', alpha=0.5, color='purple', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[0::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{phot+grism}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)
    
    ax = fig.add_subplot(133)
    
    ax.plot(zo.z_spec[2::3], zo.z_peak[2::3], marker='o', linestyle='None', alpha=0.5, color='blue', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[2::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{grism\ only}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)
    
    #
    outfile = zout.replace('.zout','_zz')
    if USE_PLOT_GUI:
        plt.savefig(outfile+'.png')
        plt.savefig(outfile+'.pdf')
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile+'.png', dpi=100, transparent=False)
        canvas.print_figure(outfile+'.pdf', dpi=100, transparent=False)
    

def run_eazy_on_all_objects(field='ERS', pipe=' > eazy.log', compress=0.75):
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')

    ######
    #logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/processed.log'

    logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/'+field+'.log'
    #field = 'ERS'

    # logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/goods-s.log'
    # field = 'GOODS-S'
    
    ######
    
    if not os.path.exists(logfile):
        fp = open(logfile,'w')
        fp.write('')
        fp.close()
        
    fp = open(logfile)
    log_lines = fp.readlines()
    fp.close()
    
    os.chdir(unicorn.GRISM_HOME+field)
    catalogs = glob.glob('HTML/*drz.cat')
    print field, len(catalogs)
    
    for catalog in catalogs:
        os.chdir(unicorn.GRISM_HOME+field)
        if not os.path.exists(catalog.replace('HTML','HTML/SED').replace('drz','match')):
            continue
        else:
            print catalog
        cat = threedhst.sex.mySexCat(catalog)
        
        pointing = os.path.basename(catalog).split('G141')[0]+'G141'
        for id in cat.id[np.cast[float](cat.FCOVER) > 0.4]:
            object = '%s_%05d' %(pointing, id)
            if (object+'\n' not in log_lines) & (os.path.exists(unicorn.GRISM_HOME+field+'/HTML/ascii/'+object+'.dat')):
                result = True
                #try:
                result = unicorn.analysis.run_eazy_fit(root=pointing, id=id, compress=compress, zmin=0.02, zmax=4, TILT_ORDER=1, pipe=pipe)
                #except:
                #    pass
                #    
                if result is False:
                    print object
                    root=pointing
                    status = os.system('rm %s_%05d' %(root, id) + '_threedhst.cat')
                    status = os.system('rm %s_%05d' %(root, id) + '.FILT.RES')
                    status = os.system('rm %s_%05d' %(root, id) + '.eazy.param')
                    status = os.system('rm templates/%s_%05d' %(root,id)+'.spectra.param')   
                    status = os.system('rm templates/%s_%05d' %(root,id)+'_spectrum.dat')
                #
                fp = open(logfile,'a')
                fp.write('%s\n' %(object))
                fp.close()
                #print unicorn.analysis.get_open_fds()

#
def get_open_fds():
    """
    Show open files because was getting "too many open files" error.
    """ 
    import fcntl, resource
    
    fds = []
    for fd in range(3,1024):
            try:
                    flags = fcntl.fcntl(fd, fcntl.F_GETFD)
            except IOError:
                    continue
			#
            fds.append(fd)
	#
    return fds

def generate_rf_color_inputs(FILTERS_RES='FILTER.RES.v8.R300', rf_filter_ids=range(135,164), filename='rf_dummy', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param'):
    """
    Make a fake catalog with all of the filters to interpolate. 
    Run with the line templates and Z_MAX = 0 to just get the template fluxes at z=0.
    """
    import threedhst.eazyPy as eazy
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
        eazy_binary = '/usr/local/bin/eazy_latest'
    else:
        eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy'
    
    #rf_filter_ids = range(135,164)
    hline = '# id '
    fluxes = ' 0 '
    for filt in rf_filter_ids:
        hline += ' F%0d E%0d' %(filt, filt)
        fluxes += ' 1.0 0.1'
    
    fp = open(filename+'.cat','w')
    fp.write(hline+'\n')
    fp.write(fluxes+'\n')
    fp.close()
    
    try:
        os.system('rm zphot.param.default')
    except:
        pass
    
    #### Generate default param file    
    os.system(eazy_binary + ' > /tmp/log')
    zp = eazy.EazyParam('zphot.param.default')
    
    zp.params['FILTERS_RES'] = FILTERS_RES
    zp.params['TEMPLATES_FILE'] = TEMPLATES_FILE
    zp.params['TEMPLATE_COMBOS'] = '1'
    zp.params['TEMP_ERR_FILE'] = 'templates/TEMPLATE_ERROR.eazy_v1.0'
    zp.params['WAVELENGTH_FILE'] = 'templates/EAZY_v1.1_lines/lambda_v1.1.def'
    zp.params['CATALOG_FILE'] = filename+'.cat'
    zp.params['Z_MIN'] = 0
    zp.params['Z_MAX'] = 4.0
    zp.params['Z_STEP'] = 0.1
    zp.params['Z_STEP_TYPE'] = 0   
    zp.params['MAIN_OUTPUT_FILE'] = filename
    zp.params['CACHE_FILE'] = filename+'.tempfilt'
    zp.params['APPLY_PRIOR'] = 'y'
    zp.params['PRIOR_FILTER'] = 153
    
    zp.write(filename+'.zphot.param')
    
    os.system(eazy_binary+' -p '+filename+'.zphot.param > /tmp/log')
    
def make_rf_flux_catalog():
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    #### UBV (Maiz-Apellaniz), ugriz (SDSS), JHK (2MASS), UV1600, UV2800
    filters = [153,154,155,156,157,158,159,160,161,162,163,218,219]
    
    unicorn.analysis.generate_rf_color_inputs(FILTERS_RES='FILTER.RES.v8.R300', rf_filter_ids=filters, filename='rf_dummy')
    
    zp_rf = eazy.EazyParam('OUTPUT/rf_dummy.param')
    head_line = '# id z_grism  DM '
    legend = ['# %s\n# \n# Only use grism+photometry redshifts\n# \n' %(time.ctime())]
    
    for filt in zp_rf.filters:
        head_line += ' L%0d' %(filt.fnumber)
        legend.append('# %0d: %s, %12.2f\n' %(filt.fnumber, filt.name, filt.lambda_c))
    
    legend.append('# \n')
    
    fp = open('full_rf_fluxes.cat','w')
    fp.write(head_line+'\n')
    fp.writelines(legend)
    
    files=glob.glob('OUTPUT/*G141*param')
    if len(files) == 0:
        status = os.system('ls OUTPUT |grep G141 |grep param > /tmp/list')
        fpf = open('/tmp/list')
        lines = fpf.readlines()
        fpf.close()
        for line in lines:
            files.append('OUTPUT/'+line[:-1])
    
    for file in files:
        print noNewLine +file
        
        root=os.path.basename(file).split('G141')[0]+'G141'
        id=int(os.path.basename(file).split('G141_')[1].split('.param')[0])
        
        try:
            zfit, DM, obs_sed_rf, filts = unicorn.analysis.get_rf_fluxes(root=root, id=id, dummy='rf_dummy', verbose=False, check_plot=False)
        except:
            zfit, DM, obs_sed_rf = np.zeros(3)-1, np.zeros(3)-1, np.zeros((len(zp_rf.filters),3))-1
            
        #### Only get the results for the first fit = spectrum + photometry
        cat_line = ' %s_%05d  %9.5f %6.2f' %(root, id, zfit[0], DM[0])
        for i in range(len(zp_rf.filters)):
            cat_line += ' %12.5e' %(obs_sed_rf[i, 0])
            
        fp.write(cat_line+'\n')
        
    fp.close()
    
    
def get_rf_fluxes(root='GOODS-S-24-G141', id=27, dummy='rf_dummy', verbose=True, check_plot=False):
    """
    Compute the rest-frame colors using the z=0 interpolated fluxes in the dummy 
    catalog and the coefficients from the actual fit.
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    zp_rf = eazy.EazyParam('OUTPUT/%s.param' %(dummy))
    zp_obj = eazy.EazyParam('OUTPUT/%s_%05d.param' %(root, id))

    tempfilt_rf, coeffs_rf, temp_seds_rf, pz_rf = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s' %(dummy), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    tempfilt_obj, coeffs_obj, temp_seds_obj, pz_obj = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    NTEMP = tempfilt_rf['NTEMP']
    ### Possible that different templates used for each
    if zp_rf.templates != zp_obj.templates:
        if verbose:
            print ' *** WARNING, different template sets used ***'
            print '\n Rest-frame:\n=========='
            for t in zp_rf.templates:
                print t
            #
            print '\n Object:\n=========='
            for t in zp_obj.templates:
                print t
            #
        if tempfilt_rf['NTEMP'] > tempfilt_obj['NTEMP']:
            NTEMP = tempfilt_obj['NTEMP']
    
    ##### Redshift from the fit
    zfit = tempfilt_obj['zgrid'][coeffs_obj['izbest']]
    DM = zfit*0
    for i, zi in enumerate(zfit):
        cc = cosmocalc.cosmocalc(zfit[i], H0=zp_obj['H0'], WM=zp_obj['OMEGA_M'])
        DM[i] = 5.*np.log10(cc['DL_Mpc']*1.e5/np.sqrt(1+zfit[i])) #### Last (1+z) term is used because f_nu is interpolated in observed-frame units.  

    ##### The following checks if a given emission line falls within DELTA_LINE Ang.
    ##### (observed frame) of the central wavelength of one of the observed 
    ##### filters.  If so, we assume that the photometry constrains the strength
    ##### of the line in the fit.  If not, remove the line from the rf_template flux.
    
    DELTA_LINE = 300 ### Angstroms, observed frame
    lc_rf = tempfilt_rf['lc'].copy()
    lc_obj = tempfilt_obj['lc'].copy()
    dlam_spec = lc_obj[-1]-lc_obj[-2]
    is_spec = np.append(np.abs(1-np.abs(lc_obj[1:]-lc_obj[0:-1])/dlam_spec) < 0.05,True)
    
    line_names = ['Ha','OIII','Hb','OII_3727']
    line_waves = [6563., 5007., 4861., 3727.]
    NLINES = len(line_names)
    
    # print DELTA_LINE/dlam_spec
    
    for i, temp in enumerate(zp_obj.templates[0:NTEMP]):
        value = 1.
        for j in range(NLINES):
            if line_names[j] in temp:
                line_obs = line_waves[j]*(1+zfit[0])
                dline = np.abs(line_obs-lc_obj)
                dline[is_spec] *= DELTA_LINE/dlam_spec
                if verbose:
                    print 'Closest filter, %8s: %7.1f A' %(line_names[j], dline.min())
                if dline.min() > DELTA_LINE:
                    if verbose:
                        print 'z=%.2f, setting %s template to zero.' %(zfit[0], line_names[j])
                    coeffs_obj['coeffs'][i,:] *= value
                    tempfilt_rf['tempfilt'][:,i,0] = (lc_rf/5500.)**2
                        
    ##### Rest-frame fnu fluxes, zeropoint from the EAZY run.
    idx=0
    obs_sed_rf = np.dot(tempfilt_rf['tempfilt'][:,0:NTEMP,0],\
                     coeffs_obj['coeffs'][0:NTEMP,:])
    
    #### Check results    
    if check_plot:
        lambdaz, temp_sed, lci, obs_sed, fobs_new, efobs_new = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        plt.plot(lci, obs_sed, marker='o', linestyle='None', alpha=0.5)
        plt.semilogx()

        lc_rf = np.arange(len(zp_rf.filters)*1.)
        for i, filt in enumerate(zp_rf.filters):
            lc_rf[i] = filt.lambda_c

        plt.plot(lc_rf*(1+zfit[0]), obs_sed_rf[:,0]/(lc_rf*(1+zfit[0])/5500.)**2, marker='o', color='red', markersize=10, linestyle='None', alpha=0.3)
        
    return zfit, DM, obs_sed_rf, zp_rf.filters
    
def make_eazy_inputs(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', check=False, bin_spec=1, spec_norm=1., zmin=None, zmax=None, zstep=0.0025, compress=1.0, TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', TILT_COEFFS=[0, 1], eazy_working_directory=None):
    
    import unicorn.analysis
    
    # root='WFC3-ERSII-G01-G141'; id=197; OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; check=False; bin_spec=1; spec_norm=1.; zmin=None; zmax=None; zstep=0.0025; compress=1.0; TEMPLATES_FILE='templates/o2_fit_lines.spectra.param'; TILT_COEFFS=[0, 1]
    from scipy import polyval
    
    unicorn.analysis.BAD_SPECTRUM = False
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    os.chdir(eazy_working_directory)
    
    if (('MARSHALL' in root) | ('UDS' in root)) & ('UDS' not in OLD_RES):
        OLD_RES='FILTER_UDS.RES'
        
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
    
    #### Dummy to get some of the output variables    
    ok = (match.rmatch < 1) & (match.logm > 10.)
    #print root, len(match.rmatch), np.max(match.logm), match.id_f140w[ok][0]
    result = False
    i=0
    while (result is False) & (i < len(match.rmatch[ok])):
        result = unicorn.analysis.specphot(id=match.id_f140w[ok][i],
            grism_root=root, SPC = SPC, 
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same',
            GET_SPEC_ONLY = True)
        i+=1
    #
    lam, spflux, sperr, lci, fobs, efobs, photom_idx = result
    
    #### Match 3D-HST object to the ancillary catalogs to get photometry
    result = unicorn.analysis.specphot(id=id,
        grism_root=root, SPC = SPC, 
        cat = cat,
        grismCat = grismCat, zout = zout, fout = fout, 
        OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
        MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
        OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
        CACHE_FILE = 'Same',
        GET_SPEC_ONLY = True)
    
    if result is not False:
        #### Object is matched in the photometry
        lam, spflux, sperr, lci, fobs, efobs, photom_idx = result
        z_spec = zout.z_spec[photom_idx]
        unicorn.analysis.HAS_PHOTOMETRY = True
        unicorn.analysis.PHOTOMETRY_ID = zout.id[photom_idx]
        
    else:
        #### No match
        sp = threedhst.catIO.Readfile('HTML/ascii/%s_%05d.dat' %(root, id))
        lam, spflux, sperr = sp.lam, (sp.flux-sp.contam)/1.e-17*(sp.lam/5500.)**2, np.sqrt(sp.error**2+(1.0*sp.contam)**2)/1.e-17*(sp.lam/5500.)**2
        photom_idx = 0
        z_spec = -1.0
        fobs = fobs*0-999
        efobs = fobs
        print 'No photometry!'
        unicorn.analysis.HAS_PHOTOMETRY = False
        unicorn.analysis.PHOTOMETRY_ID = None
        
    use = (lam > 1.1e4) & (lam < 1.65e4) & (spflux != 0.0) & np.isfinite(spflux) & np.isfinite(sperr)
    
    #### allow additional normalization term
    spflux *= spec_norm
    sperr *= spec_norm
    
    #### Apply the linear fit
    spflux *= polyval(TILT_COEFFS, lam-1.4e4)
    
    #### check plot
    if check:
        plt.semilogx([1],[1])
        plt.errorbar(lci, fobs, yerr=efobs, marker='o', linestyle='None', color='blue', ecolor='blue')
        ma = fobs.max()
        plt.plot(lam[use], spflux[use], color='red')
        plt.ylim(-0.1*ma, 1.5*ma)
        plt.xlim(3000, 9.e4)
    
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
    
    if len(lam[use]) < 10:
        unicorn.analysis.BAD_SPECTRUM = True
        return False
        
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
    
    #### Make lines slightly narrower than the direct image
    lsm -= lam[use][0]
    lsm *= compress
    lsm += lam[use][0]
    
    ###################################
    #### Make a RES file for the convolved resolution and make a catalog file with
    #### photometry and spectroscopy
    
    os.chdir(ORIG_PATH)
    
    ##### Have to fill the empty parts of the spectrum with something other than 0
    fill_value = np.median(spflux[(lam > 1.2e4) & (lam < 1.6e4) & (spflux > 0)])
    
    ##### print the spectrum to a file in the templates directory
    fp = open('templates/%s_%05d' %(root, id)+'_spectrum.dat','w')
    fp.write('50 %.5e\n%.5e %.5e\n' %(fill_value, lam[0]-dlam, fill_value))
    #spflux[~np.isfinite(spflux) | (spflux < 0)] = 1.e-8
    for i in range(len(lam)):
        if np.isfinite(spflux[i]) & (spflux[i] > 0):
            fp.write('%.5e %.5e\n' %(lam[i], spflux[i]))
    
    fp.write('%.5e %.5e\n2.e7 %.5e\n' %(lam[-1]+dlam, fill_value, fill_value))
    fp.close()
    
    #### Make the template file for this single grism spectrum
    fp = open('templates/%s_%05d' %(root, id)+'.spectra.param','w')
    fp.write('1 templates/%s_%05d_spectrum.dat 1.0 0 1.0\n' %(root, id))
    fp.close()
    
    fp = open(OLD_RES)
    res_lines = fp.readlines()
    fp.close()
    
    nfilt=0
    res_name = []
    res_lc = []
    for line in res_lines:
        if line.find('lambda_c') > 0:
            res_name.append(line.split()[1])
            nfilt += 1
    
    #### Read the parameter file associated with the zout file to get the right
    #### filters, which should be in the same order as the vectors read from the 
    #### binary files
    eazy_param = threedhst.eazyPy.EazyParam(zout.filename.replace('zout','param'))
    
    cat_head = '# phot_id id ra dec z_spec'
    cat_line = ' %d %s_%05d %f %f %8.3f' %(cat.id[photom_idx], root, id, cat.ra[photom_idx], cat.dec[photom_idx], z_spec)
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
    
    #print '\n\nMAX: %f\n\n' %np.max(fnu_spec[use])
    if np.max(fnu_spec[use]) <= 0:
        fnu_spec[use] = fnu_spec[use]*0.-100
        efnu_spec[use] = efnu_spec[use]*0.-100
        unicorn.analysis.BAD_SPECTRUM = True
        
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
    
    fp = open('%s_%05d' %(root, id) + '.FILT.RES','w')
    fp.writelines(res_lines)
    fp.close()
    
    fp = open('%s_%05d' %(root, id) + '_threedhst.cat','w')
    fp.write(cat_head+'\n')
    fp.write(cat_line+'\n')
    fp.write(no_spec+'\n')
    fp.write(no_phot+'\n')
    fp.close()
    
    #### Set some EAZY parameters
    
    eazy_param.params['READ_ZBIN'] = 'n'
    eazy_param.params['TEMPLATES_FILE'] = TEMPLATES_FILE
    eazy_param.params['FILTERS_RES'] = '%s_%05d' %(root, id) + '.FILT.RES'
    eazy_param.params['OUTPUT_DIRECTORY'] = 'OUTPUT'
    eazy_param.params['WAVELENGTH_FILE'] = 'templates/EAZY_v1.1_lines/lambda_v1.1.def'
    eazy_param.params['CATALOG_FILE'] = '%s_%05d' %(root, id) + '_threedhst.cat'
    eazy_param.params['MAIN_OUTPUT_FILE'] = '%s_%05d' %(root, id)
    eazy_param.params['CACHE_FILE'] = '%s_%05d.tempfilt' %(root, id)
    eazy_param.params['GET_ZP_OFFSETS'] = 0
    eazy_param.params['Z_STEP'] = zstep
    eazy_param.params['NMF_TOLERANCE'] = 1.e-4
    
    if zmin is None:
        zmin = zout.l99[photom_idx]*0.8
    if zmax is None:
        zmax = zout.u99[photom_idx]/0.8
        
    if zmin < 0:
        zmin = 0
    
    if zmax < zmin:
        zmax = zmin+1.e-3
        
    eazy_param.params['Z_MIN'] = zmin
    eazy_param.params['Z_MAX'] = zmax 
    
    eazy_param.params['MAGNITUDES'] = 0.0
    eazy_param.params['REST_FILTERS'] = '-,-'
    
    eazy_param.write(file='%s_%05d' %(root, id) + '.eazy.param')
    
def trim_jh_filters(input='FILTER.RES.v8.R300', output='FILTER.RES.v8.R300.trim', wlo=1.08e4, whi=1.65e4):
    """
    Need to trim H & J filters such that they don't sample the parts of the 
    grism spectra where the sensitivity blooms up.
    
    Loop through the filters in an EAZY res file.  If the filter has non-zero 
    transmission at the blue or red edges, trim it off.
    """
    import threedhst.eazyPy as eazy
    
    res = eazy.FilterFile(input)
    for i, filter in enumerate(res.filters):
        peak = filter.transmission.max()
        limits = np.interp([wlo, whi], filter.wavelength, filter.transmission)
        if np.max(limits/peak) > 1.e-1:
            print 'Trim filter: %s' %(filter.name.split()[0])
            keep = (filter.wavelength >= wlo) & (filter.wavelength <= whi)
            res.filters[i].wavelength = filter.wavelength[keep]
            res.filters[i].transmission = filter.transmission[keep]
    
    res.write(file=output)
    
def scale_to_photometry(root='GOODS-S-24-G141', id=23, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy', compress=1.0, check_results=False, spec_norm=1.0, ORDER=1, pipe=' > log', eazy_working_directory=None):
    import unicorn
    import threedhst.eazyPy as eazy
    from scipy import polyfit, polyval
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    os.chdir(eazy_working_directory)
    
    if ('MARSHALL' in root) | ('UDS' in root):
        OLD_RES='FILTER_UDS.RES'
    
    ##### Need to adjust the J/H filters in the FILTER files near the 
    ##### edges of the grism sensitivity
    if not os.path.exists(OLD_RES+'.trim'):
        unicorn.analysis.trim_jh_filters(input=OLD_RES, output=OLD_RES+'.trim')
        
    #####  Get the f_lambda fluxes with the original (not trimmed) filters
    unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=1.0, spec_norm=spec_norm, zmin=0.0000, zmax=1.e-6, compress=compress, TEMPLATES_FILE='templates/%s_%05d' %(root, id)+'.spectra.param', eazy_working_directory=eazy_working_directory)
    
    if unicorn.analysis.BAD_SPECTRUM:
        return [0,1]
    
    if unicorn.analysis.HAS_PHOTOMETRY is False:
        return [0,1]
    
    status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
    
    lambdaz, temp_sed_0, lci, obs_sed_0, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    #####  Get the f_lambda fluxes with the trimmed filters
    unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES+'.trim', bin_spec=1.0, spec_norm=spec_norm, zmin=0.0000, zmax=1.e-6, compress=compress, TEMPLATES_FILE='templates/%s_%05d' %(root, id)+'.spectra.param', eazy_working_directory=eazy_working_directory)
    
    status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
    
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    lambdaz, temp_sed, lci, obs_sed, fobs_new, efobs_new = eazy.getEazySED(2, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    if np.sum(obs_sed) == 0:
        obs_sed = obs_sed_0
        temp_sed = temp_sed_0
    
    if np.sum(obs_sed) == 0:
        print "Problem with the eazy fit to get normalization..."
        return [0, 1]
        
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        
    #### If spectrum doesn't cover both J/H, stop and return no correction
    if (lci[is_spec].min() > 1.15e4) | (lci[is_spec].max() < 1.58e4):
        print "Spectrum doesn't have full coverage: [%f, %f]" %(lci[is_spec].min(), lci[is_spec].max())
        return [0, 1]
    
    #### Find filters near J/H
    jhfilt = ~is_spec & (lci > 1.15e4) & (lci < 1.65e4) & (fobs > 0)
    
    #### No photometry:
    if fobs[~is_spec].max() <= 0:
        return [0, 1]
        
    #### Here's the simple linear fit    
    xfit = lci-1.4e4
    yfit = fobs/obs_sed
    
    if len(xfit[jhfilt]) == 0:
        print "No valid J/H filters found."
        return [0, 1]
        
    if len(xfit[jhfilt]) == 1:
        ORDER=0
        
    afit = polyfit(xfit[jhfilt], yfit[jhfilt], ORDER)
    afit[ORDER] *= 1./(coeffs['tnorm'][0]/coeffs['coeffs'][0,2])
    
    #### Reduce the effect slightly
    if ORDER == 1:
        afit[0] *= 0.8
    
    #### Diagnostic plot
    if check_results:
        if unicorn.analysis.USE_PLOT_GUI:
            fig = plt.figure(figsize=[5,5],dpi=100)
        else:
            fig = Figure(figsize=[5,5],dpi=100)

        fig.subplots_adjust(wspace=0.04,hspace=0.12,left=0.05,
                            bottom=0.15,right=0.98,top=0.98)
        #
        ax = fig.add_subplot(211)
        ax.plot(lambdaz, temp_sed)
        ax.plot(lci, fobs, marker='o', markersize=5, alpha=0.5, linestyle='None', color='yellow')
        ax.plot(lci[jhfilt], fobs[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='red')
        ax.plot(lci[jhfilt], obs_sed[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='blue')
        ax.plot(lci, obs_sed, marker='o', markersize=5, alpha=0.5, linestyle='None', color='green')
        ax.set_ylim(0, 1.1*obs_sed.max())
        ax.set_xlim(9.e3, 1.8e4)
                
    if check_results:
        print xfit[jhfilt], yfit[jhfilt]
        print afit
        
        ax.plot(lambdaz, polyval(afit, lambdaz-1.4e4), color='orange', alpha=0.4)
        
        ax = fig.add_subplot(212)
        ax.plot(lambdaz, polyval(afit, lambdaz-1.4e4)*temp_sed)
        ax.plot(lci, polyval(afit, lci-1.4e4)*fobs, marker='o', markersize=5, alpha=0.5, linestyle='None', color='yellow')
        ax.plot(lci[jhfilt], fobs[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='red')
        ax.plot(lci[jhfilt], polyval(afit, lci[jhfilt]-1.4e4)*obs_sed[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='blue')
        ax.plot(lci, polyval(afit, lci-1.4e4)*obs_sed, marker='o', markersize=5, alpha=0.5, linestyle='None', color='green')
        ax.set_ylim(0, 1.1*obs_sed.max())
        ax.set_xlim(9.e3, 1.8e4)
        outfile='scale_to_photometry.png'
        if USE_PLOT_GUI:
            plt.savefig(outfile)
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(outfile, dpi=100, transparent=False)
    
    #print afit
    
    return afit
    
def run_eazy_fit(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v9.R300', OUT_RES = 'THREEDHST.RES', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = None, zmin=None, zmax=None, zstep=0.001, compress=1.0, GET_NORM=False, COMPUTE_TILT=True, TILT_ORDER=0, clean=True, force_zrange=False, eazy_working_directory=None):
    
    # OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; TEMPLATES_FILE='templates/o2_fit_lines.spectra.param'; run=True; pipe=' > log'; bin_spec=1; spec_norm=1; eazy_binary = None; zmin=None; zmax=None; compress=1.0; GET_NORM=False; COMPUTE_TILT=True; TILT_ORDER=0; clean=True
    import matplotlib.pyplot as plt
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import unicorn.analysis
    
    t0 = time.time()
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    if (eazy_binary is None):
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            eazy_binary = '/usr/local/bin/eazy_latest'
        else:
            eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy'
    
    MAXIT = 3.
    tilt = [0,1]
    
    if ('MARSHALL' in root) | ('UDS' in root):
        OLD_RES='FILTER_UDS.RES'
    
    if run:
        ########################### Scale to photometry
        if COMPUTE_TILT:
            tilt = unicorn.analysis.scale_to_photometry(root=root, id=id, OLD_RES = OLD_RES, OUT_RES = OUT_RES, eazy_binary = eazy_binary, compress=compress, ORDER=TILT_ORDER, pipe=pipe)
        else:
            tilt = [0, 1]
        
        tnorm = time.time()
        
        if unicorn.analysis.BAD_SPECTRUM:
            return False
        
        ########################## Now run
        #### If fitting with photometry, first run with eazy line templates 
        #### and coarse sampling
        if unicorn.analysis.HAS_PHOTOMETRY:
            unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.01, compress=1.5, TILT_COEFFS=tilt, TEMPLATES_FILE='templates/eazy_v1.1_lines_suppl.spectra.param', eazy_working_directory=eazy_working_directory)

            status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
            
            unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.001, compress=compress, TILT_COEFFS=tilt, TEMPLATES_FILE=TEMPLATES_FILE, eazy_working_directory=eazy_working_directory)
            
            ztmp = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
            
            eazy_param = eazy.EazyParam('%s_%05d.eazy.param' %(root, id))
            eazy_param.params['TEMPLATES_FILE'] = TEMPLATES_FILE
            eazy_param.params['Z_STEP'] = zstep
            eazy_param.params['Z_STEP_TYPE'] = 0
                        
            if force_zrange:
                eazy_param.params['Z_MIN'] = np.max([zmin,0])
                eazy_param.params['Z_MAX'] = zmax                
            else:
                eazy_param.params['Z_MIN'] = np.max([ztmp.l99[1]-1*0.05*(1+ztmp.z_peak[1]),0])
                #
                zma = ztmp.u99[1]+1*0.05*(1+ztmp.z_peak[1])
                #### some very blue SEDs have zmin < 0, zmax~0.  Fit the full range
                if zma < 0.2:
                    zma = 3
                #   
                eazy_param.params['Z_MAX'] = zma
                print 'Refit, fine sampling: [%.2f, %.2f]' %(eazy_param.params['Z_MIN'], eazy_param.params['Z_MAX'])

            eazy_param.write(file='%s_%05d' %(root, id) + '.eazy.param')
            
        else:
            ##### No photometry found.  Look for emission lines in the spectrum, 
            ##### and if one line found, fit it assuming OIII.  If more than one line
            ##### found, fit full redshift range.
            spec = catIO.Readfile(unicorn.analysis.get_grism_path(root)+'/HTML/ascii/%s_%05d.dat' %(root, id))
            
            found_lines = threedhst.spec1d.findLines(spec, trim_abs=True)
            if found_lines is None:
                return False
                
            if len(found_lines) == 1:
                zmin = found_lines[0].wave/5007.-1
                zmax = found_lines[0].wave/4863.-1
                zmin -= 0.1*(1+zmax)
                zmax += 0.1*(1+zmax)
                
            if len(found_lines) > 1:
                wmin = 1.e5
                wmax = 0
                for line in found_lines:
                    if line.wave < wmin: 
                        wmin = line.wave
                    if line.wave > wmax:
                        wmax = line.wave
                
                zmin = wmin/5007.-1-0.1
                zmax = wmax/4863.-1+0.1
                zmin -= 0.1*(1+zmax)
                zmax += 0.1*(1+zmax)
            
            tilt = [0, 1]
            unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.002, compress=compress, TILT_COEFFS=tilt, TEMPLATES_FILE=TEMPLATES_FILE, eazy_working_directory=eazy_working_directory)
        
        status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
        
        ztmp = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
        zstep_i = zstep
        SHOW_ZOUT_FILE = 'OUTPUT/%s_%05d.zout' %(root, id)
        
        ####### 99% confidence interval is not resolved with z_step.  Shrink the step
        while (ztmp.u99[0]-ztmp.l99[0])/zstep_i <= 9.99:
            resolve_factor = (ztmp.u99[0]-ztmp.l99[0])/zstep_i
            eazy_param.params['Z_MIN'] = ztmp.l99[0]-zstep_i*10
            eazy_param.params['Z_MAX'] = ztmp.u99[0]+zstep_i*10
            eazy_param.params['MAIN_OUTPUT_FILE'] = '%s_%05d_refine' %(root, id)
                           
            zstep_i = (ztmp.u99[0]-ztmp.l99[0])/10.
            eazy_param.params['Z_STEP'] = zstep_i
            print 'N=%f, Shrink Z_STEP: %f, [%f, %f]\n' %(resolve_factor, zstep_i, eazy_param.params['Z_MIN'], eazy_param.params['Z_MAX'])
            
            eazy_param.write(file='%s_%05d_refine' %(root, id) + '.eazy.param')
            
            status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
            ztmp = catIO.Readfile('OUTPUT/%s_%05d_refine.zout' %(root, id))
            
            SHOW_ZOUT_FILE = 'OUTPUT/%s_%05d_refine.zout' %(root, id)
            
        if tilt is False:
            tilt = [0,1]
        
        if len(tilt) == 1:
            tilt = [0,tilt[0]]
        
        print 'Tilt: %e %e\n' %(tilt[0], tilt[1])
        
        fp = open('%s_%05d.tilt' %(root, id),'w')
        fp.write('%s_%05d  %.4e %.4e\n' %(root, id, tilt[0], tilt[1]))
        fp.close()
        
        if unicorn.analysis.BAD_SPECTRUM:
            return False
        
        tfit = time.time()
        
        # lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        #     eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), \
        #                       OUTPUT_DIRECTORY='OUTPUT', \
        #                       CACHE_FILE = 'Same')
        # #
        # dlam_spec = lci[-1]-lci[-2]
        # is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
            
        #### Show the results
        try:
            status = os.system('head -3 %s |tail -1' %(SHOW_ZOUT_FILE))
        except:
            pass
    else:
        tnorm = time.time()
        tfit = time.time()
                    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), \
                          OUTPUT_DIRECTORY='OUTPUT', \
                          CACHE_FILE = 'Same')
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    #### check normalization
    spec_norm = np.sum(obs_sed[is_spec]*fobs[is_spec]/efobs[is_spec]**2)/np.sum(obs_sed[is_spec]**2/efobs[is_spec]**2)

    print spec_norm    
    if GET_NORM:
        return spec_norm
        
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
    if len(fobs[is_spec & (fobs > 0)]) < 5:
        return False
        
    ymax = max(fobs[is_spec & (fobs > 0)])
    
    ## photometry
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.9, color='orange', markersize=6)
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='red', linestyle='None', markersize=6, alpha=0.2)
    
    ## best-fit SED
    ax.plot(lambdaz, temp_sed, color='red', alpha=0.6)
    
    ## Spectrum + convolved fit
    ax.plot(lci[is_spec], obs_sed[is_spec], color='red', markersize=6, alpha=0.7, linewidth=1)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.8, color='blue', linewidth=1)
        
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
    
    #### Photometry ID label
    if unicorn.analysis.PHOTOMETRY_ID is not None:
        ax.text(0.1, 0.95, '%0d' %(unicorn.analysis.PHOTOMETRY_ID), transform = ax.transAxes, horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax.set_xlim(0.95e4,1.78e4)
    # ax.set_ylim(0.8*ymin, ymax*1.1)
    ax.set_ylim(-0.1*ymax, 1.2*ymax)
    
    #################################### p(z) for combined, photometry, and spec 
    ax = fig.add_subplot(133)
    
    colors = ['purple','orange','blue']
    alpha = [0.5, 0.5, 0.2]
    zo = threedhst.catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
    zmin = 4
    zmax = 0
    ymax = 0
    for i in range(3):
        zgrid, pz = eazy.getEazyPz(i, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), 
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
    ax.xaxis.set_major_locator(MyLocator(4, prune='both'))
    
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
    
    tplot = time.time()
    
    print '%s, %.2f + %.2f + %.2f = %.2f s' %(outfile, tnorm-t0, tfit-tnorm, tplot-tfit, tplot-t0)
    
    #### Clean up temporary files
    if clean:
        status = os.system('rm %s_%05d' %(root, id) + '_threedhst.cat')
        status = os.system('rm %s_%05d' %(root, id) + '.FILT.RES')
        status = os.system('rm %s_%05d' %(root, id) + '.eazy.param')
        status = os.system('rm templates/%s_%05d' %(root, id) + '.spectra.param')
        status = os.system('rm templates/%s_%05d' %(root, id) + '_spectrum.dat')
    
    return True
#

def check_eazy_fits():
    """
    Check how many templates were run in the eazy fit in the v1.6 directory to
    assess how many objects were overwritten from the original eazy run that 
    generated the catalog redshifts.
    """
    import threedhst.eazyPy as eazy
    
    os.chdir('/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6')
    
    os.system('ls OUTPUT |grep param > /tmp/param_files')
    fp = open('/tmp/param_files')
    param_files = fp.readlines()
    fp.close()
    
    fp = open('param_log','w')
    fp.write('# object NTEMP zmin zmax temp_file\n')
    
    for file in param_files:
        print noNewLine+file[:-1]
        param = eazy.EazyParam('OUTPUT/%s' %(file[:-1]))
        fp.write('%-30s %2d %6.2f %6.2f %s\n' %(file[:-1], len(param.templates), param['Z_MIN'], param['Z_MAX'], param['TEMPLATES_FILE']))
        
    fp.close()
    
    ##### Redo fits that were overwritten
    log = catIO.Readfile('param_log')
    redo = (log.temp_file != 'templates/o2_fit_lines.spectra.param')
    
    for object in log.object[redo]:
        if object.startswith('UDF'):
            continue
        #
        object = object.split('.param')[0]
        pointing = object.split('-G141')[0]+'-G141'
        id = int(object.split('_')[1])
        result = unicorn.analysis.run_eazy_fit(root=pointing, id=id, compress=0.7, zmin=0.02, zmax=4, TILT_ORDER=1, pipe=' > log', OLD_RES = 'FILTER.RES.v9.R300', TEMPLATES_FILE='templates/o2_fit_lines.spectra.param')
        os.system('cat OUTPUT/%s.zout' %(object))
        print '\n ------ \n'
        os.system('grep %s ../FIRST_PAPER/GRISM_v1.6/full_redshift.cat' %(object))
        os.system('mv OUTPUT/%s* ../REDSHIFT_FITS_v1.6/OUTPUT/' %(object))
    
def make_all_asciifiles():
    """ 
    Make the ascii files for release v1.6
    """
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    os.chdir('/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6')
    
    fields = np.unique(phot.field)
    for field in fields:
        try:
            os.mkdir('ASCII/%s' %(field))
        except:
            pass
            
    fp = open('ASCII/failed.log','w')
    for i in range(len(zout.z_peak[0::3])):
        object = zout.id[0::3][i]
        field = phot.field[phot.idx][i]
        print noNewLine+object
        try:
            unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='./ASCII/%s' %(field))
        except:
            fp.write(object+'\n')
    #
    fp.close()
        
def make_eazy_asciifiles(object='COSMOS-8-G141_00498', eazy_output='./OUTPUT/', savepath='./ASCII/'):
    """
    Make ascii files with the best-fit eazy templates and p(z) for the 3D-HST fit.
    
    Currently only the spec + photometry combined fit is saved in the obs_sed and 
    temp_sed files.  All three fits are saved in the pz.dat file.
    
    """
    
    if not savepath.endswith('/'):
        savepath += '/'
        
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    
    eazy_param = eazy.EazyParam('%s/%s.param' %(eazy_output, object))
    for continuum_i, temp in enumerate(eazy_param.templates):
        if 'single_lines' in temp:
            break
            
    #### Photometry + spectrum
    fp = open(savepath+object+'_obs_sed.dat','w')
    fp.write('# lc fnu efnu obs_sed continuum line is_spec\n')
    
    lci = tempfilt['lc']
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)*1
    
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][0]],\
                     coeffs['coeffs'][:,0])# /(lci/5500.)**2
    
    #
    obs_sed_cont = np.dot(tempfilt['tempfilt'][:,0:continuum_i,coeffs['izbest'][0]], coeffs['coeffs'][0:continuum_i,0])# /(lci/5500.)**2
    obs_sed_line = np.dot(tempfilt['tempfilt'][:,continuum_i:,coeffs['izbest'][0]], coeffs['coeffs'][continuum_i:,0])# /(lci/5500.)**2
    
    for i in range(tempfilt['NFILT']):
        fp.write(' %10.6e %8.3e %8.3e  %8.3e  %8.3e  %8.3e %d\n' %(tempfilt['lc'][i], tempfilt['fnu'][i,0], tempfilt['efnu'][i,0], obs_sed[i], obs_sed_cont[i], obs_sed_line[i], is_spec[i]))
    
    fp.close()
        
    #### template fit
    zi = tempfilt['zgrid'][coeffs['izbest'][0]]
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,0])
    temp_sed *= (lambdaz/5500.)**2/(1+zi)**2

    temp_sed_cont = np.dot(temp_seds['temp_seds'][:,0:continuum_i],coeffs['coeffs'][0:continuum_i,0])
    temp_sed_cont *= (lambdaz/5500.)**2/(1+zi)**2

    temp_sed_line = np.dot(temp_seds['temp_seds'][:,continuum_i:],coeffs['coeffs'][continuum_i:,0])
    temp_sed_line *= (lambdaz/5500.)**2/(1+zi)**2
    
    fp = open(savepath+object+'_temp_sed.dat','w')
    fp.write('# lam fnu_temp continuum line\n')
    for i in range(temp_seds['NTEMPL']):
        fp.write(' %10.6e %8.3e %8.2e %8.2e\n' %(lambdaz[i], temp_sed[i], temp_sed_cont[i], temp_sed_line[i]))
    fp.close()
    
    #### p(z)
    
    zgrid, pz0 = eazy.getEazyPz(0, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    zgrid, pz1 = eazy.getEazyPz(1, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    zgrid, pz2 = eazy.getEazyPz(2, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    
    fp = open(savepath+object+'_pz.dat','w')
    fp.write('#  z pz_both pz_phot pz_spec\n')
    for i in range(len(zgrid)):
        fp.write('%f %.3e %.3e %.3e\n' %(zgrid[i], pz0[i], pz1[i], pz2[i]))
    
    fp.close()
            
def make_eazy_2d_continuum(root='UDF-G141', id=1279):
    """ 
    Make a 2D continuum image from the best-fit EAZY template.
    """
    from scipy import polyval
    
    FIT_DIR = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/'
    
    object = '%s_%05d' %(root, id)
    if not os.path.exists(FIT_DIR+'%s.tilt' %(object)):
        print 'Need to run_eazy_fit(root=\'%s\', id=%d)' %(root, id)
        return -1
    
    tilt = np.cast[float](np.loadtxt(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/%s.tilt' %(object), dtype=np.str)[1:])
    
    PATH = unicorn.analysis.get_grism_path(root)
    twod = pyfits.open(PATH+'HTML/images/%s_2D.fits.gz' %(object))
    thumb = pyfits.open(PATH+'HTML/images/%s_thumb.fits.gz' %(object))
    
    #### Eazy outputs
    eazy_param = eazy.EazyParam(FIT_DIR+'OUTPUT/%s_%05d.param' %(root, id))
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY=FIT_DIR+'OUTPUT', CACHE_FILE = 'Same')
    
    iz = coeffs['izbest'][0]
    z_peak = tempfilt['zgrid'][iz]
    continuum_full = np.dot(temp_seds['temp_seds'][:,0:7], coeffs['coeffs'][0:7,0])                             
    continuum_filt = np.dot(tempfilt['tempfilt'][:,0:7,iz], coeffs['coeffs'][0:7,0])                             
    
    lci = tempfilt['lc']
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    #### template interpolated at image wavelengths
    image_lambda = (np.arange(twod[1].header['NAXIS1'])+1-twod[1].header['CRPIX1'])*twod[1].header['CDELT1']+twod[1].header['CRVAL1']
    
    sens = pyfits.open(unicorn.GRISM_HOME+'CONF/WFC3.IR.G141.1st.sens.2.fits')
    sens_int = np.interp(image_lambda, sens[1].data.WAVELENGTH, sens[1].data.SENSITIVITY)
    

    yflux = np.interp(image_lambda, lci[is_spec], continuum_filt[is_spec], left=0, right=0) *3.e18*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))/image_lambda**2*(1+z_peak)
    
    yelec = yflux*sens_int
    
    #### Add the tilt computed from the fit
    yelec /= polyval(tilt, image_lambda-1.4e4)
    
    #### Normalize by the model profile
    profile = np.sum(twod[5].data, axis=1)
    profile /= np.sum(profile)
    
    #### Normalize by the object profile
    profile = np.sum(thumb[0].data, axis=1)
    profile /= np.sum(profile)
    
    #### Normalize to the spectrum itself
    mask = yelec > 0
    for i in range(len(profile)):
        profile[i] = np.sum(yelec[mask]*(twod[1].data[i,mask]-twod[4].data[i,mask]))/np.sum(yelec[mask]**2)
        
    continuum_model = np.dot(profile.reshape(-1,1), yelec.reshape(1,-1))
    
    pyfits.writeto('%s/CONTINUUM/%s_continuum.fits' %(FIT_DIR,object), continuum_model, header=twod[1].header, clobber=True)
    status = os.system('gzip -f %s/CONTINUUM/%s_continuum.fits' %(FIT_DIR,object))
    
    print '%sCONTINUUM/%s_continuum.fits' %(FIT_DIR,object)
    
    xxx = """
    id=1874
    
    im = pyfits.open('/Users/gbrammer/Downloads/UDF-G141_%05d_2D.fits.gz' %(id))
    model = pyfits.open('/research/HST/GRISM/3DHST/ANALYSIS/REDSHIFT_FITS/CONTINUUM/UDF-G141_%05d_continuum.fits.gz' %(id))
    ds9.frame(1)
    ds9.v(im[1].data-im[4].data, vmin=-0.01, vmax=0.1)
    ds9.frame(2)
    ds9.v(im[1].data-model[0].data-im[4].data, vmin=-0.01, vmax=0.1)
    
    """
def run_FAST_fit(root='COSMOS-8-G141', id=498, OLD_RES = 'FILTER.RES.v9.R300', OUT_RES = 'THREEDHST.RES', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = None, zmin=0.2, zmax=5, compress=0.8, GET_NORM=False, COMPUTE_TILT=True, TILT_ORDER=1, force_zrange=False):
    """
    Run FAST fit on the photometry + spectra.  If the catalog and template file
    for a given object isn't found, run `run_eazy_fit` first.
    
    root='COSMOS-8-G141'; id=498; OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param'; run=True; pipe=' > log'; bin_spec=1; spec_norm=1; eazy_binary = None; zmin=0.2; zmax=5; compress=0.8; GET_NORM=False; COMPUTE_TILT=True; TILT_ORDER=1; force_zrange=False
    """
    import unicorn
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FAST')
    
    object = '%s_%05d' %(root, id)
    
    zout = None
    
    if (not os.path.exists(object+'_threedhst.cat')) | (not os.path.exists(object+'.FILT.RES')):
        unicorn.analysis.run_eazy_fit(root=root, id=id, OLD_RES = OLD_RES, OUT_RES = OUT_RES, TEMPLATES_FILE=TEMPLATES_FILE, run=True, pipe=pipe, bin_spec=1, spec_norm=1, eazy_binary = eazy_binary, zmin=0.2, zmax=5, compress=0.75, GET_NORM=GET_NORM, COMPUTE_TILT=True, TILT_ORDER=1, clean=False, force_zrange=force_zrange, eazy_working_directory=unicorn.GRISM_HOME+'ANALYSIS/FAST')
    
        fp = open(object+'_threedhst.cat')
        lines = fp.readlines()
        fp.close()
        
        zout = catIO.Readfile('OUTPUT/%s.zout' %(object))
        
        #### Put z_peak for phot + grism as z_spec in the catalog    
        lines[0] = lines[0].replace(' z_spec ',' z_spec_old ').replace(' id ',' id_gris ')[:-1]+' id z_spec\n'
        for i in range(1,4):
            lines[i] = lines[i][:-1]+' %d %.4f\n' %(i,zout.z_peak[0])

        #
        fp = open(object+'_threedhst.cat','w')
        fp.writelines(lines)
        fp.close()
    
    if zout is None:
        zout = catIO.Readfile('OUTPUT/%s.zout' %(object))
    
    param = threedhst.eazyPy.EazyParam('OUTPUT/%s.param' %(object))
    
    #### Make a fast.param file
    fp = open('fast.param')  ### default parameters
    lines = fp.readlines()
    fp.close()
    
    ## change for each object
    lines[80] = 'CATALOG        = \'%s_threedhst\'\n' %(object)
    lines[81] = 'AB_ZEROPOINT   = %f\n' %(param.params['PRIOR_ABZP'])
    lines[82] = 'FILTERS_RES    = \'%s.FILT.RES\'\n' %(object)
    lines[231] = 'Z_MIN          = %f\n' %(zout.z_peak[0])
    lines[232] = 'Z_MAX          = %f\n' %(zout.z_peak[0])
    
    fp = open('%s_fast.param' %(object),'w')
    fp.writelines(lines)
    fp.close()
    
    os.system('/usr/local/share/FAST/FAST_v0.9b/fast %s_fast.param' %(object))
    
    #### Fix the FAST output file, put the header first and add the object name
    print 'Fix FAST_OUTPUT/%s_threedhst.fout' %(object)
    
    fp = open('FAST_OUTPUT/%s_threedhst.fout' %(object))
    lines = fp.readlines()
    fp.close()
    
    for i,line in enumerate(lines):
        if line.startswith('#    id'):
            start = i+1
            header = line
            
    fp = open('FAST_OUTPUT/%s_threedhst.fout' %(object),'w')
    fp.write(header.replace(' id ', ' object id'))
    fp.writelines(lines[0:start-1])
    for line in lines[start:]:
        fp.write(object + ' ' + line)
        
    fp.close()
    
    #### Show FAST results
    print lines[start]
    
    #### Make the eazy ASCII files
    print 'Make EAZY ascii files in ./ASCII ...'
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='./ASCII/')
    
class MyLocator(mticker.MaxNLocator):
    """
    Set maximum number of ticks, from
    http://matplotlib.sourceforge.net/examples/pylab_examples/finance_work2.html
    """
    def __init__(self, *args, **kwargs):
        mticker.MaxNLocator.__init__(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return mticker.MaxNLocator.__call__(self, *args, **kwargs)

#########################################
#                                       #
#                                       #
#     Fit equivalent widths             #
#                                       #
#                                       #
#########################################
def eqw_Eval(x, p):
    """ Compute the sum of the templates """
    y = np.dot(p, x)
    return y

def eqw_myfunct(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = eqw_Eval(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return [status, (y-model)/err]

def eqw_run_mpfit(fobs, efobs, fit_templates, show_results=False):
    import copy
    import unicorn.mpfit as mpfit
    import matplotlib.pyplot as plt
    
    x = fit_templates
    y = fobs
    ey = efobs
    
    #### ignore bad data
    keep = (y > -99) & (ey > 0) & np.isfinite(y) & (y != 0)    
    fa = {'x':x[:,keep], 'y':y[keep], 'err':ey[keep]}
        
    NPARAM = x.shape[0]
    
    p0 = np.ones(NPARAM,dtype='float64')  #initial conditions
    pactual = np.ones(NPARAM) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0, 0], 'limits':[0,0.]}
    parinfo=[]
    for i in range(NPARAM):
        parinfo.append(copy.deepcopy(parbase))

    for i in range(NPARAM): 
        parinfo[i]['value']=p0[i]
    
    #### Fix continuum normalization to unity
    #parinfo[0]['value'] = 1.
    #parinfo[0]['fixed'] = 1
    
    m = mpfit.mpfit(eqw_myfunct, p0, parinfo=parinfo,functkw=fa, quiet=1)
    #print m.params, m.params/m.perror
    
    ### Correlation matrix
    m.pcor = m.covar * 0.
    for i in range(NPARAM):
        for j in range(NPARAM):
            m.pcor[i,j] = m.covar[i,j]/np.sqrt(m.covar[i,i]*m.covar[j,j])
    
    ### relative error
    m.relerr = np.zeros(NPARAM)
    m.relerr[0] = m.perror[0]/m.params[0]
    
    for i in range(1, NPARAM):
        m.relerr[i] = np.sqrt((m.perror[0]/m.params[0])**2+(m.perror[i]/m.params[i])**2)
        
    if show_results:
        print 'Show'
        xarr = np.arange(x.shape[1])
        plt.errorbar(xarr, y, yerr=ey, linestyle='None', marker='o', color='blue')
        yfit = eqw_Eval(x, m.params)
        plt.plot(xarr, yfit, color='red', alpha=0.5, linewidth=3)
    
    return m
    
def equivalent_width(root='GOODS-S-24-G141', id=29):
    """ 
    Measure the line equivalent widths from the template fit.
    
    Would also be nice to add line fluxes themselves.  This should come instantly from
    the coefficients of the line templates, since the templates are normalized to 
    area unity.  Therefore the flux in "mag 25" Jy is just 'coeff/tnorm' or something
    like that.  And then f_lambda fluxes is just a unit game and is probably already
    in 'eazyPy'.
    
    Another idea is to subtract the continuum from the spectrum and integrate the 
    line directly. 
    
    USAGE:
    object, z_grism, halpha_eqw, halpha_flux, oiii_eqw, oiii_flux, hbeta_eqw, hbeta_flux = equivalent_width(root='x', id=1)
    
    """
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import cosmocalc
    import unicorn.analysis
    
    zout = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
    eazy_param = eazy.EazyParam('OUTPUT/%s_%05d.param' %(root, id))
    
    # tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'OUTPUT/line_eqw.tempfilt')
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    coeffs['coeffs'] *= 1./(1+zout.z_peak[0])**2
    
    continuum = np.dot(temp_seds['temp_seds'][:,0:6], coeffs['coeffs'][0:6,0])                             
    
    zpeak_i = tempfilt['zgrid'][coeffs['izbest'][0]]
    lci = tempfilt['lc']
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    fobs = tempfilt['fnu'][:,0]/(lci/5500.)**2
    efobs = tempfilt['efnu'][:,0]/(lci/5500.)**2
    
    idx_ha, idx_hb, idx_oiii = 6,7,8
    for i,temp in enumerate(eazy_param.templates):
        if temp == 'templates/single_lines/OIII_5006.dat':
            idx_oiii = i
        if temp == 'templates/single_lines/Hb_4861.dat':
            idx_hb = i
        if temp == 'templates/single_lines/Ha_6562.dat':
            idx_ha = i
    
    print 'IDX: ',idx_ha, idx_hb, idx_oiii, tempfilt['NTEMP']
    
    obs_sed_continuum = np.dot(tempfilt['tempfilt'][:,0:idx_ha,coeffs['izbest'][0]],coeffs['coeffs'][0:idx_ha,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_ha = np.dot(tempfilt['tempfilt'][:,idx_ha,coeffs['izbest'][0]],coeffs['coeffs'][idx_ha,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_oiii = np.dot(tempfilt['tempfilt'][:,idx_oiii,coeffs['izbest'][0]],coeffs['coeffs'][idx_oiii,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_hb = np.dot(tempfilt['tempfilt'][:,idx_hb,coeffs['izbest'][0]],coeffs['coeffs'][idx_hb,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    # obs_sed_ha /= obs_sed_ha.max()
    # obs_sed_oiii /= obs_sed_oiii.max()
    # obs_sed_hb /= obs_sed_hb.max()
    
    fit_templates = np.array([obs_sed_continuum, obs_sed_ha, obs_sed_oiii, obs_sed_hb])

    if np.std(obs_sed_ha[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[1,:] *= 0.
    #    
    if np.std(obs_sed_oiii[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[2,:] *= 0.
    #    
    if np.std(obs_sed_oiii[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[3,:] *= 0.
    #    
    
    try:
        mp = unicorn.analysis.eqw_run_mpfit(fobs[is_spec], efobs[is_spec], fit_templates[:, is_spec], show_results=False)
        relerr, perror = mp.relerr, mp.perror
    except:
        relerr, perror = np.zeros(4)-1, np.zeros(4)-1
      
    #print 'relerr', relerr, perror
      
    #mp.relerr, mp.params, mp.perror
    
    #### for testing
    plot = """
        tempfiltx, coeffsx, temp_sedsx, pzx = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        dlam_spec = lci[-1]-lci[-2]
        is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        
        # plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), continuum, color='black')
        # plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), temp_sed, color='red', alpha=0.6)
        plt.plot(lci[is_spec], fobs[is_spec], color='blue', alpha=0.4)
        plt.errorbar(lci[is_spec], fobs[is_spec], efobs[is_spec], color='blue', alpha=0.4)
        plt.plot(lci[is_spec], obs_sed_continuum[is_spec], color='black', alpha=0.4, linewidth=3)
        plt.plot(lci[is_spec], obs_sed_ha[is_spec], color='black', alpha=0.4, linewidth=3)
        # plt.plot(lci[is_spec], obs_sed_manual[is_spec], color='red', alpha=0.4, linewidth=3)
        
        plt.plot(lci[is_spec], obs_sed[is_spec], color='orange', linewidth=3, alpha=0.4)
        plt.semilogx()
        plt.xlim(3000,4.e4)
        
        #### Test getting fluxes in f_lambda
        #np.trapz(temp_seds['temp_seds'][:,6], temp_seds['templam']*coeffs['tnorm'][6])
        ha_flux_coeffs = coeffs['coeffs'][6,0]/coeffs['tnorm'][6]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2*(6563.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
        
        ha_flux_convolved = temp_seds['temp_seds'][:,6].max()*coeffs['coeffs'][6,0]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2
        
        
        ha_flux_integrated = np.trapz(halpha, temp_seds['templam']*(1+zout.z_peak[0]))
        ha_flux_integrated *= (6563.*(1+zout.z_peak[0])/5500.)**2
        ha_flux_integrated *= 10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2
        
        import cosmocalc
        cc = cosmocalc.cosmocalc(zout.z_peak[0], H0=71, WM=0.27, WV=1-0.27)
        lum_ha = ha_flux_integrated*cc['DL_cm']**2*4*np.pi
        sfr_ha = 7.9e-42*lum_ha
        
        plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), halpha)
        plt.xlim(6300*(1+zout.z_peak[0]),6900*(1+zout.z_peak[0]))
        
    """
    
    halpha = temp_seds['temp_seds'][:,idx_ha]*coeffs['coeffs'][idx_ha,0]
    halpha[halpha < 1.e-8*halpha.max()] = 0
    halpha_eqw = -np.trapz((-halpha/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    
    # fp = open('../EQW_FOR_MATTIA/test_temp.dat','w')
    # fp.write('# lam continuum line\n')
    # for i in range(1,len(temp_seds['templam'])):
    #     fp.write('%.6e %.3e %.3e\n' %(temp_seds['templam'][i], continuum[i], halpha[i]))
    # fp.close()    
    
    ## test different ways of measuring the equivalenth width from the convolved template and from the spectrum itself
    halpha_eqw_smooth = -np.trapz((-obs_sed_ha/obs_sed_continuum)[is_spec], lci[is_spec])
    use = is_spec & (obs_sed_ha > 1.e-5*obs_sed_ha.max())
    halpha_eqw_data = -np.trapz((-(tempfilt['fnu'][:,0]/(lci/5500.)**2-obs_sed_continuum)/obs_sed_continuum)[use], lci[use])
    
    # 
    # idx = np.arange(len(obs_sed_continuum))
    # fp = open('../EQW_FOR_MATTIA/test_obs.dat','w')
    # fp.write('# lam fnu continuum line\n')
    # for i in idx[is_spec]:
    #     fp.write('%.6e %.3e %.3e %.3e\n' %(lci[i], tempfilt['fnu'][i,0], obs_sed_continuum[i], obs_sed_ha[i]))
    # fp.close()
    
    halpha_flux = coeffs['coeffs'][idx_ha,0]/coeffs['tnorm'][idx_ha]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2*(6563.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    halpha_err = halpha_eqw*relerr[1]
    if perror[1] == 0:
        halpha_err = -1.
    
    #print 'relerr:' %relerr
    #print 'EQW, FLUX, ERR: %.3f %.3f %.3f' %(halpha_eqw, halpha_err, halpha_flux)
    halpha_eqw = (halpha_eqw, halpha_eqw_smooth, halpha_eqw_data)
    
    # cc = cosmocalc.cosmocalc(zout.z_peak[0], H0=71, WM=0.27, WV=1-0.27)
    # halpha_lum = halpha_flux*cc['DL_cm']**2*4*np.pi
    # halpha_sfr = 7.9e-42*halpha_lum
    
    oiii = temp_seds['temp_seds'][:,idx_oiii]*coeffs['coeffs'][idx_oiii,0]
    oiii[oiii < 1.e-8*oiii.max()] = 0
    oiii_eqw = -np.trapz((-oiii/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    oiii_flux = coeffs['coeffs'][idx_oiii,0]/coeffs['tnorm'][idx_oiii]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(5007.*(1+zout.z_peak[0]))**2*(5007.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    oiii_err = oiii_eqw*relerr[2]
    if perror[2] == 0:
        oiii_err = -1.
    
    hbeta =  temp_seds['temp_seds'][:,idx_hb]*coeffs['coeffs'][idx_hb,0]
    hbeta[hbeta < 1.e-8*hbeta.max()] = 0
    hbeta_eqw = -np.trapz((-hbeta/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    hbeta_flux = coeffs['coeffs'][idx_hb,0]/coeffs['tnorm'][idx_hb]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(4861.*(1+zout.z_peak[0]))**2*(4861.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    hbeta_err = hbeta_eqw*relerr[3]
    if perror[3] == 0:
        hbeta_err = -1.
    
    ###### MPFIT to get errors on normalizations
    
    return '%s_%05d' %(root, id), zout.z_peak[0], halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux
    
def make_o2_templates():
    """
    Leaving the OIII strength fixed allows too much freedom in the fits for galaxies at z=0.5 and only ugriz.  The code likes to put the galaxies at z>0.7 with a huge OII line that mimicks the break.  As a fix, try making a set of the `noline` templates that include the eazy OII line from v1.1
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/templates/')
    import glob
    files=glob.glob('EAZY_v1.0_lines/*nolines.dat')
    if not os.path.exists('O2_ONLY'):
        status = os.system('mkdir O2_ONLY')
    
    for file in files:
        noline = np.loadtxt(file)
        noline_wave, noline_flux = noline[:,0], noline[:,1]
        noline_keep = np.abs(noline_wave-3727) > 30
        
        hasline = np.loadtxt(file.replace('v1.0','v1.1').replace('_nolines',''))
        hasline_wave, hasline_flux = hasline[:,0], hasline[:,1]
        hasline_keep = np.abs(hasline_wave-3727) < 30
        
        xnew = np.append(noline_wave[noline_keep],hasline_wave[hasline_keep])
        ynew = np.append(noline_flux[noline_keep],hasline_flux[hasline_keep])
        
        s = np.argsort(xnew)
        xnew, ynew = xnew[s], ynew[s]
        
        fp = open('O2_ONLY/'+os.path.basename(file),'w')
        for i in range(len(xnew)):
            fp.write('%13.5e %13.5e\n' %(xnew[i], ynew[i]))
        
        fp.close()
        
    print "\n\nManually make a copy of 'fit_lines.spectra.param'\n\n"
    
def make_line_templates():
    line_wavelengths = [[6562.800], [5006.843, 4958.911], [4862.68], [3727.0], [6718.29, 6732.67], [1216.], [9068.6, 9530.6], [2799.117], [1549.48], [2326.0], [4341.68], [3889.0]]
    line_ratios = [[1], [2.98, 1], [1], [1], [1, 1], [1], [1, 2.44], [1], [1], [1], [1], [1]]
    line_names = ['Ha', 'OIII', 'Hb', 'OII', 'SII', 'Lya', 'SIII', 'MgII', 'CIV', 'CII','Hg','HeI']
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    NLINE = len(line_names)
    fp = open('templates/EAZY_v1.1_lines/lambda_v1.1.def')
    xspec = np.cast[float](fp.readlines())
    fp.close()
    
    NDEF = len(xspec)
    vel_width = 300 # km/s
    
    fp = open('templates/eazy_v1.0_nolines.spectra.param')
    spec_list = fp.readlines()
    fp.close()
    
    last = np.int(spec_list[-1][0])
    for i in range(NLINE):
        lw = line_wavelengths[i][0]
        yspec = xspec*0.
        ratio = np.array(line_ratios[i])*1./np.sum(line_ratios[i])
        for j in range(len(line_wavelengths[i])):
            lj = line_wavelengths[i][j]
            sigma = vel_width/3.e5*lj
            yspec += ratio[j]*1./np.sqrt(2*np.pi*sigma**2)*np.exp(-1.*(xspec-lj)**2/2/sigma**2)
            
        yspec /= np.trapz(yspec, xspec)
        
        #
        fp=open('templates/single_lines/'+line_names[i]+'_%0d.dat' %(lw), 'w')
        for j in range(NDEF):
            fp.write(' %.5e %.5e\n' %(xspec[j], yspec[j]+1.e-12))
        #
        fp.close()
        #
        spec_list.append('%d   templates/single_lines/%s_%0d.dat 1.0 0 1.0\n' %(last+i+1, line_names[i], lw))
    
    fp = open('templates/fit_lines.spectra.param','w')
    fp.writelines(spec_list)
    fp.close()
    
def test_equivalent_widths():
    """
    Mattia pointed out that he gets equivalent widths somewhat larger than the catalog values when he measures them with his own code.  
    """
    
    object = 'COSMOS-3-G141_00765'
    ### EW(Matt) =75.043951 \pm 3.8351794, EW(Gabe) = 47.080000 \pm 2.4630000
    unicorn.analysis.run_eazy_fit(root=object.split('_')[0], compress=0.7, TILT_ORDER=1, OLD_RES='FILTER.RES.v9.R300', zmin=0.7, zmax=1.2, id=int(object.split('_')[1]), force_zrange=True, COMPUTE_TILT=True)
    
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='../EQW_FOR_MATTIA/')
    
    obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
    
    os.chdir('../REDSHIFT_FITS_v1.6')
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='../EQW_FOR_MATTIA/')
    
    obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
    
    ### Plot:
    obs_sed = catIO.Readfile(object+'_obs_sed.dat')
    
    lci = obs_sed.lc
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)

    idx = np.arange(len(obs_sed.lc))
    fp = open('junk.txt','w')
    for i in idx[is_spec]:
        fp.write('%.1f %.3e\n' %(obs_sed.lc[i], obs_sed.obs_sed[i]*(5500./obs_sed.lc[i])**2 ))
    fp.close()
    
    temp_sed = catIO.Readfile(object+'_temp_sed.dat')
    idx = np.arange(len(temp_sed.lam))
    keep = (temp_sed.lam > 1.e4) & (temp_sed.lam < 1.6e4)
    fp = open('temp.txt','w')
    for i in idx[keep]:
        fp.write('%.3f %.3e\n' %(temp_sed.lam[i], temp_sed.fnu_temp[i]*(5500./temp_sed.lam[i])**2 ))
    fp.close()
    
    plt.plot(lci[is_spec], obs_sed.fnu[is_spec]*(5500./lci[is_spec])**2, color='black')
    plt.plot(lci[is_spec], obs_sed.obs_sed[is_spec]*(5500./lci[is_spec])**2, color='red')
    plt.plot(lci[is_spec], obs_sed.continuum[is_spec]*(5500./lci[is_spec])**2, color='blue')
    plt.plot(temp_sed.lam, temp_sed.fnu_temp*(5500./temp_sed.lam)**2, color='blue')
    plt.plot(temp_sed.lam, temp_sed.continuum*(5500./temp_sed.lam)**2, color='blue')
    
    tt = catIO.Readfile('test_temp.dat')
    plt.plot(tt.lam, tt.continuum, color='green')
    plt.plot(tt.lam, tt.line, color='green')
    
    oo = catIO.Readfile('test_obs.dat')
    plt.plot(oo.lam, oo.fnu/(oo.lam/5500.)**2, color='orange')
    plt.plot(oo.lam, oo.continuum, color='purple')
    
    plt.xlim(1.1e4,1.7e4)
    plt.ylim(0,100)
    
    np.trapz(-(obs_sed.line/obs_sed.continuum)[is_spec], obs_sed.lc[is_spec])
    keep = (temp_sed.lam > 1.2e4) & (temp_sed.lam < 1.3e4)
    np.trapz(-(temp_sed.line/temp_sed.continuum)[keep], temp_sed.lam[keep])
    
    keep = (tt.lam > 1.2e4) & (tt.lam < 1.3e4)
    keep = (tt.lam > 0.6e4) & (tt.lam < 0.7e4)
    np.trapz(-(tt.line/tt.continuum)[keep], tt.lam[keep]*(1+0.9341))
    
    ############ Mattia's list synced from the v1.6 redshift fits
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/EQW_FOR_MATTIA')
    ew = catIO.Readfile('list_for_gabe.dat')
    
    fp = open('ASCII/gabe_eqw.dat','w')
    fp.write('# id z_grism halpha_eqw halpha_eqw_err\n')
    
    for i, object in enumerate(ew.id):
        try:
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
        except:
            print 'xx %s xx' %(object)
            pass
        #
        fp.write('%s %.3f %.3f\n' %(object, halpha_eqw[0], halpha_err))
        #
        unicorn.analysis.make_eazy_asciifiles(object=object, savepath='ASCII')
    
    fp.close()
        
    new_fits = np.zeros((3, len(ew.id)))
    z_fit = np.zeros(len(ew.id))
    
    for i, object in enumerate(ew.id):
        print object
        try:
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
            z_fit[i] = z_grism
        except:
            continue
        #        
        new_fits[:,i] = np.array(halpha_eqw)
    
    #### Run splot to measure them by hand
    import iraf
    from iraf import onedspec
    for i, object in enumerate(ew.id[8:]):
        try:
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
        except:
            continue
        #
        #
        dlam_spec = lci[-1]-lci[-2]
        is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        fpf = open('SPLOT/%s_specfit.txt' %(object),'w')
        fpf.write('%s\n' %object)
        fpf.write('%d 5100.\n' %(len(lci[is_spec])))
        #
        fpd = open('SPLOT/%s_data.txt' %(object),'w')
        fpm = open('SPLOT/%s_model.txt' %(object),'w')
        for lam, flux, model, flux_err in zip(lci[is_spec], fobs[is_spec], obs_sed[is_spec], efobs[is_spec]):
            fpf.write('%.6e %.2e %.2e\n' %(lam/(1+z_grism), flux, flux_err))
            fpd.write('%.6e %.2e\n' %(lam, flux))
            fpm.write('%.6e %.2e\n' %(lam, model))
        fpf.close()
        fpd.close()
        fpm.close()
        #
        iraf.rspectext('SPLOT/%s_data.txt' %(object),'SPLOT/%s_data.fits' %(object), dtype="interp", crval1=11000.0, cdelt1 = 22)
        iraf.rspectext('SPLOT/%s_model.txt' %(object),'SPLOT/%s_model.fits' %(object), dtype="interp", crval1=11000.0, cdelt1 = 22)
        #
        print 'z=%.3f, Ha=%.1f' %(z_grism, 6563*(1+z_grism))
        print halpha_eqw, ew.ew_matt[i]
        #
        iraf.splot('SPLOT/%s_data.fits' %(object), save_file='splot.log')
        iraf.splot('SPLOT/%s_model.fits' %(object), save_file='splot.log')
    
    splot_k_data = ew.ew_gabe*0.
    splot_e_data = ew.ew_gabe*0.
    splot_k_model = ew.ew_gabe*0.
    splot_e_model = ew.ew_gabe*0.
    fp = open('splot.log.save')
    lines = fp.readlines()
    N = len(lines)/12
    for i in range(N):
        group = lines[i*12:i*12+12]
        object = group[1].split('DATA/')[1].split('_data')[0]
        splot_k_data[ew.id == object] = -float(group[3].split()[3])
        splot_e_data[ew.id == object] = -float(group[5].split()[3])
        splot_k_model[ew.id == object] = -float(group[9].split()[3])
        splot_e_model[ew.id == object] = -float(group[11].split()[3])
        
    plt.plot(ew.ew_matt, ew.ew_gabe, marker='o', linestyle='None', color='black', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[0,:], marker='o', linestyle='None', color='red', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[1,:], marker='o', linestyle='None', color='blue', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[2,:], marker='o', linestyle='None', color='green', alpha=0.3)
    
    #plt.plot(splot_k_data, splot_k_model, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_k_data, splot_k_model, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_k_data, ew.ew_gabe, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_e_data, splot_e_model, marker='o', linestyle='None', color='black', alpha=0.8)
    # plt.plot(splot_e_data, new_fits[0,:]/1.07, marker='o', linestyle='None', color='orange', alpha=0.8)
    
    #plt.plot(splot_e_model, splot_e_data, marker='o', linestyle='None', color='black', alpha=0.8)
    
    #plt.plot(splot_e_data, new_fits[2,:], marker='o', linestyle='None', color='red', alpha=0.8)

    #plt.plot(splot_e_data, splot_k_data, marker='o', linestyle='None', color='black', alpha=0.8, ms=10)

    fig = unicorn.catalogs.plot_init(square=True, aspect=1./3, left=0.12, xs=12)
        
    ax = fig.add_subplot(131)
    
    #plt.plot(splot_e_data, new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_gabe*(1+z_fit), marker='o', linestyle='None', color='blue', alpha=0.8, ms=6)
    #plt.plot(splot_e_data, new_fits[2,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_matt, marker='o', linestyle='None', color='red', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW data, splot "e"')
    plt.ylabel('EQW meas, Mattia [red], Gabe/old*(1+z) [blue]')
    
    ax = fig.add_subplot(132)
    
    plt.plot(splot_e_data, new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    #plt.plot(splot_e_data, ew.ew_gabe*(1.9), marker='o', linestyle='None', color='blue', alpha=0.8, ms=10)
    #plt.plot(splot_e_data, new_fits[2,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_matt, marker='o', linestyle='None', color='red', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW data, splot "e"')
    plt.ylabel('EQW meas, Mattia [red], Gabe/fixed [black]')
    
    ax = fig.add_subplot(133)
    plt.plot(new_fits[0,:], new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW, full template')
    plt.ylabel('EQW, convolved')
    
    #fig.savefig('eqw_comparison_old.pdf')
    
    fp = open('gabe_updated.dat','w')
    fp.write('# id z   eqw_matt eqw_gabe  splot_e splot_k  eqw_full eqw_conv\n')
    for i in range(len(ew.id)):
        fp.write('%-20s  %.3f %7.2f %7.2f  %7.2f %7.2f  %7.2f %7.2f\n' %(ew.id[i], z_fit[i],  ew.ew_matt[i], ew.ew_gabe[i], splot_e_data[i], splot_k_data[i], new_fits[0,i], new_fits[1,i]))
    fp.close()
    
    ##### Scatter between measurements
    fig = unicorn.catalogs.plot_init(square=True, aspect=1./2, left=0.12, xs=12*2./3)
        
    ax = fig.add_subplot(121)
    
    #plt.plot(splot_e_data, ew.ew_gabe*(1+z_fit), marker='o', linestyle='None', color='blue', alpha=0.8, ms=6)
    plt.plot(new_fits[1,:], ew.ew_matt, marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW, Gabe / fixed')
    plt.ylabel('EQW, Mattia')
    
    ax = fig.add_subplot(122)
    
    delta = np.log10(new_fits[1,:]/ew.ew_matt)
    keep = np.isfinite(delta)
    ax.hist(delta, bins=25, range=(-1,1))
    ax.set_xlabel(r'log EW(GABE)/EW(Mattia)')
    ax.set_ylabel(r'N')
    
    biw = threedhst.utils.biweight(delta[keep], both=True)
    nm = threedhst.utils.nmad(delta[keep])
    ax.text(-0.8,14,r'$(\mu,\ \sigma) = %.2f,\ %.2f$' %(biw[0], biw[1]))
    
    fig.savefig('eqw_scatter.pdf')

def find_brown_dwarfs():
    import glob
    import threedhst.catIO as catIO
    
    obj = catIO.Readfile('AEGIS-3-G141_00177.dat')
    

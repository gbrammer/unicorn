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

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

HAS_PHOTOMETRY = True
PHOTOMETRY_ID = None
BAD_SPECTRUM = False

SPC_FILENAME = None
SPC = None

IREF = os.getenv('iref')

try:
    flat_f140 = pyfits.open(IREF+'/uc721143i_pfl.fits')
    flat_g141 = pyfits.open(IREF+'/u4m1335mi_pfl.fits')
    flat = flat_g141[1].data[5:1019,5:1019] / flat_f140[1].data[5:1019, 5:1019]
    flat[flat <= 0] = 5
    flat[flat > 5] = 5
except:
    print '\nthreedhst.grism_sky: Flat-field files (uc721143i_pfl.fits) not found in IREF: %s\n' %(IREF)
    flat = np.ones((1014,1014))

xprofile = None
yprofile = None

try:
    not_blobs = (pyfits.open(unicorn.GRISM_HOME+'COSMOS/RAW/ibhm29wnq_flt.fits.gz')[3].data & (4+32+16+512)) == 0
except:
    not_blobs = np.ones((1014,1014))
    
def profile_msk(flt='ibhm46ioq_msk.fits', biweight=False, extension=0):
    """
    Get a cut across the columns of a FLT image, optionally masking objects and DQ 
    pixels.  
        
    If `biweight`, then the output is the biweight mean of each column.  
    Otherwise, it's just the mean.
    """
    import threedhst.grism_sky as bg
    
    im = pyfits.open(flt)
                    
    mask = (im[extension].data == 0) | (~not_blobs)
    
    shp = im[extension].data.shape
    xpix = np.arange(shp[0])
    
    N = np.ones(shp)
    
    im[extension].data[mask] = 0
    N[mask] = 0
    
    ypix = np.sum(im[extension].data, axis=0) / np.sum(N, axis=0)
    if biweight:
        for i in range(shp[0]):
            column = im[extension].data[:,i]
            ypix[i] = threedhst.utils.biweight(column[column != 0], mean=True)
    #
    bg.xprofile, bg.yprofile = xpix, ypix
    
    return xpix, ypix #, ylo, yhi

def show_profile():
    import unicorn.background as bg
    """
    Look at images collapsed along columns to separate into groups with
    similar patterns.
    """
    
    #### All fields
    
    # flt_files = []
    # for field in ['COSMOS','AEGIS','GOODS-N','GOODS-S','UDS']:
    #     fp = open('%s.G141.list' %(field))
    #     flt_files.extend(fp.readlines())
    #     fp.close()
    #
    bg_flt, bg_field, bg_val = np.loadtxt('background.%s.dat' %(filter), dtype=np.str, unpack=True)
    flt_files = bg_flt
    
    N = len(flt_files)
    profiles = np.zeros((N, 1014))
    for i,flt in enumerate(flt_files):
        print flt, i
        xi, yi = bg.profile_msk(flt=flt)
        profiles[i,:] = yi
    
    #
    norm = np.ones(N)
    test = norm > 0
    for i in range(N):
        yi = profiles[i,:]
        test[i] = True
        norm[i] = np.mean(yi[np.abs(xi-507) < 50])
        if test[i]:
            p = plt.plot(xi,yi/norm[i],color='red', alpha=0.05)
        else:
            norm[i] = 0
    
    #### Scale msk images
    for msk, ni in zip(bg_flt, norm):
        im = pyfits.open(msk)
        im[0].data *= 1./ni
        im.writeto(msk.replace('msk','msk.s'), clobber=True)
        print msk
    
    #### Get some normalization points to sort into types
    xint = np.array([50,170, 880, 940])
    yint = np.ones((N,len(xint)))
    dx = 5
    for i in range(N):
        yi = profiles[i,:]
        for j in range(len(xint)):
            yint[i,j] = np.median(yi[xint[j]-dx:xint[j]+dx])/norm[i]
    
    r1 = yint[:,1]/yint[:,0]
    r2 = yint[:,2]/yint[:,3]
    plt.hist(r1, range=(1, 1.4), bins=40, alpha=0.7)
    plt.hist(r2, range=(0.8, 1.4), bins=60, alpha=0.7)
    
    #### Loop through to make subset images
    so = np.argsort(r1)
    NSET = 5
    
    #### Get extremes
    so = np.argsort(yint[:,1]/norm[i])
    NSET, i = 25, 24
    # i = 23,24
    
    so = np.argsort(yint[:,2]/norm[i])
    NSET, i = 120, 119
    
    os.system('rm *.set*.*')

    for i in range(NSET):
        NI = len(so)/NSET
        if i == (NSET-1):
            max = len(so)
        else:
            max = (i+1)*NI
        #
        idx = so[i*NI:max]
        root = 'set%03d' %(i+1)
        unicorn.background.combine_subset(filter='G141', idx=idx, root=root)
        #
        for j in idx:
            p = plt.plot(xi, profiles[j,:]/norm[j], color='red', alpha=0.05)
        #
        xi, yi = bg.profile_msk(flt='combine.%s.%s.fits' %(filter, root))
        p = plt.plot(xi, yi, color=(0,0,(i+1.)/NSET))
    
    ##### Show earlier sky files
    sky_files = glob.glob('../CONF/sky*.fits')
    for sky_file in sky_files:
        xi, yi = bg.profile_msk(flt=sky_file)
        p = plt.plot(xi, yi, color='green')
    
    plt.xlim(-10,1024)
    plt.ylim(0.8,1.1)
    
    #### Use these in "make_imaging_flat" below
    root = 'set1'

    idx = so[1*NI:2*NI]
    root = 'set2'

    idx = so[2*NI:3*NI]
    root = 'set3'
    
    idx = so[3*NI:]
    root = 'set4'
        
    #for i in range(len(flt_files)):
    #    flt_files[i] = flt_files[i][:-1].replace('msk','flt')
    
    N = len(flt_files)
    profiles = np.zeros((N, 1014))
    for i,flt in enumerate(flt_files):
        print flt
        xi, yi = bg.profile_msk(flt=flt[:-1])
        profiles[i,:] = yi
    
    #### COSMOS
    fp = open('COSMOS.G141.list')
    flt_files = fp.readlines()
    fp.close()
    #for i in range(len(flt_files)):
    #    flt_files[i] = flt_files[i][:-1].replace('msk','flt')
    
    N = len(flt_files)
    profiles = np.zeros((N, 1014))
    for i,flt in enumerate(flt_files):
        print flt
        xi, yi = bg.profile_msk(flt=flt[:-1])
        profiles[i,:] = yi
    
    norm = np.zeros(N)
    test = norm > 0
    for i in range(N):
        yi = profiles[i,:]
        norm[i] = np.mean(yi[np.abs(xi-507) < 50])
        test[i] = np.median(yi[np.abs(xi-40) < 10]/norm[i]) < 1.95
        if test[i]:
            p = plt.plot(xi,yi/norm[i],color=(norm[i]/3.3,0,0), alpha=0.1)
        else:
            norm[i] = 0
            
    profiles_norm = profiles / np.dot(norm.reshape(N,1), np.ones((1,1014)))
    avg = np.mean(profiles_norm[norm != 0, :], axis=0)
    plt.plot(xi, avg, color='blue', alpha=0.5)
    
    # for i in range(N):
    #     yi = profiles[i,:]*1.
    #     if yi.sum() == 0:
    #         continue
    #     #
    #     yi-=0.
    #     nor = np.mean(yi[np.abs(xi-307) < 50])
    #     p = plt.plot(xi,yi/nor,color=(norm[i]/6,0,0), alpha=0.1)
    
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    plt.savefig('COSMOS_profile.png')
        
    #### GOODS-N
    fp = open('GOODS-N.G141.list')
    flt_files = fp.readlines()
    fp.close()
    
    Ng = len(flt_files)
    
    profiles_g = np.zeros((Ng, 1014))
    for i,flt in enumerate(flt_files):
        print flt
        xi, yi = bg.profile_msk(flt=flt)
        profiles_g[i,:] = yi
    
    xi = np.arange(1014)
    norm_g = np.zeros(Ng)
    test = norm_g > 0
    for i in range(Ng):
        yi = profiles_g[i,:]
        norm_g[i] = np.mean(yi[np.abs(xi-507) < 50])
        # very hi
        test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_g[i]) > 1.02
        # lo
        #test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_g[i]) < 1.01
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_g[i]) > 0.96)
        # hi
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_g[i]) < 0.96)
        #
        test[i] = True
        if test[i]:
            p = plt.plot(xi,yi/norm_g[i],color=(0,0,norm_g[i]/1.8), alpha=0.1)
        else:
            norm_g[i]*=0
    #
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    
    profiles_norm_g = profiles_g / np.dot(norm_g.reshape(Ng,1), np.ones((1,1014)))
    avg_g = np.mean(profiles_norm_g[norm_g != 0, :], axis=0)
    plt.plot(xi, avg_g, color='green', alpha=0.5)

    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    plt.savefig('GOODS-N_profile.png')
    
    #### AEGIS
    fp = open('AEGIS.G141.list')
    flt_files = fp.readlines()
    fp.close()

    Na = len(flt_files)
    
    profiles_a = np.zeros((Na, 1014))
    for i,flt in enumerate(flt_files):
        print flt
        xi, yi = bg.profile_msk(flt=flt)
        profiles_a[i,:] = yi
    
    xi = np.arange(1014)
    norm_a = np.zeros(Na)
    test = norm_a > 0
    for i in range(Na):
        yi = profiles_a[i,:]
        norm_a[i] = np.mean(yi[np.abs(xi-507) < 50])
        test[i] = True
        # very hi
        # test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_a[i]) < 1.52
        # lo
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) > 0.96)
        # hi
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) < 0.96)
        #
        if test[i]:
            p = plt.plot(xi,yi/norm_a[i],color=(0,0,norm_a[i]/1.8), alpha=0.1)
        else:
            norm_a[i]*=0
    #
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    plt.savefig('AEGIS_profile.png')
    
    #### GOODS-S
    fp = open('GOODS-S.G141.list')
    flt_files = fp.readlines()
    fp.close()

    Ngs = len(flt_files)
    
    profiles_gs = np.zeros((Ngs, 1014))
    for i,flt in enumerate(flt_files):
        print flt
        xi, yi = bg.profile_msk(flt=flt)
        profiles_gs[i,:] = yi
    
    xi = np.arange(1014)
    norm_gs = np.zeros(Ngs)
    test = norm_gs > 0
    for i in range(Ngs):
        yi = profiles_gs[i,:]
        norm_gs[i] = np.mean(yi[np.abs(xi-507) < 50])
        test[i] = True
        # very hi
        # test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_a[i]) < 1.52
        # lo
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) > 0.96)
        # hi
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) < 0.96)
        #
        if test[i]:
            p = plt.plot(xi,yi/norm_gs[i],color=(0,0,norm_gs[i]/1.8), alpha=0.1)
        else:
            norm_gs[i]*=0
    #
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    
    plt.savefig('GOODS-S_profile.png')
    
def make_g141_bg():
    """
    Make average background images with object masks
    """
    from pyraf import iraf

    os.chdir("/3DHST/Spectra/Work/Background")

    field = 'COSMOS'
    
    PATHS = []
    files = []
    
    for field in ['COSMOS','GOODS-N','GOODS-S','AEGIS','UDS']:
        info = catIO.Readfile('/3DHST/Spectra/Work/%s/PREP_FLT/files.info' %(field))
        field_files = info.file[info.filter == 'G141']
        files.extend(field_files)
        PATHS.extend(['/3DHST/Spectra/Work/%s/RAW/' %(field)] * len(info.file[info.filter == 'G141']))
        
    field = 'ALL'
        
    #files = glob.glob('ibhm*flt.seg.fits')
    #PATH = ('/3DHST/Spectra/Work/%s/RAW/' %(field))*len(files)
        
    # #### Direct flat-field
    flat = flat_g141[1].data[5:1019,5:1019] / pyfits.open('COSMOS_f140w_flat.fits')[1].data[5:-5,5:-5]
    flat[flat <= 0] = 5
    flat[flat > 5] = 5
    
    NF = len(files)
    idx = np.arange(NF)
    nxpix, nypix = 1014, 1014
    
    #nxpix, nypix = 507, 507
    
    X = np.zeros((NF, nxpix*nypix))
    
    ## Otherwise get it from "show_profile" above
    test = idx > -10
    
    for j,i in enumerate(idx):
        if ~test[i]:
            continue
        #
        fi = files[i]
        if not os.path.exists(fi.replace('flt','flt.seg')):
            continue
        #    
        if os.path.exists(fi.replace('.gz','')+'.mask.reg'):
            continue
        #
        flt = pyfits.open(PATHS[i]+files[i])
        flt[1].data *= flat
        print unicorn.noNewLine+'%d %s %s' %(i, files[i], flt[0].header['PFLTFILE'])
        #
        ### Segmentation mask
        masked = pyfits.open(fi.replace('flt','flt.seg'))[0].data == 0
        ### DQ mask, hot pixels and the "death star"
        dq_ok = (flt[3].data & (4+32+16)) == 0
        #
        ok = masked & np.isfinite(flt[1].data) & (dq_ok)
        #flt[1].data /= np.median(flt[1].data[ok])
        flt[1].data /= threedhst.utils.biweight(flt[1].data[ok], mean=True)
        flt[1].data[(ok == False)] = 0
        X[j,:] = flt[1].data[0:nypix, 0:nxpix].flatten()
        #
        #pyfits.writeto(files[i].replace('flt','msk').replace('.gz',''), flt[1].data, clobber=True, header=flt[1].header)
    
    #### Average
    #nsum = np.sum(X != 0, axis=0).reshape(1014,1014)
    #avg = np.sum(X, axis=0).reshape(1014,1014)/nsum
    
    for field in ['COSMOS','GOODS-N','GOODS-S','AEGIS','UDS']:
        info = catIO.Readfile('/3DHST/Spectra/Work/%s/PREP_FLT/files.info' %(field))
        field_files = info.file[info.filter == 'G141']
        fp = open(field+'.g141.list','w')
        for ff in field_files:
            msk = ff.replace('flt.fits.gz','msk.fits')
            if os.path.exists(msk):
                fp.write('%s\n' %(msk))
        fp.close()
        #
        iraf.imcombine ( input = '@%s.g141.list' %(field), output = 'combined_g141_%s' %(field), 
           headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
           expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
           reject = 'minmax', project = iraf.no, outtype = 'real', 
           outlimits = '', offsets = 'none', masktype = 'none', 
           maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
           weight = 'none', statsec = '', expname = '', lthreshold = 0.02, 
           hthreshold = 20.0, nlow = 3, nhigh = 3, nkeep = 1, 
           mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
           gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
           
        
    fp = open('msk_list','w')
    for file in files:
        fp.write(file+'\n')
    fp.close()
    
    iraf.imcombine ( input = '@msk_list', output = 'combine_masked', 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
       expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'real', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = 'none', statsec = '', expname = '', lthreshold = 1e-06, 
       hthreshold = 100.0, nlow = 5, nhigh = 5, nkeep = 1, 
       mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
      
    sky = pyfits.open('combine_COSMOS.fits')[0].data
    
    # #### Average
    # nsum = np.sum(X != 0, axis=0).reshape(nypix,nxpix)
    # avg = np.sum(X, axis=0).reshape(nypix,nxpix)/nsum
    #     
    # ### Fill empty pixels with no input images
    # sky = avg
    x,y = np.where((np.isfinite(sky) == False) | (sky == 0))
    NX = len(x)
    pad = 1
    for i in range(NX):
        xi = x[i]
        yi = y[i]
        sub = sky[xi-pad:xi+pad+2,yi-pad:yi+pad+2]
        if (np.sum(sub) != 0.0):
            sky[xi,yi] = np.median(sub[np.isfinite(sub)])
    
    still_bad = (np.isfinite(sky) == False) | (sky <= 0.01)
    sky[still_bad] = flat[0:nypix, 0:nxpix][still_bad]
    
    # bad_flat = (flat < 0.5)
    # sky[bad_flat] = flat[bad_flat]
        
    im_sky = pyfits.PrimaryHDU(data=sky)
    im_n = pyfits.ImageHDU(data=nsum)
    im = pyfits.HDUList([im_sky, im_n])
    im.writeto('sky.fits', clobber=True)
    
    #### for DIRECT flat
    flatim = pyfits.open('/3DHST/Spectra/Work/CONF/sky_cosmos.fits')
    flatim[0].data = sky
    flatim[1].data = sky
    #flatim[3].data[5:-5,5:-5] = nsum
    flatim.writeto('%s_g141_flat.fits' %(field), clobber=True)


def make_imaging_flat():
    """
    Make average background images with object masks
    """
    from pyraf import iraf
    
    #files = glob.glob('ibhm*flt.seg.fits')
    #PATH = ('/3DHST/Spectra/Work/%s/RAW/' %(field))*len(files)
        
    ###################### Grism sky backgrounds
    filter, flat_file = 'G141', 'u4m1335mi_pfl.fits'
    
    flat = pyfits.open(IREF+'/'+flat_file)[1].data[5:-5,5:-5] / pyfits.open(IREF+'/flat.IR_avg.fits')[1].data[5:-5,5:-5]
    flat[flat <= 0] = 5
    flat[flat > 5] = 5
    
    ##################### Direct flat-field
    filter, flat_file = 'F140W', 'uc721143i_pfl.fits'

    filter, flat_file = 'F125W', 'uc72113qi_pfl.fits'

    filter, flat_file = 'F160W', 'uc721145i_pfl.fits'

    filter, flat_file = 'F105W', 'uc72113oi_pfl.fits'
        
    flat = pyfits.open(IREF+'/'+flat_file)[1].data[5:-5,5:-5]
    
    flat[flat <= 0] = 5
    flat[flat > 5] = 5
    
    ############### 3D-HST
    os.chdir("/3DHST/Spectra/Work/Background")
    fields = ['COSMOS','GOODS-N','GOODS-S','AEGIS','UDS']
    PREP_FLT = '/3DHST/Spectra/Work/xxx/PREP_FLT/'
    RAW = '/3DHST/Spectra/Work/xxx/RAW/'
    
    ############### CANDELS
    os.chdir('/Users/gbrammer/CANDELS/Flats/')
    fields = ['GOODS-S','EGS','UDS']
    PREP_FLT = '/Users/gbrammer/CANDELS/xxx/PREP_FLT/'
    RAW = '/Users/gbrammer/CANDELS/xxx/RAW/'
        
    PATHS = []
    files = []
    file_field = []
    
    for field in fields:
        info = catIO.Readfile(PREP_FLT.replace('xxx',field)+'files.info')
        field_files = info.file[info.filter == filter]
        files.extend(field_files)
        PATHS.extend([RAW.replace('xxx',field)] * len(field_files))
        file_field.extend([field]*len(field_files))
        
    ##################
    
    
    NF = len(files)
    idx = np.arange(NF)
    
    ## Otherwise get it from "show_profile" above
    test = idx > -10
    
    fp = open('background.%s.dat' %(filter),'w')
    
    for j,i in enumerate(idx):
        if ~test[i]:
            continue
        #
        fi = files[i]
        if not os.path.exists(fi.replace('flt','flt.seg')):
            continue
        #    
        if os.path.exists(fi.replace('.gz','')+'.mask.reg'):
            continue
        #
        flt = pyfits.open(PATHS[i]+files[i])
        flt[1].data *= flat
        print unicorn.noNewLine+'%d %s %s' %(i, files[i], flt[0].header['PFLTFILE'])
        #
        ### Segmentation mask
        masked = pyfits.open(fi.replace('flt','flt.seg'))[0].data == 0
        ### DQ mask, hot pixels and the "death star"
        dq_ok = (flt[3].data & (4+32+16)) == 0
        #
        ok = masked & np.isfinite(flt[1].data) & (dq_ok)
        #flt[1].data /= np.median(flt[1].data[ok])
        level = threedhst.utils.biweight(flt[1].data[ok], mean=True)
        fp.write('%s %s %.3f\n' %(files[i].replace('flt','msk').replace('.gz',''), file_field[i], level))
        #
        #flt[1].data /= level
        #flt[1].data[(ok == False)] = 0
        #pyfits.writeto(files[i].replace('flt','msk').replace('.gz',''), flt[1].data, clobber=True, header=flt[1].header)
    
    fp.close()  ## background.dat
    
    #
    # nsum = np.sum(X != 0, axis=0).reshape(1014,1014)
    # avg = np.sum(X, axis=0).reshape(1014,1014)/nsum
    # sky = avg
    
    #### Use iraf.imcombine    
    for field in fields:
        info = catIO.Readfile(PREP_FLT.replace('xxx',field)+'files.info')
        field_files = info.file[info.filter == filter]
        if len(field_files) < 10:
            continue
        #
        fp = open('%s.%s.list' %(field, filter),'w')
        for ff in field_files:
            msk = ff.replace('flt.fits.gz','msk.fits')
            if os.path.exists(msk):
                fp.write('%s\n' %(msk))
        fp.close()
        #
        iraf.imcombine ( input = '@%s.%s.list' %(field, filter), output = 'combine.%s.%s' %(field, filter), 
           headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
           expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
           reject = 'minmax', project = iraf.no, outtype = 'real', 
           outlimits = '', offsets = 'none', masktype = 'none', 
           maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
           weight = 'none', statsec = '', expname = '', lthreshold = 1e-06, 
           hthreshold = 100.0, nlow = 5, nhigh = 5, nkeep = 1, 
           mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
           gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
    
    ##### Weight by the square of the background level (more flat signal for higher bg!)
    bg_flt, bg_field, bg = np.loadtxt('background.%s.dat' %(filter), dtype=np.str, unpack=True)
    weights = np.cast[float](bg)**2
           
    fp = open('%s.list' %(filter),'w')
    fpw = open('%s.weight' %(filter),'w')
    for msk, wht in zip(bg_flt, weights):
        if os.path.exists(msk):
            fp.write('%s\n' %(msk))
            fpw.write('%.2f\n' %(wht))
    
    fp.close()
    fpw.close()
    
    iraf.imcombine ( input = '@%s.list' %(filter), output = 'combine.%s' %(filter), 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
       expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'real', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = '@%s.weight' %(filter), statsec = '', expname = '', lthreshold = 1e-06, 
       hthreshold = 100.0, nlow = 5, nhigh = 5, nkeep = 1, 
       mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
        
    ##### Final processing
    combined_files = glob.glob('combine*%s*fits' %(filter))
    for file in combined_files:
        sky = pyfits.open(file)[0].data
        #
        ##### Fix bad pixels
        if filter != 'G141':
            ratio = sky/flat
            stats = threedhst.utils.biweight(ratio[np.isfinite(ratio)], both=True)
            sky = sky/stats[0]
            max = stats[1]*5
        else:
            max = 10
        #
        x,y = np.where((np.isfinite(sky) == False) | (sky/flat > (1+max)) | (sky == 0))
        NX = len(x)
        print '%s: N_fix = %d' %(file, NX)
        pad = 1
        for i in range(NX):
            xi = x[i]
            yi = y[i]
            sub = sky[xi-pad:xi+pad+2,yi-pad:yi+pad+2]
            if (np.sum(sub) != 0.0):
                sky[xi,yi] = np.median(sub[np.isfinite(sub)])
        #
        still_bad = (np.isfinite(sky) == False) | (sky <= 0.01)
        sky[still_bad] = flat[still_bad]
        #            
        #### for DIRECT flat
        if filter == 'G141':
            flatim = pyfits.open(unicorn.GRISM_HOME + 'CONF/sky_cosmos.fits')
            flatim[0].data = sky
            #flatim[3].data[5:-5,5:-5] = nsum
            flatim.writeto(file.replace('combine','sky'), clobber=True)
        else:
            flatim = pyfits.open(IREF+'/'+flat_file)
            flatim[1].data[5:-5,5:-5] = sky
            #flatim[3].data[5:-5,5:-5] = nsum
            flatim.writeto(file.replace('combine','flat'), clobber=True)
            
def combine_subset(filter='G141', idx=np.array([0]), root='set1', use_scaled=True):
    """
    Subset, get index array of objects to use from the "show_profile" function above
    """
    
    from pyraf import iraf
    
    bg_flt, bg_field, bg_val = np.loadtxt('background.%s.dat' %(filter), dtype=np.str, unpack=True)
    weights = np.cast[float](bg_val)**2
           
    fp = open('%s.%s.list' %(filter, root),'w')
    fpw = open('%s.%s.weight' %(filter, root),'w')
    for msk, wht in zip(bg_flt[idx], weights[idx]):
        if os.path.exists(msk):
            if use_scaled:
                img = msk.replace('msk','msk.s')
            else:
                img = msk
            fp.write('%s\n' %(img))
            fpw.write('%.4f\n' %(wht))
    #
    fp.close()
    fpw.close()
        
    iraf.imcombine ( input = '@%s.%s.list' %(filter, root), output = 'combine.%s.%s' %(filter, root), 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
       expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'real', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = '@%s.%s.weight' %(filter, root), statsec = '', expname = '', lthreshold = 1e-04, 
       hthreshold = 100.0, nlow = 2, nhigh = 2, nkeep = 1, 
       mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
            
def make_average_flat_for_grism():
    """
    Take the average of the master flats made from the CANDELS F125W and F160W images
    """
    
    os.chdir('/Users/gbrammer/CANDELS/Flats/')
    f125 = pyfits.open('flat.F125W.fits')
    f160 = pyfits.open('flat.F160W.fits')
    
    avg = f125[1].data*0.5+f160[1].data*0.5
    
    f125[1].data = avg
    f125.writeto('flat.IR_avg.fits', clobber=True)
    
    
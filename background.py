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

def make_g141_bg():
    """
    Make average background images with object masks
    """
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
        print noNewLine+'%d %s %s' %(i, files[i], flt[0].header['PFLTFILE'])
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
        
    #files = glob.glob('ibhm*flt.seg.fits')
    #PATH = ('/3DHST/Spectra/Work/%s/RAW/' %(field))*len(files)
        
    # #### Direct flat-field
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
        print noNewLine+'%d %s %s' %(i, files[i], flt[0].header['PFLTFILE'])
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
           
        
    fp = open('%s.list' %(filter),'w')
    for file in files:
        msk = file.replace('flt.fits.gz','msk.fits')
        if os.path.exists(msk):
            fp.write('%s\n' %(msk))
    fp.close()
    
    iraf.imcombine ( input = '@%s.list' %(filter), output = 'combine.%s' %(filter), 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
       expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'real', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = 'none', statsec = '', expname = '', lthreshold = 1e-06, 
       hthreshold = 100.0, nlow = 5, nhigh = 5, nkeep = 1, 
       mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
    
    #   Weight by the square of the background.  Actually helps to have higher bg levels!
    iraf.imcombine ( input = '@F125W.list', output = 'combine.weighted.F125W', 
       headers = '', bpmasks = '', rejmasks = '', nrejmasks = '', 
       expmasks = '', sigmas = '', logfile = 'STDOUT', combine = 'average', 
       reject = 'minmax', project = iraf.no, outtype = 'real', 
       outlimits = '', offsets = 'none', masktype = 'none', 
       maskvalue = '0', blank = 0.0, scale = 'none', zero = 'none', 
       weight = '@F125W.weight', statsec = '', expname = '', lthreshold = 1e-06, 
       hthreshold = 100.0, nlow = 5, nhigh = 5, nkeep = 1, 
       mclip = iraf.yes, lsigma = 3.0, hsigma = 3.0, rdnoise = '0.', 
       gain = '1.', snoise = '0.', sigscale = 0.1, pclip = -0.5)
    
    ##### Final processing
    combined_files = glob.glob('combine*%s*fits' %(filter))
    for file in combined_files:
        sky = pyfits.open(file)[0].data
        #
        ##### Fix bad pixels
        ratio = sky/flat
        stats = threedhst.utils.biweight(ratio[np.isfinite(ratio)], both=True)
        sky = sky/stats[0]
        max = stats[1]*5
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
        flatim = pyfits.open(IREF+'/'+flat_file)
        flatim[1].data[5:-5,5:-5] = sky
        #flatim[3].data[5:-5,5:-5] = nsum
        flatim.writeto(file.replace('combine','flat'), clobber=True)

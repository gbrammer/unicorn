"""
reduce_scripts.py
 
Reduction scripts for all fields starting with the v4.0 catalogs and images. Previously versions of these scripts resided in reduce.py.
 
"""
 
import os
import glob
import shutil
import time

import numpy as np
import matplotlib.pyplot as plt

USE_PLOT_GUI = False

import pyfits

import time

import threedhst
import unicorn
import unicorn.utils_c as utils_c

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import scipy.ndimage as nd

def interlace_aegis0():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time

    os.chdir(unicorn.GRISM_HOME+'AEGIS/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='AEGIS_F160W', filter='F160W', REFERENCE = '/3DHST/Ancillary/AEGIS/CANDELS/ASTRODRIZZLE/aegis-f160w-astrodrizzle-v4.0_drz_join.fits', SEGM = '/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F160W_seg.fits', Force=True)
        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = 'phot/F160W.cat'
    direct=glob.glob('AEGIS-*-F140W_asn.fits')
    direct.remove('AEGIS-20-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='AEGIS_F160W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('AEGIS-2[0-9]-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
        
    time.strftime('%X %x %Z')
    inter = glob.glob('AEGIS-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduce_scripts.fix_thumbnails(pointing=pointing)
        time.strftime('%X %x %Z')
                
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True
    
    ##### Extract and fit only spec-z objects
    models = glob.glob('AEGIS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = True)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)
                        
    ### Fit all grism redshifts down to F160/F140=24.0:

    models = glob.glob('AEGIS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = False)

    #Fit the emission lines
    models = glob.glob('AEGIS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            gris.new_fit_free_emlines(ztry=None)	

    ### Get some quality flags on the reduced spectra
    root='AEGIS'

    files = glob.glob('AEGIS*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()


def interlace_aegis1():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time

    os.chdir(unicorn.GRISM_HOME+'AEGIS/INTERLACE')

    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = 'phot/F140W.cat'
    direct=glob.glob('AEGIS-1[0-9]-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='AEGIS_F140W'

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('AEGIS-1[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            
    ##### Extract and fit only spec-z objects
    models = glob.glob('AEGIS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = True)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)

    ### Fit all grism redshifts down to F160/F140=24.0:

    models = glob.glob('AEGIS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = True)

    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

    #Fit the emission lines
    models = glob.glob('AEGIS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)	

def interlace_aegis2():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time

    os.chdir(unicorn.GRISM_HOME+'AEGIS/INTERLACE')
        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = 'phot/F140W.cat'
    direct=glob.glob('AEGIS-[2-3][0-9]-F140W_asn.fits')
    direct.remove('AEGIS-20-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='AEGIS_F140W'

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('AEGIS-[2-3][0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    ##### Extract and fit only spec-z objects
    models = glob.glob('AEGIS-[2-3][0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = True)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)

                        
    ### Fit all grism redshifts down to F160/F140=23.0:

    models = glob.glob('AEGIS-[2-3][0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = True)

    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

    #Fit the emission lines
    models = glob.glob('AEGIS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None) 

def interlace_cosmos0():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import time
    import numpy as np
    import unicorn.interlace_test as test
    import unicorn.reduce_scripts

    os.chdir(unicorn.GRISM_HOME+'COSMOS/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='COSMOS_F160W', filter='F160W', 
        REFERENCE = '/3DHST/Ancillary/COSMOS/CANDELS/ASTRODRIZZLE/cosmos-f160w-astrodrizzle-v4.0_drz_join.fits', 
        SEGM = '/3DHST/Photometry/Work/COSMOS/Sex/PSF_matched/F160W_seg.fits')
    
    NGROW=125
    pad=60
    CATALOG='phot/F160W.cat'

    direct=glob.glob('COSMOS-*-F140W_asn.fits')

    extract_limit = 35.
    skip_completed=False
    REF_ROOT='COSMOS_F160W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('COSMOS-*-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    inter = glob.glob('COSMOS-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduce_scripts.fix_thumbnails(pointing=pointing)
        time.strftime('%X %x %Z')

    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='COSMOS-11')

    skip_completed = False
    
    ##### Extract and fit only spec-z objects
    models = glob.glob('COSMOS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = True)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)
                        
    ### Fit all grism redshifts down to F160/F140=24.0:

    models = glob.glob('COSMOS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = False)

    #Fit the emission lines
    models = glob.glob('COSMOS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)
            
def interlace_cosmos1():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import time

    os.chdir(unicorn.GRISM_HOME+'COSMOS/INTERLACE_v4.0')
        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    CATALOG='phot/F160W.cat'

    extract_limit = 35.
    skip_completed=False
    REF_ROOT='COSMOS_F160W'

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('COSMOS-1[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='COSMOS-11')

    skip_completed = True

    models = glob.glob('COSMOS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
            



def interlace_cosmos2():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import time

    os.chdir(unicorn.GRISM_HOME+'COSMOS/INTERLACE_v4.0')
        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    CATALOG='phot/F160W.cat'

    extract_limit = 35.
    skip_completed=False
    REF_ROOT='COSMOS_F140W'

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('COSMOS-2[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='COSMOS-11')

    skip_completed = True

    models = glob.glob('COSMOS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('COSMOS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

    ### Get some quality flags on the reduced spectra
    root='COSMOS'

    files = glob.glob('COSMOS*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()
    
def interlace_goodsn():
    """
    Reduce the GOODS-N pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-N/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(Force=True, REF_ROOT='GOODS-N_F160W', filter='F160W', REFERENCE = '/3DHST/Ancillary/GOODS-N/CANDELS/ASTRODRIZZLE/goods-n-f160w-astrodrizzle-v4.0_drz_join.fits', SEGM = '/3DHST/Photometry/Work/GOODS-N/v4/sextr/checkimages/GOODS-N_F125W_F140W_F160W.seg.fits')
     
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Work/GOODS-N/v4/sextr/catalogs/GOODS-N_F160W.cat'
    
    direct = glob.glob('GOODS-N-*-F140W_asn.fits')

    extract_limit = 24
    skip_completed=True
    REF_ROOT='GOODS-N_F160W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        #
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)
    
    ##### Generate the spectral model and extract spectra
    inter = glob.glob('GOODS-N-4*-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)

    inter = glob.glob('GOODS-N-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduce_scripts.fix_thumbnails(pointing=pointing)
        time.strftime('%X %x %Z')

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-N-11')
        
    skip_completed = False
    

    ##### Extract and fit only spec-z objects
    models = glob.glob('GOODS-N-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = False)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)


    #Fit the redshift for objects brighter than F140W=24.
    models = glob.glob('GOODS-N-4[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = False)


    #Fit the emission lines
    models = glob.glob('GOODS-N-4[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

    ### Get some quality flags on the reduced spectra
    root='GOODS-N'
    
    files = glob.glob('GOODS-N*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])
            
    fp.close()
    fp2.close()
    
    ###################### Analysis plots ################
    ## 
    os.system('paste %s.zfit.dat %s.dqflag.dat > %s.specz.dat' %(root, root, root))
    
    zfit = catIO.Readfile('%s.specz.dat' %(root))
    
    #### v2.0 release
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-N')
    zfit = catIO.Readfile('GOODS-N.zfit.linematched.dat')
    dq = catIO.Readfile('GOODS-N.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm

    #### dq has some ****** IDL format values in max_contam
    lines = pyfits.open('GOODS-N.linefit.fits')
    cat = catIO.Readfile('goodsn_f140w_v1_det_reform.cat')
    root = 'goodsn'
    
    ####  COSMOS
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/COSMOS')
    zfit = catIO.Readfile('COSMOS.zfit.linematched.dat')
    dq = catIO.Readfile('COSMOS.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    lines = pyfits.open('COSMOS.linefit.fits')
    cat = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Catalog/3dhst.cosmos.v2.0.cat')
    cat.x_image, cat.y_image = cat.x, cat.y
    root = 'cosmos'
    fout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Fast/3dhst.cosmos.v2.0.fout')
    zout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Eazy/3dhst.cosmos.v2.0.zout')
    
    ####  GOODS-S
    os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-S')
    zfit = catIO.Readfile('GOODS-S.zfit.linematched.dat')
    dq = catIO.Readfile('GOODS-S.dqflag.linematched.dat')
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    lines = pyfits.open('GOODS-S.linefit.fits')
    cat = catIO.Readfile('GOODS-S_v2.0_PHOTOMETRY/catalog/GOODS-S_v2.0.fullz_wzp.cat')
    cat.x_image, cat.y_image = cat.ximage, cat.yimage
    root = 'goodss'
    
    zfit.mag = dq.mag
    
    dz_spec = (zfit.z_max_spec-zfit.z_spec)/(1+zfit.z_spec)
    dz_peak = (zfit.z_peak_spec-zfit.z_spec)/(1+zfit.z_spec)
    sz = np.maximum((np.abs(dz_spec)/0.001)**(0.8),8)
    
    dz_phot = (zfit.z_peak_phot - zfit.z_spec)/(1+zfit.z_spec)
    
    fig = unicorn.plotting.plot_init(square=True, xs=4, left=0.13)
    ax = fig.add_subplot(111)
    ax.scatter(zfit.z_spec, zfit.z_max_spec, color='blue', alpha=0.2, s=2)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 4)
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{gris}$')
    unicorn.plotting.savefig(fig, root+'_zphot_zspec.png')
    
    z0 = (zfit.z_spec > 0) & (zfit.z_spec < 10)
    plt.scatter(zfit.z_spec, dz_spec, color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_spec[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_z_dz_max.png')
    
    plt.scatter(zfit.z_spec, dz_peak, color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_peak[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{peak}$')
    plt.savefig(root+'_z_dz_peak.png')

    plt.scatter(zfit.z_spec, dz_phot, color='orange', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[z0], dz_phot[z0], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.ylim(-0.1, 0.1)
    plt.xlabel(r'$z_\mathrm{spec}$'); plt.ylabel(r'$\Delta z_\mathrm{phot}$')
    plt.savefig(root+'_z_dz_phot.png')
    
    zp6 = zfit.z_spec > 0.6
    keep = zp6
    
    plt.plot(dq.mag[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(18,24.5)
    plt.plot(dq.mag[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.1,0.1); plt.xlim(18,24.5)
    xm, ym, ys, N = threedhst.utils.runmed(dq.mag[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(dq.mag[keep], dz_phot[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    plt.xlabel(r'dqflag.mag'); plt.ylabel(r'$\Delta z_\mathrm{phot}$')
    plt.savefig(root+'_mag_dz_gris.png')
    keep = keep & (dq.mag < 24)
    
    plt.plot(dq.max_contam[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.02,0.02); plt.xlim(0.002, 30)
    xm, ym, ys, N = threedhst.utils.runmed(dq.max_contam[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.semilogx()
    plt.xlabel(r'dqflag.max_contam'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_max_contam_dz_gris.png')
    keep = keep & (dq.max_contam < 1)

    plt.plot(dq.int_contam[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.008,2)
    xm, ym, ys, N = threedhst.utils.runmed(dq.int_contam[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.semilogx()
    plt.xlabel(r'dqflag.int_contam'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_int_contam_dz_gris.png')
    
    keep = keep & (dq.int_contam < 0.6)
    
    plt.plot(dq.q_z[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.1,0.1); plt.xlim(0.001,100) ; plt.semilogx()
    plt.plot(dq.q_z[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.001,100) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.q_z[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    keep = keep & (zfit.q_z < 5)
    
    plt.plot(dq.f_cover[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_cover[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.1,0.1); plt.xlim(0.01,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_cover[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    plt.xlabel(r'dqflag.f_cover'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_cover_dz_gris.png')
    
    keep = keep & (dq.f_cover > 0.5)
    
    plt.plot(dq.f_flagged[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_flagged[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.01,0.01); plt.xlim(0.005,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_flagged[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')

    plt.xlabel(r'dqflag.f_flagged'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_flagged_dz_gris.png')

    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_negative = np.maximum(dq.f_negative, 1.e-2)*norm
    
    plt.plot(dq.f_negative[z0 & ~keep], dz_spec[z0 & ~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(dq.f_negative[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.01,0.01); plt.xlim(0.005,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(dq.f_negative[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')

    plt.xlabel(r'dqflag.f_negative'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_f_negative_dz_gris.png')
    
    keep = keep & (dq.f_flagged < 0.55)
    nopartial = zfit.f_flagged < 0.2
    
    #### Ha eqw
    z15 = (zfit.z_spec > 0.7) & (zfit.z_spec < 1.5)
    eqw_rf = np.maximum(lines[1].data['HALPHA_EQW'] / (1+lines[1].data['Z']), 0.2)
    plt.plot(eqw_rf[z15], dz_spec[z15], color='blue', alpha=0.4); plt.ylim(-0.05,0.05); plt.xlim(0.1,5000) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(eqw_rf[z15], dz_spec[z15], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(eqw_rf[z15], dz_phot[z15], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    plt.xlabel(r'HALPHA_EQW / (1+z)'); plt.ylabel(r'$\Delta z_\mathrm{gris}$')
    plt.savefig(root+'_ha_eqw_dz_gris_all.png')
    
    plt.xlim(10,1000)
    plt.ylim(-0.01,0.01)
    plt.savefig(root+'_ha_eqw_dz_gris.png')
        
    d99 = ((zout.u99-zout.l99)/(1+zout.z_peak))[zfit.phot_id-1]
    plt.plot(d99[~keep], dz_spec[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    plt.plot(d99[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.01,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(d99[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    #### Do the optimal selection on the full catalog
    fig = unicorn.plotting.plot_init(square=True, xs=6, aspect=0.66, left=0.08)
    ax = fig.add_subplot(111)
    has_gris = zfit.z_max_spec > 0.01
    full_selection = has_gris & (dq.mag < 24) & (dq.max_contam < 1) & (dq.f_cover > 0.5) & (dq.f_flagged < 0.55) & (dq.f_negative < 0.55)
    yh, xh = np.histogram(np.log(1+zfit.z_max_spec[full_selection]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', label='grism', color='blue')
    yh, xh = np.histogram(np.log(1+zfit.z_max_spec[has_gris]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', label='full', color='blue', alpha=0.3)
    yh, xh = np.histogram(np.log(1+zfit.z_spec[full_selection]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, 0-yh, linestyle='steps', marker='None', color='green', label='spec')
    yh, xh = np.histogram(np.log(1+zfit.z_spec[has_gris]), range=(0,np.log(5+1)), bins=np.log(5+1)/0.005)
    a = ax.plot(np.exp(xh[1:])-1, 0-yh, linestyle='steps', marker='None', color='green', label='spec/full', alpha=0.3)
    ax.legend()

    ax.set_ylim(-50,50)
    ax.set_xlim(0,3)
    
    if root=='cosmos':
        ax.set_ylim(-80,80)
        ax.set_xlim(0,3)
        
    ax.set_xlabel('z')
    unicorn.plotting.savefig(fig, root+'_full_zdist_compare_zspec.png')
    
    #### Figure looking for overdensities
    zmax = 3
    r0 = np.median(cat.ra)
    
    width = 0.005
    ztry = 0.850

    zr = (0.7,2.4)
    for count, lz in enumerate(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), width)):
        print count
        ztry = np.exp(lz)-1
        #
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        #
        fig = unicorn.plotting.plot_init(square=True, xs=5, aspect=1/0.7, NO_GUI=False)
        ax = fig.add_axes((0.05, 0.8, 0.90, 0.19))
        #
        yh, xh = np.histogram(np.log(1+zfit.z_max_spec[has_gris]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None')
        yh, xh = np.histogram(np.log(1+zfit.z_max_spec[bin]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', color='red', linewidth=2)
        a = ax.set_yticklabels([])
        #
        ax = fig.add_axes((0.05, 0.05, 0.90, 0.70))
        a = ax.plot(r0-cat.ra, cat.dec, ',', color='0.8', alpha=0.1)
        a = ax.plot(r0-cat.ra[bin], cat.dec[bin], 'o', alpha=0.4, color='red', ms=5)
        a = ax.set_yticklabels([])
        a = ax.set_xticklabels([])
        #
        unicorn.plotting.savefig(fig, 'goodsn_zhist_%04d.png' %(count+1))
    
    ### nearest neighbor density plot
    width = 0.005
    NX, NY = cat.x_image.max(), cat.y_image.max()
    fact = 100
    yi, xi = np.indices((int(NY/fact),int(NX/fact)))
        
    import scipy.spatial

    bin = cat.x_image > 0
    xy = np.array([cat.x_image[bin]/fact, cat.y_image[bin]/fact])
    tree = scipy.spatial.cKDTree(xy.T, 10)
    #
    N = 7
    mask = xi*0.
    for x in xrange(int(NX/fact)):
        #print unicorn.noNewLine+'%d' %(x)
        for y in xrange(int(NY/fact)):
            mask[y,x] = tree.query([x,y], k=N)[0][-1]
    
    mask = mask < 3
    
    ### which redshift to use
    zuse = zfit.z_max_spec
    has_gris = zuse > 0.01
    
    zuse = zout.z_peak
    has_gris = (dq.mag < 25) & (cat.use == 1)
    
    aspect, xs = 1/0.6, 5
    max_scale = np.log(0.2)
            
    if root == 'cosmos':
        aspect = 1./0.4
        xs = 5
        max_scale = np.log(0.1)
        
    zr = (0.65, 2.4)
    
    for count, lz in enumerate(np.arange(np.log(1+zr[0]), np.log(1+zr[1]), width)):
        print unicorn.noNewLine+'%d %f' %(count, lz)
        ztry = np.exp(lz)-1
        #
        bin = has_gris & (np.abs(zuse-ztry)/(1+ztry) < width)
        #
        fig = unicorn.plotting.plot_init(square=True, xs=xs, aspect=aspect, NO_GUI=True)
        ax = fig.add_axes((0.05, 0.8, 0.90, 0.19))
        #
        yh, xh = np.histogram(np.log(1+zuse[has_gris]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None')
        yh, xh = np.histogram(np.log(1+zuse[bin]), range=(0,np.log(zmax+1)), bins=np.log(zmax+1)/0.005)
        a = ax.plot(np.exp(xh[1:])-1, yh, linestyle='steps', marker='None', color='red', linewidth=2)
        a = ax.set_yticklabels([])
        a = ax.text(0.95, 0.8, 'z=%.3f' %(ztry), transform=ax.transAxes, horizontalalignment='right', fontsize=12)
        #
        ax = fig.add_axes((0.05, 0.05, 0.90, 0.70))
        #
        xy = np.array([cat.x_image[bin]/fact, cat.y_image[bin]/fact])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        #
        R = xi*0.
        for x in xrange(int(NX/fact)):
            #print unicorn.noNewLine+'%d' %(x)
            for y in xrange(int(NY/fact)):
                R[y,x] = tree.query([x,y], k=N)[0][-1]
        #
        density = N/np.pi/R**2
        #a = ax.imshow(0-np.log(density*mask), interpolation='Nearest', vmin=max_scale, vmax=6.5)
        #### Uncertainty
        # bin2 = has_gris & (np.abs(zuse-ztry)/(1+ztry) < width)
        # xy = np.array([cat.x_image[bin2]/fact, cat.y_image[bin2]/fact])
        # tree = scipy.spatial.cKDTree(xy.T, 10)
        # Nbin = bin2.sum()
        # points = np.zeros(Nbin)
        # for i in range(Nbin):
        #     points[i] = tree.query([cat.x_image[bin2][i]/fact, cat.y_image[bin2][i]/fact], k=N+1)[0][-1]
        # #
        # sigma = np.std(N/np.pi/points**2) 
        density[~mask] = -10*sigma
        #ai = ax.imshow(0-density/sigma, interpolation='Nearest', vmin=-10, vmax=0)
        #ai = ax.imshow(density/sigma, interpolation='Nearest', vmin=0, vmax=10)
        ai = ax.imshow(np.log(density/sigma), interpolation='Nearest', vmin=-6, vmax=max_scale)
        a = ax.scatter(xy[0,:], xy[1,:], alpha=0.4, color='white', s=5)
        a = ax.set_yticklabels([])
        a = ax.set_xticklabels([])
        a = ax.set_xlim(0, R.shape[1])#
        a = ax.set_ylim(0, R.shape[0])
        # ac = fig.add_axes((0.03, 0.2, 0.07, 0.5))
        # cb = plt.colorbar(ai, cax=ac)
        # cb.set_label(r'$\Sigma / \sigma$')
        #
        unicorn.plotting.savefig(fig, root+'_nn_%5.3f.png' %(ztry))
        
    #### Add zFOURGE
    centers = [[150.133, 2.302], ['10:00:15.769','+02:15:39.52'], ['10:00:18.394','+02:14:58.78'], ['10:00:23.562','+02:14:34.10']]
    xy = np.array([cat.ra, cat.dec])
    tree = scipy.spatial.cKDTree(xy.T, 10)
    for center in centers:
        if isinstance(center[0],str):
            ra = threedhst.utils.DMS2decimal(center[0], hours=True)
            dec = threedhst.utils.DMS2decimal(center[1], hours=False)
        else:
            ra, dec = center[0], center[1]
        #
        idx = tree.query([ra, dec], k=2)[1][0]
        a = ax.plot(cat.x_image[idx]/fact, cat.y_image[idx]/fact, marker='o', color='white', alpha=0.2, ms=20)
    #
    unicorn.plotting.savefig(fig, root+'_zFOURGE_%5.3f.png' %(ztry))
    
    ### 3D plot, doesn't really work
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    zsel = (zfit.z_max_spec > 0.7) & (zfit.z_max_spec < 1.5)
    ax.scatter(cat.ra[has_gris & zsel], cat.dec[has_gris & zsel], zfit.z_max_spec[has_gris & zsel])
        
    rf = catIO.Readfile('/Users/gbrammer/research/drg/PHOTZ/EAZY/GOODS_F140W/HIGHRES/OUTPUT/goodsn1.7.153-155.rf')
    
    rf_uv = catIO.Readfile('/research/HST/GRISM/3DHST/GOODS-S/FIREWORKS/fireworks.153-155.rf')
    rf_vj = catIO.Readfile('/research/HST/GRISM/3DHST/GOODS-S/FIREWORKS/fireworks.155-161.rf')
    
    uv = (-2.5*np.log10(rf_uv.l153/rf_uv.l155))[zfit.phot_id-1]
    vj = (-2.5*np.log10(rf_vj.l155/rf_vj.l161))[zfit.phot_id-1]
    
    plt.plot(uv[~keep], dz[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.01,2.5)
    plt.plot(uv[keep], dz[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.01,2.5) 
    xm, ym, ys, N = threedhst.utils.runmed(uv[keep], dz[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(uv[keep], dzphot[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    
    plt.scatter(vj[keep], uv[keep], marker='o', s=sz[keep], alpha=0.4)
    
    plt.scatter(zfit.z_max_spec[keep], uv[keep], marker='o', s=sz[keep], alpha=0.4)

    plt.scatter(zfit.z_max_spec[keep], zfit.mag[keep], marker='o', s=sz[keep], alpha=0.4)

    plt.scatter(zfit.z_max_spec[keep], cat['[24um]_totf'][zfit.phot_id-1][keep], marker='o', s=sz[keep], alpha=0.4)
    plt.ylim(0.01,1000); plt.semilogy()

    plt.scatter(zfit.z_max_spec[keep], cat.Kr50[keep], marker='o', s=sz[keep], alpha=0.4)
    
    plt.plot(np.abs(dz[~keep]), np.abs(dzphot[~keep]), color='red', alpha=0.2)
    plt.plot(np.abs(dz[keep]), np.abs(dzphot[keep]), color='blue', alpha=0.4)
    plt.plot([1.e-6,10],[1e-6,10], color='black', alpha=0.4, linewidth=2, marker='', linestyle='-')
    plt.plot([1.e-6,10],[1e-5,100], color='black', alpha=0.4, linewidth=2, marker='', linestyle='--')
    plt.loglog()
    plt.xlim(1.e-5,8); plt.ylim(1.e-5,8)
    
    delta = np.abs(zfit.z_max_spec - zfit.z_peak_phot)/(1+zfit.z_peak_phot)
    plt.plot(delta[~keep], dz_spec[~keep], color='red', alpha=0.2); plt.ylim(-0.2,0.2); plt.xlim(0.001,1.5) ; plt.semilogx()
    plt.plot(delta[keep], dz_spec[keep], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0.001,1.5) ; plt.semilogx()
    xm, ym, ys, N = threedhst.utils.runmed(delta[keep], dz_spec[keep], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    
    keep = keep & (delta < 0.1)
    
    yh, xh, qq = plt.hist(dz_spec[keep], range=(-0.01, 0.01), bins=40, alpha=0.5)
    xx = xh[:-1]/2.+xh[1:]/2.
    import gbb.pymc_gauss as gaussfit
    fit = gaussfit.init_model(xx, yh, np.maximum(np.sqrt(yh),1))
    NSAMP,NBURN=11000,1000
    fit.sample(NSAMP, NBURN)
    trace = fit.trace('eval_gaussian')[::NSAMP/25,:]
    plt.errorbar(xx, yh, np.sqrt(yh))
    for i in range(trace.shape[0]):
        plt.plot(xx, trace[i,:], color='red', alpha=0.2, linestyle='-', marker='')
        
    ###
    test = keep
    test = keep & nopartial
    plt.plot(zfit.z_spec[test], dz[test], color='blue', alpha=0.4); plt.ylim(-0.2,0.2); plt.xlim(0,4)
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[test], dz[test], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='green', linestyle='-')
    xm, ym, ys, N = threedhst.utils.runmed(zfit.z_spec[test], dzphot[test], NBIN=20, use_nmad=False, use_median=True)
    plt.plot(xm, ys, color='orange', linestyle='-')
    
    bad = keep & nopartial & (np.abs(dz) > 0.02)
    fp = open('/tmp/bad','w')
    for id in zfit.id[bad]:
        fp.write('/tmp/%s.zfit.png\n' %(id))
    fp.close()
    
    #### Local environmental density
    width = 0.02
    N = 7
    r0, d0 = np.median(cat.ra), np.median(cat.dec)
    dra = (cat.ra-r0)*3600.*np.cos(d0/360*2*np.pi)
    ddec = (cat.dec-d0)*3600.
    
    idx = np.arange(len(cat.ra))
    
    nearest = np.zeros((cat.ra.shape[0], N))
    for i in idx[has_gris]:
        ztry = zfit.z_max_spec[i]
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        xy = np.array([dra[bin], ddec[bin]])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        dist, did = tree.query([dra[i], ddec[i]], k=N+1)
        nearest[i,:] = dist[1:]
    
    zbin = (zfit.z_max_spec > 0.9) & (zfit.z_max_spec < 1.4) & (fout.lmass > 10.5)
    halpha = lines[1].data['HALPHA_FLUX']*1.
    halpha[halpha == 0] += 0.01 * (1+0.05*np.random.normal(size=(halpha == 0).sum()))
    plt.plot(1./nearest[zbin, 0]**2, halpha[zbin], alpha=0.2)
    
    #### Pairs
    min_dist = 1.5
    pairs = (zfit.z_max_spec > 0.8) & (zfit.z_max_spec < 2.5) & (fout.lmass > 10.) & (nearest[:,0] < min_dist)
    detect = pyfits.open('COSMOS_v2.0_PHOTOMETRY/Detection/COSMOS_F125W-F140W-F160W_detect.fits')
    
    ii = idx[pairs]
    ii = ii[np.argsort(zfit.z_max_spec[ii])]
    
    DX = 5/0.06
    for i in ii:
        ztry = zfit.z_max_spec[i]
        bin = has_gris & (np.abs(zfit.z_max_spec-ztry)/(1+ztry) < 2*width)
        xy = np.array([dra[bin], ddec[bin]])
        tree = scipy.spatial.cKDTree(xy.T, 10)
        dist, did = tree.query([dra[i], ddec[i]], k=N+1)
        neighbors = did[(dist > 0) & (dist < min_dist)]
        #
        xc, yc = int(np.round(cat.x_image[i])), int(np.round(cat.y_image[i]))
        subim = detect[0].data[yc-DX:yc+DX, xc-DX:xc+DX]
        plt.imshow(0-subim, vmin=-10, vmax=0.03, interpolation='Nearest')
        plt.text(DX, DX+5, '%5.3f' %(ztry), horizontalalignment='center', verticalalignment='bottom', color='red')
        for neighbor in neighbors:
            plt.plot(cat.x_image[bin][neighbor]-xc+DX, cat.y_image[bin][neighbor]-yc+DX, marker='o', ms=12, color='yellow', alpha=0.3)
            plt.text(cat.x_image[bin][neighbor]-xc+DX, cat.y_image[bin][neighbor]-yc+DX+5, '%5.3f' %(zfit.z_max_spec[bin][neighbor]), horizontalalignment='center', verticalalignment='bottom', color='green')
        #
        plt.xlim(0,2*DX)
        plt.ylim(0,2*DX)


def interlace_goodss0():
    """
    Reduce the GOODS-S pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE_v4.1')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F160W', filter='F160W', REFERENCE = '/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v4.0_drz_join.fits', SEGM = '/3DHST/Photometry/Work/GOODS-S/v4/sextr/checkimages/GOODS-S_F125W_F140W_F160W.seg.fits', Force=True)
         
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Work/GOODS-S/v4/sextr/catalogs/GOODS-S_F160W.cat'

    direct=glob.glob('GOODS-S-*-F140W_asn.fits')
            
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='GOODS-S_F160W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        #
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-3[0-9]-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)

    inter = glob.glob('GOODS-S-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduce_scripts.fix_thumbnails(pointing=pointing)
        time.strftime('%X %x %Z')


    ##### Extract and fit only spec-z objects
    models = glob.glob('GOODS-S-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = True)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)


    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-2')

    skip_completed = True

    models = glob.glob('GOODS-S-3[0-9]_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where(model.cat.mag < 23.)
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    models = glob.glob('GOODS-S-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 23.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            try: 
                gris.new_fit_free_emlines(ztry=None)
            except:
                continue

    ### Get some quality flags on the reduced spectra
    root='GOODS-S'

    files = glob.glob('GOODS-S-*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	  
    

def interlace_goodss1():
    
    """
    Reduce the GOODS-S-1* pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE_v4.1')
         
    CATALOG='/3DHST/Photometry/Work/GOODS-S/v4/sextr/catalogs/GOODS-S_F160W.cat'
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='GOODS-S_F160W'


    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-1[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-2')

    skip_completed = True

    models = glob.glob('GOODS-S-1[0-9]_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where(model.cat.mag < 23.)
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
    
def interlace_goodss2():
    
    """
    Reduce the GOODS-S-2* pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE_v4.1')
         
    CATALOG='/3DHST/Photometry/Work/GOODS-S/v4/sextr/catalogs/GOODS-S_F160W.cat'
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='GOODS-S_F160W'


    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-2[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-2')

    skip_completed = True

    models = glob.glob('GOODS-S-2[0-9]_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
        #ii = np.where(model.cat.mag < 24.)
        for id in model.cat.id:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

def interlace_ers():
    
    """
    Reduce the ERS grism pointing in GOODS-S.
    """

    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    from threedhst.prep_flt_files import process_3dhst_pair as pair

    os.chdir(unicorn.GRISM_HOME+'ERS/INTERLACE_v4.0')

    ALIGN_IMAGE = '/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f140w-astrodrizzle-v4.0_drz_sci.fits'
    FORCE = True
    
    #### Main preparation loop
    direct=glob.glob('WFC3-ERSII-G01-F*_asn.fits')
    grism = glob.glob('WFC3-ERSII-G01-G*_asn.fits')
    loop_list = range(len(direct))
    for i in loop_list:
        pointing='WFC3-ERSII-G02-F140W'
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], adjust_targname=False, ALIGN_IMAGE = ALIGN_IMAGE, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift',sky_images=['G102_master_flatcorr.v2.fits'])
            
    direct=glob.glob('WFC3-ERSII-G02-F*_asn.fits')
    grism = glob.glob('WFC3-ERSII-G02-G*_asn.fits')
    loop_list = range(len(direct))
    for i in loop_list:
        pointing='WFC3-ERSII-G02-F140W'
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], adjust_targname=False, ALIGN_IMAGE = ALIGN_IMAGE, SKIP_GRISM=False, GET_SHIFT=True, SKIP_DIRECT=False, align_geometry='rotate,shift')
    
   
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F160W', filter='F160W', REFERENCE = '/3DHST/Photometry/Release/v4.0/GOODS-S/HST_Images/goodss_3dhst.v4.0.F160W_orig_sci.fits', SEGM = '/3DHST/Photometry/Release/v4.0/GOODS-S/Detection/goodss_3dhst.v4.0.F160W_seg.fits.gz', Force=True)
         
    NGROW=300
    pad=60
    CATALOG='/3DHST/Photometry/Release/v4.0/GOODS-S/Detection/goodss_3dhst.v4.0.F160W_orig.cat'

    extract_limit = 30.
    skip_completed=False
    REF_ROOT='GOODS-S_F160W'

    direct = glob.glob('WFC3-ERSII-G0*-F*_asn.fits')
    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing = direct[i].split('_asn.fits')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing, 
            NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing, view=True, pad=pad, 
            REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True, growx=1, growy=1,
            auto_offsets=False)
        unicorn.reduce.interlace_combine(pointing, pad=pad, NGROW=NGROW,growx=1, growy=1, auto_offsets=False)
    for pointing in ['WFC3-ERSII-G01-G102','WFC3-ERSII-G02-G141']:
        unicorn.reduce.interlace_combine(pointing, pad=pad, NGROW=NGROW,growx=1, growy=1, auto_offsets=False)

    ##### Generate the spectral model
    ##### Extract all spectra 	
    redo = True
    for pointing,direct,grism in zip(['WFC3-ERSII-G01','WFC3-ERSII-G02'],['F098M','F140W'],['G102','G141']):
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=30., direct=direct, grism=grism, 
                grow_factor=1, growx=1, growy=1)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    for pointing,direct,grism in zip(['WFC3-ERSII-G01','WFC3-ERSII-G02'],['F098M','F140W'],['G102','G141']):
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=30., direct=direct, grism=grism, 
            grow_factor=1, growx=1, growy=1)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
    
def interlace_uds0():
    """
    Reduce the UDS-0* pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import time
    import unicorn.reduce_scripts       
    
    os.chdir(unicorn.GRISM_HOME+'UDS/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(Force=True, REF_ROOT='UDS_F160W', filter='F160W', REFERENCE = '/3DHST/Ancillary/UDS/CANDELS/ASTRODRIZZLE/uds-f160w-astrodrizzle-v4.0_drz_join.fits', SEGM = '/3DHST/Photometry/Work/UDS/v4/sextr/checkimages/UDS_F125W_F140W_F160W.seg.fits')
        
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Work/UDS/v4/sextr/catalogs/UDS_F160W.cat'
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='UDS_F160W'

    direct=glob.glob('UDS-*-F140W_asn.fits')

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=False, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('UDS-*-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    inter = glob.glob('UDS-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduce_scripts.fix_thumbnails(pointing=pointing)
        time.strftime('%X %x %Z')

    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-11')

    ##### Extract and fit only spec-z objects
    models = glob.glob('UDS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = False)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = True)


    #Fit the redshift for objects brighter than F140W=24.
    models = glob.glob('UDS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 23.0, model_limit=25.8, new_fit = False, skip_completed = False)


    #Fit the emission lines
    models = glob.glob('UDS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

def interlace_uds1():
    """
    Reduce the UDS-0* pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import time
    import unicorn.reduce_scripts       

    os.chdir(unicorn.GRISM_HOME+'UDS/INTERLACE_v4.0')
        
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Work/UDS/v4/sextr/catalogs/UDS_F160W.cat'
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='UDS_F160W'

    direct=glob.glob('UDS-1*-F140W_asn.fits')
    direct.remove('UDS-18-F140W_asn.fits')

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('UDS-1[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-11')

    ##### Extract and fit only spec-z objects
    models = glob.glob('UDS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = False)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = False)

    models = glob.glob('UDS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('UDS-1[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

def interlace_uds2():
    """
    Reduce the UDS-0* pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np
    import time
    import unicorn.reduce_scripts       

    #os.chdir(unicorn.GRISM_HOME+'UDS/INTERLACE_v2.1')
        
    NGROW=125
    pad=60
    CATALOG='/3DHST/Photometry/Work/UDS/v4/sextr/catalogs/UDS_F160W.cat'
    extract_limit = 35.
    skip_completed=False
    REF_ROOT='UDS_F160W'

    direct=glob.glob('UDS-2[0-9]-F140W_asn.fits')

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('UDS-2[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
    ##### Extract and fit only spec-z objects
    models = glob.glob('UDS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = False, skip_completed = False)
        unicorn.reduce_scripts.extract_spectra_spec_z(pointing=pointing, model_limit=25.8, 
            new_fit = True, skip_completed = False)

    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-11')

    skip_completed = True

    models = glob.glob('UDS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('UDS-2[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)	
            
            
def interlace_hudf():
    """
    Reduce the GOODS-S pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    #os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE')
    os.chdir(unicorn.GRISM_HOME+'HUDF/INTERLACE')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    #unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F160W', REFERENCE = 'phot/GOODS-S_F160W_sci.fits', SEGM = 'phot/GOODS-S_F160W_v1.seg.fits')
    #unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F160W', REFERENCE = 'phot/HUDF_F160W_detection.fits', SEGM = 'phot/HUDF_F160W_seg.fits')
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F140W', filter='F140W', REFERENCE = 'phot/UDF-F140W_sci_drz.fits', SEGM = 'phot/HUDF_F160W_seg.fits')
        
    NGROW=125
    pad=60
    #CATALOG='phot/GOODS-S_F160W_v1.cat'
    CATALOG='phot/convolved_h_nomean.cat'

    direct=glob.glob('GOODS-S-3[4678]-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='GOODS-S_F140W'

#    ##### Generate the interlaced images, including the "blotted" detection image
#    for i in range(len(direct)):
#        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
#        #
#        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
#        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
#        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
#        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model
    ##### Extract all spectra 	
    inter = glob.glob('GOODS-S-3[4678]*-G141_inter.fits')
    redo = False

    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]

        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            #model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
            #model.extract_spectra_and_diagnostics(MAG_LIMIT=24)

    ##### Extract and fit only mag>24 objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='GOODS-S-34')

    skip_completed = True

    models = glob.glob('GOODS-S-36_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        #
        ii = np.where(model.cat.mag < 24.)

        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print id, '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    models = glob.glob('GOODS-S-34_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)

    ### Get some quality flags on the reduced spectra
    root='HUDF'

    files = glob.glob('GOODS-S-3[4678]*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	
            
def interlace_uds18():
    """
    Reduce the COSMOS pointings on Unicorn and extract spectra, using the full
    mosaic as the detection image.
    """
    import threedhst
    import unicorn
    import glob
    import os
    import numpy as np

    os.chdir(unicorn.GRISM_HOME+'UDS/UDS-18')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='UDS-18_F140W', filter='F140W', REFERENCE = '/3DHST/Photometry/Work/UDS/UDS-18/psfmatch/run/UDS-18_F140W_sci_new/UDS_18_F140W_sci_conv.fits', SEGM = '/3DHST/Photometry/Work/UDS/UDS-18/sextr/checkimages/UDS-18_F140W.seg.fits')

        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = '/3DHST/Photometry/Work/UDS/UDS-18/sextr/catalogs/UDS-18_F140W.cat'
    direct=glob.glob('UDS-18-F140W_asn.fits')

    extract_limit = 25.8
    skip_completed=False
    REF_ROOT='UDS-18_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('UDS-18-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-18')

    skip_completed = True

    models = glob.glob('UDS-18_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.zfit.png') & skip_completed:
                continue  
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)

    #Fit the emission lines
    models = glob.glob('UDS-18_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
        ii = np.where((model.cat.mag < 25.8))
        for id in model.cat.id[ii]:
            root='%s_%05d' %(pointing, id)
            if not os.path.exists(root+'.2D.fits'):
                status = model.twod_spectrum(id)
                if not status:
                    continue
            #
            if os.path.exists(root+'.linefit.png') & skip_completed:
                continue  
            if not os.path.exists(root+'.zfit.png'):
                continue
            #
            try:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
            except:
                continue
            #
            if gris.status is False:
                continue
            #
            if gris.dr > 1:
                continue
            #
            print '\n'
            gris.new_fit_free_emlines(ztry=None)	

    ### Get some quality flags on the reduced spectra
    root='UDS-18'

    files = glob.glob('UDS-18*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.zfit.dat','w')
    first = True
    for file in files:
        print unicorn.noNewLine+file
        spec = unicorn.interlace_fit.GrismSpectrumFit(file.split('.zfit')[0], verbose=False)   
        lines = open(file).readlines()
        status = spec.stats(return_string=True)
        if status is not False:
            if first:
                fp.write(status[0])
                fp2.write(lines[0])
                first = False
            #    
            fp.write(status[1])
            fp2.write(lines[-1])

    fp.close()
    fp2.close()	  	

def interlace_sntile41():
    
    #
    ### SWarp for color image
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F160W.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_wht.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F160W_wht.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/Incoming/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F125W.fits')
    threedhst.shifts.matchImagePixels( input=glob.glob('/3DHST/Ancillary/COSMOS/ACS/*sci.fits'), matchImage='TILE41-F160W_drz.fits', match_extension=1, output='match_F814W.fits')

    ### F160W
    sub_r = pyfits.open('match_F160W.fits')[0].data*10**(-0.4*(25.96-25.96))
    ### F125W
    sub_g = pyfits.open('match_F125W.fits')[0].data*10**(-0.4*(26.25-25.96))
    ### F814W
    sub_b = pyfits.open('match_F814W.fits')[0].data*10**(-0.4*(25.94-25.96))*1.5

    Q, alpha, m0 = 4, 12, -0.005
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='TILE41.png', shape=sub_r.shape)
    Q, alpha, m0 = 3.5, 5, -0.01
    unicorn.candels.luptonRGB(sub_r, sub_g, sub_b, Q=Q, alpha=alpha, m0=m0, filename='TILE41_2.png', shape=sub_r.shape)
    
    unicorn.reduce.interlace_combine('TILE41-132-G141', NGROW=125)
    unicorn.reduce.interlace_combine('TILE41-132-F160W', NGROW=125, auto_offsets=True)
    
    model = unicorn.reduce.process_GrismModel('TILE41-132', MAG_LIMIT=25)
        
    se = threedhst.sex.SExtractor()
    
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    
    ## XXX add test for user-defined .conv file
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = 'test.cat'
    se.options['CHECKIMAGE_NAME'] = 'test_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'match_F160W_wht.fits[0]'
    se.options['FILTER']    = 'Y'
    threedhst.sex.USE_CONVFILE =  'gauss_2.0_5x5.conv'
    threedhst.sex.USE_CONVFILE =  'gauss_4.0_7x7.conv'
    se.copyConvFile()
    se.options['FILTER_NAME'] = threedhst.sex.USE_CONVFILE
    
    #### Detect thresholds (default = 1.5)
    se.options['DEBLEND_NTHRESH']    = '16' 
    se.options['DEBLEND_MINCONT']    = '0.001' 
    se.options['DETECT_MINAREA']    = '36' 
    se.options['DETECT_THRESH']    = '1.5' 
    se.options['ANALYSIS_THRESH']  = '1.5'
    se.options['MAG_ZEROPOINT'] = '%.2f' %(model.ZP)
    
    #### Run SExtractor
    status = se.sextractImage('match_F160W.fits[0]')
    threedhst.sex.sexcatRegions('test.cat', 'test.reg', format=2)
    
    unicorn.reduce.prepare_blot_reference(REF_ROOT='TILE41_F160W', filter='F160W', REFERENCE = 'match_F160W.fits', SEGM = 'test_seg.fits', Force=False, sci_extension=0)
    
    unicorn.reduce.blot_from_reference(REF_ROOT='TILE41_F160W', DRZ_ROOT = 'TILE41-132-F160W', NGROW=125, verbose=True)
    unicorn.reduce.interlace_combine_blot(root='TILE41-132-F160W', view=True, pad=60, REF_ROOT = 'TILE41_F160W', CATALOG='test.cat',  NGROW=125, verbose=True, growx=2, growy=2, auto_offsets=True, NSEGPIX=0)
    #unicorn.reduce.fill_inter_zero('TILE41-132_inter_seg.fits')
    unicorn.reduce.fill_inter_zero('TILE41-132_ref_inter.fits')
    
    unicorn.reduce.interlace_combine('TILE41-132-G141', NGROW=125, auto_offsets=True)
    unicorn.reduce.interlace_combine('TILE41-132-F160W', NGROW=125, auto_offsets=True)
    unicorn.reduce.fill_inter_zero('TILE41-132-F160W_inter.fits')
    
    model = unicorn.reduce.process_GrismModel('TILE41-132', MAG_LIMIT=25, REFINE_MAG_LIMIT=23, make_zeroth_model=False)
    
    unicorn.reduce.blot_from_reference(REF_ROOT='TILE41_F160W', DRZ_ROOT = 'TILE41-152-F140W', NGROW=125, verbose=True)
    unicorn.reduce.interlace_combine_blot(root='TILE41-152-F140W', view=True, pad=60, REF_ROOT = 'TILE41_F160W', CATALOG='test.cat',  NGROW=125, verbose=True, growx=2, growy=2, auto_offsets=True, NSEGPIX=0)
    
    unicorn.reduce.fill_inter_zero('TILE41-152_ref_inter.fits')
    
    unicorn.reduce.interlace_combine('TILE41-152-G141', NGROW=125, auto_offsets=True)
    unicorn.reduce.interlace_combine('TILE41-152-F140W', NGROW=125, auto_offsets=True)
    unicorn.reduce.fill_inter_zero('TILE41-152-F140W_inter.fits')
    
    model = unicorn.reduce.process_GrismModel('TILE41-152', MAG_LIMIT=25, REFINE_MAG_LIMIT=23)
    
    #### Seems OK, but haven't been able to check against z_specs
    
    ### Extract all objects > z>0.6, m140 > X
    cat, zout, fout = unicorn.analysis.read_catalogs('TILE41')
    matcher = catIO.CoordinateMatcher(cat)
    
    dr, idx = matcher.match_list(ra=model.ra_wcs, dec=model.dec_wcs)
    
    ok = (dr < 0.8) & ((zout.z_peak[idx] > 0.65) & (model.cat.mag < 23)) | ((zout.z_peak[idx] > 1.2) & (model.cat.mag < 24))
    
    for id in model.objects[ok]:
        try:
            model.twod_spectrum(id)
            gris = unicorn.interlace_fit.GrismSpectrumFit('TILE41-152_%05d' %(id))
            gris.fit_in_steps(zrfirst=(0.5, 4))
            #os.system('open TILE41-132_%05d*zfit*png' %(id))
        except:
            pass

def mast_header(filename='', field=None, band=None, instrument='wfc3', version='v4.0', im_type = None):

    import pyfits
    
    if field is None:
        field = filename.split('-')[0].lower()
        if field.startswith('goods'):
            field = field+filename.split('-')[1].lower()
            if band is None:
                band = filename.split('-')[2].lower()
        else:
            if band is None:
                band = filename.split('-')[1].lower()
    if im_type is None:
        im_type = filename.split('_')[-1].lower().split('.')[0]
    print 'FIELD: %s     BAND: %s    IM_TYPE: %s'%(field, band, im_type)
    
    go_ids = []
    if ((field == 'aegis') and ((band == 'f125w') or (band == 'f160w'))):
        go_ids = ['12063']
    if ((field == 'aegis') and (band == 'f140w')):
        go_ids = ['12177']
    if ((field == 'cosmos') & ((band == 'f125w') | (band == 'f160w'))):
        go_ids = ['12440','12461']
    if ((field == 'cosmos') and (band == 'f140w')):
        go_ids = ['12328']
    if ((field == 'goodsn') and ((band == 'f125w') or (band == 'f160w'))):
        go_ids = ['12443','12444','12445','12461']
    if ((field == 'goodsn') and (band == 'f140w')):
        go_ids = ['11600','12461']
    if ((field == 'goodss') and ((band == 'f125w') or (band == 'f160w'))):
        go_ids = ['11359','11563','12061','12062','12099']
    if ((field == 'goodss') and (band == 'f140w')):
        go_ids = ['11359','12177','12190']
    if ((field == 'uds') and ((band == 'f125w') or (band == 'f160w'))):
        go_ids = ['12064','12099']
    if ((field == 'uds') and (band == 'f140')):
        go_ids = ['12328']

    print go_ids

    ra_targ = 0.0
    dec_targ = 0.0
    if field == 'aegis':
        ra_targ = '14:18:36.00'
        dec_targ = '+52:39:0.00'
    if field == 'cosmos':
        ra_targ = '10:00:31.00'
        dec_targ = '+02:24:0.00'
    if field == 'goodsn':
        ra_targ = '12:35:54.98'
        dec_targ = '+62:11:51.30'
    if field == 'goodss':
        ra_targ =  '03:32:30.00'
        dec_targ = '-27:47:19.00'
    if field == 'uds':
        ra_targ = '02:17:49.00'
        dec_targ = '-05:12:02.00'
            
    outname = 'hlsp_3dhst_hst_'+instrument.lower()+'_'+field.lower()+'_'+band.lower()+'_'+version.lower()+'_'+im_type+'.fits'
    rootname = 'hlsp_3dhst_hst_'+instrument.lower()+'_'+field.lower()+'_'+band.lower()+'_'+version.lower()
    print 'WRITING OUT: %s' % (outname)
    

    file = pyfits.open(filename)
    
    old_header = file[0].header.copy()
        
    del file[0].header['HISTORY']
    
    start_key = file[0].header.keys().index('NEXTEND')
    for key in file[0].header.keys()[start_key:-1]:
        del file[0].header[key]

    
    file[0].header.update('CTYPE1', value = old_header['CTYPE1'], comment='the coordinate type for the first axis')
    file[0].header.update('CTYPE2',  value = old_header['CTYPE2'], comment=' the coordinate type for the second axis')
    file[0].header.update('CRPIX1',  value = old_header['CRPIX1'], comment=' x-coordinate of reference pixel')
    file[0].header.update('CRPIX2',  value = old_header['CRPIX2'], comment=' y-coordinate of reference pixel')
    file[0].header.update('CRVAL1',  value = old_header['CRVAL1'], comment=' first axis value at reference pixel')
    file[0].header.update('CRVAL2',  value = old_header['CRVAL2'], comment=' second axis value at reference pixel')
    file[0].header.update('ORIENTAT',  value = old_header['ORIENTAT'], comment=' position angle of image y axis (deg. e of n)')
    file[0].header.update('CDELT1',  value = old_header['CDELT1'])
    file[0].header.update('CDELT2',  value = old_header['CDELT2'])
    file[0].header.update('CD1_1',  value = old_header['CD1_1'], comment=' partial of first axis coordinate w.r.t. x')
    file[0].header.update('CD1_2',  value = old_header['CD1_2'], comment=' partial of first axis coordinate w.r.t. y')
    file[0].header.update('CD2_1',  value = old_header['CD2_1'], comment=' partial of second axis coordinate w.r.t. x')
    file[0].header.update('CD2_2',  value = old_header['CD2_2'], comment=' partial of second axis coordinate w.r.t. y')
    #file[0].header.update('RADESYS',  value = old_header['RADESYS'])
    file[0].header.update('NEXTEND',  value = old_header['NEXTEND'], comment=' Number of standard extensions')
    file[0].header.update('FILENAME',  value = outname, comment=' name of file')
    file[0].header.update('FILETYPE',  value = im_type, comment=' type of data found in data file')
    file[0].header.update('TELESCOP',  value = old_header['TELESCOP'], comment=' telescope used to acquire data')
    file[0].header.update('INSTRUME',  value = old_header['INSTRUME'], comment=' identifier for instrument used to acquire data')
    file[0].header.update('EQUINOX',  value = old_header['EQUINOX'], comment=' equinox of celestial coord. system')


    file[0].header.update('ROOTNAME',  value = rootname, comment=' rootname of the observation set')
    file[0].header.update('IMAGETYP',  value = old_header['IMAGETYP'],comment=' type of exposure identifier')
    file[0].header.update('PRIMESI' ,  value = old_header['PRIMESI'],comment=' instrument designated as prime')


    file[0].header.update('PROPOSID', value = '12177,12328', comment=' PEP proposal identifier')
    file[0].header.update('PR_INV_L', 'van Dokkum', comment=' last name of principal investigator')
    file[0].header.update('PR_INV_F', 'Pieter', comment=' first name of principal investigator')
    file[0].header.update('PR_INV_M', '', comment=' middle name initial of principal investigator')


    file[0].header.update('OBSTYPE',  value = old_header['OBSTYPE'], comment=' observation type - imaging or spectroscopic')
    file[0].header.update('OBSMODE',  value = old_header['OBSMODE'], comment=' operating mode')
    file[0].header.update('SCLAMP',  value = old_header['SCLAMP'], comment=' lamp status, NONE or name of lamp which is on')
    file[0].header.update('NRPTEXP',  value = old_header['NRPTEXP'], comment=' number of repeat exposures in set: default 1')
    file[0].header.update('SUBARRAY',  value = old_header['SUBARRAY'], comment=' data from a subarray (T) or full frame (F)')
    file[0].header.update('SUBTYPE',  value = old_header['SUBTYPE'], comment=' Size/type of IR subarray')
    file[0].header.update('DETECTOR',  value = old_header['DETECTOR'], comment=' detector in use: UVIS or IR')
    file[0].header.update('FILTER',  value = old_header['FILTER'], comment=' element selected from filter wheel')
    file[0].header.update('APERTURE',  value = old_header['APERTURE'], comment=' aperture name')
    file[0].header.update('PROPAPER',  value = old_header['PROPAPER'], comment=' proposed aperture name')
    file[0].header.update('DIRIMAGE',  value = old_header['DIRIMAGE'], comment=' direct image for grism or prism exposure')


    #file[0].header.update('DATE-OBS',  value = old_header['DATE-OBS'], comment='start of first observation')
    file[0].header.update('EXPTIME',  value = old_header['EXPTIME'], comment=' composite image')
    file[0].header.update('EXPSTART',  value = old_header['EXPSTART'], comment=' start of first observation')
    file[0].header.update('EXPEND',  value = old_header['EXPEND'],comment= ' end time of last observation')


    file[0].header.update('PHOTMODE',  value = old_header['PHOTMODE'], comment=' observation con')
    file[0].header.update('PHOTFLAM',  value = old_header['PHOTFLAM'], comment=' inverse sensitivity, ergs/cm2/Ang/electron')
    file[0].header.update('PHOTFNU',  value = old_header['PHOTFNU'], comment=' inverse sensitivity, Jy*sec/electron')
    file[0].header.update('PHOTZPT',  value = old_header['PHOTZPT'], comment=' ST magnitude zero point')
    file[0].header.update('PHOTPLAM',  value = old_header['PHOTPLAM'], comment=' Pivot wavelength (Angstroms)')
    file[0].header.update('PHOTBW',  value = old_header['PHOTBW'], comment=' RMS bandwidth of filter plus detector')
    file[0].header.update('BUNIT',  value = old_header['BUNIT'], comment=' flux unit')
    file[0].header.update('RESOLUTI', value = 0.06, comment= 'Spatial resolution')
    file[0].header.update('RESUNIT',  value = 'arcsec', comment= 'Spatial resolution units')

    file[0].header.update('TARGNAME',  field, comment = 'field name')
    file[0].header.update('RA_TARG',  ra_targ, comment = 'center of field RA [deg] (J2000)')
    file[0].header.update('DEC_TARG',  dec_targ, comment= 'center field DEC [deg] (J2000)')
    #file[0].header.update('PROPOSALID')
    file[0].header.update('HLSPLEAD',  'Ivelina Momcheva')
    file[0].header.update('PR_INV_L', 'van Dokkum')
    file[0].header.update('PR_INV_F', 'Pieter')

    file[0].header.update('HLSPNAME',  '3D-HST')
    
    file[0].header.add_blank('/ DATA DESCRIPTION KEYWORDS', after='EQUINOX', before='ROOTNAME')
    file[0].header.add_blank('', after='EQUINOX')
    file[0].header.add_blank('',before='ROOTNAME')
    file[0].header.add_blank('/ PROPOSAL INFORMATION', after='PRIMESI', before='PROPOSID')
    file[0].header.add_blank('', after='PRIMESI')
    file[0].header.add_blank('', before='PROPOSID')
    file[0].header.add_blank('/ INSTRUMENT CONFIGURATION INFORMATION', after='PR_INV_M', before='OBSTYPE')
    file[0].header.add_blank('', after='PR_INV_M')
    file[0].header.add_blank('', before='OBSTYPE')
    file[0].header.add_blank('/ DATE AND TIME KEYWORDS', after='DIRIMAGE', before='EXPTIME')         
    file[0].header.add_blank('', after='DIRIMAGE')
    file[0].header.add_blank('', before='EXPTIME')
    file[0].header.add_blank('/ DATA KEYWORDS', after='EXPEND', before='PHOTMODE')
    file[0].header.add_blank('', after='EXPEND')
    file[0].header.add_blank('', before='PHOTMODE')
    file[0].header.add_blank('/IMAGE ORIGIN KEYWORDS', after='RESUNIT', before='TARGNAME')
    file[0].header.add_blank('', after='RESUNIT')
    file[0].header.add_blank('', before='TARGNAME')
    
    
    file[0].header.update('HISTORY','')
    file[0].header.add_history('------------------------------------------------------------')
    file[0].header.add_history(' AstroDrizzle processing performed using:')
    file[0].header.add_history('     AstroDrizzle Version 1.1.9.dev23803')
    file[0].header.add_history('     Numpy Version 1.6.2')
    file[0].header.add_history('     PyFITS Version 3.2.dev2070')
    file[0].header.add_history('------------------------------------------------------------------')
    file[0].header.add_history('')
    file[0].header.add_history('3D-HST:')
    file[0].header.add_history('A Spectroscopic Galaxy Evolution Survey with the Hubble Space Telescope')
    file[0].header.add_history('==================================================================')
    file[0].header.add_history('')
    file[0].header.add_history('3D-HST (Brammer et al. 2012) is a 248-orbit Multi-Cycle Treasury')
    file[0].header.add_history('program with the Hubble Space telescope designed to study the')
    file[0].header.add_history('physical processes that shape galaxies in the distant Universe.')
    file[0].header.add_history('')
    file[0].header.add_history('Survey Principal Investigator: Pieter van Dokkum ')
    file[0].header.add_history('')
    file[0].header.add_history('Website: http://3dhst.research.yale.edu/')
    file[0].header.add_history('')
    file[0].header.add_history('Survey Data Papers')
    file[0].header.add_history('')
    file[0].header.add_history('Brammer, G. B. et al. 2012, ApJS 200, 13B')
    file[0].header.add_history('')
    file[0].header.add_history('Skelton, R. et al., 2014')
    file[0].header.add_history('')
    file[0].header.add_history('This HST image mosaic was produced by Ivelina G. Momcheva, using the ')
    file[0].header.add_history('Tweakreg and Astrodrizzle software (Gonzaga et al., 2012). The mosaic' )
    file[0].header.add_history('was made for the uses of the 3D-HST project and contains data from the' )
    file[0].header.add_history('following HST GO programs:')
    file[0].header.add_history('')
    
    for id in go_ids:
        print '  GO-%s'%(id)
        file[0].header.add_history('  GO-%s'%(id))
        
    file[0].header.add_history('')
    file[0].header.add_history('The mosaic has a scale of 0.06 arcsec/pixel, with north')
    file[0].header.add_history('toward the top in COSMOS, GOODS-N, GOODS-S and UDS. The mosaics ')
    file[0].header.add_history('in the AEGIS field are rotated 49.7 degrees.')
    file[0].header.add_history('  ')
    file[0].header.add_history('We encourage the entire community to make full use of these mosaics.')
    file[0].header.add_history('')
    
    outname = 'hlsp_3dhst_hst_'+instrument.lower()+'_'+field.lower()+'_'+band.lower()+'_'+version.lower()+'_'+im_type+'.fits'
    print 'WRITING OUT: %s' % (outname)
    file.writeto(outname, clobber=True)
    file.close()

def combined_images_script():
    
    os.chdir('/Volumes/3D-HST/MOSAICS/GOODS-N/JOIN')
    unicorn.reduce_scripts.combined_image(main_image = '../ALL_FLT/goods-n-f160w-astrodrizzle-v4.0_drz_sci.fits',fill_image = '../ALL_FLT/goods-n-f140w-astrodrizzle-v4.0_drz_sci.fits')
    
    os.chdir('/Volumes/3D-HST/MOSAICS/COSMOS/JOIN')
    unicorn.reduce_scripts.combined_image(main_image = '../ALL_FLT/cosmos-f160w-astrodrizzle-v4.0_drz_sci.fits',fill_image = '../ALL_FLT/cosmos-f140w-astrodrizzle-v4.0_drz_sci.fits')
    
    os.chdir('/Volumes/3D-HST/MOSAICS/AEGIS/JOIN')
    unicorn.reduce_scripts.combined_image(main_image = '../ALL_FLT/aegis-f160w-astrodrizzle-v4.0_drz_sci.fits',fill_image = '../ALL_FLT/aegis-f140w-astrodrizzle-v4.0_drz_sci.fits')

    os.chdir('/Volumes/3D-HST/MOSAICS/GOODS-S/JOIN')
    unicorn.reduce_scripts.combined_image(main_image = '../ALL_FLT/goods-s-f160w-astrodrizzle-v4.0_drz_sci.fits',fill_image = '../ALL_FLT/goods-s-f140w-astrodrizzle-v4.0_drz_sci.fits')

    os.chdir('/Volumes/3D-HST/MOSAICS/UDS/JOIN')
    unicorn.reduce_scripts.combined_image(main_image = '../ALL_FLT/uds-f160w-astrodrizzle-v4.0_drz_sci.fits',fill_image = '../ALL_FLT/uds-f140w-astrodrizzle-v4.0_drz_sci.fits')
   


def extract_spectra_mag_limit(pointing='UDS-10', mag_limit = 25.8, model_limit=25.8, 
    new_fit = False, skip_completed = False):
    """
    Make a model, if one does not already exist and extract all spectra down to a given
    magnitude.
    """
    
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root=pointing)
    
    model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=model_limit)
    ii = np.where((model.cat.mag < mag_limit))
    for id in model.cat.id[ii]:
        root='%s_%05d' %(pointing, id)
        if not os.path.exists(root+'.2D.fits'):
            status = model.twod_spectrum(id)
            if not status:
                continue
        #
        if ((not new_fit) & os.path.exists(root+'.zfit.png') & skip_completed) | (new_fit & os.path.exists(root+'.new_zfit.png') & skip_completed):
            continue 
        #
        try:
            if new_fit:
                gris = test.SimultaneousFit(root)
            else:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
        except:
            continue
        #
        if gris.status is False:
            continue
        #
        if gris.dr > 1:
            continue
        #
        print '\n'
        if new_fit:
            gris.read_master_templates()
            gris.new_fit_constrained()
        else:
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)    

def extract_spectra_spec_z(pointing='UDS-10', model_limit=25.8, 
    new_fit = False, skip_completed = False):
    """
    Make a model if one doesn not already exist and extract all spectra for objects 
    with ground-based spec-z's.
    """
    
    import unicorn.interlace_test as test
    
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root=pointing)
    
    id_spec_z = cat.id[np.where((cat.z_spec > 0.) & (cat.z_spec < 9.))].astype(int)
    print 'There are {0} objects in with spec_z in this field.'.format(len(id_spec_z))
    
    model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
    print "There are {0} objects with spec_z in {1}.".format(len([id for id in id_spec_z if id in model.cat.id]),pointing)
    for id in [id for id in id_spec_z if id in model.cat.id]:
        root='%s_%05d' %(pointing, id)
        if not os.path.exists(root+'.2D.fits'):
            status = model.twod_spectrum(id)
            if not status:
                continue
        #
        if ((not new_fit) & os.path.exists(root+'.zfit.png') & skip_completed) | (new_fit & os.path.exists(root+'.new_zfit.png') & skip_completed):
            continue 
        #
        try:
            if new_fit:
                gris = test.SimultaneousFit(root)
            else:
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root)
        except:
            continue
        #
        if gris.status is False:
            continue
        #
        if gris.dr > 1:
            continue
        #
        print '\n'
        if new_fit:
            gris.read_master_templates()
            gris.new_fit_constrained()
        else:
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002)
            
            
def combined_image(field='AEGIS', IMAGE_DIR=''):
    """
    Combines the F125W, F140W and F160W images, weighting them by the inverse variance and 
    scaling them to the F140W zeropoint.
    """
    
    import pyfits
    import numpy as np
    
    ZPs = {'F125W':26.25, 'F140W':26.46, 'F160W':25.96}
    
    if field == 'GOODS-N': 
        root='goodsn'
    elif field == 'GOODS-S': 
        root='goodss'
    else:
        root = field.lower()
        
    if IMAGE_DIR == '':
        IMAGE_DIR = '/3DHST/Photometry/Release/v4.0/%s/HST_Images/'%(field.upper())

    im_f160w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F160W_orig_sci.fits.gz'%(root.lower()),memmap=True)
    im_f140w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F140W_orig_sci.fits.gz'%(root.lower()),memmap=True)
    im_f125w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F125W_orig_sci.fits.gz'%(root.lower()),memmap=True)
    wht_f160w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F160W_orig_wht.fits.gz'%(root.lower()),memmap=True)
    wht_f140w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F140W_orig_wht.fits.gz'%(root.lower()),memmap=True)
    wht_f125w = pyfits.open(IMAGE_DIR+'%s_3dhst.v4.0.F125W_orig_wht.fits.gz'%(root.lower()),memmap=True)

    # weight, scale to F140W zeropoint, coadd
    f125w_zp_factor = 10**((ZPs['F140W']-ZPs['F125W'])/2.5)
    f160w_zp_factor = 10**((ZPs['F140W']-ZPs['F160W'])/2.5)
    index = (wht_f160w[0].data > 0) | (wht_f125w[0].data > 0) | (wht_f140w[0].data > 0)
    NEW_WHT = wht_f160w[0].data/np.power(f160w_zp_factor,2)+wht_f125w[0].data/np.power(f125w_zp_factor,2)+wht_f140w[0].data

    # multiply by f125w_zp_factor to scale image and divide by f125w_zp_factor**2 to scale the weight map
    NUM1 = im_f160w[0].data*wht_f160w[0].data/f160w_zp_factor
    NUM2 = im_f125w[0].data*wht_f125w[0].data/f125w_zp_factor
    NUM3 = im_f140w[0].data*wht_f140w[0].data
    
    NEW_IM = np.empty(np.shape(im_f160w[0].data))
    NEW_IM[index] = (NUM1[index]+NUM2[index]+NUM3[index])/NEW_WHT[index]
      
    im_f160w.close()
    im_f125w.close()
    wht_f160w.close()
    wht_f140w.close()
    wht_f125w.close()
    hdu = pyfits.PrimaryHDU(header=im_f140w[0].header,data=NEW_IM)
            
    hdu.writeto('%s_3dhst.v4.0.JOIN.fits'%(root), clobber=True)   
    im_f140w.close() 
    
    
def fix_thumbnails(pointing='AEGIS-1'):
    """
    Objects along the edge of the reference image did not have thumbnails extracted. 
    Rerun these objects with USE_REFERENCE_THUMB=True to get proper thumbnails for 
    them.
    """

    import pyfits
    import numpy as np
    import os
    
inter = glob.glob('COSMOS-*-G141_inter.fits')
for i in range(len(inter)):
    pointing = inter[i].split('-G141_inter')[0]
    model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
    ii=0
    for id in model.cat.id:
        if os.path.exists('{0}_{1:5d}.2D.fits'.format(pointing,id)):
            file_2d = pyfits.open('{0}_{1:5d}.2D.fits'.format(pointing,id))
            if np.max(file_2d[1].data) == 0 and np.max(file_2d[4].data) > 0:
                print '{0} {1}_{2:5d}.2D.fits'.format(ii,pointing,id)
                ii =ii+ 1
                model.twod_spectrum(id = id,USE_REFERENCE_THUMB=True)
                model.show_2d(savePNG=True)
    
    ii=0
    for i in range(len(model.cat.id)):
        id = model.cat.id[i]
        if os.path.exists('{0}_{1:5d}.2D.fits'.format(pointing,id)):
            file_2d = pyfits.open('{0}_{1:5d}.2D.fits'.format(pointing,id))
            if np.max(file_2d[1].data) > 0 and np.max(file_2d[4].data) == 0:
                print '{} {} {}_{:5d}.2D.fits'.format(ii,model.cat.MAG_AUTO[i],pointing,id)
                ii =ii + 1
                if os.path.exists('{}_{:5d}.1D.fits'.format(pointing, id)):
                    print 'Deleting {}_{:5d}.1D.fits'.format(pointing, id)
                    os.system('rm  {}_{:5d}.1D.fits'.format(pointing, id))
                if os.path.exists('{}_{:5d}.1D.png'.format(pointing, id)):
                    print 'Deleting {}_{:5d}.1D.png'.format(pointing, id)
                    os.system('rm  {}_{:5d}.1D.png'.format(pointing, id))



def delete_spectra(field = 'COSMOS'):
    

    inter = glob.glob('{}-*-G141_inter.fits'.format(field))
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        print pointing
        time.strftime('%X %x %Z')
        os.system('rm {}_*.2D.fits'.format(pointing))
        os.system('rm {}_*.1D.fits'.format(pointing))
        os.system('rm {}_*.2D.png'.format(pointing))
        os.system('rm {}_*.1D.png'.format(pointing))
        os.system('rm {}_*.zfit.fits'.format(pointing))
        os.system('rm {}_*.zfit.png'.format(pointing))
        os.system('rm {}_*.zfit.pz.fits'.format(pointing))
        os.system('rm {}_*.zfit.dat'.format(pointing))
        os.system('rm {}_*.xxx'.format(pointing))
    
    
    
    
    
    
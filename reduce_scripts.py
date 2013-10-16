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
    CATALOG = 'phot/F140W.cat'
    direct=glob.glob('AEGIS-[0-9]-F140W_asn.fits')
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
    inter = glob.glob('AEGIS-[0-9]-G141_inter.fits')
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
    
    inter = glob.glob('AEGIS-*-G141_inter.fits')
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        unicorn.reduc_scripts.fix_thumbnails(pointing=pointing)
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

    models = glob.glob('AEGIS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 24.0, model_limit=25.8, new_fit = False, skip_completed = True)

    #Fit the emission lines
    models = glob.glob('AEGIS-[0-9]_inter_model.fits')
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
            mag_limit = 24.0, model_limit=25.8, new_fit = False, skip_completed = True)

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

                        
    ### Fit all grism redshifts down to F160/F140=24.0:

    models = glob.glob('AEGIS-[2-3][0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 24.0, model_limit=25.8, new_fit = False, skip_completed = True)

    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

    #Fit the emission lines
    models = glob.glob('AEGIS-2[0-9]_inter_model.fits')
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
    inter = glob.glob('COSMOS-[0-9]-G141_inter.fits')
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

    models = glob.glob('COSMOS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        unicorn.reduce_scripts.extract_spectra_mag_limit(pointing=pointing, 
            mag_limit = 24.0, model_limit=25.8, new_fit = False, skip_completed = True)

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
    models = glob.glob('COSMOS-2[0-9]_inter_model.fits')
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
    inter = glob.glob('GOODS-N-4[0-9]-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if not (os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)

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
    models = glob.glob('GOODS-N-3[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        ii = np.where((model.cat.mag < 24.))# & (model.cat.mag > 17.)) #& (model.mag.id != 1876))
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
                gris = unicorn.interlace_fit.GrismSpectrumFit(root=root,p_flat = 0.0)
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
    models = glob.glob('GOODS-N-4[0-9]_inter_model.fits')
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

    ### Get some quality flags on the reduced spectra
    root='GOODS-N'
    
    files = glob.glob('GOODS-N-4*zfit.dat')
    fp = open(root+'.dqflag.dat','w')
    fp2 = open(root+'.new_zfit.dat','w')
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

    os.chdir(unicorn.GRISM_HOME+'GOODS-S/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='GOODS-S_F160W', filter='F160W', REFERENCE = '/3DHST/Ancillary/GOODS-S/CANDELS/ASTRODRIZZLE/goods-s-f160w-astrodrizzle-v4.0_drz_join.fits', SEGM = '/3DHST/Photometry/Work/GOODS-S/v4/sextr/checkimages/GOODS-S_F125W_F140W_F160W.seg.fits')
         
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
    inter = glob.glob('GOODS-S-*-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)


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

    models = glob.glob('GOODS-S-[0-9]_inter_model.fits')
    for file in models:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
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

    models = glob.glob('GOODS-S-[0-9]_inter_model.fits')
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
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
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
    inter = glob.glob('UDS-[0-9]-G141_inter.fits')
    redo = True
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=35.)
            
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


    #### Extract and fit objects with spec-z
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
            
            
def combined_image(main_image='F160W.fits', fill_image='F140W.fits', main_zp=25.9463, fill_zp=26.4524):
    """
    Compares the weight maps of the main image and the fill image and adds the 
    fill image to the main image where the two do not overlap. The fill image is 
    scaled to the same zeropoint
    """
    
    import pyfits
    
    im_main = pyfits.open(main_image, mmap=True)
    wht_main = pyfits.open(main_image.replace('sci.fits','wht.fits'), mmap=True)
    h_main = im_main[0].header
    data_main = im_main[0].data
    
    im_fill = pyfits.open(fill_image,   mmap = True)
    wht_fill = pyfits.open(fill_image.replace('sci.fits','wht.fits'), mmap=True)
    
    data_main[(wht_main[0].data == 0) & (wht_fill[0].data>0)] = im_fill[0].data[(wht_main[0].data == 0) & (wht_fill[0].data>0)]*(10**((main_zp-fill_zp)/2.5))
    data_main[im_main[0].data > 1e5] = 0
    
    hdu = pyfits.PrimaryHDU(header=h_main,data=data_main)
            
    hdu.writeto(main_image.split('/')[-1].replace('sci.fits','join.fits'), clobber=True)
    
def fix_thumbnails(pointing='AEGIS-1'):
    """
    Objects along the edge of the reference image did not have thumbnails extracted. 
    Rerun these objects with USE_REFERENCE_THUMB=True to get proper thumbnails for 
    them.
    """

    import pyfits
    import numpy as np
    
    model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=35.)
    ii=0
    for id in model.cat.id:
        if os.path.exists('{0}_{1:5d}.2D.fits'.format(pointing,id)):
            file_2d = pyfits.open('{0}_{1:5d}.2D.fits'.format(pointing,id))
            if np.min(file_2d[1].data) == 0:
                print '{0} {1}_{2:5d}.2D.fits'.format(ii,pointing,id)
                ii += 1
                model.twod_spectrum(id = id,USE_REFERENCE_THUMB=True)
                model.show_2d(savePNG=True)
    
    
        
    
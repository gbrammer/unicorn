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
    unicorn.reduce.prepare_blot_reference(REF_ROOT='AEGIS_F160W', filter='F160W', REFERENCE = '/3DHST/Photometry/Work/AEGIS/IMAGES/SMALL_PIXSCL/AEGIS_F160W_sci.fits', SEGM = '/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F160W_seg.fits')
        
    NGROW=125
    pad=60
    #CATALOG='phot/F160W_arjen.cat'
    #CATALOG='/3DHST/Photometry/Work/AEGIS/Sex/PSF_matched/F140W.cat'
    CATALOG = 'phot/F140W.cat'
    direct=glob.glob('AEGIS-[0-9]-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='AEGIS_F140W'

    ##### Generate the interlaced images, including the "blotted" detection image
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False).split('-F140')[0]
        unicorn.reduce.blot_from_reference(REF_ROOT=REF_ROOT, DRZ_ROOT = pointing+'-F140W', NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine_blot(root=pointing+'-F140W', view=True, pad=60, REF_ROOT=REF_ROOT, CATALOG=CATALOG,  NGROW=NGROW, verbose=True)
        unicorn.reduce.interlace_combine(pointing+'-F140W', pad=60, NGROW=NGROW)
        unicorn.reduce.interlace_combine(pointing+'-G141', pad=60, NGROW=NGROW)

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('AEGIS-[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

    models = glob.glob('AEGIS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.)
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
            gris.fit_in_steps(dzfirst=0.005, dzsecond=0.0002,p_flat=0.0)

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
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

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
    direct=glob.glob('AEGIS-2[0-9]-F140W_asn.fits')

    extract_limit = 24
    skip_completed=False
    REF_ROOT='AEGIS_F140W'

    ##### Generate the spectral model and Extract all spectra
    inter = glob.glob('AEGIS-2[0-9]-G141_inter.fits')
    redo = False
    for i in range(len(inter)):
        time.strftime('%X %x %Z')
        pointing = inter[i].split('-G141_inter')[0]
        if (not os.path.exists(pointing+'_model.fits')) | redo:
            model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=25.8)
            model.extract_spectra_and_diagnostics(MAG_LIMIT=25.8)
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='AEGIS-10')

    skip_completed = True

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

    os.chdir(unicorn.GRISM_HOME+'COSMOS/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(REF_ROOT='COSMOS_F160W', filter='F160W', 
        REFERENCE = '/3DHST/Photometry/Work/COSMOS/IMAGES/SMALL_PIXSCL_V3/COSMOS_F160W_sci.fits', 
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
    
    ### Fit just the spec_z's
    
    id_spec_z = cat.id[np.where((cat.z_spec > 0.) & (cat.z_spec < 9.))].astype(int)
    print 'There are {0} objects in with spec_z in this field.'.format(len(id_spec_z))
    
    models = glob.glob('COSMOS-[0-9]_inter_model.fits')
    for file in models[::1]:
        pointing = file.split('_inter')[0]
        model = unicorn.reduce.process_GrismModel(pointing, MAG_LIMIT=24.5)
        print "There are {0} objects with spec_z in {1}.".format(len([id for id in id_spec_z if id in model.cat.id]),pointing)
        for id in [id for id in id_spec_z if id in model.cat.id]:
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
            
    ### Fit all grism redshifts:

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
    unicorn.reduce.prepare_blot_reference(Force=True, REF_ROOT='GOODS-N_F160W', filter='F160W', REFERENCE = '/3DHST/Photometry/Work/GOODS-N/v4/images/GOODS-N_F160W_sci_sub.fits', SEGM = '/3DHST/Photometry/Work/GOODS-N/v4/sextr/checkimages/GOODS-N_F125W_F140W_F160W.seg.fits')
     
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
    inter = glob.glob('GOODS-N-1[0-9]-G141_inter.fits')
    redo = False
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
    
    os.chdir(unicorn.GRISM_HOME+'UDS/INTERLACE_v4.0')

    #### This step is needed to strip all of the excess header keywords from the mosaic for use
    #### with `blot`.
    unicorn.reduce.prepare_blot_reference(Force=True, REF_ROOT='UDS_F160W', filter='F160W', REFERENCE = '/3DHST/Photometry/Work/UDS/v4/images/UDS_F160W_sci_sub.fits', SEGM = '/3DHST/Photometry/Work/UDS/v4/sextr/checkimages/UDS_F125W_F140W_F160W.seg.fits')
        
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
            
    ##### Extract and fit only spec-z objects
    import threedhst.catIO as catIO
    cat, zout, fout = unicorn.analysis.read_catalogs(root='UDS-11')

    skip_completed = False

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

    skip_completed = True

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

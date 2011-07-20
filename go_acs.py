"""

Run the full reduction scripts on the ACS parallels

"""
import threedhst
import threedhst.prep_flt_files

import unicorn
import glob
import os
import shutil

import pyraf
from pyraf import iraf

from threedhst.prep_flt_files import process_3dhst_pair as pair

def test():
    """
    COSMOS-23 ACS overlaps with COSMOS-25 WFC3
    """
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS/PREP_FLT')
    
    root='jbhm51'
    files=glob.glob('../RAW/%s*fits*' %(root))
    for file in files:
        os.system('cp %s .' %(file))
        os.system('gunzip %s' %(os.path.basename(file)))
        print file
    
    #### Destripe + CTE correction
    threedhst.prep_flt_files.prep_acs(force=True)
    
    os.system('rm *flt.fits')
    
    asn_direct_file = 'jbhm51010_asn.fits'
    asn_grism_file = 'jbhm51020_asn.fits'
    
    #### Copy corrected FLT files to . 
    asn = threedhst.utils.ASNFile(asn_direct_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    #
    asn = threedhst.utils.ASNFile(asn_grism_file)
    for exp in asn.exposures:
        print exp
        os.system('cp ../FIXED/%s_flt.fits . ' %(exp))
    
    #### Flag CRs, subtract background and align WCS    
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    ALIGN = '/3DHST/Spectra/Work/COSMOS/MOSAIC/COSMOS-F140w_11-05-23_sci.fits'
    
    threedhst.shifts.run_tweakshifts(asn_direct_file, verbose=True)

    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)

    threedhst.shifts.refine_shifts(ROOT_DIRECT=asn_direct_file.split('_as')[0].upper(),
              ALIGN_IMAGE=ALIGN, 
              ALIGN_EXTENSION = 0,
              fitgeometry='shift', clean=True)
    
    threedhst.prep_flt_files.startMultidrizzle(asn_direct_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=False,
        updatewcs=False, clean=True, median=False)
        
    ### Grism
    threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)

    threedhst.prep_flt_files.startMultidrizzle(asn_grism_file, use_shiftfile=True,
        skysub=True,
        final_scale=0.05, pixfrac=1, driz_cr=True,
        updatewcs=True, clean=True, median=True)
    
    ## Check
    threedhst.gmap.makeImageMap(['JBHM51010_drz.fits[1]*4', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-25-F140W_drz.fits', 'JBHM51020_drz.fits[1]', unicorn.GRISM_HOME+'COSMOS/PREP_FLT/COSMOS-25-G141_drz.fits'], aper_list=[15])
    
    ###################### Run the grism
    os.system('cp *shifts.txt *tweak.fits *asn.fits ../DATA')
    
    os.chdir(unicorn.GRISM_HOME+'ACS_PARALLEL/COSMOS')
    unicorn.go_3dhst.set_parameters(direct='F814W', LIMITING_MAGNITUDE=25)

    threedhst.process_grism.set_ACS_G800L()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/JBHM51010_drz.fits'
    #threedhst.options['PREFAB_GRISM_IMAGE'] = '../PREP_FLT/JBHM51020_drz.fits'
    
    threedhst.process_grism.reduction_script(asn_grism_file='jbhm51020_asn.fits')
    
    ### SEDs unicorn.analysis.make_SED_plots(grism_root='jbhm51020')
    #os.chdir('../')
    grism_root=asn_grism_file.split('_asn')[0]
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root, cosmos=True)
    
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root.upper(), GRISM_NAME='G800L')
        
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
    
    
    unicorn.go_3dhst.clean_up()
    
    os.system('rsync -avz HTML/ ~/Sites_GLOBAL/P/GRISM_ACS/')
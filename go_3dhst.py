import threedhst
import unicorn
import glob
import os
import shutil

import pyraf
from pyraf import iraf

def redo_all_SED_plots():
    import glob
    import threedhst

    for dir in ['COSMOS','GOODS-N','AEGIS','SN-PRIMO','SN-GEORGE'][1:]:
        os.chdir(unicorn.GRISM_HOME+dir+'/DATA/')
        files=glob.glob('*G141_asn.fits')
        os.chdir('../')
        for file in files:
            print file.split('_asn')[0]
            unicorn.analysis.make_SED_plots(grism_root=file.split('_asn')[0])
            

def goods_s():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-S')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('GOODS-S-[0-9]*-G141_asn.fits')
    files=glob.glob('GOODS-S-[0-9]*-G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('GOODS-S-[0-9]*-F140W_tweak.fits'))
    for file in files:
        shutil.copy(file, '../DATA')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=24)
    
    #### Main loop for reduction
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GOODS-S-6-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='GOODS-S-6-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='GOODS-S-6-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GOODS-S-27-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='GOODS-S-27-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='GOODS-S-27-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GOODS-S-24-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='GOODS-S-24-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='GOODS-S-24-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GOODS-S-28-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='GOODS-S-28-G141_asn.fits')
    #### HUDF09, doesn't overlap with FIREWORKS
    # unicorn.analysis.make_SED_plots(grism_root='GOODS-S-27-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GOODS-S-23-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='GOODS-S-23-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='GOODS-S-23-G141')
    go.clean_up()

def aegis():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    
    import unicorn.analysis
    
    # ######################## Test!
    # go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=21)
    # 
    # threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-11-F140W_drz.fits'
    # # threedhst.options['PREFAB_GRISM_IMAGE'] = '../PREP_FLT/AEGIS-11-G141_drz.fits'
    # proc.reduction_script(asn_grism_file='AEGIS-11-G141_asn.fits')
    # ########################
    
    os.chdir(unicorn.GRISM_HOME+'AEGIS')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('AEGIS-[0-9]*-G141_asn.fits')
    files=glob.glob('AEGIS-[0-9]*-G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('AEGIS-[0-9]*-F140W_tweak.fits'))
    for file in files:
        shutil.copy(file, '../DATA')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=24)
    
    #### Main loop for reduction
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-4-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-4-G141_asn.fits')
    unicorn.analysis.aegis_SED_plots(grism_root='AEGIS-4-G141')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-5-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-5-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-5-G141')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-11-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-11-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-11-G141')
    go.clean_up()
 
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-2-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-2-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-2-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-1-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-1-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-1-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-3-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-3-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-3-G141')
    go.clean_up()

    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/AEGIS-9-F140W_drz.fits'
    proc.reduction_script(asn_grism_file='AEGIS-9-G141_asn.fits')
    unicorn.analysis.make_SED_plots(grism_root='AEGIS-9-G141')
    go.clean_up()
    
def cosmos():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'COSMOS')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('COSMOS-[0-9]*-G141_asn.fits')
    files=glob.glob('COSMOS-[0-9]*-G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('COSMOS-[0-9]*-F140W_tweak.fits'))
    for file in files:
        shutil.copy(file, '../DATA')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=24)
    
    #threedhst.options['DRZRESOLA'] = '100.0'
    
    grism_asn = grism_asn
    
    #### Main loop for reduction
    for i in range(len(grism_asn))[16:]:
        asn = grism_asn[i]
        threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/' +  asn.replace('G141_asn','F140W_drz')
        # threedhst.options['PIXFRAC'] = 0.8
        # threedhst.options['DRZRESOLA'] = '35'
        # threedhst.options['DRZSCALE'] = '0.10'
        #### Images for a better fluxcube
        root=asn.replace('_asn.fits','')
        threedhst.options['OTHER_BANDS'] = []
        # for wave in [1.1e4,1.25e4,1.6e4]:
        #     out = unicorn.analysis.make_fluximage(grism_root=root,
        #                wavelength=wave)
        #     threedhst.options['OTHER_BANDS'].append([os.path.basename(out), 'F%03dW' %(wave/100), wave/10., 26.46])
        proc.reduction_script(asn_grism_file=asn)
        unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()

def goodsn():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'GOODS-N')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('GOODS-N-[0-9]*-G141_asn.fits')
    files=glob.glob('GOODS-N-[0-9]*-G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('GOODS-N-[0-9]*-F140W_tweak.fits'))
    for file in files:
        shutil.copy(file, '../DATA')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=24)
            
    #### Main loop for reduction
    for i in range(len(grism_asn))[11:]:
        asn=grism_asn[i]
        threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/' +  asn.replace('G141_asn','F140W_drz')
        proc.reduction_script(asn_grism_file=asn)
        unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()
   
#
def sn_primo():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'SN-PRIMO')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('PRIMO-1???-G141_asn.fits')
    files=glob.glob('PRIMO-1???-G141_shifts.txt')
    files.extend(grism_asn)
    for file in files:
        shutil.copy(file,'../DATA/')
    
    try:
        iraf.imcopy('PRIMO_F125W_drz.fits[1]', '../DATA/f125w.fits')
    except:
        os.remove('../DATA/f125w.fits')
        iraf.imcopy('PRIMO_F125W_drz.fits[1]', '../DATA/f125w.fits')
    
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F160W', LIMITING_MAGNITUDE=24.3)
    
    #### Main loop for reduction
    for i, asn in enumerate(grism_asn):
        threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/PRIMO_F160W_drz.fits'
        threedhst.options['OTHER_BANDS'] = [['f125w.fits', 'F125W' , 1248.6, 26.25]]
        proc.reduction_script(asn_grism_file=asn)
        unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()

def sn_george():
    """

    """
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'SN-GEORGE')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('GEORGE-G141_asn.fits')
    files=glob.glob('GEORGE-G141_shifts.txt')
    files.extend(grism_asn)
    for file in files:
        shutil.copy(file,'../DATA/')
    
    try:
        iraf.imcopy('GEORGE_F125W_drz.fits[1]', '../DATA/f125w.fits')
    except:
        os.remove('f125w.fits')
        iraf.imcopy('GEORGE_F125W_drz.fits[1]', '../DATA/f125w.fits')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F160W', LIMITING_MAGNITUDE=24.3)
    
    #### Main loop for reduction
    for i, asn in enumerate(grism_asn):
        threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/GEORGE_F160W_drz.fits'
        threedhst.options['OTHER_BANDS'] = [['f125w.fits', 'F125W' , 1248.6, 26.25]]
        proc.reduction_script(asn_grism_file=asn)
        unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()

def sn_marshall():
    """

    """
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis

    os.chdir(unicorn.GRISM_HOME+'SN-MARSHALL')

    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('MARSHALL-2??-G141_asn.fits')
    files=glob.glob('MARSHALL-2*-G141_shifts.txt')
    files.extend(grism_asn)
    for file in files:
        shutil.copy(file,'../DATA/')

    try:
        iraf.imcopy('MARSHALL-F125W_drz.fits[1]', '../DATA/f125w_sci.fits')
    except:
        os.remove('../DATA/f125w_sci.fits')
        iraf.imcopy('MARSHALL-F125W_drz.fits[1]', '../DATA/f125w_sci.fits')
    
    os.chdir('../')

    #### Initialize parameters
    go.set_parameters(direct='F160W', LIMITING_MAGNITUDE=25)
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/MARSHALL-F160W_drz.fits'
    threedhst.options['OTHER_BANDS'] = [['f125w_sci.fits', 'F125W' , 1248.6, 26.25]]

    #### Main loop for reduction
    for i in range(len(grism_asn)):
        asn = grism_asn[i]
        proc.reduction_script(asn_grism_file=asn)
        unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()
    
#
def daddi():
    """

    """
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis

    os.chdir(unicorn.GRISM_HOME+'DADDI')

    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('HIGHZ-CLUSTER-?-G141_asn.fits')
    files=glob.glob('HIGHZ-CLUSTER-?-G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('*tweak.fits'))
    for file in files:
        shutil.copy(file,'../DATA/')

    os.chdir('../')

    #### Initialize parameters
    go.set_parameters(direct='F140W', LIMITING_MAGNITUDE=25)
    threedhst.options['OTHER_BANDS'] = []

    #### Main loop for reduction
    for i in range(len(grism_asn)):
        asn = grism_asn[i]
        threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/'+asn.replace('asn','drz').replace('G141', 'F140W')
        proc.reduction_script(asn_grism_file=asn)
        #unicorn.analysis.make_SED_plots(grism_root=asn.split('_asn.fits')[0])
        go.clean_up()

def stanford():
    import unicorn.go_3dhst as go
    import threedhst.process_grism as proc
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'STANFORD')
    
    #### Copy necessary files from PREP_FLT to DATA
    os.chdir('PREP_FLT')
    grism_asn  = glob.glob('ISCS*G141_asn.fits')
    files=glob.glob('ISCS*G141_shifts.txt')
    files.extend(grism_asn)
    files.extend(glob.glob('ISCS*F160W_tweak.fits'))
    for file in files:
        shutil.copy(file, '../DATA')
    os.chdir('../')
    
    #### Initialize parameters
    go.set_parameters(direct='F160W', LIMITING_MAGNITUDE=23)
    
    #### Main loop for reduction
    # threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1425.3+3250-F160W_drz.fits'
    # proc.reduction_script(asn_grism_file='ISCSJ1425.3+3250-G141_asn.fits')
    # go.clean_up()
    # 
    # threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1426.5+3339-F160W_drz.fits'
    # proc.reduction_script(asn_grism_file='ISCSJ1426.5+3339-G141_asn.fits')
    # go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1429.2+3357-F160W_drz.fits'
    proc.reduction_script(asn_grism_file='ISCSJ1429.2+3357-G141_asn.fits')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1429.3+3437-F160W_drz.fits'
    proc.reduction_script(asn_grism_file='ISCSJ1429.3+3437-G141_asn.fits')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1431.1+3459-F160W_drz.fits'
    proc.reduction_script(asn_grism_file='ISCSJ1431.1+3459-G141_asn.fits')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1432.4+3250-F160W_drz.fits'
    proc.reduction_script(asn_grism_file='ISCSJ1432.4+3250-G141_asn.fits')
    go.clean_up()
    
    threedhst.options['PREFAB_DIRECT_IMAGE'] = '../PREP_FLT/ISCSJ1434.5+3427-F160W_drz.fits'
    proc.reduction_script(asn_grism_file='ISCSJ1434.5+3427-G141_asn.fits')
    go.clean_up()
    
    
def set_parameters(direct='F140W', LIMITING_MAGNITUDE=25):
    
    threedhst.defaultOptions()
    
    threedhst.options['DETECT_THRESH'] = 2
    threedhst.options['ANALYSIS_THRESH'] = 1
    threedhst.options['LIMITING_MAGNITUDE'] = LIMITING_MAGNITUDE

    threedhst.options['FULL_EXTRACTION_GEOMETRY'] = False

    #### Already processed background and shifts
    threedhst.options['OTHER_BANDS'] = []
    threedhst.options['PATH_TO_RAW'] = '../PREP_FLT/'
    threedhst.options['SKY_BACKGROUND'] = None
    threedhst.options['MAKE_SHIFTFILES'] = False
    threedhst.options['ALIGN_IMAGE'] = None

    threedhst.options['DRZRESOLA'] = '22.0'
    threedhst.options['PIX_FRAC'] = '0.8'
    threedhst.options['DRZSCALE'] = '0.06'

    threedhst.options['AXE_EDGES'] = "250,0,0,0"
    threedhst.options['USE_TAXE'] = True

    #### Use F140W as detection image
    threedhst.options['MAG_ZEROPOINT'] = 26.46
    threedhst.options['FILTWAVE'] = 1392.
        
    #### Use F125W as detection image
    if direct == 'F125W':
        threedhst.options['MAG_ZEROPOINT'] = 26.25
        threedhst.options['FILTWAVE'] = 1248.6
    
    #### Use F160W as detection image
    if direct == 'F160W':
        threedhst.options['MAG_ZEROPOINT'] = 25.96
        threedhst.options['FILTWAVE'] = 1537.
    
def clean_up():
    files=glob.glob('OUTPUT_G141/*')
    files.extend(glob.glob('DATA/*FLX*'))
    files.extend(glob.glob('PREP_FLT/*FLX*'))
    files.extend(glob.glob('DRIZZLE_G141/*mef*'))
    files.extend(glob.glob('DRIZZLE_G141/*PET*'))
    files.extend(glob.glob('DATA/threedhst_auto*swarp'))
    files.extend(glob.glob('DATA/*coeffs?.dat'))
    files.extend(glob.glob('HTML/ascii/*.FITS'))

    for file in files:
        os.remove(file)
       
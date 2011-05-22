#!/usr/bin/env python
# encoding: utf-8
"""
prepare.py

Scripts to "prepare" the FLT files for the 3D-HST reduction:

Direct images:
    1) Align to WCS reference
    2) Subtract sky backgrounds
    3) Run multidrizzle to flag cosmic rays
    4) Make mosaics combining pointings for a given field
    
Grism images:
    1) Copy shifts from direct images
    2) Divide by the F140W flat
    3) Subtract sky background
    4) Run multidrizzle to flag cosmic rays
    5) ...
    
$URL$
$Author$
$Date$

"""
__version__ = " $Rev$"

def COSMOS(FORCE=False):
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    import glob
    import os
    
    os.chdir('/3DHST/Spectra/Work/COSMOS/PREP_FLT')
    ALIGN = '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits'
    #ALIGN = '/3DHST/Ancillary/COSMOS/NMBS/COSMOS-1_K_sci.fits'
    
    #### Direct images only
    direct=glob.glob('*30_asn.fits')
    grism = glob.glob('*40_asn.fits')
    for i in range(len(direct)):
        pointing=threedhst.prep_flt_files.make_targname_asn(direct[i], newfile=False)
        if (not os.path.exists(pointing)) | FORCE:
            pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=True)
        
    # # Check alignment
    # threedhst.gmap.makeImageMap(['COSMOS-27-F140W_drz.fits', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.01'], aper_list=[16])
    # 
    # threedhst.gmap.makeImageMap(['COSMOS-27-F140W_drz.fits','/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'], aper_list=[16])

    #### Direct mosaic
    direct_files = glob.glob('COSMOS-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.06
    #SCALE = 0.5
    NX, NY = int(9355*0.06/SCALE), int(11501*0.06/SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.6, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12356, dec=2.3608425,
             final_outnx = NX, final_outny=NY)
    #
    #### Temporary release of the direct mosaic
    # iraf.imcopy('COSMOS-F140W_drz.fits[1]','COSMOS-F140w_11-05-20_sci.fits')
    # iraf.imcopy('COSMOS-F140W_drz.fits[2]','COSMOS-F140w_11-05-20_wht.fits')
    # !tar czvf COSMOS-F140w_11-05-20.tar.gz COSMOS-*-F140W_shifts.txt COSMOS-*-F140W_tweak.fits COSMOS-*-F140W_asn.fits COSMOS-F140W_shifts.txt COSMOS-F140W_asn.fits
    # !mv COSMOS-F140w_11-05-20* ../MOSAIC
    
    threedhst.gmap.makeImageMap(['COSMOS-F140W_drz.fits', '/3DHST/Ancillary/COSMOS/WIRDS/WIRDS_Ks_100028+021230_T0002.fits[0]*0.01','/3DHST/Ancillary/COSMOS/ACS/acs_I_030mas_077_sci.fits[0]*3'], aper_list=[13,14,15,16,17], tileroot=['F140W','WIRDS-K','F814W'])
    
    
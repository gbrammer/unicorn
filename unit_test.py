"""

Start-to-finish test from preparation steps to fitting redshifts.

"""

def misc():
    """
    Make cutouts of the large 3D-HST mosaic / segm. image for reference.
    
    Put in the "Catalog" directory
    
    """
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/GOODS-S/UnitTest/Misc')
    iraf.imcopy('/Users/brammer/3DHST/Spectra/Work/ASTRODRIZZLE_FLT/Catalog/goodss_3dhst.v4.0.IR_orig_sci.fits[5976:11034, 9396:14135]', 'goodss_3dhst.v4.0.IR_cutout_sci.fits')
    iraf.imcopy('/Users/brammer/3DHST/Spectra/Work/ASTRODRIZZLE_FLT/Catalog/GOODS-S_IR.seg.fits[5976:11034, 9396:14135]', 'goodss_3dhst.v4.0.IR_cutout_seg.fits')
    
    os.system('cp /Users/brammer/3DHST/Spectra/Work/ASTRODRIZZLE_FLT/Catalog/GOODS-S_IR.cat goodss_3dhst.v4.0.IR_cutout.cat')
        
    os.system('cp /Users/brammer/3DHST/Spectra/Work/ASTRODRIZZLE_FLT/Catalog/goodss_radec.dat .')
    
def run_wfc3(run_prep=True, run_interlace=True, run_redshifts=True, align_threshold=6):
    """
    Run end-to-end WFC3 test
    """
    import os
    import glob
    import time
    
    import astropy.io.fits as pyfits
    from astropy.table import Table as table
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init
    
    import unicorn
    import unicorn.interlace_test as test
    
    #### Preparation steps
    #run_prep = True
    
    root = 'goodss-34'
    
    if run_prep:
        
        #### Make ASN files
        unicorn.candels.make_asn_files(uniquename=False)
        
        #### Run main preparation wrapper to align to reference and 
        #### subtract backgrounds
        init.prep_direct_grism_pair(direct_asn=root+'-F140W_asn.fits', grism_asn=root+'-G141_asn.fits', radec='../Catalog/goodss_radec.dat', raw_path='../RAW/', mask_grow=8, scattered_light=False, final_scale=None, skip_direct=False, ACS=False, align_threshold=align_threshold)
        
        #### test
        # if True:
        #     import threedhst.prep_flt_astrodrizzle as init
        #     init.subtract_flt_background(root=root+'-F140W', scattered_light=False)
                    
    
    #### Interlaced images + reference
    #run_interlace = True    
    if run_interlace:
        
        #### Interlace images themselves
        NGROWX = 180
        NGROWY = 30
        grow, auto_off = 2, False
        #grow, auto_off = 1, True
        
        unicorn.reduce.interlace_combine(root=root+'-F140W', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)
        unicorn.reduce.interlace_combine(root=root+'-G141', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)
        
        #### Interlaced reference
        adriz_blot = unicorn.reduce.adriz_blot_from_reference
        adriz_blot(pointing=root+'-F140W', pad=60, NGROWX=NGROWX, NGROWY=NGROWY, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0, ref_image='../Catalog/goodss_3dhst.v4.0.IR_cutout_sci.fits', ref_ext=0, ref_filter='F140W', seg_image='../Catalog/goodss_3dhst.v4.0.IR_cutout_seg.fits', cat_file='../Catalog/goodss_3dhst.v4.0.IR_cutout.cat')
        
        #### Make a model to 25 AB
        model = unicorn.reduce.process_GrismModel(root=root, grow_factor=grow, growx=grow, growy=grow, MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, model_list=None, direct='F140W', grism='G141', BEAMS=['A', 'B', 'C', 'D', 'E'], align_reference=True)
        
        #### Remove residual backgroudn with model mask
        model.refine_mask_background(threshold=0.002, grow_mask=14, update=True, resid_threshold=4, clip_left=640, save_figure=True, interlace=True)
        
    #run_redshifts = True
    
    if run_redshifts:
        
        model = unicorn.reduce.process_GrismModel(root=root, grow_factor=2, growx=2, growy=2, MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, model_list=None, direct='F140W', grism='G141', BEAMS=['A', 'B', 'C', 'D', 'E'])
        
        #### Some objects with z_spec, em lines; v4.0 IDs
        ids = [25778, 26087, 26550, 28052, 28870, 29122, 30520, 30662]
        for id in ids:
            if os.path.exists('%s_%05d.2D.fits' %(model.baseroot, id)):
                continue
            #
            model.twod_spectrum(id=id, grow=1, miny=-36, maxy=None, CONTAMINATING_MAGLIMIT=23, refine=False, verbose=False, force_refine_nearby=False, USE_REFERENCE_THUMB=True, USE_FLUX_RADIUS_SCALE=3, BIG_THUMB=False, extract_1d=True)
            #
            ### Original fit
            gris = unicorn.interlace_fit.GrismSpectrumFit(root='%s_%05d' %(model.baseroot, id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.4, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True, dr_match=1.0)
            gris.fit_in_steps(dzfirst=0.003, zrfirst=(0.6, 3.3), dzsecond=0.0002, save=True, make_plot=True, skip_second=False, oned=False)
            ### New fit
            gris = test.SimultaneousFit(root='%s_%05d' %(model.baseroot, id), FIGURE_FORMAT='png', lowz_thresh=0.4)
            gris.new_fit_constrained(zrfirst=[0.6, 3.3], dzfirst=0.005, dzsecond=0.0005, make_plot=True, ignore_photometry=False, ignore_spectrum=False, refit_norm=False, get_tilt=True, faint_limit=24)
        
def run_acs(run_prep=True, run_interlace=True, run_redshifts=True):
    import glob
    
    import threedhst
    import threedhst.prep_flt_astrodrizzle as init

    import unicorn
    import unicorn.interlace_acs
    
    ### Make ACS associations
    files=glob.glob('../RAW/j*asn.fits')
    for file in files:
        pointing=threedhst.prep_flt_files.make_targname_asn(file, newfile=True, field='goodss-34', ext='flc')
    
    if run_prep:
        ### Alignment and background subtraction
        init.prep_direct_grism_pair(direct_asn='goodss-34-03-F814W_asn.fits', grism_asn='goodss-34-03-G800L_asn.fits', radec='../Catalog/goodss_radec.dat', ACS=True)
        
    if run_interlace:
        #### Interlace images
        interlace_acs = unicorn.interlace_acs.interlace_combine_acs
        for chip in [1,2]:
            interlace_acs(root='goodss-34-03-F814W', view=False, chip=chip, filter='F814W', outroot='goodss-34-03', center=True, pad=100, growx=1, growy=1)
            interlace_acs(root='goodss-34-03-G800L', view=False, chip=chip, filter='G800L', outroot='goodss-34-03', center=True, pad=100, growx=1, growy=1)
        
        #### Interlaced reference, self as reference image
        adriz_blot = unicorn.reduce.adriz_blot_from_reference
        adriz_blot(pointing='goodss-34-03-F814W', ACS=True, pad=100, NGROW=0, growx=1, growy=1, auto_offsets=False, ref_exp=0, ref_image='goodss-34-03-F814W_drc_sci.fits', ref_ext=0, ref_filter='F814W', seg_image='Catalog/goodss_3dhst.v4.0.IR_cutout_seg.fits', cat_file='Catalog/goodss_3dhst.v4.0.IR_cutout.cat')
        
        ### test blot from UDF reference [works]
        #adriz_blot(pointing='goodss-34-03-F814W', ACS=True, pad=100, NGROW=0, growx=1, growy=1, auto_offsets=False, ref_exp=0, ref_image='hlsp_xdf_hst_acswfc-60mas_hudf_f814w_v1_sci.fits', ref_ext=0, ref_filter='F140W', seg_image='Catalog/goodss_3dhst.v4.0.IR_cutout_seg.fits', cat_file='Catalog/goodss_3dhst.v4.0.IR_cutout.cat')
        
    if run_redshifts:
        import unicorn.interlace_test
        ids = {1:[27549, 30367, 30861], 2:[]}
        for chip in ids.keys():
            model = unicorn.reduce.process_GrismModel(root='goodss-34-03-chip%d' %(chip), grow_factor=1, growx=1, growy=1, MAG_LIMIT=21, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, model_list=None, direct='F814W', grism='G800L', BEAMS=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
            # 
            # model = unicorn.reduce.GrismModel('goodss-34-03-chip%d' %(chip), grow_factor=1, growx=1, growy=1, MAG_LIMIT=28, use_segm=False, grism='G800L', direct='F814W')
            # MLIM = [21,21]
            # model.compute_full_model(refine=False, MAG_LIMIT=MLIM[0], save_pickle=False, model_slope=0, BEAMS=['A','B','C','D','E','F','G'])   
            # ### For the brighter galaxies, refine the model with the observed spectrum         
            # model.compute_full_model(refine=True, MAG_LIMIT=MLIM[1], save_pickle=True, model_slope=0, BEAMS=['A','B','C','D','E','F','G'])
            #
            for id in ids[chip]:
                model.twod_spectrum(id=id, grow=1, miny=-40, maxy=None, CONTAMINATING_MAGLIMIT=23, refine=False, verbose=False, force_refine_nearby=False, USE_REFERENCE_THUMB=True, USE_FLUX_RADIUS_SCALE=3, BIG_THUMB=False, extract_1d=True)
                model.show_2d(savePNG=True, verbose=False)
                #
                gris = unicorn.interlace_fit.GrismSpectrumFit(root='%s_%05d' %(model.baseroot, id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.05, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True, dr_match=1.0)
                #
                gris = unicorn.interlace_test.SimultaneousFit(root='%s_%05d' %(model.baseroot, id), FIGURE_FORMAT='png', verbose=True, lowz_thresh=0.05, fix_direct_thumbnail=True, RELEASE=False, OUTPUT_PATH='./', BASE_PATH='./', skip_photometric=False, p_flat=0.0001, use_mag_prior=True, dr_match=1.0)
                gris.new_fit_constrained(zrfirst=[0.05, 0.8], dzfirst=0.005, dzsecond=0.0005, make_plot=True, ignore_photometry=False, ignore_spectrum=False, refit_norm=False, get_tilt=True, faint_limit=24)
                
        
        
        
    
    
    
    
        
        
        
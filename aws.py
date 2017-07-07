"""
Testing Amazon Web Services for the 3D-HST reductions.

"""

import threedhst
import unicorn
import os
import glob
import multiprocessing


def run_test():

    import unicorn.unit_test
    import threedhst

    os.system('mkdir /home/ubuntu/AWS_TEST/REDUCE')
    os.chdir('/home/ubuntu/AWS_TEST/REDUCE')
    os.system('rsync -av ../RAW/*asn* .')
    os.system('sh '+threedhst.__path__[0]+'/bin/flt_info.sh')
    
    unicorn.unit_test.run_wfc3(run_prep=True, run_interlace=False, run_redshifts=False, align_threshold=6)
    
    unicorn.unit_test.run_wfc3(run_prep=False, run_interlace=True, run_redshifts=False)

    unicorn.unit_test.run_wfc3(run_prep=False, run_interlace=False, run_redshifts=True)


def multi_test(field='goodss', n_proc=28):
    
    """
    from multiprocessing import Pool
    p = Pool(8) 
    pixVals = p.map(compute_pixel_value,neighborInds.size)    
    """
    
    from multiprocessing import Pool
    import glob
    roots = [r.split('-F140W')[0] for r in glob.glob('{}-*F140W_asn.fits'.format(field))]
    roots.sort()
        
    p1 = Pool(processes=n_proc)
    p1.map(interlace_func, roots)
    p1.close()
    
    for root in roots:
        adriz_func(root)
        model_func(root)
        
    p2 = Pool(processes=n_proc)
    p2.map(reduce_func, roots)
    p2.close()

def interlace_func(root):
   
    NGROWX = 200
    NGROWY = 30 
    grow= 2
    
    if root in ['aegis-01','goodsn-12', 'goodsn-13', 'goodsn-21', 'goodsn-22', 'goodsn-24', 'goodsn-31','goodsn-111','goodsn-114','goodsn-123','goodss-15']:
        auto_off = True
    else:
        auto_off = False
        
    unicorn.reduce.interlace_combine(root=root+'-F140W', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)
    unicorn.reduce.interlace_combine(root=root+'-G141', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)

def adriz_func(root):
    
    NGROWX = 200
    NGROWY = 30 
    grow = 2
    
    if root in ['aegis-01','goodsn-12', 'goodsn-13', 'goodsn-21', 'goodsn-22', 'goodsn-24', 'goodsn-31','goodsn-111','goodsn-114','goodsn-123','goodss-15']:
        auto_off = True
    else:
        auto_off = False

    if root.startswith('aegis'):
        CATALOG = '../REF/aegis_3dhst.v4.0.IR_orig.cat'
        REF_IMAGE = '../REF/aegis_3dhst.v4.0.IR_orig_sci.fits'
        SEG_IMAGE = '../REF/F160W_seg.fits'

    if root.startswith('cosmos'):
        CATALOG = '../REF/cosmos_3dhst.v4.0.IR_orig.cat'
        REF_IMAGE = '../REF/cosmos_3dhst.v4.0.IR_orig_sci.fits'
        SEG_IMAGE = '../REF/F160W_seg.fits'

    if root.startswith('goodsn'):
        CATALOG = '../REF/GOODS-N_IR.cat'
        REF_IMAGE = '../REF/goodsn_3dhst.v4.0.IR_orig_sci.fits'
        SEG_IMAGE = '../REF/goodsn_3dhst.v4.0.F160W_seg.fits'
        
    if root.startswith('goodss'):
        CATALOG = '../REF/GOODS-S_IR.cat'
        REF_IMAGE = '../REF/goodss_3dhst.v4.0.IR_orig_sci.fits'
        SEG_IMAGE = '../REF/GOODS-S_IR.seg.fits'
        
    if root.startswith('uds'):
        CATALOG = '../REF/UDS_IR.cat'
        REF_IMAGE = '../REF/uds_3dhst.v4.0.IR_orig_sci.fits'
        SEG_IMAGE = '../REF/UDS_F125W_F140W_F160W.seg.fits'
            
    #### Interlaced reference
    adriz_blot = unicorn.reduce.adriz_blot_from_reference
    adriz_blot(pointing=root+'-F140W', pad=60, NGROWX=NGROWX, NGROWY=NGROWY, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0, ref_image=REF_IMAGE, ref_ext=0, ref_filter='F140W', seg_image=SEG_IMAGE, cat_file=CATALOG)
      
def model_func(root):
    
    if root == 'goodsn-111':
        align_reference = False
    else:
        align_reference = True
    
    m0 = unicorn.reduce.GrismModel(root=root)
    model_list = m0.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
    model = unicorn.reduce.process_GrismModel(root=root, model_list=model_list, grow_factor=2, growx=2, growy=2, MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141', BEAMS=['A', 'B', 'C', 'D', 'E'], align_reference=align_reference)
    if not os.path.exists(os.path.basename(model.root) + '-G141_maskbg.dat'):
        model.refine_mask_background(threshold=0.002, grow_mask=14, update=True, resid_threshold=4, clip_left=640, save_figure=True, interlace=True)
    

def reduce_func(root):
      
    import unicorn
    import numpy as np
    import unicorn.interlace_test as test
    import os
      
    if root == 'goodsn-111':
        align_reference = False
    else:
        align_reference = True  
      
    model = unicorn.reduce.process_GrismModel(root=root, align_reference=align_reference)
    model.mask_zeroth()
     
    model.extract_spectra_and_diagnostics(MAG_LIMIT=26.)    

    cat, zout, fout = unicorn.analysis.read_catalogs(root=root)
    
    ii = np.where(model.cat.mag < 26.)
    for id in model.cat.id[ii]:
        obj_root='{}-G141_{:05d}'.format(root, id)
        if not os.path.exists(obj_root+'.2D.fits'):
            status = model.twod_spectrum(id)
            if not status:
                continue
        if not os.path.exists('{}-G141-big_{:05d}.2D.fits'.format(root, id)):
            status_big = model.twod_spectrum(id, verbose=False, miny=80, USE_FLUX_RADIUS_SCALE=3,
                USE_REFERENCE_THUMB=True, 
                BIG_THUMB=True, extract_1d=False)                                
        try:
            gris = test.SimultaneousFit(obj_root,lowz_thresh=0.01, FIGURE_FORMAT='png') 
        except:
            continue
        #
        if gris.status is False:
            continue
        #
        print '\n'
        if not os.path.exists(obj_root+'.new_zfit.pz.fits'):
            try:
                gris.new_fit_constrained()
                gris.new_save_results()
                gris.make_2d_model()
            except:
                continue
        if os.path.exists('{}-G141-big_{:05d}.2D.fits'.format(root, id)) and not os.path.exists('{}-G141-big_{:05d}.new_zfit.fits'.format(root, id)):
            os.system('cp {}.1D.fits {}-G141-big_{:05d}.1D.fits'.format(obj_root, root, id))
            os.system('cp {}.new_zfit.pz.fits {}-G141-big_{:05d}.new_zfit.pz.fits'.format(obj_root,root, id))
            gris_big = test.SimultaneousFit('{}-G141-big_{:05d}'.format(root, id))
            gris_big.make_2d_model(base='new_zfit', write_fits = True)
            os.system('rm {0}-G141-big_{1:05d}.1D.fits {0}-G141-big_{1:05d}.new_zfit.pz.fits'.format(root, id))
        if not os.path.exists(obj_root+'.linefit.fits'):
            try:
                gris.new_fit_free_emlines(ztry=None, NSTEP=600)
            except:
                continue
                
    ii = np.where(model.cat.mag >= 26.)
    for id in model.cat.id[ii]:
        obj_root='{}-G141_{:05d}'.format(root, id)
        status = model.twod_spectrum(id, miny=40)
        if not status:
            continue
        try:
            gris = test.SimultaneousFit(obj_root,lowz_thresh=0.01, FIGURE_FORMAT='png') 
        except:
            continue
        #
        if gris.status is False:
            continue
        #
        print '\n'
        if not os.path.exists(obj_root+'.new_zfit.pz.fits'):
            try:
                gris.new_fit_constrained()
                gris.new_save_results()
                gris.make_2d_model()
            except:
                continue
        if not os.path.exists(obj_root+'.linefit.fits'):
            try:
                gris.new_fit_free_emlines(ztry=None, NSTEP=600)
            except:
                continue
                
                
def extract_func(root):
    
    import unicorn
    import numpy as np
    import unicorn.interlace_test as test
    import os
      
    if root == 'goodsn-111':
        align_reference = False
    else:
        align_reference = True        
      
    model = unicorn.reduce.process_GrismModel(root=root, align_reference=align_reference)
    model.mask_zeroth()
     
    model.extract_spectra_and_diagnostics(MAG_LIMIT=30.)                
        
def time_test():
    
    import unicorn.interlace_test as test
    
    obj_root = 'goodss-34-G141_25778'
    gris = test.SimultaneousFit(obj_root,lowz_thresh=0.01, FIGURE_FORMAT='png')
    gris.new_fit_constrained()
    
def fit_timer(n_repeats=3, n_runs=1):
    """
    hyperion: [64.93909096717834, 62.57095193862915, 67.32873106002808] 64.9462579886
    unicorn: [97.09501099586487, 93.94038319587708, 93.99588394165039] 95.0104260445
    aws: [97.73548483848572, 97.6720449924469, 97.60361194610596] 97.6703805923
    """
    
    import timeit
    import numpy as np
    
    t = timeit.Timer("unicorn.aws.time_test()","import unicorn")
    times = t.repeat(n_repeats, n_runs)
    print times, np.mean(times)

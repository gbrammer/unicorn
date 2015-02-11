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


def multi_test(field='goodss', n_proc=32):
    
    """
    from multiprocessing import Pool
    p = Pool(8) 
    pixVals = p.map(compute_pixel_value,neighborInds.size)    
    """
    
    from multiprocessing import Pool
    import glob
    roots = [r.split('-F140W')[0] for r in glob.glob('{}-*F140W_asn.fits'.format(field))]

    for root in roots:
        interlace_func(root)
        
    p = Pool(processes=n_proc)
    p.map(reduce_func, roots)


def interlace_func(root):

    import unicorn
    import numpy as np
    import unicorn.interlace_test as test
    import os
   
    NGROWX = 200
    NGROWY = 30 
    grow, auto_off = 2, False

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
        
    unicorn.reduce.interlace_combine(root=root+'-F140W', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)
    unicorn.reduce.interlace_combine(root=root+'-G141', view=False, use_error=True, make_undistorted=False, pad=60, NGROWX=NGROWX, NGROWY=NGROWY, ddx=0, ddy=0, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0)
        
    #### Interlaced reference
    adriz_blot = unicorn.reduce.adriz_blot_from_reference
    adriz_blot(pointing=root+'-F140W', pad=60, NGROWX=NGROWX, NGROWY=NGROWY, growx=grow, growy=grow, auto_offsets=auto_off, ref_exp=0, ref_image=REF_IMAGE, ref_ext=0, ref_filter='F140W', seg_image=SEG_IMAGE, cat_file=CATALOG)
      
def reduce_func(root):
      
    m0 = unicorn.reduce.GrismModel(root=root)
    model_list = m0.get_eazy_templates(dr_min=0.5, MAG_LIMIT=25)
    model = unicorn.reduce.process_GrismModel(root=root, model_list=model_list, grow_factor=2, growx=2, growy=2, MAG_LIMIT=25, REFINE_MAG_LIMIT=21, make_zeroth_model=False, use_segm=False, model_slope=0, direct='F140W', grism='G141', BEAMS=['A', 'B', 'C', 'D', 'E'])
    model.mask_zeroth()
     
    model.extract_spectra_and_diagnostics(MAG_LIMIT=24.)    

    cat, zout, fout = unicorn.analysis.read_catalogs(root=root)
    
    ii = np.where(model.cat.mag < 24.)
    for id in model.cat.id[ii]:
        obj_root='%s-G141_%05d' %(root, id)
        if os.path.exists(obj_root+'.new_zfit.pz.fits'):
            continue
        if not os.path.exists(obj_root+'.2D.fits'):
            status = model.twod_spectrum(id)
            print status
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
        try:
            gris.new_fit_constrained()
        except:
            continue
        
def time_test():
    
    import unicorn.interlace_test as test
    
    obj_root = 'goodss-34-G141_25778'
    gris = test.SimultaneousFit(obj_root,lowz_thresh=0.01, FIGURE_FORMAT='png')
    gris.new_fit_constrained()
    
def fit_timer(n_repeats=3, n_runs=1):
    """
    My desktop: ~??s
    """
    
    import timeit
    import numpy as np
    
    t = timeit.Timer("unicorn.aws.time_test()","import unicorn")
    times = t.repeat(n_repeats, n_runs)
    print times, np.mean(times)

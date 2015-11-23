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

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

HAS_PHOTOMETRY = True
PHOTOMETRY_ID = None
BAD_SPECTRUM = False

SPC_FILENAME = None
SPC = None
 
def get_grism_path(root):
    """ 
    Given a rootname for a grism pointing, get the path to the 
    grism files
    """
    PATH = './'
    if root.startswith('COSMOS'):
        PATH = unicorn.GRISM_HOME+'COSMOS/'
    if root.startswith('AEGIS'):
        PATH = unicorn.GRISM_HOME+'AEGIS/'
    if root.startswith('GOODS-N'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('UDS'):
        PATH = unicorn.GRISM_HOME+'UDS/'
    if root.startswith('GN20'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('G850.1'):
        PATH = unicorn.GRISM_HOME+'GOODS-N/'
    if root.startswith('GOODS-S'):
        PATH = unicorn.GRISM_HOME+'GOODS-S/'
    if root.startswith('UDF'):
        PATH = unicorn.GRISM_HOME+'UDF/'
    if root.startswith('WFC3-ERS'):
        PATH = unicorn.GRISM_HOME+'ERS/'
    if root.startswith('MARSHALL'):
        PATH = unicorn.GRISM_HOME+'SN-MARSHALL/'
    if root.startswith('PRIMO'):
        PATH = unicorn.GRISM_HOME+'SN-PRIMO/'
    if root.startswith('GEORGE'):
        PATH = unicorn.GRISM_HOME+'SN-GEORGE/'
    if root.startswith('TILE41'):
        PATH = unicorn.GRISM_HOME+'SN-TILE41/'
    if root.startswith('EGS1'):
        PATH = unicorn.GRISM_HOME+'COOPER/'
    if root.startswith('COLFAX'):
        PATH = unicorn.GRISM_HOME+'SN-COLFAX/'
    for fi in ['HS0105', 'HS1603', 'Q0100', 'Q0142', 'Q0207', 'Q0449', 'Q0821', 'Q1009', 'Q1217', 'Q1442', 'Q1549', 'Q1623', 'Q1700', 'Q2206', 'Q2343']:
        if fi in root:
            PATH =  unicorn.GRISM_HOME+'Erb/'
            
    return PATH
    
def read_catalogs(root='', cosmos=False, aegis=False, goodsn=False, cdfs=False, ecdfs=False, uds=False, udf=False, aegis_wirds=False):
    """
    
    Read photometry, redshift, SPS catalogs for a given field.
    
    Usage: cat, zout, fout = read_catalogs(cosmos=True)
    
    """        
    import threedhst.catIO as catIO
    import numpy as np
    
    KTOT_COL = None
    ### Optionally supply a grism basename and determine the field from that
    if root.upper().startswith('UDS'):
        uds=True
    if root.upper().startswith('COSMOS'):
        cosmos=True
    if root.upper().startswith('AEGIS'):
        aegis=True
    if root.startswith('GOODS-N') | root.startswith('goodsn'):
        goodsn=True
    if root.startswith('GN20') | root.startswith('G850.1'):
        goodsn=True
    if root.startswith('GOODS-S') | root.startswith('goodss') | root.startswith('gs'):
        cdfs=True

    if root.startswith('CDFS'):
        cdfs=True
    
    #### M. Cooper
    if root.startswith('EGS1'):
        aegis=True
    
    #### G. Barro
    if 'gdn' in root.lower():
        goodsn=True
    #
    if 'GN-z10' in root:
        goodsn=True
    
    ### CLEAR
    if root.startswith('GD'):
        goodsn=True
    
    if root.startswith('GS'):
        cdfs=True
    
        
    if root.startswith('GOODS-S-34') | root.startswith('GOODS-S-36') | root.startswith('GOODS-S-37') | root.startswith('GOODS-S-38') | root.startswith('UDF') | root.startswith('HUDF'):
        udf=True
    
    #### Supernova pointings    
    if root.startswith('MARSHALL'):
        uds=True
    if root.startswith('PRIMO'):
        cdfs=True
    if root.startswith('GEORGE'):
        cdfs=True
    if root.startswith('TILE41'):
        cosmos=True
    if root.startswith('COLFAX'):
        goodsn=True
    
    if 'FIGS-GS' in root:
        cdfs=True
    
    if 'FIGS-GN' in root:
        goodsn=True
        
    #### ERS
    if root.startswith('WFC3-ERS'):
        cdfs=True
    
    #### Papovich cluster
    if root.startswith('IRC0218'):
        uds=True
        
    aegis_wirds=False
    #if root.startswith('AEGIS-11') | root.startswith('AEGIS-2-') | root.startswith('AEGIS-1-') | root.startswith('AEGIS-6-'):
    #    aegis=False
    #    aegis_wirds=True
        
    uds_cluster = False
    if root.startswith('IRC0218'):
        uds=True

    erb = False
    for fi in ['HS0105', 'HS1603', 'Q0100', 'Q0142', 'Q0207', 'Q0449', 'Q0821', 'Q1009', 'Q1217', 'Q1442', 'Q1549', 'Q1623', 'Q1700', 'Q2206', 'Q2343']:
        if fi in root:
            erb = True
    
    CAT_FILE = None
    ZOUT_FILE = None
    FOUT_FILE = None
    
    MAGS = False
    
    #### Define paths
    if cosmos:
        
        ### NMBS
        GRISM_PATH=unicorn.GRISM_HOME+'COSMOS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'        
        CAT_FILE = CAT_PATH + '../cosmos-1.v4.6.cat'
        ZOUT_FILE = CAT_PATH + 'cosmos-1.v4.6.zout'
        FOUT_FILE = CAT_PATH+'../cosmos-1.bc03.v4.6.fout'
        KTOT_COL = 'ktot'
                
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            CAT_PATH = '/3DHST/Ancillary/COSMOS/NMBS/Photometry/'
            CAT_FILE = CAT_PATH + 'cosmos-1.deblend.v5.1.cat'
            ZOUT_FILE = CAT_PATH + '/cosmos-1.deblend.redshifts/cosmos-1.deblend.v5.1.zout'
            FOUT_FILE = CAT_PATH+'cosmos-1.deblend.sps/cosmos-1.bc03.del.deblend.v5.1.fout'
            KTOT_COL = 'K'

        # For the v4.1 reduction on Unicorn
        
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('hyp'):
            CAT_PATH = '/3DHST/Photometry/Release/v4.1/COSMOS/'            
            CAT_FILE = CAT_PATH + 'Catalog/cosmos_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/cosmos_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH+'Fast/cosmos_3dhst.v4.1.fout'
            KTOT_COL = 'f_k'
                     
        #### ZFOURGE
        if unicorn.hostname().startswith('850dhcp'):
            CAT_PATH = unicorn.GRISM_HOME+'COSMOS/ZFOURGE/gabe/'        
            CAT_FILE = CAT_PATH + 'cosmos-1.v0.7.2.cat'
            ZOUT_FILE = CAT_PATH + 'cosmos-1.v0.7.2.zout'
            FOUT_FILE = CAT_PATH+'cosmos-1.v0.7.2.fout'
            KTOT_COL = 'Kstot'
        
        #### Gabe's laptop
        if ('brammer' in unicorn.hostname().lower()) | ('gabriel' in unicorn.hostname().lower()):
            # CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.0/'
            # CAT_FILE = CAT_PATH + 'Eazy/EazyRerun/cosmos_3dhst.v4.0.cat'
            # ZOUT_FILE = CAT_PATH + 'Eazy/EazyRerun/OUTPUT/cosmos.dusty3.zout'
            # FOUT_FILE = CAT_PATH + 'Fast/cosmos_3dhst.v4.0.fout'
            # KTOT_COL = 'f_F160W'
            
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/cosmos_3dhst.v4.1.cats/'
            CAT_FILE = CAT_PATH + 'Catalog/cosmos_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/cosmos_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH + 'Fast/cosmos_3dhst.v4.1.fout'
            KTOT_COL = 'f_F160W'
            
            ### Use dusty-template fit
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun'
            ZOUT_FILE = CAT_PATH + '/OUTPUT/cosmos_3dhst.v4.1.dusty.zout'
            
    if aegis:
        GRISM_PATH=unicorn.GRISM_HOME+'AEGIS/'
        
        #### NMBS
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v5.1/'
        CAT_FILE = CAT_PATH + 'aegis-n2.deblend.v5.1.cat'
        ZOUT_FILE = CAT_PATH + '/aegis-n2.deblend.redshifts/aegis-n2.deblend.v5.1.zout'
        FOUT_FILE = CAT_PATH+'aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout'
        KTOT_COL = 'K'

        CAT_PATH = os.getenv('THREEDHST')+'/Spectra/Release/v4.0/'
        CAT_FILE = CAT_PATH + 'Catalog/aegis_3dhst.v4.0.cat'
        ZOUT_FILE = CAT_PATH + 'Eazy/aegis_3dhst.v4.0.zout'
        FOUT_FILE = CAT_PATH+'Fast/aegis_3dhst.v4.0.fout'
        KTOT_COL = 'f_F160W'
        
        # For the v4.1 reduction
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('hyp'):
            CAT_PATH = '/3DHST/Photometry/Release/v4.1/AEGIS/'
            CAT_FILE = CAT_PATH + 'Catalog/aegis_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/aegis_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH+'Fast/aegis_3dhst.v4.1.fout'
            KTOT_COL = 'f_k'
        
        #### Gabe's laptop
        if ('brammer' in unicorn.hostname().lower()) | ('gabriel' in unicorn.hostname().lower()):
            # CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.0/'
            # CAT_FILE = CAT_PATH + 'Eazy/EazyRerun/aegis_3dhst.v4.0.cat'
            # ZOUT_FILE = CAT_PATH + 'Eazy/EazyRerun/OUTPUT/aegis.dusty3.zout'
            # FOUT_FILE = CAT_PATH + 'Fast/aegis_3dhst.v4.0.fout'
            # KTOT_COL = 'f_F160W'
            
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/aegis_3dhst.v4.1.cats/'
            CAT_FILE = CAT_PATH + 'Catalog/aegis_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/aegis_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH + 'Fast/aegis_3dhst.v4.1.fout'
            KTOT_COL = 'f_F160W'
            
            ### Use dusty-template fit
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun'
            ZOUT_FILE = CAT_PATH + '/OUTPUT/aegis_3dhst.v4.1.dusty.zout'
            
    if aegis_wirds:
        GRISM_PATH=unicorn.GRISM_HOME+'AEGIS/'
        CAT_PATH = '/Users/gbrammer/research/drg/PHOTZ/EAZY/WIRDS/FAST/'
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            CAT_PATH = '/3DHST/Ancillary/AEGIS/WIRDS/Photometry/FAST/'
        #
        ZOUT_FILE = CAT_PATH + '../EAZY/OUTPUT/egs_candels.zout'
        FOUT_FILE = CAT_PATH+'egs_candels.fout'
        CAT_FILE = CAT_PATH + 'egs_candels.cat'
        KTOT_COL = 'ks'
        CAT_FILE = CAT_PATH + '../WIRDS_D3-95_Ks_ugrizJHKs_141927+524056_T0002.cat.candels'
        KTOT_COL = 'kstot'
        MAGS = True
        
        CAT_FILE = '/Users/brammer/3DHST/Photometry/Work/Test_v4.0/Catalogs/aegis_3dhst.v4.0.nzpcat.bright'
        ZOUT_FILE = '/Users/brammer/3DHST/Photometry/Work/Test_v4.0/ZPs/OUTPUT/aegis.zout'
        FOUT_FILE = '/Users/brammer/3DHST/Photometry/Work/Test_v4.0/ZPs/OUTPUT/aegis.fout'
        KTOT_COL = 'f_K'
        MAGS = False
        
    if goodsn:        
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('hyp'):
            CAT_PATH = '/3DHST/Photometry/Release/v4.1/GOODS-N/'
            CAT_FILE = CAT_PATH+'Catalog/goodsn_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH+'Eazy/goodsn_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH+'Fast/goodsn_3dhst.v4.1.fout'
            KTOT_COL = 'f_ks'
        
        #### Gabe's laptop
        if ('brammer' in unicorn.hostname().lower()) | ('gabriel' in unicorn.hostname().lower()):
            # CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.0/'
            # CAT_FILE = CAT_PATH + 'Eazy/EazyRerun/goodsn_3dhst.v4.0.cat'
            # ZOUT_FILE = CAT_PATH + 'Eazy/EazyRerun/OUTPUT/goodsn.dusty3.zout'
            # FOUT_FILE = CAT_PATH + 'Fast/goodsn_3dhst.v4.0.fout'
            # KTOT_COL = 'f_F160W'
        
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/goodsn_3dhst.v4.1.cats/'
            CAT_FILE = CAT_PATH + 'Catalog/goodsn_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/goodsn_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH + 'Fast/goodsn_3dhst.v4.1.fout'
            KTOT_COL = 'f_F160W'

            ### Use dusty-template fit
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun'
            ZOUT_FILE = CAT_PATH + '/OUTPUT/goodsn_3dhst.v4.1.dusty.zout'
            
    # #### Force UDF to use GOODS-S
    # if udf:
    #     cdfs = True
    #     udf = False
          
    # if udf and (unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp')):
    #     GRISM_PATH=unicorn.GRISM_HOME+'UDF/'
    #     CAT_PATH = '/3DHST/Ancillary/GOODS-S/HUDF09/CAT/'
    # 
    #     CAT_FILE = CAT_PATH+'hudf09.cat'
    #     ZOUT_FILE = CAT_PATH+'hudf09_noU.zout'
    #     FOUT_FILE = CAT_PATH + 'hudf09.fout'
    #     KTOT_COL = 'f_kv2_tot'
    # 
    # if udf and (unicorn.hostname().startswith('hyp')):
    #     GRISM_PATH=unicorn.GRISM_HOME+'HUDF/'
    #     CAT_PATH = '/3DHST/Spectra/Work/HUDF/INTERLACE/phot/'
    #     print CAT_PATH
    # 
    #     CAT_FILE = CAT_PATH+'hudfcat_mar3012.cat'
    #     ZOUT_FILE = CAT_PATH+'hudfcat_mar3012.zout'
    #     FOUT_FILE = CAT_PATH + 'hudfcat_mar3012.fout'
    #     KTOT_COL = 'f_kv2_tot'
    
    if udf | cdfs:
        
        CAT_PATH = '/3DHST/Photometry/Release/v4.1/GOODS-S/'
        CAT_FILE = CAT_PATH + 'Catalog/goodss_3dhst.v4.1.cat'
        ZOUT_FILE = CAT_PATH + 'Eazy/goodss_3dhst.v4.1.zout'
        FOUT_FILE = CAT_PATH + 'Fast/goodss_3dhst.v4.1.fout'
        KTOT_COL = 'f_ks'
        
        #### Gabe's laptop
        if ('brammer' in unicorn.hostname().lower()) | ('gabriel' in unicorn.hostname().lower()):
            # CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.0/'
            # CAT_FILE = CAT_PATH + 'Eazy/EazyRerun/goodss_3dhst.v4.0.cat'
            # ZOUT_FILE = CAT_PATH + 'Eazy/EazyRerun/OUTPUT/goodss.dusty3.zout'
            # FOUT_FILE = CAT_PATH + 'Fast/goodss_3dhst.v4.0.fout'
            # KTOT_COL = 'f_F160W'

            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/goodss_3dhst.v4.1.cats/'
            CAT_FILE = CAT_PATH + 'Catalog/goodss_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/goodss_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH + 'Fast/goodss_3dhst.v4.1.fout'
            KTOT_COL = 'f_F160W'
            
            ### Use dusty-template fit
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun'
            ZOUT_FILE = CAT_PATH + '/OUTPUT/goodss_3dhst.v4.1.dusty.zout'
            
    if uds:
        #for the v4.1 reduction:
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('hyp'):
            CAT_PATH = '/3DHST/Photometry/Release/v4.1/UDS/'
            CAT_FILE = CAT_PATH+'Catalog/uds_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH+'Eazy/uds_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH +'Fast/uds_3dhst.v4.1.fout'
            KTOT_COL = 'f_k'
                
        #### Gabe's laptop
        if ('brammer' in unicorn.hostname().lower()) | ('gabriel' in unicorn.hostname().lower()):
            # CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.0/'            
            # CAT_FILE = CAT_PATH + 'Eazy/EazyRerun/uds_3dhst.v4.0.cat'
            # ZOUT_FILE = CAT_PATH + 'Eazy/EazyRerun/OUTPUT/uds.dusty3.zout'
            # FOUT_FILE = CAT_PATH + 'Fast/uds_3dhst.v4.0.fout'
            # KTOT_COL = 'f_F160W'
            
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/uds_3dhst.v4.1.cats/'
            CAT_FILE = CAT_PATH + 'Catalog/uds_3dhst.v4.1.cat'
            ZOUT_FILE = CAT_PATH + 'Eazy/uds_3dhst.v4.1.zout'
            FOUT_FILE = CAT_PATH + 'Fast/uds_3dhst.v4.1.fout'
            KTOT_COL = 'f_F160W'

            ### Use dusty-template fit
            CAT_PATH = '/Users/brammer/3DHST/Spectra/Release/v4.1/EazyRerun'
            ZOUT_FILE = CAT_PATH + '/OUTPUT/uds_3dhst.v4.1.dusty.zout'
            
    #### UltraVISTA in COSMOS
    if 'UltraV' in root:
        CAT_PATH = '/Users/brammer/3DHST/Ancillary/COSMOS/UltraVISTA/'
        CAT_FILE = CAT_PATH + 'UVISTA_final_v4.1_head.cat'
        ZOUT_FILE = CAT_PATH + 'UVISTA_final_v4.1.zout'
        ZOUT_FILE = CAT_PATH + 'OUTPUT/uvista.full.dusty3.zout'
        FOUT_FILE = CAT_PATH+'UVISTA_final_BC03_v4.1.fout'
        KTOT_COL = 'Ks_tot'
                
    #### D. Erb's grism program
    if erb:
        GRISM_PATH=unicorn.GRISM_HOME+'Erb/'
        # CAT_PATH = GRISM_PATH+'Catalogs/Combined/EAZY/FAST/'
        # CAT_FILE = CAT_PATH+'erb.cat'
        # ZOUT_FILE = CAT_PATH+'erb.zout'
        # FOUT_FILE = CAT_PATH + 'erb.fout'
        CAT_PATH = GRISM_PATH+'Catalogs/Combined/EAZY/'
        CAT_FILE = CAT_PATH+'reddy_ugrJH.cat'
        ZOUT_FILE = CAT_PATH+'OUTPUT/reddy.zout'
        FOUT_FILE = CAT_PATH + 'reddy.fout'
        
        # CAT_PATH = GRISM_PATH+'Catalogs/Combined/EAZY/'
        # CAT_FILE = CAT_PATH+'erb2.cat'
        # ZOUT_FILE = CAT_PATH+'OUTPUT/erb2.zout'
        # FOUT_FILE = CAT_PATH + 'erb2.fout'
        KTOT_COL = 'f_h'
    
    # #### Frontier fields
    # if ('2744' in root.split('_')[0]) & root.upper().startswith('A'):
    #     GRISM_PATH = '/Users/brammer/Research/HST/FrontierFields/HLSP/Catalog/'
    #     CAT_PATH = '/Users/brammer/Research/HST/FrontierFields/HLSP/Catalog/'
    #     CAT_FILE = CAT_PATH + 'a2744_Ks_v0.1.cat'
    #     ZOUT_FILE = CAT_PATH + 'OUTPUT/a2744.zout'
    #     FOUT_FILE = CAT_PATH + 'OUTPUT/a2744.fout'
    #     KTOT_COL = 'f_ks_iso'
    # 
    #     CAT_PATH = '/Users/brammer/Research/HST/FrontierFields/Mosaic/CombinedCatalog/EAZY/'
    #     CAT_FILE = CAT_PATH + 'ff_Kscl.cat'
    #     ZOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_withK.zout'
    #     FOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_withK.fout'
    #     KTOT_COL = 'f_Ks'
        
            
    #### Supernova cluster pointing
    if ('MACSJ1720' in root) | ('PANCHA' in root):
        CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/Perlmutter/Eazy/'
        #CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_macs1720_cat.txt.flux_AB25'
        #ZOUT_FILE = CAT_PATH+'OUTPUT/m1720.zout'
        #FOUT_FILE = CAT_PATH+'OUTPUT/m1720.fout'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_ir_macs1720_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m1720ir.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m1720ir.fout'
        KTOT_COL = 'f160w_flux'

    #### CLASH / GLASS / FF
    if ('2744' in root.split('_')[0]) & root.upper().startswith('A'):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'abell2744.cat.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/a2744.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/a2744.fout'
        KTOT_COL = 'f160w_flux'
        
        CAT_PATH = '/Users/brammer/Research/HST/FrontierFields/Mosaic/CombinedCatalog/EAZY/'
        CAT_FILE = CAT_PATH + 'ff_Kscl.cat'
        hff='with'
        #hff='no'
        ZOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_%sK.zout' %(hff)
        FOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_%sK.fout' %(hff)
        KTOT_COL = 'f_Ks'
        
    if ('0717' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_macs0717_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m0717.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m0717.fout'
        KTOT_COL = 'f160w_flux'
    
    if ('0416' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_macs0416_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m0416.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m0416.fout'
        KTOT_COL = 'f160w_flux'
        
        CAT_PATH = '/Users/brammer/Research/HST/FrontierFields/Mosaic/CombinedCatalog/EAZY/'
        CAT_FILE = CAT_PATH + 'ff_Kscl.cat'
        hff='with'
        #hff='no'
        ZOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_%sK.zout' %(hff)
        FOUT_FILE = CAT_PATH + 'OUTPUT/hawkiff_%sK.fout' %(hff)
        KTOT_COL = 'f_Ks'

        CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/GLASS/MACS0416/Catalog/'
        CAT_FILE = CAT_PATH + 'macs0416_HST_vX.cat'
        ZOUT_FILE = CAT_PATH + 'EAZY/OUTPUT/m0416_test.zout'
        FOUT_FILE = CAT_PATH + 'EAZY/OUTPUT/m0416_test.fout'
        KTOT_COL = 'flux_f160w'
        
    if ('1423' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp-clash_hst_acs-ir_macs1423_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m1423.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m1423.fout'
        KTOT_COL = 'f160w_flux'
    
    if ('1347' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_rxj1347_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/r1347.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/r1347.fout'
        KTOT_COL = 'f160w_flux'
    #
    if ('2248' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        #CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_rxj2248_cat.txt.flux_AB25'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_ir_rxj2248_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/rxj2248.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/rxj2248.fout'
        KTOT_COL = 'f160w_flux'
     
    if ('1149' in root.split('_')[0]) | ('REFSDAL' in root.upper()):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_macs1149_cat.txt.flux_hasH' #'.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m1149.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m1149.fout'
        KTOT_COL = 'f160w_flux'
        # 
        # CAT_PATH = '/Users/brammer/Research/HST/FrontierFields/HLSP/M1149_cluster/'
        # CAT_FILE = CAT_PATH+'/macs1149_HST_vX.cat' #'.flux_AB25'
        # ZOUT_FILE = CAT_PATH+'EAZY/OUTPUT/macs1149_ff_noTE.zout'
        # FOUT_FILE = CAT_PATH+'EAZY/OUTPUT/macs1149_ff.fout'
        # KTOT_COL = 'flux_f160w'
        
    #
    if ('2129' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_ir_macs2129_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m2129.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m2129.fout'
        KTOT_COL = 'f160w_flux'
    #
    if ('0647' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_ir_macs0647_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m0647.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m0647.fout'
        KTOT_COL = 'f160w_flux'
    #
    if ('0744' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/Research/HST/CLASH/Eazy/'
        #hlsp_clash_hst_acs-ir_macs0744_cat.txt.flux_AB25
        CAT_FILE = CAT_PATH+'hlsp_clash_hst_acs-ir_macs0744_cat.txt.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/m0744.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/m0744.fout'
        KTOT_COL = 'f160w_flux'
    #
    if ('0327' in root.split('_')[0]):
        CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/Rigby/EAZY/'
        CAT_FILE = CAT_PATH+'RCS0327-2015.225-HST.cat' #'.flux_AB25'
        ZOUT_FILE = CAT_PATH+'OUTPUT/rcs0327.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/rcs0327.fout'
        KTOT_COL = 'flux_F160W'
    
    #### D. Newman's grism program
    if ('JKCS041' in root):
        CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/Newman/WIRDS/'
        CAT_FILE = CAT_PATH+'WIRDS_D1_KSN7.cat'
        ZOUT_FILE = CAT_PATH+'OUTPUT/wirds.zout'
        FOUT_FILE = CAT_PATH+'WIRDS_D1_KSN7_sci.cat'
        KTOT_COL = 'f_ks'
    
    #### Muzzin SpARCS clusters
    if 'SPARCS' in root.upper():
        # cluster = '-'.join(root.split('-')[:-1])
        # CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/Muzzin/Catalogs/'
        # CAT_FILE = CAT_PATH+'%s_v2.0.cat' %(cluster)
        # ZOUT_FILE = CAT_PATH+'OUTPUT/%s_v2.0.zout' %(cluster)
        # FOUT_FILE = CAT_PATH+'OUTPUT/%s_v2.0.fout' %(cluster)
        # KTOT_COL = 'K'
        cluster = '-'.join(root.split('-')[:2]).replace('SPARCS','SpARCS')
        CAT_PATH = '/Users/brammer/3DHST/Spectra/Work/Muzzin/Catalogs/RvdB/'
        CAT_FILE = CAT_PATH+'%s_total_F140Wsel.cat' %(cluster)
        ZOUT_FILE = CAT_PATH+'EAZY/%s_IRAC_HST/OUTPUT/photz.zout' %(cluster)
        FOUT_FILE = CAT_PATH+'EAZY/%s_IRAC_HST/OUTPUT/photz.fout' %(cluster)
        KTOT_COL = 'K'
        
    if ecdfs | root.startswith('ECDFS'):
        CAT_PATH = '/Users/brammer/3DHST/Ancillary/GOODS-S/Cardamone/'
        CAT_FILE = CAT_PATH+'ECDFS_BVRdet_Subaru_v1_hasz.cat'
        ZOUT_FILE = CAT_PATH+'OUTPUT/cardamone_full.zout'
        FOUT_FILE = CAT_PATH+'OUTPUT/cardamone_full.fout'
        KTOT_COL = 'f_K'

    if KTOT_COL is None:
        """
        All catalog flags are False
        """
        return None, None, None
    
    print CAT_FILE
    
    #### Read the catalogs
    #cat = catIO.ReadASCIICat(CAT_FILE)
    cat = catIO.Readfile(CAT_FILE, force_lowercase=False); cat.names = cat.columns
    #print KTOT_COL
    if 'field' in cat.names:
        ktot_field = cat[KTOT_COL.lower()]
    else:
        ktot_field = cat.field(KTOT_COL)
        
    if MAGS:
        cat.kmag = ktot_field
    else:
        cat.kmag = 25.-2.5*np.log10(ktot_field)
    
    cat.MAGS_UNIT = MAGS
    
    if 'star_flag' in cat.names:
        cat.star_flag = cat.field('star_flag')
    else:
        cat.star_flag = cat.kmag*0
    
    #### Different column names
    #if cdfs:
    #    cat.ra = cat.field('RA')
    #    cat.dec = cat.field('DEC')
    #    cat.id = np.cast[int](cat.field('ID'))
    
    cat.filename=CAT_FILE
    
    if ZOUT_FILE:
        zout = catIO.Readfile(ZOUT_FILE, force_lowercase=False)
        zout.filename=ZOUT_FILE
    else:
        zout = None

    if FOUT_FILE:
        fout = catIO.Readfile(FOUT_FILE, force_lowercase=False)
        fout.filename=FOUT_FILE
    else:
        fout = None
    
    return cat, zout, fout

def read_grism_files(root='COSMOS-3-G141', BASE_PATH='', GRISM_NAME='G141'):
    """
    Read root+'_drz.cat' and the associated SPC file.
    """
    import threedhst
    import unicorn.analysis
    
    grismCat, SPC = None, None
    
    if not BASE_PATH:
        BASE_PATH = get_grism_path(root)
        
    ##### catalog
    grismCat = threedhst.sex.mySexCat(BASE_PATH+'DATA/'+root+'_drz.cat')
    for col in grismCat.column_names:
        if col.startswith('MAG_F'):
            grismCat.MAG = grismCat[col]
            grismCat.DETECT_FILTER = col
            break
            
    ##### SPC file        
    if root+'_2_opt.SPC.fits' == unicorn.analysis.SPC_FILENAME:
        SPC = unicorn.analysis.SPC
    else:
        try:
            try:
                unicorn.analysis.SPC.fits.close()
            except:
                pass
            SPC = threedhst.plotting.SPCFile(root+'_2_opt.SPC.fits',
                axe_drizzle_dir=BASE_PATH+'DRIZZLE_'+GRISM_NAME)
            unicorn.analysis.SPC = SPC
            unicorn.analysis.SPC_FILENAME = SPC.filename
        except:
            SPC = None
            unicorn.analysis.SPC_FILENAME = None
            unicorn.analysis.SPC = None
    
    return grismCat, SPC
    
def make_SED_plots(grism_root='COSMOS-3-G141'):
    import unicorn.analysis
    
    PATH = unicorn.analysis.get_grism_path(grism_root)
    print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root)
    
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root)
        
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
    
def convolveWithThumb(id, lambdaz, temp_sed, SPC, oned=True, xint=None, verbose=False):
    """ 
    
    Convolve the best-fit eazy spectrum with the object shape
    
    """
    from scipy import convolve as conv
    
    thumb_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_thumb.fits.gz'
    try:
        thumb = pyfits.open(thumb_file)
    except:
        return lambdaz, temp_sed
    
    size = thumb[0].data.shape
    DLAM = np.sqrt(thumb[0].header['CD1_1']**2+thumb[0].header['CD1_2']**2)*3600./0.128254*46.5
    
    twod_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_2D.fits.gz'
    twod = pyfits.open(twod_file)
    model1D = np.matrix(twod[5].data.sum(axis=1))
    model2D = np.array(np.dot(np.transpose(model1D),np.ones((1,size[0]))))
    
    if oned:
        profile = np.sum(thumb[0].data*model2D,axis=0)
        profile /= profile.sum()
        if verbose: 
            xprof = np.arange(len(profile))
            print 'Profile center:', np.trapz(xprof, xprof*profile)/np.trapz(xprof, profile), np.mean(xprof)
    else:
        profile = thumb[0].data*model2D
        #for i in range(size[0]):
        #    profile[i,:] /= profile[i,:].sum()
        profile /= profile.sum()
    
    profile = profile[::-1]
    
    LSM = size[0]/2
    xmin = 3000
    xmax = 2.4e4
    
    q = np.where((lambdaz > (xmin-LSM)) & (lambdaz < (xmax+LSM)))[0]
    lambdaz = lambdaz[q]
    temp_sed = temp_sed[q]
    
    ### convolve with some factor of the object size
    #LSM = np.sqrt(LSM**2+(0.1*R[grism_idx]/0.128*46.5)**2)
    
    temp_sed_sm = temp_sed*1.
    # for i in range(len(lambdaz)):
    #     temp_sed_sm[i] = np.trapz(1./np.sqrt(2*np.pi*50**2)*np.exp(-0.5*(lambdaz-lambdaz[i])**2/50**2)*temp_sed, lambdaz)
    
    if verbose:
        print 'DLAM: ', DLAM
        
    #
    if xint is None:
        xint = np.arange(xmin-LSM*DLAM,xmax+LSM*DLAM,DLAM)
    
    yint = np.interp(xint, lambdaz, temp_sed_sm)
    
    # #### Original units
    yint_x = temp_sed/np.sum(temp_sed)
    dl = lambdaz[1]-lambdaz[0]
    #### Convolve with a gaussian
    xgauss = np.arange(10*46/dl)*dl-5*46
    ygauss = np.exp(-1*xgauss**2/2/((35/dl)**2))
    ygauss /= np.sum(ygauss)
    yintc = conv(yint_x, ygauss, mode='same')
    
    xprof = np.arange(0,len(profile))*21.75
    xprof_int = np.arange(0,len(profile)*21.75,dl)
    prof_int = np.interp(xprof_int, xprof, profile)
    prof_int /= np.sum(prof_int)
    full_profile = conv(yintc, prof_int, mode='same')
    
    if len(full_profile) == len(xgauss):
        x_full = xgauss
    else:
        x_full = xprof_int
        
    x_full -= x_full.mean()
    
    plt.rcParams['font.size'] = 9
    
    fig = Figure(figsize=[4,3], dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(xgauss, ygauss, color='blue')
    ax.plot(x_full, full_profile, color='red')
    ax.plot(xprof_int-xprof_int.mean(), prof_int, color='orange')
    ax.set_title('%s-G141_%05d' %(SPC.filename.split('-G141')[0], id))
    ax.set_xlabel(r'$\Delta\lambda$')
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure('/tmp/profile.png', dpi=100, transparent=False)
    
    x0 = lambdaz[temp_sed == temp_sed.max()][0]
    return x_full+x0, full_profile
    
    #### Idea: convolve to the size of the gaussian, don't need whole spectrum!
    
    #### Convolve with a gaussian
    xgauss = np.arange(20)*DLAM-10*DLAM
    ygauss = np.exp(-1*xgauss**2/2/35**2)
    ygauss /= np.sum(ygauss)
    yint = conv(yint, ygauss, mode='same')
    
    if oned:
        temp_sed_conv = conv(yint, profile, mode='same')
        # temp_sed_conv = yint*0.
        # for i in range(LSM,len(xint)-LSM):
        #     temp_sed_conv[i] = np.sum(yint[i-LSM:i+LSM]*profile)
    else:
        NX, NY = len(yint), size[0]
        temp_sed_conv = np.zeros((NY,NX))
        for i in range(size[0]):
            temp_sed_conv[i,:] = conv(yint, profile[i,:].flatten(), mode='same') #np.dot(yint[i-LSM:i+LSM],profile).reshape((NY,))
            
    
    thumb.close()
    
    return xint, temp_sed_conv
    
def specphot(id=69, grism_root='ibhm45030',
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
    CACHE_FILE = 'Same', Verbose=False,
    SPC = None, cat=None, grismCat = None,
    zout = None, fout = None, OUT_PATH='/tmp/', OUT_FILE_FORMAT=True,
    OUT_FILE='junk.png', GET_SPEC_ONLY=False):
    """
specphot(id)
    
    Get photometry/SED fit and compare G141 spectrum
    """
    #import scipy.interpolate as interpol
       
    ### 69, 54!
    
    xxx = """
    id=199
    grism_root='ibhm48'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    CACHE_FILE = 'Same'
    """
    
    #### Get G141 spectrum
    if Verbose:
        print 'Read SPC'
    
    if SPC is None:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                    axe_drizzle_dir='DRIZZLE_G141')
                    
    spec = SPC.getSpec(id)
    if spec is False:
        return False
        
    xmin = 3000
    xmax = 2.4e4
    
    lam = spec.field('LAMBDA')
    flux = spec.field('FLUX')
    ffix = flux-spec.field('CONTAM')
    ferr = spec.field('FERROR') #*0.06/0.128254
        
    if Verbose:
        print 'Read grism catalog'
        
    #### Read the grism catalog and get coords of desired object
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    #### Source size
    R = np.sqrt(np.cast[float](grismCat.A_IMAGE)*np.cast[float](grismCat.B_IMAGE))
    grism_idx = np.where(grismCat.id == id)[0][0]
    
    Rmatch = R[grism_idx]*1.
    
    #print 'R=%f"' %(Rmatch)
    ra0 = grismCat.ra[grismCat.id == id][0]
    de0 = grismCat.dec[grismCat.id == id][0]
    
    #### Read EAZY outputs and get info for desired object
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
    
    
    photom_idx = np.where(dr == np.min(dr))[0][0]
    
    drMatch = dr[photom_idx]*1.
    #print 'dr = %7.2f\n' %(drMatch)
    #print drMatch, np.min(dr)
    
    if drMatch > 2:
        return False
        
    if Verbose:
        print 'Read zout'
    if zout is None:    
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
        
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    if Verbose:
        print 'Read binaries'
        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(photom_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE)
         
    try:
        lambdaz, temp_sed_sm = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
    except:
        temp_sed_sm = temp_sed*1.
        
    if Verbose: 
        print 'Normalize spectrum'
        
    #### Normalize G141 spectrum
    #interp = interpol.interp1d(lambdaz, temp_sed_sm, kind='linear')

    q = np.where((lam > 1.08e4) & (lam < 1.68e4) & (flux > 0))[0]
    #### G102
    if lam.min() < 9000:
        q = np.where((lam > 0.8e4) & (lam < 1.13e4) & (flux > 0))[0]
    
    #### ACS G800L
    if lam.min() < 5000:
        q = np.where((lam > 0.55e4) & (lam < 1.0e4) & (flux > 0))[0]
        
    if len(q) == 0:
        return False

    yint = np.interp(lam[q], lambdaz, temp_sed_sm)
        
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    if np.isnan(anorm):
        anorm = 1.0
    total_err = np.sqrt((ferr)**2+(1.0*spec.field('CONTAM'))**2)*anorm
    
    if GET_SPEC_ONLY:
        if drMatch > 1:
            return False
        else:
            return lam, ffix*anorm, total_err, lci, fobs, efobs, photom_idx
         
    if Verbose:
        print 'Start plot'
        
    #### Make the plot
    threedhst.plotting.defaultPlotParameters()
    
    xs=5.8
    ys = xs/4.8*3.2
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[xs,ys],dpi=100)
    else:
        fig = Figure(figsize=[xs,ys], dpi=100)
    
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.13*4.8/xs, bottom=0.15*4.8/xs,right=1.-0.02*4.8/xs,top=1-0.10*4.8/xs)
    
    ax = fig.add_subplot(111)
    
    ymax = np.max((ffix[q])*anorm)
    
    if Verbose:
        print 'Make the plot'
        
    ax.plot(lambdaz, temp_sed_sm, color='red')
    # plt.errorbar(lam[q], ffix[q]*anorm, yerr=ferr[q]*anorm, color='blue', alpha=0.8)
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.2, linewidth=1)
    
    #### Show own extraction
    sp1d = threedhst.spec1d.extract1D(id, root=grism_root, path='./HTML', show=False, out2d=False)
    lam = sp1d['lam']
    flux = sp1d['flux']
    ffix = sp1d['flux']-sp1d['contam'] #-sp1d['background']
    ferr = sp1d['error']
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    ax.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.6, linewidth=1)
    
    #### Show photometry + eazy template
    ax.errorbar(lci, fobs, yerr=efobs, color='orange', marker='o', markersize=10, linestyle='None', alpha=0.4)
    ax.plot(lambdaz, temp_sed_sm, color='red', alpha=0.4)

    ax.set_ylabel(r'$f_{\lambda}$')
    
    if plt.rcParams['text.usetex']:
        ax.set_xlabel(r'$\lambda$ [\AA]')
        ax.set_title('%s: \#%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
    else:
        ax.set_xlabel(r'$\lambda$ [$\AA$]')
        ax.set_title('%s: #%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[photom_idx]))
        
    #kmag = 25-2.5*np.log10(cat.ktot[photom_idx])
    kmag = cat.kmag[photom_idx]
    
    ##### Labels
    label = 'ID='+r'%s   K=%4.1f  $\log M$=%4.1f' %(np.int(cat.id[photom_idx]),
        kmag, fout.field('lmass')[photom_idx])
        
    ax.text(5e3,1.08*ymax, label, horizontalalignment='left',
      verticalalignment='bottom')
    
    
    label = 'R=%4.1f"' %(drMatch)
    if drMatch > 1.1:
        label_color = 'red'
    else:
        label_color = 'black'
    ax.text(2.2e4,1.08*ymax, label, horizontalalignment='right',
      color=label_color, verticalalignment='bottom')
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.1*ymax,1.2*ymax)
    
    if Verbose:
        print 'Save the plot'
    
    if OUT_FILE_FORMAT:
        out_file = '%s_%05d_SED.png' %(grism_root, id)
    else:
        out_file = OUT_FILE
        
    if USE_PLOT_GUI:
        fig.savefig(OUT_PATH+'/'+out_file,dpi=100,transparent=False)
        plt.close()
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(OUT_PATH+'/'+out_file, dpi=100, transparent=False)
    
    print unicorn.noNewLine+OUT_PATH+'/'+out_file
    
    if Verbose:
        print 'Close the plot window'
        
    
def match_grism_to_phot(grism_root='ibhm45', MAIN_OUTPUT_FILE = 'cosmos-1.v4.6', OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/', CACHE_FILE = 'Same', Verbose=False, SPC = None, cat=None, grismCat = None, zout = None, fout = None, OUTPUT='/tmp/match.cat'):

    """
    Generate the matched catalogs of 3D-HST objects and objects in ancillary catalogs
    """
    import unicorn.analysis
    # 
    # MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    # OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    # CACHE_FILE = 'Same'
    
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    if zout is None:
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    #rfUV = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.153-155.rf')
    
    z_peak = zout.field('z_peak')
    # kmag = 25-2.5*np.log10(cat.field('ktot'))
    kmag = cat.kmag
    
    #uv = -2.5*np.log10(rfUV.field('L153')/rfUV.field('L155'))
    
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    drs = grismCat.ra*1.
    ids = np.cast[int](grismCat.ra*1.)
    
    cosfact = np.cos(np.median(grismCat.dec)/360*2*np.pi)
    
    for i in range(len(drs)):
        dr = np.sqrt((cat.ra-grismCat.ra[i])**2*cosfact**2+
                     (cat.dec-grismCat.dec[i])**2)*3600
        mat = np.where(dr == np.min(dr))[0][0]
        drs[i] = dr[mat]
        ids[i] = mat
    
    fp = open(OUTPUT,'w')
    fp.write("# id_f140w  mag_f140w  id_phot  mag_Ktot  Rmatch  z_peak  logM star_flag fcontam\n# Rmatch = match distance in arcsec (should be < 1)\n# %s\n" %(MAIN_OUTPUT_FILE))
    for i in range(len(drs)):
        j = ids[i]
        line = "%5d %6.2f %8d %6.2f %6.2f %6.2f %6.2f %d %5s\n" %(grismCat.id[i],
                  np.float(grismCat.MAG[i]), np.int(cat.id[j]),
                  kmag[j], drs[i], 
                  z_peak[j], fout.field('lmass')[j], 
                  np.int(cat.star_flag[j]), grismCat['FCONTAM'][i])
        if grismCat.id[i] in SPC._ext_map:
            fp.write(line)
    fp.close()

def make_multiple_fluximage(grism_root='COSMOS-3-G141'):
    import unicorn.analysis
    waves = [1.1e4,1.25e4,1.6e4]
    for wave in waves:
        unicorn.analysis.make_fluximage(grism_root=grism_root, wavelength=wave)

def make_fluximage(grism_root='COSMOS-3-G141', wavelength=1.1e4, direct_image=None, match_toler=1, verbose=True):
    """
    1) Read a SExtractor 3DHST catalog
    2) Open the corresponding direct image
    3) Match each object in the 3DHST catalog to an external catalog, which
       has redshifts / SED fits
    4) Get the ('wavelength' - detection band) color for each object and 
       scale its segmentation image accordingly
    5) If not match within 'match_toler', just use the detected flux
    6) Write out the new scaled image
    
    (2011-10-11, this is old code where I tried to make a fake fluxcube with observed
    colors from photometric catalogs.  Never really used and deprecated with the
    availability of CANDELS imaging.)
    """
    
    from pyraf import iraf
    
    out_image = 'DATA/'+grism_root.replace('G141','f%03d' %(wavelength/100))+'.fits'
    
    ##### Get the path and read the catalogs
    PATH = unicorn.analysis.get_grism_path(grism_root)
    PWD=os.getcwd()
    print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=grism_root)
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    ## read grism outputs
    grismCat, SPC = unicorn.analysis.read_grism_files(root=grism_root)
    
    detect_wlen = np.float(grismCat.DETECT_FILTER.strip('MAG_FW'))*10
    
    #### Detection and segmentation images
    if direct_image is None:
        direct_image = glob.glob('PREP_FLT/'+grism_root.replace('G141','*')+'_drz.fits')[0]
    
    seg_file = grismCat.filename.replace('drz.cat','seg.fits')
    seg_file = threedhst.utils.find_fits_gz(seg_file)
    
    direct = pyfits.open(direct_image)
    seg = pyfits.open(seg_file)
    
    #### Loop through objects in the catalog
    cosfact = np.cos(np.median(grismCat.dec)/360*2*np.pi)
    xint = np.array([wavelength, detect_wlen])

    #### If a Grism SPC file exists, only use IDs defined there
    #### Otherwise use all objects in the SEx. catalog
    if SPC is not None:
        ids = SPC._ext_map
    else:
        ids = grismCat.id
        
    for j, id in enumerate(ids):
        progress = '%2d' %(np.int(j*100./len(ids))) + '%'
        print unicorn.noNewLine+out_image+':  '+progress
            
        i = np.arange(grismCat.nrows)[grismCat.id == id][0]
           
        dr = np.sqrt((cat.ra-grismCat.ra[i])**2*cosfact**2+
                         (cat.dec-grismCat.dec[i])**2)*3600
        mat = np.where(dr == np.min(dr))[0][0]
        
        scale = 1.
        
        if dr[mat] < match_toler:
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
                eazy.getEazySED(mat, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                     OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                     CACHE_FILE = 'Same')
                        
            #### Smooth to grism resolution
            try:
                lambdaz, temp_sed = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
            except:
                pass
                # lambdaz2, temp_sed2 = unicorn.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC, oned=False)
            
            yint = np.interp(xint, lambdaz, temp_sed)
            scale = yint[0]/yint[1]    
            direct[1].data[seg[0].data == id] *= scale
    
    ### f_nu
    direct[1].data *= wavelength**2/detect_wlen**2
    
    ### Write the image, but keep just the SCI extension
    direct[1].writeto('/tmp/fluximage.fits', clobber=True)
    try:
        os.remove(out_image)
    except:
        pass
    
    os.chdir(PWD)
    iraf.imcopy('/tmp/fluximage.fits[1]',out_image)
    
    return out_image
        
        
def show_massive_galaxies(masslim=10.5, maglim=23.5, zrange=(0,5), 
    use_kmag=False, contam=0.5, coverage=0.9, skip_goodsn=False):        
    """
    Make a webpage showing objects selected on mass, redshift, contamination, mag...
    """
    
    if unicorn.hostname().startswith('unicorn'):
        os.chdir('/Library/WebServer/Documents/P/GRISM_v1.7/ANALYSIS')
        scripts="../scripts"
        files = glob.glob('../SED/*match.cat')
        matches = []
        for file in files:
            if skip_goodsn:
                if not 'GOODS-N' in file:
                    matches.append(file)
            else:
                matches.append(file)
    else:
        os.chdir(unicorn.GRISM_HOME+'ANALYSIS')
        scripts="http://localhost/~gbrammer/COSMOS/scripts"

        matches = []
        matches.extend(glob.glob('../AEGIS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../COSMOS/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../GOODS-S/HTML/SED/*match.cat'))
        if not skip_goodsn:
            matches.extend(glob.glob('../GOODS-N/HTML/SED/*match.cat'))
        
        # matches.extend(glob.glob('../SN-GEORGE/HTML/SED/*match.cat'))
        # matches.extend(glob.glob('../SN-PRIMO/HTML/SED/*match.cat'))    
        # matches.extend(glob.glob('../SN-MARSHALL/HTML/SED/*match.cat'))
        matches.extend(glob.glob('../ERS/HTML/SED/*match.cat'))
    
    print matches
    
    fplist = open('massive.dat','w')
    fplist.write('# ID   ra   dec  z lmass\n')
    fplist.write('# mag > %.1f, mass > %.1f, z=(%.1f,%.1f), contam=%.2f\n' %(maglim, masslim, zrange[0], zrange[1], contam))
    
    fpreg = open('massive.reg', 'w')
    fpreg.write('fk5\n')
    
    fp = open('massive.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="%s/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="%s/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="%s/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[2,2]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                        6: {
                                sorter: false
                        },
                        7: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <p>
        mag_F140W < %.2f
        <br>log M/Msun > %.2f
        <br>%.2f < z < %.2f
        <br> contam(1.4mum) < %.2f
        <br> coverage > %.2f
    </p>
    
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> Mag_WFC3 </th>
        <th> z </th>
        <th> logM </th>
        <th> Thumb </th>
        <th> 2D </th>
        <th> 1D </th>
        <th> SED </th>
    </thead>
    <tbody>
    """ %(scripts, scripts, scripts,
          maglim, masslim, zrange[0], zrange[1], contam, coverage))
    
    NUSE = 0
    for match in matches:
        root_path = os.path.dirname(match).split('SED')[0]
        root = os.path.basename(match).split('_match.cat')[0]
        
        c = catIO.Readfile(match)
        xml_file = root_path+root+'.xml'
        xml = catIO.readMarkerXML(xml_file)
        #print c.keys()
        
        #### SExtractor + contamination catalog
        sex = threedhst.sex.mySexCat(root_path+root+'_drz.cat')
        mat = c.id_f140w*0
        for i in range(len(c.id_f140w)):
            mat[i] = np.where(np.cast[int](sex.NUMBER) == c.id_f140w[i])[0][0]
        
        try:
            fcover = np.cast[float](sex.FCOVER)[mat]
        except:
            fcover = c.id_f140w*0.+1
            
        select_mag = c.mag_f140w
        if use_kmag:
            #print c.keys()
            select_mag = c.mag_ktot
        
        #### Had used wrong zeropoint in FAST    
        # if ('ERS' in root_path) | ('GOODS-S' in root_path):
        #     c.logm -= 0*-0.4*(23.86-25)
            
        use = (c.star_flag < 1) & (c.logm > masslim) & (c.rmatch < 1) & (select_mag < maglim) & (c.z_peak > zrange[0]) & (c.z_peak < zrange[1]) & (fcover > coverage)
        
        if 'fcontam' in c.keys():
            if contam < 0:
                use = use & (c.fcontam > -contam)
            else:
                use = use & (c.fcontam < contam)

        use = np.where(use == True)[0]
        use = use[np.argsort(c.logm[use])]
        NUSE += len(use)
        
        for i in use:
            ID = c.id_f140w[i]
            
            file="%s_%05d" %(root, ID)
            
            print "%4d %4d %s" %(ID, len(xml), file)
            
            fplist.write("%15s %14.6f %14.6f %8.3f %5.2f\n" %(file, xml[ID].ra, xml[ID].dec, c.z_peak[i], c.logm[i]))
            
            fpreg.write('circle(%14.6f,%14.6f,1") # text={%s} color=magenta  \n' %(xml[ID].ra, xml[ID].dec, file))
            
            fp.write("""
        <tr>
            <td> %s<br> %13.6f %13.6f <br> %4.1f </td>
            <td> %5.1f </td>
            <td> %5.2f </td>
            <td> %5.2f </td>
            <td> <img src=%s/images/%s_thumb.png height=180px> </td>
            <td> <img src=%s/images/%s_2D.png height=180px> </td>
            <td> <img src=%s/images/%s_1D.png height=180px> </td>
            <td> <img src=%s/SED/%s_SED.png height=180px> </td>
        </tr>
        """ %(file, xml[ID].ra, xml[ID].dec, select_mag[i],
              c.mag_f140w[i], c.z_peak[i], c.logm[i],
              root_path, file,
              root_path, file,
              root_path, file,
              root_path, file))
            
    fp.write("</tbody></table></body></html>")
    fp.close()
    fplist.close()
    fpreg.close()
    
    print 'N = %d' %NUSE
        
def field_area(ra_in, dec_in):
    """
    Given a list of ra / dec coordinates, compute the field area covered by those
    objects.
    """
    from shapely.geometry import MultiPoint, Point, MultiPolygon
    import numpy as np
    import matplotlib.pyplot as plt
    from descartes import PolygonPatch
    from shapely.ops import cascaded_union
    
    ra = ra_in*np.cos(np.mean(dec_in)/360.*2*np.pi)
    
    coords = []
    for i in xrange(len(ra)):
        coords.append((ra[i], dec_in[i]))
    
    env = MultiPoint(coords).envelope 
    env_area = env.area
    
    hull = MultiPoint(coords).convex_hull 
    hull_area = hull.area
    
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    
    ax.plot(ra, dec_in, marker='.', linestyle='None', alpha=0.2, color='blue')
    x_e, y_e = env.boundary.xy
    ax.plot(x_e, y_e, color='red')

    x_h, y_h = hull.boundary.xy
    ax.plot(x_h, y_h, color='yellow')
        
    # plt.show()
    
    return hull_area


def make_full_selection(zmin=None, zmax=None):
    import unicorn
    
    ########## Full selections to get everything
    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=23.5, zrange=(0.2,3.5))
    out='full_faint.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=22, zrange=(0.2,3.5),  use_kmag=False, contam=0.5, coverage=0.8)
    out='full_bright.html'
    unicorn.analysis.run_eazy_products_on_html(out)

    unicorn.analysis.show_massive_galaxies(masslim=10.49, maglim=24, zrange=(0.2,3.5),  use_kmag=False, contam=0.5, coverage=0.8)
    out='full_massive.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    ########## Bright galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=21., zrange=(0.7,2.8),  use_kmag=False, contam=0.05, coverage=0.9)
    out='full_z_F140W_lt_21.html'

    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=22., zrange=(1.0,1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='full_z_F140W_lt_22.html'
    
    ########## Massive galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=22., zrange=(0.7,2.8),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_11_F140W_22.html'

    ########## Massive z>1 galaxies
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=24., zrange=(1, 1.5),  use_kmag=False, contam=0.2, coverage=0.9)
    out='mass_11_z_1_1.5.html'
    unicorn.analysis.run_eazy_products_on_html(out)

    unicorn.analysis.show_massive_galaxies(masslim=10.49, maglim=23., zrange=(1, 1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_10.5_z_1_1.5.html'

    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=23., zrange=(1.5, 2.5),  use_kmag=False, contam=0.2, coverage=0.9)
    out='mass_11_z_1.5_2.5.html'
    
    unicorn.analysis.show_massive_galaxies(masslim=10., maglim=22., zrange=(1.0,1.5),  use_kmag=False, contam=0.05, coverage=0.9)
    out='mass_10_mag_22_z_1.0_1.5.html'

    unicorn.analysis.show_massive_galaxies(masslim=8., maglim=21.5, zrange=(0.8,1.4),  use_kmag=False, contam=0.05, coverage=0.95, skip_goodsn=False)
    unicorn.analysis.show_massive_galaxies(masslim=10.99, maglim=22.0, zrange=(0.8,1.5),  use_kmag=False, contam=0.05, coverage=0.95, skip_goodsn=False)
    out='test.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
    unicorn.analysis.show_massive_galaxies(masslim=9, maglim=22, zrange=(1.1,4),  use_kmag=False, contam=0.1, coverage=0.8)
    out='test.html'
    unicorn.analysis.run_eazy_products_on_html(out)
    
def run_eazy_products_on_html(out):
    import unicorn
    import os
    import shutil
    
    if out is not 'massive.html':
        shutil.move('massive.html', out)
        
    shutil.copy('/Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS/'+out, unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    #out='full_faint.html'
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    
    unicorn.analysis.process_eazy_redshifts(html=out, zmin=0.2, zmax=3.8, compress=0.8)
    unicorn.analysis.show_triple_zphot_zspec(zout=out.replace('html','zout'), zmin=0, zmax=2.5)
    
    #status = os.system('rsync -avz *.png *.html *.pdf ~/Sites_GLOBAL/P/GRISM_v1.6/EAZY/')
    status = os.system('rsync -avz *zz.png *.html *.pdf ~/Sites_GLOBAL/P/GRISM_v1.6/EAZY/')
    
def process_eazy_redshifts(html='massive.html', zmin=None, zmax=None, compress=1.0):
    import unicorn
    
    fp = open(html)
    lines = fp.readlines()
    fp.close()
    
    # example_zout = glob.glob('OUTPUT/*.zout')[0]
    # os.system('head -1 ' + example_zout + ' > '+html.replace('html','zout'))
    fp = open(html.replace('html','zout'),'w')
    fp.write('# id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68  l95 u95  l99 u99  nfilt q_z z_peak peak_prob z_mc\n')
    fp.close()
    
    for i,line in enumerate(lines):
        failed = False
        
        if 'G141' in line and 'png' not in line:
            object = line.split('<br>')[0].split('<td>')[1].strip()
            root=object.split('G141_')[0]+'G141'
            id = int(object.split('G141_')[1])
            
            if not os.path.exists('OUTPUT/%s.zout' %(object)):
                try:
                    print 'MISSING! '+object
                    #unicorn.analysis.run_eazy_fit(root=root, id=id, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = '/usr/local/bin/eazy_latest', zmin=zmin, zmax=zmax, compress=compress, COMPUTE_TILT=True, TILT_ORDER=1)
                    #os.system('cat OUTPUT/threedhst.zout > OUTPUT/%s.zout' %(object))
                except:
                    status = os.system('echo "" > OUTPUT/%s.zout' %(object))
                    failed = True
                    pass
            
            if not failed:
                status = os.system('grep -v "#" OUTPUT/%s.zout >> %s' %(object, html.replace('html','zout')))
        
        if ('/SED/' in line) & (not failed):
            lines[i] = line.replace('..//SED/','./').replace('SED.png height=180','eazy.png height=250')
        
        if "<table id=" in line:
            insert_line = i
    
    #### photz-specz plot
    unicorn.analysis.show_triple_zphot_zspec(zout=html.replace('html','zout'), zmin=0, zmax=2.5)
    lines.insert(insert_line, '<a href=%s_zz.pdf><img src=%s_zz.png></a>\n' %(html.replace('.html',''), html.replace('.html','')))
    
    fp = open(html.replace('.html','_eazy.html'),'w')
    fp.writelines(lines)
    fp.close()
    
def show_triple_zphot_zspec(zout='massive.zout', zmin=0, zmax=4):
    import threedhst.catIO as catIO
    zo = catIO.Readfile(zout)
    dz = (zo.z_peak-zo.z_spec)/(1+zo.z_spec)
    dz = dz[zo.z_spec > 0]
    
    threedhst.utils.biweight(dz[0::3])
    
    ### setup
    plt.rcParams['font.size'] = 10
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[7,2.6],dpi=100)
    else:
        fig = Figure(figsize=[7,2.6],dpi=100)
    
    #
    #fig = plt.figure(figsize=[8, 2.6],dpi=100)
    
    fig.subplots_adjust(wspace=0.27,hspace=0.12,left=0.07,
                        bottom=0.15,right=0.98,top=0.98)
    
    #################################### Broad-band + spectrum
    ax = fig.add_subplot(131)
    
    point_size = 4
    
    ax.plot(zo.z_spec[1::3], zo.z_peak[1::3], marker='o', linestyle='None', alpha=0.5, color='orange', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[1::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.text(0.5, 0.9, r'$N=$'+'%0d' %(len(dz[1::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{phot}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)

    ax = fig.add_subplot(132)
    
    ax.plot(zo.z_spec[0::3], zo.z_peak[0::3], marker='o', linestyle='None', alpha=0.5, color='purple', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[0::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{phot+grism}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)
    
    ax = fig.add_subplot(133)
    
    ax.plot(zo.z_spec[2::3], zo.z_peak[2::3], marker='o', linestyle='None', alpha=0.5, color='blue', markersize=point_size)
    ax.plot([0,4],[0,4], color='black', alpha=0.1)
    ax.text(0.5, 0.05, r'$\sigma=$'+'%5.3f' %(threedhst.utils.biweight(dz[2::3])), transform = ax.transAxes, horizontalalignment='center')
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    ax.set_ylabel(r'$z_\mathrm{grism\ only}$')
    ax.set_xlim(zmin, zmax)
    ax.set_ylim(zmin, zmax)
    
    #
    outfile = zout.replace('.zout','_zz')
    if USE_PLOT_GUI:
        plt.savefig(outfile+'.png')
        plt.savefig(outfile+'.pdf')
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile+'.png', dpi=100, transparent=False)
        canvas.print_figure(outfile+'.pdf', dpi=100, transparent=False)
    
def eazy_lists():
    """
    Various sublists to re-run eazy
    """
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    #### Objects with spec-z
    maglim = 24
    qzmax = 0.2
    contam_max = 0.05
    
    keep = (phot.mag_f1392w[phot.idx] < maglim) & (phot.fcontam[phot.idx] < contam_max) & (zout.q_z[0::3] < qzmax) & (phot.fcover[phot.idx] > 0.9) & (mcat.logm[mcat.idx] > 0) & (mcat.rmatch[mcat.idx] < 0.5) & (zsp.zspec[zsp.mat_idx] > 0) & (zsp.dr < 1)
    keep = keep & (zout.q_z[0::3] != zout.q_z[2::3])
    
    #### Original templates, individual lines, no tilt
    unicorn.analysis.run_eazy_on_list(ids=zout.id[0::3][keep], TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', TILT_ORDER=1, SCALE_SPEC_ERROR=2)
    unicorn.catalogs.make_full_redshift_catalog()
    os.system('mv full_redshift.cat full_redshift_scaleSpecErr2_yesTilt.cat')
    
    unicorn.analysis.run_eazy_on_list(ids=zout.id[0::3][keep], TEMPLATES_FILE='templates/fixed_lines_suppl.spectra.param', TILT_ORDER=0, SCALE_SPEC_ERROR=2)
    unicorn.catalogs.make_full_redshift_catalog()
    os.system('mv full_redshift.cat full_redshift_fixed_scaleSpecErr2_noTilt.cat')
    
    ### single
    ii = keep & (zout.id[0::3] == 'GOODS-N-14-G141_00526')
    unicorn.analysis.run_eazy_on_list(ids=zout.id[0::3][ii], TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', TILT_ORDER=0, SCALE_SPEC_ERROR=3.)
    
def run_eazy_on_list(ids = ['COSMOS-20-G141_01097'], compress=0.75, pipe=' > eazy.log', COMPUTE_TILT=True, TILT_ORDER=1, TEMPLATES_FILE='templates/fixed_lines_suppl.spectra.param', SCALE_SPEC_ERROR=1):
    """
    Run the eazy redshift code on a list of objects with id name like:
    [POINTING]-G141_[000ID].
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    for id in ids:
        pointing = id.split('G141')[0]+'G141'
        number = int(id.split('G141_')[1])
        result = unicorn.analysis.run_eazy_fit(root=pointing, id=number, compress=compress, zmin=0.1, zmax=4, TILT_ORDER=TILT_ORDER, pipe=pipe, force_zrange=True, COMPUTE_TILT=COMPUTE_TILT, TEMPLATES_FILE=TEMPLATES_FILE, zstep=0.025, SCALE_SPEC_ERROR=SCALE_SPEC_ERROR)
        
def run_eazy_on_all_objects(field='ERS', pipe=' > eazy.log', compress=0.75):
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')

    ######
    #logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/processed.log'

    logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/'+field+'.log'
    #field = 'ERS'

    # logfile = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/goods-s.log'
    # field = 'GOODS-S'
    
    ######
    
    if not os.path.exists(logfile):
        fp = open(logfile,'w')
        fp.write('')
        fp.close()
        
    fp = open(logfile)
    log_lines = fp.readlines()
    fp.close()
    
    os.chdir(unicorn.GRISM_HOME+field)
    catalogs = glob.glob('HTML/*drz.cat')
    print field, len(catalogs)
    
    for catalog in catalogs:
        os.chdir(unicorn.GRISM_HOME+field)
        if not os.path.exists(catalog.replace('HTML','HTML/SED').replace('drz','match')):
            continue
        else:
            print catalog
        cat = threedhst.sex.mySexCat(catalog)
        
        pointing = os.path.basename(catalog).split('G141')[0]+'G141'
        for id in cat.id[np.cast[float](cat.FCOVER) > 0.4]:
            object = '%s_%05d' %(pointing, id)
            if (object+'\n' not in log_lines) & (os.path.exists(unicorn.GRISM_HOME+field+'/HTML/ascii/'+object+'.dat')):
                result = True
                #try:
                result = unicorn.analysis.run_eazy_fit(root=pointing, id=id, compress=compress, zmin=0.1, zmax=4, TILT_ORDER=1, pipe=pipe, force_zrange=True, COMPUTE_TILT=True, TEMPLATES_FILE='templates/fixed_lines_suppl.spectra.param', zstep=0.025)
                #except:
                #    pass
                #    
                if result is False:
                    print object
                    root=pointing
                    status = os.system('rm %s_%05d' %(root, id) + '_threedhst.cat')
                    status = os.system('rm %s_%05d' %(root, id) + '.FILT.RES')
                    status = os.system('rm %s_%05d' %(root, id) + '.eazy.param')
                    status = os.system('rm templates/%s_%05d' %(root,id)+'.spectra.param')   
                    status = os.system('rm templates/%s_%05d' %(root,id)+'_spectrum.dat')
                #
                fp = open(logfile,'a')
                fp.write('%s\n' %(object))
                fp.close()
                #print unicorn.analysis.get_open_fds()

#
def get_open_fds():
    """
    Show open files because was getting "too many open files" error.
    """ 
    import fcntl, resource
    
    fds = []
    for fd in range(3,1024):
            try:
                    flags = fcntl.fcntl(fd, fcntl.F_GETFD)
            except IOError:
                    continue
			#
            fds.append(fd)
	#
    return fds

def generate_rf_color_inputs(FILTERS_RES='FILTER.RES.v8.R300', rf_filter_ids=range(135,164), filename='rf_dummy', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param'):
    """
    Make a fake catalog with all of the filters to interpolate. 
    Run with the line templates and Z_MAX = 0 to just get the template fluxes at z=0.
    """
    import threedhst.eazyPy as eazy
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
        eazy_binary = '/usr/local/bin/eazy_latest'
    else:
        eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy'
    
    #rf_filter_ids = range(135,164)
    hline = '# id '
    fluxes = ' 0 '
    for filt in rf_filter_ids:
        hline += ' F%0d E%0d' %(filt, filt)
        fluxes += ' 1.0 0.1'
    
    fp = open(filename+'.cat','w')
    fp.write(hline+'\n')
    fp.write(fluxes+'\n')
    fp.close()
    
    try:
        os.system('rm zphot.param.default')
    except:
        pass
    
    #### Generate default param file    
    os.system(eazy_binary + ' > /tmp/log')
    zp = eazy.EazyParam('zphot.param.default')
    
    zp.params['FILTERS_RES'] = FILTERS_RES
    zp.params['TEMPLATES_FILE'] = TEMPLATES_FILE
    zp.params['TEMPLATE_COMBOS'] = '1'
    zp.params['TEMP_ERR_FILE'] = 'templates/TEMPLATE_ERROR.eazy_v1.0'
    zp.params['WAVELENGTH_FILE'] = 'templates/EAZY_v1.1_lines/lambda_v1.1.def'
    zp.params['CATALOG_FILE'] = filename+'.cat'
    zp.params['Z_MIN'] = 0
    zp.params['Z_MAX'] = 4.0
    zp.params['Z_STEP'] = 0.1
    zp.params['Z_STEP_TYPE'] = 0   
    zp.params['MAIN_OUTPUT_FILE'] = filename
    zp.params['CACHE_FILE'] = filename+'.tempfilt'
    zp.params['APPLY_PRIOR'] = 'y'
    zp.params['PRIOR_FILTER'] = 153
    
    zp.write(filename+'.zphot.param')
    
    os.system(eazy_binary+' -p '+filename+'.zphot.param > /tmp/log')
    
def make_rf_flux_catalog():
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    #### UBV (Maiz-Apellaniz), ugriz (SDSS), JHK (2MASS), UV1600, UV2800
    filters = [153,154,155,156,157,158,159,160,161,162,163,218,219]
    
    unicorn.analysis.generate_rf_color_inputs(FILTERS_RES='FILTER.RES.v8.R300', rf_filter_ids=filters, filename='rf_dummy')
    
    zp_rf = eazy.EazyParam('OUTPUT/rf_dummy.param')
    head_line = '# id z_grism  DM '
    legend = ['# %s\n# \n# Only use grism+photometry redshifts\n# \n' %(time.ctime())]
    
    for filt in zp_rf.filters:
        head_line += ' L%0d' %(filt.fnumber)
        legend.append('# %0d: %s, %12.2f\n' %(filt.fnumber, filt.name, filt.lambda_c))
    
    legend.append('# \n')
    
    fp = open('full_rf_fluxes.cat','w')
    fp.write(head_line+'\n')
    fp.writelines(legend)
    
    files=glob.glob('OUTPUT/*G141*param')
    if len(files) == 0:
        status = os.system('ls OUTPUT |grep G141 |grep param > /tmp/list')
        fpf = open('/tmp/list')
        lines = fpf.readlines()
        fpf.close()
        for line in lines:
            files.append('OUTPUT/'+line[:-1])
    
    for file in files:
        print unicorn.noNewLine +file
        
        root=os.path.basename(file).split('G141')[0]+'G141'
        id=int(os.path.basename(file).split('G141_')[1].split('.param')[0])
        
        try:
            zfit, DM, obs_sed_rf, filts = unicorn.analysis.get_rf_fluxes(root=root, id=id, dummy='rf_dummy', verbose=False, check_plot=False)
        except:
            zfit, DM, obs_sed_rf = np.zeros(3)-1, np.zeros(3)-1, np.zeros((len(zp_rf.filters),3))-1
            
        #### Only get the results for the first fit = spectrum + photometry
        cat_line = ' %s_%05d  %9.5f %6.2f' %(root, id, zfit[0], DM[0])
        for i in range(len(zp_rf.filters)):
            cat_line += ' %12.5e' %(obs_sed_rf[i, 0])
            
        fp.write(cat_line+'\n')
        
    fp.close()
    
    
def get_rf_fluxes(root='GOODS-S-24-G141', id=27, dummy='rf_dummy', verbose=True, check_plot=False):
    """
    Compute the rest-frame colors using the z=0 interpolated fluxes in the dummy 
    catalog and the coefficients from the actual fit.
    """
    import cosmocalc
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    zp_rf = eazy.EazyParam('OUTPUT/%s.param' %(dummy))
    zp_obj = eazy.EazyParam('OUTPUT/%s_%05d.param' %(root, id))

    tempfilt_rf, coeffs_rf, temp_seds_rf, pz_rf = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s' %(dummy), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    tempfilt_obj, coeffs_obj, temp_seds_obj, pz_obj = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    NTEMP = tempfilt_rf['NTEMP']
    ### Possible that different templates used for each
    if zp_rf.templates != zp_obj.templates:
        if verbose:
            print ' *** WARNING, different template sets used ***'
            print '\n Rest-frame:\n=========='
            for t in zp_rf.templates:
                print t
            #
            print '\n Object:\n=========='
            for t in zp_obj.templates:
                print t
            #
        if tempfilt_rf['NTEMP'] > tempfilt_obj['NTEMP']:
            NTEMP = tempfilt_obj['NTEMP']
    
    ##### Redshift from the fit
    zfit = tempfilt_obj['zgrid'][coeffs_obj['izbest']]
    DM = zfit*0
    for i, zi in enumerate(zfit):
        cc = cosmocalc.cosmocalc(zfit[i], H0=zp_obj['H0'], WM=zp_obj['OMEGA_M'])
        DM[i] = 5.*np.log10(cc['DL_Mpc']*1.e5/np.sqrt(1+zfit[i])) #### Last (1+z) term is used because f_nu is interpolated in observed-frame units.  

    ##### The following checks if a given emission line falls within DELTA_LINE Ang.
    ##### (observed frame) of the central wavelength of one of the observed 
    ##### filters.  If so, we assume that the photometry constrains the strength
    ##### of the line in the fit.  If not, remove the line from the rf_template flux.
    
    DELTA_LINE = 300 ### Angstroms, observed frame
    lc_rf = tempfilt_rf['lc'].copy()
    lc_obj = tempfilt_obj['lc'].copy()
    dlam_spec = lc_obj[-1]-lc_obj[-2]
    is_spec = np.append(np.abs(1-np.abs(lc_obj[1:]-lc_obj[0:-1])/dlam_spec) < 0.05,True)
    
    line_names = ['Ha','OIII','Hb','OII_3727']
    line_waves = [6563., 5007., 4861., 3727.]
    NLINES = len(line_names)
    
    # print DELTA_LINE/dlam_spec
    
    for i, temp in enumerate(zp_obj.templates[0:NTEMP]):
        value = 1.
        for j in range(NLINES):
            if line_names[j] in temp:
                line_obs = line_waves[j]*(1+zfit[0])
                dline = np.abs(line_obs-lc_obj)
                dline[is_spec] *= DELTA_LINE/dlam_spec
                if verbose:
                    print 'Closest filter, %8s: %7.1f A' %(line_names[j], dline.min())
                if dline.min() > DELTA_LINE:
                    if verbose:
                        print 'z=%.2f, setting %s template to zero.' %(zfit[0], line_names[j])
                    coeffs_obj['coeffs'][i,:] *= value
                    tempfilt_rf['tempfilt'][:,i,0] = (lc_rf/5500.)**2
                        
    ##### Rest-frame fnu fluxes, zeropoint from the EAZY run.
    idx=0
    obs_sed_rf = np.dot(tempfilt_rf['tempfilt'][:,0:NTEMP,0],\
                     coeffs_obj['coeffs'][0:NTEMP,:])
    
    #### Check results    
    if check_plot:
        lambdaz, temp_sed, lci, obs_sed, fobs_new, efobs_new = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        plt.plot(lci, obs_sed, marker='o', linestyle='None', alpha=0.5)
        plt.semilogx()

        lc_rf = np.arange(len(zp_rf.filters)*1.)
        for i, filt in enumerate(zp_rf.filters):
            lc_rf[i] = filt.lambda_c

        plt.plot(lc_rf*(1+zfit[0]), obs_sed_rf[:,0]/(lc_rf*(1+zfit[0])/5500.)**2, marker='o', color='red', markersize=10, linestyle='None', alpha=0.3)
        
    return zfit, DM, obs_sed_rf, zp_rf.filters
    
def make_eazy_inputs(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', check=False, bin_spec=1, spec_norm=1., zmin=None, zmax=None, zstep=0.0025, compress=1.0, TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', TILT_COEFFS=[0, 1], eazy_working_directory=None, SCALE_SPEC_ERROR=1):
    
    import unicorn.analysis
    
    # root='WFC3-ERSII-G01-G141'; id=197; OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; check=False; bin_spec=1; spec_norm=1.; zmin=None; zmax=None; zstep=0.0025; compress=1.0; TEMPLATES_FILE='templates/o2_fit_lines.spectra.param'; TILT_COEFFS=[0, 1]
    from scipy import polyval
    
    unicorn.analysis.BAD_SPECTRUM = False
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    os.chdir(eazy_working_directory)
    
    if (('MARSHALL' in root) | ('UDS' in root)) & ('UDS' not in OLD_RES):
        OLD_RES='FILTER_UDS.RES'
        
    ORIG_PATH = os.getcwd()
    
    PATH = unicorn.analysis.get_grism_path(root)
    #print PATH
    os.chdir(PATH)
    
    ## read photometric, redshift, SPS catalogs
    cat, zout, fout = unicorn.analysis.read_catalogs(root=root)
        
    ## path where other eazy outputs live
    OUTPUT_DIRECTORY = os.path.dirname(zout.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zout.filename).split('.zout')[0]
    
    grismCat, SPC = unicorn.analysis.read_grism_files(root=root)
    
    match = threedhst.catIO.Readfile('HTML/SED/'+root+'_match.cat')
    
    #### Dummy to get some of the output variables    
    ok = (match.rmatch < 1) & (match.logm > 10.)
    #print root, len(match.rmatch), np.max(match.logm), match.id_f140w[ok][0]
    result = False
    i=0
    while (result is False) & (i < len(match.rmatch[ok])):
        result = unicorn.analysis.specphot(id=match.id_f140w[ok][i],
            grism_root=root, SPC = SPC, 
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same',
            GET_SPEC_ONLY = True)
        i+=1
    #
    lam, spflux, sperr, lci, fobs, efobs, photom_idx = result
    
    #### Match 3D-HST object to the ancillary catalogs to get photometry
    result = unicorn.analysis.specphot(id=id,
        grism_root=root, SPC = SPC, 
        cat = cat,
        grismCat = grismCat, zout = zout, fout = fout, 
        OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
        MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
        OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
        CACHE_FILE = 'Same',
        GET_SPEC_ONLY = True)
    
    if result is not False:
        #### Object is matched in the photometry
        lam, spflux, sperr, lci, fobs, efobs, photom_idx = result
        z_spec = zout.z_spec[photom_idx]
        unicorn.analysis.HAS_PHOTOMETRY = True
        unicorn.analysis.PHOTOMETRY_ID = zout.id[photom_idx]
        
    else:
        #### No match
        sp = threedhst.catIO.Readfile('HTML/ascii/%s_%05d.dat' %(root, id))
        lam, spflux, sperr = sp.lam, (sp.flux-sp.contam)/1.e-17, np.sqrt(sp.error**2+(1.0*sp.contam)**2)/1.e-17
        photom_idx = 0
        z_spec = -1.0
        fobs = fobs*0-999
        efobs = fobs
        print 'No photometry!'
        print TILT_COEFFS
        unicorn.analysis.HAS_PHOTOMETRY = False
        unicorn.analysis.PHOTOMETRY_ID = None
        
    use = (lam > 1.1e4) & (lam < 1.65e4) & (spflux != 0.0) & np.isfinite(spflux) & np.isfinite(sperr)
    
    #### Scale the errors of the spectrum
    #sperr *= 5
    x0, sig = 1.1e4, 1000
    # err_scale = np.exp(-(lam-x0)**2/2/sig**2)
    # err_scale /= err_scale.max()
    # err_scale[lam <= x0] = 1
    # err_scale = err_scale*3+1
    #plt.plot(lam, err_scale)
    sperr*=SCALE_SPEC_ERROR
    
    #### allow additional normalization term
    spflux *= spec_norm
    sperr *= spec_norm
    
    #### Apply the linear fit
    spflux *= polyval(TILT_COEFFS, lam-1.4e4)
    
    #### check plot
    if check:
        plt.semilogx([1],[1])
        plt.errorbar(lci, fobs, yerr=efobs, marker='o', linestyle='None', color='blue', ecolor='blue')
        ma = fobs.max()
        plt.plot(lam[use], spflux[use], color='red')
        plt.ylim(-0.1*ma, 1.5*ma)
        plt.xlim(3000, 9.e4)
    
    #### Rebin if bin_spec is specified
    bin_spec = int(bin_spec)
    if bin_spec > 1:
        NUSE = len(lam[use])
        lbin = np.zeros(NUSE/bin_spec)
        fbin = np.zeros(NUSE/bin_spec)
        ebin = np.zeros(NUSE/bin_spec)
        for i in range(0,NUSE,bin_spec):
            try:
                lbin[i/bin_spec] = np.mean(lam[use][i:i+bin_spec])
                ebin[i/bin_spec] = 1./np.sum(1./sperr[use][i:i+bin_spec]**2)
                fbin[i/bin_spec] = np.sum(spflux[use][i:i+bin_spec]/sperr[use][i:i+bin_spec]**2)*ebin[i/bin_spec]
            except:
                pass
        
        use = lbin > 0
        lam = lbin; spflux = fbin; sperr = ebin
    
    if len(lam[use]) < 10:
        unicorn.analysis.BAD_SPECTRUM = True
        return False
        
    #### Take a wavelength bin and convolve it with the thumbnail
    dlam = (lam[use][1]-lam[use][0])*bin_spec
    NBIN = 30
    DX = 100
    xarr = (np.arange(NBIN*2)-NBIN)*1./NBIN*3*dlam
    yarr = xarr*0
    yarr[np.abs(xarr) < dlam/2.] = 1
    
    os.chdir(PATH)
    lsm, ysm = unicorn.analysis.convolveWithThumb(id, xarr+lam[use][0], yarr, SPC, verbose=False)
    #### Recenter the convolved response
    #xoff = lam[use][0]-np.trapz(lsm, ysm*lsm)/np.trapz(lsm, ysm)
    #print 'XOFF: %f %f' %(xoff, lam[use][0])
    #lsm, ysm = unicorn.analysis.convolveWithThumb(id, xarr+lam[use][0]+xoff, yarr, SPC, verbose=False)
    #lsm += xoff
    keep = ysm > (1.e-4*np.max(ysm))
    
    #### Make lines slightly narrower than the direct image
    lsm -= lam[use][0]
    lsm *= compress
    lsm += lam[use][0]
    
    ###################################
    #### Make a RES file for the convolved resolution and make a catalog file with
    #### photometry and spectroscopy
    
    os.chdir(ORIG_PATH)
    
    ##### Have to fill the empty parts of the spectrum with something other than 0
    fill_value = np.median(spflux[(lam > 1.2e4) & (lam < 1.6e4) & (spflux > 0)])
    
    ##### print the spectrum to a file in the templates directory
    fp = open('templates/%s_%05d' %(root, id)+'_spectrum.dat','w')
    fp.write('50 %.5e\n%.5e %.5e\n' %(fill_value, lam[0]-dlam, fill_value))
    #spflux[~np.isfinite(spflux) | (spflux < 0)] = 1.e-8
    for i in range(len(lam)):
        if np.isfinite(spflux[i]) & (spflux[i] > 0):
            fp.write('%.5e %.5e\n' %(lam[i], spflux[i]))
    
    fp.write('%.5e %.5e\n2.e7 %.5e\n' %(lam[-1]+dlam, fill_value, fill_value))
    fp.close()
    
    #### Make the template file for this single grism spectrum
    fp = open('templates/%s_%05d' %(root, id)+'.spectra.param','w')
    fp.write('1 templates/%s_%05d_spectrum.dat 1.0 0 1.0\n' %(root, id))
    fp.close()
    
    fp = open(OLD_RES)
    res_lines = fp.readlines()
    fp.close()
    
    nfilt=0
    res_name = []
    res_lc = []
    for line in res_lines:
        if line.find('lambda_c') > 0:
            res_name.append(line.split()[1])
            nfilt += 1
    
    #### Read the parameter file associated with the zout file to get the right
    #### filters, which should be in the same order as the vectors read from the 
    #### binary files
    eazy_param = threedhst.eazyPy.EazyParam(zout.filename.replace('zout','param'))
    
    cat_head = '# phot_id id ra dec z_spec'
    cat_line = ' %d %s_%05d %f %f %8.3f' %(cat.id[photom_idx], root, id, cat.ra[photom_idx], cat.dec[photom_idx], z_spec)
    no_spec = cat_line+''
    no_phot = cat_line+''
    for i in range(len(lci)):
        fnumber = eazy_param.filters[i].fnumber
        #print res_lc[fnumber-1], lci[i]
        cat_head += ' F%0d E%0d' %(fnumber, fnumber)
        fnui = (lci[i]/5500.)**2
        cat_line += ' %9.3e %9.3e' %(fobs[i]*fnui, efobs[i]*fnui)
        no_spec += ' %9.3e %9.3e' %(fobs[i]*fnui, efobs[i]*fnui)
        no_phot += ' -99 -99' 
    
    NKEEP = len(lsm[keep])    
    fnu_spec = spflux * (lam/5500.)**2
    efnu_spec = sperr * (lam/5500.)**2
    
    #print '\n\nMAX: %f\n\n' %np.max(fnu_spec[use])
    if np.max(fnu_spec[use]) <= 0:
        fnu_spec[use] = fnu_spec[use]*0.-100
        efnu_spec[use] = efnu_spec[use]*0.-100
        unicorn.analysis.BAD_SPECTRUM = True
        
    for i, l0 in enumerate(lam[use]):
        ### CAT
        cat_head += ' F%0d E%0d' %(nfilt+i+1, nfilt+i+1)
        cat_line += ' %9.3e %9.3e' %(fnu_spec[use][i], efnu_spec[use][i]/3)
        no_spec += ' -99  -99'
        no_phot += ' %9.3e %9.3e' %(fnu_spec[use][i], efnu_spec[use][i]/3)
        #
        ### RES
        lsm0 = lsm[keep]-lam[use][0]+l0
        ysmk = ysm[keep]
        res_lines.append('%d spec lambda_c=%f\n' %(NKEEP, l0))
        for j in range(NKEEP):
            res_lines.append(' %5i %13.5e %13.5e\n' %(j+1, lsm0[j], ysmk[j]))
    
    fp = open('%s_%05d' %(root, id) + '.FILT.RES','w')
    fp.writelines(res_lines)
    fp.close()
    
    fp = open('%s_%05d' %(root, id) + '_threedhst.cat','w')
    fp.write(cat_head+'\n')
    fp.write(cat_line+'\n')
    fp.write(no_spec+'\n')
    fp.write(no_phot+'\n')
    fp.close()
    
    #### Set some EAZY parameters
    
    eazy_param.params['READ_ZBIN'] = 'n'
    eazy_param.params['TEMPLATES_FILE'] = TEMPLATES_FILE
    eazy_param.params['FILTERS_RES'] = '%s_%05d' %(root, id) + '.FILT.RES'
    eazy_param.params['OUTPUT_DIRECTORY'] = 'OUTPUT'
    #eazy_param.params['WAVELENGTH_FILE'] = 'templates/EAZY_v1.1_lines/lambda_v1.1.def'
    eazy_param.params['WAVELENGTH_FILE'] = 'templates/dobos11/lambda_sdss.def'
    eazy_param.params['CATALOG_FILE'] = '%s_%05d' %(root, id) + '_threedhst.cat'
    eazy_param.params['MAIN_OUTPUT_FILE'] = '%s_%05d' %(root, id)
    eazy_param.params['CACHE_FILE'] = '%s_%05d.tempfilt' %(root, id)
    eazy_param.params['GET_ZP_OFFSETS'] = 0
    eazy_param.params['Z_STEP'] = zstep
    eazy_param.params['NMF_TOLERANCE'] = 1.e-4
    
    if zmin is None:
        zmin = zout.l99[photom_idx]*0.8
    if zmax is None:
        zmax = zout.u99[photom_idx]/0.8
        
    if zmin < 0:
        zmin = 0
    
    if zmax < zmin:
        zmax = zmin+1.e-3
        
    eazy_param.params['Z_MIN'] = zmin
    eazy_param.params['Z_MAX'] = zmax 
    
    eazy_param.params['MAGNITUDES'] = 0.0
    eazy_param.params['REST_FILTERS'] = '-,-'
    
    eazy_param.write(file='%s_%05d' %(root, id) + '.eazy.param')
    
def trim_jh_filters(input='FILTER.RES.v8.R300', output='FILTER.RES.v8.R300.trim', wlo=1.08e4, whi=1.65e4):
    """
    Need to trim H & J filters such that they don't sample the parts of the 
    grism spectra where the sensitivity blooms up.
    
    Loop through the filters in an EAZY res file.  If the filter has non-zero 
    transmission at the blue or red edges, trim it off.
    """
    import threedhst.eazyPy as eazy
    
    res = eazy.FilterFile(input)
    for i, filter in enumerate(res.filters):
        peak = filter.transmission.max()
        limits = np.interp([wlo, whi], filter.wavelength, filter.transmission)
        if np.max(limits/peak) > 1.e-1:
            print 'Trim filter: %s' %(filter.name.split()[0])
            keep = (filter.wavelength >= wlo) & (filter.wavelength <= whi)
            res.filters[i].wavelength = filter.wavelength[keep]
            res.filters[i].transmission = filter.transmission[keep]
    
    res.write(file=output)
    
def scale_to_photometry(root='GOODS-S-24-G141', id=23, OLD_RES = 'FILTER.RES.v8.R300', OUT_RES = 'THREEDHST.RES', eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy', compress=1.0, check_results=False, spec_norm=1.0, ORDER=1, pipe=' > log', eazy_working_directory=None):
    import unicorn
    import threedhst.eazyPy as eazy
    from scipy import polyfit, polyval
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    os.chdir(eazy_working_directory)
    
    if ('MARSHALL' in root) | ('UDS' in root):
        OLD_RES='FILTER_UDS.RES'
    
    ##### Need to adjust the J/H filters in the FILTER files near the 
    ##### edges of the grism sensitivity
    if not os.path.exists(OLD_RES+'.trim'):
        unicorn.analysis.trim_jh_filters(input=OLD_RES, output=OLD_RES+'.trim')
        
    #####  Get the f_lambda fluxes with the original (not trimmed) filters
    unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=1.0, spec_norm=spec_norm, zmin=0.0000, zmax=1.e-6, compress=compress, TEMPLATES_FILE='templates/%s_%05d' %(root, id)+'.spectra.param', eazy_working_directory=eazy_working_directory)
    
    if unicorn.analysis.BAD_SPECTRUM:
        return [0,1]
    
    if unicorn.analysis.HAS_PHOTOMETRY is False:
        return [0,1]
    
    status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
    
    lambdaz, temp_sed_0, lci, obs_sed_0, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    #####  Get the f_lambda fluxes with the trimmed filters
    unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES+'.trim', bin_spec=1.0, spec_norm=spec_norm, zmin=0.0000, zmax=1.e-6, compress=compress, TEMPLATES_FILE='templates/%s_%05d' %(root, id)+'.spectra.param', eazy_working_directory=eazy_working_directory)
    
    status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
    
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT',CACHE_FILE = 'Same')
    
    lambdaz, temp_sed, lci, obs_sed, fobs_new, efobs_new = eazy.getEazySED(2, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    if np.sum(obs_sed) == 0:
        obs_sed = obs_sed_0
        temp_sed = temp_sed_0
    
    if np.sum(obs_sed) == 0:
        print "Problem with the eazy fit to get normalization..."
        return [0, 1]
        
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        
    #### If spectrum doesn't cover both J/H, stop and return no correction
    if (lci[is_spec].min() > 1.15e4) | (lci[is_spec].max() < 1.58e4):
        print "Spectrum doesn't have full coverage: [%f, %f]" %(lci[is_spec].min(), lci[is_spec].max())
        return [0, 1]
    
    #### Find filters near J/H
    jhfilt = ~is_spec & (lci > 1.15e4) & (lci < 1.65e4) & (fobs > 0)
    
    #### No photometry:
    if fobs[~is_spec].max() <= 0:
        return [0, 1]
        
    #### Here's the simple linear fit    
    xfit = lci-1.4e4
    yfit = fobs/obs_sed
    
    if len(xfit[jhfilt]) == 0:
        print "No valid J/H filters found."
        return [0, 1]
        
    if len(xfit[jhfilt]) == 1:
        ORDER=0
        
    afit = polyfit(xfit[jhfilt], yfit[jhfilt], ORDER)
    afit[ORDER] *= 1./(coeffs['tnorm'][0]/coeffs['coeffs'][0,2])
    
    #### Reduce the effect slightly
    if ORDER == 1:
        afit[0] *= 0.8
    
    #### Diagnostic plot
    if check_results:
        if unicorn.analysis.USE_PLOT_GUI:
            fig = plt.figure(figsize=[5,5],dpi=100)
        else:
            fig = Figure(figsize=[5,5],dpi=100)

        fig.subplots_adjust(wspace=0.04,hspace=0.12,left=0.05,
                            bottom=0.15,right=0.98,top=0.98)
        #
        ax = fig.add_subplot(211)
        ax.plot(lambdaz, temp_sed)
        ax.plot(lci, fobs, marker='o', markersize=5, alpha=0.5, linestyle='None', color='yellow')
        ax.plot(lci[jhfilt], fobs[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='red')
        ax.plot(lci[jhfilt], obs_sed[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='blue')
        ax.plot(lci, obs_sed, marker='o', markersize=5, alpha=0.5, linestyle='None', color='green')
        ax.set_ylim(0, 1.1*obs_sed.max())
        ax.set_xlim(9.e3, 1.8e4)
                
    if check_results:
        print xfit[jhfilt], yfit[jhfilt]
        print afit
        
        ax.plot(lambdaz, polyval(afit, lambdaz-1.4e4), color='orange', alpha=0.4)
        
        ax = fig.add_subplot(212)
        ax.plot(lambdaz, polyval(afit, lambdaz-1.4e4)*temp_sed)
        ax.plot(lci, polyval(afit, lci-1.4e4)*fobs, marker='o', markersize=5, alpha=0.5, linestyle='None', color='yellow')
        ax.plot(lci[jhfilt], fobs[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='red')
        ax.plot(lci[jhfilt], polyval(afit, lci[jhfilt]-1.4e4)*obs_sed[jhfilt], marker='o', markersize=10, alpha=0.5, linestyle='None', color='blue')
        ax.plot(lci, polyval(afit, lci-1.4e4)*obs_sed, marker='o', markersize=5, alpha=0.5, linestyle='None', color='green')
        ax.set_ylim(0, 1.1*obs_sed.max())
        ax.set_xlim(9.e3, 1.8e4)
        outfile='scale_to_photometry.png'
        if USE_PLOT_GUI:
            plt.savefig(outfile)
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(outfile, dpi=100, transparent=False)
    
    #print afit
    
    return afit
    
def run_eazy_fit(root='COSMOS-23-G141', id=39, OLD_RES = 'FILTER.RES.v9.R300', OUT_RES = 'THREEDHST.RES', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = None, zmin=None, zmax=None, zstep=0.001, compress=1.0, GET_NORM=False, COMPUTE_TILT=True, TILT_ORDER=0, clean=True, force_zrange=False, eazy_working_directory=None, SCALE_SPEC_ERROR=1.):
    
    # OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; TEMPLATES_FILE='templates/o2_fit_lines.spectra.param'; run=True; pipe=' > log'; bin_spec=1; spec_norm=1; eazy_binary = None; zmin=None; zmax=None; compress=1.0; GET_NORM=False; COMPUTE_TILT=True; TILT_ORDER=0; clean=True
    import matplotlib.pyplot as plt
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import unicorn.analysis
    
    t0 = time.time()
    
    if eazy_working_directory is None:
        eazy_working_directory = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS'
        
    if (eazy_binary is None):
        if unicorn.hostname().startswith('uni') | unicorn.hostname().startswith('850dhcp'):
            eazy_binary = '/usr/local/bin/eazy_latest'
        else:
            eazy_binary = '/research/drg/PHOTZ/EAZY/code/SVN/src/eazy'
    
    MAXIT = 3.
    tilt = [0,1]
    
    if ('MARSHALL' in root) | ('UDS' in root):
        OLD_RES='FILTER_UDS.RES'
    
    if run:
        ########################### Scale to photometry
        if COMPUTE_TILT:
            tilt = unicorn.analysis.scale_to_photometry(root=root, id=id, OLD_RES = OLD_RES, OUT_RES = OUT_RES, eazy_binary = eazy_binary, compress=compress, ORDER=TILT_ORDER, pipe=pipe)
        else:
            tilt = [0, 1]
        
        tnorm = time.time()
        
        if unicorn.analysis.BAD_SPECTRUM:
            return False
        
        ########################## Now run
        #### If fitting with photometry
        if unicorn.analysis.HAS_PHOTOMETRY:
            ### first run with eazy line templates 
            ### and coarse sampling, broaden the compression to hopefully catch a line
            if not force_zrange:
                unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.025, compress=3, TILT_COEFFS=tilt, TEMPLATES_FILE='templates/eazy_v1.1_lines_suppl.spectra.param', eazy_working_directory=eazy_working_directory, SCALE_SPEC_ERROR=SCALE_SPEC_ERROR)
            
            #os.system('grep Z_ %s_%05d.eazy.param' %(root, id))
                status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
                ztmp = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))

            ### Remake the input files for the desired compression
            unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.001, compress=compress, TILT_COEFFS=tilt, TEMPLATES_FILE=TEMPLATES_FILE, eazy_working_directory=eazy_working_directory, SCALE_SPEC_ERROR=SCALE_SPEC_ERROR)
            
            
            eazy_param = eazy.EazyParam('%s_%05d.eazy.param' %(root, id))
            eazy_param.params['TEMPLATES_FILE'] = TEMPLATES_FILE
            eazy_param.params['Z_STEP'] = zstep
            eazy_param.params['Z_STEP_TYPE'] = 0
                        
            if force_zrange:
                #### z range set by hand
                eazy_param.params['Z_MIN'] = np.max([zmin,0])
                eazy_param.params['Z_MAX'] = zmax                
            else:
                #### z range set from 99 % confidence interval of first coarse fit
                zmi = np.max([ztmp.l99[1]-1*0.05*(1+ztmp.z_peak[1]),0])
                zma = ztmp.u99[1]+1*0.05*(1+ztmp.z_peak[1])
                #### some very blue SEDs have zmin < 0, zmax~0.  Fit the full range to 
                #### allow LBG fit
                if zma < 0.2:
                    zma = 3
                #   
                eazy_param.params['Z_MIN'] = zmi
                eazy_param.params['Z_MAX'] = zma
                print 'Refit, fine sampling: [%.2f, %.2f]' %(eazy_param.params['Z_MIN'], eazy_param.params['Z_MAX'])

            eazy_param.write(file='%s_%05d' %(root, id) + '.eazy.param')
            
        else:
            ##### No photometry found.  Look for emission lines in the spectrum, 
            ##### and if one line found, fit it assuming OIII.  If more than one line
            ##### found, fit full redshift range.
            spec = catIO.Readfile(unicorn.analysis.get_grism_path(root)+'/HTML/ascii/%s_%05d.dat' %(root, id))
            
            eazy_param = eazy.EazyParam('%s_%05d.eazy.param' %(root, id))
            eazy_param.params['TEMPLATES_FILE'] = TEMPLATES_FILE
            eazy_param.params['Z_STEP'] = zstep
            eazy_param.params['Z_STEP_TYPE'] = 0
            
            found_lines = threedhst.spec1d.findLines(spec, trim_abs=True)
            if found_lines is None:
                return False
                
            if len(found_lines) == 1:
                zmin = found_lines[0].wave/5007.-1
                zmax = found_lines[0].wave/4863.-1
                zmin -= 0.1*(1+zmax)
                zmax += 0.1*(1+zmax)
                
            if len(found_lines) > 1:
                wmin = 1.e5
                wmax = 0
                for line in found_lines:
                    if line.wave < wmin: 
                        wmin = line.wave
                    if line.wave > wmax:
                        wmax = line.wave
                
                zmin = wmin/5007.-1-0.1
                zmax = wmax/4863.-1+0.1
                zmin -= 0.1*(1+zmax)
                zmax += 0.1*(1+zmax)
            
            tilt = [0, 1]
            unicorn.analysis.make_eazy_inputs(root=root, id=id, OLD_RES = OLD_RES, bin_spec=bin_spec, spec_norm=spec_norm, zmin=zmin, zmax=zmax, zstep=0.002, compress=compress, TILT_COEFFS=tilt, TEMPLATES_FILE=TEMPLATES_FILE, eazy_working_directory=eazy_working_directory, SCALE_SPEC_ERROR=SCALE_SPEC_ERROR)
        
        status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
        
        ztmp = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
        zstep_i = zstep
        SHOW_ZOUT_FILE = 'OUTPUT/%s_%05d.zout' %(root, id)
        
        ####### 99% confidence interval is not resolved with z_step.  Shrink the step
        resolve_factor = (ztmp.u95[0]-ztmp.l95[0])/zstep_i
        while (resolve_factor <= 9.99) & (resolve_factor != 0):
            eazy_param.params['Z_MIN'] = ztmp.l99[0]-zstep_i*5
            eazy_param.params['Z_MAX'] = ztmp.u99[0]+zstep_i*5
            eazy_param.params['MAIN_OUTPUT_FILE'] = '%s_%05d_refine' %(root, id)
            eazy_param.params['CACHE_FILE'] = '%s_%05d_refine.tempfilt' %(root, id)
                           
            zstep_i = (ztmp.u95[0]-ztmp.l95[0])/10.
            eazy_param.params['Z_STEP'] = zstep_i
            print 'N=%d, Shrink Z_STEP: %f, [%f, %f]\n' %(resolve_factor, zstep_i, eazy_param.params['Z_MIN'], eazy_param.params['Z_MAX'])
            
            eazy_param.write(file='%s_%05d' %(root, id) + '.eazy.param')
            
            status = os.system(eazy_binary + ' -p '+'%s_%05d' %(root, id)+'.eazy.param '+pipe)
            ztmp = catIO.Readfile('OUTPUT/%s_%05d_refine.zout' %(root, id))
            resolve_factor = (ztmp.u95[0]-ztmp.l95[0])/zstep_i
            
            SHOW_ZOUT_FILE = 'OUTPUT/%s_%05d_refine.zout' %(root, id)
            
        if tilt is False:
            tilt = [0,1]
        
        if len(tilt) == 1:
            tilt = [0,tilt[0]]
        
        print 'Tilt: %e %e\n' %(tilt[0], tilt[1])
        
        fp = open('%s_%05d.tilt' %(root, id),'w')
        fp.write('%s_%05d  %.4e %.4e\n' %(root, id, tilt[0], tilt[1]))
        fp.close()
        
        if unicorn.analysis.BAD_SPECTRUM:
            return False
        
        tfit = time.time()
        
        # lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        #     eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), \
        #                       OUTPUT_DIRECTORY='OUTPUT', \
        #                       CACHE_FILE = 'Same')
        # #
        # dlam_spec = lci[-1]-lci[-2]
        # is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
            
        #### Show the results
        try:
            status = os.system('head -3 %s |tail -1' %(SHOW_ZOUT_FILE))
        except:
            pass
    else:
        tnorm = time.time()
        tfit = time.time()
                        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(0, MAIN_OUTPUT_FILE=os.path.basename(SHOW_ZOUT_FILE).replace('.zout',''), \
                          OUTPUT_DIRECTORY='OUTPUT', \
                          CACHE_FILE = 'Same')
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    #### check normalization
    spec_norm = np.sum(obs_sed[is_spec]*fobs[is_spec]/efobs[is_spec]**2)/np.sum(obs_sed[is_spec]**2/efobs[is_spec]**2)

    print spec_norm    
    if GET_NORM:
        return spec_norm
        
    ################
    ### Make plot
    ################
    
    ### setup
    plt.rcParams['font.size'] = 10
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[7,3],dpi=100)
    else:
        fig = Figure(figsize=[7,3],dpi=100)
        
    fig.subplots_adjust(wspace=0.04,hspace=0.12,left=0.05,
                        bottom=0.15,right=0.98,top=0.98)
    
    #################################### Broad-band + spectrum
    ax = fig.add_subplot(131)
    
    ax.semilogx([1],[1])
    ymax = max(fobs)
    if len(fobs[is_spec & (fobs > 0)]) < 5:
        return False
        
    ymax = max(fobs[is_spec & (fobs > 0)])
    
    ## photometry
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.9, color='orange', markersize=6)
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='red', linestyle='None', markersize=6, alpha=0.2)
    
    ## best-fit SED
    ax.plot(lambdaz, temp_sed, color='red', alpha=0.6)
    
    ## Spectrum + convolved fit
    ax.plot(lci[is_spec], obs_sed[is_spec], color='red', markersize=6, alpha=0.7, linewidth=1)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.8, color='blue', linewidth=1)
        
    ax.set_yticklabels([])
    ax.set_ylabel(r'$f_\lambda$')
    ax.set_xlabel(r'$\lambda$')
    xtick = ax.set_xticks(np.array([0.5, 1., 2, 4])*1.e4)
    ax.set_xticklabels(np.array([0.5, 1., 2, 4]))
    
    ax.set_xlim(3000,9.e4)
    ax.set_ylim(-0.1*ymax, 1.2*ymax)
    
    ymax = max(fobs[is_spec])
    ymin = min(fobs[is_spec])
    
    ################################## Spectrum only
    ax = fig.add_subplot(132)
    ## photometry
    ax.errorbar(lci[~is_spec], fobs[~is_spec], efobs[~is_spec], marker='o', linestyle='None', alpha=0.9, color='orange', markersize=10)
    ax.plot(lci[~is_spec], obs_sed[~is_spec], marker='o', color='red', linestyle='None', markersize=10, alpha=0.2)
    
    ## best-fit template
    ax.plot(lambdaz, temp_sed, color='red', alpha=0.2)
    
    ## spectrum
    ax.plot(lci[is_spec], obs_sed[is_spec], color='red', markersize=6, alpha=0.7, linewidth=2)
    ax.plot(lci[is_spec], fobs[is_spec], marker='None', alpha=0.6, color='blue', linewidth=2)
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.array([1.,1.2, 1.4, 1.6])*1.e4)
    ax.set_xticklabels(np.array([1.,1.2, 1.4, 1.6]))
    ax.set_xlabel(r'$\lambda$')
    
    #### Photometry ID label
    if unicorn.analysis.PHOTOMETRY_ID is not None:
        ax.text(0.1, 0.95, '%0d' %(unicorn.analysis.PHOTOMETRY_ID), transform = ax.transAxes, horizontalalignment='left', verticalalignment='top', fontsize=8)

    ax.set_xlim(0.95e4,1.78e4)
    # ax.set_ylim(0.8*ymin, ymax*1.1)
    ax.set_ylim(-0.1*ymax, 1.2*ymax)
    
    #################################### p(z) for combined, photometry, and spec 
    ax = fig.add_subplot(133)
    
    colors = ['purple','orange','blue']
    alpha = [0.5, 0.5, 0.2]
    zo = threedhst.catIO.Readfile(SHOW_ZOUT_FILE)
    zmin = 4
    zmax = 0
    ymax = 0
    for i in range(3):
        zgrid, pz = eazy.getEazyPz(i, MAIN_OUTPUT_FILE= os.path.basename(SHOW_ZOUT_FILE).replace('.zout',''), 
                          OUTPUT_DIRECTORY='./OUTPUT', 
                          CACHE_FILE='Same')
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        ax.fill_between(zgrid, pz, pz*0., color=colors[i], alpha=alpha[i], edgecolor=colors[i])
        #
        if pz.max() > ymax:
            ymax = pz.max()
        #
        if zgrid[pz > 1.e-6].min() < zmin:
            zmin = zgrid[pz > 1.e-6].min()
        #
        if zgrid[pz > 1.e-6].max() > zmax:
            zmax = zgrid[pz > 1.e-6].max()
    
    ax.plot(zo.z_spec[0]*np.array([1,1]),[0,1.e4], color='green', linewidth=1)
    
    ax.set_yticklabels([])
    ax.set_xlabel(r'$z$')
    ax.xaxis.set_major_locator(MyLocator(4, prune='both'))
    
    ### Plot labels
    ax.text(0.5, 0.9, root+'_%05d' %(id), transform = ax.transAxes, horizontalalignment='center')
    ax.text(0.95, 0.8, r'$z_\mathrm{phot}=$'+'%5.3f' %(zo.z_peak[1]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
    ax.text(0.95, 0.7, r'$z_\mathrm{gris}=$'+'%5.3f' %(zo.z_peak[0]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
    if zo.z_spec[0] > 0:
        ax.text(0.95, 0.6, r'$z_\mathrm{spec}=$'+'%5.3f' %(zo.z_spec[0]), transform = ax.transAxes, horizontalalignment='right', fontsize=9)
        
    ax.set_xlim(zmin-0.2, zmax+0.2)
    ax.set_xlim(zgrid.min(), zgrid.max())
    ax.set_ylim(0,1.1*ymax)
    
    #### Save an image
    outfile = '%s_%05d_eazy.png' %(root, id)
    if USE_PLOT_GUI:
        plt.savefig(outfile)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(outfile, dpi=100, transparent=False)
    
    tplot = time.time()
    
    print '%s, %.2f + %.2f + %.2f = %.2f s' %(outfile, tnorm-t0, tfit-tnorm, tplot-tfit, tplot-t0)
    
    #### Clean up temporary files
    if clean:
        status = os.system('rm %s_%05d' %(root, id) + '_threedhst.cat')
        status = os.system('rm %s_%05d' %(root, id) + '.FILT.RES')
        status = os.system('rm %s_%05d' %(root, id) + '.eazy.param')
        status = os.system('rm templates/%s_%05d' %(root, id) + '.spectra.param')
        status = os.system('rm templates/%s_%05d' %(root, id) + '_spectrum.dat')
    
    return True
#

def check_eazy_fits():
    """
    Check how many templates were run in the eazy fit in the v1.6 directory to
    assess how many objects were overwritten from the original eazy run that 
    generated the catalog redshifts.
    """
    import threedhst.eazyPy as eazy
    
    os.chdir('/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6')
    
    os.system('ls OUTPUT |grep param > /tmp/param_files')
    fp = open('/tmp/param_files')
    param_files = fp.readlines()
    fp.close()
    
    fp = open('param_log','w')
    fp.write('# object NTEMP zmin zmax temp_file\n')
    
    for file in param_files:
        print unicorn.noNewLine+file[:-1]
        param = eazy.EazyParam('OUTPUT/%s' %(file[:-1]))
        fp.write('%-30s %2d %6.2f %6.2f %s\n' %(file[:-1], len(param.templates), param['Z_MIN'], param['Z_MAX'], param['TEMPLATES_FILE']))
        
    fp.close()
    
    ##### Redo fits that were overwritten
    log = catIO.Readfile('param_log')
    redo = (log.temp_file != 'templates/o2_fit_lines.spectra.param')
    
    for object in log.object[redo]:
        if object.startswith('UDF'):
            continue
        #
        object = object.split('.param')[0]
        pointing = object.split('-G141')[0]+'-G141'
        id = int(object.split('_')[1])
        result = unicorn.analysis.run_eazy_fit(root=pointing, id=id, compress=0.7, zmin=0.02, zmax=4, TILT_ORDER=1, pipe=' > log', OLD_RES = 'FILTER.RES.v9.R300', TEMPLATES_FILE='templates/o2_fit_lines.spectra.param')
        os.system('cat OUTPUT/%s.zout' %(object))
        print '\n ------ \n'
        os.system('grep %s ../FIRST_PAPER/GRISM_v1.6/full_redshift.cat' %(object))
        os.system('mv OUTPUT/%s* ../REDSHIFT_FITS_v1.6/OUTPUT/' %(object))
    
def make_all_asciifiles():
    """ 
    Make the ascii files for release v1.6
    """
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    os.chdir('/3DHST/Spectra/Work/ANALYSIS/REDSHIFT_FITS_v1.6')
    
    fields = np.unique(phot.field)
    for field in fields:
        try:
            os.mkdir('ASCII/%s' %(field))
        except:
            pass
            
    fp = open('ASCII/failed.log','w')
    for i in range(len(zout.z_peak[0::3])):
        object = zout.id[0::3][i]
        field = phot.field[phot.idx][i]
        print unicorn.noNewLine+object
        try:
            unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='./ASCII/%s' %(field))
        except:
            fp.write(object+'\n')
    #
    fp.close()
        
def make_eazy_asciifiles(object='COSMOS-8-G141_00498', eazy_output='./OUTPUT/', savepath='./ASCII/'):
    """
    Make ascii files with the best-fit eazy templates and p(z) for the 3D-HST fit.
    
    Currently only the spec + photometry combined fit is saved in the obs_sed and 
    temp_sed files.  All three fits are saved in the pz.dat file.
    
    """
    
    if not savepath.endswith('/'):
        savepath += '/'
        
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    
    eazy_param = eazy.EazyParam('%s/%s.param' %(eazy_output, object))
    for continuum_i, temp in enumerate(eazy_param.templates):
        if 'single_lines' in temp:
            break
            
    #### Photometry + spectrum
    fp = open(savepath+object+'_obs_sed.dat','w')
    fp.write('# lc fnu efnu obs_sed continuum line is_spec\n')
    
    lci = tempfilt['lc']
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)*1
    
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][0]],\
                     coeffs['coeffs'][:,0])# /(lci/5500.)**2
    
    #
    obs_sed_cont = np.dot(tempfilt['tempfilt'][:,0:continuum_i,coeffs['izbest'][0]], coeffs['coeffs'][0:continuum_i,0])# /(lci/5500.)**2
    obs_sed_line = np.dot(tempfilt['tempfilt'][:,continuum_i:,coeffs['izbest'][0]], coeffs['coeffs'][continuum_i:,0])# /(lci/5500.)**2
    
    for i in range(tempfilt['NFILT']):
        fp.write(' %10.6e %8.3e %8.3e  %8.3e  %8.3e  %8.3e %d\n' %(tempfilt['lc'][i], tempfilt['fnu'][i,0], tempfilt['efnu'][i,0], obs_sed[i], obs_sed_cont[i], obs_sed_line[i], is_spec[i]))
    
    fp.close()
        
    #### template fit
    zi = tempfilt['zgrid'][coeffs['izbest'][0]]
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,0])
    temp_sed *= (lambdaz/5500.)**2/(1+zi)**2

    temp_sed_cont = np.dot(temp_seds['temp_seds'][:,0:continuum_i],coeffs['coeffs'][0:continuum_i,0])
    temp_sed_cont *= (lambdaz/5500.)**2/(1+zi)**2

    temp_sed_line = np.dot(temp_seds['temp_seds'][:,continuum_i:],coeffs['coeffs'][continuum_i:,0])
    temp_sed_line *= (lambdaz/5500.)**2/(1+zi)**2
    
    fp = open(savepath+object+'_temp_sed.dat','w')
    fp.write('# lam fnu_temp continuum line\n')
    for i in range(temp_seds['NTEMPL']):
        fp.write(' %10.6e %8.3e %8.2e %8.2e\n' %(lambdaz[i], temp_sed[i], temp_sed_cont[i], temp_sed_line[i]))
    fp.close()
    
    #### p(z)
    
    zgrid, pz0 = eazy.getEazyPz(0, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    zgrid, pz1 = eazy.getEazyPz(1, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    zgrid, pz2 = eazy.getEazyPz(2, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY=eazy_output, CACHE_FILE = 'Same')
    
    fp = open(savepath+object+'_pz.dat','w')
    fp.write('#  z pz_both pz_phot pz_spec\n')
    for i in range(len(zgrid)):
        fp.write('%f %.3e %.3e %.3e\n' %(zgrid[i], pz0[i], pz1[i], pz2[i]))
    
    fp.close()
            
def make_eazy_2d_continuum(root='UDF-G141', id=1279):
    """ 
    Make a 2D continuum image from the best-fit EAZY template.
    """
    from scipy import polyval
    
    FIT_DIR = unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/'
    
    object = '%s_%05d' %(root, id)
    if not os.path.exists(FIT_DIR+'%s.tilt' %(object)):
        print 'Need to run_eazy_fit(root=\'%s\', id=%d)' %(root, id)
        return -1
    
    tilt = np.cast[float](np.loadtxt(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/%s.tilt' %(object), dtype=np.str)[1:])
    
    PATH = unicorn.analysis.get_grism_path(root)
    twod = pyfits.open(PATH+'HTML/images/%s_2D.fits.gz' %(object))
    thumb = pyfits.open(PATH+'HTML/images/%s_thumb.fits.gz' %(object))
    
    #### Eazy outputs
    eazy_param = eazy.EazyParam(FIT_DIR+'OUTPUT/%s_%05d.param' %(root, id))
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY=FIT_DIR+'OUTPUT', CACHE_FILE = 'Same')
    
    iz = coeffs['izbest'][0]
    z_peak = tempfilt['zgrid'][iz]
    continuum_full = np.dot(temp_seds['temp_seds'][:,0:7], coeffs['coeffs'][0:7,0])                             
    continuum_filt = np.dot(tempfilt['tempfilt'][:,0:7,iz], coeffs['coeffs'][0:7,0])                             
    
    lci = tempfilt['lc']
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    #### template interpolated at image wavelengths
    image_lambda = (np.arange(twod[1].header['NAXIS1'])+1-twod[1].header['CRPIX1'])*twod[1].header['CDELT1']+twod[1].header['CRVAL1']
    
    sens = pyfits.open(unicorn.GRISM_HOME+'CONF/WFC3.IR.G141.1st.sens.2.fits')
    sens_int = np.interp(image_lambda, sens[1].data.WAVELENGTH, sens[1].data.SENSITIVITY)
    

    yflux = np.interp(image_lambda, lci[is_spec], continuum_filt[is_spec], left=0, right=0) *3.e18*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))/image_lambda**2*(1+z_peak)
    
    yelec = yflux*sens_int
    
    #### Add the tilt computed from the fit
    yelec /= polyval(tilt, image_lambda-1.4e4)
    
    #### Normalize by the model profile
    profile = np.sum(twod[5].data, axis=1)
    profile /= np.sum(profile)
    
    #### Normalize by the object profile
    profile = np.sum(thumb[0].data, axis=1)
    profile /= np.sum(profile)
    
    #### Normalize to the spectrum itself
    mask = yelec > 0
    for i in range(len(profile)):
        profile[i] = np.sum(yelec[mask]*(twod[1].data[i,mask]-twod[4].data[i,mask]))/np.sum(yelec[mask]**2)
        
    continuum_model = np.dot(profile.reshape(-1,1), yelec.reshape(1,-1))
    
    pyfits.writeto('%s/CONTINUUM/%s_continuum.fits' %(FIT_DIR,object), continuum_model, header=twod[1].header, clobber=True)
    status = os.system('gzip -f %s/CONTINUUM/%s_continuum.fits' %(FIT_DIR,object))
    
    print '%sCONTINUUM/%s_continuum.fits' %(FIT_DIR,object)
    
    xxx = """
    id=1874
    
    im = pyfits.open('/Users/gbrammer/Downloads/UDF-G141_%05d_2D.fits.gz' %(id))
    model = pyfits.open('/research/HST/GRISM/3DHST/ANALYSIS/REDSHIFT_FITS/CONTINUUM/UDF-G141_%05d_continuum.fits.gz' %(id))
    ds9.frame(1)
    ds9.v(im[1].data-im[4].data, vmin=-0.01, vmax=0.1)
    ds9.frame(2)
    ds9.v(im[1].data-model[0].data-im[4].data, vmin=-0.01, vmax=0.1)
    
    """
    
def run_FAST_fit(root='COSMOS-8-G141', id=498, OLD_RES = 'FILTER.RES.v9.R300', OUT_RES = 'THREEDHST.RES', TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param', run=True, pipe=' > log', bin_spec=1, spec_norm=1, eazy_binary = None, zmin=0.2, zmax=5, compress=0.8, GET_NORM=False, COMPUTE_TILT=True, TILT_ORDER=1, force_zrange=False):
    """
    Run FAST fit on the photometry + spectra.  If the catalog and template file
    for a given object isn't found, run `run_eazy_fit` first.
    
    root='COSMOS-8-G141'; id=498; OLD_RES = 'FILTER.RES.v8.R300'; OUT_RES = 'THREEDHST.RES'; TEMPLATES_FILE='templates/o2_fit_lines_suppl.spectra.param'; run=True; pipe=' > log'; bin_spec=1; spec_norm=1; eazy_binary = None; zmin=0.2; zmax=5; compress=0.8; GET_NORM=False; COMPUTE_TILT=True; TILT_ORDER=1; force_zrange=False
    """
    import unicorn
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FAST')
    
    object = '%s_%05d' %(root, id)
    
    zout = None
    
    if (not os.path.exists(object+'_threedhst.cat')) | (not os.path.exists(object+'.FILT.RES')):
        unicorn.analysis.run_eazy_fit(root=root, id=id, OLD_RES = OLD_RES, OUT_RES = OUT_RES, TEMPLATES_FILE=TEMPLATES_FILE, run=True, pipe=pipe, bin_spec=1, spec_norm=1, eazy_binary = eazy_binary, zmin=0.2, zmax=5, compress=0.75, GET_NORM=GET_NORM, COMPUTE_TILT=True, TILT_ORDER=1, clean=False, force_zrange=force_zrange, eazy_working_directory=unicorn.GRISM_HOME+'ANALYSIS/FAST')
    
        fp = open(object+'_threedhst.cat')
        lines = fp.readlines()
        fp.close()
        
        zout = catIO.Readfile('OUTPUT/%s.zout' %(object))
        
        #### Put z_peak for phot + grism as z_spec in the catalog    
        lines[0] = lines[0].replace(' z_spec ',' z_spec_old ').replace(' id ',' id_gris ')[:-1]+' id z_spec\n'
        for i in range(1,4):
            lines[i] = lines[i][:-1]+' %d %.4f\n' %(i,zout.z_peak[0])

        #
        fp = open(object+'_threedhst.cat','w')
        fp.writelines(lines)
        fp.close()
    
    if zout is None:
        zout = catIO.Readfile('OUTPUT/%s.zout' %(object))
    
    param = threedhst.eazyPy.EazyParam('OUTPUT/%s.param' %(object))
    
    #### Make a fast.param file
    fp = open('fast.param')  ### default parameters
    lines = fp.readlines()
    fp.close()
    
    ## change for each object
    lines[80] = 'CATALOG        = \'%s_threedhst\'\n' %(object)
    lines[81] = 'AB_ZEROPOINT   = %f\n' %(param.params['PRIOR_ABZP'])
    lines[82] = 'FILTERS_RES    = \'%s.FILT.RES\'\n' %(object)
    lines[231] = 'Z_MIN          = %f\n' %(zout.z_peak[0])
    lines[232] = 'Z_MAX          = %f\n' %(zout.z_peak[0])
    
    fp = open('%s_fast.param' %(object),'w')
    fp.writelines(lines)
    fp.close()
    
    os.system('/usr/local/share/FAST/FAST_v0.9b/fast %s_fast.param' %(object))
    
    #### Fix the FAST output file, put the header first and add the object name
    print 'Fix FAST_OUTPUT/%s_threedhst.fout' %(object)
    
    fp = open('FAST_OUTPUT/%s_threedhst.fout' %(object))
    lines = fp.readlines()
    fp.close()
    
    for i,line in enumerate(lines):
        if line.startswith('#    id'):
            start = i+1
            header = line
            
    fp = open('FAST_OUTPUT/%s_threedhst.fout' %(object),'w')
    fp.write(header.replace(' id ', ' object id'))
    fp.writelines(lines[0:start-1])
    for line in lines[start:]:
        fp.write(object + ' ' + line)
        
    fp.close()
    
    #### Show FAST results
    print lines[start]
    
    #### Make the eazy ASCII files
    print 'Make EAZY ascii files in ./ASCII ...'
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='./ASCII/')
    
class MyLocator(mticker.MaxNLocator):
    """
    Set maximum number of ticks, from
    http://matplotlib.sourceforge.net/examples/pylab_examples/finance_work2.html
    """
    def __init__(self, *args, **kwargs):
        mticker.MaxNLocator.__init__(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return mticker.MaxNLocator.__call__(self, *args, **kwargs)

#########################################
#                                       #
#                                       #
#     Fit equivalent widths             #
#                                       #
#                                       #
#########################################
def eqw_Eval(x, p):
    """ Compute the sum of the templates """
    y = np.dot(p, x)
    return y

def eqw_myfunct(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = eqw_Eval(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return [status, (y-model)/err]

def eqw_run_mpfit(fobs, efobs, fit_templates, show_results=False):
    import copy
    import unicorn.mpfit as mpfit
    import matplotlib.pyplot as plt
    
    x = fit_templates
    y = fobs
    ey = efobs
    
    #### ignore bad data
    keep = (y > -99) & (ey > 0) & np.isfinite(y) & (y != 0)    
    fa = {'x':x[:,keep], 'y':y[keep], 'err':ey[keep]}
        
    NPARAM = x.shape[0]
    
    p0 = np.ones(NPARAM,dtype='float64')  #initial conditions
    pactual = np.ones(NPARAM) #actual values used to make data
    parbase={'value':0., 'fixed':0, 'limited':[0, 0], 'limits':[0,0.]}
    parinfo=[]
    for i in range(NPARAM):
        parinfo.append(copy.deepcopy(parbase))

    for i in range(NPARAM): 
        parinfo[i]['value']=p0[i]
    
    #### Fix continuum normalization to unity
    #parinfo[0]['value'] = 1.
    #parinfo[0]['fixed'] = 1
    
    m = mpfit.mpfit(eqw_myfunct, p0, parinfo=parinfo,functkw=fa, quiet=1)
    #print m.params, m.params/m.perror
    
    ### Correlation matrix
    m.pcor = m.covar * 0.
    for i in range(NPARAM):
        for j in range(NPARAM):
            m.pcor[i,j] = m.covar[i,j]/np.sqrt(m.covar[i,i]*m.covar[j,j])
    
    ### relative error
    m.relerr = np.zeros(NPARAM)
    m.relerr[0] = m.perror[0]/m.params[0]
    
    for i in range(1, NPARAM):
        m.relerr[i] = np.sqrt((m.perror[0]/m.params[0])**2+(m.perror[i]/m.params[i])**2)
        
    if show_results:
        print 'Show'
        xarr = np.arange(x.shape[1])
        plt.errorbar(xarr, y, yerr=ey, linestyle='None', marker='o', color='blue')
        yfit = eqw_Eval(x, m.params)
        plt.plot(xarr, yfit, color='red', alpha=0.5, linewidth=3)
    
    return m
    
def equivalent_width(root='GOODS-S-24-G141', id=29):
    """ 
    Measure the line equivalent widths from the template fit.
    
    Would also be nice to add line fluxes themselves.  This should come instantly from
    the coefficients of the line templates, since the templates are normalized to 
    area unity.  Therefore the flux in "mag 25" Jy is just 'coeff/tnorm' or something
    like that.  And then f_lambda fluxes is just a unit game and is probably already
    in 'eazyPy'.
    
    Another idea is to subtract the continuum from the spectrum and integrate the 
    line directly. 
    
    USAGE:
    object, z_grism, halpha_eqw, halpha_flux, oiii_eqw, oiii_flux, hbeta_eqw, hbeta_flux = equivalent_width(root='x', id=1)
    
    """
    import threedhst.eazyPy as eazy
    import threedhst.catIO as catIO
    import cosmocalc
    import unicorn.analysis
    
    zout = catIO.Readfile('OUTPUT/%s_%05d.zout' %(root, id))
    eazy_param = eazy.EazyParam('OUTPUT/%s_%05d.param' %(root, id))
    
    # tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'OUTPUT/line_eqw.tempfilt')
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
    
    coeffs['coeffs'] *= 1./(1+zout.z_peak[0])**2
    
    continuum = np.dot(temp_seds['temp_seds'][:,0:6], coeffs['coeffs'][0:6,0])                             
    
    zpeak_i = tempfilt['zgrid'][coeffs['izbest'][0]]
    lci = tempfilt['lc']
    
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
    
    fobs = tempfilt['fnu'][:,0]/(lci/5500.)**2
    efobs = tempfilt['efnu'][:,0]/(lci/5500.)**2
    
    idx_ha, idx_hb, idx_oiii = 6,7,8
    for i,temp in enumerate(eazy_param.templates):
        if temp == 'templates/single_lines/OIII_5006.dat':
            idx_oiii = i
        if temp == 'templates/single_lines/Hb_4861.dat':
            idx_hb = i
        if temp == 'templates/single_lines/Ha_6562.dat':
            idx_ha = i
    
    print 'IDX: ',idx_ha, idx_hb, idx_oiii, tempfilt['NTEMP']
    
    obs_sed_continuum = np.dot(tempfilt['tempfilt'][:,0:idx_ha,coeffs['izbest'][0]],coeffs['coeffs'][0:idx_ha,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_ha = np.dot(tempfilt['tempfilt'][:,idx_ha,coeffs['izbest'][0]],coeffs['coeffs'][idx_ha,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_oiii = np.dot(tempfilt['tempfilt'][:,idx_oiii,coeffs['izbest'][0]],coeffs['coeffs'][idx_oiii,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    obs_sed_hb = np.dot(tempfilt['tempfilt'][:,idx_hb,coeffs['izbest'][0]],coeffs['coeffs'][idx_hb,0])/(lci/5500.)**2*(1+zpeak_i)**2
    
    # obs_sed_ha /= obs_sed_ha.max()
    # obs_sed_oiii /= obs_sed_oiii.max()
    # obs_sed_hb /= obs_sed_hb.max()
    
    fit_templates = np.array([obs_sed_continuum, obs_sed_ha, obs_sed_oiii, obs_sed_hb])

    if np.std(obs_sed_ha[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[1,:] *= 0.
    #    
    if np.std(obs_sed_oiii[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[2,:] *= 0.
    #    
    if np.std(obs_sed_oiii[is_spec]/obs_sed_ha.max()) < 1.e-4:
        fit_templates[3,:] *= 0.
    #    
    
    try:
        mp = unicorn.analysis.eqw_run_mpfit(fobs[is_spec], efobs[is_spec], fit_templates[:, is_spec], show_results=False)
        relerr, perror = mp.relerr, mp.perror
    except:
        relerr, perror = np.zeros(4)-1, np.zeros(4)-1
      
    #print 'relerr', relerr, perror
      
    #mp.relerr, mp.params, mp.perror
    
    #### for testing
    plot = """
        tempfiltx, coeffsx, temp_sedsx, pzx = eazy.readEazyBinary(MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE='%s_%05d' %(root, id), OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
        dlam_spec = lci[-1]-lci[-2]
        is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        
        # plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), continuum, color='black')
        # plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), temp_sed, color='red', alpha=0.6)
        plt.plot(lci[is_spec], fobs[is_spec], color='blue', alpha=0.4)
        plt.errorbar(lci[is_spec], fobs[is_spec], efobs[is_spec], color='blue', alpha=0.4)
        plt.plot(lci[is_spec], obs_sed_continuum[is_spec], color='black', alpha=0.4, linewidth=3)
        plt.plot(lci[is_spec], obs_sed_ha[is_spec], color='black', alpha=0.4, linewidth=3)
        # plt.plot(lci[is_spec], obs_sed_manual[is_spec], color='red', alpha=0.4, linewidth=3)
        
        plt.plot(lci[is_spec], obs_sed[is_spec], color='orange', linewidth=3, alpha=0.4)
        plt.semilogx()
        plt.xlim(3000,4.e4)
        
        #### Test getting fluxes in f_lambda
        #np.trapz(temp_seds['temp_seds'][:,6], temp_seds['templam']*coeffs['tnorm'][6])
        ha_flux_coeffs = coeffs['coeffs'][6,0]/coeffs['tnorm'][6]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2*(6563.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
        
        ha_flux_convolved = temp_seds['temp_seds'][:,6].max()*coeffs['coeffs'][6,0]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2
        
        
        ha_flux_integrated = np.trapz(halpha, temp_seds['templam']*(1+zout.z_peak[0]))
        ha_flux_integrated *= (6563.*(1+zout.z_peak[0])/5500.)**2
        ha_flux_integrated *= 10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2
        
        import cosmocalc
        cc = cosmocalc.cosmocalc(zout.z_peak[0], H0=71, WM=0.27, WV=1-0.27)
        lum_ha = ha_flux_integrated*cc['DL_cm']**2*4*np.pi
        sfr_ha = 7.9e-42*lum_ha
        
        plt.plot(temp_seds['templam']*(1+zout.z_peak[0]), halpha)
        plt.xlim(6300*(1+zout.z_peak[0]),6900*(1+zout.z_peak[0]))
        
    """
    
    halpha = temp_seds['temp_seds'][:,idx_ha]*coeffs['coeffs'][idx_ha,0]
    halpha[halpha < 1.e-8*halpha.max()] = 0
    halpha_eqw = -np.trapz((-halpha/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    
    # fp = open('../EQW_FOR_MATTIA/test_temp.dat','w')
    # fp.write('# lam continuum line\n')
    # for i in range(1,len(temp_seds['templam'])):
    #     fp.write('%.6e %.3e %.3e\n' %(temp_seds['templam'][i], continuum[i], halpha[i]))
    # fp.close()    
    
    ## test different ways of measuring the equivalenth width from the convolved template and from the spectrum itself
    halpha_eqw_smooth = -np.trapz((-obs_sed_ha/obs_sed_continuum)[is_spec], lci[is_spec])
    use = is_spec & (obs_sed_ha > 1.e-5*obs_sed_ha.max())
    halpha_eqw_data = -np.trapz((-(tempfilt['fnu'][:,0]/(lci/5500.)**2-obs_sed_continuum)/obs_sed_continuum)[use], lci[use])
    
    # 
    # idx = np.arange(len(obs_sed_continuum))
    # fp = open('../EQW_FOR_MATTIA/test_obs.dat','w')
    # fp.write('# lam fnu continuum line\n')
    # for i in idx[is_spec]:
    #     fp.write('%.6e %.3e %.3e %.3e\n' %(lci[i], tempfilt['fnu'][i,0], obs_sed_continuum[i], obs_sed_ha[i]))
    # fp.close()
    
    halpha_flux = coeffs['coeffs'][idx_ha,0]/coeffs['tnorm'][idx_ha]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(6563.*(1+zout.z_peak[0]))**2*(6563.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    halpha_err = halpha_eqw*relerr[1]
    if perror[1] == 0:
        halpha_err = -1.
    
    #print 'relerr:' %relerr
    #print 'EQW, FLUX, ERR: %.3f %.3f %.3f' %(halpha_eqw, halpha_err, halpha_flux)
    halpha_eqw = (halpha_eqw, halpha_eqw_smooth, halpha_eqw_data)
    
    # cc = cosmocalc.cosmocalc(zout.z_peak[0], H0=71, WM=0.27, WV=1-0.27)
    # halpha_lum = halpha_flux*cc['DL_cm']**2*4*np.pi
    # halpha_sfr = 7.9e-42*halpha_lum
    
    oiii = temp_seds['temp_seds'][:,idx_oiii]*coeffs['coeffs'][idx_oiii,0]
    oiii[oiii < 1.e-8*oiii.max()] = 0
    oiii_eqw = -np.trapz((-oiii/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    oiii_flux = coeffs['coeffs'][idx_oiii,0]/coeffs['tnorm'][idx_oiii]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(5007.*(1+zout.z_peak[0]))**2*(5007.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    oiii_err = oiii_eqw*relerr[2]
    if perror[2] == 0:
        oiii_err = -1.
    
    hbeta =  temp_seds['temp_seds'][:,idx_hb]*coeffs['coeffs'][idx_hb,0]
    hbeta[hbeta < 1.e-8*hbeta.max()] = 0
    hbeta_eqw = -np.trapz((-hbeta/continuum)[1:-1], temp_seds['templam'][1:-1]*(1+zpeak_i))
    hbeta_flux = coeffs['coeffs'][idx_hb,0]/coeffs['tnorm'][idx_hb]*10**(-0.4*(eazy_param.params['PRIOR_ABZP']+48.6))*3.e18/(4861.*(1+zout.z_peak[0]))**2*(4861.*(1+zout.z_peak[0])/5500.)**2*(1+zout.z_peak[0])
    hbeta_err = hbeta_eqw*relerr[3]
    if perror[3] == 0:
        hbeta_err = -1.
    
    ###### MPFIT to get errors on normalizations
    
    return '%s_%05d' %(root, id), zout.z_peak[0], halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux
    
def make_o2_templates():
    """
    Leaving the OIII strength fixed allows too much freedom in the fits for galaxies at z=0.5 and only ugriz.  The code likes to put the galaxies at z>0.7 with a huge OII line that mimicks the break.  As a fix, try making a set of the `noline` templates that include the eazy OII line from v1.1
    """
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/templates/')
    import glob
    files=glob.glob('EAZY_v1.0_lines/*nolines.dat')
    if not os.path.exists('O2_ONLY'):
        status = os.system('mkdir O2_ONLY')
    
    for file in files:
        noline = np.loadtxt(file)
        noline_wave, noline_flux = noline[:,0], noline[:,1]
        noline_keep = np.abs(noline_wave-3727) > 30
        
        hasline = np.loadtxt(file.replace('v1.0','v1.1').replace('_nolines',''))
        hasline_wave, hasline_flux = hasline[:,0], hasline[:,1]
        hasline_keep = np.abs(hasline_wave-3727) < 30
        
        xnew = np.append(noline_wave[noline_keep],hasline_wave[hasline_keep])
        ynew = np.append(noline_flux[noline_keep],hasline_flux[hasline_keep])
        
        s = np.argsort(xnew)
        xnew, ynew = xnew[s], ynew[s]
        
        fp = open('O2_ONLY/'+os.path.basename(file),'w')
        for i in range(len(xnew)):
            fp.write('%13.5e %13.5e\n' %(xnew[i], ynew[i]))
        
        fp.close()
        
    print "\n\nManually make a copy of 'fit_lines.spectra.param'\n\n"
    
def make_line_templates():
    line_wavelengths = [[6562.800], [5006.843, 4958.911], [4862.68], [3727.0], [6718.29, 6732.67], [1216.], [9068.6, 9530.6], [2799.117], [1549.48], [2326.0], [4341.68], [3889.0]]
    line_ratios = [[1], [2.98, 1], [1], [1], [1, 1], [1], [1, 2.44], [1], [1], [1], [1], [1]]
    line_names = ['Ha', 'OIII', 'Hb', 'OII', 'SII', 'Lya', 'SIII', 'MgII', 'CIV', 'CII','Hg','HeI']
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    
    NLINE = len(line_names)
    #fp = open('templates/EAZY_v1.1_lines/lambda_v1.1.def')
    fp = open('templates/dobos11/lambda_sdss.def')  ### higher resolution
    xspec = np.cast[float](fp.readlines())
    fp.close()
    
    NDEF = len(xspec)
    vel_width = 150 # km/s
    
    fp = open('templates/eazy_v1.0_nolines.spectra.param')
    spec_list = fp.readlines()
    fp.close()
    
    last = np.int(spec_list[-1][0])
    for i in range(NLINE):
        lw = line_wavelengths[i][0]
        yspec = xspec*0.
        ratio = np.array(line_ratios[i])*1./np.sum(line_ratios[i])
        for j in range(len(line_wavelengths[i])):
            lj = line_wavelengths[i][j]
            sigma = vel_width/3.e5*lj
            yspec += ratio[j]*1./np.sqrt(2*np.pi*sigma**2)*np.exp(-1.*(xspec-lj)**2/2/sigma**2)
            
        yspec /= np.trapz(yspec, xspec)
        
        #
        fp=open('templates/single_lines/'+line_names[i]+'_%0d.dat' %(lw), 'w')
        for j in range(NDEF):
            fp.write(' %.5e %.5e\n' %(xspec[j], yspec[j]+1.e-12))
        #
        fp.close()
        #
        spec_list.append('%d   templates/single_lines/%s_%0d.dat 1.0 0 1.0\n' %(last+i+1, line_names[i], lw))
    
    fp = open('templates/fit_lines.spectra.param','w')
    fp.writelines(spec_list)
    fp.close()
    
def test_equivalent_widths():
    """
    Mattia pointed out that he gets equivalent widths somewhat larger than the catalog values when he measures them with his own code.  
    """
    
    object = 'COSMOS-3-G141_00765'
    ### EW(Matt) =75.043951 \pm 3.8351794, EW(Gabe) = 47.080000 \pm 2.4630000
    unicorn.analysis.run_eazy_fit(root=object.split('_')[0], compress=0.7, TILT_ORDER=1, OLD_RES='FILTER.RES.v9.R300', zmin=0.7, zmax=1.2, id=int(object.split('_')[1]), force_zrange=True, COMPUTE_TILT=True)
    
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='../EQW_FOR_MATTIA/')
    
    obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
    
    os.chdir('../REDSHIFT_FITS_v1.6')
    unicorn.analysis.make_eazy_asciifiles(object=object, eazy_output='./OUTPUT/', savepath='../EQW_FOR_MATTIA/')
    
    obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
    
    ### Plot:
    obs_sed = catIO.Readfile(object+'_obs_sed.dat')
    
    lci = obs_sed.lc
    dlam_spec = lci[-1]-lci[-2]
    is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)

    idx = np.arange(len(obs_sed.lc))
    fp = open('junk.txt','w')
    for i in idx[is_spec]:
        fp.write('%.1f %.3e\n' %(obs_sed.lc[i], obs_sed.obs_sed[i]*(5500./obs_sed.lc[i])**2 ))
    fp.close()
    
    temp_sed = catIO.Readfile(object+'_temp_sed.dat')
    idx = np.arange(len(temp_sed.lam))
    keep = (temp_sed.lam > 1.e4) & (temp_sed.lam < 1.6e4)
    fp = open('temp.txt','w')
    for i in idx[keep]:
        fp.write('%.3f %.3e\n' %(temp_sed.lam[i], temp_sed.fnu_temp[i]*(5500./temp_sed.lam[i])**2 ))
    fp.close()
    
    plt.plot(lci[is_spec], obs_sed.fnu[is_spec]*(5500./lci[is_spec])**2, color='black')
    plt.plot(lci[is_spec], obs_sed.obs_sed[is_spec]*(5500./lci[is_spec])**2, color='red')
    plt.plot(lci[is_spec], obs_sed.continuum[is_spec]*(5500./lci[is_spec])**2, color='blue')
    plt.plot(temp_sed.lam, temp_sed.fnu_temp*(5500./temp_sed.lam)**2, color='blue')
    plt.plot(temp_sed.lam, temp_sed.continuum*(5500./temp_sed.lam)**2, color='blue')
    
    tt = catIO.Readfile('test_temp.dat')
    plt.plot(tt.lam, tt.continuum, color='green')
    plt.plot(tt.lam, tt.line, color='green')
    
    oo = catIO.Readfile('test_obs.dat')
    plt.plot(oo.lam, oo.fnu/(oo.lam/5500.)**2, color='orange')
    plt.plot(oo.lam, oo.continuum, color='purple')
    
    plt.xlim(1.1e4,1.7e4)
    plt.ylim(0,100)
    
    np.trapz(-(obs_sed.line/obs_sed.continuum)[is_spec], obs_sed.lc[is_spec])
    keep = (temp_sed.lam > 1.2e4) & (temp_sed.lam < 1.3e4)
    np.trapz(-(temp_sed.line/temp_sed.continuum)[keep], temp_sed.lam[keep])
    
    keep = (tt.lam > 1.2e4) & (tt.lam < 1.3e4)
    keep = (tt.lam > 0.6e4) & (tt.lam < 0.7e4)
    np.trapz(-(tt.line/tt.continuum)[keep], tt.lam[keep]*(1+0.9341))
    
    ############ Mattia's list synced from the v1.6 redshift fits
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/EQW_FOR_MATTIA')
    ew = catIO.Readfile('list_for_gabe.dat')
    
    fp = open('ASCII/gabe_eqw.dat','w')
    fp.write('# id z_grism halpha_eqw halpha_eqw_err\n')
    
    for i, object in enumerate(ew.id):
        try:
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
        except:
            print 'xx %s xx' %(object)
            pass
        #
        fp.write('%s %.3f %.3f\n' %(object, halpha_eqw[0], halpha_err))
        #
        unicorn.analysis.make_eazy_asciifiles(object=object, savepath='ASCII')
    
    fp.close()
        
    new_fits = np.zeros((3, len(ew.id)))
    z_fit = np.zeros(len(ew.id))
    
    for i, object in enumerate(ew.id):
        print object
        try:
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
            z_fit[i] = z_grism
        except:
            continue
        #        
        new_fits[:,i] = np.array(halpha_eqw)
    
    #### Run splot to measure them by hand
    import iraf
    from iraf import onedspec
    for i, object in enumerate(ew.id[8:]):
        try:
            lambdaz, temp_sed, lci, obs_sed, fobs, efobs = eazy.getEazySED(0, MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE = 'Same')
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=object.split('_')[0], id=int(object.split('_')[1]))
        except:
            continue
        #
        #
        dlam_spec = lci[-1]-lci[-2]
        is_spec = np.append(np.abs(1-np.abs(lci[1:]-lci[0:-1])/dlam_spec) < 0.05,True)
        fpf = open('SPLOT/%s_specfit.txt' %(object),'w')
        fpf.write('%s\n' %object)
        fpf.write('%d 5100.\n' %(len(lci[is_spec])))
        #
        fpd = open('SPLOT/%s_data.txt' %(object),'w')
        fpm = open('SPLOT/%s_model.txt' %(object),'w')
        for lam, flux, model, flux_err in zip(lci[is_spec], fobs[is_spec], obs_sed[is_spec], efobs[is_spec]):
            fpf.write('%.6e %.2e %.2e\n' %(lam/(1+z_grism), flux, flux_err))
            fpd.write('%.6e %.2e\n' %(lam, flux))
            fpm.write('%.6e %.2e\n' %(lam, model))
        fpf.close()
        fpd.close()
        fpm.close()
        #
        iraf.rspectext('SPLOT/%s_data.txt' %(object),'SPLOT/%s_data.fits' %(object), dtype="interp", crval1=11000.0, cdelt1 = 22)
        iraf.rspectext('SPLOT/%s_model.txt' %(object),'SPLOT/%s_model.fits' %(object), dtype="interp", crval1=11000.0, cdelt1 = 22)
        #
        print 'z=%.3f, Ha=%.1f' %(z_grism, 6563*(1+z_grism))
        print halpha_eqw, ew.ew_matt[i]
        #
        iraf.splot('SPLOT/%s_data.fits' %(object), save_file='splot.log')
        iraf.splot('SPLOT/%s_model.fits' %(object), save_file='splot.log')
    
    splot_k_data = ew.ew_gabe*0.
    splot_e_data = ew.ew_gabe*0.
    splot_k_model = ew.ew_gabe*0.
    splot_e_model = ew.ew_gabe*0.
    fp = open('splot.log.save')
    lines = fp.readlines()
    N = len(lines)/12
    for i in range(N):
        group = lines[i*12:i*12+12]
        object = group[1].split('DATA/')[1].split('_data')[0]
        splot_k_data[ew.id == object] = -float(group[3].split()[3])
        splot_e_data[ew.id == object] = -float(group[5].split()[3])
        splot_k_model[ew.id == object] = -float(group[9].split()[3])
        splot_e_model[ew.id == object] = -float(group[11].split()[3])
        
    plt.plot(ew.ew_matt, ew.ew_gabe, marker='o', linestyle='None', color='black', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[0,:], marker='o', linestyle='None', color='red', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[1,:], marker='o', linestyle='None', color='blue', alpha=0.3)
    plt.plot(ew.ew_matt, new_fits[2,:], marker='o', linestyle='None', color='green', alpha=0.3)
    
    #plt.plot(splot_k_data, splot_k_model, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_k_data, splot_k_model, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_k_data, ew.ew_gabe, marker='o', linestyle='None', color='black', alpha=0.8)
    plt.plot(splot_e_data, splot_e_model, marker='o', linestyle='None', color='black', alpha=0.8)
    # plt.plot(splot_e_data, new_fits[0,:]/1.07, marker='o', linestyle='None', color='orange', alpha=0.8)
    
    #plt.plot(splot_e_model, splot_e_data, marker='o', linestyle='None', color='black', alpha=0.8)
    
    #plt.plot(splot_e_data, new_fits[2,:], marker='o', linestyle='None', color='red', alpha=0.8)

    #plt.plot(splot_e_data, splot_k_data, marker='o', linestyle='None', color='black', alpha=0.8, ms=10)

    fig = unicorn.catalogs.plot_init(square=True, aspect=1./3, left=0.12, xs=12)
        
    ax = fig.add_subplot(131)
    
    #plt.plot(splot_e_data, new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_gabe*(1+z_fit), marker='o', linestyle='None', color='blue', alpha=0.8, ms=6)
    #plt.plot(splot_e_data, new_fits[2,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_matt, marker='o', linestyle='None', color='red', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW data, splot "e"')
    plt.ylabel('EQW meas, Mattia [red], Gabe/old*(1+z) [blue]')
    
    ax = fig.add_subplot(132)
    
    plt.plot(splot_e_data, new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    #plt.plot(splot_e_data, ew.ew_gabe*(1.9), marker='o', linestyle='None', color='blue', alpha=0.8, ms=10)
    #plt.plot(splot_e_data, new_fits[2,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=10)
    plt.plot(splot_e_data, ew.ew_matt, marker='o', linestyle='None', color='red', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW data, splot "e"')
    plt.ylabel('EQW meas, Mattia [red], Gabe/fixed [black]')
    
    ax = fig.add_subplot(133)
    plt.plot(new_fits[0,:], new_fits[1,:] , marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW, full template')
    plt.ylabel('EQW, convolved')
    
    #fig.savefig('eqw_comparison_old.pdf')
    
    fp = open('gabe_updated.dat','w')
    fp.write('# id z   eqw_matt eqw_gabe  splot_e splot_k  eqw_full eqw_conv\n')
    for i in range(len(ew.id)):
        fp.write('%-20s  %.3f %7.2f %7.2f  %7.2f %7.2f  %7.2f %7.2f\n' %(ew.id[i], z_fit[i],  ew.ew_matt[i], ew.ew_gabe[i], splot_e_data[i], splot_k_data[i], new_fits[0,i], new_fits[1,i]))
    fp.close()
    
    ##### Scatter between measurements
    fig = unicorn.catalogs.plot_init(square=True, aspect=1./2, left=0.12, xs=12*2./3)
        
    ax = fig.add_subplot(121)
    
    #plt.plot(splot_e_data, ew.ew_gabe*(1+z_fit), marker='o', linestyle='None', color='blue', alpha=0.8, ms=6)
    plt.plot(new_fits[1,:], ew.ew_matt, marker='o', linestyle='None', color='black', alpha=0.8, ms=6)
    
    plt.plot([1,400],[1,400], color='orange')
    plt.semilogy()
    plt.semilogx()
    plt.xlim(10,600)
    plt.ylim(10,600)
    plt.xlabel('EQW, Gabe / fixed')
    plt.ylabel('EQW, Mattia')
    
    ax = fig.add_subplot(122)
    
    delta = np.log10(new_fits[1,:]/ew.ew_matt)
    keep = np.isfinite(delta)
    ax.hist(delta, bins=25, range=(-1,1))
    ax.set_xlabel(r'log EW(GABE)/EW(Mattia)')
    ax.set_ylabel(r'N')
    
    biw = threedhst.utils.biweight(delta[keep], both=True)
    nm = threedhst.utils.nmad(delta[keep])
    ax.text(-0.8,14,r'$(\mu,\ \sigma) = %.2f,\ %.2f$' %(biw[0], biw[1]))
    
    fig.savefig('eqw_scatter.pdf')

def find_brown_dwarfs():
    import glob
    import threedhst.catIO as catIO
    
    obj = catIO.Readfile('AEGIS-3-G141_00177.dat')
    
#
def old_analysis_plots():
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

"""
Work with v2.x catalogs
"""

import os
from astropy.io import fits as pyfits
import numpy as np
import glob
import shutil
import re
import time

import matplotlib.pyplot as plt

# from matplotlib.figure import Figure
# from matplotlib.backends.backend_agg import FigureCanvasAgg

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

PATH_TO_CAT = unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/'

zout = None
fout = None
cat = None
zfit = None
lines = None
dq = None
selection_params = None
field = None
rf = None
irf = None
    
class SpeczCatalog():
    """
    Container for spectroscopic redshifts
    """
    def __init__(self, force_new=False):
        self.tree = None
        self.sproot = unicorn.GRISM_HOME+'/ANALYSIS/SPECTROSCOPIC_REDSHIFTS/'
        
        has_pickle = os.path.exists(os.path.join(self.sproot, 'full_zspec.pkl'))
        if (not has_pickle) | (force_new):
            self.read_catalogs()
            self.save_pickle()
        else:
            self.read_pickle()
        
        self.init_tree()
        
    def init_tree(self):
        """
        Set up the tree for retrieving nearest zSpec objects from
        test coordinate positions.
        
        The match is done using the scipy.spatial.cKDTree function, which is
        apparently some 80 times faster than doing it "by hand" with numpy.
        """
        import scipy.spatial
        
        cosd = self.ra*np.cos(self.dec/360*2*np.pi)
        self.xy = np.array([cosd, self.dec]).T
        self.tree = scipy.spatial.cKDTree(self.xy, 10)
        
    def test_by_hand(self, ra, dec, N=1):
        """
        Test to compare speed of scipy.spatial.cKDTree to numpy for 
        finding N nearest neighbors.
        """
        xy_test = [ra*np.cos(dec/360.*2*np.pi), dec]
        dr = np.sqrt(np.sum((self.xy-xy_test)**2, axis=1))
        so = np.argsort(dr)
        return dr[so[0:N]]*3600, so[0:N]
        
    def find_nearest(self, ra, dec, N=1):
        """
        Find N nearest neighbors to (ra, dec) in the zSpec catalogs.  
        
        Example: 
            
            >>> dist, ids = zsp.find_nearest(zsp.ra[100], zsp.dec[100], N=5)
            >>> print dist, ids
            (array([  0.        ,  12.96253365,  17.63697491,  29.72497372,  31.16232403]), array([100,  86, 119, 116,  80], dtype=int32))
            
        """
        if self.tree is None:
            self.init_tree()
            
        xy_test = [ra*np.cos(dec/360.*2*np.pi), dec]
        dist, ids = self.tree.query(xy_test, k=N)
        return dist*3600, ids
    #
    def match_list(self, ra=[], dec=[], N=1, MATCH_SELF=False, verbose=True):
        """
        Make a full matched list, input 'ra' and 'dec' are
        arrays
        
        If MATCH_SELF, find nearest matches *within* the self catalog
        """
        noNewLine = '\x1b[1A\x1b[1M'
        
        if MATCH_SELF:
            ra = self.ra
            dec = self.ra
        
        Nlist = len(ra)
        dr_match = ra*0.
        id_match = np.cast[int](ra*0)
        
        for i in range(Nlist):
            if verbose:
                print noNewLine+'%d of %d' %(i+1, Nlist)
            
            dist, ids = self.find_nearest(ra[i], dec[i], N=1+N)
            dr_match[i] = dist[N-1+MATCH_SELF]
            id_match[i] = ids[N-1+MATCH_SELF]
        
        self.dr_zsp = dr_match
        self.id_zsp = id_match
        return dr_match, id_match
    
    def flag_multiples(self, mag_list, toler=1):
        """
        Find multiple objects matched to the same zsp object.
        Make flags `mult_match` with the number matched and 
        `brightest`=1 for the brightest galaxy of multiple matches.
        """
        self.mult_match = self.id_zsp*0
        self.brightest = self.id_zsp*0+1
        for ix, i in enumerate(self.id_zsp):
            mat = (self.id_zsp == i) & (self.dr_zsp < toler)
            self.mult_match[ix] = mat.sum()
            if self.mult_match[ix] > 1:
                self.brightest[ix] = (mag_list[ix] < mag_list[mat].min()+0.3) | (not np.isfinite(mag_list[ix]))
        
    def save_pickle(self):
        import cPickle as pickle
        print 'Save to %s' %(os.path.join(self.sproot, 'full_zspec.pkl'))
        fp =  open(os.path.join(self.sproot, 'full_zspec.pkl'), 'wb')
        pickle.dump(self, fp)
        fp.close()
    
    def read_pickle(self):
        import cPickle as pickle
        print 'Read from %s' %(os.path.join(self.sproot, 'full_zspec.pkl'))
        fp =  open(os.path.join(self.sproot, 'full_zspec.pkl'), 'rb')
        zsp = pickle.load(fp)
        fp.close()
        
        self.ra = zsp.ra
        self.dec = zsp.dec
        self.zspec = zsp.zspec
        self.catqual = zsp.catqual
        self.source = zsp.source
        self.catid = zsp.catid
        
    def read_catalogs(self):
        sproot = self.sproot
        #os.chdir(sproot)

        #### Make simple arrays, and *only include good redshifts*
        self.ra = []
        self.dec = []
        self.zspec = []
        self.catqual = []
        self.source = []
        self.catid = []

        ################### AEGIS
        #### DEEP2
        print '\nAEGIS: DEEP2\n'
        deep2 = pyfits.open(sproot+'DEEP2+3/zcat.dr3.v1_0.uniq.fits.gz')[1].data
        id = np.cast[float](deep2.field('OBJNO'))
        ra = np.cast[float](deep2.field('RA'))
        dec = np.cast[float](deep2.field('DEC'))
        dq = np.cast[float](deep2.field('ZQUALITY'))
        z_spec = np.cast[float](deep2.field('Z'))

        keep = dq >= 3
        N = keep.sum()

        self.ra.extend(ra[keep])
        self.dec.extend(dec[keep])
        self.zspec.extend(z_spec[keep])
        self.catqual.extend(dq[keep])
        self.source.extend(['DEEP2']*N)
        self.catid.extend(id[keep])

        #### Steidel
        print '\nSteidel 2003 LBGs\n'
        stei = unicorn.catalogs.read_steidel(verbose=False)
        N = len(stei['ra'])

        self.ra.extend(stei['ra'])
        self.dec.extend(stei['dec'])
        self.zspec.extend(stei['zuse'])
        self.catqual.extend([1]*N)
        self.source.extend(['Steidel03']*N)
        self.catid.extend(np.arange(N))

        ################### COSMOS
        #### zCOSMOS
        print '\nCOSMOS: zCOSMOS-bright\n'
        zc = catIO.Readfile(sproot+'/COSMOS/zCOSMOS_VIMOS_BRIGHT_DR2_TABLE.tab')
        keep = (np.floor(zc.cc) == 3) | (np.floor(zc.cc) == 4) | (zc.cc == 2.5) | (np.floor(zc.cc) == 13) | (np.floor(zc.cc) == 14) | (np.floor(zc.cc) == 23) | (np.floor(zc.cc) == 24)
        N = keep.sum()

        self.ra.extend(zc.ra[keep])
        self.dec.extend(zc.dec[keep])
        self.zspec.extend(zc.z[keep])
        self.catqual.extend(zc.cc[keep])
        self.source.extend(['zCOSMOS-Bright']*N)
        self.catid.extend(zc.name[keep])

        #### Brusa AGN compilation
        print '\nCOSMOS: Brusa\n'
        br = catIO.Readfile(sproot+'COSMOS/brusa_apj341700t2_mrt.txt')
        keep = (br.iflag == 1) & (br.zspec > 0) & (br.r_zspec != 5)
        N = keep.sum()

        self.ra.extend(br.ra[keep])
        self.dec.extend(br.dec[keep])
        self.zspec.extend(br.zspec[keep])
        self.catqual.extend(br.obj_class[keep])
        self.source.extend(['Brusa']*N)
        self.catid.extend(br.iid[keep])
        
        #### Onodera passive galaxies
        # http://iopscience.iop.org/0004-637X/755/1/26/fulltext/apj_755_1_26.html
        print '\nCOSMOS: Onodera\n'
        on = catIO.Readfile(sproot+'COSMOS/onodera.txt')
        keep = on.qual >= 3
        N = keep.sum()

        self.ra.extend(on.ra[keep])
        self.dec.extend(on.dec[keep])
        self.zspec.extend(on.zsp[keep])
        self.catqual.extend(on.qual[keep])
        self.source.extend(['Onodera']*N)
        self.catid.extend(on.id[keep])
        
        
        ################### GOODS-N
        #### DEEP3
        print '\nGOODS-N: DEEP3\n'
        d3 = pyfits.open(sproot+'DEEP2+3/cooper.2011apjs.goodsn.fits.gz')[1].data
        id = np.cast[float](d3.field('ID'))
        ra = np.cast[float](d3.field('RA'))
        dec = np.cast[float](d3.field('DEC'))
        dq = np.cast[float](d3.field('ZQUALITY'))
        z_spec = np.cast[float](d3.field('Z'))

        keep = dq >= 3
        N = keep.sum()

        self.ra.extend(ra[keep])
        self.dec.extend(dec[keep])
        self.zspec.extend(z_spec[keep])
        self.catqual.extend(dq[keep])
        self.source.extend(['DEEP3']*N)
        self.catid.extend(id[keep])
        
        #### Wirth (2004 TKRS)
        print '\nGOODS-N: Wirth TKRS\n'
        fp = open(sproot+'GOODS-North/wirth_TKRS_catalog.dat')
        lines = fp.readlines()
        fp.close()
        for line in lines:
            zsp = line[45:54]
            if zsp.strip() == '':
                continue
            #
            qual = line.split()
            lspl = line.split()
            qual = lspl[11]
            if int(qual) < 3:
                continue
            #
            #if 'J' not in lspl[-1]:
            #    continue
            #
            self.ra.append(threedhst.utils.DMS2decimal( ':'.join(lspl[1:4]), hours=True))
            self.dec.append(threedhst.utils.DMS2decimal( ':'.join(lspl[4:7]), hours=False))
            #self.ra.append(threedhst.utils.DMS2decimal(lspl[-1][1:10], hours=True))
            #self.dec.append(threedhst.utils.DMS2decimal(lspl[-1][11:19], hours=False))            
            self.zspec.append(np.float(lspl[10]))
            self.catqual.append(np.float(lspl[11]))
            self.source.append('WirthTKRS')
            self.catid.append(lspl[0])
        
        #### Yoshikawa (2010 ApJ, 718, 112) MOIRCS
        print '\nGOODS:N: Yoshikawa, z=2 (MOIRCS)\n'
        fp = open(sproot+'GOODS-North/apj338806t1_ascii.txt')
        lines = fp.readlines()
        fp.close()
        for line in lines:
            if not line.startswith('MODS'):
                continue
            #
            lspl = line.split()
            if lspl[14] != '...':
                self.ra.append(threedhst.utils.DMS2decimal( ':'.join(lspl[1:4]), hours=True))
                self.dec.append(threedhst.utils.DMS2decimal( ':'.join(lspl[4:7]), hours=False))
                self.zspec.append(np.float(lspl[14]))
                self.catqual.append(1)
                self.source.append('Yoshikawa')
                self.catid.append(lspl[0])
        
        #### Reddy (2006) LBGs, BM/BX
        print '\nGOODS:N: Reddy (2006)\n'
        fp = open(sproot+'GOODS-North/reddy_65294.tb2.txt')
        lines = fp.readlines()
        fp.close()
        for line in lines[1:]:
            lspl = line.split()
            zem = lspl[7]
            zabs = lspl[8]
            type = lspl[9]
            if type == '\ldots':
                continue
            #
            if (':' not in zem) & (':' not in zabs):
                if zem == '\ldots':
                    zsp = float(zabs)
                    type='abs_%s' %(type)
                else:
                    zsp = float(zem)
                    type='em_%s' %(type)
                #
                self.ra.append(threedhst.utils.DMS2decimal( ':'.join(lspl[1:4]), hours=True))
                self.dec.append(threedhst.utils.DMS2decimal( ':'.join(lspl[4:7]), hours=False))
                self.zspec.append(zsp)
                self.catqual.append(type)
                self.source.append('Reddy06')
                self.catid.append(lspl[0])
        
        #### Barger 2008
        print '\nGOODS-N: Barger\n'
        fp = open(sproot+'GOODS-North/table4.dat')
        lines = fp.readlines()
        fp.close()
        for line in lines:
            flag = line[65]
            if flag != 's':
                lspl = line.split()
                self.ra.append(np.float(lspl[0]))
                self.dec.append(np.float(lspl[1]))
                self.zspec.append(np.float(lspl[8]))
                self.catqual.append(1)
                self.source.append('Barger')
                self.catid.append('xx')
                
        ################### GOODS-S
        #### Fireworks
        print '\nGOODS-S: FIREWORKS\n'
        fw = catIO.Readfile(sproot+'GOODS-South/FIREWORKS_redshift.cat')
        fw_phot = catIO.Readfile(sproot+'GOODS-South/FIREWORKS_phot.cat')
        keep = (fw.zsp_qual == 1)
        N = keep.sum()
        
        self.ra.extend(fw_phot.ra[keep])
        self.dec.extend(fw_phot.dec[keep])
        self.zspec.extend(fw.zsp[keep])
        self.catqual.extend(fw.zsp_qual[keep])
        self.source.extend(['FIREWORKS']*N)
        self.catid.extend(fw.id[keep])

        #### Ballestra 2010
        # print '\nGOODS-S: Ballestra\n'
        # mr = catIO.Readfile(sproot+'GOODS-South/Ballestra_mrv2.dat')
        # keep = mr.zqual == 'A'
        # N = keep.sum()
        # 
        # self.ra.extend(mr.ra_vim[keep])
        # self.dec.extend(mr.dec_vim[keep])
        # self.zspec.extend(mr.zspec[keep])
        # self.catqual.extend(mr.zqual[keep])
        # self.source.extend(['Ballestra-mr']*N)
        # self.catid.extend(mr.source_id[keep])
        # 
        # lr = catIO.Readfile(sproot+'GOODS-South/Ballestra_lrv2.dat')
        # keep = lr.zqual == 'A'
        # N = keep.sum()
        # 
        # self.ra.extend(lr.ra_vim[keep])
        # self.dec.extend(lr.dec_vim[keep])
        # self.zspec.extend(lr.zspec[keep])
        # self.catqual.extend(lr.zqual[keep])
        # self.source.extend(['Ballestra-lr']*N)
        # self.catid.extend(lr.source_id[keep])
        
        #### ACES ECDFS survey, Cooper et al. 
        ## http://mur.ps.uci.edu/~cooper/ACES/zcatalog.html
        print '\nGOODS-S: ACES (Cooper)\n'
        aces = pyfits.open(sproot+'GOODS-South/zcat.ACES.2012jun04.fits.gz')[1].data
        keep = aces.ZQUALITY >= 3
        N = keep.sum()
        
        self.ra.extend(aces.RA[keep])
        self.dec.extend(aces.DEC[keep])
        self.zspec.extend(aces.Z[keep])
        self.catqual.extend(aces.ZQUALITY[keep])
        self.source.extend(['ACES']*N)
        self.catid.extend(aces.OBJNAME[keep])
        
        #### Full ESO compilation
        ## http://www.eso.org/sci/activities/garching/projects/goods/MasterSpectroscopy.html
        ## There is a lot of overlap with the Wuyts FIREWORKS catalog.
        print '\nGOODS-S: Full ESO compilation \n'
        eso = catIO.Readfile(sproot+'GOODS-South/MASTERCAT_v3.0.dat.fix')
        keep = eso.zspec < -100

        best_flags = [[0,'A'], [14,'A'], [1,'3'], [1,'4'], [2,'3.0'], [2,'2.0'], [8,'3'], [15,'1'], [16,'2'], [17,'2']]
        for flag in best_flags:
            keep = keep | ((eso.survey == flag[0]) & (eso.flag == flag[1]))
            
        no_flags = [3,4,5,6,7,10,11,13,18]
        for flag in no_flags:
            keep = keep | (eso.survey == flag)
        
        N = keep.sum()
        
        self.ra.extend(eso.ra[keep])
        self.dec.extend(eso.dec[keep])
        self.zspec.extend(eso.zspec[keep])
        self.catqual.extend(eso.flag[keep])
        ESO_ID = ['FORS2_vanzella','VVDS_v1', 'Szokoly', 'Croom', 'Dickinson', 'vderWel', 'Bunker', 'Stanway', 'LCIRS', 'x', 'Christiani', 'Strolger', 'x',  'Stanway04b', 'VIMOS', 'K20', 'IMAGES', 'Silverman', 'GMASS']
        for survey in eso.survey[keep]:
            self.source.append('ESO_%s' %(ESO_ID[survey]))

        self.catid.extend(eso.oid[keep])
        
        #### Stern (in prep) from Teplitz et al. 2011
        ### http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/141/1
        print '\nGOODS-N/S Stern (in prep) from Teplitz. 2011\n'
        stern = pyfits.open(sproot+'GOODS-South/teplitz_stern.fits')[1].data
        
        keep = stern['zsp'] > -0.1
        N = keep.sum()
        
        self.ra.extend(stern['RAJ2000'][keep])
        self.dec.extend(stern['DEJ2000'][keep])
        self.zspec.extend(stern['zsp'][keep])
        self.catqual.extend(np.array([1]*N))
        self.source.extend(['Stern']*N)
        self.catid.extend(stern['ID'][keep])
        
        ################### UDS
        #### UDS website compilation
        print '\nUDS\n'
        uds = pyfits.open(sproot+'UDS/UDS_redshifts_18Oct2010.fits')[1].data
        id = np.cast[str](uds.field('sxdfID'))
        zspec = np.cast[float](uds.field('z'))
        zq = np.cast[str](uds.field('quality'))
        zsou = np.cast[str](uds.field('specRun'))

        aki = zsou == 'Akiyama_FOCAS'
        keep = ((aki & (zq == 'A')) | (~aki)) & (zspec > 0)

        self.zspec.extend(zspec[keep])
        self.catqual.extend(list(zq[keep]))
        self.catid.extend(id[keep])

        idx = np.arange(len(zsou))
        ra_sxdf = np.cast[str](uds.field('sxdfRA'))
        dec_sxdf = np.cast[str](uds.field('sxdfDEC'))
        
        #### Convert sexagesimal coords to degrees
        for ii in idx[keep]:
            self.source.append('UDS-%s' %(zsou[ii]))
            self.ra.append(np.sum(np.cast[float](ra_sxdf[ii].split(':')) * np.array([1.,1./60,1./60**2]))*360./24)
            self.dec.append(-1*np.sum(np.abs(np.cast[float](dec_sxdf[ii].split(':'))) * np.array([1.,1./60,1./60**2])))
        
        ### Done Reading
        self.ra, self.dec, self.zspec = np.array(self.ra), np.array(self.dec), np.array(self.zspec)
        self.source = np.array(self.source)
        self.catqual = np.array(self.catqual)
        self.catid = np.array(self.catid)
        
def read_catalogs(field='COSMOS', force=False, v20=False):
    """ 
    Read all of the catalogs and run the matching between them.
    """
    import unicorn.catalogs2 as cat2
    
    #### Extract field from pointing ID if specified
    spl = field.split('-')
    if (len(spl) > 1) & (field[-1].isdigit()):
        field = '-'.join(spl[:-1])
        
    if (cat2.zout is not None) & (not force) & (field == cat2.field):
        print 'Looks like %s catalogs already read in.  To redo, use `read_catalogs(force=True)`.' %(field)
        return True
    
    PATH = '/3DHST/Spectra/Release/v2.1/'
    PATH = '/Users/brammer/3DHST/Spectra/Release/v2.1/'
    root = field.lower().replace('-','')
    
    cat = catIO.Readfile(PATH+'%s/%s_CATALOGS/%s_3dhst.v2.1.cat' %(field, field, root))
    zout = catIO.Readfile(PATH+'%s/%s_CATALOGS/%s_3dhst.v2.1.zout' %(field, field, root))
    fout = catIO.Readfile(PATH+'%s/%s_CATALOGS/%s_3dhst.v2.1.fout' %(field, field, root))
    try:
        rf = catIO.Readfile(PATH+'%s/%s_RF/%s.v2.1.rf.eazy.cat' %(field, field, root))
    except:
        rf = None
        
    #irf = catIO.Readfile(PATH+'%s/%s_RF/%s.v2.1.rf.interrest.rfout' %(field, field, root))

    zfit = catIO.Readfile(PATH+'%s/%s_CATALOGS/%s.zfit.linematched.dat' %(field, field, field))
    lines = pyfits.open(PATH+'%s/%s_CATALOGS/%s.linefit.fits' %(field, field, field))[1].data
    dq = catIO.Readfile(PATH+'%s/%s_CATALOGS/%s.dqflag.linematched.dat' %(field, field, field))
    
    ### Galfit files
    gfit_file = PATH+'%s/%s-GALFIT_v2.1/%s_F140W.cat' %(field, field, field.replace('-',''))
    if os.path.exists(gfit_file):
        gfit = catIO.Readfile(gfit_file)
    else:
        ### Empty catalog for fields that don't yet have GALFIT outputs
        columns = ['f', 'mag', 'dmag', 're', 'dre', 'n', 'dn', 'q', 'dq', 'pa', 'dpa', 'sn']
        gfit = catIO.EmptyCat()
        gfit['number'] = cat.id
        gfit['ra'] = cat.ra
        gfit['dec'] = cat.dec
        empty = cat.id*0-999
        for column in columns:
            gfit[column] = empty
        
        #gfit.parse()
        
    #### A few extra things    
    cat.m140 = 25-2.5*np.log10(cat.f_f140w)
    icolumn = {}
    icolumn['AEGIS'] = 'f_f814w'
    icolumn['COSMOS'] = 'f_f814w'
    icolumn['GOODS-N'] = 'f_f775w'
    icolumn['GOODS-S'] = 'f_f775w'
    icolumn['UDS'] = 'f_f814w'
    cat.imag = 25-2.5*np.log10(cat[icolumn[field]])
    
    cat.x_image, cat.y_image = cat.x, cat.y
    
    zfit.mag = dq.mag
    
    zfit.pointing = np.array([p.split('_')[0] for p in zfit.spec_id])
    
    ### Add some small scatter to the f_flagged numbers
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    
    cat2.field = field
    cat2.cat = cat
    cat2.zout = zout
    cat2.fout = fout
    cat2.zfit = zfit
    cat2.lines = lines
    cat2.dq = dq
    cat2.galfit = gfit
    cat2.root = root
    cat2.rf = rf
    cat2.irf = irf
    
    cat2.matcher = catIO.CoordinateMatcher(cat)
    
def read_catalogs_old(field='COSMOS', force=False, return_dir=True):
    """ 
    Read all of the catalogs and run the matching between them.
    """
    import unicorn.catalogs2 as cat2
    
    original_directory = os.getcwd()
    
    if (cat2.zout is not None) & (not force) & (field == cat2.field):
        print 'Looks like catalogs already read in.  To redo, use `read_catalogs(force=True)`.'
        return True
    
    if field == 'GOODS-N':
        if unicorn.hostname().startswith('uni'):
            os.chdir('/3DHST/Photometry/Release/v2.0/GOODS-N')
        else:
            os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-N')
            
        zfit = catIO.Readfile('GOODS-N.zfit.linematched.dat')
        dq = catIO.Readfile('GOODS-N.dqflag.linematched.dat')
        cat = catIO.Readfile('GOODS-N_v2.0_PHOTOMETRY/catalogs/goodsn_v1.8.fullz_wzp.cat')
        cat.x_image, cat.y_image = cat.ximage, cat.yimage
        zout = catIO.Readfile('GOODS-N_v2.0_PHOTOMETRY/goodsn_v1.8_eazy/goodsn_v1.8.fullz_wzp.zout')
        fout = catIO.Readfile('GOODS-N_v2.0_PHOTOMETRY/goodsn_v1.8_fast/goodsn_v1.8.fullz_wzp.fout')
        #### dq has some ****** IDL format values in max_contam
        lines = pyfits.open('GOODS-N.linefit.fits')
        #cat = catIO.Readfile('goodsn_f140w_v1_det_reform.cat')
        root = 'goodsn'
        cat.m140 = 25-2.5*np.log10(cat.f_f140)
    
    if field == 'COSMOS':
        ####  COSMOS
        if unicorn.hostname().startswith('uni'):
            ### still not working because v2.0 release directory is 
            ### [correctly] not writeable and Readfile can't make FITS file
            os.chdir('/3DHST/Photometry/Release/v2.0/COSMOS')
            gris_path = '/3DHST/Spectra/Release/v2.0/COSMOS/'
            zfit = catIO.Readfile(gris_path+'COSMOS.zfit.linematched.dat')
            dq = catIO.Readfile(gris_path+'COSMOS.dqflag.linematched.dat')
            lines = pyfits.open(gris_path+'COSMOS.linefit.fits')
            cat = catIO.Readfile('Catalog/3dhst.cosmos.v2.0.cat')
            fout = catIO.Readfile('Fast/3dhst.cosmos.v2.0.fout')
            zout = catIO.Readfile('Eazy/3dhst.cosmos.v2.0.zout')
        else:
            os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/COSMOS')
            zfit = catIO.Readfile('COSMOS.zfit.linematched.dat')
            dq = catIO.Readfile('COSMOS.dqflag.linematched.dat')
            lines = pyfits.open('COSMOS.linefit.fits')
            cat = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Catalog/3dhst.cosmos.v2.0.cat')
            fout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Fast/3dhst.cosmos.v2.0.fout')
            zout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Eazy/3dhst.cosmos.v2.0.zout')
        
        
        cat.x_image, cat.y_image = cat.x, cat.y
        root = 'cosmos'
        cat.m140 = 25-2.5*np.log10(cat.f140w)
        
    if field == 'GOODS-S':
        ####  GOODS-S
        if unicorn.hostname().startswith('uni'):
            os.chdir('/3DHST/Photometry/Release/v2.0/GOODS-S')
        else:
            os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/GOODS-S')
        
        zfit = catIO.Readfile('GOODS-S.zfit.linematched.dat')
        dq = catIO.Readfile('GOODS-S.dqflag.linematched.dat')
        lines = pyfits.open('GOODS-S.linefit.fits')
        cat = catIO.Readfile('GOODS-S_v2.0_PHOTOMETRY/catalog/GOODS-S_v2.0.fullz_wzp.cat')
        cat.x_image, cat.y_image = cat.ximage, cat.yimage
        fout = catIO.Readfile('GOODS-S_v2.0_PHOTOMETRY/GOODS-S_v2.0_fast/GOODS-S_v2.0.fullz_wzp.fout')
        zout = catIO.Readfile('GOODS-S_v2.0_PHOTOMETRY/GOODS-S_v2.0_eazy/GOODS-S_v2.0.fullz_wzp.zout')
        root = 'goodss'
        cat.m140 = 25-2.5*np.log10(cat.f_f140)
        
    #### A few extra things    
    zfit.mag = dq.mag
    norm = 1-np.random.normal(size=dq.f_flagged.shape)*0.05
    dq.f_flagged *= norm
    
    cat2.field = field
    cat2.cat = cat
    cat2.zout = zout
    cat2.fout = fout
    cat2.zfit = zfit
    cat2.lines = lines
    cat2.dq = dq
    cat2.root = root

    if return_dir:
        os.chdir(original_directory)
        
def view_selection(sel, OUTPUT='/tmp/selection.html', verbose=True, extra={}):
    import unicorn.catalogs2 as cat2
    
    scripts="http://localhost/~gbrammer/UDS/scripts"
    
    scripts="http:///localhost/~gbrammer/threedhst"
    
    threedhst.plotting.makeCSS(path=os.path.dirname(OUTPUT), title_size=18)
    
    PATH_TO_V2 = unicorn.GRISM_HOME+'../Release/v2.1/'
    
    if not os.path.exists(os.path.dirname(OUTPUT)):
        os.path.mkdir(os.path.dirname(OUTPUT))
    
    #
    extra_header = ''
    for key in extra.keys():
        extra_header += '<th> %s </th>' %(key)
        
    zrange = [np.min(cat2.zfit.z_max_spec[sel]), np.max(cat2.zfit.z_max_spec[sel])]
    maglim = np.max(cat2.cat.m140[sel])
    masslim = np.min(cat2.fout.lmass[sel])
    contam = -99 
    coverage = -99
    
    fplist = open(OUTPUT.replace('html','dat'),'w')
    fplist.write('# ID   ra   dec  z lmass\n')
    
    fpreg = open(OUTPUT.replace('html','reg'), 'w')
    fpreg.write('fk5\n')
    
    fp = open(OUTPUT,'w')
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
                        %d: {
                                sorter: false
                        },
                        %d: {
                                sorter: false
                        },
                        %d: {
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
        <th> m140 </th>
        <th> z </th>
        <th> logM </th>
        %s<th> RGB </th>
        <th> ZFIT </th>
        <th> 2D </th>
        <th> LINE </th>
    </thead>
    <tbody>
    """ %(scripts, scripts, scripts, 4+len(extra.keys()), 5+len(extra.keys()), 6+len(extra.keys()),
          maglim, masslim, zrange[0], zrange[1], contam, coverage, extra_header))
    
    NUSE = 0
    idx = np.arange(len(cat2.zfit.spec_id))[sel]
    
    for i in idx:            
        spec_id = cat2.zfit.spec_id[i]
        pointing = spec_id.split('_')[0]
        
        root_path = os.path.join(PATH_TO_V2, cat2.field, cat2.field+'-WFC3_v2.1_SPECTRA', pointing)
        rgb_path = os.path.join(PATH_TO_V2, cat2.field, cat2.field+'-WFC3_v2.1_SPECTRA/'+cat2.field+'_RGB/06.0')
        
        #print root_path
        
        fplist.write("%15s %14.6f %14.6f %8.3f %5.2f\n" %(spec_id, cat2.cat.ra[i], cat2.cat.dec[i], cat2.zfit.z_max_spec[i], cat2.fout.lmass[i]))
    
        ra_hex = threedhst.utils.decimal2HMS(cat2.cat.ra[i], hours=True)
        dec_hex = threedhst.utils.decimal2HMS(cat2.cat.dec[i], hours=False)
        coords = '%s %s' %(ra_hex, dec_hex)
        format='+'.join(coords.replace('+','%2B').replace('-','%2D').split())
        vizier="http://vizier.u-strasbg.fr/viz-bin/VizieR?-c=%s&-c.rs=%.1f" %(format, 1)
        
        fpreg.write('circle(%14.6f,%14.6f,1") # text={%s} color=magenta  \n' %(cat2.cat.ra[i], cat2.cat.dec[i], spec_id))

        fp.write("""
        <tr>
            <td> <a href="%s">%s</a><br> %13.6f  %13.6f <br> %s <br> %s </td>
            <td> %5.1f </td>
            <td> %5.2f </td>
            <td> %5.2f </td>
        """ %(vizier, spec_id, cat2.cat.ra[i], cat2.cat.dec[i],
              threedhst.utils.decimal2HMS(cat2.cat.ra[i], hours=True), 
              threedhst.utils.decimal2HMS(cat2.cat.dec[i], hours=False),
              cat2.cat.m140[i], cat2.zfit.z_max_spec[i],
                  cat2.fout.lmass[i]))
        
        for key in extra.keys():
            fp.write("            <td> %s </td>\n" %(extra[key][i]))
            
        fp.write("""            <td> <img src=%s/%s_rgb_06.0.png height=180px> </td>
            <td> <img src=%s/ZFIT/PNG/%s.zfit.png height=180px> </td>
            <td>  </td>
            <td> <img src=%s/LINE/PNG/%s.linefit.png height=180px> </td>
        </tr>""" %(rgb_path, spec_id, root_path, spec_id,
                   root_path, spec_id))

    
    fp.write("</tbody></table></body></html>")
    fp.close()
    fplist.close()
    fpreg.close()
    
    print 'N = %d' %(len(idx))
    if verbose:
        print OUTPUT
    if verbose > 1:
        os.system('open %s' %(OUTPUT))
        

        
def redshift_outliers():
    import unicorn.catalogs2 as cat2
    
    cat2.read_catalogs(field='COSMOS')

    cat2.read_catalogs(field='GOODS-S')
    
    #### Find cases within the photometric catalog that have nearby matches
    cat_coords = catIO.CoordinateMatcher(cat2.cat)
    dr_self, id_self = cat_coords.match_list(MATCH_SELF=True, N=1)
    dmag = cat2.cat.m140-cat2.cat.m140[id_self]
    
    has_brighter_neighbor =  (dr_self < 2) & (dmag > 0.7)
    neighbor_flag = (dr_self < 1) | has_brighter_neighbor
        
    zsp = cat2.SpeczCatalog()
    brusa = zsp.source[id_zsp] == 'Brusa'
    
    #### Match to spec catalog
    dr_zsp, id_zsp = zsp.match_list(ra=cat2.cat.ra, dec=cat2.cat.dec)
    has_spec = (dr_zsp < 1) & (cat2.zfit.z_max_spec > 0) & (zsp.zspec[id_zsp] > 0)
    ### zsp in original catalog
    cat_has_zsp = (cat2.zfit.z_spec > 0)
    
    dz = (cat2.zfit.z_max_spec-zsp.zspec[id_zsp]) / (1+zsp.zspec[id_zsp])
        
    plt.scatter(zsp.zspec[id_zsp][has_spec & ~brusa & cat_has_zsp], dz[has_spec & ~brusa & cat_has_zsp], color='blue', alpha=0.5)
    plt.scatter(zsp.zspec[id_zsp][has_spec & ~brusa & ~cat_has_zsp], dz[has_spec & ~brusa & ~cat_has_zsp], color='green', alpha=0.5)
    plt.scatter(zsp.zspec[id_zsp][has_spec & brusa], dz[has_spec & brusa], color='red', alpha=0.5)
    
    zp7 = zsp.zspec[id_zsp] > 0.7
    
    #### Outliers
    bad = np.abs(dz) > 0.05
    outliers = has_spec & bad & zp7 & ~neighbor_flag & ~brusa
    ids = cat2.zfit.spec_id[has_spec & outliers]
    bad_list = ''
    for id in ids:
        pointing = id.split('_')[0]
        bad_list += ' %s/ZFIT/PNG/%s.zfit*png' %(pointing, id)
        
    
    plt.scatter(cat2.zfit.z_spec[has_spec], cat2.zfit.z_max_spec[has_spec], alpha=0.5, color='blue')
    plt.scatter(zsp.zspec[id_zsp][has_spec], cat2.zfit.z_max_spec[has_spec], alpha=0.2, color='red')
    plt.xlim(0,2); plt.ylim(0,2)
    
    dz = (cat2.zfit.z_max_spec - cat2.zfit.z_spec) / (1+cat2.zfit.z_spec)
    
    plt.scatter(cat2.zfit.z_spec,  cat2.zfit.z_max_spec, alpha=0.5)
    plt.xlim(0,6); plt.ylim(0,6)
    
    sproot = unicorn.GRISM_HOME+'/ANALYSIS/SPECTROSCOPIC_REDSHIFTS/'
    print '\nGOODS-S: FIREWORKS\n'
    fw = catIO.Readfile(sproot+'GOODS-South/FIREWORKS_redshift.cat')
    fw_phot = catIO.Readfile(sproot+'GOODS-South/FIREWORKS_phot.cat')
    keep = (fw.zsp_qual == 1)
    N = keep.sum()
    
    cosd_gris = cat2.cat.ra*np.cos(cat2.cat.dec/360*2*np.pi)
    cosd_spec = fw_phot.ra*np.cos(fw_phot.dec/360*2*np.pi)
    
    k_c1 = -2.5*np.log10(cat2.cat.f_f160/cat2.cat.f_irac1)
    
    import scipy.spatial

    xy_gris = np.array([cosd_gris, cat2.cat.dec])
    xy_spec = np.array([cosd_spec, fw_phot.dec])
    tree_gris = scipy.spatial.cKDTree(xy_gris.T, 10)
    tree_spec = scipy.spatial.cKDTree(xy_spec.T, 10)
    
    zsp = copy.deepcopy(cat2.zout)
    zsp.dr_gris = np.zeros(zsp.N)-1
    zsp.mag_gris = np.zeros(zsp.N)-1
    zsp.idx_gris = np.zeros(zsp.N)-1
    zsp.dr_spec = np.zeros(zsp.N)-1
    zsp.z_spec = np.zeros(zsp.N)-1
    zsp.z_qual = np.zeros(zsp.N)-1
    
    for i in range(cat.N):
        print unicorn.noNewLine+'%d' %(i)
        point_gris = [cosd_gris[i], cat2.cat.dec[i]]
        dr_gris = tree_gris.query(point_gris, k=2)
        zsp.dr_gris[i] = dr_gris[0][1]*3600.
        zsp.idx_gris[i] = dr_gris[1][1]
        zsp.mag_gris[i] = cat2.cat.m140[dr_gris[1][1]]
        #
        dr_spec = tree_spec.query(point_gris, k=1)
        zsp.dr_spec[i] = dr_spec[0]*3600.
        zsp.z_spec[i] = fw.zsp[dr_spec[1]]
        zsp.z_qual[i] = fw.zsp_qual[dr_spec[1]]
    
    #
    dz = (cat2.zfit.z_max_spec - zsp.z_spec) / (1+zsp.z_spec)
    mat = (zsp.dr_spec < 1) & (zsp.z_spec > 0) & (zsp.z_qual > 0.9) & (zout.z_spec > 0) & (zsp.z_spec > 0.7) & (cat2.zfit.z_max_spec > 0)
    
    neighbor = (zsp.dr_gris < 2) & ((cat2.cat.m140-zsp.mag_gris) > 0)
    
    plt.scatter(zsp.z_spec[mat & neighbor], dz[mat & neighbor])
    
    #### make csv for web plot
    has_spec = (zsp.z_spec > 0) & (zsp.dr_spec < 1) & (cat2.cat.m140 < 23) & (cat2.zfit.z_max_spec > 0)
    data = np.array([cat2.zfit.phot_id[has_spec], zsp.z_spec[has_spec], cat2.zfit.z_max_spec[has_spec], cat2.cat.m140[has_spec], cat2.fout.lmass[has_spec]]).T
    np.savetxt(cat2.root+'.csv', data, fmt='%d,%.3f,%.3f,%.2f,%.2f')

def broad_band_eqw(flux, mag):
    """
    Use equation from the 3D-HST paper to compute equivalent widths based 
    on the emission line flux and the broad-band F140W magnitude.  
    
    Assumes:
    
    1) Flat continuum in Flam across F140W bandpass
    2) Line flux is the total line flux, in uits of 10^-17 erg/s/cm2
    """
    eqw = 8.78e-8 * (flux/5) / (10**(-0.4*mag)-2.21e-11*flux/5)
    return eqw
     
def high_eqw():
    
    import unicorn.catalogs2 as cat2
    import gbb
    
    cat2.read_catalogs(field='GOODS-S')
    
    zmax = cat2.zfit.z_max_spec
    
    sn_OIII_eqw = cat2.lines[1].data.OIII_EQW/cat2.lines[1].data.OIII_EQW_ERR
    sn_OII_eqw = cat2.lines[1].data.OII_EQW/cat2.lines[1].data.OII_EQW_ERR
    sn_HA_eqw = cat2.lines[1].data.HALPHA_EQW/cat2.lines[1].data.HALPHA_EQW_ERR
    sn_HB_eqw = cat2.lines[1].data.HBETA_EQW/cat2.lines[1].data.HBETA_EQW_ERR

    sn_OII = cat2.lines[1].data.OII_FLUX/cat2.lines[1].data.OII_FLUX_ERR
    sn_OIII = cat2.lines[1].data.OIII_FLUX/cat2.lines[1].data.OIII_FLUX_ERR
    sn_HA = cat2.lines[1].data.HALPHA_FLUX/cat2.lines[1].data.HALPHA_FLUX_ERR
    sn_HB = cat2.lines[1].data.HBETA_FLUX/cat2.lines[1].data.HBETA_FLUX_ERR
    
    OII_eqw = cat2.lines[1].data.OII_EQW/(1+zmax)
    OIII_eqw = cat2.lines[1].data.OIII_EQW/(1+zmax)
    HA_eqw = cat2.lines[1].data.HALPHA_EQW/(1+zmax)
    HB_eqw = cat2.lines[1].data.HBETA_EQW/(1+zmax)
    
    #### Compute eqw from line fluxes and magnitudes
    OII_mask = (zmax > 1.1e4/3727-1) & (zmax < 1.63e4/3727-1)
    OIII_mask = (zmax > 1.1e4/5007-1) & (zmax < 1.63e4/5007-1)
    HA_mask = (zmax > 1.1e4/6563.-1) & (zmax < 1.63e4/6563.-1)
    HB_mask = (zmax > 1.1e4/4861-1) & (zmax < 1.63e4/4861-1)
    
    flux_sum = cat2.lines[1].data.OIII_FLUX*OIII_mask + cat2.lines[1].data.HBETA_FLUX*HB_mask + cat2.lines[1].data.HALPHA_FLUX*HA_mask + cat2.lines[1].data.OII_FLUX*OII_mask
    
    mag_eqw = cat2.broad_band_eqw(flux_sum, cat2.cat.m140)/(1+zmax)
    mag_oiii_eqw = cat2.lines[1].data.OIII_FLUX/flux_sum*mag_eqw
    mag_halpha_eqw = cat2.lines[1].data.HALPHA_FLUX/flux_sum*mag_eqw
    mag_hbeta_eqw = cat2.lines[1].data.HBETA_FLUX/flux_sum*mag_eqw
    
    #### Compare measured and "mag" equivalent widths
    plt.plot(OIII_eqw[(sn_OIII > 2)], mag_oiii_eqw[(sn_OIII > 2)], marker='.', alpha=0.3, linestyle='None', label='OIII')
    plt.loglog(); plt.xlim(5,10000); plt.ylim(5,10000)
    plt.plot([1,10000],[1,1.e4])

    plt.plot(HA_eqw[(sn_HA > 2)], mag_halpha_eqw[(sn_HA > 2)], alpha=0.3, marker='.', linestyle='None', label=r'H$\alpha$')
    #plt.errorbar(HA_eqw[(sn_HA > 2)], mag_halpha_eqw[(sn_HA > 2)], xerr=cat2.lines[1].data.HALPHA_EQW_ERR[sn_HA > 2], alpha=0.3, marker='.', linestyle='None', label=r'H$\alpha$')
    plt.loglog(); plt.xlim(5,10000); plt.ylim(5,10000)
    plt.plot([1,10000],[1,1.e4])
    
    plt.legend(loc='upper left')
    
    #### THere is an offset in the measured and computed eqw, but
    #### probably color dependent since comiputed assumes flat continuum
    delta_HA = HA_eqw/mag_halpha_eqw
    delta_OIII = OIII_eqw/mag_oiii_eqw
    
    if cat2.field == 'GOODS-N':
        JH = -2.5*np.log10(cat2.cat.f_j/cat2.cat.f_h)
        JH = -2.5*np.log10(cat2.cat.f_acsz/cat2.cat.f_f140)
    #
    if cat2.field == 'GOODS-N':
        JH = -2.5*np.log10(cat2.cat.f_acsz/cat2.cat.f_f140)
    
    plt.plot(JH[(sn_HA > 2)], delta_HA[(sn_HA > 2)], alpha=0.5, marker='o', linestyle='None')
    plt.plot(JH[(sn_OIII > 2)], delta_OIII[(sn_OIII > 2)], alpha=0.5, marker='o', linestyle='None')
    plt.plot([-10,10000],[0,0])
    plt.semilogy(); plt.xlim(-1,1); plt.ylim(0.1,10)
    
    plt.plot(cat2.cat.m140[sn_OIII > 2], mag_oiii_eqw[sn_OIII > 2], 'o', linestyle='None', alpha=0.5)
    
    #### Select samples
    extra = {}
    extra['OIII'] = gbb.utils.formatted_array(mag_oiii_eqw,'%.f')
    extra['OIIIb'] = gbb.utils.formatted_array(OIII_eqw, '%.f')
    cat2.view_selection((mag_oiii_eqw > 500) & (sn_OIII > 2), extra=extra)
    
    mass = cat2.fout.lmass > 11.4
    cat2.view_selection(mass & (cat2.zfit.z_max_spec > 1.1), extra=extra)
    
    extra = {}
    extra[r'H&alpha;'] = mag_halpha_eqw
    cat2.view_selection((mag_halpha_eqw > 200) & (sn_HA > 3) & (zmax > 0.70), extra=extra)
    
    extra[r'H&alpha;'] = cat2.lines[1].data.HALPHA_FLUX
    cat2.view_selection((cat2.lines[1].data.HALPHA_FLUX > 50) & (zmax > 0.70), extra=extra)
    
    O2 = cat2.lines[1].data.OII_FLUX
    sn_O2 = cat2.lines[1].data.OII_FLUX/cat2.lines[1].data.OII_FLUX_ERR
    extra={}
    extra[r'OII'] = cat2.lines[1].data.OII_FLUX
    cat2.view_selection((sn_O2 > 3) & (zmax > 0.70), extra=extra)

    O3 = cat2.lines[1].data.OIII_FLUX
    sn_O3 = cat2.lines[1].data.OIII_FLUX/cat2.lines[1].data.OIII_FLUX_ERR
    extra={}
    extra[r'OIII'] = O3
    cat2.view_selection((sn_O3 > 20) & (zmax > 0.70), extra=extra)
    
    ### Metallicity
    O2 = cat2.lines[1].data.OII_FLUX
    O3 = cat2.lines[1].data.OIII_FLUX
    hb = cat2.lines[1].data.HBETA_FLUX
    sn_O2 = cat2.lines[1].data.OII_FLUX/cat2.lines[1].data.OII_FLUX_ERR
    sn_O3 = cat2.lines[1].data.OIII_FLUX/cat2.lines[1].data.OIII_FLUX_ERR
    sn_hb = cat2.lines[1].data.HBETA_FLUX/cat2.lines[1].data.HBETA_FLUX_ERR

    R23 = (O2+O3)/hb
    logOH_lo = gbb.utils.log_OH(o2=O2, o3=O3, hbeta=hb, track=-1)
    logOH_hi = gbb.utils.log_OH(o2=O2, o3=O3, hbeta=hb, track=1)
    O_ratio = O3/O2
    logOH = logOH_lo
    logOH[np.log10(O_ratio) < 0.45] = logOH_hi[np.log10(O_ratio) < 0.45]
    
    extra={}
    extra[r'Z'] = logOH
    extra[r'log R23'] = np.log10(R23)
    cat2.view_selection((sn_O3 > 3) & (sn_O2 > 2) & (sn_hb > 2), extra=extra)
    
    
    sel = (sn_O3 > 3) & (sn_O2 > 2) & (sn_hb > 2)
    plt.plot(cat2.fout.lmass[sel], logOH[sel], linestyle='None', marker='o', alpha=0.5)
    
    xwuyts = [9.5,8.4,9.9,9.0,9.2,9.2,8.7,8.6,9.4,9.9]
    ywuyts = [8.887,8.42, 8.11, 8.66, 8.37, 8.11, 7.67, 8.60, 8.65, 8.67]
    plt.scatter(xwuyts, ywuyts, color='red', s=30)
    
    # Maiolino 2008 polynomical fits to M-Z relation
    xm08 = np.arange(9,11.3,0.1)
    coeffs = [[11.18,9.04], [11.57, 9.04], [12.38, 8.99], [12.76, 8.79]]
    colors = ['black']*4
    mass_lim = [9,9,9.5,9.5]
    for i in range(4):
        m0, k0 = coeffs[i]
        ym08 = -0.0864*(xm08-m0)**2+k0
        plt.plot(xm08[xm08 > mass_lim[i]], ym08[xm08 > mass_lim[i]], color=colors[i], linewidth=1, alpha=0.8, zorder=1)
        if i == 2:
            plt.plot(xm08[xm08 > mass_lim[i]], ym08[xm08 > mass_lim[i]], color='orange', linewidth=2, alpha=0.8, zorder=-1)
    
    plt.xlabel(r'$\log M/M_\odot$')   
    plt.ylabel(r'12 + log O/H')

def cdfs_lowz_for_s3_lines():
    """
    Added [SIII] lines at ~9000A to the emission line templates.  Check line ratios.
    
    1) Check objects with H-alpha lines at z~0.7 where could also have [SIII]
    2) Check X-ray sources in CDFS with 0.22 < z_phot < 0.66

    """    
    import unicorn.catalogs2 as cat2
    cat2.read_catalogs(field='GOODS-S')
    
    sn = cat2.lines.SIII_FLUX/cat2.lines.SIII_FLUX_ERR > 2
    
    cat2.view_selection(sn)
    
    import astropy.io.ascii
    from astropy import coordinates as coord
    table = astropy.io.ascii.read('/3DHST/Ancillary/GOODS-S/CDFS/paper100-cdfs-xue-table3.txt')
    ra, dec = [], []
    for i in range(len(table)):
        rstr = '%0dh%02dm%02ds' %(table['RAh'][i], table['RAm'][i], table['RAs'][i])
        destr = '-%0dd%02dm%02s' %(table['DEd'][i], table['DEm'][i], table['DEs'][i]*0.99)
        rd = coord.ICRSCoordinates('%s %s' %(rstr, destr))
        ra.append(rd.ra.degrees)
        dec.append(rd.dec.degrees)
    
    ra, dec = np.array(ra), np.array(dec)
    m = catIO.CoordinateMatcher(cat2.cat)
    dr, ix = m.match_list(ra, dec)
    ok = (dr < 4) & (cat2.zfit.z_max_spec[ix] > 0)
    cat2.view_selection(ix[ok])
    
    ok = (cat2.zfit.z_max_spec > 0) & (cat2.zout.z_spec > 0.3) & (cat2.zout.z_spec < 0.7) & (cat2.cat.m140 < 21)
    cat2.view_selection(ok)
    
    obj = 'GOODS-S-23_28004'; ii = np.arange(cat2.cat.N)[cat2.zfit.spec_id == obj][0];  model = unicorn.reduce.process_GrismModel(cat2.zfit.pointing[ii]); x = catIO.CoordinateMatcher({'ra': model.ra_wcs, 'dec':model.dec_wcs}); dr, ix = x.find_nearest(cat2.cat.ra[ii], cat2.cat.dec[ii]); id = model.objects[ix]; 
    model.twod_spectrum(id, miny=25, USE_REFERENCE_THUMB=True); gris = unicorn.interlace_fit.GrismSpectrumFit('%s_%05d' %(model.root, id)); gris.fit_in_steps()
    
    obj = 'GOODS-S-23_28004'
    gris = unicorn.interlace_fit.GrismSpectrumFit(obj, RELEASE=True, lowz_thresh=0.05); gris.fit_in_steps(); os.system('open refit/%s*fit*png' %(obj))
    
    # GOODS-S-10_08509
    # GOODS-S-7_11762
    # GOODS-S-7_13537
    
def brown_dwarf():
    
    import unicorn.catalogs2 as cat2
    
    cat2.read_catalogs(field='GOODS-S')
    known = cat2.zfit.spec_id == 'GOODS-S-4_11332'

    cat2.read_catalogs(field='COSMOS')
    det = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Detection/F140W_psfmatched.cat')
    cat2.cat.flux_radius = det.flux_radius
    known = cat2.zfit.spec_id == 'GOODS-S-4_11332'
    xx = -2.5*np.log10(cat2.cat.v/cat2.cat.uvista_j)
    yy = -2.5*np.log10(cat2.cat.uvista_j/cat2.cat.uvista_ks)
    cat_coords = catIO.CoordinateMatcher(cat2.cat)
    dr, id = cat_coords.find_nearest(150.093115, 2.331402, N=2)
    known = id
    
    
    cat2.read_catalogs(field='GOODS-N')
    det = catIO.Readfile('goodsn_f140w_v1_det_reform.cat')
    cat2.cat.flux_radius = det.flux_radius
    cat_coords = catIO.CoordinateMatcher(cat2.cat)
    dr, id = cat_coords.find_nearest(189.006231, 62.170721, N=2)
    dr, id = cat_coords.find_nearest(189.223984, 62.188266, N=2)
    
    known = id
    
    has_spec = (cat2.zfit.z_max_spec > 0) 
    
    plt.scatter(cat2.cat.m140[has_spec], cat2.cat.flux_radius[has_spec], alpha=0.5)
    
    # 'GOODS-S', GOODS-N
    xx = -2.5*np.log10(cat2.cat.f_acsv/cat2.cat.f_j)
    yy = -2.5*np.log10(cat2.cat.f_j/cat2.cat.f_ks)
    
    compact = cat2.cat.flux_radius < 6
    
    plt.scatter(xx[compact & (cat2.cat.m140 < 26)], yy[compact & (cat2.cat.m140 < 26)], alpha=0.1, s=10, color='orange')
    plt.scatter(xx[has_spec & ~compact], yy[has_spec & ~compact], alpha=0.1, color='green')
    plt.scatter(xx[has_spec & compact], yy[has_spec & compact], alpha=0.5)
    
    plt.scatter(xx[known], yy[known], color='red', alpha=0.5, s=30)
    
    #plt.scatter(5.06020194, -0.65126636, s=30, color='red', alpha=0.5)
    
    selection = has_spec & (xx > 4.) & (yy < 0.1)
    cat2.zfit.spec_id[selection]
    
    ids = cat2.zfit.spec_id[selection]
    bad_list = ''
    for id in ids:
        pointing = id.split('_')[0]
        bad_list += ' %s/ZFIT/PNG/%s.zfit*png' %(pointing, id)
    
    ### Do the fits
    os.chdir('/Volumes/WD_1/3D-HST/v2.0/%s' %(cat2.field))
    try:
        os.mkdir('BD')
    except:
        pass
        
    os.chdir('BD')
    bd = unicorn.brown_dwarf.BD_fit()
    for id in ids:
        pointing = id.split('_')[0]
        print ''
        bd.fit('../%s/1D/ASCII/%s.1D.ascii' %(pointing, id), chi2_limit=100, max_contam=2, trim_mtype=False)
    
def match_v17_v20():
    
    import unicorn.catalogs2 as cat20
    import unicorn.catalogs as cat17
    
    field = 'GOODS-N'
    cat20.read_catalogs(field=field)
    
    cat17.read_catalogs()
    
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit, zsp
    
    keep17 = cat17.run_selection(zmin=0.7, zmax=2.3, fcontam=0.1, qzmin=0., qzmax=0.1, dr=1.0, has_zspec=True, fcovermin=0.9, fcovermax=1.0, massmin=8.5, massmax=15, magmin=0, magmax=23.8)
    
    sn_ha = (lines.halpha_eqw/lines.halpha_eqw_err)[lines.idx]
    #keep17 = keep17 & (sn_ha > 5) & (sn_ha < 7) & (phot.field[phot.idx] == field)
    
    ### Run the matcher on the full selection, then do second selection
    ### based on the v2.0 catalog
    dz17 = (zout.z_peak[::3][keep17]- cat20.zfit.z_spec[ix20])/(1+cat20.zfit.z_spec[ix20])
    dz20 = (cat20.zfit.z_max_spec[ix20]- cat20.zfit.z_spec[ix20])/(1+cat20.zfit.z_spec[ix20])
    ### better in v2.0
    keep17[keep17] = (cat20.zfit.z_spec[ix20] > 0) & (np.abs(dz17) > 0.03) & (np.abs(dz20) < 0.01) & (np.array(dr20) < 1)
    ### better in v1.7
    #keep17[keep17] = (cat20.zfit.z_spec[ix20] > 0) & (np.abs(dz20) > 0.03) & (np.abs(dz17) < 0.01) & (np.array(dr20) < 1)
    
    match20 = catIO.CoordinateMatcher(cat20.cat)
    
    idx_full = np.arange(len(keep17))
    dr20, ix20 = [], []
    for ix17 in idx_full[keep17]:
        pix17 = phot.idx[ix17]
        dr, ix = match20.find_nearest(phot.x_world[pix17], phot.y_world[pix17], N=1)
        dr20.append(dr)
        ix20.append(ix)
    
    plt.scatter(zout.z_peak[::3][keep17], cat20.zfit.z_max_spec[ix20], alpha=0.5)
    
    extra = {}
    
    z17 = np.zeros(cat20.zfit.N)
    z17[ix20] = zout.z_peak[::3][keep17]
    extra['z_v1.7'] = z17
    
    dz = np.zeros(cat20.zfit.N) 
    dz[ix20] = (zout.z_peak[::3][keep17]- cat20.zfit.z_max_spec[ix20])/(1+cat20.zfit.z_max_spec[ix20])
    extra['dz'] = gbb.utils.formatted_array(dz, fmt='%.3f')
    
    id_v17 = [' ']*cat20.zfit.N
    for i in range(len(ix20)):
        url = '<img src=http://unicorn.astro.yale.edu/P/GRISM_v1.7/EAZY/%s_eazy.png height=180px>' %(zout.id[::3][keep17][i])
        id_v17[ix20[i]] = url
    extra['ID v1.7'] = id_v17
    
    cat20.view_selection(ix20, extra=extra)
    
    phot.x_world[phot.idx]
    
    PATH_TO_V2 = unicorn.GRISM_HOME+'../Release/v2.0/'

    spec_id = cat2.zfit.spec_id[i]
    pointing = spec_id.split('_')[0]
    root_path = os.path.join(PATH_TO_V2, cat2.field, pointing)
    
def redshifts_v21(field='GOODS-S'):
    """
    Run some checks on the v2.1 catalog photo-zs
    """
    import unicorn.catalogs2 as cat2
    
    #field='GOODS-S'
    
    cat2.read_catalogs(field)

    zmax = {}
    zmax['AEGIS'] = 4
    zmax['COSMOS'] = 4
    zmax['GOODS-N'] = 7
    zmax['GOODS-S'] = 7
    zmax['UDS'] = 7

    #### Match to spec catalog
    zsp = cat2.SpeczCatalog(force_new=False)    
    dr_zsp, id_zsp = zsp.match_list(ra=cat2.cat.ra, dec=cat2.cat.dec)
    
    ### reject multiple matches
    print 'Flag multiple zsp matches'
    zsp.flag_multiples(cat2.cat.imag, toler=1.)
            
    #### Selection
    test = (dr_zsp < 1) & (cat2.zout.z_peak > 0) & (zsp.zspec[id_zsp] > 0) & (zsp.brightest == 1) #& (zsp.source[id_zsp] == 'ACES')
    
    #### GMASS
    #test = test & (zsp.source[id_zsp] == 'ESO_18')
    #plt.scatter(zsp.zspec[id_zsp][test], zout.z_peak[test], alpha=0.5)
    
    dz_phot = -(zsp.zspec[id_zsp]-cat2.zout.z_peak)/(1+zsp.zspec[id_zsp])
    dz_gris = -(zsp.zspec[id_zsp]-cat2.zfit.z_max_spec)/(1+zsp.zspec[id_zsp])
    threedhst.utils.nmad(dz_phot[test])
    threedhst.utils.nmad(dz_gris[test])
    
    #### Show the bad examples
    # bad = (np.abs(dz) > 0.15) & test & (cat2.zfit.z_max_spec > 0)
    extra = {}
    extra['dz'] = dz_gris
    extra['zsp'] = zsp.zspec[id_zsp]
    cat2.view_selection(bad, extra=extra)

    #### Colors for DQ flags
    colors = ['blue','green','red','orange','yellow','purple','cyan','magenta']
    
    surveys = np.unique(zsp.source[id_zsp][test])
    
    for survey in surveys:
        sub = test & (zsp.source[id_zsp] == survey)
        if sub.sum() > 0:
            quals = np.unique(zsp.catqual[id_zsp][sub])
            print survey
            fig = unicorn.plotting.plot_init(square=True, xs=9, aspect=1.2, left=0.15, bottom=0.1, top=0.08, hspace=0.0)
            ### Make plots for different redshift measurements
            for subplot, redshift in zip([221,222], [cat2.zout.z_peak, cat2.zfit.z_max_spec]):
                sub_z = test & (zsp.source[id_zsp] == survey) & (redshift > 0)
                dz = (zsp.zspec[id_zsp]-redshift)/(1+zsp.zspec[id_zsp])
                #
                #### Difference, large scale
                ax = fig.add_subplot(subplot+200)
                for iq, qual in enumerate(quals):
                    sub_qual = sub_z & (zsp.catqual[id_zsp] == qual)
                    cat2.zphot_zspec(zsp.zspec[id_zsp][sub_qual], redshift[sub_qual], alpha=0.4, label='%s %4d %6.4f' %(str(qual), sub_qual.sum(), threedhst.utils.nmad(dz[sub_qual])), ms=3, marker='s', verbose=1, ax=ax, zmax=zmax[field], diff = True, color=colors[iq])
                #
                p = ax.plot([0,np.log10(1+zmax[field])], [0,0], color='red', alpha=0.5)
                p = ax.set_ylim([-0.9,0.9]); p = ax.set_xlabel(''); p = ax.set_xticklabels([])
                #
                p = plt.title('%s, %s (%d) $\delta z$=%0.4f' %(field.lower().replace('-',''), survey, sub_z.sum(), threedhst.utils.nmad(dz[sub_z])))
                ####
                #### Difference, small scale
                ax = fig.add_subplot(subplot+202)
                for iq, qual in enumerate(quals):
                    sub_qual = sub_z & (zsp.catqual[id_zsp] == qual)
                    cat2.zphot_zspec(zsp.zspec[id_zsp][sub_qual], redshift[sub_qual], alpha=0.4, label='%s %4d %6.4f' %(str(qual), sub_qual.sum(), threedhst.utils.nmad(dz[sub_qual])), ms=3, marker='s', verbose=1, ax=ax, zmax=zmax[field], diff = True, color=colors[iq])
                #
                p = ax.plot([0,np.log10(1+zmax[field])], [0,0], color='red', alpha=0.5)
                p = ax.set_ylim([-0.038,0.038]); p = ax.set_xlabel(''); p = ax.set_xticklabels([])
                #
                #### One-to-one
                ax = fig.add_subplot(subplot+2)
                for iq, qual in enumerate(quals):
                    sub_qual = sub & (zsp.catqual[id_zsp] == qual) & (redshift > 0)
                    cat2.zphot_zspec(zsp.zspec[id_zsp][sub_qual], redshift[sub_qual], alpha=0.4, label='%s %4d %6.4f %4d %6.4f' %(str(qual), sub_qual.sum(), threedhst.utils.nmad(dz[sub_qual]), (sub_qual & (redshift > 0.7)).sum(), threedhst.utils.nmad(dz[sub_qual & (redshift > 0.7)])), ms=3, marker='s', verbose=1, ax=ax, zmax=zmax[field], diff = False, color=colors[iq])
                #
                p = ax.legend(prop=dict(size=8))
                ax.plot([0,np.log10(1+zmax[field])], [0,np.log10(1+zmax[field])], color='red', alpha=0.5)
                #
            #
            p = plt.savefig('%s_%s.png' %(field.lower().replace('-',''), survey.strip()))
            p = plt.close()
    #
    
    gris = test & (cat2.zfit.z_max_spec > 0)
    gris = gris & (zsp.zspec[id_zsp] > 0.68)
    
    fig = cat2.zphot_zspec(zsp.zspec[id_zsp][gris], cat2.zout.z_peak[gris], alpha=0.2, label='z_peak', ms=3, color='black', marker='s', verbose=2)
    fig.savefig(field.lower().replace('-','')+'_zpeak.png')
    plt.close()
    
    fig = cat2.zphot_zspec(zsp.zspec[id_zsp][gris], cat2.zfit.z_max_spec[gris], alpha=0.2, label='z_max_spec', ms=3, color='black', marker='s', verbose=2)
    fig.savefig(field.lower().replace('-','')+'_z_max_spec.png')
    plt.close()
    
def zphot_zspec(zspec, zphot, zmax=6, ax=None, marker='o', ylabel='phot', verbose=True, diff=False, *args, **kwargs):
    """
    Make a zphot_zspec plot with "pseudo-log" scaling
    """
    ret = False
    if ax is None:
        fig = unicorn.plotting.plot_init(square=True, xs=5, aspect=1, left=0.1)
        ax = fig.add_subplot(111)
        ret = True
    
    if diff:
        ydata = (zphot-zspec) / (1+zspec)
    else:
        ydata = np.log10(1+zphot)
        
    ax.plot(np.log10(1+zspec), ydata, linestyle='None', marker=marker, *args, **kwargs)
    
    xlab = np.arange(0,zmax+0.1,0.5)
    tlab = np.array(['%d' %(x) for x in xlab])
    tlab[1::2] = ''
    ax.set_xticks(np.log10(1+xlab))
    ax.set_xticklabels(tlab)
    ax.set_xlim(0, np.log10(1+zmax))
    ax.set_xlabel(r'$z_\mathrm{spec}$')
    
    if diff:
        ax.set_ylabel(r'$\Delta z/(1+z)$')
        zarr = np.arange(0,zmax,0.1)
        for ztest in range(zmax):
            ax.plot(np.log10(1+zarr), (ztest-zarr)/(1+zarr), color='0.8', alpha=0.05)
            
    else:
        ax.set_yticks(np.log10(1+xlab))
        ax.set_ylim(0, np.log10(1+zmax))
        ax.set_yticklabels(tlab) 
        ax.set_ylabel(r'$z_\mathrm{%s}$' %(ylabel))
        
    if verbose:
        dz = (zphot-zspec)/(1+zspec)
        if verbose > 1:
            ax.text(0.5,0.95,r'$N=%d$, $\delta z=%.4f$' %(len(dz), threedhst.utils.nmad(dz)), ha='center', va='top', transform=ax.transAxes)
        else:
            print 'N=%d\n NMAD(dz)=%.4f' %(len(dz), threedhst.utils.nmad(dz))
            
    if ret:
        ax.plot([0,np.log10(1+zmax)],[0,np.log10(1+zmax)], color='red', alpha=0.5)
        return fig
     
def check_rf_colors():
    """
    Kate noticed something with the rest-frame colors in that there is
    a cloud of apparently quiescent galaxies that scatter to redder U-V
    colors in InterRest compared to EAZY.
    """   
    
    ids = ['AEGIS-26_00339', 
    'AEGIS-6_02236', 
    'AEGIS-6_03994', 
    'AEGIS-8_06138', 
    'AEGIS-21_03851',
    'AEGIS-9_05683', 
    'AEGIS-14_06299',
    'AEGIS-27_12031', 
    'COSMOS-25_21952',
    'COSMOS-11_29935', 
    'COSMOS-12_09739'] 
    
    os.chdir('/Volumes/Crucial/3DHST/Spectra/Release/v2.1/Tests')
    for id in ids:
        gris = unicorn.interlace_fit.GrismSpectrumFit(id, RELEASE=True)
        cat2.show_rf_color(gris)
        plt.savefig(id+'_rf.png')
        plt.close()
        
def show_rf_color(gris):
    """
    gris is a unicorn.interlace_fit.GrismSpectrumFit object.
    """
    import unicorn.catalogs2 as cat2
    cat2.read_catalogs(gris.field)
    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = gris.eazy_fit
    
    plt.plot(lambdaz, temp_sed, color='blue', alpha=0.8)
    plt.semilogx()
    plt.xlim(3000, 9.e4)
    
    ok = efobs > 0
    plt.errorbar(lci[ok], fobs[ok], efobs[ok], linestyle='None', color='black', ecolor='black', alpha=.8, marker='o')
    
    irf_u = cat2.irf.rf_f6[gris.ix]
    irf_v = cat2.irf.rf_f8[gris.ix]
    
    rf_u = cat2.rf.l135[gris.ix]
    rf_v = cat2.rf.l137[gris.ix]
    
    irf_uv = -2.5*np.log10(irf_u/irf_v)
    rf_uv = -2.5*np.log10(rf_u/rf_v)
    
    red = (1+cat2.irf.redshift[gris.ix])
    lc = np.array([3.586e3, 5.4776e3])
    flam_factor = 10**(-0.4*(25+48.6))*3.e18/1.e-17/(lc*red)**2
    
    plt.plot(lc*red, np.array([rf_u, rf_v])*flam_factor, color='blue', linewidth=2, marker='o', ms=10, alpha=0.8, label='Eazy: U-V=%.2f' %(rf_uv))
    plt.plot(lc*red, np.array([irf_u, irf_v])*flam_factor, color='red', linewidth=2, marker='o', ms=10, alpha=0.8, label='Irest: U-V=%.2f' %(irf_uv))
    plt.legend(prop=dict(size=10))
    
    plt.ylim((fobs[ok]-3*efobs[ok]).min(), (fobs[ok]+3*efobs[ok]).max())
    plt.title(gris.grism_id)
    
class SpecObject():
    """
    Data holder for a general 3D-HST object.  Store spectrum values
    as well as photometric catalog information
    """
    def __init__(self, spec_id='AEGIS-26_00339', version=2.1, path=None):
        import unicorn
        import unicorn.catalogs2 as cat2
        
        if path is None:
            path = '/3DHST/Spectra/Release/v2.1/'
        
        self.spec_id = spec_id
        self.field = '-'.join(spec_id.split('-')[:-1])
        self.pointing = spec_id.split('_')[0]
        id = int(spec_id.split('_')[-1])
        self.ix = np.where(cat2.cat.id == id)[0][0]
        
        self.spec_path = path+'%s/%s-WFC3_v2.1_SPECTRA/%s/' %(self.field, self.field, self.pointing)
        
        self.gris = unicorn.interlace_fit.GrismSpectrumFit(self.spec_id, RELEASE=True)
        
        # lambdaz, temp_sed, lci, obs_sed, fobs, efobs = self.eazy_fit
        
def extract_wavelet_lines():
    import unicorn.catalogs2 as cat2
    
    mag = 25-2.5*np.log10(cat2.cat.f_f140w)
    keep = mag < 25.8
    
    try:
        os.chdir('/3DHST/Spectra/Release/v2.1/%s/%s_WAVELET/' %(cat2.field, cat2.field))
    except:
        os.mkdir('/3DHST/Spectra/Release/v2.1/%s/%s_WAVELET/' %(cat2.field, cat2.field))
        os.chdir('/3DHST/Spectra/Release/v2.1/%s/%s_WAVELET/' %(cat2.field, cat2.field))
    
    BASE_PATH = '/3DHST/Spectra/Release/v2.1/%s/%s-WFC3_v2.1_SPECTRA/' %(cat2.field, cat2.field)
    
    for i in np.arange(cat2.cat.N)[keep]:
        print unicorn.noNewLine+'%d' %(cat2.cat.id[i])
        pointings = threedhst.utils.which_3dhst_pointing(cat2.cat.ra[i], cat2.cat.dec[i])
        for point in pointings:
            if point.startswith('PRIM'):
                continue
            #
            p = point.replace('SOUTH','S').replace('GOODSN','GOODS-N')
            object = '%s_%05d' %(p, cat2.cat.id[i])
            oned_file = '%s/%s/1D/FITS/%s.1D.fits' %(BASE_PATH, p, object)
            if not os.path.exists(oned_file):
                continue
            #
            oned = unicorn.reduce.Interlace1D(oned_file, PNG=False)
            lines = oned.find_em_lines()
            #### Check for UDFj redshift
            for i in range(len(lines['wave'])):
                if (np.abs(lines['wave'][i] - 1.599e4) < 100) & (lines['sn'][i] > 3):
                    print oned_file.replace('FITS','PNG').replace('fits','png')
                    
    l0 = 1.599e4
    files=glob.glob('*wavelet.dat')
    for file in files:
        spl = open(file).readlines()[1].split()
        waves = np.cast[float](spl[4::6])
        sn = np.cast[float](spl[7::6])
        count = 0
        for i in range(len(waves)):
            if (waves[i] > 1.2e4) & (sn[i] > 1):
                count = count + 1
                # if (np.abs(waves[i] - l0) < 50) & (sn[i] > 1):
            #     print spl[0], spl[3], waves[i], sn[i]
        if (float(spl[3]) > 24) & (count > 1):
            print spl[0], spl[3], waves, sn   
        
    
        
def make_web_browser(pointing='UDS-17'):
    """
    Modify the gmap webpages to pull up the nearest object with a spectrum + fit.
    """
    import json
    import unicorn.catalogs2 as cat2
    
    cat2.read_catalogs(pointing)
    #this = cat2.zfit.pointing == pointing
    this = cat2.zfit.pointing != '00000'
    
    data_list = []
    idx = np.arange(cat2.zfit.N)[this]
    for i in idx:
        data = {}
        data['id'] = cat2.zfit.spec_id[i]
        data['ra'] = float('%.6f' %(cat2.cat.ra[i]))
        data['dec'] = float('%.6f' %(cat2.cat.dec[i]))
        data['z_gris'] = float('%.3f' %(cat2.zfit.z_max_spec[i]))
        data['z_spec'] = float('%.3f' %(cat2.zfit.z_spec[i]))
        data['mass'] = float('%.3f' %(cat2.fout.lmass[i]))
        data['m140'] = float('%.3f' %(cat2.cat.m140[i]))
        #data['r50'] = float('%.3f' %(cat2.cat.f160W_flux_radius[i]))
        data_list.append(data)
        
    for dir in ['HTML','HTML/tiles','HTML/scripts']:
        try:
            os.mkdir(dir)
        except:
            pass
        
    fp = open('HTML/%s.json' %(pointing),'w')
    json.dump(data_list, fp)
    fp.close()
    
    line = open('HTML/%s.json' %(pointing)).readlines()[0].replace('NaN', '-99')
    fp = open('HTML/%s.json' %(pointing),'w')
    fp.write(line)
    fp.close()
    
    # threedhst.gmap.makeImageMap(['%s-F140W_drz.fits' %(pointing),'%s-G141_drz.fits[1]*4' %(pointing)], aper_list=[15,16], zmin=-0.15, zmax=1.5, path='./HTML/', tileroot=['F140W','G141'], invert=True)
    
    #### Insert these bits into the HTML file
    ## <!-- xxx -->
    script = """
    <script type="text/javascript"> 
    
    var objects = [];
    var Nobj = 0;
    $.getJSON("%s.json", null, 
        function(json) {
            $.each(json, function(i, item){
              objects.push(item);
              Nobj += 1;
            });
        }
    );
    
    ///// Add listener for "s" key to open/close the spectrum popup
    var ShowSpectrum = 0;
    $(document).keypress(function(event) {
        // alert('charCode: '+event.charCode);
        // 's'
        if (event.charCode == 115) {
            //alert('Next!');
            if (ShowSpectrum == 1) {
                ShowSpectrum = 0;
                $('#spectrum').hide();
            } else {
                nearest_3dhst();
                ShowSpectrum = 1;
            }
        };
    });
    
    function nearest_3dhst() {
        var mapcenter = map.getCenter();
        var dec = mapcenter.lat()+centerLat;                       
        var ra = ((360-mapcenter.lng()/Math.cos(centerLat/360.*2*3.14159)+offset-centerLng));
        var match = 0;
        var matched = [];
        var drmin = 1000.;
        $.each(objects, function(i, object) {
            var dr = Math.sqrt(Math.pow((ra-object.ra)*Math.cos(dec/360.*2*3.14159),2) + Math.pow(dec-object.dec,2))*3600.;
            //alert(objects[i].id + ',' + dr + ',' + drmin);
            if ( dr < drmin ) {
                drmin = dr*1.;
                match = i*1;
            };
            if ( dr < 3 ) { 
                matched.push(i);
            }
        });
        $("#spectrum").show();
        var popup = "<img id=\\"zfit\\" src=\\"../ZFIT/PNG/" + objects[match].id + ".zfit.png\\" /> <br>";
        popup += "<img src=../RGB/" + objects[match].id + "_rgb_03.0.png onMouseOver=\\"$('#zfit').attr('src', '../ZFIT/PNG/" + objects[match].id + ".zfit.png');\\" margin-left='10px'/>";
        $.each(matched, function(i, im) {
            if ( im != match ) {
                popup += "<img src=../RGB/" + objects[im].id + "_rgb_03.0.png onMouseOver=\\"$('#zfit').attr('src', '../ZFIT/PNG/" + objects[im].id + ".zfit.png');\\" margin-left='10px'/>";
            }
        });
        $("#spectrum").html(popup);
        return [match, drmin];
    }; 
    
    </script>
    
    <!-- xxx -->
    
    """ %(pointing)
    
    ### </div> <!-- content -->
    spectrum_div = """
    <div id="spectrum" onClick="$(this).hide();"> </div>
    
    </div> <!-- content -->
    """
    
    lines = open('HTML/map.html').readlines()
    for i, line in enumerate(lines):
        if '<!-- xxx -->' in line:
            lines[i] = script
        if '<!-- content -->' in line:
            lines[i] = spectrum_div
    
    fp = open('HTML/map.html','w')
    fp.writelines(lines)
    fp.close()
     
    ### style.css / header    
    css = """
    #spectrum {
        max-width:90%;
        display:none;
        position:absolute;
        top:10px;
        left:10px;
        border:1px solid red;
        background: white;
    }
    """
    
    lines = open('HTML/scripts/style.css').readlines()
    if 'spectrum' not in lines[1]:
        lines[0] = css
    
    fp = open('HTML/scripts/style.css','w')
    fp.writelines(lines)
    fp.close()
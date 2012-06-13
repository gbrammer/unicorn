"""
Work with v2.0 catalogs
"""

import os
import pyfits
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

import cosmocalc

PATH_TO_CAT = unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/'

try:
    zout
except: 
    zout = None
    fout = None
    cat = None
    zfit = None
    lines = None
    dq = None
    selection_params = None
    field = None

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
        
        return dr_match, id_match
    
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
        N = len(keep[keep])

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
        N = len(keep[keep])

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
        N = len(keep[keep])

        self.ra.extend(br.ra[keep])
        self.dec.extend(br.dec[keep])
        self.zspec.extend(br.zspec[keep])
        self.catqual.extend(br.obj_class[keep])
        self.source.extend(['Brusa']*N)
        self.catid.extend(br.iid[keep])

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
        N = len(keep[keep])

        self.ra.extend(ra[keep])
        self.dec.extend(dec[keep])
        self.zspec.extend(z_spec[keep])
        self.catqual.extend(dq[keep])
        self.source.extend(['DEEP2']*N)
        self.catid.extend(id[keep])

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
        N = len(keep[keep])

        self.ra.extend(fw_phot.ra[keep])
        self.dec.extend(fw_phot.dec[keep])
        self.zspec.extend(fw.zsp[keep])
        self.catqual.extend(fw.zsp_qual[keep])
        self.source.extend(['FIREWORKS']*N)
        self.catid.extend(fw.id[keep])

        #### Ballestra 2010
        print '\nGOODS-S: Ballestra\n'
        mr = catIO.Readfile(sproot+'GOODS-South/Ballestra_mrv2.dat')
        keep = mr.zqual == 'A'
        N = len(keep[keep])

        self.ra.extend(mr.ra_vim[keep])
        self.dec.extend(mr.dec_vim[keep])
        self.zspec.extend(mr.zspec[keep])
        self.catqual.extend(mr.zqual[keep])
        self.source.extend(['Ballestra-mr']*N)
        self.catid.extend(mr.source_id[keep])

        lr = catIO.Readfile(sproot+'GOODS-South/Ballestra_lrv2.dat')
        keep = lr.zqual == 'A'
        N = len(keep[keep])

        self.ra.extend(lr.ra_vim[keep])
        self.dec.extend(lr.dec_vim[keep])
        self.zspec.extend(lr.zspec[keep])
        self.catqual.extend(lr.zqual[keep])
        self.source.extend(['Ballestra-lr']*N)
        self.catid.extend(lr.source_id[keep])

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
        
        
def read_catalogs(field='COSMOS', force=False, return_dir=True):
    """ 
    Read all of the catalogs and run the matching between them.
    """
    import unicorn.catalogs2 as cat2
    
    original_directory = os.getcwd()
    
    if (cat2.zout is not None) & (not force) & (field == cat2.field):
        print 'Looks like catalogs already read in.  To redo, use `read_catalogs(force=True)`.'
        return True
    
    if field == 'GOODS-N':
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
        os.chdir('/research/HST/GRISM/3DHST/RELEASE_v2.0/COSMOS')
        zfit = catIO.Readfile('COSMOS.zfit.linematched.dat')
        dq = catIO.Readfile('COSMOS.dqflag.linematched.dat')
        lines = pyfits.open('COSMOS.linefit.fits')
        cat = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Catalog/3dhst.cosmos.v2.0.cat')
        
        cat.x_image, cat.y_image = cat.x, cat.y
        root = 'cosmos'
        fout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Fast/3dhst.cosmos.v2.0.fout')
        zout = catIO.Readfile('COSMOS_v2.0_PHOTOMETRY/Eazy/3dhst.cosmos.v2.0.zout')
        cat.m140 = 25-2.5*np.log10(cat.f140w)
        
    if field == 'GOODS-S':
        ####  GOODS-S
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
    
    PATH_TO_V2 = unicorn.GRISM_HOME+'../Release/v2.0/'
    
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
        %s<th> ZFIT </th>
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
        
        root_path = os.path.join(PATH_TO_V2, cat2.field, pointing)
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
            <td> <a href="%s">%s</a><br> %13.6f <br> %13.6f </td>
            <td> %5.1f </td>
            <td> %5.2f </td>
            <td> %5.2f </td>
        """ %(vizier, spec_id, cat2.cat.ra[i],
                  cat2.cat.dec[i], cat2.cat.m140[i], cat2.zfit.z_max_spec[i],
                  cat2.fout.lmass[i]))
        
        for key in extra.keys():
            fp.write("            <td> %s </td>\n" %(extra[key][i]))
            
        fp.write("""            <td> <img src=%s/ZFIT/PNG/%s.zfit.png height=180px> </td>
            <td>  </td>
            <td> <img src=%s/LINE/PNG/%s.linefit.png height=180px> </td>
        </tr>""" %(root_path, spec_id,
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
    
    brusa = zsp.source[id_zsp] == 'Brusa'
    
    zsp = cat2.SpeczCatalog()
    
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
    N = len(keep[keep])
    
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
    
    cat2.read_catalogs(field='COSMOS')
    
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
    
    
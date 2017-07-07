"""
Processing scripts for version 4.1.4
"""

def make_zbest_catalogs():
    """
    Make z_best = z_grism / z_peak catalogs
    """
    from threedhst import catIO
    
    os.chdir('/Users/brammer/3DHST/Spectra/Release/v4.1.4/RestFrameColors')
    
    PATH_TO_NEW = '/Users/brammer/3DHST/Spectra/Release/v4.1.4'
    PATH_TO_PHOT = '/Users/brammer/3DHST/Spectra/Release/v4.1'
    
    PATH_TO_NEW = '/Users/brammer/3DHST/Spectra/Release/v4.1.5/RF_Colors'
    
    filters = []
    filters.extend(['153', '154', '155']) # Maiz UBV
    filters.extend(['161', '162', '163']) # 2MASS JHK
    filters.extend(['156', '157', '158', '159', '160']) # SDSS ugriz
    filters.extend(['135','136','137','138','139']) # Bessel UX BVRI
    filters.extend(['270','271','272','273','274','275']) # Tophat filters
    
    #### Add z_gris to catalog column to use as z_spec
    field='aegis'
    
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        print field
        
        old_cat = catIO.Table('%s/%s_3dhst.v4.1.cats/Catalog/%s_3dhst.v4.1.cat' %(PATH_TO_PHOT, field, field))
    
        #### v4.1.4
        # old_zout = catIO.Table('%s/%s_3dhst.v4.1.cats/Eazy/%s_3dhst.v4.1.zout' %(PATH_TO_PHOT, field, field))
        # zfit = catIO.Table('%s/%s-WFC3_v4.1.4_SPECTRA/%s-catalogs/%s.new_zfit.linematched.v4.1.4.dat' %(PATH_TO_NEW, field.upper(), field, field))
        #     
        # z_best = old_zout['z_peak']
        # z_type = np.ones(len(z_best), dtype=int)
        # has_grism = zfit['z_max_grism'] > 0
        # z_best[has_grism] = zfit['z_max_grism'][has_grism]
        # z_type[has_grism] = 2
        
        #### v4.1.5
        # zfit = catIO.Table('%s/%s.zbest.v4.1.5.fits' %(PATH_TO_NEW, field))
        # z_best = zfit['z_best']
        # z_type = zfit['z_best_s']
        # 
        # z_best_column = catIO.table_base.Column(data=z_best, name='z_best')
        # old_cat.add_column(z_best_column)
        # 
        # z_type_column = catIO.table_base.Column(data=z_type, name='z_type')
        # old_cat.add_column(z_type_column)
        # old_cat.write(os.path.basename(old_cat.filename).replace('.cat', '.zbest.cat'), format='ascii.commented_header')
    
        #### v4.1.5 zbest + zgrism
        zfit = catIO.Table('%s/%s_3dhst.zfit.linematched.v4.1.5.fits' %(PATH_TO_NEW, field))
        old_cat.add_column(zfit['z_max_grism'])
        old_cat.add_column(zfit['z_best'])
        old_cat.add_column(zfit['z_best_s'])
        
        old_cat.write(os.path.basename(old_cat.filename).replace('.cat', '.zfit.cat'), format='ascii.commented_header')
        
    
    
def go():
    import unicorn.v414
    
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        unicorn.v414.fit_zbest_RF_colors(field=field)

    #### Parse into "master" files
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        
        filters = []
        filters.extend(['153', '154', '155']) # Maiz UBV
        filters.extend(['161', '162', '163']) # 2MASS JHK
        filters.extend(['156', '157', '158', '159', '160']) # SDSS ugriz
        filters.extend(['135','136','137','138','139']) # Bessel UX BVRI
        filters.extend(['270','271','272','273','274','275']) # Tophat filters
        
        # cat = catIO.Table('%s_3dhst.v4.1.zbest.cat' %(field))
        # rf0 = catIO.Table('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, filters[0]))
        cat = catIO.Table('%s_3dhst.v4.1.zfit.cat' %(field))
        rf0 = catIO.Table('OUTPUT/%s_3dhst.v4.1.zfit.%s.rf' %(field, filters[0]))
        for c in ['nfilt_fit', 'chi2_fit', 'L%s' %(filters[0])]:
            rf0.remove_column(c)
            
        # rf0.rename_column('z', 'z_best')
        # rf0.add_column(cat['z_type'], index=2)
        # rf0.add_column(cat['z_spec'], index=3)
        
        rf0.rename_column('z', 'z_max_grism')
        
        for f in filters:
            print '%s, filter: %s' %(field, f)
            #rf = catIO.Table('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, f))
            rf = catIO.Table('OUTPUT/%s_3dhst.v4.1.zfit.%s.rf' %(field, f))
            rf.rename_column('nfilt_fit', 'nfilt%s' %(f))
            rf0.add_columns([rf['L%s' %(f)], rf['nfilt%s' %(f)]])
        
        #rf0.write('%s_3dhst.v4.1.4.zbest.rf' %(field), format='ascii.commented_header')
        rf0.write('%s_3dhst.v4.1.5.z_max_grism.rf' %(field), format='ascii.commented_header')
    
    #
    descrip = ['#\n']
    for f in filters:
        #print '%s, filter: %s' %(field, f)
        #fp = open('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, f))
        fp = open('OUTPUT/%s_3dhst.v4.1.zfit.%s.rf' %(field, f))
        for i in range(3):
            line = fp.readline()
        
        descrip.append(line)
    
    #descrip.append('#\n# Rest-frame colors computed using templates in tweak_cosmos_v4.1/spectra.param\n#\n# z_type: 1 (phot), 2 (grism)\n#\n')
    descrip.append('#\n# Rest-frame colors computed using templates in tweak_cosmos_v4.1/spectra.param\n#\n# z_type: 0 (star), 1 (spec), 2 (grism), 3(phot)\n#\n')
    
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        #lines = open('%s_3dhst.v4.1.4.zbest.rf' %(field)).readlines()
        lines = open('%s_3dhst.v4.1.5.z_max_grism.rf' %(field)).readlines()
        for line in descrip[::-1]:
            lines.insert(1, line)
            
        #fp = open('%s_3dhst.v4.1.4.zbest.rf' %(field), 'w')
        fp = open('%s_3dhst.v4.1.5.z_max_grism.rf' %(field), 'w')
        fp.writelines(lines)
        fp.close()
    
    ###
    # for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
    #     print 'mv %s_3dhst.v4.1.4.zbest.rf %s_3dhst.v4.1.5.zbest.rf' %(field, field)
        
def fit_zbest_RF_colors(field='aegis'):
    """
    Use z_grism as the redshift for fitting EAZY RF colors
    """
    import os
    import threedhst.eazyPy as eazy
    
    param = eazy.EazyParam('%s_3dhst.v4.1.param' %(field))
    
    param['FILTERS_RES'] = 'FILTER.RES.latest'
    param['CATALOG_FILE'] = '%s_3dhst.v4.1.zfit.cat' %(field)
    param['OUTPUT_DIRECTORY'] = 'OUTPUT'
    
    #### Use z_spec / z_best
    param['BINARY_OUTPUT'] = 'n'
    param['FIX_ZSPEC'] = 'y'
    
    #### Dusty template fit
    param['TEMPLATES_FILE'] = 'tweak_cosmos_v4.1/spectra.param'
    param['CACHE_FILE'] = '%s_3dhst.v4.1.zbest.tempfilt' %(field)
    #param['MAIN_OUTPUT_FILE'] = field+'_3dhst.v4.1.zbest'
    
    param['MAIN_OUTPUT_FILE'] = field+'_3dhst.v4.1.zfit'
    
    #### Baseline
    param['REST_FILTERS'] = '---'
    param['RF_ERRORS'] = 'n'
    param.write('/tmp/zphot.%s.param' %(field))
    
    ### translate with z_best = z_spec
    lines = open('%s_3dhst.v4.1.translate' %(field)).readlines()
    
    ### z_best
    #lines.append('z_spec z_spec_old\nz_best z_spec\n')
    
    ### z_max_grism
    lines.append('z_spec z_spec_old\nz_max_grism z_spec\n')
    
    fp = open('/tmp/zphot.%s.translate' %(field), 'w')
    fp.writelines(lines)
    fp.close()
    
    os.system('eazy -p /tmp/zphot.%s.param -t /tmp/zphot.%s.translate' %(field, field))
    
    #### Loop over filters
    filters = []
    filters.extend(['153', '154', '155']) # Maiz UBV
    filters.extend(['161', '162', '163']) # 2MASS JHK
    filters.extend(['156', '157', '158', '159', '160']) # SDSS ugriz
    filters.extend(['135','136','137','138','139']) # Bessel UX BVRI
    filters.extend(['270','271','272','273','274','275']) # Tophat filters
    
    for rf in filters:
        param['REST_FILTERS'] = rf
        param['READ_ZBIN'] = 'n'
        param['RF_ERRORS'] = 'n'
        param['RF_PADDING'] = '%d' %(1500-500*(',' in rf))
        param.write('/tmp/zphot.%s.param' %(field))
        os.system('eazy -p /tmp/zphot.%s.param -t /tmp/zphot.%s.translate' %(field, field))
        
def make_master_catalog():
    """
    Combine all columns into a master FITS file
    """
    from collections import OrderedDict
    import numpy as np
    from astropy.table import Table, Column
    from threedhst import catIO
    
    #### Redshifts
    zout = catIO.Table('3dhst.dusty.v4.1.zout', format='ascii.commented_header')
    
    #### FAST
    fout = catIO.Table('3dhst.v4.1.fout', format='ascii.commented_header')
    
    #### RF colors
    rf = catIO.Table('3dhst.v4.1.master.RF')
    rf.add_column(Column(data=-2.5*np.log10(rf['L153']/rf['L155']), name='UmV'))
    rf.add_column(Column(data=-2.5*np.log10(rf['L153']/rf['L154']), name='UmB'))
    rf.add_column(Column(data=-2.5*np.log10(rf['L155']/rf['L161']), name='VmJ'))
    
    #### Grism fits
    zfit = catIO.Table('3dhst.new_zfit.v4.1.4.dat')
    
    #### Grism DQ
    dq = catIO.Table('3dhst.dqflag.v4.1.4.dat')
    
    #### Line fits
    files = glob.glob('*-WFC3_v4.1.4_SPECTRA/*catalogs/*fits')
    lf = Table.read(files[0])
    line_dict = OrderedDict()
    field = np.ones(len(lf), dtype=int)
    for line in ['OIII', 'OII', 'OIIIX', 'HALPHA', 'HBETA', 'SII']:
        for attr in ['FLUX', 'FLUX_ERR', 'SCALE', 'EQW', 'EQW_ERR']:
            col = '%s_%s' %(line, attr)
            print col
            line_dict[col] = lf[col].data
            
    for i, file in enumerate(files[1:]):
        lf = Table.read(file)
        field = np.append(field, np.ones(len(lf), dtype=int)*(i+2))
        
        for line in ['OIII', 'OII', 'OIIIX', 'HALPHA', 'HBETA', 'SII']:
            for attr in ['FLUX', 'FLUX_ERR', 'SCALE', 'EQW', 'EQW_ERR']:
                col = '%s_%s' %(line, attr)
                print col
                line_dict[col] = np.append(line_dict[col], lf[col].data)
        
    
    #### Photometric catalog
    phot = catIO.Table('../v4.1/3dhst_master.phot.v4.1/3dhst_master.phot.v4.1.cat')
    phot.add_column(Column(data=25-2.5*np.log10(phot['f_F160W']), name='hmag'))

    ok = np.isfinite(phot['hmag']) & (zout['z_peak'] >= 0) & np.isfinite(fout['lmass'])
    phot.add_column(Column(data=phot['use_phot']*ok, name='use_all'))
    
    for table, tname in zip([zout, fout, rf, zfit, dq], ['zout', 'fout','rf','zfit','dq']):
        for col in table.columns:
            if col in phot.columns:
                table.rename_column(col, '%s_%s' %(tname, col))
                col = '%s_%s' %(tname, col)
                
            print '%s %s' %(tname, col)    
            phot.add_column(table[col])
            
    for col in line_dict.keys():
        print col
        phot.add_column(Column(data=line_dict[col], name=col))
        
    #### Fraction of old + young dusty templates
    import research.v4 as cat2
    import research.dusty
    research.dusty.dusty_fraction()
    
    phot.add_column(Column(data=cat2.dusty_old, name='dusty_old'))
    phot.add_column(Column(data=cat2.dusty_young, name='dusty_young'))
    phot.add_column(Column(data=cat2.dusty_total, name='dusty_total'))
    
    phot.add_column(Column(data=cat2.eazyMass, name='eazy_mass'))
    phot.add_column(Column(data=cat2.MLv, name='eazy_MLv'))
    phot.add_column(Column(data=cat2.Lv, name='eazy_Lv'))
    phot.add_column(Column(data=cat2.template_Halpha, name='eazy_Halpha'))
    phot.add_column(Column(data=cat2.template_O2, name='eazy_OII'))
    phot.add_column(Column(data=cat2.template_O3, name='eazy_OIII'))
    
    field_id = np.ones(len(phot), dtype=int)
    for fi, field in enumerate(['AEGIS', 'COSMOS', 'GOODS-N', 'GOODS-S', 'UDS']):
        field_id[phot['field'] == field] = fi+1
    
    phot.add_column(Column(data=field_id, name='field_id'), index=2)
    
    phot.write('3dhst.v4.1.4.full.v2.fits')
    
    #full = catIO.Table('3dhst.v4.1.4.full.v1.fits')
    
    #############
    full = catIO.Table('../../v4.1.4/3dhst.v4.1.4.full.v2.fits')
    for line in ['OIII', 'OIIIX', 'OII', 'HALPHA', 'HBETA', 'SII']:
        for ext in ['SCALE', 'FLUX', 'FLUX_ERR', 'EQW', 'EQW_ERR']:
            full.remove_column('%s_%s' %(line, ext))
            
    ### update with v4.1.5 redshifts
    ### cat aegis.new_zfit.linematched.v4.1.5.dat cosmos.new_zfit.linematched.v4.1.5.dat
    fields = ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']
    ### columns that need "zfit" before the name
    prepend_columns = ['z_spec', 'l95', 'l68', 'u68', 'u95']
    for field in fields:
        print field
        #zfit = catIO.Table('Files/%s.new_zfit.linematched.v4.1.5.dat' %(field))
        zfit = catIO.Table('Files/%s.new_zfit.linematched_all.v4.1.5.fits' %(field))
        ix = full['field'] == field.upper().replace('GOODS','GOODS-')
        for col in zfit.columns:
            print col
            outcol = col
            if col in prepend_columns:
                outcol = 'zfit_%s' %(col)
            
            full[outcol][ix] = zfit[col]
    
    #### linefit
    tt = catIO.Table('Files/aegis.linefit.linematched_all.v4.1.5.fits')
    cols = list(tt.columns[3:])
    cols[0] = 'linefit_z'
    
    for col in cols:
        data = np.zeros(len(full))-99
        full.add_column(catIO.Column(data=data, name=col))
        
    for field in fields:
        print field
        linefit = catIO.Table('Files/%s.linefit.linematched_all.v4.1.5.fits' %(field))
        ix = full['field'] == field.upper().replace('GOODS','GOODS-')
        linefit.rename_column('z', 'linefit_z')
        for col in cols:
            full[col][ix] = linefit[col]
    
    full.write('3dhst.v4.1.5.full.v1.fits')

    full = catIO.Table('3dhst.v4.1.5.full.v1.fits')
#
class SpecSorted(object):
    def __init__(self, spec_rest, wave, skip, selection, sort, sort_name):
        self.spec_rest = spec_rest*1
        self.wave = wave*1
        self.sh = self.spec_rest.shape
        self.skip = skip
        self.selection = selection
        self.N = self.selection.sum()
        self.sort = sort
        self.sort_name = sort_name
        
def master_spec():
    full = catIO.Table('3dhst.v4.1.4.full.v1.fits')
    lran = np.arange(1.e4,1.7e4,23)
    spec = np.zeros((len(full), len(lran)))
    contam = np.zeros((len(full), len(lran)))
    err = np.zeros((len(full), len(lran)))
    
    N = len(full)
    PATH = '/Volumes/WD5_ESOMAC/3DHST/Spectra/Release/v4.1.4'
    for i in range(N):
        if full['z_max_grism'][i] >= 0:
            obj = full['spec_id'][i]
            field = obj.split('-')[0]
            pointing = obj.split('-G1')[0]
            try:
                print obj
                oned = unicorn.reduce.Interlace1D('%s/%s-wfc3-spectra_v4.1.4/%s/1D/FITS/%s.1D.fits' %(PATH, field, pointing, obj))
                spec[i,:] = np.interp(lran, oned.lam, oned.flux/oned.sens, left=0, right=0)
                contam[i,:] = np.interp(lran, oned.lam, oned.contam/oned.sens, left=0, right=0)
                err[i,:] = np.interp(lran, oned.lam, oned.error/oned.sens, left=0, right=0)
            except:
                print '(failed)'
                pass
    
    hdu = [pyfits.PrimaryHDU()]
    hdu.append(pyfits.ImageHDU(data=lran, name='WAVE'))
    hdu.append(pyfits.ImageHDU(data=spec, name='FLUX'))
    hdu.append(pyfits.ImageHDU(data=err, name='ERR'))
    hdu.append(pyfits.ImageHDU(data=contam, name='CONTAM'))
    
    pyfits.HDUList(hdu).writeto('master_spec_v4.1.4.fits')
    
    #### v4.1.5
    # os.chdir("/Users/brammer/3DHST/Spectra/Release/v4.1.5/Catalogs")
    # full = catIO.Table('3dhst.v4.1.5.full.v1.fits')

    os.chdir("/Users/brammer/3DHST/Spectra/Release/v4.1.5/FullRelease")
    full = catIO.Table('3dhst.v4.1.5.master.fits')

    lran = np.arange(1.e4,1.7e4,23)
    spec = np.zeros((len(full), len(lran)))
    contam = np.zeros((len(full), len(lran)))
    err = np.zeros((len(full), len(lran)))
    
    N = len(full)
    # PATH = '/Volumes/KEYTAR/INTERLACE_v4.1.5'
    # for i in range(N):
    #     #for i in range(162462, N):
    #     if full['z_max_grism'][i] >= 0:
    #         obj = full['spec_id'][i]
    #         field = obj.split('-')[0]
    #         pointing = obj.split('-G1')[0]
    #         #break
    #         print "%d %s" %(i, obj)
    #         oned_file = '%s/%s_INTERLACE_v4.1.5/%s.1D.fits' %(PATH, field.upper(), obj)
    #         if os.path.exists(oned_file):
    #             oned = unicorn.reduce.Interlace1D(oned_file)
    #             #print len(oned_file)
    #             spec[i,:] = np.interp(lran, oned.lam, oned.flux/oned.sens, left=0, right=0)
    #             contam[i,:] = np.interp(lran, oned.lam, oned.contam/oned.sens, left=0, right=0)
    #             err[i,:] = np.interp(lran, oned.lam, oned.error/oned.sens, left=0, right=0)
    #         else:
    #             print '%s (not found)' %(oned_file)
    
    
    N = len(full)
    wave = np.zeros((N, 312), dtype=np.float32)
    spec = np.zeros((N, 312), dtype=np.float32)
    contam = np.zeros((N, 312), dtype=np.float32)
    err = np.zeros((N, 312), dtype=np.float32)
    continuum = np.zeros((N, 312), dtype=np.float32)
        
    PATH = '/Volumes/KEYTAR/INTERLACE_v4.1.5'
    #for i in range(N):
    for i in range(77708, N):
        if full['npoint'][i] > 0:
            obj = full['grism_id'][i]
            if obj == '00000':
                obj = '%s-%s-G141_%05d' %(full['field'][i], full['pointings'][i].split(',')[0], full['phot_id'][i])
                #print 'x ' + obj
            else:
                pass
                #print obj
                
            field = obj.split('-')[0]
            pointing = obj.split('-G1')[0]
            #break
            print "%d %s" %(i, obj)
            oned_file = '%s/%s_INTERLACE_v4.1.5/%s.1D.fits' %(PATH, field.upper(), obj)
            if os.path.exists(oned_file):
                oned = unicorn.reduce.Interlace1D(oned_file)
                nl = len(oned.lam)
                wave[i,:nl] = oned.lam*1
                spec[i,:nl] = oned.flux/oned.sens
                contam[i,:nl] = oned.contam/oned.sens
                err[i,:nl] = oned.error/oned.sens
                zf_file = oned_file.replace('1D','new_zfit')
                if os.path.exists(zf_file):
                    zf = pyfits.open(zf_file)
                    continuum[i,:nl] = zf['CONT1D'].data/oned.sens
                
                os.remove(os.path.basename(oned_file).replace('.fits', '.png'))
            else:
                print '%s (not found)' %(oned_file)
                
    #
    hdu = [pyfits.PrimaryHDU()]
    #hdu.append(pyfits.ImageHDU(data=lran, name='WAVE'))
    hdu.append(pyfits.ImageHDU(data=wave, name='WAVE'))
    hdu.append(pyfits.ImageHDU(data=spec, name='FLUX'))
    hdu.append(pyfits.ImageHDU(data=err, name='ERR'))
    hdu.append(pyfits.ImageHDU(data=contam, name='CONTAM'))
    hdu.append(pyfits.ImageHDU(data=continuum, name='CONTIN'))
    pyfits.HDUList(hdu).writeto('master_spec_v4.1.5.v2.fits', clobber=True)
            
    #####
    
    img = pyfits.open('/Users/brammer/3DHST/Spectra/Release/v4.1.4/master_spec_v4.1.4.fits')


    from threedhst import catIO
    import numpy as np
    import astropy.io.fits as pyfits
    
    full = catIO.Table('/Users/brammer/3DHST/Spectra/Release/v4.1.5/Catalogs/3dhst.v4.1.5.full.v1.fits')    
    img = pyfits.open('/Users/brammer/3DHST/Spectra/Release/v4.1.5/Catalogs/master_spec_v4.1.5.fits')
    
    full = catIO.Table('/Users/brammer/3DHST/Spectra/Release/v4.1.5/FullRelease/3dhst.v4.1.5.master.fits')
    img = pyfits.open('/Users/brammer/3DHST/Spectra/Release/v4.1.5/FullRelease/master_spec_v4.1.5.v2.fits')
    phot = catIO.Table('/Users/brammer/3DHST/Spectra/Release/v4.1.5/Catalogs/3dhst.v4.1.5.full.v1.fits')
    
    full['UmV'] = -2.5*np.log10(full['L153']/full['L155'])
    full['VmJ'] = -2.5*np.log10(full['L161']/full['L161'])
    
    full['hmag'] = full['jh_mag']
    full['star_flag'] = (full['z_best_s'] == 0)*1
    full['z_peak'] = full['z_peak_phot']
    
    lran = img['WAVE'].data
    spec = img['FLUX'].data
    err = img['ERR'].data
    contam = img['CONTAM'].data
    contin = img['CONTIN'].data
    
    dz = np.abs(full['z_max_grism'] - full['z_peak']) / (1+full['z_peak'])
    
    ### Galaxies
    npix = (img['FLUX'].data != 0).sum(axis=1)
    ncontam = (img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1)
    
    #ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1)
    ok = (npix > 80) & (ncontam < 190) & (full['hmag'] < 25) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1) & (full['z_max_grism'] < 3.3) & (dz < 0.2) ### full plot

    ok = (npix > 80) & (ncontam < 190) & (full['hmag'] > 26) & (full['star_flag'] != 1)  ### faint zphot
    
    ok = (npix > 100) & (ncontam < 150) & (full['hmag'] < 26) & (full['hmag'] > .24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1) & (full['z_max_grism'] < 3.3)
    so = np.argsort(full['z_max_grism'][ok])
    skip = 15
    
    ok = (npix > 80) & (ncontam < 190) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.605) & (full['z_max_grism'] < 3.3) & (full['star_flag'] != 1) & (full['lmass'] > 9)
    ok = (npix > 10) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.605) & (full['z_max_grism'] < 3.3) & (full['star_flag'] != 1) #& (full['lmass'] > 9)
    
    ok = (npix > 80) & (ncontam < 190) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.15) & (full['z_max_grism'] < 3.3) & (full['star_flag'] != 1) #& (full['lmass'] > 9)
    
    skip = ok.sum()/305
    xx = full['z_max_grism'][ok]
    xx = full['z_peak_phot'][ok]
    so = np.argsort(xx)
    
    ### UVJ
    xuvj = [0, 0.475806, 0.876573, 1.297379, 1.678108, 1.711505]
    yuvj = [1.39, 1.39, 1.41, 1.69, 2.16, 2.62]
    yint_uvj = np.interp(full['VmJ'], xuvj, yuvj)
    is_uvj = full['UmV'] > yint_uvj
    
    ok = (npix > 80) & (ncontam < 190) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1) & (full['z_max_grism'] < 2.8) & (dz < 0.2) & (full['lmass'] > 9) ### full plot
    ok &= ~is_uvj
    skip = ok.sum()/305
    so = np.argsort(full['lmass'][ok])
    xx = full['lmass'][ok]
    
    #### Continuum
    line_SN = npix < 0
    no_line_SN = npix < 0
    line_list = ['OII', 'OIII', 'HALPHA']
    line_list = ['OII', 'OIII', 'Ha', 'Hb', 'SIII']
    for line in line_list:
        line_SN |= (full['%s_FLUX' %(line)]/full['%s_FLUX_ERR' %(line)] > 3) & (full['%s_SCALE' %(line)] > 0)
        #no_line_SN[full['%s_SCALE' %(line)] != 0] &= (full['%s_FLUX' %(line)]/full['%s_FLUX_ERR' %(line)] < 1)[full['%s_SCALE' %(line)] != 0]
        no_line_SN |= (full['%s_FLUX' %(line)]/full['%s_FLUX_ERR' %(line)] > 2) & (full['%s_SCALE' %(line)] > 0)
        #no_line_SN[full['%s_SCALE' %(line)] != 0] &= (full['%s_EQW' %(line)] - full['%s_EQW_ERR' %(line)] < 50)[full['%s_SCALE' %(line)] != 0]
        
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & line_SN
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & line_SN
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 25) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & ~line_SN
    so = np.argsort(full['z_max_grism'][ok])
    #skip = 5
    skip = 18
    skip = ok.sum()/305

    ### No line
    #ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1) & no_line_SN
    #ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 23) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & no_line_SN & (~line_SN)
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & ~no_line_SN #& (~line_SN)
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 26) & (full['hmag'] > 24) & (full['z_max_grism'] > 0.605) & (full['star_flag'] != 1) & ~no_line_SN #& (~line_SN)
    so = np.argsort(full['z_max_grism'][ok])
    #so = np.argsort(full['hmag'][ok])
    skip = 18
    skip = ok.sum()/305
    
    ## sort by stellar mass
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 2.6)
    ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 2.6)
    #ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 26) & (full['hmag'] > 24) & (full['z_max_grism'] > 0.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 2.6)
    #ok = (npix > 80) & (ncontam < 150) & (full['hmag'] < 26) & (full['z_max_grism'] > 0.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 3.6)
    
    ### High-z
    ok = (npix > 80) & (ncontam < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 1.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 3.2)
    ### lowz
    ok = (npix > 80) & (ncontam < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.18) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 0.9)
    
    
    so = np.argsort(full['lmass'][ok])
    xx = full['lmass'][ok]
    #skip = ok.sum()/200

    so = np.argsort(full['UmV'][ok])
    xx = full['UmV'][ok]

    so = np.argsort(full['lssfr'][ok])
    so = np.argsort(full['fout_Av'][ok])
    xx = full['fout_Av'][ok] + np.random.rand(ok.sum())*0.1
    so = np.argsort(xx)

    so = np.argsort(full['eazy_MLv'][ok])
    xx = full['eazy_MLv'][ok]

    ok &= full['galfit_J_f'] <= 1
    so = np.argsort(full['galfit_J_re'][ok])
    xx = full['galfit_J_re'][ok]

    so = np.argsort(full['Ha_EQW'][ok])
    xx = full['Ha_EQW'][ok]

    haO3 = (full['Ha_FLUX']/full['OIII_FLUX'])[ok]
    
    hbO3 = (full['Hb_FLUX']/full['OIII_FLUX'])[ok]
    hbO3[full['Hb_SCALE'][ok] < 0] = -100
    so = np.argsort(hbO3)
    xx = hbO3
    
    zerr = np.abs((full['zfit_u95']-full['zfit_l95'])/(1+full['z_max_grism']))[ok]
    zerr = np.abs((full['zfit_u68']-full['zfit_l68'])/(1+full['z_max_grism']))[ok]
    so = np.argsort(zerr)
    xx = zerr
    
    so = np.argsort(full['z_peak'][ok])
    xx = full['z_peak'][ok]
    
    so = np.argsort(full['OIIIx_FLUX'][ok])
    xx = full['OIIIx_FLUX'][ok]

    so = np.argsort(full['NeIII_FLUX'][ok])
    xx = full['NeIII_FLUX'][ok]
    
    #skip = ok.sum()/200
    skip = ok.sum()/305
    
    #### stars
    imj = -2.5*np.log10(full['f_F814W']/full['f_F125W'])
    imj = -2.5*np.log10(full['f_F125W']/full['f_F140W'])
    ok = (npix > 150) & ((img['CONTAM'].data/img['FLUX'].data > 3).sum(axis=1) < 60) & (full['hmag'] < 24) & (full['star_flag'] == 1) & np.isfinite(imj)
    so = np.argsort(imj[ok])
    skip = 5
    #avg = np.median(smooth[0:100,:], axis=0)
    
    ########## Copy from here
    
    sub = (img['FLUX'].data-img['CONTAM'].data)[ok,:][so,:]
    sub -= img['CONTIN'].data[ok,:][so,:]
    #scl = 10**(0.4*(full['hmag'][ok][so]-21))
    h2 = 25-2.5*np.log10(phot['f_F140W'][ok][so])
    scl = 10**(0.4*(h2-21))
     
    import numpy.ma as ma
    mask = np.isfinite(sub) & (sub > 0)
    sub = ma.masked_array(sub, mask=(~mask))
    
    sub = (sub.T*scl).T
    smooth = sub*1.
    for i in range(ok.sum())[::skip]:
        #smooth[i,:] = np.median(sub[i:i+skip,:], axis=0)
        smooth[i,:] = ma.median(sub[i:i+skip,:], axis=0)
    
    smooth = smooth[::skip]
    avg = np.median(smooth, axis=0)
        
    zi = full['z_max_grism'][ok][so][::skip]
    lrest = np.arange(2500, 1.4e4, 10)

    lrest = np.arange(2500, 1.7e4/(1+0.605), 10)  #### limit z range to h-alpha
    lrest = np.linspace(2500, 1.7e4/(1+0.605), smooth.shape[1])  #### limit z range to h-alpha

    sh = smooth.shape
    spec_rest = np.zeros((sh[0], len(lrest)))
    for i in range(sh[0]):
        #spec_rest[i,:] = np.interp(lrest, lran[30:]/(1+zi[i]), smooth[i,30:]/avg[30:], left=0, right=0)
        spec_rest[i,:] = unicorn.utils_c.interp_conserve_c(lrest, np.cast[np.double](lran[100,33:])/(1+zi[i]), smooth[i,33:]/avg[33:], left=0, right=0)
    
    spec_rest[spec_rest == 0] = 100
    
    # spec_rest_line = SpecSorted(spec_rest, lrest, skip, ok, zi, 'z')
    # 
    # spec_rest_cont = SpecSorted(spec_rest, lrest, skip, ok, zi, 'z')
            
    ###### To here
    
    ### full rest, for when not sorted by redshift
    #g=1
    lrest = np.linspace(3300, 8000, smooth.shape[1]*g)  #### limit z range to h-alpha
    #lrest = np.linspace(3000, 10000, smooth.shape[1]*g)  #### limit z range to h-alpha
    #lrest = 10**np.linspace(np.log10(3000), np.log10(13200), smooth.shape[1]*g)  #### limit z range to h-alpha
    #lrest = 10**np.linspace(np.log10(3300), np.log10(10000), smooth.shape[1]*g)  #### limit z range to h-alpha
    #lrest = np.linspace(3300, 13000, smooth.shape[0])  #### limit z range to h-alpha
    lami = lran[ok,:][so,:]
    zzi = full['z_max_grism'][ok][so]
    #zzi = full['z_peak'][ok][so]
    sh = sub.shape
    spec_rest_full = np.zeros((sh[0], len(lrest)))
    for i in range(sh[0]):
        print i
        #spec_rest_full[i,:] = np.interp(lrest, lran[30:-5]/(1+zzi[i]), sub[i,30:-5]/avg[30:-5], left=0, right=0)
        sl = slice(40,-10)
        sl = (lami[i,:] > 1.14e4) & (lami[i,:] < 1.58e4)
        
        lrest_i = lami[i,sl]/(1+zzi[i])
        if (lrest_i[0] > lrest[-1]) | (lrest_i[-1] < lrest[0]):
            print 'No overlap (%.2f)' %(zzi[i])
            continue
        #
        spec_rest_full[i,:] = unicorn.utils_c.interp_conserve_c(lrest, np.cast[np.float64](lrest_i), np.cast[np.float64](sub[i,sl]/avg[sl]), left=0, right=0)
    
    spec_rest_full[~np.isfinite(spec_rest_full)] = 0
    smooth_full = spec_rest_full*1.
    for i in range(ok.sum())[::skip]:
        print i
        sub_region = np.ma.masked_array(spec_rest_full[i:i+skip,:], mask=spec_rest_full[i:i+skip,:] == 0)
        smooth_full[i,:] = np.ma.median(sub_region, axis=0)
    
    smooth_full = smooth_full[::skip]
    smooth_full[smooth_full == 0] = 100
    
    #####
    
    spec_rest_z = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'z')

    spec_rest_uvj_q = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'lmass UVJ q')
    spec_rest_uvj_sf = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'lmass UVJ sf')
    
    spec_rest_z_phot = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'z zphot')
    
    spec_rest_mass_zphot = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'lmass zphot')
    
    spec_rest_mass = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'lmass')
    spec_rest_av = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'Av')
    spec_rest_UmV = SpecSorted(smooth_full, lrest, skip, ok, xx[so][::skip], r'UV')
    
    ### Stack
    # spec_rest[~np.isfinite(spec_rest) | (spec_rest > 99)] = 0
    # sum_rest = np.sum(spec_rest, axis=0)
    # N = np.sum(spec_rest != 0, axis=0)    
    # spec_rest[spec_rest == 0] = 100
    smooth_full[~np.isfinite(smooth_full) | (smooth_full > 99)] = 0
    sum_rest = np.sum(smooth_full, axis=0)
    N = np.sum(smooth_full != 0, axis=0)    
    smooth_full[smooth_full == 0] = 100
        
    #### Show lines and continuum
    import seaborn as sns
    sns.set(style="ticks")
    plt.rcParams['xtick.direction'] = u'in'
    plt.rcParams['ytick.direction'] = u'in'
    #sns.set_style("darkgrid", {"axes.facecolor": ".9"})
    plt.ioff()
    show_label=True
    fig = unicorn.plotting.plot_init(aspect=0.55, xs=8, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)
    #fig = unicorn.plotting.plot_init(aspect=1/0.55, xs=6, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)

    #for subplot, obj in zip([121, 122], [spec_rest_line, spec_rest_cont]):
    #for subplot, obj in zip([121, 122], [spec_rest_UmV, spec_rest_av]):
    for subplot, obj in zip([121, 122], [spec_rest_mass, spec_rest_av]):
        #for subplot, obj in zip([121, 122], [spec_rest_mass, spec_rest_mass_zphot]):
        #for subplot, obj in zip([121, 122], [spec_rest_uvj_sf, spec_rest_uvj_q]):
        #for subplot, obj in zip([111], [spec_rest_z_phot]):
        zi = obj.sort
        sh = obj.sh
        spec_rest = obj.spec_rest*1
        N = obj.N
        
        ax = fig.add_subplot(subplot) ### run twice
        #ax = fig.add_subplot(122) ### run twice
    
        #ax.imshow(255-unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
        if (obj.sort_name in ['lmass', 'UV']) & show_label:
            #dy=6
            spec_rest[sh[0]/18/3+dy-2:sh[0]/18+dy,:] += 100
            for line in [[3727, '[OII]'], [4341, r'H$\gamma$'], [4980, r'H$\beta$ [OIII]'], [6640, r'H$\alpha$ [SII]']]:
                xpix = np.interp(line[0], obj.wave, np.arange(len(obj.wave)))
                ax.text(xpix, sh[0]/45.+dy, line[1], ha='center', va='baseline', zorder=10, fontsize=8)
                
        ax.imshow(255-unicorn.candels.clipLog(spec_rest, lexp=1000, cmap=[9.62222, 0.90352], scale=[-0.02, 2]), origin='lower', aspect='auto', interpolation='Nearest')
    
        if 'z' in obj.sort_name:
            if subplot == 111:
                xlam = np.arange(0.3,1.31, 0.1)            
                zpix = np.arange(0.25,3.1,0.25)
            else:
                xlam = np.arange(0.3,1.01, 0.1)
                zpix = np.arange(0.75,3.1,0.25)
            
            if 'zphot' in obj.sort_name:
                ax.set_ylabel(r'$z_\mathrm{phot}$')
            else:
                ax.set_ylabel(r'$z$')
                
            zpixstr = list(zpix)
            for i in range(len(zpix))[::2]:
                zpixstr[i] = ''
        
            if subplot == 122:
                zpixstr[-1] = ''
        
        elif 'lmass' in obj.sort_name:
            xlam = np.arange(0.4,0.81, 0.1)
            ax.set_ylabel(r'$\log\ M/M_\odot$')
            zpix = np.arange(9,11.5,0.25)
            zpixstr = list(zpix)
            for i in range(len(zpix))[1::2]:
                zpixstr[i] = ''
            
            if 'zphot' in obj.sort_name:
                ax.text(0.99, 0.99, r'$z_\mathrm{phot}$', ha='right', va='top', fontsize=14, backgroundcolor='white', transform=ax.transAxes)
                show_label=False
                
        elif obj.sort_name == 'Av':
            xlam = np.arange(0.4,0.81, 0.1)
            ax.set_ylabel(r'$A_V$')
            zpix = np.arange(0,3.6,0.5)
            zpixstr = list(zpix)
            zpixstr[-3] = ''
            zpixstr[-1] = ''
        
        elif obj.sort_name == 'UV':
            xlam = np.arange(0.4,0.81, 0.1)
            ax.set_ylabel(r'$(U-V)_\mathrm{rest\ frame}$')
            zpix = np.arange(0.5,2.1, 0.5)
            zpixstr = list(zpix)
            # zpixstr[-3] = ''
            # zpixstr[-1] = ''
        
        ypix = np.interp(zpix, zi, np.arange(sh[0]))
        ax.set_yticks(ypix); ax.set_yticklabels(zpixstr)
        ax.set_ylim(0, sh[0])
    
        ax2 = ax.twinx()
        #ax2.imshow(unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
        ynticks = np.arange(0, N, 1000)
        ynv = ynticks/obj.skip
        ax2.set_yticks(ynv, minor=True)
        #
        ynticks = np.arange(0, N, 5000)
        ynv = ynticks/obj.skip
        ax2.set_yticks(ynv, minor=False)
        #
        ax2.set_yticklabels([]) #ynticks)
        ax2.set_ylim(0, sh[0])
        #ax2.set_ylabel(r'$N$')
        
        if subplot == 111:
            ax2.set_ylabel(r'$N\times1000$')
            keep_label = [0,2,5,10,20]
            if ok.sum() > 30000:
                keep_label.append(30)
            
            if ok.sum() > 40000:
                keep_label.append(40)
            
            yy = ['']*len(ynv)
            for i in range(len(yy)):
                if ynticks[i]/1000 in keep_label:
                    yy[i] = r'%d' %(ynticks[i]/1000)
            
            yy[0] = '0'
            
            ax2.set_yticklabels(yy)
                     
        #xlam = [0.25,0.5,0.75,1,1.25]
        xpix = np.interp(xlam, obj.wave/1.e4, np.arange(len(obj.wave)))
        ax.set_xticks(xpix); ax.set_xticklabels(xlam)
        ax2.set_xticks(xpix); ax2.set_xticklabels(xlam)
        ax.set_xlabel(r'$\lambda_\mathrm{rest}$ / $\mu\mathrm{m}$')
    
    if 'z' not in obj.sort_name:
        for ax in fig.axes:
            ax.tick_params(axis='both', color='0.9', width=1.5, which='both')
        #
    else:
        if show_label:
            lines = [[3969, r'Balmer/4000\AA'+'\nbreak (abs.)', -1], [3727, '[OII]', -1], 
                     [4341, r'H$\gamma$', -1], [4861, r'H$\beta$', -1], [5007, r'[OIII]', -1], [6563, r'H$\alpha$+[NII]', -1], [6724, '[SII]', -1], [9068, '[SIII]', 1], [9530, '[SIII]', 1], [10830, 'HeI', 1]]
            for line in lines:
                l0, llabel, top = line
                if top < 0:
                    yl = np.interp(1.1e4/l0-1, zzi, np.arange(ok.sum()))-0.01*ok.sum()
                    hal = 'right'
                    val = 'top'
                else:
                    yl = np.interp(1.68e4/l0-1, zzi, np.arange(ok.sum()))+0.01*ok.sum()
                    hal = 'left'
                    val = 'bottom'
                    
                xl = np.interp(l0, obj.wave, np.arange(len(obj.wave)))
                ax.text(xl, yl/obj.skip, llabel, ha=hal, va=val, zorder=1000, fontsize=9, rotation=45)
                print xl, yl, llabel
                
    fig.tight_layout(pad=0.3)
    
    unicorn.plotting.savefig(fig, '3dhst_allspecs_x.pdf', dpi=300)
    plt.close()

def twod_headers():
    """
    Pull information from 2D.fits headers, e.g. for getting FLT positions
    """
    os.chdir("/Users/brammer/3DHST/Spectra/Release/v4.1.5/FullRelease")
    full = catIO.Table('3dhst.v4.1.5.master.fits')

    N = len(full)
    
    PATH = '/Volumes/KEYTAR/Release/v4.1.5/Spectra'
    
    idx = np.arange(N)
    #ok = (full['z_spec'] > 0) & (full['npoint'] > 0)
    ok = (full['npoint'] > 0)
    
    import collections
    header_data = collections.OrderedDict()
    header_data['obj'] = ['']*N
    
    twod0 = pyfits.open('/Volumes/KEYTAR/Release/v4.1.5/Spectra/aegis-02/2D/FITS/aegis-02-G141_07223.2D.fits')
    for ext in [0, 'SCI']:
        h = twod0[ext].header
        for key in h:
            header_data[key] = [h[key]*0]*N
        
    for i in range(191370, N):
        print i
        if not ok[i]:
            continue
        
        #for i in range(77708, N):
        obj = full['grism_id'][i]
        if obj == '00000':
            obj = '%s-%s-G141_%05d' %(full['field'][i], full['pointings'][i].split(',')[0], full['phot_id'][i])
            #print 'x ' + obj
        else:
            pass
            #print obj
            
        field = obj.split('-')[0]
        pointing = obj.split('-G1')[0]
        #break
        print "%d %s" %(i, obj)
        twod_file = '%s/%s/2D/FITS/%s.2D.fits' %(PATH, pointing, obj)
        if not os.path.exists(twod_file):
            continue
        
        twod = pyfits.open(twod_file)
        if 'SCI' not in twod:
            print 'Skip'
            continue
        
        header_data['obj'][i] = obj
        
        for ext in [0, 'SCI']:
            h = twod[ext].header
            for key in h:
                header_data[key][i] = h[key]
    
    tab = catIO.table_base()
    for col in header_data.keys():
        print col
        tab.add_column(catIO.Column(data=header_data[col], name=col))
    
    try:
        os.remove('3dhst.v4.1.5.2D_header.fits')
    except:
        pass
    
    tab.write('3dhst.v4.1.5.2D_header.fits')
                            
def flt_redshift_residuals():
    """
    Compute FLT redshift residuals
    """
    import collections
    
    full = catIO.Table('3dhst.v4.1.5.master.fits')
    twod = catIO.Table('3dhst.v4.1.5.2D_header.fits')
    
    ha_sn = full['Ha_FLUX']/full['Ha_FLUX_ERR']
    o3_sn = full['OIII_FLUX']/full['OIII_FLUX_ERR']
    line_snlim = 2

    dz = (full['z_max_grism'] - full['z_spec'])/(1+full['z_spec'])
    
    x_flt = (twod['X_PIX'] - twod['PAD']/2 - twod['NGROWX']*twod['GROWX']) / twod['GROWX']
    y_flt = (twod['Y_PIX'] - twod['PAD']/2 - twod['NGROWY']*twod['GROWY']) / twod['GROWY']

    #### Get field-dependent parameters for test27s.gbb
    lx, dldp = get_dldp_wavelengths(dldp_field=None, x_flt=x_flt+0, y_flt=y_flt+0, dx=x_flt*0)

    xspec_flt = (twod['XGRISM0'] - twod['PAD']/2 - twod['NGROWX']*twod['GROWX']) / twod['GROWX']
    yspec_flt = (twod['YGRISM0'] - twod['PAD']/2 - twod['NGROWY']*twod['GROWY']) / twod['GROWY']
        
    #### Pixel values of emission lines 
    lrest = collections.OrderedDict()
    xpix = collections.OrderedDict()
    #xpix2 = collections.OrderedDict()
    dxpix = collections.OrderedDict()
    #dxpix2 = collections.OrderedDict()
    
    for col in ['z_spec', 'z_max_grism', 'z_grism_l68', 'z_grism_u68']:    
        print col
        lrest[col] = full[col]*0
        lrest[col][ha_sn > line_snlim] = 6563.*(1+full[col][ha_sn > line_snlim])
        lrest[col][o3_sn > line_snlim] = 5007.*(1+full[col][o3_sn > line_snlim])
    
        #xpix[col] = (lrest[col] - twod['CRVAL1'])/(twod['CD1_1'])        
        #dxpix[col] = xpix[col]/twod['GROWX'] + xspec_flt - x_flt
    
        ### need quadratic xpix for test27s
        dxpix[col] = (-dldp[1] + np.sqrt(dldp[1]**2-4*dldp[2]*(dldp[0]-lrest[col])))/(2*dldp[2])/np.sqrt(1+dldp[-1]**2)
        xpix[col] = (dxpix[col] + x_flt - xspec_flt)*twod['GROWX']
    
    dpix = (xpix['z_max_grism']-xpix['z_spec'])/twod['GROWX']
    
    #### Selection
    ok = (full['z_spec'] > 0) & (full['grism_id'] == twod['obj'])
    ok &= (ha_sn > line_snlim) | (o3_sn > line_snlim)
    #ok &= np.abs(dz) < 0.01
    ok &= np.abs(dpix) < 2
    
    ### GOODS-S plot
    gs = ok & (full['dec'] < -26)
    #plt.xlim(53.26, 52.98); plt.ylim(-27.97, -27.66)
    
    mos = pyfits.open('/Users/brammer/3DHST/Ancillary/Mosaics/goodss_3dhst.v4.0.F160W_orig_sci.fits')
    
    reg_file = '/Users/brammer/3DHST/Spectra/Work/FIGS/12177.reg'
    r = pyregion.open(reg_file).as_imagecoord(header=mos[0].header)
    patches = r.get_mpl_patches_texts()

    fig = plt.figure(figsize=[9,4])
    ax = fig.add_subplot(121)
    ax.scatter(full['x'][gs], full['y'][gs], color='k', alpha=0.5)
    #
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.set_xlabel('RA'); ax.set_ylabel('Dec')
    
    ax2 = fig.add_subplot(122)
    ax2.scatter(x_flt[gs], y_flt[gs], color='black', alpha=0.5)
    ax2.set_xlim(-200,1014); ax2.set_ylim(-10,1024)
    ax2.plot([0,1014,1014,0,0], [0,0,1014,1014,0], color='red', linewidth=2)
    ax2.set_xlabel(r'$x$ FLT'); ax2.set_ylabel(r'$y$ FLT')
    fig.tight_layout(pad=1)

    for patch in patches[0]:
        p = ax.add_patch(patch)
        p.set_color('orange')
        p.set_alpha(0.1)
        p.set_fc('None')

    fig.savefig('../Wavecal/GS_layout.png') #, dpi=300)
    
    go_plot=True
    if go_plot:
        #### Show residuals in 2D
        fig = plt.figure()
        plt.set_cmap('cubehelix')
        sc = plt.scatter(x_flt[ok], y_flt[ok], c=dpix[ok], vmin=-1.5, vmax=1.5, alpha=0.7, marker='s', s=80, edgecolor='None')
        cb = plt.colorbar(sc)
        cb.set_label(r'$\Delta$ pix')
    
        plt.plot([0,0,1014,1014,0], [0,1014,1014,0,0], color='black', linestyle='--')
        plt.xlim(-200,1064); plt.ylim(-50,1064)
        plt.xlabel('X_FLT'); plt.ylabel('Y_FLT')
        fig.tight_layout(pad=0.5)
    
        #### Show residuals
        fig = plt.figure()
        plt.scatter(lrest['z_spec'][ok], (xpix['z_spec'][ok]-xpix['z_max_grism'][ok])/twod['GROWX'][ok], alpha=0.3, marker='s', color='black', s=30)
        xm, ym, ys, n = threedhst.utils.runmed(lrest['z_spec'][ok], (xpix['z_spec'][ok]-xpix['z_max_grism'][ok])/twod['GROWX'][ok], NBIN=20, use_median=True, use_nmad=True)
        plt.plot(xm, ym, color='red', linewidth=4, alpha=0.6)
        plt.fill_between(xm, ym+ys, ym-ys, color='red', alpha=0.2)
        plt.xlabel(r'$\lambda$'); plt.ylabel(r'$\Delta$ pix')
        fig.tight_layout(pad=0.5)
    
    ### Test optimization, black points should be zero residuals [OK]
    # ltest, dldp_ok = get_dldp_wavelengths(dldp_field=None, x_flt=x_flt[ok]+0, y_flt=y_flt[ok]+0, dx=dxpix['z_spec'][ok])
    # plt.scatter(lrest['z_spec'][ok], ltest-lrest['z_spec'][ok], alpha=0.1, color='black')
    # plt.scatter(lrest['z_spec'][ok], ltest-lrest['z_max_grism'][ok], color='red', alpha=0.1)
    # 
    # print dldp_ok[1]/twod['CD1_1'][ok]/twod['GROWX'][ok]
    # print (dldp_ok[0]+(xspec_flt-x_flt)[ok]*dldp_ok[1])/twod['CRVAL1'][ok]
    
    ###### Fit it!
    init =  [8987, 0   , 0   , 0   , 0    , 0    , 44.9, 0   , 0   , 0   , 0   , 0   ]
    delta = [30  , 1e-1, 1e-1, 1e-3, 1.e-3, 1.e-3,  2  , 1e-2, 1e-2, 1e-6, 1e-6, 1e-6]
    delta = [30  , 0e-1, 0e-1, 0e-3, 0.e-3, 0.e-3,  2  , 0e-2, 0e-2, 0e-6, 0e-6, 0e-6]
    delta = [30  , 1e-1, 1e-1, 0e-3, 0.e-3, 0.e-3,  2  , 1e-2, 1e-2, 0e-6, 0e-6, 0e-6]
    
    init = [  8.98763363e+03,   2.46977484e-02,  -1.66954148e-02, 5.27231086e-05,   5.84052102e-05,   9.30612950e-06,
              4.48741889e+01,   1.05867596e-03,   2.83833876e-03, -1.20922575e-06,  -6.34223284e-07,  1.20093841e-07]
    npar = 2
    
    ### 2nd order
    init.extend([0.003727603229696248, 5.351182862289981e-6, -3.158204223825075e-6, -1.0285571344257502e-8, 3.4669363123239304e-9, 5.786658518854924e-10])
    delta.extend([1.e-2, 1.e-6, 1.e-6, 1.e-8, 1.e-8, 1.e-9])
    npar = 3
    
    init = np.array([0.]*20)
    init[0:6] += [  8.98763363e+03,   2.46977484e-02,  -1.66954148e-02, 5.27231086e-05,   5.84052102e-05,   9.30612950e-06]
    init[10:16] += [4.48741889e+01,   1.05867596e-03,   2.83833876e-03, -1.20922575e-06,  -6.34223284e-07,   1.20093841e-07]
    delta = np.array([1e-1]*20)
    delta[0] = 30; delta[10] = 2
    delta[1:3] = 1e-1; delta[11:13] = 1e-2
    delta[3:6] = 1e-3; delta[13:16] = 1.e-6
    delta[6:10] = 5e-8; delta[16:] = 5e-12
    npar = 2
    
    NWALKERS, NSTEP = 100, 2000
    p0 = [init+np.random.normal(size=len(init))*delta for i in range(NWALKERS)]
    
    NP = len(init)
    obj_fun = unicorn.v414._objective_dldp
    obj_args = [npar, x_flt[ok], y_flt[ok], lrest['z_spec'][ok], dxpix['z_max_grism'][ok], dxpix['z_grism_l68'][ok]] 
    
    unicorn.v414._objective_dldp(init, *obj_args) #x_flt, y_flt, lrest, dx, dxlo)
        
    NTHREADS=8
    ndim = len(init)
        
    sampler = emcee.EnsembleSampler(NWALKERS, ndim, obj_fun, threads=NTHREADS, args=obj_args)
    result = sampler.run_mcmc(p0, NSTEP)
    
    NP = len(init)/npar    
    param_names = ['x%d_%d' %(j, i+1) for j in range(npar) for i in range(NP)]
    chain = unicorn.interlace_fit.emceeChain(chain=sampler.chain, param_names = param_names)
    dldp_field = np.reshape(chain.map, (npar, -1)) #[chain.map[0:NP/2], chain.map[NP/2:]]
    
    if False:
        ## test41
        dldp_field = []
        dldp_field.append([8951.386205717843, 0.08044032819916265, -0.009279698766495334, 0.000021856641668116504, -0.000011048008881387708, 0.00003352712538187608])
        dldp_field.append([44.97227893276267, 0.0004927891511929662, 0.0035782416625653765, -9.175233345083485e-7, 2.2355060371418054e-7, -9.258690000316504e-7])
        chain_label = 'test41'
        
        ### test27s
        dldp_field = []
        dldp_field.append([9028.930809875766, 0.011337227492123389, -0.079859567190042, 0.00003518029786572532, 0.0000681593925920102, 0.000060520907385702346]) 
        dldp_field.append([44.260335031281336, -0.0004031703834933002, 0.0032719446699151777, 1.1145275992836486e-6, -6.014104797258286e-7, -9.856215339544352e-8])
        dldp_field.append([0.003727603229696248, 5.351182862289981e-6, -3.158204223825075e-6, -1.0285571344257502e-8, 3.4669363123239304e-9, 5.786658518854924e-10])
        chain_label = 'test27s'
        
    lobs, dldpx = unicorn.v414.get_dldp_wavelengths(dldp_field=dldp_field, x_flt=x_flt[ok], y_flt=y_flt[ok], dx=dxpix['z_max_grism'][ok])
    
    #plt.set_cmap('cubehelix'); plt.close()
    fig = plt.figure(figsize=[11,4])
    ax = fig.add_subplot(121)
    
    sc = ax.scatter(x_flt[ok], y_flt[ok], c=lobs-lrest['z_spec'][ok], vmin=-80, vmax=80, alpha=0.7, marker='s', s=40, edgecolor='None')
    cb = plt.colorbar(sc)
    cb.set_label(r'$\Delta\lambda$')

    ax.plot([0,0,1014,1014,0], [0,1014,1014,0,0], color='black', linestyle='--')
    ax.set_xlim(-200,1064); ax.set_ylim(-50,1064)
    ax.set_xlabel('X_FLT'); ax.set_ylabel('Y_FLT')
    #fig.tight_layout(pad=0.5)
    
    ax = fig.add_subplot(122)
    ax.scatter(lrest['z_spec'][ok]/1.e4, lobs-lrest['z_spec'][ok], alpha=0.3, marker='s', color='black', s=30)
    xm, ym, ys, n = threedhst.utils.runmed(lrest['z_spec'][ok], lobs-lrest['z_spec'][ok], NBIN=20, use_median=True, use_nmad=True)
    ax.plot(xm/1.e4, ym, color='red', linewidth=4, alpha=0.6)
    ax.fill_between(xm/1.e4, ym*0+ys, ym*0-ys, color='red', alpha=0.2)
    ax.set_xlabel(r'$\lambda$'); ax.set_ylabel(r'$\Delta\lambda$')
    ax.set_ylim(-80,80)
    ax.grid(alpha=0.6)
    
    fig.tight_layout(pad=0.5)
    fig.savefig('chain_%s.png' %(chain_label))
    
    yp, xp = np.indices((128,128))
    ddx,lr = yp*0+100, [13400,13900]
    ddx,lr = yp*0+60, [11600,12000]
    ddx,lr = yp*0+150, [15650,16300]
    lgrid, dldpx = unicorn.v414.get_dldp_wavelengths(dldp_field=dldp_field, x_flt=xp.flatten()*8, y_flt=yp.flatten()*8, dx=ddx.flatten())
    lgrid1, dldpx = unicorn.v414.get_dldp_wavelengths(dldp_field=dldp_field, x_flt=xp.flatten()*8, y_flt=yp.flatten()*8, dx=ddx.flatten()+1)
    
    fig = plt.figure(figsize=[11,4])
    ax = fig.add_subplot(121)
    ims = ax.imshow(lgrid.reshape(yp.shape), interpolation='Nearest', origin='lower', vmin=lr[0], vmax=lr[1])
    levels = np.linspace(lr[0], lr[1], 7)
    ax.contour(lgrid.reshape(yp.shape), level=levels, colors='white', alpha=0.8, linewidth=2)
    ax.contour(lgrid.reshape(yp.shape), level=levels, colors='black', alpha=0.6, linewidth=1.)
    ax.text(1.05, 0.98, chain_label, ha='left', va='top', transform=ax.transAxes)
    
    ax.set_xticklabels([]); ax.set_yticklabels([])
    cb = plt.colorbar(ims, shrink=0.8)
    cb.set_label(r'$\lambda,\ \Delta x=%d$' %(ddx[0,0]))
    
    ax = fig.add_subplot(122)
    ims = ax.imshow((lgrid1-lgrid).reshape(yp.shape), interpolation='Nearest', origin='lower', vmin=43, vmax=48)
    levels = np.linspace(43,48, 7)
    ax.contour(lgrid.reshape(yp.shape), level=levels, colors='white', alpha=0.8, linewidth=2)
    ax.contour(lgrid.reshape(yp.shape), level=levels, colors='black', alpha=0.6, linewidth=1)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    cb = plt.colorbar(ims, shrink=0.8)
    cb.set_label(r'$\Delta\lambda/$pix')
    
    fig.tight_layout(pad=0.3)
    fig.savefig('wavepars_%s.png' %(chain_label))
    
    #### Show chain RMS
    chain = unicorn.interlace_fit.emceeChain(file='chain_order4_2K.fits') 
    NDRAW = 200
    draws = chain.draw_random(N=NDRAW)
    
    yp, xp = np.indices((128,128))
    ddx,lr = yp*0+150, [15650,16300]
    
    out = np.zeros((NDRAW, yp.size))
    NP = draws.shape[1]
    for i in range(NDRAW):
        print i
        dldp_field = [draws[i,0:NP/2], draws[i,NP/2:]]
        lgrid, dldpx = unicorn.v414.get_dldp_wavelengths(dldp_field=dldp_field, x_flt=xp.flatten()*8, y_flt=yp.flatten()*8, dx=ddx.flatten())
        out[i,:] = lgrid*1
    
    rms = np.std(out, axis=0).reshape(yp.shape)

    fig = plt.figure(figsize=[11,4])
    
    ax = fig.add_subplot(221)
    chain.show_chain('x0_1', ax=ax, color='black', alpha=0.1)
    ax.set_ylabel('DLDP_A_0 [0]'); ax.set_xticklabels([]); ax.set_xlabel('')

    ax = fig.add_subplot(223)
    chain.show_chain('x1_1', ax=ax, color='black', alpha=0.1)
    ax.set_ylabel('DLDP_A_1 [0]')
    
    ax = fig.add_subplot(122)
    ims = ax.imshow(rms, interpolation='Nearest', origin='lower', vmin=0, vmax=10)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    cb = plt.colorbar(ims, shrink=0.8)
    cb.set_label(r'$\sigma(\lambda)$')
    fig.tight_layout(pad=0.3)
    
    fig.savefig('wavelength_uncertainty.png')
    
def _objective_dldp(params, npar, x_flt, y_flt, lrest, dx, dxlo):
    """
    Objective function for EMCEE fitting
    """
    import numpy as np
    #NP = len(params)
    #dldp_field = [params[0:NP/2], params[NP/2:]]
    dldp_field = np.reshape(params, (npar, -1))
    #print 'SHAPE:', dldp_field.shape
    
    lfit, dldp = get_dldp_wavelengths(dldp_field=dldp_field, x_flt=x_flt, y_flt=y_flt, dx=dx)
    lfit_lo, dldp = get_dldp_wavelengths(dldp_field=dldp_field, x_flt=x_flt, y_flt=y_flt, dx=dxlo)
    
    lnprob = -0.5*np.sum((lfit-lrest)**2/(lfit-lfit_lo)**2)
    pstr = ' '.join(['%11.3e' %(p) for p in params])
    print '%s %14.2f' %(pstr, lnprob)
    
    return lnprob
    
def get_dldp_wavelengths(dldp_field=None, x_flt=None, y_flt=None, dx=None):
    """
    Function to compute the pixel wavelengths for a given FLT position + dx & configuration polynomials
    """
    import numpy as np
    
    ### test41.gbb
    # if dldp_field is None:
    #     dldp_field = []
    #     dldp_field.append([8951.386205717843, 0.08044032819916265, -0.009279698766495334, 0.000021856641668116504, -0.000011048008881387708, 0.00003352712538187608])
    #     dldp_field.append([44.97227893276267, 0.0004927891511929662, 0.0035782416625653765, -9.175233345083485e-7, 2.2355060371418054e-7, -9.258690000316504e-7])
        
    ### test41.gbb
    #dydx_1_field = [0.010205281672977665, -6.06056923866002e-6, -3.2485600412356953e-6, 4.2363866304617406e-10, 1.230956851333159e-8, 1.6123073931033502e-9]
    
    ## test27s.gbb
    if dldp_field is None:
        dldp_field = []
        dldp_field.append([9028.930809875766, 0.011337227492123389, -0.079859567190042, 0.00003518029786572532, 0.0000681593925920102, 0.000060520907385702346]) 
        dldp_field.append([44.260335031281336, -0.0004031703834933002, 0.0032719446699151777, 1.1145275992836486e-6, -6.014104797258286e-7, -9.856215339544352e-8])
        dldp_field.append([0.003727603229696248, 5.351182862289981e-6, -3.158204223825075e-6, -1.0285571344257502e-8, 3.4669363123239304e-9, 5.786658518854924e-10])
            
    ### test27s.gbb
    dydx_1_field = [0.010205281672977665, -6.06056923866002e-6, -3.2485600412356953e-6, 4.2363866304617406e-10, 1.230956851333159e-8, 1.6123073931033502e-9]
    
    order = int(-1+np.sqrt(1+8*len(dldp_field[0])))/2
    
    xy = []
    for p in range(order):
        for px in range(p+1):
            xy.append(x_flt**(p-px)*y_flt**(px))

    ## Evaluate the polynomial, allowing for N-dimensional inputs
    dydx_1 = np.sum((np.array(xy).T[:,:len(dydx_1_field)]*dydx_1_field).T, axis=0)
    dp = np.sqrt(1+dydx_1**2)*dx
    #dp = dx
    
    lam = dp*0.
    dldp = []
    for i in range(len(dldp_field)):
        dldp_i = np.sum((np.array(xy).T*dldp_field[i]).T, axis=0)
        lam += dldp_i*dp**i
        dldp.append(dldp_i)
        
    dldp.append(dydx_1)
    
    return lam, dldp
    
    
def pure_line_flux():
    """
    Calculate lien flux for pure emission line at H140=26
    """
    waves = np.linspace(1.2e4,1.6e4,100)
    fluxes = waves*0.
    for i in range(len(waves)):
        print waves[i]
        G = S.GaussianSource(1.e-17, waves[i], 10)
        sp = G.renorm(26, 'ABMag', S.ObsBandpass('wfc3,ir,f140w'))
        fluxes[i] = np.trapz(sp.flux, sp.wave)
        
    ### EW limit
    line = S.GaussianSource(5.e-17, 1.4e4, 10)
    cont = S.FlatSpectrum(26, fluxunits='ABMag')
    
    line_obs = S.Observation(line, S.ObsBandpass('wfc3,ir,f140w'))
    cont_obs = S.Observation(cont, S.ObsBandpass('wfc3,ir,f140w'))
    
    scl = (cont_obs.countrate() - line_obs.countrate()) / (cont_obs.countrate())
    sp_line = cont*scl+line
    sp_cont = cont*scl+line*0
    sp_line.convert('flam')
    sp_cont.convert('flam')
    EW = np.trapz((sp_line.flux-sp_cont.flux)/(sp_cont.flux), sp_cont.wave)
    
def line_sensitivities():
    """
    Plots like in Brammer 2012
    """
    
    import pysynphot as S
    g141 = S.ObsBandpass('wfc3,ir,g141')
    
    line_id, line_wave = 'Ha', 6563.

    #line_id, line_wave = 'OIII', 5007.
    #line_id, line_wave = 'OII', 3727.
    
    lx = (1+full['z_max_grism'])*line_wave
    ok = (full['%s_SCALE' %(line_id)] > 0) & (full['star_flag'] == 0) & (full['hmag'] < 24) & (lx > 1.1e4) & (lx < 1.65e4) & (full['%s_FLUX_ERR' %(line_id)] > 0) #(full['z_max_grism'] > 0.68) & (full['z_max_grism'] < 1.53)
    
    pow = 2
    
    inv_through = 1/g141.throughput**pow
    inv_through /= np.interp(1.5e4, g141.wave, inv_through)

    ### Throughput / z / wavelength
    interp_z = np.interp(full['z_max_grism'], g141.wave/line_wave-1, inv_through)
    
    #xm, ym, ys, nn = threedhst.utils.runmed((1+full['z_max_grism'][ok])*line_wave/1.e4, full['%s_FLUX_ERR' %(line_id)][ok]/interp_size[ok])
    #cc = polyfit(xm-1.49, ym, 2)    
    #interp_z = np.interp((1+full['z_max_grism'])*line_wave/1.e4, xm, ym/ym.min())
    
    xsize = np.arange(0,1000)
    ysize = (xsize/4.)    
    size_key = 'flux_radius'
    #size_key = 'kron_radius'
    interp_size = np.interp(full[size_key], xsize, ysize)
    
    ### size
    okx = ok & (np.abs(full[size_key]-5) < 0.2)
    size_scl = np.median(full['%s_FLUX_ERR' %(line_id)][okx]/interp_z[okx]) 
    size_scl = 0.689
    
    plt.ioff()
    fig = unicorn.plotting.plot_init(aspect=0.5, xs=8, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)

    ax = fig.add_subplot(121)
    aa = ax.scatter((1+full['z_max_grism'][ok])*line_wave/1.e4, full['%s_FLUX_ERR' %(line_id)][ok]/interp_size[ok], alpha=0.1, color='black', marker='.')
    #ax.scatter(full['z_max_grism'][ok], full['%s_FLUX_ERR' %(line_id)][ok], c=full[size_key][ok], vmin=0, vmax=10, alpha=0.2)
    #ax.scatter((1+full['z_max_grism'][ok])*line_wave/1.e4, full['%s_FLUX_ERR' %(line_id)][ok]/interp_size[ok]/interp_z[ok], c=full[size_key][ok], vmin=0, vmax=10, alpha=0.2)
    
    #aa = ax.scatter([0],[0], c=[0], vmin=0, vmax=10, alpha=0.8)
    
    #ax.plot(xm, ym, color='orange', linewidth=3)
    
    ax.plot(g141.wave/1.e4, inv_through*size_scl, color='red', linewidth=3, alpha=0.5, label=r'$\sigma\propto$ (G141 throughput)$^{-2}$')
    ax.legend(loc='upper left', fontsize=9)
    
    ax.set_xlim(1.08, 1.7)
    ax.set_ylim(0,3)
    ax.set_xlabel(r'Line $\lambda$')
    ax.set_ylabel(r'Line uncertainty, 1$\sigma$, $R=4$ pix') #, ($\times 10^{-17}\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^2$)')
    
    # cb = plt.colorbar(aa)
    # cb.set_label(size_key.replace('_','\_'))
    
    ax = fig.add_subplot(122)
    aa = ax.scatter(full[size_key][ok], full['%s_FLUX_ERR' %(line_id)][ok]/interp_z[ok], alpha=0.1, color='black', marker='.')
    #ax.scatter(full[size_key][ok], full['%s_FLUX_ERR' %(line_id)][ok]/interp_z[ok]/interp_size[ok], c=full['z_max_grism'][ok], vmin=0.7, vmax=1.5, alpha=0.2)
    #aa = ax.scatter([0],[0], c=[0], vmin=1.1, vmax=1.6, alpha=0.8)
    
    ax.plot(xsize, ysize*size_scl, color='red', linewidth=3, alpha=0.5, label=r'$\sigma\propto R$')
    ax.legend(loc='upper left', fontsize=9)
        
    ax.set_xlim(2.8,15)
    ax.semilogx()
    ax.set_ylim(0,3)
    ax.set_xlabel(r'$R = $ %s, pix' %(size_key.replace('_','\_')))
    
    #ax.set_ylabel(r'Line uncertainty, 1$\sigma$ ($\times 10^{-17}\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^2$)')
    ax.set_xticks([5,10,15])
    ax.set_xticklabels([5,10,15])
    ax.set_ylabel(r'Line uncertainty, 1$\sigma$, $\lambda=1.5~\mu\mathrm{m}$') #, ($\times 10^{-17}\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^2$)')
    
    # cb = plt.colorbar(aa)
    # cb.set_label(r'Line $\lambda$')
    
    fig.tight_layout(pad=0.5)
    
    fig.savefig('3dhst_line_sensitivity.pdf')
    
    
def get_full_catalog():
    """
    Helper function to read the master catalog
    
    Example workflow:
    
    ur_setup common ssbx
    cd /Library/WebServer/Documents/P/GRISM_v4.1.4/Selections
    ipython # start python
    
    import unicorn.v414
    full = unicorn.v414.get_full_catalog()
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['OIII_EQW'] > 2000) & (full['OIII_EQW']/full['OIII_EQW_ERR'] > 3)
    
    unicorn.v414.make_selection_webpage(full, selection, output='highEW_OIII_GBr.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism'])
    
    #### COSMOS z_spec
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['field'] == 'COSMOS') & (full['z_spec'] > 0)
    
    unicorn.v414.make_selection_webpage(full, selection, output='cosmos_zspec_GBr.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'z_spec'])

    #### z_spec outliers
    import numpy as np
    dz = (full['z_max_grism']-full['z_spec']) / (1+full['z_spec'])
    from astropy.table import Column
    full.add_column(Column(data=dz, name='dz', format='%.3f'))
    full.add_column(Column(data=np.abs(dz), name='abs(dz)', format='%.3f'))
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['z_spec'] > 0) & (np.abs(dz) > 0.08)

    #selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['z_spec'] > 0) & (np.abs(dz) > 0.05) & (np.abs(dz) < 0.1)
    
    unicorn.v414.make_selection_webpage(full, selection, output='zspec_outliers_GBr.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'z_spec', 'abs(dz)'])
    
    #### stars
    imj = -2.5*np.log10(full['f_F814W']/full['f_F125W'])
    full.add_column(Column(data=imj, name='i-J', format='%.3f'))
    jmh = -2.5*np.log10(full['f_F125W']/full['f_F140W'])
    full.add_column(Column(data=jmh, name='J-H', format='%.3f'))

    selection = (full['z_max_grism'] > 0) & (full['star_flag'] == 1) & (full['hmag'] < 24) & (full['hmag'] > 21.5) & np.isfinite(imj) #& (imj > 1) 

    selection = (full['z_max_grism'] > 0) & (full['star_flag'] == 1) & (full['hmag'] < 24) #& (jmh > 1)

    unicorn.v414.make_selection_webpage(full, selection, output='red_stars_GBr.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'i-J', 'J-H'])
    """
    from threedhst import catIO
    PATH = '/Library/WebServer/Documents/P/GRISM_v4.1.4'
    full = catIO.Table('%s/RELEASE/3dhst.v4.1.4.full.v1.fits' %(PATH))
    return full
    
def select_list(full, list_file='/tmp/objects_list'):
    """ 
    Make a boolean selection that returns true for objects listed 
    in an ASCII file
    
    reload(unicorn.v414)
    selection = unicorn.v414.select_list(full, list_file='/tmp/objects_list')
    
    """
    objects = list(np.loadtxt(list_file, dtype=str))
    selection = full['id'] < 0
    for object in objects:
        selection |= (full['spec_id'] == object)
    
    return selection
    
def parse_inspection(file='Inspect/inspect_3dhst_ivo_8.1.info.fits', 
    columns=['contam','comment'], output=''):
    
    import os
    import copy
    
    import numpy as np
    from astropy.table import Table, Column
    from threedhst import catIO
    
    if not output:
        output = os.path.split(file)[-1].replace('fits','html')
                
    inspect = catIO.Table(file)
    sub = np.where((inspect['comment'] != '---') | (inspect['contam'] > 0) | (inspect['zp'] > 0))[0]
    s_inspect = inspect[sub]
    
    BASE = 'http://unicorn.astro.yale.edu/P/GRISM_v4.1.4/HTML/PNG'
    
    zfit, zfit2, linefit, thumb = [], [], [], []
    for filename in s_inspect['images']:
        id = filename.split('.new_zfit.png')[0]
        field = id.split('-')[0]
        pointing = id.split('-G141')[0]
        objid = id.split('G141_')[1]
        pngdir = '%s/%s/%s/' %(BASE, field, pointing)
        zfit.append('<img src=%s/%s.new_zfit.png height=200>' %(pngdir, id))
        zfit2.append('<img src=%s/%s.new_zfit.2D.png height=200>' %(pngdir, id))
        linefit.append('<img src=%s/%s.linefit.png height=200>' %(pngdir, id))
        thumb.append('<img src=%s/RGB/%s_%s_vJH_6.png height=180>' %(BASE, field, objid))

    s_inspect.add_column(Column(data=np.array(zfit), name='zfitpng'))
    s_inspect.add_column(Column(data=np.array(zfit2), name='zfit2png'))
    s_inspect.add_column(Column(data=np.array(linefit), name='linefit'))
    s_inspect.add_column(Column(data=np.array(thumb), name='thumb'))
    
    show_cols = copy.copy(columns)
    show_cols.extend(['thumb', 'zfitpng', 'zfit2png', 'linefit'])
    show = s_inspect[show_cols]
        
    show.write(output, format='ascii.html')
    os.system('perl -pi -e "s/\&lt;/</g" %s' %(output))
    os.system('perl -pi -e "s/\&gt;/>/g" %s' %(output))
    
    ##### Modify HTML table
    
    sorter = """
    <link rel="stylesheet" href="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/table_style.css" type="text/css" id="" media="print, projection, screen" /> 

    <script type="text/javascript" src="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/jquery-1.4.2.min.js"></script>

    <script type="text/javascript" src="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/jquery.tablesorter.min.js"></script>

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
    
    """ %(len(show_cols)-3, len(show_cols)-2, len(show_cols)-1)
    
    lines = open(output).readlines()
    thead_start = True
    for i, line in enumerate(lines):
        if "<head>" in line:
            lines.insert(i+1, sorter)
        
        if "<table>" in line:
            lines[i] = " <table id=\"myTable\" cellspacing=\"1\" class=\"tablesorter\">\n"
        
        if ("<th>" in line) & thead_start:
            lines[i-1] = " <thead>\n"
            thead_start = False
            j = i
            while "<th>" in lines[j]:
                j += 1
            
            lines[j] = " </thead>\n <tbody>\n"
        
        if "</table>" in lines[i]:
            lines[i] = "    </tbody>\n  </thead>\n"
        
    fp = open(output,'w')
    fp.writelines(lines)
    fp.close()
    
def random_objects(N=200):
    full = catIO.Table('../../3dhst.v4.1.4.full.v1.fits')
    has_spec = full['z_max_grism'] > 0.01
    objects = full['spec_id'][has_spec]
    idx = np.argsort(np.random.normal(size=has_spec.sum()))
    
    for i in range(N):
        obj = objects[idx][i]
        os.system('ln -s /Volumes/WD5_ESOMAC/3DHST/Spectra/Release/v4.1.4/*4/*/ZFIT/PNG/%s*png .' %(obj))
      
  
        
def make_selection_webpage(full, selection, output='test.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism'], local=False):
    """
    Make an HTML table with image links
    
    full = catIO.Table('3dhst.v4.1.4.full.v1.fits')

    cd "/Users/brammer/3DHST/Spectra/Release/v4.1.4/Selection"
    full = catIO.Table('../3dhst.v4.1.4.full.v1.fits')
    
    ### z_spec outliers
    dz = (full['z_max_grism'] - full['z_spec'])/(1+full['z_spec'])
    selection = (np.abs(dz) > 0.1) & (full['z_spec'] > 0.5) & (full['z_max_grism'] > 0.4)
    unicorn.v414.make_selection_webpage(full, selection, output='zspec_outliers.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'z_spec'], local=True)

    ### z_phot outliers
    dz = (full['z_max_grism'] - full['z_peak'])/(1+full['z_peak'])
    selection = (np.abs(dz) > 0.1) & (full['z_peak'] > 0.5) & (full['z_max_grism'] > 0.4) & (full['hmag'] < 24) & (full['star_flag'] == 0) & (full['near_star'] == 0)
    unicorn.v414.make_selection_webpage(full, selection, output='zpeak_outliers.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'z_peak'], local=True)
    
    ### H-alpha and [OIII]
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['HALPHA_EQW']/full['HALPHA_EQW_ERR'] > 3) &  (full['OIII_EQW']/full['OIII_EQW_ERR'] > 3)
    unicorn.v414.make_selection_webpage(full, selection, output='test.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'HALPHA_EQW', 'OIII_EQW'])
    
    selection = (full['z_max_grism'] > 1.5) & (full['use_all'] > 0) & (full['hmag'] < 22) & (full['lmass'] > 10.9) 
    unicorn.v414.make_selection_webpage(full, selection, output='test.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism', 'lmass', 'OIII_EQW'])
    
    
    ### Example, high EW OIII
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['OIII_EQW'] > 2000) & (full['OIII_EQW']/full['OIII_EQW_ERR'] > 3)
    
    unicorn.v414.make_selection_webpage(full, selection, output='test.html')
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['HALPHA_EQW'] > 1000) & (full['HALPHA_EQW']/full['HALPHA_EQW_ERR'] > 3)
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 22) & (full['HALPHA_FLUX']/full['HALPHA_FLUX_ERR'] > 3) & (full['OIII_FLUX']/full['OIII_FLUX_ERR'] > 3)

    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['HALPHA_FLUX']/full['HALPHA_FLUX_ERR'] > 3) 
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['OIIIX_FLUX']/full['OIIIX_FLUX_ERR'] > 2) 
    
    
    selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 24) & (full['z_peak'] > 3.5)
    
    
    """
    import os
    import copy
    
    import numpy as np
    from astropy.table import Table, Column
    
    #### e.g.
    #### selection = (full['z_max_grism'] > 0) & (full['use_all'] > 0) & (full['hmag'] < 22) & (full['OIII_EQW'] > 500)
    
    sub = full[selection]
    
    if local:
        BASE = 'file:///Users/brammer/3DHST/Spectra/Release/v4.1.4/'
    else:
        BASE = 'http://unicorn.astro.yale.edu/P/GRISM_v4.1.4/HTML/PNG'
    
    zfit, zfit2, linefit, thumb = [], [], [], []
    for id in sub['spec_id']:
        field = id.split('-')[0]
        pointing = id.split('-G141')[0]
        objid = id.split('G141_')[1]
        if local:
            zfitdir = '%s/%s-wfc3-spectra_v4.1.4/%s/ZFIT/PNG' %(BASE, field, pointing)
            lfitdir = '%s/%s-wfc3-spectra_v4.1.4/%s/LINE/PNG' %(BASE, field, pointing)
        else:
            zfitdir = '%s/%s/%s/' %(BASE, field, pointing)
            lfitdir = '%s/%s/%s/' %(BASE, field, pointing)
            
        zfit.append('<img src=%s/%s.new_zfit.png height=200>' %(zfitdir, id))
        zfit2.append('<img src=%s/%s.new_zfit.2D.png height=200>' %(zfitdir, id))
        linefit.append('<img src=%s/%s.linefit.png height=200>' %(lfitdir, id))
        thumb.append('<img src=%s/RGB/%s_%s_vJH_6.png height=180>' %(BASE, field, objid))
        
    sub.add_column(Column(data=np.array(zfit), name='zfitpng'))
    sub.add_column(Column(data=np.array(zfit2), name='zfit2png'))
    sub.add_column(Column(data=np.array(linefit), name='linefit'))
    sub.add_column(Column(data=np.array(thumb), name='thumb'))
    
    
    show_cols = copy.copy(columns)
    show_cols.extend(['thumb', 'zfitpng', 'zfit2png', 'linefit'])
    show = sub[show_cols]
    
    show['ra'].format = '%.5f'
    show['dec'].format = '%.5f'
    show['hmag'].format = '%4.2f'
    show['z_max_grism'].format = '%6.4f'
    
    show.write(output, format='ascii.html')
    os.system('perl -pi -e "s/\&lt;/</g" %s' %(output))
    os.system('perl -pi -e "s/\&gt;/>/g" %s' %(output))
    
    ##### Modify HTML table
    
    sorter = """
    <link rel="stylesheet" href="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/table_style.css" type="text/css" id="" media="print, projection, screen" /> 

    <script type="text/javascript" src="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/jquery-1.4.2.min.js"></script>

    <script type="text/javascript" src="http://monoceros.astro.yale.edu/RELEASE_V3.0/Spectra/UDF/Web/scripts/jquery.tablesorter.min.js"></script>

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
    
    """ %(len(show_cols)-3, len(show_cols)-2, len(show_cols)-1)
    
    lines = open(output).readlines()
    thead_start = True
    for i, line in enumerate(lines):
        if "<head>" in line:
            lines.insert(i+1, sorter)
        
        if "<table>" in line:
            lines[i] = " <table id=\"myTable\" cellspacing=\"1\" class=\"tablesorter\">\n"
        
        if ("<th>" in line) & thead_start:
            lines[i-1] = " <thead>\n"
            thead_start = False
            j = i
            while "<th>" in lines[j]:
                j += 1
            
            lines[j] = " </thead>\n <tbody>\n"
        
        if "</table>" in lines[i]:
            lines[i] = "    </tbody>\n  </thead>\n"
        
    fp = open(output,'w')
    fp.writelines(lines)
    fp.close()
   
def web_duplicates():
    
    #### add 0000 column
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        dups_list = np.loadtxt('%s-WFC3_v4.1.4_SPECTRA/%s-catalogs/%s.duplicates_zfit.v4.1.4.dat' %(field.upper(), field, field), dtype=str)
        sh = dups_list.shape
        idx = np.arange(sh[0])
        #
        zfit = catIO.Table('%s-WFC3_v4.1.4_SPECTRA/%s-catalogs/%s.new_zfit.v4.1.4.dat' %(field.upper(), field, field))
        #
        ok = dups_list[:,2] != '00000'
        for i in idx[ok]:
            print dups_list[i,2]
            objid = dups_list[i,2].split('-G141_')[1]
            dups = dups_list[i,1:][dups_list[i,1:] != '00000']
            fig = unicorn.plotting.plot_init(xs=5, square=True)
            ax = fig.add_subplot(111)
            
            zmin, zmax = 4, 0
            for j, dup in enumerate(dups):
                ix = zfit['spec_id'] == dup
                ran = np.array([zfit['l95'][ix][0], zfit['u95'][ix][0]])
                ax.scatter(zfit['z_max_grism'][ix][0], 0.5+j, color='black', marker='o')
                ax.fill_between(ran, [j,j],[j+1,j+1], color='black', alpha=0.5)
                ax.text(0.02, j/7.+0.5/7., dup, ha='left', va='center', fontsize=7, transform=ax.transAxes)
                zmin = np.minimum(zmin, ran[0])
                zmax = np.maximum(zmax, ran[1])
            #
            ax.set_ylim([0,7]); ax.set_xlabel('z')
            if zfit['z_spec'][ix][0] > 0:
                ax.plot(zfit['z_spec'][ix][0]*np.ones(2), [0,7], color='red')
                zmin = np.minimum(zmin, zfit['z_spec'][ix][0])
                zmax = np.maximum(zmax, zfit['z_spec'][ix][0])
            #
            ax.set_xlim(zmin-0.1, zmax+0.1)
            fig.tight_layout()
            unicorn.plotting.savefig(fig, 'Duplicates/%s_%s.dup.png' %(field, objid))
            plt.close()
        
        #
        from astropy.table import Table, Column    
        tab = Table()
        for column in ['id', 'Ndup', 'zavg', 'zlo', 'zhi']:
            tab.add_column(Column(name=column, dtype=float))
        
        tab.add_column(Column(name='dupPNG', dtype=str))
        for i in range(7):
            tab.add_column(Column(name='z%d' %(i+1), dtype=float))
        
        for i in range(7):
            tab.add_column(Column(name='id%d' %(i+1), dtype=float))
            tab.add_column(Column(name='2D_%d' %(i+1), dtype=str))
            
        
            
    for field, add in zip(['aegis', 'cosmos', 'goodsn', 'goodss', 'uds'], [3,3,4,0,3]):
        
        ix, id1, id2, id3, id4 = np.loadtxt('3dhst.duplicates_zfit.v4.1.4.dat', unpack=True, dtype=str)
    
        dups = catIO.Table('3dhst.duplicates_zfit.v4.1.4.dat')
    
        all_zfit = catIO.Table('3dhst.new_zfit_full.v4.1.4.dat')
    
        ok = dups['id']
    
def random_assignments():
    """
    Assign spectra to each person that will do inspections
    
    dirs=`ls |grep p |grep _ |grep -v unix |grep -v tar`
    for dir in $dirs; do
        for ext in linefit tilt 1D; do
            echo $dir $ext
            rm ${dir}/*${ext}.png
        done
    done
    
    for dir in $dirs; do
        echo $dir
        rm ${dir}/*[0-9].2D.png
    done
    
    """
    files=glob.glob('p[0-9]*unix')
    for file in files:
        dir=file.split('.unix')[0]
        list = open(file).readlines()
        for line in list:
            obj = line.split('/')[3].split('.')[0]
            field = obj.split('-')[0]
            id = obj.split('_')[1]
            rgb = 'RGB/%s_%s_vJH_6.png' %(field, id)
            print 'cp %s %s' %(rgb, dir)
            os.system('cp %s %s' %(rgb, dir))
            
def redshift_plot():
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1)
    
    ok_zsp = ok & (full['z_spec'] > 0.)
    
    line = ((full['OIII_EQW']/full['OIII_EQW_ERR'] > 1) | (full['OII_EQW']/full['OII_EQW_ERR'] > 1) | (full['HALPHA_EQW']/full['HALPHA_EQW_ERR'] > 1))
    
    dz = (full['z_max_grism']-full['z_spec'])/(1+full['z_spec'])
    
    plt.ioff()
    
    fig = unicorn.plotting.plot_init(aspect=0.55, square=True, xs=7, use_tex=True)
    
    ax = fig.add_subplot(121)
    ax.scatter(full['z_spec'][ok_zsp & ~line], full['z_max_grism'][ok_zsp & ~line], alpha=0.2, color='black')
    ax.scatter(full['z_spec'][ok_zsp & line], full['z_max_grism'][ok_zsp & line], alpha=0.2, color='red')
    ax.set_xlim(0,3.4); ax.set_ylim(0,3.4)
    ax.set_xlabel(r'$z_\mathrm{spec}$'); ax.set_ylabel(r'$z_\mathrm{grism}$')
    
    ax = fig.add_subplot(122)
    ax.hist(dz[ok_zsp & ~line], range=[-0.03, 0.03], bins=100, color='black', alpha=0.3, normed=True, label='NO lines, $\sigma$=%0.3f' %(threedhst.utils.nmad(dz[ok_zsp & ~line])))
    ax.hist(dz[ok_zsp & line], range=[-0.03, 0.03], bins=100, color='red', alpha=0.3, normed=True, label='Lines, $\sigma$=%0.3f' %(threedhst.utils.nmad(dz[ok_zsp & line])))
    ax.legend(loc='upper right', fontsize=8, ncol=1)
    ax.set_xlabel(r'$\Delta z/(1+z)$'); ax.set_ylabel(r'$N$')
    ax.set_xlim(-0.03, 0.03)
    
    fig.tight_layout()
    
    unicorn.plotting.savefig(fig, 'master_deltaz.pdf')
#

fakeLine = [1.4e4/(1+0.85)][::-1]
def test_line_fit():
    """
    Add a strong HeI line to an observed spectrum and try to recover the line 
    flux.
    """
    import pysynphot as s
    import astropy.io.fits as pyfits
    import threedhst
    from threedhst import catIO
    import unicorn
    import unicorn.interlace_test as test
    import os
    
    hizels = catIO.Table('Orig/match_0.84.cat')
    hizels = catIO.Table('Orig/match_1.47.cat')
    
    gris = test.SimultaneousFit('Orig/cosmos-23-G141_11986', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-02-G141_03882', fast=False, lowz_thresh=0.01)
    #gris = test.SimultaneousFit('Orig/cosmos-28-G141_23070', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-28-G141_23019', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-02-G141_04752', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-01-G141_19894', fast=False, lowz_thresh=0.01)
    #gris = test.SimultaneousFit('Orig/cosmos-27-G141_19894', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-05-G141_07983', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/uds-17-G141_19754', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/uds-24-G141_31562', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/uds-13-G141_27783', fast=False, lowz_thresh=0.01)
    
    id=int(gris.root.split('_')[-1])
    hizels_flux = 10**(hizels['logFlux'][hizels['id'] == id][0])
    
    lines = [5877.2] # HeI
    lines = [6718.29, 6732.67] # SII
    #lines = [6564.61] # Ha
    lines = unicorn.interlace_fit.fakeLine
    spec = s.FlatSpectrum(0, fluxunits='flam')
    for l0 in lines:
        #spec += s.GaussianSource(10.e-16/len(lines), l0*(1+gris.z_max_spec), 2.)
        spec += s.GaussianSource(hizels_flux, l0*(1+0.85), 2.)
    
    gris.twod.compute_model(spec.wave, spec.flux/1.e-17/gris.twod.total_flux)
    threedhst.showMessage('%f %f.  HiZELS:%.1f' %(gris.twod.model.sum(), s.Observation(spec, s.ObsBandpass('wfc3,ir,g141')).countrate(), hizels_flux/1.e-17))
    
    gris.twod.im['SCI'].data += gris.twod.model
    gris.twod.im.writeto(os.path.basename(gris.twod.file), clobber=True)

    os.system('cp %s .' %(gris.twod.file.replace('2D.','1D.')))
    os.system('cp %s .' %(gris.twod.file.replace('2D.','new_zfit.pz.')))
    os.system('cp %s .' %(gris.twod.file.replace('2D.','new_zfit.')))
    twod = unicorn.reduce.Interlace2D(os.path.basename(gris.twod.file))
    twod.oned_spectrum()
    gris = test.SimultaneousFit(gris.root, fast=False, lowz_thresh=0.01)
    
    gris.new_fit_free_emlines(NTHREADS=1, NSTEP=400)
    
    #### Try aperture flux
    gris = test.SimultaneousFit('Orig/uds-13-G141_27783', fast=False, lowz_thresh=0.01)
    gris = test.SimultaneousFit('Orig/cosmos-02-G141_03882', fast=False, lowz_thresh=0.01)
    f = pyfits.open(gris.twod.file.replace('.2D', '.new_zfit'))
    chain = unicorn.interlace_fit.emceeChain(file=os.path.basename(gris.twod.file.replace('.2D', '.linefit')))
    
    import scipy.ndimage as nd
    kern = nd.maximum_filter((gris.twod.im['DSEG'].data == gris.twod.id)*1., size=3)
    #kern /= kern.sum()
        
    var2D = gris.twod.im['WHT'].data*1.
    var2D[gris.twod.im['WHT'].data > 1] = 0.
    var = nd.convolve((var2D/(gris.twod.im['SENS'].data*4.))**2, kern)
    
    flux = nd.convolve((gris.twod.cleaned - f['CONT2D'].data*1)*(var2D > 0)/(gris.twod.im['SENS'].data*4.), kern)
    xf, yf = gris.twod.trace_extract(flux, width=0, dy=0)
    plt.plot(xf, yf * threedhst.utils.diff(xf), color='black', alpha=0.7)

    flux = nd.convolve((f['LINE2D'].data+f['CONT2D'].data*0)*(var2D > 0)/(gris.twod.im['SENS'].data*4.), kern)
    xf, yf = gris.twod.trace_extract(flux, width=0, dy=0)
    plt.plot(xf, yf * threedhst.utils.diff(xf), color='red', alpha=0.5)
    
    xvar, yvar = gris.twod.trace_extract(var, width=0, dy=0)
    plt.plot(xvar, np.sqrt(yvar) * threedhst.utils.diff(xf), color='orange', alpha=0.5)
    
    plt.fill_between(xf, yf*0.+np.exp(chain.stats['Ha']['q16'])*(1+gris.z_max_spec), yf*0.+np.exp(chain.stats['Ha']['q84'])*(1+gris.z_max_spec), color='blue', linewidth=2, alpha=0.5)
    
    plt.ylim(-5,30)
    
def check_twod_flux():
    """
    Do full PySynphot test to check line fits
    """
    import numpy as np
    import matplotlib.pyplot as plt

    import pysynphot as s
    import unicorn
    
    flux, z = 31.3e-17, 2.1947
    #### 1D model
    spec = s.GaussianSource(flux*1./3.98, 4959*(1+z), 10) + s.GaussianSource(flux*2.98/3.98, 5007*(1+z), 10)
    check_flux = np.trapz(spec.flux, spec.wave)
    #### Is this what we think
    print 'Model flux=%.1e (ratio=%.2f)' %(np.trapz(spec.flux, spec.wave), np.trapz(spec.flux, spec.wave)/flux)
    
    #### model 2D spectrum
    twod = unicorn.reduce.Interlace2D('cosmos-23-G141_10280.2D.fits')
    twod.compute_model(spec.wave, spec.flux/1.e-17/twod.total_flux)
    sens_scaled = twod.im['SENS'].data*4 # this will be fixed eventually
    flux_density = twod.model / sens_scaled
    dlam = np.diff(twod.im['WAVE'].data) # this missing an element
    dlam = np.append(dlam[0], dlam)      # add the element for the first entry
    integrated = np.sum(flux_density*dlam)*1.e-17
    
    ### Check 2D integrated flux
    print '2D flux = %.1e (ratio=%.2f)' %(integrated, integrated/flux)
    
    ### Check that looks right in the image
    fig = plt.figure(figsize=[4,4])
    ax = fig.add_subplot(311); ax.imshow(0-twod.model, vmin=-0.1, vmax=0.01, interpolation='Nearest')
    ax.text(0.95, 0.2, 'Model', transform=ax.transAxes, color='red', ha='right', va='top', fontsize=12)
    ax.text(0.05, 0.9, twod.file.split('.2D')[0], ha='left', va='top', color='black', fontsize=13, transform=ax.transAxes)
    ax.text(0.05, 0.15, r'[OIII] flux = %.1e erg/s/cm2/A' %(flux), ha='left', va='top', color='black', fontsize=10, transform=ax.transAxes)
    
    ax = fig.add_subplot(312); ax.imshow(0-twod.cleaned, vmin=-0.1, vmax=0.01, interpolation='Nearest')
    ax.text(0.95, 0.2, 'Spectrum', transform=ax.transAxes, color='red', ha='right', va='top', fontsize=12)

    ax = fig.add_subplot(313); ax.imshow(0-(twod.cleaned - twod.model), vmin=-0.1, vmax=0.01, interpolation='Nearest')
    ax.text(0.95, 0.2, 'Diff.', transform=ax.transAxes, color='red', ha='right', va='top', fontsize=12)
    
    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xlim(207, 285)
        
    fig.tight_layout(pad=0.1)
    fig.savefig('line_flux_2D.pdf')
    
def check_cutoff_pz():
    """
    Some objects have grism p(z) that looks cut off with respect 
    to photoz
    """
    gris = test.SimultaneousFit('goodss-29-G141_41381', fast=False, lowz_thresh=0.01)

def join_v415_tables():
    import astropy.table
    import os
    import numpy as np
    
    tt = astropy.table.Table()
    
    os.chdir("/Users/brammer/3DHST/Spectra/Release/v4.1.5/FullRelease")
    
    linefit = []
    zfit = []
    rf = []
    sfr = []
    fout = []
    
    dups = []
    
    for ifield, field in enumerate(['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']):
        print field
        
        out = np.loadtxt('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.duplicates_2d.dat' %(field, field), dtype='a', unpack=True)
        pointings = []
        npoint = []
        for i in range(len(out[0])):
            pi = []
            npi = 0
            for j in range(1,len(out)):
                if 'G141' in out[j,i]:
                    npi+=1
                    pi.append(out[j,i].split('-')[1])
            
            npoint.append(npi)    
            if len(pi) == 0:
                pointings.append('--')
            else:
                pointings.append(','.join(pi))
        
        ti = astropy.table.Table()
        ti.add_column(astropy.table.Column(data=pointings, name='pointings'))
        ti.add_column(astropy.table.Column(data=npoint, name='npoint'))
        dups.append(ti)
             
        t = tt.read('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.linefit.linematched.fits' %(field, field))
        linefit.append(t)
        
        t = tt.read('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.zfit.linematched.fits' %(field, field))
        t.add_column(astropy.table.Column(data=[field]*len(t), name='field'), index=2)
        t.add_column(astropy.table.Column(data=[ifield+1]*len(t), name='ifield'), index=2)
        zfit.append(t)
        
        
        t = tt.read('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.zbest.rf' %(field, field), format='ascii.commented_header')
        rf.append(t)
        
        t = tt.read('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.zbest.sfr' %(field, field), format='ascii.commented_header')
        sfr.append(t)
        
        #t = tt.read('%s_3dhst_v4.1.5_catalogs/%s_3dhst.v4.1.5.zbest.fout' %(field, field), format='ascii.commented_header')
        t = tt.read('NewFAST_Oct8/%s_3dhst.v4.1.5.zbest.fout' %(field), format='ascii.commented_header')
        fout.append(t)
        
    #
    zfit_full = astropy.table.vstack(zfit)
    linefit_full = astropy.table.vstack(linefit)
    rf_full = astropy.table.vstack(rf)
    sfr_full = astropy.table.vstack(sfr)
    fout_full = astropy.table.vstack(fout)
    dups_full = astropy.table.vstack(dups)
    
    # count = np.arange(len(zfit_full), dtype=int)+1
    # ccol = astropy.table.Column(data=count, name='counter')
    # for tab in [linefit_full, zfit_full, rf_full, sfr_full, fout_full, dups_full]:
    #     tab.add_column(ccol)
     
    # linefit_full.rename_column('number', 'phot_id')
    # rf_full.rename_column('id', 'phot_id')
    # 
    # sfr_full.rename_column('flag', 'sfr_flag')
    # sfr_full.rename_column('id', 'phot_id')
    # sfr_full.rename_column('beta', 'sfr_beta')
    # 
    # fout_full.rename_column('id', 'phot_id')
    # fout_full.rename_column('z', 'z_best')
    
    #### rename columns
    for tab, name in zip([linefit_full, rf_full, sfr_full, fout_full], ['linefit', 'rf', 'sfr', 'fout']):
        for c in tab.keys():
            if c in zfit_full.keys():
                tab.rename_column(c, '%s_%s' %(name, c))
                
    full = astropy.table.hstack([zfit_full, linefit_full, rf_full, sfr_full, fout_full])
    full.add_column(dups_full['pointings'], index=4)
    full.add_column(dups_full['npoint'], index=5)
    
    phot = tt.read('../../v4.1.5/Catalogs/3dhst.v4.1.5.full.v1.fits')
    for c in ['y','x','dec','ra']:
        print c
        full.add_column(phot[c], index=1)
    
    full.write('3dhst.v4.1.5.master.fits')
    
    ## Add Galfit
    import astropy.table
    tt = catIO.Table('3dhst.v4.1.5.master.fits')
    
    for filt, band in zip(['f125w', 'f160w'], ['J', 'H']):
        gf = []
        for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
            gf.append(catIO.Table('Galfit/%s/%s_3dhst.v4.1_%s.galfit' %(field, field, filt)))
        
        gfull = astropy.table.vstack(gf)
        gfull.remove_columns(['NUMBER', 'RA', 'DEC'])
        for c in gfull.colnames:
            gfull.rename_column(c, 'galfit_%s_%s' %(band, c))
        
        tt = astropy.table.hstack([tt, gfull])
    
    tt.write('3dhst.v4.1.5.master.galfit.fits')
        
    #### Demo, plot UVJ diagram color-coded by H-alpha equivalent width
    import astropy.table
    tt = astropy.table.Table()
    full = tt.read('3dhst.v4.1.5.master.fits')
    
    UV = -2.5*np.log10(full['L153']/full['L155'])
    UB = -2.5*np.log10(full['L153']/full['L154'])
    VJ = -2.5*np.log10(full['L155']/full['L161'])
    #VJ = -2.5*np.log10(full['L155']/full['L163'])
    HK = -2.5*np.log10(full['L162']/full['L163'])
    
    #UV = -2.5*np.log10(full['L272']/full['L154']) # 2200 / B
    #VJ = -2.5*np.log10(full['L154']/full['L160']) # B / z
    
    #### H-alpha
    ok = (full['z_best_s'] > 0) & (full['z_max_grism'] > 0.72) & (full['z_max_grism'] < 1.58) & (full['Ha_SCALE'] != -99) & (full['lmass'] > 9.5)
    line = 'Ha'

    ok = (full['z_best_s'] > 0) & (full['z_max_grism'] > 1.25) & (full['z_max_grism'] < 2.3) & (full['OIII_SCALE'] != -99) & (full['lmass'] > 9.2)
    line = 'OIII'

    # ok = (full['z_best_s'] > 0) & (full['z_max_grism'] > 1.92) & (full['z_max_grism'] < 3.4) & (full['OII_SCALE'] != -99) & (full['lmass'] > 10)
    # line = 'OII'
    
    eqw = np.log10(full['%s_EQW' %(line)]/(1+full['z_max_grism']))
    eqw[(full['%s_EQW' %(line)] < 0) < 0] = -1
    vmi = [-1,3]
    
    eqw = np.log10(full['SII_FLUX']/full['Ha_FLUX'])
    eqw[(full['SII_FLUX'] < 0) < 0] = -3
    vmi = [-1.,0.0]
    ok &= (full['SII_FLUX']/full['SII_FLUX_ERR'] > 2) & (full['Ha_FLUX']/full['Ha_FLUX_ERR'] > 2)
    
    #eqw[(full['%s_EQW' %(line)] < 0) | (full['%s_EQW' %(line)]-full['%s_EQW_ERR' %(line)] < 0)] = -1
    
    ####
    plt.set_cmap('coolwarm_r'); plt.close()

    fig = plt.figure(figsize=[7.5,5]); ax = fig.add_subplot(111)
    
    ax.scatter(VJ[ok], UV[ok], c=eqw[ok], vmin=vmi[0],vmax=vmi[1], alpha=0.4)
    sc = ax.scatter(VJ[ok][0]-10, UV[ok][0]-10, c=eqw[ok][0], vmin=vmi[0],vmax=vmi[1], alpha=0.8)
    ax.set_xlim(-0.5,2.5); ax.set_ylim(-0.1,2.5)
    #ax.set_xlim(-0.5,2.7); ax.set_ylim(-0.1,5.) 
    ax.set_xlabel(r'$V-J$'); ax.set_ylabel(r'$U-V$')
    
    if False:
        ### Check H-K color test
        ok = (full['z_max_grism'] > 0.5) & (full['z_max_grism'] < 2.) & (full['lmass'] > 9.5)
        
        fig = plt.figure(figsize=[7.5,5]); ax = fig.add_subplot(111)
        
        sc = ax.scatter(VJ[ok], UB[ok], c=HK[ok], vmin=vm[0],vmax=vm[1], alpha=0.4)
        #sc = ax.scatter(VJ[ok][0]-10, UV[ok][0]-10, c=eqw[ok][0], vmin=vmi[0],vmax=vmi[1], alpha=0.8)
        ax.set_xlim(-0.5,2.5); ax.set_ylim(-0.1,2.5)
        ax.set_xlabel(r'$V-J$'); ax.set_ylabel(r'$U-V$')
        plt.colorbar(sc)
        
    cb = plt.colorbar(sc); cb.set_label(r'$\log_{10}$ %s EQW' %(line))
    ax.text(0.02,0.98, r'$%.1f < z < %.1f$' %(full['z_max_grism'][ok].min(), full['z_max_grism'][ok].max()) + '\nlog M/'+r'M$_\odot$'+' > 9.2\n' + r'$N$=%d' %(ok.sum()), backgroundcolor='white', ha='left', va='top', transform=ax.transAxes)
    
    fig.tight_layout(pad=0.1)
    
    fig.savefig('3dhst_v4.1.5_UVJ_%s_demo.png' %(line))
    
    #### Dust correction
    import research.dusty
    old = catIO.Table('../..//v4.1.4/3dhst.v4.1.4.full.v2.fits')
    Av_eazy = np.interp(old['dusty_total'], [0,1], [0,2])
    
    dUV, dVJ = research.dusty.get_calzetti_correction(UV, VJ, full['Av'])
    #dUV, dVJ = research.dusty.get_calzetti_correction(UV, VJ, Av_eazy)
    #plt.set_cmap('coolwarm_r'); plt.close()
    fig = plt.figure(figsize=[7.5,5]); ax = fig.add_subplot(111)
    
    ax.scatter((VJ-dVJ)[ok], (UV-dVJ)[ok], c=eqw[ok], vmin=vmi[0],vmax=vmi[1], alpha=0.4)
    sc = ax.scatter(VJ[ok][0]-10, UV[ok][0]-10, c=eqw[ok][0], vmin=vmi[0],vmax=vmi[1], alpha=0.8)
    ax.set_xlim(-0.5,2.5); ax.set_ylim(-0.1,2.5)
    ax.set_xlabel(r'$V-J$'); ax.set_ylabel(r'$U-V$')
    
    cb = plt.colorbar(sc); cb.set_label(r'$\log_{10}$ %s EQW' %(line))
    ax.text(0.02,0.98, r'$%.1f < z < %.1f$' %(full['z_max_grism'][ok].min(), full['z_max_grism'][ok].max()) + '\nlog M/'+r'M$_\odot$'+' > 9.2\n' + r'$N$=%d' %(ok.sum()), backgroundcolor='white', ha='left', va='top', transform=ax.transAxes)
    
    fig.tight_layout(pad=0.1)
    
    ### Hb/OIII
    ok = (full['z_best_s'] > 0) & (full['Hb_SCALE'] != -99) & (full['OIII_SCALE'] != -99) & (full['lmass'] > 9.2)
    
    ok &= (full['Hb_FLUX']/full['Hb_FLUX_ERR'] > 1)
    
    HbO3 = full['Hb_FLUX']/full['OIII_FLUX']

    line = 'OIII'
    eqw = np.log10(full['%s_EQW' %(line)]/(1+full['z_max_grism']))
    eqw[(full['%s_EQW' %(line)] < 0) < 0] = -1
    
    err = np.sqrt((full['Hb_FLUX_ERR']/full['Hb_FLUX'])**2 + (full['OIII_FLUX_ERR']/full['OIII_FLUX'])**2)
    
    plt.errorbar(full['lmass'][ok], 1/HbO3[ok], (1/HbO3*err)[ok], ecolor='0.6', alpha=0.5, linestyle='None', marker='None', zorder=-1)
    plt.scatter(full['lmass'][ok], 1/HbO3[ok], c=eqw[ok], vmin=-1, vmax=3, alpha=0.5, zorder=1)
    plt.xlim(9.2,11.5); plt.ylim(-0.5,20)

    #### R23
    R23 = (full['OII_FLUX']+full['OIII_FLUX'])/(full['Hb_FLUX'])
    err_R23 = np.sqrt((full['Hb_FLUX_ERR']/full['Hb_FLUX'])**2 + (full['OII_FLUX_ERR']**2+full['OIII_FLUX_ERR']**2)/(full['OII_FLUX']+full['OIII_FLUX'])**2)
    
    O32 = (full['OIII_FLUX']/full['OII_FLUX'])
    x = np.log10(R23)
    y = np.log10(O32)
    logOH = 12 - 4.944 + 0.767*x + 0.602*x**2 - y*(0.29+0.332*x-0.331*x**2)
    
    ok = (full['z_best_s'] > 0) & (full['Hb_SCALE'] != -99) & (full['OIII_SCALE'] != -99) & (full['lmass'] > 9.2) & (full['OII_SCALE'] != -99)
    ok &= (full['OII_FLUX']/full['OII_FLUX_ERR'] > 1) & (full['OIII_FLUX']/full['OIII_FLUX_ERR'] > 1) & (full['Hb_FLUX']/full['Hb_FLUX_ERR'] > 1) & (O32 > 0)

    o4363 = (full['OIIIx_FLUX']/full['OIIIx_FLUX_ERR'] > 2)
    
    plt.errorbar(full['lmass'][ok], R23[ok], (R23*err_R23)[ok], ecolor='0.6', alpha=0.3, linestyle='None', marker='None', zorder=-1)
    plt.scatter(full['lmass'][ok & ~o4363], R23[ok & ~o4363], c=eqw[ok & ~o4363], vmin=-1, vmax=3, alpha=0.5, zorder=1)
    plt.scatter(full['lmass'][ok & o4363], R23[ok & o4363], c=eqw[ok & o4363], vmin=-1, vmax=3, alpha=0.8, zorder=1, marker='s', s=60)
    plt.xlim(9.2,11.5); plt.ylim(-0.5,20)

    
    #### Compare EAZY masses
    dz = (full['z']-full['z_peak'])/(1+full['z_peak'])
    ok = (full['hmag'] < 25) & (full['use_phot'] == 1) & (np.abs(dz) < 0.03) & (full['lmass'] > 8) & np.isfinite(full['lmass']) & np.isfinite(full['eazy_mass'])

    plt.set_cmap('gray_r')
    plt.close(); astroML.plotting.scatter_contour(full['lmass'][ok], full['eazy_mass'][ok]-full['lmass'][ok], levels=10, threshold=10, plot_args={'linestyle':'None', 'marker':'o', 'alpha':0.1, 'color':'black'}, histogram2d_args={'bins':[30,30], 'range':[[8,12],[-1.5,1.5]]}, log_counts=False)

    so = np.argsort(full['lmass'][ok])[::-1]
    xm, ym, ys, n = threedhst.utils.runmed(full['lmass'][ok][so], (full['eazy_mass']-full['lmass'])[ok][so], NBIN=40, use_median=True, use_nmad=False)
    #plt.errorbar(xm, ym, ys, linestyle='None', marker='o', color='red', ecolor='red')
    plt.plot(xm, ym, color='red', alpha=0.6, linewidth=3)
    plt.plot(xm, ym+ys, color='red', alpha=0.5, linewidth=1)
    plt.plot(xm, ym-ys, color='red', alpha=0.5, linewidth=1)
    plt.fill_between(xm, ym+ys, ym-ys, color='red', alpha=0.1, linewidth=3)

    plt.ylim(-1.5,1.5); plt.xlim(8,12)
    plt.grid()
    plt.xlabel('FAST lmass'); plt.ylabel(r'$\Delta$Mass (EAZY $-$ FAST)')
    plt.tight_layout(pad=0.18)
    plt.text(11.6,-1.4, '3D-HST v4.1, H<25, log M > 8', ha='right', va='bottom', backgroundcolor='white')
    plt.savefig('compare_FAST_EAZY_masses.pdf')
    
    
    
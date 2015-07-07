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
    
    filters = []
    filters.extend(['153', '154', '155']) # Maiz UBV
    filters.extend(['161', '162', '163']) # 2MASS JHK
    filters.extend(['156', '157', '158', '159', '160']) # SDSS ugriz
    filters.extend(['135','136','137','138','139']) # Bessel UX BVRI
    filters.extend(['270','271','272','273','274','275']) # Tophat filters
    
    #### Add z_gris to catalog column to use as z_spec
    field='aegis'
    
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        
        zfit = catIO.Table('%s/%s-WFC3_v4.1.4_SPECTRA/%s-catalogs/%s.new_zfit.linematched.v4.1.4.dat' %(PATH_TO_NEW, field.upper(), field, field))
    
        old_cat = catIO.Table('%s/%s_3dhst.v4.1.cats/Catalog/%s_3dhst.v4.1.cat' %(PATH_TO_PHOT, field, field))
        old_zout = catIO.Table('%s/%s_3dhst.v4.1.cats/Eazy/%s_3dhst.v4.1.zout' %(PATH_TO_PHOT, field, field))
    
        z_best = old_zout['z_peak']
        z_type = np.ones(len(z_best), dtype=int)
        
        has_grism = zfit['z_max_grism'] > 0
        z_best[has_grism] = zfit['z_max_grism'][has_grism]
        z_type[has_grism] = 2
        
        z_best_column = catIO.table_base.Column(data=z_best, name='z_best')
        old_cat.add_column(z_best_column)

        z_type_column = catIO.table_base.Column(data=z_type, name='z_type')
        old_cat.add_column(z_type_column)
    
        old_cat.write(os.path.basename(old_cat.filename).replace('.cat', '.zbest.cat'), format='ascii.commented_header')
    
    
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
        
        cat = catIO.Table('%s_3dhst.v4.1.zbest.cat' %(field))
        rf0 = catIO.Table('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, filters[0]))
        for c in ['nfilt_fit', 'chi2_fit', 'L%s' %(filters[0])]:
            rf0.remove_column(c)
            
        rf0.rename_column('z', 'z_best')
        rf0.add_column(cat['z_type'], index=2)
        rf0.add_column(cat['z_spec'], index=3)
        
        for f in filters:
            print '%s, filter: %s' %(field, f)
            rf = catIO.Table('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, f))
            rf.rename_column('nfilt_fit', 'nfilt%s' %(f))
            rf0.add_columns([rf['L%s' %(f)], rf['nfilt%s' %(f)]])
        
        rf0.write('%s_3dhst.v4.1.4.zbest.rf' %(field), format='ascii.commented_header')
    
    #
    descrip = ['#\n']
    for f in filters:
        #print '%s, filter: %s' %(field, f)
        fp = open('OUTPUT/%s_3dhst.v4.1.zbest.%s.rf' %(field, f))
        for i in range(3):
            line = fp.readline()
        
        descrip.append(line)
    
    descrip.append('#\n# Rest-frame colors computed using templates in tweak_cosmos_v4.1/spectra.param\n#\n# z_type: 1 (phot), 2 (grism)\n#\n')
    
    for field in ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']:
        lines = open('%s_3dhst.v4.1.4.zbest.rf' %(field)).readlines()
        for line in descrip[::-1]:
            lines.insert(1, line)
            
        fp = open('%s_3dhst.v4.1.4.zbest.rf' %(field), 'w')
        fp.writelines(lines)
        fp.close()
        
def fit_zbest_RF_colors(field='aegis'):
    """
    Use z_grism as the redshift for fitting EAZY RF colors
    """
    import os
    import threedhst.eazyPy as eazy
    
    param = eazy.EazyParam('%s_3dhst.v4.1.param' %(field))
    
    param['FILTERS_RES'] = 'FILTER.RES.latest'
    param['CATALOG_FILE'] = '%s_3dhst.v4.1.zbest.cat' %(field)
    param['OUTPUT_DIRECTORY'] = 'OUTPUT'
    
    #### Use z_spec / z_best
    param['BINARY_OUTPUT'] = 'n'
    param['FIX_ZSPEC'] = 'y'
    
    #### Dusty template fit
    param['TEMPLATES_FILE'] = 'tweak_cosmos_v4.1/spectra.param'
    param['CACHE_FILE'] = '%s_3dhst.v4.1.zbest.tempfilt' %(field)
    param['MAIN_OUTPUT_FILE'] = field+'_3dhst.v4.1.zbest'
    
    #### Baseline
    param['REST_FILTERS'] = '---'
    param['RF_ERRORS'] = 'n'
    param.write('/tmp/zphot.%s.param' %(field))
    
    ### translate with z_best = z_spec
    lines = open('%s_3dhst.v4.1.translate' %(field)).readlines()
    lines.append('z_spec z_spec_old\nz_best z_spec\n')
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
    from astropy.table import Table, Column
    
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
    
    ### Galaxies
    npix = (img['FLUX'].data != 0).sum(axis=1)
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1)
    so = np.argsort(full['z_max_grism'][ok])
    skip = 15
    
    #### Continuum
    line_SN = npix < 0
    for line in ['OII', 'OIII', 'HALPHA']:
        line_SN |= (full['%s_FLUX' %(line)]/full['%s_FLUX_ERR' %(line)] > 2) & (full['%s_FLUX' %(line)] > 0)
    
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.15) & (full['star_flag'] != 1) & line_SN
    so = np.argsort(full['z_max_grism'][ok])
    skip = 5
    
    ## sort by stellar mass
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 2.3)
    
    ### High-z
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 1.7) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 3.2)
    ### lowz
    ok = (npix > 80) & ((img['CONTAM'].data/img['FLUX'].data > 4).sum(axis=1) < 100) & (full['hmag'] < 24) & (full['z_max_grism'] > 0.18) & (full['star_flag'] != 1) &  (full['z_max_grism'] < 0.9)

    so = np.argsort(full['lmass'][ok])
    skip = 21
    
    #### stars
    imj = -2.5*np.log10(full['f_F814W']/full['f_F125W'])
    imj = -2.5*np.log10(full['f_F125W']/full['f_F140W'])
    ok = (npix > 150) & ((img['CONTAM'].data/img['FLUX'].data > 3).sum(axis=1) < 60) & (full['hmag'] < 24) & (full['star_flag'] == 1) & np.isfinite(imj)
    so = np.argsort(imj[ok])
    skip = 5
    #avg = np.median(smooth[0:100,:], axis=0)
    
    sub = (img['FLUX'].data-img['CONTAM'].data)[ok,:][so,:]
    scl = 10**(0.4*(full['hmag'][ok][so]-21))
    sub = (sub.T*scl).T
    
    smooth = sub*1.
    for i in range(ok.sum())[::skip]:
        smooth[i,:] = np.median(sub[i:i+skip,:], axis=0)
    
    smooth = smooth[::skip]
    avg = np.median(smooth, axis=0)
    
    zi = full['z_max_grism'][ok][so][::skip]
    lrest = np.arange(2500, 1.3e4, 10)
    lrest = np.arange(2500, 1.4e4, 10)
    sh = smooth.shape
    spec_rest = np.zeros((sh[0], len(lrest)))
    for i in range(sh[0]):
        spec_rest[i,:] = np.interp(lrest, lran[30:]/(1+zi[i]), smooth[i,30:]/avg[30:], left=0, right=0)
    
    ### full rest, for when not sorted by redshift
    zzi = full['z_max_grism'][ok][so]
    sh = sub.shape
    spec_rest_full = np.zeros((sh[0], len(lrest)))
    for i in range(sh[0]):
        spec_rest_full[i,:] = np.interp(lrest, lran[30:-5]/(1+zzi[i]), sub[i,30:-5]/avg[30:-5], left=0, right=0)
    
    spec_rest_full[~np.isfinite(spec_rest_full)] = 0
    smooth_full = spec_rest_full*1.
    for i in range(ok.sum())[::skip]:
        print i
        sub_region = np.ma.masked_array(spec_rest_full[i:i+skip,:], mask=spec_rest_full[i:i+skip,:] == 0)
        smooth_full[i,:] = np.ma.median(sub_region, axis=0)
    
    smooth_full = smooth_full[::skip]
    
    ### Stack
    spec_rest[~np.isfinite(spec_rest)] = 0
    sum_rest = np.sum(spec_rest, axis=0)
    N = np.sum(spec_rest != 0, axis=0)
    
    spec_rest[spec_rest == 0] = 100
    
    import seaborn as sns
    sns.set(style="ticks")
    plt.ioff()
    fig = unicorn.plotting.plot_init(aspect=sh[0]*1./sh[1]*0.5, xs=6, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)
    
    ax = fig.add_subplot(111)
    ax.imshow(255-unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
    
    zpix = np.arange(0.25,3.1,0.25)
    zpixstr = list(zpix)
    for i in range(len(zpix))[::2]:
        zpixstr[i] = ''
    
    ypix = np.interp(zpix, zi, np.arange(sh[0]))
    ax.set_yticks(ypix); ax.set_yticklabels(zpixstr)
    ax.set_ylim(0, sh[0])
    
    ax2 = ax.twinx()
    #ax2.imshow(unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
    ynticks = np.arange(0, ok.sum(), 5000)
    ynv = ynticks/skip
    ax2.set_yticks(ynv); ax2.set_yticklabels(ynticks)
    ax2.set_ylim(0, sh[0])
    ax2.set_ylabel(r'$N$')
    ax.set_ylabel(r'$z$')
    
    xlam = [0.25,0.5,0.75,1,1.25]
    xpix = np.interp(xlam, lrest/1.e4, np.arange(len(lrest)))
    ax.set_xticks(xpix); ax.set_xticklabels(xlam)
    ax2.set_xticks(xpix); ax2.set_xticklabels(xlam)
    ax.set_xlabel(r'$\lambda_\mathrm{rest}$ / $\mu\mathrm{m}$')
    
    fig.tight_layout()
    
    unicorn.plotting.savefig(fig, '3dhst_allspecs.pdf')
    
    
    ##### Horizontal
    import seaborn as sns
    sns.set(style="ticks")
    plt.ioff()
    fig = unicorn.plotting.plot_init(aspect=1./(sh[0]*1./sh[1]*0.5), xs=10, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)
    
    ax = fig.add_subplot(111)
    ax.imshow(255-unicorn.candels.clipLog(spec_rest.T, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
    
    zpix = np.arange(0.25,3.1,0.25)
    zpixstr = list(zpix)
    for i in range(len(zpix))[::2]:
        zpixstr[i] = ''
    
    ypix = np.interp(zpix, zi, np.arange(sh[0]))
    ax.set_xticks(ypix); ax.set_xticklabels(zpixstr)
    ax.set_xlim(0, sh[0])
    
    ax2 = ax.twiny()
    #ax2.imshow(unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
    ynticks = np.arange(0, ok.sum(), 5000)
    ynv = ynticks/skip
    ax2.set_xticks(ynv); ax2.set_xticklabels(ynticks)
    ax2.set_xlim(0, sh[0])
    ax2.set_xlabel(r'$N$')
    ax.set_xlabel(r'$z$')
    
    xlam = [0.25,0.5,0.75,1,1.25]
    xpix = np.interp(xlam, lrest/1.e4, np.arange(len(lrest)))
    ax.set_yticks(xpix); ax.set_yticklabels(xlam)
    ax2.set_yticks(xpix); ax2.set_yticklabels(xlam)
    ax.set_ylabel(r'$\lambda_\mathrm{rest}$ / $\mu\mathrm{m}$')
    
    fig.tight_layout()
    
    unicorn.plotting.savefig(fig, '3dhst_allspecs_horiz.pdf')
    
    #### Sort stellar mass
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as font_manager

    path = '/Library/Fonts/Microsoft/Calibri.ttf'
    prop = font_manager.FontProperties(fname=path)
    mpl.rcParams['font.family'] = prop.get_name()
    mpl.rcParams['font.style'] = 'normal'
    mpl.rcParams['font.weight'] = 'normal'
    
    sh = smooth_full.shape
    import seaborn as sns
    sns.set(style="ticks")
    plt.ioff()
    smooth_full[smooth_full == 0] = 100
    fig = unicorn.plotting.plot_init(aspect=1./1.3, xs=6, top=0.01, right=0.1, bottom=0.1, left=0.1, square=True, use_tex=True)
    
    ax = fig.add_subplot(111)
    ax.imshow(255-unicorn.candels.clipLog(smooth_full.T, lexp=1e5, cmap=[6.47, 0.88235], scale=[0,5]), origin='lower', aspect='auto')
    
    mass = full['lmass'][ok][so][::skip]
    zpix = np.arange(9,11.1,0.5)
    zpixstr = list(zpix)
    
    ypix = np.interp(zpix, mass, np.arange(sh[0]))
    ax.set_xticks(ypix); ax.set_xticklabels(zpixstr)
    ax.set_xlim(0, sh[0])
    
    ax2 = ax.twiny()
    #ax2.imshow(unicorn.candels.clipLog(spec_rest, lexp=1e5, cmap=[6, 0.921569], scale=[0.6,1.5]), origin='lower', aspect='auto')
    ynticks = np.arange(0, ok.sum(), 2000)
    ynv = ynticks/skip
    ax2.set_xticks(ynv); ax2.set_xticklabels(ynticks)
    ax2.set_xlim(0, sh[0])
    ax2.set_xlabel(r'$N$')
    ax.set_xlabel(r'$\log\ M/M_\odot$')
    
    xlam = [0.25,0.5,0.75,1,1.25]
    xlam = np.arange(4000,8000,1000)/1.e4
    xpix = np.interp(xlam, lrest/1.e4, np.arange(len(lrest)))
    ax.set_yticks(xpix); ax.set_yticklabels(np.cast[int](xlam*1e4))
    ax2.set_yticks(xpix); ax2.set_yticklabels(np.cast[int](xlam*1e4))
    ax.set_ylabel(r'$\lambda_\mathrm{rest}$ / \AA')
    ax.set_ylim(np.interp([0.33, 0.72], lrest/1.e4, np.arange(len(lrest))))
    fig.tight_layout()
    
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    ax.yaxis.set_minor_locator(MultipleLocator(np.diff(xpix)[0]/2.))
    ax2.yaxis.set_minor_locator(MultipleLocator(np.diff(xpix)[0]/2.))
    
    unicorn.plotting.savefig(fig, '3dhst_allspecs_mass.pdf')
    
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

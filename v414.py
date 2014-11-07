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
        
    
    phot.write('3dhst.v4.1.4.full.v1.fits')
    
    full = catIO.Table('3dhst.v4.1.4.full.v1.fits')
    

def make_selection_webpage(full, selection, output='test.html', columns=['spec_id', 'ra', 'dec', 'hmag', 'z_max_grism']):
    """
    Make an HTML table with image links
    
    full = catIO.Table('3dhst.v4.1.4.full.v1.fits')
    
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
    
    BASE = 'http://unicorn.astro.yale.edu/P/GRISM_v4.1.4/HTML/PNG'
    
    zfit, zfit2, linefit = [], [], []
    for id in sub['spec_id']:
        field = id.split('-')[0]
        pointing = id.split('-G141')[0]
        pngdir = '%s/%s/%s/' %(BASE, field, pointing)
        zfit.append('<img src=%s/%s.new_zfit.png height=200>' %(pngdir, id))
        zfit2.append('<img src=%s/%s.new_zfit.2D.png height=200>' %(pngdir, id))
        linefit.append('<img src=%s/%s.linefit.png height=200>' %(pngdir, id))
        
    sub.add_column(Column(data=np.array(zfit), name='zfitpng'))
    sub.add_column(Column(data=np.array(zfit2), name='zfit2png'))
    sub.add_column(Column(data=np.array(linefit), name='linefit'))
    
    
    show_cols = copy.copy(columns)
    show_cols.extend(['zfitpng', 'zfit2png', 'linefit'])
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
    
    
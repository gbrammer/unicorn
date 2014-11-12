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
    phot.add_column(Column(data=cat2.template_O2, name='eazy_OIII'))
    
    field_id = np.ones(len(phot), dtype=int)
    for fi, field in enumerate(['AEGIS', 'COSMOS', 'GOODS-N', 'GOODS-S', 'UDS']):
        field_id[phot['field'] == field] = fi+1
    
    phot.add_column(Column(data=field_id, name='field_id'), index=2)
    
    phot.write('3dhst.v4.1.4.full.v2.fits')
    
    #full = catIO.Table('3dhst.v4.1.4.full.v1.fits')
    
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
    
def random_objects(N=200):
    full = catIO.Table('../../3dhst.v4.1.4.full.v1.fits')
    has_spec = full['z_max_grism'] > 0.01
    objects = full['spec_id'][has_spec]
    idx = np.argsort(np.random.normal(size=has_spec.sum()))
    
    for i in range(N):
        obj = objects[idx][i]
        os.system('ln -s /Volumes/WD5_ESOMAC/3DHST/Spectra/Release/v4.1.4/*4/*/ZFIT/PNG/%s*png .' %(obj))
        
        
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
    
    zfit, zfit2, linefit, thumb = [], [], [], []
    for id in sub['spec_id']:
        field = id.split('-')[0]
        pointing = id.split('-G141')[0]
        objid = id.split('G141_')[1]
        pngdir = '%s/%s/%s/' %(BASE, field, pointing)
        zfit.append('<img src=%s/%s.new_zfit.png height=200>' %(pngdir, id))
        zfit2.append('<img src=%s/%s.new_zfit.2D.png height=200>' %(pngdir, id))
        linefit.append('<img src=%s/%s.linefit.png height=200>' %(pngdir, id))
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
            
            
    
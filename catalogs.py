import os
import pyfits
import numpy as np
import glob
import shutil
import re
import time

import matplotlib.pyplot as plt

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from pyraf import iraf
from iraf import iraf

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO
import unicorn

import cosmocalc

noNewLine = '\x1b[1A\x1b[1M'

zout = None
phot = None
mcat = None
lines = None
rest = None
gfit = None
selection_params = None

## from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit

def read_catalogs(force=False):
    """ 
    Read all of the catalogs and run the matching between them.
    """
    import unicorn.catalogs
    from unicorn.catalogs import match_string_arrays
    
    if (unicorn.catalogs.rest is not None) & (not force):
        print 'Looks like catalogs already read in.  To redo, use `read_catalogs(force=True)`.'
        return True
        
    PATH_TO_CAT = unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/'
        
    ######################################
    #### Redshifts, master ID list
    ######################################
    
    print noNewLine+'Reading Redshifts...'
    
    zout = catIO.Readfile(PATH_TO_CAT+'full_redshift.cat')
    
    pointing = []
    field = []
    for obj in zout.id[0::3]:
        point = obj.split('-G141')[0]
        pointing.append(point)
        field.append(re.split('-[1-9]',point)[0])
    
    pointing = np.array(pointing)
    field = np.array(field)
    
    zout.pointing = pointing
    zout.field = field
    
    unicorn.catalogs.zout = zout
    
    ######################################
    #### Full photometry
    ######################################
    
    print noNewLine+'Reading SExtractor photometry...'
    
    phot = catIO.Readfile(PATH_TO_CAT+'full_sextractor.cat')

    if not os.path.exists(PATH_TO_CAT+'phot.npz'):
        found, idx = match_string_arrays(zout.id[0::3], phot.id)
        np.savez(PATH_TO_CAT+'phot.npz', found=found, idx=idx)
    else:
        npz = np.load(PATH_TO_CAT+'phot.npz')
        found, idx = npz['found'], npz['idx']

    #found, idx = match_string_arrays(zout.id[0::3], phot.id)
    phot.idx = idx
    
    unicorn.catalogs.phot = phot
    
    ######################################
    #### Matches to photometric catalogs
    ######################################
    
    print noNewLine+'Reading external catalog matches...'
    
    mcat = catIO.Readfile(PATH_TO_CAT+'full_match.cat')
    if not os.path.exists(PATH_TO_CAT+'mcat.npz'):
        found, idx = match_string_arrays(zout.id[0::3], mcat.id_f140w)
        np.savez(PATH_TO_CAT+'mcat.npz', found=found, idx=idx)
    else:
        npz = np.load(PATH_TO_CAT+'mcat.npz')
        found, idx = npz['found'], npz['idx']
        
    mcat.idx = idx
    unicorn.catalogs.mcat = mcat
    
    ######################################
    #### Emission line fits
    ######################################
    
    print noNewLine+'Reading Emission lines...'
    
    lines = catIO.Readfile(PATH_TO_CAT+'full_emission_lines.cat')
    if not os.path.exists(PATH_TO_CAT+'lines.npz'):
        found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
        np.savez(PATH_TO_CAT+'lines.npz', found=found_lines, idx=idx_lines)
    else:
        npz = np.load(PATH_TO_CAT+'lines.npz')
        found_lines, idx_lines = npz['found'], npz['idx']
    
    lines.idx = idx_lines
    unicorn.catalogs.lines = lines
    
    ######################################
    #### Rest-frame colors
    ######################################
    
    print noNewLine+'Reading Rest-frame colors...'
    
    rest = catIO.Readfile(PATH_TO_CAT+'full_rf_fluxes.cat.partial')
    if not os.path.exists(PATH_TO_CAT+'rest.npz'):
        found_rest, idx_rest = match_string_arrays(zout.id[0::3], rest.id)
        np.savez(PATH_TO_CAT+'rest.npz', found=found_rest, idx=idx_rest)
    else:
        npz = np.load(PATH_TO_CAT+'rest.npz')
        found_rest, idx_rest = npz['found'], npz['idx']
    
    rest.idx = idx_rest
    unicorn.catalogs.rest = rest
    
    ##################################### 
    #### Galfit
    #####################################
    
    print noNewLine+'Reading GALFIT...'
    
    gfit = catIO.Readfile(PATH_TO_CAT+'full_galfit.cat')
    if not os.path.exists(PATH_TO_CAT+'gfit.npz'):
        found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.id, zout.id[0::3])
        found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.id)
        np.savez(PATH_TO_CAT+'gfit.npz', found=found_gfit, idx=idx_gfit, found_rev=found_gfit_rev, idx_rev=idx_gfit_rev)
    else:
        npz = np.load(PATH_TO_CAT+'gfit.npz')
        found_gfit, idx_gfit, found_gfit_rev, idx_gfit_rev = npz['found'], npz['idx'], npz['found_rev'], npz['idx_rev']

    # found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.id, zout.id[0::3])
    # found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.id)
    
    #### Make a grid of redshifts to compute the plate scale and then interpolate it.
    zgrid = np.arange(100)/100.*4+1./100
    scale = zgrid*0.
    for i in range(100):
        cc = cosmocalc.cosmocalc(zgrid[i])
        scale[i] = cc['PS_kpc']
    
    gfit.r_e_kpc = gfit.r_e*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_err = gfit.r_e_err*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(np.abs(gfit.ba)) 
    gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(np.abs(gfit.ba)) 
    
    gfit.idx = idx_gfit
    unicorn.catalogs.gfit = gfit

class selectionParams():
    def __init__(self, zmin=1.0, zmax=1.5, fcontam=0.2, qzmin=0., qzmax=0.4, dr=1.0, has_zspec=False, fcovermin=0.9, fcovermax=1.0, massmin=7, massmax=15, magmin=0, magmax=30):
        self.zmin=1.0
        self.zmax=1.5
        self.fcontam=0.2
        self.qzmin=0.
        self.qzmax=0.4
        self.dr=1.0
        self.has_zspec=False
        self.fcovermin=0.9
        self.fcovermax=1.0
        self.massmin=11
        self.massmax=15
        self.magmin=0
        self.magmax=30
        
def run_selection(zmin=1.0, zmax=1.5, fcontam=0.2, qzmin=0., qzmax=0.4, dr=1.0, has_zspec=False, fcovermin=0.9, fcovermax=1.0, massmin=7, massmax=15, magmin=0, magmax=30):
    """
    Run a selection on the 3D-HST catalogs
    """
    import unicorn.catalogs
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    dr_match = mcat.rmatch[mcat.idx] < dr
    
    zrange = (zout.z_peak[0::3] >= zmin) & (zout.z_peak[0::3] <= zmax)
    
    keep = dr_match & zrange & (mcat.fcontam[mcat.idx] <= fcontam) & (zout.q_z[0::3] >= qzmin)  & (zout.q_z[0::3] <= qzmax) 
    
    keep = keep & (phot.fcover[phot.idx] >= fcovermin) & (phot.fcover[phot.idx] <= fcovermax)
    
    keep = keep & (mcat.logm[mcat.idx] >= massmin) & (mcat.logm[mcat.idx] <= massmax)

    keep = keep & (phot.mag_f1392w[phot.idx] >= magmin) & (phot.mag_f1392w[phot.idx] <= magmax)
    
    unicorn.catalogs.selection_params = selectionParams(zmin=zmin, zmax=zmax, fcontam=fcontam, qzmin=qzmin, qzmax=qzmax, dr=dr, has_zspec=has_zspec, fcovermin=fcovermin, fcovermax=fcovermax, massmin=massmin, massmax=massmax)
    
    return keep

def make_selection_catalog(selection, filename='selection.cat', make_html=True):
    """ 
    Make a single catalog for a given `selection` that combines the outputs
    of the separate 3D-HST catalogs
    """
    import unicorn
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    UV = -2.5*np.log10(rest.l153/rest.l155)
    VJ = -2.5*np.log10(rest.l155/rest.l161)
    
    int_idx = np.arange(len(phot.idx))[selection]
    #print gfit.id[gfit.idx][selection][0:10]
    #print gfit.id[gfit.idx][int_idx][0:10]
    
    fp = open(filename,'w')
    fp.write('#   object   z_fast   z_gris   z_phot2   z_spec   Q_z   mag_f140w   fcontam   lmass    Av   Umv   VmJ      eqwHa   eqwHa_err   r_e_pix  r_e_pix_err   r_e_kpc   r_e_circ   sersic_n  sersic_n_err  ba \n# eqw: observed frame\n')
    
    for ii in int_idx:
        object = gfit.id[gfit.idx][ii]
        print noNewLine+object
        #### Match to FAST catalog to get Av
        cat, zout_phot, fout = unicorn.analysis.read_catalogs(root=object)
        fmat = np.where(fout.id == mcat.id_phot[mcat.idx][ii])[0][0]
        #fout.lmass[fmat], mcat.logm[mcat.idx][ii]
        
        string_line = ' %-30s %8.4f %8.4f %8.4f %8.4f  %.1e   %6.3f  %4.2f  %7.3f  %7.3f  %8.3f %8.3f   %15.2e  %15.2e  %8.2f %8.2f  %8.2f  %8.2f  %8.1f  %8.1f  %8.2f' %(object, fout.z[fmat], zout.z_peak[0::3][ii], zout.z_peak[1::3][ii], zout.z_spec[0::3][ii], zout.q_z[0::3][ii], phot.mag_f1392w[phot.idx][ii], phot.fcontam[phot.idx][ii],  fout.lmass[fmat], fout.Av[fmat], UV[rest.idx][ii], VJ[rest.idx][ii], lines.halpha_eqw[lines.idx][ii], lines.halpha_eqw_err[lines.idx][ii], gfit.r_e[gfit.idx][ii], gfit.r_e_err[gfit.idx][ii], gfit.r_e_kpc[gfit.idx][ii], gfit.r_e_kpc_circ[gfit.idx][ii], gfit.n[gfit.idx][ii], gfit.n_err[gfit.idx][ii], gfit.ba[gfit.idx][ii])
        
        fp.write(string_line+'\n')
    
    fp.close()
    
    if make_html:
        make_selection_html(catalog_file=filename)
        
def make_selection_html(catalog_file='selection.cat'):
    """
    Make a webpage for a given selection.
    """
    import unicorn.catalogs
    from unicorn.catalogs import selection_params as par
    
    if par is None:
        par = unicorn.catalogs.selectionParams()
        
    cat = catIO.Readfile(catalog_file)
    
    head="""
    <html> 
    <head> 
    <link rel="stylesheet" href="../scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="../scripts/jquery-1.4.2.min.js"></script> 
    
    <script type="text/javascript" src="../scripts/jquery.tablesorter.min.js"></script> 
    
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
                        8: {
                                sorter: false
                        },
                        9: {
                                sorter: false
                        },
                }
        });        
    });
    </script> 
    
    </head> 
    <body> 
    <p> 
        %.2f < <b>mag_F140W</b> < %.2f
        <br> %.2f < <b>log M/Msun</b> < %.2f
        <br> %.2f < <b>z_gris</b> < %.2f
        <br> <b>fcontam</b> < %.2f
        <br> %.2f < <b>fcover</b> < %.2f
    </p> 
    <h2> N=%d, <a href=./%s>catalog</a> </h2>
    
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
        <th> EAZY </th> 
        <th> GALFIT </th> 
    </thead> 
    <tbody> 
    """ %(par.magmin, par.magmax, par.massmin, par.massmax, par.zmin, par.zmax, par.fcontam, par.fcovermin, par.fcovermax, cat.N, catalog_file)
    
    lines = [head]
    for i in range(cat.N):
        line="""
        <tr> 
            <td> %s </td> 
            <td> %.2f </td> 
            <td> %.3f </td> 
            <td> %.2f </td> 
            <td> <img src=../images/%s_thumb.png height=180px> </td> 
            <td> <img src=../images/%s_2D.png height=180px> </td> 
            <td> <img src=../images/%s_1D.png height=180px> </td> 
            <td> <img src=../SED/%s_SED.png height=180px> </td> 
            <td> <img src=../EAZY/%s_eazy.png height=180px> </td> 
            <td> <img src=../GALFIT/%s_galfit.png height=180px> </td> 
        </tr> """ %(cat.object[i], cat.mag_f140w[i], cat.z_gris[i], cat.lmass[i], cat.object[i], cat.object[i], cat.object[i], cat.object[i], cat.object[i], cat.object[i])
        
        lines.append(line)
    
    lines.append("        </tbody></table></body></html>")
    
    fp = open(catalog_file.replace('.cat','.html'),'w')
    fp.writelines(lines)
    fp.close()
    
    base = catalog_file.split('.cat')[0]
    
    print '! rsync -avz %s.cat %s.html ~/Sites_GLOBAL/P/GRISM_v1.6/ANALYSIS/' %(base, base)

def make_selection_for_Pieters_paper():
    import unicorn.catalogs
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    
    unicorn.catalogs.read_catalogs()
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    keep = unicorn.catalogs.run_selection(zmin=1.0, zmax=1.5, fcontam=0.2, qzmin=0., qzmax=0.1, dr=1.0, has_zspec=False, fcovermin=0.9, fcovermax=1.0, massmin=11, massmax=15, magmin=0, magmax=30)
    
    print len(keep[keep])
    
    unicorn.catalogs.make_object_tarfiles(zout.id[0::3][keep], thumbs=False)
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    os.system('mv ../spec.tar.gz for_pieter_Aug10.spec.tar.gz')
    
    unicorn.catalogs.make_selection_catalog(keep, filename='for_pieter_Aug10.cat', make_html=True)

    os.system('rsync -avz for_pieter_Aug10.* ~/Sites_GLOBAL/P/GRISM_v1.6/ANALYSIS/')
    
def match_string_arrays(target=['b','a','d'], source=['a','b','c']):
    """
    
    Given a list of strings in `target`, return the integer indices of
    `source` where the strings match.  Fill with -1 if no match so that 
    the output array has the same size as `target`.
    
    Example:
    >>> target=np.array(['b','a','d'])
    >>> source=np.array(['a','b','c'])
    >>> found_in_source, indices = matche_string_arrays(target, source)
    >>> print target[found_in_source], source[indices]
    
    """
    
    indices = np.zeros(len(target), dtype='int')-1
    source_arr = np.array(source)
    for i,targ in enumerate(target):
        match = np.where(source_arr == targ)[0]
        if len(match) > 0:
            indices[i] = int(match[0])
    
    found_in_source = indices >= 0
    
    return found_in_source, indices
    
def test_plots():
    import unicorn.catalogs
    from unicorn.catalogs import match_string_arrays
    import cosmocalc
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    
    ######################################
    #### Full photometry
    ######################################
    
    phot = catIO.Readfile('full_sextractor.cat')
    #phot_goods = catIO.Readfile('goodsn_sextractor.cat')
    
    ##### Marijn's test of size vs mag
    ours = (phot.field == 'AEGIS') | (phot.field == 'COSMOS') | (phot.field == 'GOODS-S')
    theirs = (phot.field == 'GOODS-N') 
    candels = (phot.field == 'PRIMO') | (phot.field == 'GEORGE') | (phot.field == 'MARSHALL')
    
    unicorn.catalogs.plot_init()
    
    plt.plot(phot.mag_f1392w[ours], phot.flux_radius[ours], marker='.', alpha=0.1, linestyle='None', color='red', markersize=3)
    plt.plot(phot.mag_f1392w[theirs], phot.flux_radius[theirs], marker='.', alpha=0.1, linestyle='None', color='blue', markersize=3)
    plt.plot(phot.mag_f1392w[candels], phot.flux_radius[candels], marker='.', alpha=0.5, linestyle='None', color='green', markersize=3)

    plt.xlim(14,25.7)
    plt.ylim(0,20)
    plt.xlabel(r'$\mathrm{MAG\_AUTO_{F140W}}$')
    plt.ylabel(r'$\mathrm{FLUX\_RADIUS}\ (0.06^{\prime\prime}\ \mathrm{pix})$')
    
    plt.text(15,17,'3D-HST',color='red')
    plt.text(15,15,'GOODS-N',color='blue')

    plt.savefig('mag_vs_SEx_size.pdf')
    plt.savefig('mag_vs_SEx_size.png')
    
    #### S/N vs mag
    sn = 2.5/np.log(10)/phot.magerr_aper
    #sn = phot.flux_auto/phot.fluxerr_auto
    
    ok = np.isfinite(sn)
    unicorn.catalogs.plot_init()
    plt.plot(phot.mag_f1392w[ours & ok], sn[ours & ok], marker='.', alpha=0.1, linestyle='None', color='red', markersize=3)
    plt.plot(phot.mag_f1392w[theirs & ok], sn[theirs & ok], marker='.', alpha=0.1, linestyle='None', color='blue', markersize=3)
    
    plt.semilogy()
    plt.xlim(18,25.7)
    plt.ylim(3,1000)
    plt.xlabel(r'$\mathrm{MAG\_AUTO_{F140W}}$')
    plt.ylabel(r'$\mathrm{S/N}\ (1^{\prime\prime}\ \mathrm{aper})$')
    
    plt.savefig('mag_vs_SN.pdf')
    plt.savefig('mag_vs_SN.png')
    
    ######################################
    #### Redshifts and masses
    ######################################
    
    zout = catIO.Readfile('full_redshift.cat')
    unicorn.catalogs.zout = zout
    
    pointing = []
    field = []
    for obj in zout.id[0::3]:
        point = obj.split('-G141')[0]
        pointing.append(point)
        field.append(re.split('-[1-9]',point)[0])
    
    pointing = np.array(pointing)
    field = np.array(field)
    
    mcat = catIO.Readfile('full_match.cat')
        
    zsp = zout.z_spec[0::3] > 0
    dz = (zout.z_peak-zout.z_spec)/(1+zout.z_spec)
    found, idx = match_string_arrays(zout.id[0::3], mcat.id_f140w)
    #idx = idx[idx >= 0]
    zsp = zsp & found

    ######################################
    #### Emission line fits
    ######################################
    
    lines = catIO.Readfile('full_emission_lines.cat')
    found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
    #lines.halpha_eqw[idx_lines] /= (1+zout.z_peak[0::3][found_lines])
    red_limit = 15
    red_sequence = lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3]) < red_limit

    ######################################
    #### Rest-frame colors
    ######################################
    
    rest = catIO.Readfile('full_rf_fluxes.cat')
    found_rest, idx_rest = match_string_arrays(zout.id[0::3], rest.id)
    
    ### plt.plot(zout.z_peak[0::3], rest.z_grism[idx_rest], marker='o', linestyle='None', color='blue', alpha=0.5); plt.xlim(0,5); plt.ylim(0,5)
    
    ##################################### 
    #### Galfit
    #####################################
    gfit = catIO.Readfile('full_galfit.cat')
    found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.id, zout.id[0::3])
    found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.id)
    
    #### Make a grid of redshifts to compute the plate scale and then interpolate it.
    zgrid = np.arange(100)/100.*4+1./100
    scale = zgrid*0.
    for i in range(100):
        cc = cosmocalc.cosmocalc(zgrid[i])
        scale[i] = cc['PS_kpc']
    
    gfit.r_e_kpc = gfit.r_e*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_err = gfit.r_e_err*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(np.abs(gfit.ba)) 
    gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(np.abs(gfit.ba)) 
    
    #############
    ####   Compare galfit mags to SExtractor
    ##############
    phot = catIO.Readfile('full_sextractor.cat')
    gfit = catIO.Readfile('full_galfit.cat')
    found_gfit, idx_gfit = match_string_arrays(gfit.id, phot.id)
    dmag = gfit.mag-phot.mag_f1392w[idx_gfit]
    
    f140 = (phot.field == 'COSMOS') | (phot.field == 'AEGIS') | (phot.field == 'GOODS-S') | (phot.field == 'GOODS-N')
    
    ok = (gfit.n > 0) & (gfit.mag_err < 0.5)
    plt.plot(gfit.mag, gfit.mag_err, marker='o', alpha=0.1, linestyle='None')

    plt.plot(gfit.mag[ok & f140[idx_gfit]], dmag[ok & f140[idx_gfit]], marker='o', alpha=0.05, linestyle='None')
    plt.xlim(14,25)
    plt.ylim(-2,2)
    
    ##### Selection slices
    dr = mcat.rmatch[idx] < 1
    zrange = (zout.z_peak[0::3] > 0.2)
    
    keep = dr & zrange & (mcat.logm[idx] > 10.98) & (mcat.fcontam[idx] < 0.5)
    keep = dr & zrange & (mcat.mag_f140w[idx] < 24) & (mcat.fcontam[idx] < 0.5)
    
    #### refine
    # zrange = (zout.z_peak[0::3] > 2) & (zout.z_peak[0::3] < 3)
    # zrange = (zout.z_peak[0::3] > 1.5) & (zout.z_peak[0::3] < 2.)
    # zrange = (zout.z_peak[0::3] > 0.4) & (zout.z_peak[0::3] < 1.0)
    zrange = (zout.z_peak[0::3] > 1.0) & (zout.z_peak[0::3] < 1.5)
    keep = dr & zrange & (mcat.fcontam[idx] < 0.2) & (zout.q_z[0::3] < 0.4) 
    #keep = dr & zrange & (mcat.fcontam[idx] < 0.02) & (zout.q_z[0::3] < 0.1) 
    #keep = dr & zrange & (mcat.fcontam[idx] < 0.2) & (zout.q_z[0::3] < 1) 
    
    ### copy to macbook
    copy = keep & (mcat.logm[idx] > 8.)
    # unicorn.catalogs.make_object_tarfiles(zout.id[0::3][copy], thumbs=False)
    
    ##### Q_z vs dz
    plt.semilogx(zout.q_z[0::3][zsp & dr], dz[0::3][zsp & dr], marker='o', linestyle='None', color='red', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(1.e-3, 10)
    plt.xlabel(r'$Q_z$')
    plt.ylabel(r'$\Delta z$')

    ##### z vs dz
    qz = zout.q_z < 0.1
    eqw = lines.halpha_eqw[idx_lines] < 30
    
    main = zsp & dr
    test = eqw[main]
    
    plt.plot(zout.z_spec[0::3][main][test], dz[0::3][main][test], marker='o', linestyle='None', color='red', alpha=0.05)
    plt.plot(zout.z_spec[0::3][main][~test], dz[0::3][main][~test], marker='o', linestyle='None', color='blue', alpha=0.05)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0, 4)
    plt.xlabel(r'$z_\mathrm{spec}$')
    plt.ylabel(r'$\Delta z$')
    
    ################################################
    #### Real figure, redshift errors
    ################################################
    fig = unicorn.catalogs.plot_init()
    
    #fig = plt.figure(figsize=(5.5,5), dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.15,
                        bottom=0.10,right=0.99,top=0.97)

    stats = abs(dz[0::3]) < 0.06 
    ## zphot
    ax = fig.add_subplot(211)
    
    dzmcat = (mcat.z_peak[idx]-zout.z_spec[0::3])/(1+zout.z_spec[0::3])
    plt.plot(zout.z_spec[0::3][zsp & keep], dzmcat[zsp & keep], marker='o', linestyle='None', color='orange', alpha=0.4)
    hist = np.histogram(dzmcat[zsp & keep], bins=100, range=(-0.05, 0.05))
    plt.plot(hist[0]*1./np.sum(hist[0])+0.85, hist[1][1:], linestyle='steps', alpha=0.8, color='orange')

    plt.ylim(-0.059, 0.059)
    plt.xlim(0.85, 1.7)
    plt.xlabel(r'$z_\mathrm{spec}$')
    plt.ylabel(r'$\Delta z$')
    plt.text(1.5,-0.05,r'$\sigma=%.3f$' %(threedhst.utils.biweight(dzmcat[zsp & keep & stats])), fontsize=16)
    
    ## gspec
    ax = fig.add_subplot(212)
    plt.plot(zout.z_spec[0::3][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='purple', alpha=0.4)
    hist = np.histogram(dz[0::3][zsp & keep], bins=100, range=(-0.05, 0.05))
    plt.plot(hist[0]*1./np.sum(hist[0])+0.85, hist[1][1:], linestyle='steps', alpha=0.8, color='purple')

    plt.ylim(-0.059, 0.059)
    plt.xlim(0.85, 1.7)
    plt.xlabel(r'$z_\mathrm{spec}$')
    plt.ylabel(r'$\Delta z$')
    plt.text(1.5,-0.05,r'$\sigma=%.3f$' %(threedhst.utils.biweight(dz[0::3][zsp & keep & stats])), fontsize=16)
    plt.text(0.9,-0.05,r'$N=%0d$' %(len(dz[0::3][zsp & keep & stats])), fontsize=16)

    plt.savefig('redshift_errors.pdf')
    plt.savefig('redshift_errors.png')
    
    ##### F140W_mag vs dz
    plt.plot(mcat.mag_f140w[idx][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(15, 24)
    plt.xlabel(r'$m_\mathrm{f140w}$')
    plt.ylabel(r'$\Delta z$')
        
    ##### Contamination vs dz
    plt.plot(mcat.fcontam[idx][zsp & keep]+1.e-2, dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0, 0.5)
    plt.xlabel(r'$f_\mathrm{contam}$')
    plt.ylabel(r'$\Delta z$')
    
    ##### mass vs dz
    plt.plot(mcat.logm[idx][zsp & keep], dz[0::3][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(8, 12)
    plt.xlabel(r'$\log M/M_\odot$')
    plt.ylabel(r'$\Delta z$')

    ################################################
    ##### Redshift distribution
    ################################################
    unicorn.catalogs.plot_init()
    
    hist = np.histogram(zout.z_spec[0::3][keep], bins=200, range=(1,1.5))
    plt.plot((hist[1][:-1]+hist[1][1:])/2., 0-hist[0], linestyle='steps', color='red', alpha=0.8)
    hist = np.histogram(zout.z_peak[0::3][keep], bins=200, range=(1,1.5))
    plt.plot((hist[1][:-1]+hist[1][1:])/2., hist[0], linestyle='steps', color='black', alpha=0.8)
    hist_wide = np.histogram(zout.z_peak[0::3][keep], bins=6, range=(0.9,1.5))
    plt.plot(hist_wide[1][1:], hist_wide[0]/20., linestyle='steps', color='black', alpha=0.5, linewidth=4)
    #plt.plot((hist_wide[1][:-1]+hist_wide[1][1:])/2., hist_wide[0]/20., linestyle='steps', color='black', alpha=0.5, linewidth=4)
    plt.xlim(1,1.5)
    plt.ylim(-10,30)
    plt.xlabel(r'$z_\mathrm{grism}\ (\Delta z=0.0025)$')
    plt.ylabel(r'$N\ \ (N_\mathrm{tot}=%d)$' %len(zout.z_peak[0::3][keep]))
   
    peak = np.abs(zout.z_peak[0::3][keep]-1.020) < 0.005
    zout.id[0::3][keep][peak]
    
    plt.savefig('redshift_dist.pdf')
    plt.savefig('redshift_dist.png')
    
    ##### line strength vs dz
    
    haw = lines.oiii_eqw
    xmax = 2500
    NBAD = len(haw[haw > xmax])
    haw[haw > xmax] = xmax + np.random.randn(NBAD)*10.
    
    plt.semilogx(haw[idx_lines][zsp & dr], dz[0::3][zsp & dr], marker='o', linestyle='None', color='blue', alpha=0.4)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0.1, xmax+100)
    plt.xlabel(r'$\mathrm{H}\alpha\ \mathrm{EQW}$')
    plt.ylabel(r'$\Delta z$')
    
    ################################################
    ##### plot mass vs equivalent width
    ################################################
    unicorn.catalogs.plot_init()
    
    plt.semilogy(mcat.logm[idx][~zsp & keep], lines.halpha_eqw[idx_lines][~zsp & keep]/(1+zout.z_peak[0::3][~zsp & keep]), marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & keep], lines.halpha_eqw[idx_lines][zsp & keep]/(1+zout.z_peak[0::3][zsp & keep]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    # oiii = keep & (zout.z_peak[0::3] > 1.2) & (lines.oiii_eqw[idx_lines] > 25)
    # plt.semilogy(mcat.logm[idx][oiii], lines.halpha_eqw[idx_lines][oiii], marker='s', linestyle='None', color='red', alpha=0.2, markersize=8)
    
    ## red sequence
    plt.plot([0,100],[red_limit,red_limit], linestyle='--', alpha=0.5, color='green', linewidth=4)

    plt.text(9.5,1,r'has $z_\mathrm{spec}$', color='blue', fontsize=18)
    
    ### For SDSS, read in catalogs and get selection from unicorn.paper1
    # sel, NOBJ = p1.sdss_selection(zmax=0.2)
    # hist, xedge, yedge = np.histogram2d(p1.sdss_logm.AVG[sel], np.log10(-p1.sdss_line.H_ALPHA_EQW[sel]/(1+p1.sdss_info.Z[sel])), bins=100, range=[[9,11.6], [-1, 2.7]])
    # Vbins = [2, 4, 8, 16, 32, 64, 128, 256, 512, 4096]
    # values =   1.-np.arange(len(Vbins))*1./len(Vbins)
    # Vcolors = []
    # for i in range(len(Vbins)):
    #     Vcolors.append('%f' %(values[i]))
    # 
    # unicorn.catalogs.plot_init()
    # plt.contourf(xedge[1:], 10**yedge[1:], hist.transpose(), Vbins, colors=Vcolors, alpha=1.0, linethick=2)
    # plt.semilogy(mcat.logm[idx][keep], lines.halpha_eqw[idx_lines][keep]/(1+zout.z_peak[0::3][keep]), marker='o', linestyle='None', color='red', alpha=0.5, markersize=5)
    
    plt.ylim(0.1, 500)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{H}\alpha\ \mathrm{eqw}$')

    # plt.savefig('ha_mass_sdss.pdf')
    # plt.savefig('ha_mass_sdss.png')
    
    plt.savefig('ha_mass_zspec.pdf')
    plt.savefig('ha_mass_zspec.png')
    
    ################################################
    ####### H-alpha equivalent width signal to noise
    ################################################
    fig = unicorn.catalogs.plot_init()
    ax = fig.add_subplot(111)
    
    sn_halpha = (lines.halpha_eqw/lines.halpha_eqw_err)
    ### scale by object size
    #sn_halpha[idx_lines] *= gfit.r_e[idx_gfit]/10.
    ax.plot([0.01,3000],[3,3], color='orange',alpha=0.2, linewidth=3)
    #sn_halpha = lines.halpha_eqw_err
    #plt.plot([0.03,3000],[0.01,1000], color='orange',alpha=0.2, linewidth=3)
    
    keep = keep & (lines.halpha_eqw_err[idx_lines] > 0)

    ax.semilogx(lines.halpha_eqw[idx_lines][~zsp & keep]/(1+0*zout.z_peak[0::3][~zsp & keep]), sn_halpha[idx_lines][~zsp & keep], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    ax.semilogx(lines.halpha_eqw[idx_lines][zsp & keep]/(1+0*zout.z_peak[0::3][zsp & keep]), sn_halpha[idx_lines][zsp & keep], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=6)
    
    xm, ym, ys, ns = threedhst.utils.runmed(lines.halpha_eqw[idx_lines][keep]/(1+0*zout.z_peak[0::3][keep]), sn_halpha[idx_lines][keep], NBIN=20)
    ax.plot(xm, ym, color='orange', alpha=0.5, linewidth=8)
    
    ax.semilogy()
    
    ax.set_xticklabels(['1','10','100'])
    ax.set_xticks([1,10,100])

    ax.set_yticklabels(['0.3','1','3','10','30','100'])
    ax.set_yticks([0.3,1,3,10,30,100])
    
    ax.plot([10,10], [0.01,3000], color='orange',alpha=0.2, linewidth=3)
    ax.set_xlim(0.3, 800)
    ax.set_ylim(0.1, 200)
    ax.set_ylabel(r'$\mathrm{eqw\ signal-to-noise}$')
    ax.set_xlabel(r'$\mathrm{H}\alpha\ \mathrm{eqw\ (observed\ frame)}$')
    
    fig.savefig('ha_eqw_errors.pdf')
    fig.savefig('ha_eqw_errors.png')
    
    ####### OIII emitters
    oiii = keep & (zout.z_peak[0::3] > 1.2)
    plt.semilogy(mcat.logm[idx][oiii], lines.halpha_eqw[idx_lines][oiii], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & oiii], lines.halpha_eqw[idx_lines][zsp & oiii], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    plt.ylim(0.01, 500)
    plt.xlim(8.5, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{H}\alpha\ \mathrm{EQW}$')

    plt.semilogy(mcat.logm[idx][oiii], lines.oiii_eqw[idx_lines][oiii], marker='o', linestyle='None', color='black', alpha=0.2, markersize=4)
    plt.semilogy(mcat.logm[idx][zsp & oiii], lines.oiii_eqw[idx_lines][zsp & oiii], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    
    plt.ylim(0.01, 500)
    plt.xlim(8.5, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\mathrm{OIII}\ \mathrm{eqw}$')
    
    ################################################
    #### Mass - size
    ################################################
    unicorn.catalogs.plot_init(square=True)
    
    plt.semilogy(mcat.logm[idx][keep & ~red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.1, markersize=8)
    plt.semilogy(mcat.logm[idx][keep & red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.1, markersize=8)
    
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep], gfit.r_e_kpc_circ[idx_gfit][keep], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='black', ecolor='black', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)
    
    vd08 = catIO.Readfile('vd08.dat')
    plt.plot(vd08.logm, 10**vd08.log_re, marker='o', color='orange', alpha=0.7, markersize=12, linestyle='None')

    vd08 = catIO.Readfile('vd08_sdss.dat')
    plt.plot(vd08.logm, 10**vd08.log_re, marker='None', color='green', alpha=0.3, linewidth=14)
    
    plt.ylim(0.2, 20)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$r_\mathrm{e}\ \mathrm{(kpc)}$')
    plt.savefig('mass_size.pdf')
    plt.savefig('mass_size.png')
    
    ################################################
    #### Mass - n
    ################################################
    fig = unicorn.catalogs.plot_init()
    ax = fig.add_subplot(111)
    
    ax.semilogy(mcat.logm[idx][keep & ~red_sequence], gfit.n[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.1, markersize=8)
    ax.semilogy(mcat.logm[idx][keep & red_sequence], gfit.n[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.1, markersize=8)

    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep], gfit.n[idx_gfit][keep], NBIN=10)
    ax.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='black', ecolor='black', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.n[idx_gfit][keep & ~red_sequence], NBIN=10)
    ax.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.n[idx_gfit][keep & red_sequence], NBIN=10)
    ax.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)

    ax.plot([0,100],[1,1], color='blue', linestyle='--', alpha=0.4, linewidth=4)
    ax.plot([0,100],[4,4], color='red', linestyle='--', alpha=0.4, linewidth=4)
    ax.text(11.5,1.1,r'$n=1$', fontsize=14, horizontalalignment='right')
    ax.text(11.5,4.4,r'$n=4$', fontsize=14, horizontalalignment='right')
    
    ax.set_yticklabels(['1','4','10'])
    ax.set_yticks([1,4,10])
    
    ax.set_ylim(0.2, 21)
    ax.set_xlim(9, 11.6)
    ax.set_xlabel(r'$\log\ M/M_\odot$')
    ax.set_ylabel(r'$n$')
    
    fig.savefig('mass_n.pdf')
    fig.savefig('mass_n.png')
    
    ################################################
    #### Mass - b/a
    ################################################
    unicorn.catalogs.plot_init()
    
    plt.plot(mcat.logm[idx][keep & ~red_sequence], gfit.ba[idx_gfit][keep & ~red_sequence], marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(mcat.logm[idx][keep & red_sequence], gfit.ba[idx_gfit][keep & red_sequence], marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)

    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & ~red_sequence], gfit.ba[idx_gfit][keep & ~red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='blue', ecolor='blue', alpha=0.6, marker='o', markersize=12)
    xm, ym, ys, ns = threedhst.utils.runmed(mcat.logm[idx][keep & red_sequence], gfit.ba[idx_gfit][keep & red_sequence], NBIN=10)
    plt.errorbar(xm, ym, yerr=ys/np.sqrt(ns), color='red', ecolor='red', alpha=0.6, marker='o', markersize=12)

    plt.ylim(0., 1)
    plt.xlim(9, 11.6)
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$b/a$')
    
    plt.savefig('mass_ba.pdf')
    plt.savefig('mass_ba.png')
        
    ################################################
    ######### Composite spectra
    ################################################
    eqw = lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3])
    
    mass = keep & (mcat.logm[idx] > 10.)
    
    eqws = eqw[mass]
    eqws.sort()
    quartiles = eqws[np.cast[int](np.round(np.array([0.25,0.5,0.75])*len(eqws)))]
    thirds = eqws[np.cast[int](np.round(np.array([1./3, 2./3])*len(eqws)))]
    
    high_eqw = (eqw > 40) & mass
    med_eqw = (eqw > 10) & (eqw <= 40) & mass
    low_eqw = (eqw <= 10) & mass
        
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    bins = quartiles
    colors = ['red','orange','green','blue']
    
    lnorm = 5900
    #### Lower limit
    eqw_bin = (eqw < bins[0]) & mass
    unicorn.catalogs.composite_spectra(lines.id[idx_lines][eqw_bin], color=colors[0], alpha=0.02, lnorm=lnorm, NITER=3, show_lines=True)
    plt.text(7700,1.6,r'$N_\mathrm{EW<%.1f}=%d$' %(bins[0], len(lines.id[idx_lines][eqw_bin])), color=colors[0], horizontalalignment='right')
    # plt.ylim(0.7,2)
    # plt.xlim(4300, 7900)
    
    #### Middle bins
    for i in range(0,len(bins)-1):
        eqw_bin = (eqw >= bins[i]) & (eqw < bins[i+1]) & mass
        unicorn.catalogs.composite_spectra(lines.id[idx_lines][eqw_bin], color=colors[i+1], alpha=0.02, lnorm=lnorm, NITER=3, show_lines=True)
        plt.text(7700,1.6+(i+1)*0.1,r'$N_\mathrm{%.1f<EW<%.1f}=%d$' %(bins[i], bins[i+1], len(lines.id[idx_lines][eqw_bin])), color=colors[i+1], horizontalalignment='right')
    
    #### Upper bin
    eqw_bin = (eqw > bins[-1]) & mass
    unicorn.catalogs.composite_spectra(lines.id[idx_lines][eqw_bin], color=colors[-1], alpha=0.02, lnorm=lnorm, NITER=3, show_lines=True)
    plt.text(7700,1.6+(i+2)*0.1,r'$N_\mathrm{EW>%.1f}=%d$' %(bins[-1], len(lines.id[idx_lines][eqw_bin])), color=colors[-1], horizontalalignment='right')
    
    plt.text(7700,1.5,r'$\log\ M/M_\odot > %.1f$' %(np.min(mcat.logm[idx][mass])), color='black', horizontalalignment='right')
    
    plt.ylim(0.7,2)
    plt.xlim(4500, 7900)
    
    plt.xlabel(r'$\lambda_\mathrm{rest}$')
    plt.ylabel(r'$f_\lambda$')
    
    plt.savefig('composite_spectra.pdf')
    plt.savefig('composite_spectra.png')
    
    ################################################
    #### Eq-W vs size
    ################################################
    eqw = lines.halpha_eqw[idx_lines]/(1+zout.z_peak[0::3])
    eqw[eqw < 0.2] = 0.2
    mass_limit = 10.6
    massive = mcat.logm[idx] > mass_limit
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    plt.plot(mcat.logm[idx][keep & ~red_sequence & massive], np.log10(eqw[keep & ~red_sequence & massive]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(mcat.logm[idx][keep & red_sequence & massive], np.log10(eqw[keep & red_sequence & massive]), marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)
    
    plt.plot([-10,-10],[-10,-10])
    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.xlim(9, 11.6)
    x = mcat.logm[idx][keep & massive]
    plt.ylim(-1, np.log10(200))
    y = np.log10(eqw[keep & massive])
    
    plt.savefig('eqw_mass9_morph_A.pdf')
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    
    plt.plot(np.log10(eqw[keep & ~red_sequence & massive]), np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & ~red_sequence & massive]), marker='o', linestyle='None', color='blue', alpha=0.2, markersize=8)
    plt.plot(np.log10(eqw[keep & red_sequence & massive]), np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & red_sequence & massive]), marker='o', linestyle='None', color='red', alpha=0.2, markersize=8)

    plt.xlim(-1, np.log10(200))
    plt.ylim(-0.2, 1.5)
    plt.xlabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.ylabel(r'$\log\ r_\mathrm{e}\ (\mathrm{kpc})$')
    y = np.log10(gfit.r_e_kpc_circ[idx_gfit][keep & massive])
    x = np.log10(eqw[keep & massive])
    
    plt.savefig('eqw_size_morph_A.pdf')
        
    objects = mcat.id_f140w[idx][keep & massive]
    xrange = plt.xlim()
    yrange = plt.ylim()
    scale = 1./10**(-0.4*(mcat.mag_f140w[idx][keep & massive]-22))
    
    #scale = scale*0+1
    
    im = unicorn.catalogs.fill_image(objects, x, y, scale=scale, xrange=xrange, yrange=yrange, NX=1024*1., NY=1024*1.)
    
    unicorn.catalogs.plot_init(square=True, xs=7)
    plt.imshow(0-im, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), aspect='auto', interpolation='nearest', vmin=-0.5, vmax=0.05)

    plt.xlabel(r'$\log\ \mathrm{EQW\ H}\alpha$')
    plt.ylabel(r'$\log\ r_\mathrm{e}\ (\mathrm{kpc})$')
    plt.savefig('eqw_size_morph_B.pdf')

    plt.xlabel(r'$\log\ M/M_\odot$')
    plt.ylabel(r'$\log\ \mathrm{EQW\ H}\alpha$')    
    plt.savefig('eqw_mass9_morph_B.pdf')

def uvj_test():
    """
    ################################################
    #### UVJ, ha_lines
    ################################################
    """
    import unicorn.catalogs
    
    unicorn.catalogs.read_catalogs()
    
    from unicorn.catalogs import zout, phot, mcat, lines, rest, gfit
    
    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/GRISM_v1.6/')
    
    #### Selection
    keep = unicorn.catalogs.run_selection(zmin=1.0, zmax=1.5, fcontam=0.25, qz=0.4, dr=1.0)
    
    #### Only massive
    keep = keep & (mcat.logm[mcat.idx] >= 11)
    
    #### Figure
    fig = unicorn.catalogs.plot_init(square=True)
    
    UV = -2.5*np.log10(rest.l153/rest.l155)
    VJ = -2.5*np.log10(rest.l155/rest.l161)
    
    keep_rest = keep & (rest.idx > -1) 
    detected_lines = (lines.halpha_eqw[lines.idx]/(1+zout.z_peak[0::3]) > 5) & (lines.halpha_eqw[lines.idx]/lines.halpha_eqw_err[lines.idx] > 3)
    
    plt.plot(VJ[rest.idx][keep_rest & detected_lines], UV[rest.idx][keep_rest & detected_lines], color='blue', marker='o', linestyle='None', alpha=0.2, markersize=6)
    plt.plot(VJ[rest.idx][keep_rest & ~detected_lines], UV[rest.idx][keep_rest & ~detected_lines], color='red', marker='o', linestyle='None', alpha=0.4, markersize=6)
    
    #### Compare colors of fit templates
    tempfilt, coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='rf_dummy', OUTPUT_DIRECTORY=unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/OUTPUT', CACHE_FILE = 'Same')
    
    uflux = tempfilt['tempfilt'][0,0:6,0]*1.
    vflux = tempfilt['tempfilt'][2,0:6,0]*1.
    jflux = tempfilt['tempfilt'][8,0:6,0]*1.

    uflux /= vflux
    jflux /= vflux
    vflux /= vflux
    
    uv_temp = -2.5*np.log10(uflux/vflux)
    vj_temp = -2.5*np.log10(vflux/jflux)
    plt.plot(vj_temp, uv_temp, marker='o', markersize=12, alpha=0.5, linestyle='None', color='black')
    
    s = np.argsort(vj_temp)
    s = np.append(s, s[0])
    
    for i in range(len(s)-1):
        uv_int = []
        vj_int = []
        for f in np.arange(0,1.01,0.01):
            uv_int.append( -2.5*np.log10((uflux[s][i]*f+uflux[s][i+1]*(1-f))/(vflux[s][i]*f+vflux[s][i+1]*(1-f))))
            vj_int.append( -2.5*np.log10((vflux[s][i]*f+vflux[s][i+1]*(1-f))/(jflux[s][i]*f+jflux[s][i+1]*(1-f))))
        #
        plt.plot(vj_int, uv_int, color='black', alpha=0.1, marker='None', linewidth=3)
    
    
    # dusty = (VJ[idx_rest] > 1.5) & (UV[idx_rest] > 1.5)
    # plt.plot(VJ[idx_rest][keep & zsp & dusty], UV[idx_rest][keep & zsp & dusty], color='green', marker='x', linestyle='None', alpha=0.8, markersize=6)
    # print rest.id[idx_rest][keep & zsp & dusty]
    
    plt.text(0,2.1,r'$1 < z < 1.5$', fontsize=13)
    plt.text(0,1.9,r'Split: EQW H$\alpha=5\AA$', fontsize=11)
    plt.xlim(-0.35,2.15)
    plt.ylim(0,2.5)
    plt.xlabel(r'$(V-J)$')
    plt.ylabel(r'$(U-V)$')
    
    fig.savefig('UVJ_v1.0.png')
    
def fill_image(objects, x, y, scale=None, xrange=(0,1), yrange=(0,1), NX=1024, NY=1024):
    
    im = np.zeros((NY, NX))
    NOBJ = len(objects)
    dx = xrange[1]-xrange[0]
    dy = yrange[1]-yrange[0]
    
    if scale is None:
        scale = np.ones(NOBJ)
        
    for i in range(NOBJ):
        print noNewLine+objects[i]
        thu = pyfits.open('DATA/'+objects[i]+'_thumb.fits.gz')
        shp = thu[0].data.shape
        
        xi = np.round((x[i]-xrange[0])/dx*NX)-shp[1]/2
        yi = np.round((y[i]-yrange[0])/dy*NY)-shp[0]/2
        
        x0, x1 = 0, shp[1]
        y0, y1 = 0, shp[0]
        
        if (xi+x1/2 < 0) | (xi > NX) | (yi+y1/2 < 0) | (yi > NY):
            print '\n', x[i], y[i], shp, '\n'
            continue
            
        if xi < 0:
            x0 = -xi
            xi = 0
        if (xi+shp[1]) > NX:
            x1 = NX-xi
            
        if yi < 0:
            y0 = -yi
            yi = 0
        if (yi+shp[0]) > NY:
            y1 = NY-yi
        
        print x0,x1,y0,y1,shp, xi, yi, scale[i]
        
        thu = thu[0].data[y0:y1, x0:x1]*scale[i] #/thu[0].data.sum()
        shp = thu.shape
        im[yi:yi+shp[0], xi:xi+shp[1]] += thu 
    
    return im
    
def composite_spectra(objects, color='red', alpha=0.1, lnorm=8.e3, NITER=3, show_lines=False):
    import unicorn.catalogs
    
    zout = unicorn.catalogs.zout #catIO.Readfile('full_redshift.cat')
    
    #objects = ['GOODS-S-24-G141_00459','GOODS-N-27-G141_00447']
    xm = np.array([3000,8000])
    ym = np.array([1,1])
    
    ### Iterate on normalizing to the average spectrum.
    for iter in range(1,NITER+1):
      avgx = []
      avgy = []
      for object in objects:
        print noNewLine+object
        zi = zout.z_peak[0::3][zout.id[0::3] == object][0]
        #
        try:
            data = np.loadtxt('DATA/'+object+'_obs_sed.dat')
        except:
            continue
        #           
        lc = data[:,0]
        fnu = data[:,1]
        efnu = data[:,2]
        #
        dlam_spec = lc[-1]-lc[-2]
        is_spec = np.append(np.abs(1-np.abs(lc[1:]-lc[0:-1])/dlam_spec) < 0.05,True)
        #
        flam = fnu/(lc/5500.)**2    
        s = np.argsort(lc)
        fint = np.interp(lnorm, lc[s]/(1+zi), flam[s])
        
        mask = is_spec & (np.abs(lc/(1+zi)-6563.) > 300)

        #mask = is_spec & (np.abs(lc/(1+zi)-lnorm) < 200)

        ymint = np.interp(lc[mask]/(1+zi), xm, ym)
        fnorm = np.sum(flam[mask]*ymint*efnu[mask]**2)/np.sum(ymint**2*efnu[mask]**2)
        
        #fnorm = fint
        
        flam /= fnorm
        eflam = efnu/(lc/5500.)**2/fnorm
        #
        if iter == NITER:
            #plt.plot(lc[is_spec]/(1+zi), flam[is_spec], color=color, alpha=alpha)
            #plt.plot(lc[~is_spec]/(1+zi), flam[~is_spec], color=color, alpha=alpha, marker='o', markersize=5, linestyle='None')
            pass
            
        avgx.extend(lc[is_spec]/(1+zi))
        avgy.extend(flam[is_spec])
    
      avgx = np.array(avgx)
      avgy = np.array(avgy)
      xm, ym, ys, N = threedhst.utils.runmed(avgx, avgy, NBIN=100)
      if iter == NITER:
          plt.plot(xm, ym, color='white', alpha=0.5, linewidth=4) #, drawstyle='steps-mid')
          plt.plot(xm, ym, color=color, alpha=0.8, linewidth=2) #, drawstyle='steps-mid')
    
    if show_lines:
        for line in [4861, 4959, 5007, 5178, 5891, 6563, 6585, 6718, 6731]:
            plt.plot(line*np.array([1,1]), [0,10], color='black', linestyle='--', alpha=0.5)
            
def plot_init(square=True, xs=6, aspect=1):
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = ['Times']
    plt.rcParams['patch.edgecolor'] = 'None'
    plt.rcParams['font.size'] = 10

    if square:
        #xs=5
        lrbt = np.array([0.22,0.02,0.11,0.02])*5./xs     
        ys = (1-lrbt[1]-lrbt[0])/(1-lrbt[3]-lrbt[2])*xs*aspect
        lrbt[[2,3]] /= aspect
        fig = plt.figure(figsize=(xs,ys), dpi=100)
        #fig.subplots_adjust(left=0.13,bottom=0.10,right=0.98,top=0.98)
        fig.subplots_adjust(left=lrbt[0],bottom=lrbt[2],right=1-lrbt[1],top=1-lrbt[3])
        # plt.plot([0,2])
        # plt.xlabel('x')
        # plt.ylabel('y')        
    else:
        fig = plt.figure(figsize=(7,5), dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.10,
                        bottom=0.10,right=0.99,top=0.97)        
    return fig
    
def make_object_tarfiles(objects, thumbs=False):
    """
    Make a script to get spectra and thumbnails from unicorn.
    """
    
    ###### Thumbnails
    if thumbs:
        os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/images')
        line = 'tar czvf thumbs.tar.gz'
        for object in objects:
            line += ' %s_thumb.fits.gz' %(object)
        
            os.system(line)
        
        os.system('mv thumbs.tar.gz /3DHST/Spectra/Work/ANALYSIS/FIRST_PAPER/')
    
    
    #### Get photometry + scaled spectra from the tempfilt files
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/FIRST_PAPER/GRISM_v1.6')
    
    tarline = 'tar czvf spec.tar.gz'
    for object in objects:
        print noNewLine+object
        #
        tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=object, OUTPUT_DIRECTORY='../../REDSHIFT_FITS/OUTPUT', CACHE_FILE = 'Same')
        
        #### Photometry + spectrum
        fp = open(object+'_obs_sed.dat','w')
        fp.write('# lc fnu efnu obs_sed \n')

        obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][0]],\
                         coeffs['coeffs'][:,0])# /(lci/5500.)**2
        
        for i in range(tempfilt['NFILT']):
            fp.write(' %8.4e %8.4e %8.4e  %8.4e\n' %(tempfilt['lc'][i], tempfilt['fnu'][i,0], tempfilt['efnu'][i,0], obs_sed[i]))
        
        fp.close()
        
        tarline += ' %s_obs_sed.dat' %(object)
        
        #### template fit
        zi = tempfilt['zgrid'][coeffs['izbest'][0]]
        lambdaz = temp_seds['templam']*(1+zi)
        temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,0])
        temp_sed *= (lambdaz/5500.)**2/(1+zi)**2
        
        fp = open(object+'_temp_sed.dat','w')
        fp.write('# lam fnu_temp\n')
        for i in range(temp_seds['NTEMPL']):
            fp.write(' %8.4e %8.4e\n' %(lambdaz[i], temp_sed[i]))
        fp.close()
        tarline += ' %s_temp_sed.dat' %(object)
        
        #
    os.system(tarline)
    
    print 'scp $UNICORN:'+unicorn.GRISM_HOME+'ANALYSIS/FIRST_PAPER/GRISM_V1.6/spec*tar.gz .'
    
def combine_sextractor_catalogs():
    os.chdir('/Users/gbrammer/Sites_GLOBAL/P/GRISM_v1.6/')

    files=glob.glob('*drz.cat')
    
    fp = open(files[0])
    lines = fp.readlines()
    fp.close()
    
    line = lines[0]
    columns = []
    i=0
    while line.startswith('#'):
        columns.append(line.split()[2])
        i+=1
        line = lines[i]
    
    header = '# ID FIELD POINTING '+' '.join(columns)+'\n# '+time.ctime()+'\n'
    out_lines = []
    for file in files:
        print noNewLine+file
        if file.startswith('GN20'):
            continue
        ### GOODS-N catalogs don't have the MAG_APER columns
        # if file.startswith('GOODS-N'):
        #     continue
        #
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        pointing = file.split('-G141')[0]
        field = re.split('-[1-9]',pointing)[0]
        if len(pointing.split(field+'-')) == 1:
            pointing_number == '1'
        else:
            pointing_number = pointing.split(field+'-')[1]
        #
        for line in lines:
            if not line.startswith('#'):
                number = int(line.split()[0])
                out_lines.append('%s-G141_%05d %s %3s ' %(pointing, number, field, pointing_number)+line)
    
    fp = open('ANALYSIS/full_sextractor.cat','w')
    fp.write(header)
    fp.writelines(out_lines)
    fp.close()
    
    ########## GOODS-N separate for now because doesn't have same columns
    # files=glob.glob('GOODS-N*drz.cat')
    # 
    # fp = open(files[0])
    # lines = fp.readlines()
    # fp.close()
    # 
    # line = lines[0]
    # columns = []
    # i=0
    # while line.startswith('#'):
    #     columns.append(line.split()[2])
    #     i+=1
    #     line = lines[i]
    # 
    # header = '# ID FIELD POINTING '+' '.join(columns)+'\n'
    # out_lines = []
    # for file in files:
    #     print noNewLine+file
    #     if file.startswith('GN20'):
    #         continue
    #     fp = open(file)
    #     lines = fp.readlines()
    #     fp.close()
    #     pointing = file.split('-G141')[0]
    #     field = re.split('-[1-9]',pointing)[0]
    #     if len(pointing.split(field+'-')) == 1:
    #         pointing_number = '1'
    #     else:
    #         pointing_number = pointing.split(field+'-')[1]
    #     #
    #     for line in lines:
    #         if not line.startswith('#'):
    #             number = int(line.split()[0])
    #             out_lines.append('%s-G141_%05d %s %3s ' %(pointing, number, field, pointing_number)+line)
    # 
    # fp = open('ANALYSIS/goodsn_sextractor.cat','w')
    # fp.write(header)
    # fp.writelines(out_lines)
    # fp.close()
    
    status = os.system('gzip ANALYSIS/*_sextractor.cat')
    
def make_full_redshift_catalog():
    """
    Cat all individual zout files into a single file
    """
    import copy
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS/')
    files = glob.glob('OUTPUT/*G141*zout')
    if len(files) == 0:
        os.system("ls OUTPUT/ |grep G141 |grep zout |awk '{print \"x\" $1 }' |sed \"s/x/OUTPUT\//\" > files.list")
        files = np.loadtxt('files.list', dtype=np.str)
        
    fp = open(files[0])
    lines = fp.readlines()
    lines[0]+='# '+time.ctime()+'\n'
    
    specphot_lines = copy.deepcopy(lines[0:3])
    phot_lines = copy.deepcopy(lines[0:3])
    spec_lines = copy.deepcopy(lines[0:3])
    
    fp.close()
    for file in files[1:]:
        print noNewLine+file
        fp = open(file)
        flines = fp.readlines()[2:]
        if len(flines) == 3:
            lines.extend(flines)
            specphot_lines.append(flines[0])
            phot_lines.append(flines[1])
            spec_lines.append(flines[2])
        #
        fp.close()
        
    fp = open('full_redshift.cat','w')
    fp.writelines(lines)
    fp.close()

    fp = open('specphot_redshift.cat','w')
    fp.writelines(specphot_lines)
    fp.close()

    fp = open('phot_redshift.cat','w')
    fp.writelines(phot_lines)
    fp.close()

    fp = open('spec_redshift.cat','w')
    fp.writelines(spec_lines)
    fp.close()
    
    if unicorn.hostname().startswith('uni'):
        status = os.system('cp *_redshift.cat /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
        status = os.system('gzip /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS/*_redshift.cat')
    
def combine_matched_catalogs():
    """
    Combine all of the match catalogs in the HTML/SED directories, adding
    the correct object id with the full pointing name.
    """
    
    os.chdir('/Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
    files = glob.glob('../SED/*match.cat')
    fp = open(files[0])
    full_lines = fp.readlines()[0:2]
    full_lines[0]+='# '+time.ctime()+'\n'
    fp.close()
    for file in files:
        pointing = os.path.basename(file).split('_match')[0]
        print noNewLine+pointing
        fp = open(file)
        lines = fp.readlines()[3:]
        fp.close()
        for line in lines:
            spl = line.split()
            id = int(spl[0])
            object = "%s_%05d  " %(pointing, id)
            full_lines.append(object+'  '.join(spl[1:])+'\n')
    
    fp = open('full_match.cat','w')
    fp.writelines(full_lines)
    fp.close()
    
    status = os.system('gzip full_match.cat')

def eqw_catalog():
    """ 
    Make a full catalog of the line fluxes / eq. widths
    """
    import unicorn.analysis
    
    os.chdir(unicorn.GRISM_HOME+'ANALYSIS/REDSHIFT_FITS')
    files=glob.glob('OUTPUT/*G141*.coeff')
    if len(files) == 0:
        os.system("ls OUTPUT/ |grep G141 |grep coeff |awk '{print \"x\" $1 }' |sed \"s/x/OUTPUT\//\" > files.list")
        files = np.loadtxt('files.list', dtype=np.str)

    lines = ['# id  z_grism halpha_eqw  halpha_eqw_err  halpha_flux   oiii_eqw oiii_eqw_err  oiii_flux   hbeta_eqw  hbeta_eqw_err  hbeta_flux\n# '+time.ctime()+'\n']
    for file in files:
        object=os.path.basename(file).split('.coeff')[0]
        print noNewLine+object
        root=object.split('G141')[0]+'G141'
        id = int(object.split('G141_')[1])
        #
        try:
            obj, z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = unicorn.analysis.equivalent_width(root=root, id=id)
        except:
            print '\n\nFail!\n\n'
            z_grism, halpha_eqw, halpha_err, halpha_flux, oiii_eqw, oiii_err, oiii_flux, hbeta_eqw, hbeta_err, hbeta_flux = -1,-1,-1,-1,-1,-1,-1,-1,-1,-1
        #
        lines.append('%s %8.3f %16.3e %11.3e %11.3e %16.3e %11.3e %11.3e %16.3e %11.3e %11.3e\n' %(object, z_grism,halpha_eqw,halpha_err,halpha_flux,oiii_eqw,oiii_err,oiii_flux,hbeta_eqw,hbeta_err,hbeta_flux))
    
    fp = open('full_emission_lines.cat','w')
    fp.writelines(lines)
    fp.close()
     
    status = os.system('cp full_emission_lines.cat /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS')
    status = os.system('gzip /Library/WebServer/Documents/P/GRISM_v1.6/ANALYSIS/full_emission_lines.cat')

def show_acs_spectra():
    import unicorn.catalogs
    from unicorn.catalogs import match_string_arrays
    import cosmocalc

    os.chdir(unicorn.GRISM_HOME+'/ANALYSIS/FIRST_PAPER/')

    ######################################
    #### Redshifts and masses
    ######################################
    
    zout = catIO.Readfile('full_redshift.cat')
    pointing = []
    field = []
    for obj in zout.id[0::3]:
        point = obj.split('-G141')[0]
        pointing.append(point)
        field.append(re.split('-[1-9]',point)[0])
    
    pointing = np.array(pointing)
    field = np.array(field)
    
    # mcat = catIO.Readfile('full_match.cat')
    #     
    # zsp = zout.z_spec[0::3] > 0
    # dz = (zout.z_peak-zout.z_spec)/(1+zout.z_spec)
    # found, idx = match_string_arrays(zout.id[0::3], mcat.id_f140w)
    # idx = idx[idx >= 0]
    # zsp = zsp & found
    # 
    # ######################################
    # #### Emission line fits
    # ######################################
    # 
    # lines = catIO.Readfile('full_emission_lines.cat')
    # found_lines, idx_lines = match_string_arrays(zout.id[0::3], lines.id)
    # #lines.halpha_eqw[idx_lines] /= (1+zout.z_peak[0::3][found_lines])
    # 
    # ##################################### 
    # #### Galfit
    # #####################################
    # gfit = catIO.Readfile('full_galfit.cat')
    # found_gfit_rev, idx_gfit_rev = match_string_arrays(gfit.object, zout.id[0::3])
    # found_gfit, idx_gfit = match_string_arrays(zout.id[0::3], gfit.object)
    # 
    # #### Make a grid of redshifts to compute the plate scale and then interpolate it.
    # zgrid = np.arange(100)/100.*4+1./100
    # scale = zgrid*0.
    # for i in range(100):
    #     cc = cosmocalc.cosmocalc(zgrid[i])
    #     scale[i] = cc['PS_kpc']
    # 
    # gfit.r_e_kpc = gfit.r_e*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    # gfit.r_e_kpc_err = gfit.r_e_err*0.06*np.interp(zout.z_peak[0::3][idx_gfit_rev], zgrid, scale)
    # gfit.r_e_kpc_circ = gfit.r_e_kpc * np.sqrt(1./np.abs(gfit.ba)) 
    # gfit.r_e_kpc_circ_err = gfit.r_e_kpc_err * np.sqrt(1./np.abs(gfit.ba)) 
    # 
    # ##### Selection slices
    # dr = mcat.rmatch[idx] < 1
    # zrange = (zout.z_peak[0::3] > 0.2)
    # 
    # keep = dr & zrange & (mcat.logm[idx] > 10.98) & (mcat.fcontam[idx] < 0.5)
    # keep = dr & zrange & (mcat.mag_f140w[idx] < 24) & (mcat.fcontam[idx] < 0.5)
    # 
    # #### refine
    # # zrange = (zout.z_peak[0::3] > 2) & (zout.z_peak[0::3] < 3)
    # # zrange = (zout.z_peak[0::3] > 1.5) & (zout.z_peak[0::3] < 2.)
    # # zrange = (zout.z_peak[0::3] > 0.4) & (zout.z_peak[0::3] < 1.0)
    # zrange = (zout.z_peak[0::3] > .0) & (zout.z_peak[0::3] < 1.5)
    # keep = dr & zrange & (mcat.fcontam[idx] < 0.2) & (zout.q_z[0::3] < 0.1) 
    
    #### Photometry with coords
    phot = catIO.Readfile('full_sextractor.cat')
    found_phot, idx_phot = match_string_arrays(zout.id[0::3], phot.id)
    
    #### ACS
    root = 'jbhm54020'
    acs = threedhst.sex.mySexCat('ACS/'+root+'_drz.cat')
    acs.ra = np.cast[float](acs.X_WORLD)
    acs.dec = np.cast[float](acs.Y_WORLD)

    acs_match = []
    acs_dr = np.arange(len(acs.NUMBER))+1.
    
    for i in range(len(acs.NUMBER)):
        dr = np.sqrt((acs.ra[i]-phot.x_world)**2*np.cos(acs.dec[i]/360.*2*np.pi)**2+(acs.dec[i]-phot.y_world)**2)*3600.
        ma = dr == dr.min()
        acs_match.append(phot.id[ma][0])
        acs_dr[i] = dr[ma][0]
    
    acs_match = np.array(acs_match)
    acs.mag = np.cast[float](acs.MAG_F806W)
    matches = (acs_dr < 1.0) & (acs.mag < 28)
    idx = np.arange(len(acs.mag))[matches] 
    
    cat, zcat, fout = unicorn.analysis.read_catalogs(root='COSMOS')
    
    OUTPUT_DIRECTORY = os.path.dirname(zcat.filename)
    MAIN_OUTPUT_FILE = os.path.basename(zcat.filename).split('.zout')[0]
    
    for i0, i in enumerate(idx):
        i = idx[i0]
        ra0, de0 = acs.ra[i], acs.dec[i]
        dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
        photom_idx = np.where(dr == np.min(dr))[0][0]
        drMatch = dr[photom_idx]*1.
        #if drMatch > 0.5: 
        #    continue
        
        ####
        acs_file = '%s_%05d' %(root, acs.id[i])
        wfc_file = acs_match[matches][i0]
        
        if not os.path.exists(acs_file+'.dat'):
            os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_ACS/ascii/%s.dat' %(acs_file))
        
        if not os.path.exists(wfc_file+'.dat'):    
            os.system('wget http://3dhst:getspecs@unicorn.astro.yale.edu/P/GRISM_v1.5/ascii/%s.dat' %(wfc_file))

        if (not os.path.exists(acs_file+'.dat')) | (not os.path.exists(wfc_file+'.dat')):
            continue
            
        lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(photom_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = 'Same')
    
        acs_spec = catIO.Readfile(acs_file+'.dat')
        wfc_spec = catIO.Readfile(wfc_file+'.dat')

        acs_spec.flux -= acs_spec.contam
        wfc_spec.flux -= wfc_spec.contam
        
        zi = zout.z_peak[::3][zout.id[::3] == wfc_file]
        
        if len(zi) == 0:
            continue
        else:
            zi = zi[0]
        
        keep_acs = (acs_spec.lam > 6000) & (acs_spec.lam < 9000) & (acs_spec.flux > 0)
        mask_acs = (np.abs(acs_spec.lam/(1+zi)-6563.) > 500) & (np.abs(acs_spec.lam/(1+zi)-5007.) > 500)
        
        keep_wfc = (wfc_spec.lam > 1.15e4) & (wfc_spec.lam < 1.62e4) & (wfc_spec.flux > 0)
        mask_wfc = (np.abs(wfc_spec.lam/(1+zi)-6563.) > 500) & (np.abs(wfc_spec.lam/(1+zi))-5007. > 500)
        
        yint_acs = np.interp(acs_spec.lam, lambdaz, temp_sed)
        norm_acs = np.sum((acs_spec.flux*yint_acs)[keep_acs & mask_acs])/np.sum(yint_acs[keep_acs & mask_acs]**2)

        yint_wfc = np.interp(wfc_spec.lam, lambdaz, temp_sed)
        norm_wfc = np.sum((wfc_spec.flux*yint_wfc)[keep_wfc & mask_wfc])/np.sum(yint_wfc[keep_wfc & mask_wfc]**2)
        
        xg = np.arange(1000)-500
        yg = np.exp(-xg**2/2/100**2)
        yg /= np.sum(yg)
        ll = np.arange(3.8e4)+2000
        yy = np.interp(ll, lambdaz, temp_sed)
        
        tconv = np.convolve(yy, yg, mode='same')
        
        fig = plt.figure(figsize=[5,3.7],dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.1,left=0.1,
                            bottom=0.15,right=0.98,top=0.98)
        
        ax = fig.add_subplot(111)
        
        zcorr = (1+zi)
        zcorr = 1.
        
        #plt.plot(lambdaz, temp_sed, color='black', alpha=0.4)
        ax.plot(ll/zcorr, tconv, color='black', alpha=0.3, linewidth=2)
        #plt.plot(ll, tconv, color='orange', alpha=0.4)
        ax.plot(wfc_spec.lam[keep_wfc]/zcorr, wfc_spec.flux[keep_wfc]/norm_wfc, color='red', linewidth=3, alpha=0.8)
        ax.plot(acs_spec.lam[keep_acs]/zcorr, acs_spec.flux[keep_acs]/norm_acs, color='blue', linewidth=3, alpha=0.8)
        ax.errorbar(lci/zcorr, fobs, yerr=efobs, marker='o', markersize=10, color='orange', ecolor='orange', linestyle='None', alpha=0.3)
        ax.text(0.95, 0.85, r'$z=%.2f$' %(zi), ha='right', transform = ax.transAxes, fontsize=20)
        ax.text(0.95, 0.1, r'$i_{814}=%5.2f$' %(acs.mag[i]), ha='right', transform = ax.transAxes, fontsize=16)
        ax.text(0.95, 0.2, r'$H_{140}=%5.2f$' %(phot.mag_f1392w[phot.id == wfc_file][0]), ha='right', transform = ax.transAxes, fontsize=16)
        
        #ax.set_xlabel(r'$\lambda_\mathrm{rest}\ [\AA]$')
        ax.set_xlabel(r'$\lambda_\mathrm{obs}\ [\AA]$')
        
        #from pylab import *
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter

        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2500)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)

        #for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)

        ax.set_xlim(3000,1.9e4)
        ymax_acs = np.max(acs_spec.flux[keep_acs]/norm_acs)
        ymax_wfc = np.max(wfc_spec.flux[keep_wfc]/norm_wfc)
        
        ymax = np.max([ymax_acs, ymax_wfc])
        
        ax.set_ylim(-0.1*ymax, 1.5*ymax)
        
        plt.savefig('%04.2f_' %(zi) + acs_file+'_example.pdf')
        plt.close()
        

def gen_rgb_colortable(Vbins = range(10), rgb = (1,0,0), reverse=False):
    values = np.arange(len(Vbins))*1./len(Vbins)
    if reverse:
        values = 1.-values
    Vcolors = []
    for i in range(len(Vbins)):
        Vcolors.append((rgb[0]*values[i], rgb[1]*values[i], rgb[2]*values[i]))
    #
    return Vcolors

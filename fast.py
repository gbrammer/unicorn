def go_all():
    import glob
    
    force=False
    files=glob.glob('FAST_OUTPUT/*fout')
    for file in files: 
        object = os.path.basename(file).split('_threedhst')[0]
        if (not os.path.exists(object+'_fast.png')) | force:
            check_fast(object=object)
        
def check_fast(object='AEGIS-1-G141_00497'):

    if object.startswith('GOODS-S') | object.startswith('WFC3') | object.startswith('GEORGE') | object.startswith('PRIMO'):
        abzp=23.86
    else:
        abzp = 25

    obs_sed = catIO.Readfile('ASCII/%s_obs_sed.dat' %(object))
    temp_sed = catIO.Readfile('ASCII/%s_temp_sed.dat' %(object))

    lc = obs_sed.lc
    dlam_spec = lc[-1]-lc[-2]
    is_spec = np.append(np.abs(1-np.abs(lc[1:]-lc[0:-1])/dlam_spec) < 0.05,True)

    obs_convert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2/10.**-19
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5, aspect=2./3, left=0.12)
    
    plt.plot(lc[~is_spec], obs_sed.fnu[~is_spec]*obs_convert[~is_spec], marker='o', color='orange', linestyle='None', markersize=15, alpha=0.7)
    plt.plot(lc[is_spec], obs_sed.fnu[is_spec]*obs_convert[is_spec], color='blue', alpha=0.5)
    plt.plot(lc[~is_spec], obs_sed.obs_sed[~is_spec]*obs_convert[~is_spec], color='red', alpha=0.7, marker='o', linestyle='None', markersize=8)

    temp_convert = 10**(-0.4*(abzp+48.6))*3.e18/temp_sed.lam**2/10.**-19
    flam_temp = temp_sed.fnu_temp*temp_convert

    plt.plot(temp_sed.lam, flam_temp, color='red', alpha=0.3)
    #fast_norm = 1./np.interp(2.2e4, wfast, tfast)*np.interp(2.2e4, temp_sed.lam, flam_temp)

    wfast, tfast = np.loadtxt('FAST_OUTPUT/BEST_FITS/%s_threedhst_1.fit' %(object), skiprows=1, unpack=True)
    plt.plot(wfast, tfast, color='green', alpha=0.5)

    wfast, tfast = np.loadtxt('FAST_OUTPUT/BEST_FITS/%s_threedhst_2.fit' %(object), skiprows=1, unpack=True)
    plt.plot(wfast, tfast, color='purple', alpha=0.5)
    
    #plt.semilogx()
    plt.xlim(2000,2.3e4)

    ymax = np.max(obs_sed.fnu*obs_convert)
    plt.ylim(-0.1*ymax,1.3*ymax)

    fout = catIO.Readfile('FAST_OUTPUT/%s_threedhst.fout' %(object))
    
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$f_\lambda\ (10^{-19}$)')
    xtext = 2e4
    xal = 'right'
    
    plt.text(xtext,0.4*ymax,r'$z_\mathrm{gris}=%.4f$' %(fout.z[0]), horizontalalignment=xal)

    plt.text(3000,1.1*ymax, object)
    
    plt.text(xtext,0.3*ymax,r'log M: $%.1f^{\ %.1f}_{\ %.1f}$   $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.lmass[0], fout.u68_lmass[0], fout.l68_lmass[0], fout.lmass[1], fout.u68_lmass[1], fout.l68_lmass[1]), horizontalalignment=xal)

    plt.text(xtext,0.2*ymax,r'$A_V$: $%.1f^{\ %.1f}_{\ %.1f}$   $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.av[0], fout.u68_av[0], fout.l68_av[0], fout.av[1], fout.u68_av[1], fout.l68_av[1]), horizontalalignment=xal)
    
    plt.text(xtext,0.1*ymax,r'log $\tau$: $%.1f^{\ %.1f}_{\ %.1f}$    $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.ltau[0], fout.u68_ltau[0], fout.l68_ltau[0], fout.ltau[1], fout.u68_ltau[1], fout.l68_ltau[1]), horizontalalignment=xal)
    
    plt.text(xtext,0.0*ymax,r'log Age: $%.1f^{\ %.1f}_{\ %.1f}$    $%.1f^{\ %.1f}_{\ %.1f}$' %(fout.lage[0], fout.u68_lage[0], fout.l68_lage[0], fout.lage[1], fout.u68_lage[1], fout.l68_lage[1]), horizontalalignment=xal)
    
    fig.savefig(object+'_fast.png')
    plt.close()
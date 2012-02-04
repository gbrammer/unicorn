"""
Test suite for Cython utilities.
"""
import numpy as np

def interp_conserve(x, xp, fp, left=0., right=0.):
    """
    Interpolate `xp`,`yp` array to the output x array, conserving flux.  
    `xp` can be irregularly spaced.
    """
    midpoint = (x[1:]-x[:-1])/2.+x[:-1]
    midpoint = np.append(midpoint, np.array([x[0],x[-1]]))
    midpoint = midpoint[np.argsort(midpoint)]
    int_midpoint = np.interp(midpoint, xp, fp, left=left, right=right)
    int_midpoint[midpoint > xp.max()] = 0.
    int_midpoint[midpoint < xp.min()] = 0.
    
    fullx = np.append(xp, midpoint)
    fully = np.append(fp, int_midpoint)
    
    so = np.argsort(fullx)
    fullx, fully = fullx[so], fully[so]
    
    outy = x*0.
    dx = midpoint[1:]-midpoint[:-1]
    for i in range(len(x)):
        bin = (fullx >= midpoint[i]) & (fullx <= midpoint[i+1])
        outy[i] = np.trapz(fully[bin], fullx[bin])/dx[i]
    
    return outy

def test_nmf():
    
    import interp_c
    import numpy as np
    import time
    import matplotlib.pyplot as plt
    
    x = np.arange(0,16,0.5)*np.pi*2
    NFILT = len(x)
    NTEMP = 7
    
    coeffs = np.random.random(NTEMP)
    
    y = x*0.
    templates = np.zeros((NTEMP, NFILT))
    norms = np.zeros(NTEMP)
    for i in range(0,NTEMP):
        templates[i,:] = np.sin(x/((i+1)*1.6))+1.1
        norms[i] = np.trapz(templates[i,:], x)
        
    y = np.dot(coeffs.reshape((1,NTEMP)), templates).flatten()
    
    err = y*0.+np.median(y)/10.
    
    yerr = y+np.random.normal(size=NFILT)*err
    
    t0 = time.time()
    
    amatrix = interp_c.prepare_nmf_amatrix(err**2, templates)
    t1= time.time()
    coeffs_fit = interp_c.run_nmf(yerr, err**2, templates, amatrix, toler=1.e-5)
    t2 = time.time()
    
    print 'Prepare: %.4f' %(t1-t0)
    print 'Fit    : %.4f' %(t2-t1)
    
    yfit = np.dot(coeffs_fit.reshape((1,-1)), templates).flatten()
    
    fig = plt.figure()
    
    ax = fig.add_subplot(311)
    ax.plot(x, y/np.median(y), color='blue')
    ax.errorbar(x, yerr/np.median(y), err/np.median(y), color='blue', marker='o', linestyle='None')
    ax.plot(x, yfit/np.median(y), color='red')
    ax.set_ylabel('"Flux"')
        
    ax = fig.add_subplot(312)
    ax.plot(x, y-yfit, color='blue')
    ax.errorbar(x, yerr-yfit, err, color='blue', marker='o', linestyle='None')
    ax.plot(x, yfit-yfit, color='red')
    ax.set_ylabel(r'$\Delta$(obs-fit)')
    chi2 = np.sum((yerr-yfit)**2/err**2)/(len(x)-1)
    ax.text(0.1,0.8,r'$\chi^2_\nu$=%.3f' %(chi2), transform=ax.transAxes)
    
    ax = fig.add_subplot(313)
    ax.plot(np.log10(coeffs/coeffs_fit), color='orange')
    ax.set_ylabel(r'$\log$($\Delta$coeff)')
    
    #### Is sum of normalizations conserved?
    norm_input = np.sum(norms*coeffs)
    norm_fit = np.sum(norms*coeffs_fit)
    int_fit = np.trapz(yfit, x)
    print 'Norm_in: %.2f, Norm_fit: %.2f, trapz_fit: %.2f' %(norm_input, norm_fit, int_fit)
    
    fig.savefig('/tmp/nmf.png')
    
def test():

    import interp_c
    import time
    import scipy
    import threedhst
    import numpy as np
    
    N = int(1.e6)
    
    xfull = np.arange(0,N+1,1)*1.
    #yfull = np.sin(xfull/(N/1239.)*2*np.pi)+1
    yfull = np.sin(xfull/np.pi/2/20)+0.2

    # coeffs = np.random.random(size=12)*5
    # yfull = scipy.polyval(coeffs, xfull)
    
    xint = np.arange(0,N+1,N/100)*1.
    
    tstart = time.time()
    denom = np.trapz(yfull,xfull)
    
    tstart = time.time()
    yint_0 = np.interp(xint, xfull, yfull)
    t0 = time.time()
    print 'Linear           : %.3f   (%.4e)' %(t0-tstart, np.trapz(yint_0, xint)/denom-1)

    yint_x = interp_c.interp_c(xint, xfull, yfull)
    tx = time.time()
    print 'Linear(c)        : %.3f   (%.4e)' %(tx-t0, np.trapz(yint_x, xint)/denom-1)
    
    xreverse = xint[::-1]
    yint_y = interp_c.interp_c(xreverse, xfull, yfull, assume_sorted=0)
    ty = time.time()
    print 'Linear(c) rev    : %.3f   (%.4e)' %(ty-tx, np.trapz(yint_y, xint)/denom-1)
    
    yint_1 = threedhst.utils.interp_conserve(xint, xfull, yfull)
    t1 = time.time()
    print 'Conserve         : %.3f   (%.4e)' %(t1-ty, np.trapz(yint_1, xint)/denom-1)
    
    yint_2 = interp_c.interp_conserve(xint, xfull, yfull)
    t2 = time.time()
    print 'Conserve (Cython): %.3f   (%.4e)' %(t2-t1, np.trapz(yint_2, xint)/denom-1)

    yint_3 = interp_c.interp_conserve_c(xint, xfull, yfull)
    t3 = time.time()
    print 'Conserve (more c): %.3f   (%.4e)' %(t3-t2, np.trapz(yint_3, xint)/denom-1)

    yint_4 = threedhst.utils.interp_conserve_c(xint, xfull, yfull)
    t4 = time.time()
    print 'Inline c         : %.3f   (%.4e)' %(t4-t3, np.trapz(yint_4, xint)/denom-1)

    #### Test interpolation
    threedhst.showMessage('Interpolation')
    
    #### Faster while n(int)/n(full) < 1./50
    
    xint = xfull[1000:-1000:40]

    tstart = time.time()    
    yint = np.interp(xint, xfull, yfull, left=0., right=0.)
    t0 = time.time()
    print 'Python         : %.4f' %(t0-tstart)

    yint1 = interp_c.interp_c(xint, xfull, yfull, extrapolate=0.)
    t1 = time.time()
    print 'Cython rewrite : %.4f   (%.2e)' %(t1-t0, np.sum((yint1-yint)**2))
    
    #### Test midpoint definition --- slices work better than by hand

    threedhst.showMessage('Midpoint')
    xmid = xfull
    
    tstart = time.time()
    midpoint = (xmid[1:]+xmid[:-1])/2.
    midpoint = np.append(midpoint, np.array([xmid[0],xmid[-1]]))
    midpoint = midpoint[np.argsort(midpoint)]
    t0 = time.time()
    print 'Python      :  %.3f  %.2e'   %(t0-tstart, np.sum((midpoint-midpoint)**2))

    midpoint_c1 = interp_c.midpoint(xmid)
    t1 = time.time()
    print 'Cython      :  %.3f  %.2e'   %(t1-t0, np.sum((midpoint_c1-midpoint)**2))

    midpoint_c2 = interp_c.midpoint_c(xmid, N+1)
    t2 = time.time()
    print 'Cython (opt):  %.3f  %.2e'   %(t2-t1, np.sum((midpoint_c2-midpoint)**2))
    
    
if __name__ == "__main__":
    test()
    
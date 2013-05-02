"""
Pythonic utilities ported to C [Cython] for speedup.
"""
import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t

import cython
    
@cython.boundscheck(False)
def interp_c(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] xp, np.ndarray[DTYPE_t, ndim=1] fp, double extrapolate=0., short int assume_sorted=1):
    """
    interp_c(x, xp, fp, extrapolate=0., assume_sorted=0)
    
    Fast interpolation: [`xp`, `fp`] interpolated at `x`.
    
    Extrapolated values are set to `extrapolate`.
    
    The default `assume_sorted`=1 assumes that the `x` array is sorted and single-
    valued, providing a significant gain in speed. (xp is always assumed to be sorted)
    
    """
    cdef unsigned long i, j, N, Np
    cdef DTYPE_t x1,x2,y1,y2,out
    cdef DTYPE_t fout, xval, xmin
    
    N, Np = len(x), len(xp)
    cdef np.ndarray[DTYPE_t, ndim=1] f = np.zeros(N)

    i=0
    j=0
    ### Handle left extrapolation
    xmin = xp[0]    
    if assume_sorted == 1:
        while x[j] < xmin: 
            f[j] = extrapolate
            j+=1
        
    while j < N:
        xval = x[j]
        if assume_sorted == 0:
            if x[j] < xmin:
                f[j] = extrapolate
                j+=1
                continue
            else:
                i=0
        while (xp[i] < xval) & (i < Np-1): i+=1;
        if i == (Np-1):
            if x[j] != xp[i]:
                f[j] = extrapolate
            else:
                f[j] = fp[i]
            j+=1
            continue   
        #### x[i] is now greater than xval because the 
        #### expression (x[i]<xval) is false, assuming
        #### that xval < max(x).
        
        x1 = xp[i];
        x2 = xp[i+1];
        y1 = fp[i];
        y2 = fp[i+1];
        out = ((y2-y1)/(x2-x1))*(xval-x1)+y1;
        f[j] = out
        j+=1
                
    return f
    
#
def interp_conserve(x, xp, fp, left=0., right=0.):
    """
    Interpolate `xp`,`yp` array to the output x array, conserving flux.  
    `xp` can be irregularly spaced.
    """
    cdef np.ndarray fullx, fully, so, outy, dx
    cdef long N, i
    
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

@cython.boundscheck(False)
def interp_conserve_c(np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] tlam, np.ndarray[DTYPE_t, ndim=1] tf, double left=0, double right=0):
    """
    interp_conserve_c(x, xp, fp, left=0, right=0)
    
    Interpolate `xp`,`yp` array to the output x array, conserving flux.  
    `xp` can be irregularly spaced.
    """
    cdef np.ndarray[DTYPE_t, ndim=1] templmid
    cdef np.ndarray[DTYPE_t, ndim=1] tempfmid
    cdef np.ndarray[DTYPE_t, ndim=1] outy
    cdef unsigned long i,k,istart,ntlam,NTEMPL
    cdef DTYPE_t h, numsum
    
    # templmid = (x[1:]+x[:-1])/2. #2.+x[:-1]
    # templmid = np.append(templmid, np.array([x[0], x[-1]]))
    # templmid = templmid[np.argsort(templmid)]
    NTEMPL = len(x)
    ntlam = len(tlam)

    templmid = midpoint_c(x, NTEMPL)
    #tempfmid = np.interp(templmid, tlam, tf, left=left, right=right)
    tempfmid = interp_c(templmid, tlam, tf, extrapolate=0.)
    
    outy = np.zeros(NTEMPL, dtype=DTYPE)

    ###### Rebin template grid to master wavelength grid, conserving template flux
    i=0
    for k in range(NTEMPL):
        numsum=0.;

        #### Go to where tlam is greater than the first midpoint
        while (tlam[i] < templmid[k]) & (i < ntlam): i+=1;
        istart=i;

        ####### First point
        if tlam[i] < templmid[k+1]: 
            h = tlam[i]-templmid[k];
            numsum+=h*(tf[i]+tempfmid[k]);
            i+=1;

        if i==0: i+=1;

        ####### Template points between master grid points
        while (tlam[i] < templmid[k+1]) & (i < ntlam):
            h = tlam[i]-tlam[i-1];
            numsum+=h*(tf[i]+tf[i-1]);
            i+=1;

        #### If no template points between master grid points, then just use interpolated midpoints
        if i == istart:
            h = templmid[k+1]-templmid[k];
            numsum=h*(tempfmid[k+1]+tempfmid[k]);
        else:  
            ##### Last point              
            if (templmid[k+1] == tlam[i]) & (i < ntlam):
                h = tlam[i]-tlam[i-1];
                numsum+=h*(tf[i]+tf[i-1]);
            else:
                i-=1;
                h = templmid[k+1]-tlam[i];
                numsum+=h*(tempfmid[k+1]+tf[i]);

        outy[k] = numsum*0.5/(templmid[k+1]-templmid[k]);
    
    return outy
    
def midpoint(x):
    mp = (x[1:]+x[:-1])/2.
    mp = np.append(mp, np.array([x[0],x[-1]]))
    mp = mp[np.argsort(mp)]
    return mp

def midpoint_c(np.ndarray[DTYPE_t, ndim=1] x, long N):
    cdef long i
    cdef DTYPE_t xi,xi1
    # N = len(x)
    cdef np.ndarray[DTYPE_t, ndim=1] midpoint = np.zeros(N+1, dtype=DTYPE)
    midpoint[0] = x[0]
    midpoint[N] = x[N-1]
    xi1 = x[0]
    for i in range(1, N):
        xi = x[i]
        midpoint[i] = 0.5*xi+0.5*xi1
        xi1 = xi
        
    return midpoint
    
@cython.boundscheck(False)
def prepare_nmf_amatrix(np.ndarray[DTYPE_t, ndim=1] variance, np.ndarray[DTYPE_t, ndim=2] templates):
    """
    prepare_nmf_amatrix(variance, templates)
    
    Generate the "A" matrix needed for the NMF fit, which is essentially 
    T.transpose() dot T.  This function is separated from the main fitting routine 
    because it does not depend on the actual measured "flux", which the user
    might want to vary independent of the variance
    
    `templates` (T) is a 2D matrix of size (NTEMPLATE, NBAND) in terms of photo-z
    fitting ofphotometric bands.
    
    `variance` is an array with size (NBAND) representing the *measured* variance in 
    each band.  
    
    --- Cythonified from eazy/getphotz.c (G. Brammer et al. 2008) ---
    
    """
    
    cdef np.ndarray[DTYPE_t, ndim=2] amatrix
    cdef unsigned int i,j,k,NTEMP,NFILT
    
    NTEMP, NFILT = np.shape(templates)
    amatrix = np.zeros((NTEMP,NTEMP))

    for i in range(NTEMP):
        for j in range(NTEMP):
            amatrix[i,j] = 0.
            for k in range(NFILT):
                amatrix[i,j]+=templates[i,k]*templates[j,k]/variance[k]
            
    return amatrix
    
@cython.boundscheck(False)
def run_nmf(np.ndarray[DTYPE_t, ndim=1] flux, np.ndarray[DTYPE_t, ndim=1] variance, np.ndarray[DTYPE_t, ndim=2] templates, np.ndarray[DTYPE_t, ndim=2] amatrix, double toler=1.e-4, long MAXITER=100000, init_coeffs=1, verbose=False):
    """
    run_nmf(flux, variance, templates, amatrix, toler=1.e-4, MAXITER=100000, verbose=False)
    
    Run the "NMF" fit to determine the non-negative coefficients of the `templates`
    matrix that best-fit the observed `flux` and `variance` arrays.
    
    `amatrix` is generated with the `prepare_nmf_amatrix` function.
    
    e.g.
    
    >>> coeffs = run_nmf(flux, variance, templates, amatrix)
    >>> flux_fit = np.dot(coeffs.reshape((1,-1)), templates).flatten()
    >>> chi2 = np.sum((flux-flux_fit)**2/variance)
    
    --- Cythonified from eazy/getphotz.c (G. Brammer et al. 2008) ---
    """
    cdef unsigned long i,j,k,itcount,NTEMP
    cdef double tolnum,toldenom,tol
    cdef double vold,av
    cdef np.ndarray[DTYPE_t, ndim=1] bvector
    
    cdef extern from "math.h":
        double fabs(double)

    NTEMP, NFILT = np.shape(templates)
            
    #### Make Bvector
    bvector = np.zeros(NTEMP, dtype=DTYPE)
    for i in range(NTEMP):
        bvector[i] = 0.;
        for k in range(NFILT):
            bvector[i]+=flux[k]*templates[i,k]/variance[k];
    
    #### Fit coefficients
    cdef np.ndarray[DTYPE_t, ndim=1] coeffs = np.ones(NTEMP, dtype=DTYPE)*init_coeffs
    tol = 100

    itcount=0;
    while (tol>toler) & (itcount<MAXITER):
        tolnum=0.
        toldenom = 0.
        tol=0
        for i in range(NTEMP):
            vold = coeffs[i];
            av = 0.;
            for j in range(NTEMP): av+=amatrix[i,j]*coeffs[j];
            #### Update coeffs in place      
            coeffs[i]*=bvector[i]/av;
            #tolnum+=np.abs(coeffs[i]-vold);
            tolnum+=fabs(coeffs[i]-vold);
            toldenom+=vold;

        tol = tolnum/toldenom;
        
        if verbose & 2:
            print 'Iter #%d, tol=%.2e' %(itcount, tol)
        
        itcount+=1
    
    if verbose & 1:
        print 'Iter #%d, tol=%.2e' %(itcount, tol)

    return coeffs
    
#    
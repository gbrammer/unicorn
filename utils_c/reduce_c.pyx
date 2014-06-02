"""
Pythonic utilities ported to C [Cython] for speedup.
"""

from __future__ import division

import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t
ctypedef np.int_t INT_t

import cython

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.embedsignature(True)
def disperse_grism_object(np.ndarray[DTYPE_t, ndim=2] flux, unsigned int id, np.ndarray[UINT_t, ndim=2] segm, np.ndarray[INT_t, ndim=1] xord, np.ndarray[INT_t, ndim=1] yord, np.ndarray[DTYPE_t, ndim=1] ford, np.ndarray[DTYPE_t, ndim=2] object, xpix=None, ypix=None):
    """
    disperse_grism_object(flux, id, segm, xord, yord, ford, object)
    
    The orders of the grism spectrum are defined in xord, yord, ford.  For a given object, `id`, loop through the pixels within the `segm` segmentation image scaling the dispersed spectrum by the `flux` value in each pixel.  
    
    The result is updated in place to the `object` image.
    
    """
    cdef unsigned int NX, NY, x, y, Nord, iord, x0, x1, y0, y1
    cdef int xxi, yyi
    #cdef np.ndarray[DTYPE_t, ndim=2] object
    cdef double flux_i
            
    NY, NX = np.shape(segm)
    Nord = np.shape(xord)[0]
    
    if xpix is not None:
        x0, x1 = xpix
    else:
        x0, x1 = 0, NX
        
    if ypix is not None:
        y0, y1 = ypix
    else:
        y0, y1 = 0, NY
    
    #object = np.zeros((NY,NX), dtype=DTYPE)
    
    for x in range(x0, x1):
        for y in range(y0, y1):
            if segm[y,x] == id:
                flux_i = flux[y,x]
                for iord in range(Nord):
                    xxi = x+xord[iord]
                    yyi = y+yord[iord]
                    if (xxi >= 0) & (xxi < NX) & (yyi >= 0) & (yyi < NY):
                        object[yyi,xxi] += ford[iord]*flux_i
                        
    return True

@cython.boundscheck(False)
@cython.embedsignature(True)
def total_flux(np.ndarray[DTYPE_t, ndim=2] flux, np.ndarray[UINT_t, ndim=2] segm):
    """
    total_flux(flux, segm)
    
    `segm` is a segmentation image with intensities defined in the `flux` image`.  
    
    Return an array, `total_flux`, where
        total_flux[i] = sum(flux[where(segm == i)])
    
    Note that there will be zero-valued entries in the output arrays where
    there were no pixels with that index number in the segmentation image.
    """
    cdef unsigned int NX, NY, x, y
    
    cdef unsigned int id_max
    cdef np.ndarray[DTYPE_t, ndim=1] total_fluxes
    
    id_max = segm.max()
    total_fluxes = np.zeros(id_max+1, dtype=np.double)
    
    NY, NX = np.shape(segm)        
    
    for x in range(NX):
        for y in range(NY):
            total_fluxes[segm[y,x]] += flux[y,x]
    
    return total_fluxes

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.embedsignature(True)
def get_model_ratio(np.ndarray[DTYPE_t, ndim=2] object, np.ndarray[DTYPE_t, ndim=2] model, np.ndarray[DTYPE_t, ndim=2] observed):
    
    cdef unsigned int NX, NY, ix, iy, x
    cdef np.ndarray[DTYPE_t, ndim=1] ratio_extract
    cdef double wht_sum, wht_i, obj_i, obj_sum, observed_i, model_i
    
    NY, NX = np.shape(object)
    ratio_extract = np.zeros(NX, dtype=DTYPE)
    
    for x in range(NX):
        wht_sum = 0.
        obj_sum = 0.
        for y in range(NY):
            obj_i = object[y,x]    
            observed_i = observed[y,x]
            model_i = model[y,x]
            if (obj_i > 0) & ((observed_i-model_i) > 0):
                wht_i = obj_i**4
                obj_sum += (observed_i-model_i)/obj_i*wht_i
                wht_sum += wht_i
        
        if wht_sum > 0:
            ratio_extract[x] = obj_sum/wht_sum
        
    return ratio_extract
    
#
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.embedsignature(True)
def get_model_ratio_optimal(np.ndarray[DTYPE_t, ndim=2] object, np.ndarray[DTYPE_t, ndim=2] model, np.ndarray[DTYPE_t, ndim=2] observed, np.ndarray[DTYPE_t, ndim=2] error):
    
    cdef unsigned int NX, NY, ix, iy
    cdef np.ndarray[DTYPE_t, ndim=1] ratio_extract
    cdef double obj_i, observed_i, var_i, wht_i
    cdef double data_num_sum, obj_num_sum, denom_sum, prof_sum
    
    NY, NX = np.shape(object)
    ratio_extract = np.zeros(NX, dtype=DTYPE)
        
    for x in range(NX):
        wht_sum = 0.
        obj_sum = 0.
        data_num_sum = 0.
        obj_num_sum = 0.
        #denom_sum = 0.
        #prof_sum = 0.
        for y in range(NY):
            obj_i = object[y,x]    
            observed_i = observed[y,x]
            var_i = error[y,x]**2
            if (obj_i > 0) & (observed_i != 0) & (var_i > 0):
                #wht_i = obj_i**5/var_i
                data_num_sum += obj_i*(observed_i-model[y,x])/var_i
                obj_num_sum += obj_i*obj_i/var_i
                #denom_sum += obj_i**2/var_i
                #prof_sum += obj_i
        
        if obj_num_sum > 0:
            #ratio_extract[x] = (num_sum) / (denom_sum)
            ratio_extract[x] = (data_num_sum) / (obj_num_sum)
        
    return ratio_extract

#### Compute IGM absorption factors as in EAZY, Hyperz, etc.
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.embedsignature(True)
def getigmfactors(np.ndarray[DTYPE_t, ndim=1] ztarg, np.ndarray[DTYPE_t, ndim=1] da, np.ndarray[DTYPE_t, ndim=1] db):
    """
    Get Madau IGM da/db extinction factors
    """
    
    cdef double l2=1216.0, l3=1026.0, l4=973.0, l5=950.0
    cdef double a2=3.6e-3, a3=1.7e-3, a4=1.2e-3, a5=9.3e-4
    cdef double a,b,c,d, madau_sum
    cdef double daz, dbz, zi
    
    cdef double lam1, lam2
    cdef double dl
    
    #### Extend further down lyman series
    cdef int i, N
    cdef np.ndarray[DTYPE_t, ndim=1] ll, aa
    ll = np.zeros(16, dtype=np.double)
    aa = np.zeros(16, dtype=np.double)
    
    ll[0] = 1216.0
    aa[0] =   3.6e-03
    ll[1] = 1026.0
    aa[1] =   1.7e-03
    ll[2] = 972.8
    aa[2] =   1.2e-03
    ll[3] = 950.0
    aa[3] =   9.4e-04
    ll[4] = 938.1
    aa[4] =   8.2e-04
    ll[5] = 931.0
    aa[5] =   7.5e-04
    ll[6] = 926.5
    aa[6] =   7.1e-04
    ll[7] = 923.4
    aa[7] =   6.8e-04
    ll[8] = 921.2
    aa[8] =   6.6e-04
    ll[9] = 919.6
    aa[9] =   6.4e-04
    ll[10] = 918.4
    aa[10] =   6.3e-04
    ll[11] = 917.5
    aa[11] =   6.2e-04
    ll[12] = 916.7
    aa[12] =   6.1e-04
    ll[13] = 916.1
    aa[13] =   6.0e-04
    ll[14] = 915.6
    aa[14] =   6.0e-04
    ll[15] = 915.2
    aa[15] =   6.0e-04
    
   
    N = ztarg.shape[0]
    #da = np.zeros(N, dtype=np.double)
    #db = np.zeros(N, dtype=np.double)
    
    for i in range(N):
        #
        #print i
        zi = ztarg[i]
        if (zi <= 0.0):
            da[i] = 1
            db[i] = 1
            #continue
        #    
        lam1 =  (1050.0*(1+zi)+0.5)
        lam2 =  (1170.0*(1+zi)+0.5)
        #
        daz = 0.0
        #for (dl=lam1;dl<=lam2; dl+=1/4.) {
        for dl in np.arange(lam1, lam2, 1./4):
            a = exp(-a2*(dl/l2)**3.46)
            daz += a
        #
        da[i] = 1.0-daz/((lam2-lam1)*4.)
        #
        lam1  = (915.0*(1+zi)+0.5)
        lam2 =  (1015.0*(1+zi)+0.5)
        #
        dbz = 0
        #for (dl=lam1;dl<=lam2; dl+=1/4.)  {
        for dl in np.arange(lam1, lam2, 1./4):
            # a=a3*(dl/l3)**3.46
            # b=a4*(dl/l4)**3.46
            # c=a5*(dl/l5)**3.46
            # d=exp(-(a+b+c))
            madau_sum = 0.
            #for (i=1;i<16;++i)
            for j in range(1, 16):
                madau_sum += aa[j]*(dl/ll[j])**3.46
            #
            dbz += exp(-1*madau_sum)
        #
        db[i] = 1.0-dbz/((lam2-lam1)*4.)
    
    return 0

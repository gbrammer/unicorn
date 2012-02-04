"""
Pythonic utilities ported to C [Cython] for speedup.
"""
import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t

import cython

def disperse_grism_object(np.ndarray[DTYPE_t, ndim=2] flux, unsigned int id, np.ndarray[UINT_t, ndim=2] segm, np.ndarray[UINT_t, ndim=1] xord, np.ndarray[UINT_t, ndim=1] yord, np.ndarray[DTYPE_t, ndim=1] ford, np.ndarray[DTYPE_t, ndim=2] object):
    """
    disperse_grism_object(flux, id, segm, xord, yord, ford, object)
    
    The orders of the grism spectrum are defined in xord, yord, ford.  For a given object, `id`, loop through the pixels within the `segm` segmentation image scaling the dispersed spectrum by the `flux` value in each pixel.  
    
    The result is updated in place to the `object` image.
    
    """
    cdef unsigned int NX, NY, x, y, Nord, iord, xxi, yyi
    #cdef np.ndarray[DTYPE_t, ndim=2] object
    cdef double flux_i
    
    NY, NX = np.shape(segm)
    Nord = np.shape(xord)[0]
    
    #object = np.zeros((NY,NX), dtype=DTYPE)
    
    for x in range(NX):
        for y in range(NY):
            if segm[y,x] == id:
                flux_i = flux[y,x]*10.
                for iord in range(Nord):
                    xxi = x+xord[iord]
                    yyi = y+yord[iord]
                    if (xxi >= 0) & (xxi < NX) & (yyi >= 0) & (yyi < NY):
                        object[yyi,xxi] += ford[iord]*flux_i
                        
    return True

@cython.boundscheck(False)
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
def get_model_ratio(np.ndarray[DTYPE_t, ndim=2] object, np.ndarray[DTYPE_t, ndim=2] model, np.ndarray[DTYPE_t, ndim=2] observed):
    
    cdef unsigned int NX, NY, ix, iy
    cdef np.ndarray[DTYPE_t, ndim=1] ratio_extract
    cdef double wht_sum, wht_i, obj_i, obj_sum, observed_i
    
    NY, NX = np.shape(object)
    ratio_extract = np.zeros(NX, dtype=DTYPE)
    
    for x in range(NX):
        wht_sum = 0.
        obj_sum = 0.
        for y in range(NY):
            obj_i = object[y,x]    
            observed_i = observed[y,x]
            if (obj_i > 0) & (observed_i != 0):
                wht_i = obj_i**4
                obj_sum += (observed_i-model[y,x])/obj_i*wht_i
                wht_sum += wht_i
        
        if wht_sum > 0:
            ratio_extract[x] = obj_sum/wht_sum
        
    return ratio_extract
    

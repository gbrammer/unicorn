"""
Compute dust reddening.  Currently just implemented Calzetti (2000) law
"""
import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.uint_t UINT_t
ctypedef np.int_t INT_t

import cython

###### Exactly from Hyperz
def redden(np.ndarray[DTYPE_t, ndim=1] wl,np.ndarray[DTYPE_t, ndim=1] fl, DTYPE_t av, rv=4.05):
    """
    Apply a Calzetti et al. (2000) reddening curve to a spectrum (`wl`, `fl`)
    with A_v=`av`.  The function is translated directly from the FORTRAN code
    of HYPERZ (Bolzonella et al. 2000).
    """
    cdef unsigned int ik, NWL
    cdef double p11,f11,p12,ff12,slope1,ff99,ff100,slope2
    cdef double ala, p, ff, alambda
    cdef double Rv
    
    Rv = rv
    NWL = len(wl)
    
    yob = np.zeros(NWL, dtype=DTYPE)
    
    #rv=4.05
    p11=1./0.11
    ff11=2.659*(-2.156+1.509*p11-0.198*p11**2+0.011*p11**3)+Rv
    p12=1./0.12
    ff12=2.659*(-2.156+1.509*p12-0.198*p12**2+0.011*p12**3)+Rv
    slope1=(ff12-ff11)/100.
    ff99=2.659*(-1.857+1.040/2.19)+Rv
    ff100=2.659*(-1.857+1.040/2.2)+Rv
    slope2=(ff100-ff99)/100.
    
    for ik in range(NWL):
        ala=wl[ik]*1.e-4 #### wavelength in microns
        p=1./ala
        if (ala >= 0.63) & (ala <= 2.2):
           ff=2.659*(-1.857+1.040*p)+Rv
        if (ala >= 0.12) & (ala < 0.63):
           ff=2.659*(-2.156+1.509*p-0.198*p**2+0.011*p**3)+Rv
        if (ala < 0.12):
           ff=ff11+(wl[ik]-1100.)*slope1
        if (ala > 2.2):
           ff=ff99+(wl[ik]-21900.)*slope2

        alambda=av*ff/Rv
        if (alambda < 0.):
            alambda=0.
        
        yob[ik]=fl[ik]*(10**(-0.4*alambda))

    return yob

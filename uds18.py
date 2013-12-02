"""
Custom fits for the lens in UDS-18
"""
import unicorn
import pyfits
import emcee

z_lens = 1.634
z_arc = 2.25

def test():
    
    twod_model = unicorn.reduce.Interlace2D('model_31684.2D.fits')
    twod_orig = unicorn.reduce.Interlace2D('orig_31684.2D.fits')
    twod_line = unicorn.reduce.Interlace2D('acs_31684.2D.fits')
    
    tx, ty_hi = np.loadtxt('SF0_0.emline.hiOIII.txt', unpack=True)
    #tx, ty_lo = np.loadtxt('SF0_0.emline.loOIII.txt', unpack=True)
    
    ty = ty_hi
    
    status = twod_line.compute_model(tx*(1+z_arc), ty*3)
    
    zgrid = np.arange(2.25,2.30,0.001)
    chi2 = zgrid*0.
    
    for i in range(len(zgrid)):
        status = twod_line.compute_model(tx*(1+zgrid[i]), ty*2.5)
        mask = (twod_line.model > 1.e-5) & (twod_line.im['WHT'].data > 0)
        resid = (twod_line.im['SCI'].data-twod_line.im['CONTAM'].data - twod_line.model)
        chi2[i] = np.sum(resid[mask]**2/(twod_line.im['WHT'].data[mask]**2))
        
    z_best = zgrid[chi2 == chi2.min()][0]
    status = twod_line.compute_model(tx*(1+z_best), ty*2.5)
    
    status = twod_model.compute_model()
    
def _objective_full_fit(params, **data):
    """
    Get log(chi2/2) for some set of input parameters
    """
    
    
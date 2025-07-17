def test_placeholder():
    assert 1 + 1 == 2

def test_imports():
    #Test all import statements that occur in this codebase.
    from starrotator import StarRotator
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.special import sph_harm
    import requests
    import shutil
    import urllib.request as request
    from contextlib import closing
    import os.path
    import starrotator.lib.test as test
    import starrotator.lib.operations as ops
    import starrotator.lib.stellar_spectrum
    import starrotator.lib.vgrid as vgrid
    import starrotator.lib.integrate_depr
    import starrotator.lib.util as util
    from starrotator.lib.util import gaussian
    import astropy.io.fits as fits
    import astropy.units as u
    import scipy.interpolate as interp
    import astropy.constants as consts
    import numpy as np
    import copy
    import imp
    import sys
    from importlib.resources import files
    from functools import partial
    import jax
    from jax import jit, lax
    import jax.numpy as jnp
    from starrotator.lib.constants import rad_in_deg, R_sun, d_in_seconds, c as c_light


def test_demo_files():
    from importlib.resources import files
    import starrotator.lib.util as util
    data_path = files("starrotator.data").joinpath("demo_system.txt")
    obs_path = files("starrotator.data").joinpath("demo_observations.txt")
    paramdict = util.read_into_dictionary(data_path)
    util.check_integrity_input(paramdict)



def test_integrators():
    """This tests the consistency of the v1 and v2 integrators."""
    import numpy as np
    from starrotator.lib.util import gaussian
    from starrotator.lib.vgrid import calc_vel_stellar, calc_flux_stellar
    from starrotator.lib.integrate import sum_stellar_spectrum_v1, sum_stellar_spectrum_v2

    Nwl = 5000
    wl = np.linspace(399,401,Nwl)
    fx = gaussian(wl,-0.5,400,0.01,1.0) 

    N = 100
    x = np.linspace(-1,1,N)
    y = x*1.0
    a1,a2 = 0.2,0.2
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0

    flux_disk = calc_flux_stellar(x,y,a1,a2)
    vel_disk  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)

    F1 = sum_stellar_spectrum_v1(wl,fx,x,vel_eq,i_stellar,a1,a2)
    F2 = sum_stellar_spectrum_v2(wl,fx,vel_disk,flux_disk)

    maxdiff = np.max(np.abs((F1-F2)/F1))
    
    assert(maxdiff < 3e-4 * (100/N)**2 )
    #The relative error between these two methods may be more than 
    # 3e-4 for any wavelength point when N = 100. Error goes as the
    # square of N. And I suspect that the 2D integrator is actually
    # the source of the inaccuracy, for low N.






#This breaks until integration, dynamics and new ways of calling these, are refactored. 
#I.e. until a long while from now. But we're getting there!
"""
def test_StarRotator():
    from starrotator import StarRotator
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.special import sph_harm
    import requests
    import shutil
    import urllib.request as request
    from contextlib import closing
    import os.path
    import starrotator.lib.test as test
    import starrotator.lib.operations as ops
    import starrotator.lib.stellar_spectrum
    import starrotator.lib.vgrid as vgrid
    import starrotator.lib.integrate_depr
    import astropy.io.fits as fits
    import astropy.units as u
    import scipy.interpolate as interp
    import astropy.constants as consts
    import numpy as np
    import copy
    import imp
    import sys
    from importlib.resources import files




    KELT2 = StarRotator(586.0,592.0,50.0)
    wl = KELT2.wl
    F_out = KELT2.stellar_spectrum
    spectra = KELT2.spectra
    residuals = KELT2.residual
    KELT2.convolve_spectral_resolution()

    assert KELT2.status == 'success'


    #Test that multiple-convolution is detected and blocked:
    try:
        KELT2.convolve_spectral_resolution()
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Trying convolution twice in a row should be caught.")
    error_trigger=0

    #Testing another grid and another vartype.
    KELT3 = StarRotator(586,592.0,13)
    KELT3.convolve_spectral_resolution()

    assert KELT3.status == 'success'

    #Now test the whole thing with a dictionary as input. First test that the input is well
    #tested:
    in_dict = {'lala':1.0}
    try:
        KELT4 = StarRotator(586,592.0,13,input=in_dict)
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Wrong input dictionary not caught.")
    error_trigger=0


    in_dict = {'veq':114000.0,
    'stelinc':90.0,
    'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
    'u1':0.93,'u2':-0.23,'R':115000.,'mus':0,'model':'PHOENIX','sma_Rs':3.153,
            'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
            'phases':[-0.02,-0.01,0.0,0.01,0.02]}
    KELT5 = StarRotator(586,592.0,13,input=in_dict)

    assert KELT5.status == 'success'

    #
    pysme_error = 0
    try:
        from pysme.sme import SME_Structure as SME_Struct
        from pysme.abund import Abund
        from pysme.synthesize import synthesize_spectrum
        from pysme.solve import solve
        from pysme.linelist.vald import ValdFile
    except:
        print('WARNING: PYSME cannot be imported. PSME functionality will not be available.')
        n_warnings+=1
        pysme_error = 1


    if pysme_error == 0:
        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':0,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02]}#Without setting abund and grid_model keywords.
        try:
            KELT9 = StarRotator(586,592.0,13,input=in_dict)
            error_trigger=1
        except:
            pass
        if error_trigger==1:
            raise Exception("ERROR: Wrong input dictionary not caught.")
        error_trigger=0

        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':5,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02],'grid_model':'atlas12.sav','abund':[],
        'linelist_path':'input/demo_linelist.dat'}
        KELT9 = StarRotator(586,592.0,13,input=in_dict)


        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':6000.0,'FeH':0.3,'logg':4.2,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':5,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02],'grid_model':'marcs2014.sav','abund':[],
        'linelist_path':'input/demo_linelist.dat'}
        KELT9 = StarRotator(586,592.0,13,input=in_dict)

    print('')
    print('')
    print('')
    print('Tests complete.')
    print(f'{n_warnings} warnings triggered.')

    """
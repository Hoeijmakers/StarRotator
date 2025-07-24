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
    import matplotlib.pyplot as plt

    Nwl = 1000
    wl = np.linspace(399,401,Nwl)
    fx = gaussian(wl,-0.5,400,0.01,1.0) 

    N = 200
    x = np.linspace(-1,1,N)
    dx = x[1]-x[0]
    y = x*1.0
    dy = y[1]-y[0]
    a1,a2 = 0.2,0.3
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0

    flux_disk = calc_flux_stellar(x,y,a1,a2,norm=False)
    vel_disk  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)

    F1 = sum_stellar_spectrum_v1(wl,fx,x,vel_eq,i_stellar,a1,a2)*dx*2#Because it integrates only a semicircle.
    F2 = sum_stellar_spectrum_v2(wl,fx,vel_disk,flux_disk)*dx*dy
    maxdiff = np.max(np.abs((F1-F2)/F1))
    
    # print(maxdiff)
    # plt.plot(F1)
    # plt.plot(F2)
    # plt.show()
    assert(maxdiff < 7e-4)
    #The relative error between these two methods may be more than 
    # 7e-4 for any wavelength point when N = 200. Error goes as the
    # square of N. And I suspect that the 2D integrator is actually
    # the source of the inaccuracy, for low N.

    flux_disk = calc_flux_stellar(x,y,a1,a2,norm=True)
    vel_disk  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)

    F1 = sum_stellar_spectrum_v1(wl,fx,x,vel_eq,i_stellar,a1,a2,norm=True)
    F2 = sum_stellar_spectrum_v2(wl,fx,vel_disk,flux_disk)
    maxdiff = np.max(np.abs((F1-F2)/F1))

    assert(maxdiff < 1.2e-4)
    #The relative error shrinks when the fluxes are brute-force normalised.





def test_analytical_limb_darkening():
    from starrotator.lib.operations import vert_int_q_ld, vert_int_q_ld_bounded, circ_int_q_ld
    from starrotator.lib.vgrid import calc_flux_stellar
    import numpy as np
    import matplotlib.pyplot as plt
    N = 1000
    x = np.linspace(-1,1,N)
    dx = x[1]-x[0]
    y = np.linspace(0,1,int(N/2))
    dy = y[1]-y[0]

    a1 = 0.2
    a2 = 0.4

    half_disk = calc_flux_stellar(x,y,a1,a2,norm=False)
    I1 = np.nansum(half_disk,axis=0)*dy
    I2 = vert_int_q_ld(x,a1,a2)
    I3 = vert_int_q_ld_bounded(x,0.0,np.sqrt(1-x**2),a1,a2)
    I4 = circ_int_q_ld(a1,a2)
    # plt.plot(x,I1) # As long as these two do not agree there is some form of error or typo. Need to redo this from scratch or stop caring.
    # plt.plot(x,I2,'--')
    # plt.title('Wow!')
    # plt.show()
    median_error_1 = np.nanmedian(np.abs(I2-I1)/I2)
    median_error_2 = np.nanmedian(np.abs(I3-I2)/I3)
    median_error_3 = np.nanmedian(np.abs(I4-np.nansum(I2)*dx*2)/I4)
    # print(median_error_1)
    # print(median_error_2)
    # print(median_error_3)

    # The error between the analytical column-wise vertical integral and the numerical
    # column-wise integral:
    assert(median_error_1 < 2e-3)

    # The error between the analytical column-wise vertical integral from y=0 to y=disk-edge
    # and the analytical column-wise vertical integral with the indefinite integral explicitly
    # filled in:
    assert(median_error_2 < 1e-9)

    # The error between the column-wise vertial integration over a half disk times two, 
    # and the analytical disk-integral:
    assert(median_error_3 < 3e-5)



    # Test that if the edges are avoided by a minimal amount, the indefinite integral version
    # of the code does not hit singularities.
    x = np.linspace(-1,1,200)*0.99999
    I = vert_int_q_ld_bounded(x,0,np.sqrt(1-x**2)*0.99999,a1,a2)
    assert(np.sum(np.isnan(I))==0)

    # And that if the edges of x=-1 and x=1 are hit, there are precisely 2 nans.
    x = np.linspace(-1,1,200)
    I = vert_int_q_ld_bounded(x,0,np.sqrt(1-x**2)*0.99999,a1,a2)
    assert(np.sum(np.isnan(I))==2)

    # Test that the analytical full-circle integration is normalised the same way as the 
    # edge-bounded indefinite-integral version.
    X = np.linspace(-1,1,800)*0.99999
    dX = X[1]-X[0]
    I1 = circ_int_q_ld(0.0,0.0)
    I2 = np.sum(vert_int_q_ld_bounded(X,-np.sqrt(1-X**2)*0.99999,np.sqrt(1-X**2)*0.99999,0,0)*dX)
    assert(np.abs(I1-I2)/I1 < 1e-4)




def test_analytical_integration():
    from starrotator.lib.operations import circ_int_q_ld
    import numpy as np
    from starrotator.lib.vgrid import calc_flux_stellar
    # This tests that the numerical integration of the limb darkened disk approaches
    # the analytical result.
    N = 100
    x = np.linspace(-1,1,N)
    y = x*1.0
    a1,a2 = 0.2,0.3
    flux_disk_2 = calc_flux_stellar(x,y,a1,a2,norm=False)
    dx = x[1]-x[0]
    dy = dx*1.0
    int_numerical = np.nansum(flux_disk_2*dx*dy)
    int_theoretical = circ_int_q_ld(a1,a2)
    numerical_error = (int_numerical-int_theoretical) / int_theoretical
    assert(np.abs(numerical_error) < 0.003)
    # For N=100, the relative error on the total integral is 3e-3 and 
    # # this depends a bit on the limb darkening parameters.



def test_hidden_flux():
    import numpy as np
    from starrotator.lib.integrate import create_hidden_grid_array, sum_hidden_spectrum_v2, sum_stellar_spectrum_v1, sum_hidden_spectrum_v1
    from starrotator.lib.vgrid import calc_vel_stellar
    from starrotator.lib.vgrid import calc_flux_stellar
    import jax.numpy as jnp
    import matplotlib.pyplot as plt
    #This tests that a planet with a radius of 0.3 Rstar transiting a star without limb darkening
    #creates a transit depth of 0.09.
    N = 400
    x = jnp.linspace(-1,1,N)
    y = x*1.0
    dx = x[1]-x[0]

    xp = jnp.array([-0.3,-0.29,-0.28-0.27])
    yp = xp*0.5


    a1,a2 = 0.0,0.0
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0
    Rp = 0.3
    wl = np.linspace(399,401,1000)
    fx = np.ones_like(wl)

    flux_disk_total = calc_flux_stellar(jnp.array(x),jnp.array(y),a1,a2,norm=False) *dx**2
    flux_array,vel_array,dxR = create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=200)
    F_out = sum_hidden_spectrum_v2(wl,fx,vel_array,flux_array,batched=True) *dxR**2 / jnp.nansum(flux_disk_total)
    mean_error_per_wl = np.mean(np.abs(F_out-0.09))   
    assert(mean_error_per_wl < 1e-4)


    #Final testing: 
    xp = jnp.array([-0.2,0.0,0.2])
    yp = xp*0.0

    F_in_v1 = sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=500)

    x = jnp.linspace(-1,1,500)
    y = x*1.0
    dx = x[1]-x[0]

    flux_grid_array,vel_grid_array,dxR = create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=300)
    F_in_v2 = sum_hidden_spectrum_v2(wl,fx,vel_grid_array,flux_grid_array,batched=True) *dxR**2 

    mean_error = np.mean((F_in_v1[0]-F_in_v2[0])/F_in_v1[0])
    assert(mean_error < 6e-4)
    # plt.plot(wl,(F_in_v1[0]-F_in_v2[0])/F_in_v1[0])
    # plt.show()




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
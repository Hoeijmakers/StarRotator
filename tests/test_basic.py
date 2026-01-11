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
    # import starrotator.lib.test as test
    import starrotator.lib.operations as ops
    import starrotator.lib.stellar_spectrum
    import starrotator.lib.vgrid as vgrid
    import starrotator.lib.integrate_depr
    import starrotator.lib.util as util
    import starrotator.lib.dynamics as dynamics
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
    from starrotator.lib.operations import circ_int_q_ld
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

    # flux_disk = calc_flux_stellar(x,y,a1,a2,norm=False)
    # vel_disk  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)

    F1 = sum_stellar_spectrum_v1(wl,fx,vel_eq,i_stellar,a1,a2,N=N)
    F2 = sum_stellar_spectrum_v2(wl,fx,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=N,batched=True)
    maxdiff = np.max(np.abs((F1-F2)/F1))
    # print(circ_int_q_ld(a1,a2))
    # print(maxdiff)
    # plt.plot(F1)
    # plt.plot(F2)
    # plt.show()
    assert(maxdiff < 7e-4)
    #The relative error between these two methods may not be more than 
    # 7e-4 for any wavelength point when N = 200. Error goes as the
    # square of N. And I suspect that the 2D sum is actually
    # the source of the inaccuracy, for low N.






def test_analytical_limb_darkening():
    from starrotator.lib.operations import vert_int_q_ld, vert_int_q_ld_bounded, circ_int_q_ld
    from starrotator.lib.vgrid import calc_flux_stellar
    import numpy as np
    # import matplotlib.pyplot as plt
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
    # The passing of this test indicates that the analytical computation is correct,
    # and should be used - no need to do numerical integration in the v1 case.


def test_mu_integration_v1():
    """This tests that integrating the whole disk using v1 and v1_mu produces the same output for 
    zero limb darkening on v1 and equal spectra on v1_mu (fake mu-dependence)."""
    from starrotator.lib.integrate import sum_stellar_spectrum_v1_mu, sum_stellar_spectrum_v1
    import numpy as np
    import jax.numpy as jnp
    from starrotator.lib.util import gaussian
    import matplotlib.pyplot as plt

    wl = np.linspace(500,505,1000)
    mu = jnp.linspace(0,1,20)

    fx = wl*0.0+1.0
    fx = gaussian(wl,-0.5,jnp.mean(wl),0.02,cont=1.0)
    fx_array = wl[None,:]*0.0+ mu[:,None]*0.0+fx

    vsini= 100
    fx_sum = sum_stellar_spectrum_v1_mu(wl,fx_array,vsini,90.0,mu,N=201)
    fx_v1 = sum_stellar_spectrum_v1(wl,fx,vsini,90.0,0.0,0.0,N=201)

    # plt.plot(fx_sum)
    # plt.plot(fx_v1)
    # plt.show()

    max_rel_error = np.max(np.abs((fx_sum-fx_v1)/fx_v1))
    assert max_rel_error < 1e-6



def test_hidden_flux():
    import numpy as np
    from starrotator.lib.integrate import sum_hidden_spectrum_v2, sum_hidden_spectrum_v1
    from starrotator.lib.vgrid import calc_vel_stellar
    from starrotator.lib.vgrid import calc_flux_stellar
    import jax.numpy as jnp
    import matplotlib.pyplot as plt
    #This tests that a planet with a radius of 0.3 Rstar transiting a star without limb darkening
    #creates a transit depth of 0.09.

    xp = jnp.array([-0.3,-0.29,-0.28,-0.27])
    yp = xp*0.5


    a1,a2 = 0.0,0.0
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0
    Rp = 0.3
    wl = np.linspace(399,401,1000)
    fx = np.ones_like(wl)


    F_out = sum_hidden_spectrum_v2(wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=200,batched=True)



    mean_error_per_wl = np.mean(np.abs(F_out-0.09))   
    assert(mean_error_per_wl < 1.5e-4)


    #Final testing to ensure that v1 and v2 produce essentially the same output (and that means a transit depth of 0.9): 
    xp = jnp.array([-0.2,0.0,0.2])
    yp = xp*0.0
    F_in_v1 = sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=500)
    F_in_v2 = sum_hidden_spectrum_v2(wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=300,batched=True)

    mean_rel_error = np.mean((F_in_v1[0]-F_in_v2[0])/F_in_v1[0])

    # print(mean_rel_error)
    # plt.plot(wl,F_in_v1[0])
    # plt.plot(wl,F_in_v2[0])
    # plt.show()

    assert(mean_rel_error < 8e-4)


def test_mu_interpolation():
    """This tests that the interpolation of a spectrum between two mu angles works as intended."""
    import numpy as np
    from starrotator.lib.integrate import interp_fx_array
    from starrotator.lib.vgrid import calc_vel_stellar
    from starrotator.lib.vgrid import calc_flux_stellar
    import jax.numpy as jnp
    import matplotlib.pyplot as plt
    from starrotator.lib.util import gaussian
    import sys
    wl = np.linspace(500,505,1001)
    mu = jnp.array([0.0,0.5,1.0]) # Choose 3 mu angles, 0, half and full.
    fx = wl*0.0+1.0
    fx = gaussian(wl,-0.5,jnp.mean(wl),0.07,cont=0.0)


    fx_array = wl[None,:]*0.0+ mu[:,None]*fx + 1.0 # Scale the linedepth by the mu angle. So the linedepth increases linearly.

    mu_to_test = jnp.array([0.0,0.25,0.5,0.75,1.0,1.2])

    fx_interpolated = interp_fx_array(mu,fx_array,mu_to_test) #Interpolate at a range of values:
    #0 should return a line with depth 0.0
    #0.25 should return a line with depth 25% of that of fx, etc.
    #1.2 maxes out, so it should return a line depth of 100% of fx. Not extrapolate or fail.

    linedepth_in = np.min(fx)
    linedepths = np.min(fx_interpolated-1.0,axis=1)

    for i in range(len(mu_to_test)):
        depth_error = linedepths[i] - linedepth_in * mu_to_test[i]
        if 0.0<=mu_to_test[i]<=1.0:
            assert(depth_error < 1e-6)
        if mu_to_test[i]>1:
            assert(depth_error-(mu_to_test[i]-1.0) < 1e-6)




def test_hidden_flux_mu_v1():
    """This tests that the flux of the obscured part of the disk works in the case of a mu-dependent spectrum (v1_mu),
    by comparing to the mu-independent method (v1) in case that the input spectrum is equal for all mu (fake mu dependence)."""
    import numpy as np
    from starrotator.lib.integrate import create_hidden_grid_array, sum_hidden_spectrum_v2, sum_stellar_spectrum_v1, sum_hidden_spectrum_v1, sum_hidden_spectrum_v1_mu
    from starrotator.lib.vgrid import calc_vel_stellar
    from starrotator.lib.vgrid import calc_flux_stellar
    import jax.numpy as jnp
    import matplotlib.pyplot as plt
    from starrotator.lib.util import gaussian
    import sys
    import pdb
    # This tests that the mu-dependent hidden-flux works, and returns the same answer as mu-independent, if all spectra are the same.
    # It also tests that the transit depth of a planet with R=0.1 is 1%.
    wl = np.linspace(500,505,1000)
    mu = jnp.linspace(0,1,20)

    fx = wl*0.0+1.0
    fx = gaussian(wl,-0.5,jnp.mean(wl),0.02,cont=1.0)
    fx_array = wl[None,:]*0.0+ mu[:,None]*0.0+fx #All spectra are the same.


    a1,a2 = 0.0,0.0
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0
    Rp = 0.1

    N = 300


    xp = jnp.array([-0.7,-0.3,-0.29,-0.28,-0.27,0.0,0.1,0.15,0.2,0.7])
    yp = xp*0.0 + 0.2

    F_in_v1 = sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=N)
    F_in_v2 = sum_hidden_spectrum_v1_mu(wl,fx_array,xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)


    for i in range(len(xp)):
        mean_error = np.mean(np.abs(F_in_v1[i]-F_in_v2[i])/F_in_v1[i])
        assert(mean_error < 1e-8)




    # Testing a more complicated mu-dependence:

    xp = jnp.array([-0.9,-0.45,0.0,0.45,0.9])
    yp = xp*0.0

 

    fx_array_2 = []
    for i in range(len(mu)):
        fx_array_2.append(gaussian(wl,-0.5*(1-mu[i]),jnp.mean(wl),0.02,cont=1.0+mu[i]))
    fx_array_2=np.array(fx_array_2)


    # plt.figure()
    # for i in range(len(fx_array_2)):
    #     plt.plot(wl,fx_array_2[i])
    # plt.title('Input spectrum arrays.')
    # plt.show()

    F_in_v3 = sum_hidden_spectrum_v1_mu(wl,jnp.array(fx_array_2),xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)


    linedepths = F_in_v3[:,0] - jnp.min(F_in_v3,axis=1)
    contlevel = F_in_v3[:,0]

    assert(linedepths[0]>linedepths[1])
    assert(linedepths[2] < 1e-7) # The central position is at disk center, meaning mu=1.0, meaning line amplitude zero.
    assert(np.abs(linedepths[0]-linedepths[-1])<1e-5) # Symmetry.
    assert(contlevel[2]/contlevel[0] - 2/(1+jnp.sqrt(1-xp[0]**2)) < 1e-4)
    # plt.figure()
    # for i in range(len(F_in_v4)):
    #     plt.plot(wl,F_in_v4[i])
    # plt.title('Hidden spectrum arrays.')
    # plt.show()

    # Testing that xp and yp are interchangeable as far as mu is concerned:
    xp = jnp.array([-0.9,-0.9/jnp.sqrt(2),0.0,0.0,0.9/jnp.sqrt(2),0.9])
    yp = jnp.array([ 0.0,0.9/jnp.sqrt(2),0.9,-0.9,0.9/jnp.sqrt(2),0.0])

    # rp = jnp.sqrt(xp**2+yp**2)
    # print(rp)
    # mup = jnp.sqrt(1-rp**2)
    # print(mup)
    F_in_v4 = sum_hidden_spectrum_v1_mu(wl,jnp.array(fx_array_2),xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)
    linedepths = F_in_v4[:,0] - jnp.min(F_in_v4,axis=1)
    contlevel = F_in_v4[:,0]

    # print(linedepths)
    # print(contlevel)
    assert(np.mean(np.abs(linedepths-linedepths[0])/linedepths[0]) < 1e-2)
    assert(np.mean(np.abs(contlevel-contlevel[0])/contlevel[0]) < 1e-2)
    # plt.figure()
    # for i in range(len(F_in_v4)):
    #     plt.plot(wl,F_in_v4[i])
    # plt.title('Hidden spectrum arrays.')
    # plt.show()


    # Testing that travelling over the disk changes the velocity, or not:
    xp = jnp.linspace(-1,1,10)
    yp = xp*0.0
    F_in_v5 = sum_hidden_spectrum_v1_mu(wl,jnp.array(fx_array_2),xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)
    linepos = jnp.argmin(F_in_v5,axis=1)
    delta_index = linepos[1:]-linepos[0:-1]
    assert(jnp.min(delta_index)>0) #Make sure that the line always moves to longer wavelengths.

    # When transiting vertically:
    yp = jnp.linspace(-1,1,10)
    xp = yp*0.0
    F_in_v5 = sum_hidden_spectrum_v1_mu(wl,jnp.array(fx_array_2),xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)
    linepos = jnp.argmin(F_in_v5,axis=1)
    delta_index = linepos[1:]-linepos[0:-1]
    assert(jnp.max(jnp.abs(delta_index))==0) #Make sure that the line is constant when moving parallel to the star's axis.


    N=1000
    xp = jnp.array([-1,0.0])
    yp = jnp.array([-1,0.0])
    F_in_v6 = sum_hidden_spectrum_v1_mu(wl,jnp.array(fx_array_2),xp,yp,Rp,vel_eq,i_stellar,mu,N=N,small_planet=True)
    assert(np.max(F_in_v6[0]) < 1e-9) # Outside of the disk, we are 0.0
    # print(Rp**2)
    # print(np.mean(np.abs(F_in_v6[1]-(2*Rp**2))/(2*Rp**2)))
    assert(np.mean(np.abs(F_in_v6[1]-(2*Rp**2))/(2*Rp**2)) <  5e-5) # At disk center, we are 2x Rp**2 in transit depth.




def test_cache():
    """This tests the handling of cache configuration and creation."""
    import starrotator.lib.util as ut
    import copy
    from pathlib import Path
    # Define a temporary testing config file to not overwrite existing configs.
    current_configfile = copy.deepcopy(ut.CONFIG_FILE)
    ut.CONFIG_FILE = ut.CONFIG_DIR / "config_testing.json"


    # Test that a default cache path can be defined and that nothing exists there.
    P = ut.get_default_cache_dir()
    # assert(P.exists() == False) # This would fail if starrotator is not freshly installed.

    # Test that load_config (which tries to read a config file that does not exist)
    # returns the default cache path.
    config = ut.load_config()
    assert("cache_dir" in config)
    assert(config["cache_dir"] == str(P))

    # Test that get_cache_dir does the same thing, in absence of a config file it returns
    # the default cache dir and makes it.
    P2 = ut.get_cache_dir()
    assert(str(P) == str(P2))
    assert(P.exists()==True)


    # Thest that a dummy config file can be created. Permission errors are caught and should not make the test fail.
    ut.save_config({'cache_dir' : 'void'}) 
    config = ut.load_config()
    assert("cache_dir" in config)
    assert(config["cache_dir"] == 'void')

    ut.save_default_config()
    config = ut.load_config()
    assert("cache_dir" in config)
    assert(config["cache_dir"] == str(P))

    P3 = ut.get_cache_dir()
    assert(str(P3) == str(P))


    # Test that set_cache_dir works to update the config file.
    ut.set_cache_dir('void')
    config = ut.load_config()
    assert("cache_dir" in config)
    assert(config["cache_dir"] == 'void')

    # Remove garbage.
    Path(ut.CONFIG_FILE).unlink()
    ut.CONFIG_FILE = copy.deepcopy(current_configfile)





def test_orbit():
    """This tests Keplerian orbits computed using Jaxoplanet."""
    import numpy as np
    from starrotator.lib.dynamics import orbit_euclidian
    phase = np.linspace(-0.5,0.5,1000)
    x,y,z = orbit_euclidian(phase, a = 5.0, m = 0.0, P = 4.0, e = 0.5, omega = 0.0,i=90.0)
    num_error = np.mean(np.abs(y))
    assert num_error < 1e-6
    assert np.abs(np.max(x-7.5)) < 1e-4


def test_fast_doppler():
    from starrotator.lib.integrate import sum_stellar_spectrum_v1, sum_stellar_spectrum_v2, sum_stellar_spectrum_v1_mu, sum_hidden_spectrum_v2, sum_hidden_spectrum_v1, sum_hidden_spectrum_v1_mu
    import timeit
    import numpy as np
    import jax.numpy as jnp
    from starrotator.lib.util import gaussian
    from starrotator.lib.dynamics import doppler_shift,doppler_factor,doppler_shift_dlogl
    from starrotator.lib.operations import vert_int_q_ld, circ_int_q_ld, constant_velocity_grid

    from jax import jit
    import matplotlib.pyplot as plt    
    xp = jnp.array([-0.3,-0.29,-0.28,-0.27,0.0,0.1])
    yp = xp*0.5
    a1,a2 = 0.0,0.0
    i_stellar,vel_eq,diff_rot_rate = 90.0,100.0,0.0
    Rp = 0.1
    wlmin,wlmax,Nwl = 463,467,10000
    wl = np.linspace(wlmin,wlmax,Nwl)
    mu = jnp.linspace(0,1,10)
    N = 300
    fx = gaussian(wl,-0.5,jnp.mean(wl),0.02,cont=1.0)
    #Define a constant velocity wavelength grid as well of the same size and limits as the linear one:
    wlc = np.mean(np.array([wlmin,wlmax]))
    wl2,dlogl = constant_velocity_grid(wlmin,wlmax,Nwl)
    fx2 = gaussian(wl2,-0.5,wlc,0.02,cont=1.0)
    fx2_array = wl[None,:]*0.0+ mu[:,None]*0.0+fx2


    x = jnp.linspace(-1,1,N)
    v_axis = x*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.

    fx_shifted_1 = doppler_shift(wl2,fx2,v_axis)
    fx_shifted_2 = doppler_shift_dlogl(dlogl,fx2,v_axis)
    mean_difference = np.mean(np.abs((fx_shifted_2[0]-fx_shifted_1[0])/fx_shifted_2[0])) # mean relative difference.
    assert(mean_difference < 3e-5)
    mean_difference = np.mean(np.abs((fx_shifted_2[-1]-fx_shifted_1[-1])/fx_shifted_2[-1]))
    assert(mean_difference < 3e-5)


    integrated_1 = sum_stellar_spectrum_v1(wl2,fx2,vel_eq,i_stellar,a1,a2,N=400,constant_dlogl=False)
    integrated_2 = sum_stellar_spectrum_v1(dlogl,fx2,vel_eq,i_stellar,a1,a2,N=400,constant_dlogl=True)
    mean_difference = np.mean(np.abs((integrated_1 - integrated_2)/integrated_1))

    assert(mean_difference < 1e-5)


    #Test 2D input to dopplershift functions.
    vaxis_2d = np.vstack([v_axis,v_axis])
    fx_shifted_1 = doppler_shift(wl2,fx2,vaxis_2d)
    fx_shifted_2 = doppler_shift_dlogl(dlogl,fx2,vaxis_2d)
    assert np.shape(fx_shifted_1)[0] == np.shape(fx_shifted_2)[0]
    assert np.shape(fx_shifted_1)[1] == np.shape(fx_shifted_2)[1]
    assert np.shape(fx_shifted_1)[2] == np.shape(fx_shifted_2)[2]
    assert np.shape(fx_shifted_1)[0] == len(v_axis)
    assert np.shape(fx_shifted_1)[1] == 2
    assert np.shape(fx_shifted_1)[2] == len(fx2)
    mean_difference = np.mean(np.abs((fx_shifted_1[0,0,:]-fx_shifted_2[0,0,:])/fx_shifted_1[0,0,:]))
    assert(mean_difference < 3e-5)


    integrated_1 = sum_stellar_spectrum_v2(wl2,fx2,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=200,batched=True,constant_dlogl=False)
    integrated_2 = sum_stellar_spectrum_v2(dlogl,fx2,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=200,batched=True,constant_dlogl=True)
    mean_difference = np.mean(np.abs((integrated_1 - integrated_2)/integrated_1))
    assert(mean_difference < 1e-5)


    integrated_1 = sum_hidden_spectrum_v1(wl2,fx2,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100,constant_dlogl=False)
    integrated_2 = sum_hidden_spectrum_v1(dlogl,fx2,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100,constant_dlogl=True)
    mean_difference = np.mean(np.abs((integrated_1 - integrated_2)/integrated_1))
    assert(mean_difference < 4e-5)


    integrated_1 = sum_hidden_spectrum_v1_mu(wl2,fx2_array,xp,yp,Rp,vel_eq,i_stellar,mu,N=100,constant_dlogl=False)
    integrated_2 = sum_hidden_spectrum_v1_mu(dlogl,fx2_array,xp,yp,Rp,vel_eq,i_stellar,mu,N=100,constant_dlogl=True)
    mean_difference = np.mean(np.abs((integrated_1 - integrated_2)/integrated_1))
    assert(mean_difference < 4e-5)


    integrated_1 = sum_hidden_spectrum_v2(wl2,fx2,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=100,batched=True,constant_dlogl=False)
    integrated_2 = sum_hidden_spectrum_v2(dlogl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=100,batched=True,constant_dlogl=True)
    mean_difference = np.mean(np.abs((integrated_1 - integrated_2)/integrated_1))
    assert(mean_difference < 2e-3) # This is a large relative error. Why?






def test_default_computation():
    from starrotator import StarRotator
    KELT9 = StarRotator(500,502,100)
    assert KELT9.status == 'success computing spectra'
    KELT9.drr = 0.1
    KELT9.status = ''
    KELT9.compute_spectrum()
    assert KELT9.status == 'success computing spectra'



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
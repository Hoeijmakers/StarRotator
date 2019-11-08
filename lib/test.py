#######################
#all testing functions
#
#
#######################



def test_parser(args):
    if args.wave_end <= args.wave_start:
        raise ValueError("End of wavelength range < start of wavelength range.")
    if args.velStar > 1000.:
        raise ValueError("Please give the equatorial velocity in meters per second.")
    if not -90. <= args.stelinc <= 90.:
        raise ValueError("The stellar inclination has to be [-90.,90.].")
    if not 0. <= args.orbinc <= 90.:
        raise ValueError("The orbital inclination has to be [0.,90.].")
    if not 5 <= args.grid_size <= 1000:
        raise ValueError("Grid size has to be [5,1000].")


def test_library():
    import lib.operations
    import lib.stellar_spectrum

def test_ld():
    import lib.operations as ops

    w1 = ops.limb_darkening(0,0,0)
    w2 = ops.limb_darkening(0,1,0)
    w3 = ops.limb_darkening(0,0,1)

    if w1 != 1 or w2 !=0 or w3 != 0:
        print('ERROR: limb darkening does not pass.')


def test_integrators():
    """This tests the consistency of the fast and the slow surface integrators."""
    import numpy as np
    import lib.vgrid as vgrid
    import lib.plotting as plt
    import lib.operations as ops
    import lib.stellar_spectrum as spectrum
    import lib.integrate as integrate
    import time
    import matplotlib.pyplot as pl
    import sys
    import pdb
    T = 10000.0
    logg = 4.5
    u1 = 0.387
    u2 = 0.178
    wlmin = 588.2
    wlmax = 590.0
    #two arrays for the x and y axis
    x = np.linspace(-1,1,num=2*25) #in units of stellar radius
    y = np.linspace(-1,1,num=2*25) #in units of stellar radius

    #calculate the velocity grid
    vel_grid = vgrid.calc_vel_stellar(x,y,90.0,100000.0,0.0,0.0)
    flux_grid = vgrid.calc_flux_stellar(x,y,u1,u2)
    wl,fx = spectrum.read_spectrum(T,logg)
    wlF,F = integrate.build_spectrum_fast(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid)


    #Test the building of the local spectrum in-between.
    wlFp,Fp,fluxp,mask = integrate.build_local_spectrum_fast(0.4,0.0,0.15,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid)
    wlFp2,Fp2,fluxp2,mask2 = integrate.build_local_spectrum_slow(0.4,0.0,0.15,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid)
    epsilon = np.abs(np.nansum(Fp-Fp2)/np.nansum(Fp))#Test that they are consistent within error.
    if epsilon > 1e-8:
        raise Exception("Error: The integrated difference of the local spectra built with the slow and fast integrators differs by more than 1e-8 times the integrated spectrum.")

    #
    # pdb.set_trace()
    # pl.plot(wlF,F-np.nanmax(F))
    # pl.plot(wlF,(F-Fp)-np.nanmax(F-Fp))
    # pl.show()

    #Test the slow star integrator and compare.
    wlF2,F2 = integrate.build_spectrum_slow(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid)
    epsilon = np.abs(np.nansum(F-F2)/np.nansum(F))
    if epsilon > 1e-8:
        raise Exception("Error: The integrated difference of the disk-integrated spectra built with the slow and fast integrators differs by more than 1e-8 times the integrated spectrum.")
    # pl.plot(wl,fx)
    # pl.plot(wlF,F)
    # pl.plot(wlF2,F2)
    # pl.xlim((550.0,560.0))
    # pl.show()


def test_plotting():
    import numpy as np
    import imp
    #We create a fake velocity map with a hole in it.
    x = np.arange(-200,200,1)
    y = np.arange(-200,200,1)
    z = np.zeros((len(x),len(y)))
    xp = 50
    yp = 20#Fiducial planet coordinates.
    for i in x:
        for j in y:
            if np.sqrt(x[i]**2 + y[j]**2) <= 190:
                z[j,i] = x[i]/np.e
            else:
                z[j,i] = np.nan
            if np.sqrt((x[i]-xp)**2 + (y[j]-yp)**2) <= 10:
                z[j,i] = np.nan

    import lib.plotting
    imp.reload(lib.plotting)
    lib.plotting.plot_star_2D(x,y,z,quantities=('x','y','v'),units=('','','km/s'),noshow=True)

def test_plot_3D():
    import lib.plotting
    lib.plotting.plot_star_3D()

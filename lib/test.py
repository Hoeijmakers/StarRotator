#######################
#all testing functions
#
#
#######################


def test_parser(args):
    if args.wave_end <= args.wave_start:
        raise ValueError("End of wavelength range < start of wavelength range.")
    if args.velStar < 500.:
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

def investigate_SPECTRUM():
    """

    PHOENIX UNIT:  'erg/s/cm2/cm'
    SPECTRUM UNIT: 'erg/s/cm2/A'
    """
    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    import lib.vgrid as vgrid
    import lib.stellar_spectrum as spectrum
    import lib.integrate as integrate
    import pdb
    import sys
    wave_start = 580.0
    wave_end = 600.0
    velStar = 110000.0
    stelinc = 90.0
    drr = 0.0
    pob = 0.0
    grid_size = 100
    T = 5000.0
    logg = 4.5
    Z = 0.0
    # u1 = 0.387
    # u2 = 0.178
    u1 = 0.93
    u2 = -0.23
    mus = 1.0#Normal operation without CLV.
    # mus = np.linspace(0.0,1.0,3)#Uncomment this to run in CLV mode with SPECTRUM.


    test_KURUCZ()
    wl,fx_list = spectrum.compute_spectrum(T,logg,Z,mus,wave_start,wave_end,mode='anf')
    wl3,fx_list3 = spectrum.compute_spectrum(T,logg,Z,1.0,wave_start,wave_end,mode='anM')
    wl*=u.nm
    fx_list*=u.erg/u.cm/u.cm/u.AA/u.s
    fx_list3*=u.erg/u.cm/u.cm/u.AA/u.s
    print(fx_list)
    wl2,fx2 = spectrum.read_spectrum(T,logg,metallicity=Z)
    wl2*=u.nm
    fx2*=u.erg/u.cm/u.cm/u.cm/u.s
    plt.plot(wl,fx_list.to(u.erg/u.cm/u.cm/u.cm/u.s),label='SPECTRUM (disk-integrated)')
    plt.plot(wl2,fx2.to(u.erg/u.cm/u.cm/u.cm/u.s),label='PHOENIX (disk-integrated)')
    plt.plot(wl3,fx_list3.to(u.erg/u.cm/u.cm/u.cm/u.s),label='SPECTRUM (center of disk)')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absolute flux (erg/cm^2/cm/s)')
    plt.xlim(wave_start,wave_end)
    plt.legend()
    plt.show()
    sys.exit()






    #Two arrays for the x and y axes
    x = np.linspace(-1,1,num=2*grid_size) #in units of stellar radius
    y = np.linspace(-1,1,num=2*grid_size) #in units of stellar radius
    #Calculate the velocity and flux grids
    vel_grid = vgrid.calc_vel_stellar(x,y,stelinc,velStar,drr, pob)
    flux_grid = vgrid.calc_flux_stellar(x,y,u1,u2)
    # for fx in fx_list:
    #     plt.plot(wl,fx)
    # plt.show()
    # pdb.set_trace()
    wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,mus,wave_start,wave_end,x,y,vel_grid)
    wlF2,F2 = integrate.build_spectrum_fast(wl2,fx2,wave_start,wave_end,x,y,vel_grid,flux_grid)
    wlp,Fp,flux,mask = integrate.build_local_spectrum_limb_resolved(-0.3,0.0,0.1,wl,fx_list,mus,wave_start,wave_end,x,y,vel_grid)
    wlp2,Fp2,flux2,mask2 = integrate.build_local_spectrum_fast(-0.3,0.0,0.1,wl2,fx2,wave_start,wave_end,x,y,vel_grid,flux_grid)
    #This overplots non-rotating SPECTRUM and PHOENIX spectra, normalised.
    # plt.plot(wl,fx_list[-1]/max(fx_list[-1]),label='SPECTRUM')
    # plt.plot(wl2,fx2/5e15,label='PHOENIX')
    # plt.xlabel('Wavelength (nm)')
    # plt.ylabel('Max-normalised flux')
    # plt.title('T = %s K, log(g) = %s' % (T,logg))
    # plt.legend()
    # plt.show()
    # pdb.set_trace()
    plt.plot(wlF,F/max(F),color='skyblue',alpha=0.5)
    plt.plot(wlF,(F-Fp)/np.nanmax(F-Fp),color='skyblue',label='SPECTRUM')
    plt.plot(wlF2,F2/max(F2),color='red',alpha=0.5)
    plt.plot(wlF2,(F2-Fp2)/np.nanmax(F2-Fp2),color='red',label='PHOENIX')
    plt.legend()
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Max-normalised flux')
    plt.title('T = %s K, log(g) = %s, vsini = 110km/s' % (T,logg))
    plt.show()



def test_KURUCZ():
    """This tests that the Kurucz models used to run SPECTRUM are in place,
    and asks the user to attempt to automatically download them, otherwise.
    This test is used to make sure that compute_spectrum() has the necessary files
    to run.

    After running all of this, all the kurucz models should be in place, and
    compute_spectrum() should be able to call them safely, as long as it is pointed
    to the correct kurucz_root folder ('data/KURUCZ' by default)."""
    import os
    from os import path
    import lib.stellar_spectrum
    import lib.test as test
    import sys
    kurucz_root = 'data/KURUCZ/'#THIS IS THE FIRST PLACE WHERE THIS ROOT FOLDER IS DEFINED. THE ONLY OTHER IS IN STELLAR_SPECTRUM.PY.
    #THESE TWO SHOULD MATCH.

    if path.isdir(kurucz_root) != True:
        os.mkdir(kurucz_root)
        print('WARNING: Kurucz model folder does not exist at '+kurucz_root+'.')
        print('Proceed with automatic download from http://kurucz.harvard.edu?')
        val = input("Proceed with automatic download from http://kurucz.harvard.edu? (y/n) ")
        if val == 'Y' or val == 'y':
            lib.stellar_spectrum.download_kurucz(kurucz_root)
        else:
            print('Aborting')
            sys.exit()

    #This hardcodes the models that should be in the kurucz supermodels.
    #If the models are present, but these are wrong, compute_spectrum() will raise
    # an error.
    T_a,logg_a,Z_a = lib.stellar_spectrum.available_kurucz_models()
    #We now need to construct the filename of the supermodel, using the provided metallicity.
    for Z in Z_a:
        mpath = lib.stellar_spectrum.construct_kurucz_path(Z,kurucz_root)
        if path.isfile(mpath) != True:
            print('WARNING: Kurucz model '+mpath+' does not exist.')
            val = input("Proceed with automatic download from http://kurucz.harvard.edu? (y/n) ")
            if val == 'Y' or val == 'y':
                lib.stellar_spectrum.download_kurucz(kurucz_root)
            else:
                print('Aborting')
                sys.exit()



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






#I copied the following functions from my own project that I use to test variables.
def nantest(var,varname=''):
    import numpy as np
    if np.isnan(var).any()  == True:
        raise ValueError("Variable %s contains NaNs." % varname)
    if np.isinf(var).any()  == True:
        raise ValueError("Variable %s contains in-finite values." % varname)

def postest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) <= 0:
        raise ValueError('Variable %s is not strictly positive' % varname)

def notnegativetest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) < 0:
        raise ValueError('Variable %s is negative.' % varname)

def file_exists(file,varname=''):
    """This program tests if a file exists, and prints an error otherwise."""
    from os import path
    import sys
    typetest(file,str,varname='file in test.file_exists')
    if path.isfile(file) != True:
        print('Error: File %s does not exist.' % varname)
        sys.exit()

def dir_exists(dir,varname=''):
    """This program tests if a directory exists, and prints an error otherwise."""
    from os import path
    import sys
    typetest(dir,str,varname='dir in test.dir_exists')
    if path.isdir(dir) != True:
        print('Error: Directory %s does not exist.' % varname)
        sys.exit()

def typetest(var,vartype,varname=''):
    """This program tests the type of var which has the name varname against
    the type vartype (or the types in list vartype), and raises an exception if
    either varname is not a string or if type(var) is not equal to (any element in)
    vartype.

    Example:
    a = 'ohai'
    utils.typtest('a',a,str)"""
    if isinstance(varname,str) != True:
        print(varname)
        raise Exception("Input error in typetest: varname should be of type string.")

    #Test two cases, if vartype is a list, we test each element.
    if isinstance(vartype,list) == True:
        error = 1
        msg = "Type error: Variable %s should be one of these types: " % varname
        for type in vartype:
            msg+=' %s' % type
            if isinstance(var,type) == True:
                error = 0
        if error == 1:
            raise Exception(msg)
    else:
        #Otherwise, we do the classical typetest.
        if isinstance(var,vartype) != True:
            raise Exception("Type error: Variable %s should be  of type %s." % (varname,vartype))

def typetest_array(var,vartype,varname=''):
    """This program tests the type of the elements in the array or list var which has the
    name varname, against the type vartype, and raises an exception if either
    varname is not a string, type(var) is not equal to numpy.array or list, or the elements of
    var are not ALL of a type equal to vartype.
    """
    #NEED TO FIX: MAKE SURE THAT A RANGE OF TYPES CAN BE TESTED FOR, SUCH AS
    #float, np.float32, np.float64... should all pass as a float.
    import numpy as np
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")
    if (isinstance(var,list) != True) and (isinstance(var,np.ndarray) != True):
        raise Exception("Input error in typetest_array: %s should be of class list or numpy array." % varname)
    for i in range(0,len(var)):
        typetest(var[i],vartype,varname='element %s of %s' % (i,varname))

def dimtest(var,sizes,varname=''):
    """This program tests the dimensions and shape of the input array var.
    Sizes is the number of elements on each axis.
    The program uses the above type tests to make sure that the input is ok.
    If an element in sizes is set to zero, that dimension is not checked against.
    Example:
    import numpy as np
    a=[[1,2,3],[4,3,9]]
    b=np.array(a)
    dimtest(a,[2,3])
    dimtest(a,[3,10])
    """
    import numpy as np
    typetest(sizes,list,varname='sizes in dimtest')
    typetest_array(sizes,int,varname='sizes in dimtest')

    ndim=len(sizes)

    #First check that the number of dimensions is correct
    if np.ndim(var) != ndim:
        raise Exception("Dimension Error:  Variable %s ndim = %s but was required to be %s." % (varname,np.ndim(var),ndim))

    #Then check that each of the axes match.
    sizes_var=np.shape(var)
    for i in range(0,len(sizes)):
        if sizes[i] < 0:
            raise Exception("Sizes in dimtest was not set correctly. It contains negative values. (%s)" % sizes(i))
        if sizes[i] > 0:
            if sizes[i] != sizes_var[i]:
                raise Exception("Dimension Error: Axis %s of variable %s contains %s elements, but %s were required." % (i,varname,sizes_var[i],sizes[i]))

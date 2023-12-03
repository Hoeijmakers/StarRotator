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

def test_smoothing():
    from StarRotator import StarRotator
    import matplotlib.pyplot as plt
    import scipy.interpolate as interp
    import numpy as np
    from lib.integrate import statusbar as statusbar
    import sys
    import copy
    KELT9 = StarRotator(586.0,592.0,100.0)
    N = len(KELT9.spectra)
    F_unsmooth_1 = KELT9.spectra[0]


    wl_low = copy.deepcopy(KELT9.wl)
    wl_high = np.linspace(np.min(KELT9.wl),np.max(KELT9.wl),num=len(KELT9.wl)*30.0)
    fx_high = np.zeros((N,len(wl_high)))



    for i in range(N):
        fx_high[i]=interp.interp1d(KELT9.wl,KELT9.spectra[i])(wl_high)
        statusbar(i,N)


    Res_1 = KELT9.residuals()
    KELT9.spectral_resolution(20000.0)
    Res_2 = KELT9.residuals()
    KELT9.wl = copy.deepcopy(wl_high)
    KELT9.spectra = copy.deepcopy(fx_high)
    KELT9.stellar_spectrum = fx_high[0]
    KELT9.spectra_smooth = copy.deepcopy(fx_high)
    KELT9.spectral_resolution(20000.0)
    # F_smooth = KELT9.spectra[0]
    # Res_2 = KELT9.residuals()
    Res_3 = KELT9.residuals()



    # plt.plot(KELT9.wl,F_unsmooth_1)
    # plt.plot(KELT9.wl,F_smooth)
    # plt.show()

    fig = plt.figure(figsize=(12,6))
    plt.plot(wl_low,Res_1[10],label='No blurring')
    plt.plot(wl_low,Res_2[10],label='Blurring, normal sampling')
    plt.plot(wl_high,Res_3[10],'--',label='Blurring, 30x oversampling')
    plt.xlim((588.8,590.2))
    plt.ylabel('In/out Residual')
    plt.xlabel('Wavelength (nm)')
    plt.title('With oversampling')
    plt.legend()
    plt.show()



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

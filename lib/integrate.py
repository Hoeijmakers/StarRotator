def statusbar(i,x):
    """A beautiful percentage indicator for for-loops.
        Parameters
        ----------
        i : float
            The iteration number of the loop.

        x : float or list
            The list that is being iterated over in the loop, or its number of elements.
            The former is useful if you are doing something like 'for i in list'.
        """
    if type(x) == int:
        print('  '+f"{i/(float(x))*100:.1f} %", end="\r")
    else:
        print('  '+f"{i/(len(x))*100:.1f} %", end="\r")#Statusbar.

def input_tests_local(xp,yp,RpRs):
    """This wraps input tests for the local integrator functions, as these are the same in both cases."""
    import lib.test as test
    #xp and yp should be ints or floats, and may not be NaN.
    test.typetest(xp,[int,float],varname='xp in build_local_spectrum_fast')
    test.typetest(yp,[int,float],varname='xp in build_local_spectrum_fast')
    test.nantest(xp,varname='xp in build_local_spectrum_fast')
    test.nantest(yp,varname='yp in build_local_spectrum_fast')

    #RpRs should be a non-negative float, and may not be NaN:
    test.typetest(RpRs,float,varname='RpRs in build_local_spectrum_fast')
    test.nantest(RpRs,varname='RpRs in build_local_spectrum_fast')
    test.notnegativetest(RpRs,varname='RpRs in build_local_spectrum_fast')


def input_tests_global(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid,fname=''):
    """This wraps the input tests for all the integrator functions, as these are all the same."""
    import lib.test as test
    import numpy as np
    import lib.operations as ops
    #wl and fx should be numpy arrays with the same length.
    test.typetest(wl,np.ndarray,varname='wl in '+fname)
    test.typetest(fx,np.ndarray,varname='fx in '+fname)
    test.dimtest(fx,[len(wl)],varname='fx in '+fname)

    #wlmin and wlmax should be floats or ints, and wlmin < wlmax:
    test.typetest(wlmin,[int,float],varname='wlmin in '+fname)
    test.typetest(wlmax,[int,float],varname='wlmax in '+fname)
    if wlmax <= wlmin:
        raise ValueError('Wlmax (%s) in '+fname+' should be smaller than wlmin (%s)' % (wlmax,wlmin))

    #x and y should be numpy arrays, contain no NaNs and their dimensions should match
    #the velocity and flux grids, which should be 2D numpy arrays.
    test.typetest(x,np.ndarray,varname='x in '+fname)
    test.typetest(y,np.ndarray,varname='y in '+fname)
    test.nantest(x,varname='x in '+fname)
    test.nantest(y,varname='x in '+fname)
    test.dimtest(vel_grid,[len(y),len(x)],varname='y in '+fname)
    test.dimtest(flux_grid,[len(y),len(x)],varname='y in '+fname)
    test.typetest(vel_grid,np.ndarray,varname='vel_grid in '+fname)
    test.typetest(flux_grid,np.ndarray,varname='flux_grid in '+fname)

    dvm = np.nanmax(np.abs(vel_grid))
    if wlmin/ops.doppler(dvm) <= np.nanmin(wl):
        raise Exception("Value error in "+fname+": wlmin (%s) is smaller than min(wl) (%s) after accounting for the maximum velocity shift (%s)" % (wlmin/ops.doppler(dvm),min(wl),dvm))
    if wlmax*ops.doppler(dvm) >= np.nanmax(wl):
        raise Exception("Value error in "+fname+": wlmax (%s) is greater than max(wl) (%s) after accounting for the maximum velocity shift (%s)" % (wlmax*ops.doppler(dvm),max(wl),dvm))



    #fname should be a string.
    test.typetest(fname,str,varname='fname in input_tests_global')

def build_local_spectrum_fast(xp,yp,zp,RpRs,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the fast and easy way of building the missing local spectrum that
    is blocked by the planet, assuming that the projected spin axis is vertical
    and differential rotation is zero, i.e. that all rows contain the same velocity.
        Parameters
        ----------
        xp : int,float
            The x position of the planet in units of x (see below).

        yp : int,float
            The y position of the planet in units of y (see below).

        RpRs : float
            The radius of the planet in units of stellar radii.

        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        wlmin: float
            The minimum wavelength to be considered, in units of wl.

        wlmax: float
            The maximum wavelength to be considered, in units of wl.

        x,y:
            The x and y axis of the velocity grid, typically in units of stellar
            radius.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

        flux_grid: 2D np.array()
            The broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum,
        flux: float
            The lightcurve of the remainder of the stellar disk.
        mask:
            The inverse mask that shows the stellar disk with the planet masked out.

    """
    import numpy as np
    import lib.operations as ops
    import warnings
    import pdb
    import matplotlib.pyplot as plt

    #Standard tests on input:
    input_tests_local(xp,yp,RpRs)
    input_tests_global(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid,fname='build_local_spectrum_fast')

    #We start by creating a mask that selects just the area of the star that is
    #covered by the planet, as well as its inverse which we like to return for
    #plotting purposes.
    x_full = np.tile(x,(len(y),1))-xp
    y_full = np.tile(y,(len(x),1)).T-yp
    d = np.sqrt(x_full**2 + y_full**2)
    di = d*1.0
    d[d > RpRs] = np.nan#Mask out everything but the location of the planet.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if np.nanmin(d) <= RpRs:
            di[d <= RpRs] = np.nan#out the location of the planet. Inverse of d.
    if zp < 0.0:
        d*=np.nan#Force the planet to not be in front of the disk if the z-coordinate is such that it is behind the star.
        di[np.isnan(di)]=1.0#started as 1.0, set back to 1.0
    mask = flux_grid*(d*0.0+1.0)#Set that to 1.0 and multiply with flux grid. Nansum coming!
    mask_i = flux_grid*(di*0.0+1.0)#Inverse of mask.
    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))


    F = 0#output
    flux = np.nansum(mask,axis = 0)#This is the sum of the flux grid.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        v=np.nanmedian(vel_grid,axis = 0)
    #Here we only loop over columns. The summation in the vertical direction is contained
    #in the flux variable.
    for i in range(len(x)):
        if np.isnan(v[i]) == False:
            F+=ops.shift(wlc,wlc_wide,fxc_wide,v[i])*flux[i]

    return(wlc,F,np.nansum(mask_i),(di*0.0)+1.0)


def build_local_spectrum_slow(xp,yp,zp,RpRs,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the brute-force (slow) way of building the missing local spectrum that
    is blocked by the planet, applicable to any velocity/flux grid. This is used when differential
    rotation is a thing.

        Parameters
        ----------
        xp : float
            The x position of the planet in units of x (see below).

        yp : float
            The y position of the planet in units of y (see below).

        RpRs : float
            The radius of the planet in units of stellar radii.

        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        wlmin: float
            The minimum wavelength to be considered, in units of wl.

        wlmax: float
            The maximum wavelength to be considered, in units of wl.

        x,y:
            The x and y axis of the velocity grid, typically in units of stellar
            radius.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

        flux_grid: 2D np.array()
            The broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum,
        flux: float
            The lightcurve of the remainder of the stellar disk.
        mask:
            The inverse mask that shows the stellar disk with the planet masked out.

    """
    import numpy as np
    import lib.operations as ops
    import warnings
    #Standard tests on input:
    input_tests_local(xp,yp,RpRs)
    input_tests_global(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid,fname='build_local_spectrum_slow')


    #We start by creating a mask that selects just the area of the star that is
    #covered by the planet, as well as its inverse which we like to return for
    #plotting purposes.
    x_full = np.tile(x,(len(y),1))-xp
    y_full = np.tile(y,(len(x),1)).T-yp
    d = np.sqrt(x_full**2 + y_full**2)
    if zp < 0: d+=np.inf
    di = d*1.0
    d[d > RpRs] = np.nan#Mask out everything but the location of the planet.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if np.nanmin(d) <= RpRs:
            di[d <= RpRs] = np.nan#out the location of the planet. Inverse of d.
    if zp < 0.0:
        d*=np.nan#Force the planet to not be in front of the disk if the z-coordinate is such that it is behind the star.
        di[np.isnan(di)]=1.0#started as 1.0, set back to 1.0
    mask = flux_grid*(d*0.0+1.0)#Set that to 1.0 and multiply with flux grid. Nansum coming!
    mask_i = flux_grid*(di*0.0+1.0)#Inverse of mask.

    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    flux = np.nansum(mask,axis = 0)#This is the sum of the flux grid.

    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(mask[j,i]) == False:
                F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
        # statusbar(i,len(x))
    return(wlc,F,np.nansum(mask_i),(di*0.0)+1.0)






def build_spectrum_fast(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the fast and easy way of building the stellar spectrum assuming
    that the projected spin axis is vertical and differential rotation is zero,
    i.e. that all rows contain the same velocity.
        Parameters
        ----------
        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        wlmin: float
            The minimum wavelength to be considered, in units of wl.

        wlmax: float
            The maximum wavelength to be considered, in units of wl.

        x,y:
            The x and y axis of the velocity grid, typically in units of stellar
            radius.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

        flux_grid: 2D np.array()
            The broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum.
    """
    import numpy as np
    import lib.operations as ops
    import warnings
    import pdb
    #Standard tests on input:
    input_tests_global(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid,fname='build_spectrum_fast')

    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))
    F = 0#output
    flux = np.nansum(flux_grid,axis = 0)#This is the sum of the flux grid.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        v=np.nanmedian(vel_grid,axis = 0)
    v[(np.isnan(v))] = 0.0
    for i in range(len(x)):
        if np.isnan(v[i]) == False:
            F+=ops.shift(wlc,wlc_wide,fxc_wide,v[i])*flux[i]
        statusbar(i,len(x))
    return(wlc,F)






def build_spectrum_slow(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the default, brute-force way of integrating the spectrum.
        Parameters
        ----------
        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        wlmin: float
            The minimum wavelength to be considered, in units of wl.

        wlmax: float
            The maximum wavelength to be considered, in units of wl.

        x,y:
            The x and y axis of the velocity grid, typically in units of stellar
            radius.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

        flux_grid: 2D np.array()
            The broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.
    """

    import numpy as np
    import lib.operations as ops
    import lib.stellar_spectrum as spectrum
    input_tests_global(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid,fname='build_spectrum_slow')
    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    # start = time.time()
    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(vel_grid[j,i]) == False:
                F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
        statusbar(i,len(x))
    # print(time.time()-start)
    return(wlc,F)



def build_spectrum_limb_resolved(wl,fx_list,mu_list,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """WRITE THIS.
        Parameters
        ----------

    """
    #I copy paste as much as possible from above. The roles of wlc and wlc_wide have
    #changed because wl,fx is already narrow by construction. Therefore in order to
    #crop the spectrum with margin, I actually need to crop the wl axis inwards.
    #To do this, I needed to convert clip_spectrum to crop_spectrum.
    import numpy as np
    import lib.operations as ops
    import lib.stellar_spectrum as spectrum
    import sys
    import lib.test as test
    wlc_wide = wl

    #Standard tests on input
    test.typetest(mu_list,np.ndarray,varname='mu_list in build_spectrum_limb_resolved')
    test.dimtest(mu_list,[len(fx_list)],varname='mu_list in build_spectrum_limb_resolved')
    input_tests_global(wlc_wide,fx_list[0],wlmin,wlmax,x,y,vel_grid,vel_grid,fname='build_spectrum_limb_resolved')
    wlc,fxc = ops.crop_spectrum(wl,fx_list[0],1.0*np.nanmax(np.abs(vel_grid)))#I do this for only one spectrum because I only care about wlc
    F = 0#output
    # start = time.time()
    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(vel_grid[j,i]) == False:
                mu = 1 - np.sqrt(x[i]**2 + y[j]**2) # mu angle corresponding to the pixel
                diff = abs(mu_list - mu) # chooses nearest mu
                index = np.argmin(diff)
                if mu_list[index] > 0:
                    fxc_wide = fx_list[index]
                    # print(i,j,x[i],y[j],mu,mu_list[index])
                    F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
        statusbar(i,len(x))
    # print(time.time()-start)
    return(wlc,F)


def build_local_spectrum_limb_resolved(xp,yp,zp,RpRs,wl,fx_list,mu_list,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is used when SPECTRUM was called to provide spectra wth mu angles, in which
    case the flux map doesn't exist and the computation is different from the
    disk-averaged PHOENIX spectra.

    WRITE THIS MORE

        Parameters
        ----------
        xp : float
            The x position of the planet in units of x (see below).

        yp : float
            The y position of the planet in units of y (see below).

        RpRs : float
            The radius of the planet in units of stellar radii.

        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        wlmin: float
            The minimum wavelength to be considered, in units of wl.

        wlmax: float
            The maximum wavelength to be considered, in units of wl.

        x,y:
            The x and y axis of the velocity grid, typically in units of stellar
            radius.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.


        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum,
        flux: float
            The lightcurve of the remainder of the stellar disk.
        mask:
            The inverse mask that shows the stellar disk with the planet masked out.

    """
    import numpy as np
    import lib.operations as ops
    import lib.test as test
    import warnings
    import sys
    import pdb
    wlc_wide = wl
    #Standard tests on input:
    input_tests_local(xp,yp,RpRs)
    input_tests_global(wl,fx_list[0],wlmin,wlmax,x,y,vel_grid,vel_grid,fname='build_local_spectrum_limb_resolved')
    test.typetest(mu_list,np.ndarray,varname='mu_list in build_local_spectrum_limb_resolved')
    test.dimtest(mu_list,[len(fx_list)],varname='mu_list in build_local_spectrum_limb_resolved')

    #We start by creating a mask that selects just the area of the star that is
    #covered by the planet, as well as its inverse which we like to return for
    #plotting purposes.
    x_full = np.tile(x,(len(y),1))-xp
    y_full = np.tile(y,(len(x),1)).T-yp
    d = np.sqrt(x_full**2 + y_full**2)
    if zp < 0: d+=np.inf
    di = d*1.0
    d[d > RpRs] = np.nan#Mask out everything but the location of the planet.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if np.nanmin(d) <= RpRs:
            di[d <= RpRs] = np.nan#out the location of the planet. Inverse of d.
    if zp < 0.0:
        d*=np.nan#Force the planet to not be in front of the disk if the z-coordinate is such that it is behind the star.
        di[np.isnan(di)]=1.0#started as 1.0, set back to 1.0
    mask = flux_grid*(0.0*vel_grid+1.0)*(d*0.0+1.0)#Set that to 1.0 and multiply with flux grid. Nansum coming!
    mask_i = flux_grid*(0.0*vel_grid+1.0)*(di*0.0+1.0)#Inverse of mask.
#MOVE ALL OF THIS INTO A WRAPPER? I THINK THIS IS NOW COPYPASTED 3 TIMES?
#ACTUALLY, NOW I AM USING VEL GRID INS TEAD OF FLUX GRID. THAT MEANS THAT THE
#FLUX IS NO LONGER CALCULATED PROPERLY. NO LIGHTCURVE!
    # wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))
    wlc,fxc = ops.crop_spectrum(wl,fx_list[0],1.0*np.nanmax(np.abs(vel_grid)))#I do this for only one spectrum because I only care about wlc


    F = 0#output
    # start = time.time()

    # pdb.set_trace()
    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(mask[j,i]) == False:
                mu = 1 - np.sqrt(x[i]**2 + y[j]**2)
                diff = np.abs(mu_list - mu)
                index = np.argmin(diff)
                if mu_list[index] > 0:
                    fxc_wide = fx_list[index]
                    # print(i,j,x[i],y[j],mu,mu_list[index])
                    F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
    # print(time.time()-start)
    return(wlc,F,np.nansum(mask_i),(di*0.0)+1.0)



    # F = 0#output
    # flux = np.nansum(mask,axis = 0)#This is the sum of the flux grid.
    # for i in range(len(x)):
    #     for j in range(len(y)):
    #         if np.isnan(mask[j,i]) == False:
    #             F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
    #     # statusbar(i,len(x))
    # return(wlc,F,np.nansum(mask_i),(di*0.0)+1.0)

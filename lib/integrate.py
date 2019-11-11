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
        print('  '+f"{i/(float(x)-1)*100:.1f} %", end="\r")
    else:
        print('  '+f"{i/(len(x)-1)*100:.1f} %", end="\r")#Statusbar.



def build_local_spectrum_fast(xp,yp,RpRs,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the fast and easy way of building the missing local spectrum that
    is blocked by the planet, assuming that the projected spin axis is vertical
    and differential rotation is zero, i.e. that all rows contain the same velocity.
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

    #We start by creating a mask that selects just the area of the star that is
    #covered by the planet, as well as its inverse which we like to return for
    #plotting purposes.
    x_full = np.tile(x,(len(x),1))-xp
    y_full = np.tile(y,(len(y),1)).T-yp
    d = np.sqrt(x_full**2 + y_full**2)
    di = d*1.0
    d[d > RpRs] = np.nan#Mask out everything but the location of the planet.
    di[d <= RpRs] = np.nan#out the location of the planet. Inverse of d.
    mask = flux_grid*(d*0.0+1.0)#Set that to 1.0 and multiply with flux grid. Nansum coming!
    mask_i = flux_grid*(di*0.0+1.0)#Inverse of mask.


    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    flux = np.nansum(mask,axis = 0)#This is the sum of the flux grid.
    v=np.nanmedian(vel_grid,axis = 0)
    # v[(np.isnan(v))] = 0.0
    #Here we only loop over columns. The summation in the vertical direction is contained
    #in the flux variable.
    for i in range(len(x)):
        if np.isnan(v[i]) == False:
            F+=ops.shift(wlc,wlc_wide,fxc_wide,v[i])*flux[i]
    return(wlc,F,np.nansum(mask_i),(di*0.0)+1.0)


def build_local_spectrum_slow(xp,yp,RpRs,wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the brute-force (slow) way of building the missing local spectrum that
    is blocked by the planet, applicable to any velocity/flux grid.

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

    #We start by creating a mask that selects just the area of the star that is
    #covered by the planet, as well as its inverse which we like to return for
    #plotting purposes.
    x_full = np.tile(x,(len(x),1))-xp
    y_full = np.tile(y,(len(y),1)).T-yp
    d = np.sqrt(x_full**2 + y_full**2)
    di = d*1.0
    d[d > RpRs] = np.nan#Mask out everything but the location of the planet.
    di[d <= RpRs] = np.nan#out the location of the planet. Inverse of d.
    mask = flux_grid*(d*0.0+1.0)#Set that to 1.0 and multiply with flux grid. Nansum coming!
    mask_i = flux_grid*(di*0.0+1.0)#Inverse of mask.

    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    flux = np.nansum(mask,axis = 0)#This is the sum of the flux grid.

    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(mask[j,i]) == False:
                F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
        statusbar(i,len(x))
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
    import time
    import matplotlib.pyplot as plt
    import lib.plotting as plotting

    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    start = time.time()
    flux = np.nansum(flux_grid,axis = 0)#This is the sum of the flux grid.
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
    import time
    import matplotlib.pyplot as plt

    wlc,fxc,wlc_wide,fxc_wide = ops.clip_spectrum(wl,fx,wlmin,wlmax,pad=2.0*np.nanmax(np.abs(vel_grid)))

    F = 0#output
    start = time.time()
    for i in range(len(x)):
        for j in range(len(y)):
            if np.isnan(vel_grid[j,i]) == False:
                F+=ops.shift(wlc,wlc_wide,fxc_wide,vel_grid[j,i])*flux_grid[j,i]
        statusbar(i,len(x))
    print(time.time()-start)
    return(wlc,F)

def statusbar(i,x):
    if type(x) == int:
        print('  '+f"{i/(float(x)-1)*100:.1f} %", end="\r")
    else:
        print('  '+f"{i/(len(x)-1)*100:.1f} %", end="\r")#Statusbar.






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

        args:
            The arguments class passed to main.py.

        x,y:
            The x and y axis of the velocity grid.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

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
        F+=ops.shift(wlc,wlc_wide,fxc_wide,v[i])*flux[i]
    statusbar(i,len(x))
    print(time.time()-start)

    return(wlc,F)






def build_spectrum_slow(wl,fx,wlmin,wlmax,x,y,vel_grid,flux_grid):
    """This is the default, brute-force way of integrating the spectrum.
        Parameters
        ----------
        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
            The stellar model flux.

        args:
            The arguments class passed to main.py.

        x,y:
            The x and y axis of the velocity grid.

        vel_grid: 2D np.array()
            The velocity grid of the stellar disk. Values outside of the disk are
            to be set to NaN.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum.
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

import jax
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.dynamics import doppler_shift
from starrotator.lib.operations import vert_int_q_ld, circ_int_q_ld, vert_int_q_ld_bounded
from starrotator.lib.vgrid import calc_vel_stellar
from starrotator.lib.vgrid import calc_flux_stellar



@partial(jax.jit,static_argnames=["N","norm"])
def sum_stellar_spectrum_v1(wl,fx,vel_eq,i_stellar,a1,a2,N=400,norm=False):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a SIMPLE VELOCITY GRID (that means: without drr)
        and a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation is 
        modelled by broad-band limb darkening. This version of the integration is fastest because
        it uses the axial symmetry of the disk, so it only computes the velocity of the disk along
        the projected equator and the limb darkening function is integrated analytically vertically.

        This is the fastest possible approximation.

        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths.

        fx : array-like
            The stellar model fluxes corresponding to wl.

        x : array-like
            The horizontal grid axis. Equivalent to x,y in the 2D grid computation of the other integration versions.

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees.

        a1 : float
            Linear limb darkening coefficient.

        a2 : float
            Quadratic limb darkening coefficient.

        norm : bool
            If set to true, the weights are normalised such that all fluxes sum to 1.0

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
    """
    x = jnp.linspace(-1,1,N)
    dx = x[1]-x[0]
    v_axis = x*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.
    weights = vert_int_q_ld(x,a1,a2)

    fx_shifted = doppler_shift(wl,fx,v_axis).T*weights

    if norm:
        F_out = jnp.nansum(fx_shifted,axis=1)/jnp.nansum(weights)
    else:
        F_out = jnp.nansum(fx_shifted,axis=1)*2*dx
    return(F_out)#Factor of 2 because this is calculating a semi-circle.





@partial(jax.jit,static_argnames=['batched'])
def sum_stellar_spectrum_v2(wl,fx,vel_grid,flux_grid,batched=True):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a FREE VELOCITY GRID (that means: with drr)
        but a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation that is 
        available, is modelled by a broad-band, quadratic limb darkening model.
    
        In this scenario we cannot assume that any velocity in the grid repeats itself when drr is present, and there is
        probably no symmetry. This means that we need to interpolate the spectrum for all velocities and
        sum over the whole grid. We can do this fully vectorised, or batched (row-by-row in the grid) to 
        not load this huge array into memory. Full vectorization may be a tiny bit faster but may only really be useful on
        big servers with sufficient RAM, so batching=True is the default.

        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths. One-dimensional.

        fx : array-like
            The stellar model fluxes corresponding to wl. One-dimensional.

        vel_grid : array-like
            The 2D velocity grid of the stellar disk. Generally this is a circular map in a square matrix.
            Values outside of the disk (in the corners of the square) are assumed to be set to NaN.

        flux_grid : array-like
            The 2D broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        batched : bool
            Compute the sum row-by-row or entirely at once (warning: very memory hungry).

        Returns
        -------
        fx : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
    """

    if batched: #We use lax.scan to sum over the grid row-by-row to reduce memory load.
        def scan_fn(carry, inputs):
            vel_row, flux_row = inputs  # Each is a row of len(y)
            fx_shifted = doppler_shift(wl, fx, vel_row).T * flux_row  # shape (len(wl), len(y))
            partial_sum = jnp.nansum(fx_shifted, axis=1)  # shape (N_wl,)
            new_total = carry + partial_sum
            return new_total, None  # No output needs to be collected.
        # 
        # Initial spectrum is all zeros (shape of output)
        init_spectrum = jnp.zeros_like(fx)
        # 
        # Scan over axis 0 (rows)
        F_out, _ = lax.scan(scan_fn, init_spectrum, (vel_grid, flux_grid))
        return(F_out)
    else: # If a LOT of memory is available, then we may not need to batch and do it directly.
        # This is very memory-hungry because the Nx n Ny x Nwl array gets done all at once.
        # You very likely want to set batched=True instead....
        fx_shifted = doppler_shift(wl,fx,jnp.ravel(vel_grid)).T*jnp.ravel(flux_grid)
        F_out = jnp.nansum(fx_shifted,axis=1)
        return(F_out)
    


# Note that this approximation may be very good even for the case of mu-dependent spectra without drr 
# because the range of mu-angles covered by the planet is small. Mu dependence then still comes in via
# the change in xp,yp but that can be addressed by either:
# -writing the function such that accepts an fx-array and a mu-array. The switching logic may make this hard to jit.
# -calling this function in a jitted loop or vmap over different fx-es, corresponding to the mu's the planet traverses.
@partial(jax.jit,static_argnames=['N'])
def sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100):
    """This builds the stellar spectrum that is hidden behind the planet by doppler-shifting 
        and interpolating the input spectrum, summing only over the range of coordinates that is
        obscured by the planet disk.

        This version of the integration assumes a SIMPLE VELOCITY GRID (that means: without drr)
        and a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation is 
        modelled by broad-band limb darkening. This version of the integration is fastest because
        it uses the axial symmetry of the disk, so it only computes the velocity of the disk along
        the projected equator, and the limb darkening function is integrated analytically.

        This is the fastest possible approximation.

        This function integrates in a manner that is equivalent to sum_stellar_spectrum_v1.

        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths.

        fx : array-like
            The stellar model fluxes corresponding to wl.

        xp : array-like
            The horizontal location of the planet, either as a float or as a 1D array.

        yp : array-like
            The vertical location of the planet, either as a float or as a 1D array.

        Rp : float
            The radois of the planet in units of stellar radii.

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees.

        a1 : float
            Linear limb darkening coefficient.

        a2 : float
            Quadratic limb darkening coefficient.

        N : int
            Number of grid-points in the x direction onto which the planet is simulated. Note that the resolution of the planet is independent
            of the resolution of the grid of the star. Note that as this calculation is batched over planet positions, cases where
            you have many planet positions and high N will result in high memory load if batched is set to false.

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is a 2D array as xp,yp are 1D). This is scaled to an arbitrary number, but this scaling corresponds
            to the scaling of the full disk. Normalisation is therefore achieved by dividing by the integral of the full
            disk, that depends only on a1 and a2.
    """
    x = jnp.linspace(-1,1,N)
    dxR = (x[1]-x[0])*Rp
    x_array = x*Rp + xp[:,None]
    y_max = yp[:,None] + jnp.sqrt(1-x**2)*Rp
    y_min = yp[:,None] - jnp.sqrt(1-x**2)*Rp
    y_star_max = jnp.sqrt(1-x_array**2)*0.999999
    y_star_min = -jnp.sqrt(1-x_array**2)*0.999999
    y_final_max = jnp.clip(y_max,y_star_min,y_star_max)
    y_final_min = jnp.clip(y_min,y_star_min,y_star_max)
    v_axes = x_array*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.
    weights = vert_int_q_ld_bounded(x_array,y_final_min,y_final_max,a1,a2)
    fx_shifted = doppler_shift(wl,fx,v_axes).T*weights # /jnp.nansum(weights)
    F_out = jnp.nansum(fx_shifted,axis=2).T * dxR #Multiply with the differential for normalization.
    return(F_out)




@partial(jax.jit,static_argnames=['batched'])
def sum_hidden_spectrum_v2(wl,fx,vel_grid_array,flux_grid_array,batched=True):
    """This builds the stellar spectrum that is hidden behind the planet by doppler-shifting 
        and interpolating the input spectrum, summing only over the range of coordinates that is
        obscured by the planet disk.

        This version of the integration assumes a FREE VELOCITY GRID (that means: with drr)
        and a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation is 
        modelled by broad-band limb darkening. 

        This function integrates in a manner that is equivalent to sum_stellar_spectrum_v2.

        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths.

        fx : array-like
            The stellar model fluxes corresponding to wl.

        vel_grid : array-like
            An array of 2D velocity grids of the stellar disk that are hidden behind the planet.
            Generally this is a circular map in a square matrix.
            Values outside of each disk (in the corners of the square) are assumed to be set to NaN,
            as well as values that are outside of the stellar disk. Use the function 
            create_hidden_grid_array() to make this out of known planet-star parameters.

        flux_grid : array-like
            An array of the 2D broad-band flux maps of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        batched : bool
            Compute the integration phase-by-phase (True) or all at once (warning: very memory hungry).


        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is generally a 2D array). Note that the result is unnormalised and therefore depends on the
            normalization of the flux grid array and of the input spectrum.
    """

    def scan_fn(carry, inputs):
        vel_i, flux_i = inputs  # Each is a planetary flux grid.
        sum_i = sum_stellar_spectrum_v2(wl,fx,vel_i,flux_i,batched=batched)
        carry += 1
        return carry, sum_i  # No output needs to be collected.

    _, F_out = lax.scan(scan_fn, 0, (vel_grid_array.T, flux_grid_array.T))

    #F_out_norm = F_out * dxR**2 / int_analytical 

    return(F_out)





@partial(jax.jit,static_argnames=['N'])
def create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=50):
    """This builds arrays of the flux grid and velocity grid of the stellar disk that are
        hidden behind the planet.

        This version creates a FREE VELOCITY GRID (that means: with possible drr), and the 
        possibility to add center-to-limb variation is modelled by broad-band limb darkening. 

        The output are 3D arrays of flux and velocity grids, one for each phase of the planet,
        with areas that are outside the planet or stellar disk set to NaN. These are designed
        to be passed into sum_hidden_spectrum_v2.

        Parameters
        ----------
        xp : array-like
            The horizontal location of the planet, either as a float or as a 1D array.

        yp : array-like
            The vertical location of the planet, either as a float or as a 1D array.

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees.

        diff_rot_rate : float
            If you set this to zero, v1 may be the more suitable integration method.

        a1 : float
            Linear limb darkening coefficient.

        a2 : float
            Quadratic limb darkening coefficient.

        N : int
            Number of grid-points onto which the planet is simulated. Note that the resolution of 
            the planet is independent of the resolution of the grid of the star. Also note that 
            as this calculation is batched over planet positions, cases where you have many planet 
            positions and high N will result in high memory load if batched is set to false.


        Returns
        -------
        flux_grid_array : array
            3D array of flux grids that are hidden by the planet, meant as input for sum_hidden_spectrum_v2
            
        vel_grid_array : array
            3D array of velocity grids that are hidden by the planet, meant as input for sum_hidden_spectrum_v2

        differential : float
            The size of the differential (dx*R) that may be used to renormalise the flux array, since it has an 
            aritrary size compared to that of the stellar disk.
    """

    x = jnp.linspace(-1,1,N)
    dx = x[1]-x[0]

    x_array = x*Rp + xp[:,None] # Rescaling the grid to only encompass the planet (much finer than star typically)
    y_array = x*Rp + yp[:,None] # Note that we do a circle in a square so y=x.
    dxR = dx*Rp # Rescaling the size of the differential.
    d = jnp.sqrt(x[None,:]**2 + x[:,None]**2)
    mask = jnp.where(d > 1,jnp.nan,d)*0.0+1.0

    # int_analytical = circ_int_q_ld(a1,a2)

    flux_disk_array  = calc_flux_stellar(x_array.T,y_array.T,a1,a2,norm=False) #This is really, really awesome. #Need to figure out how to scale these slices using the total integral of the LD profile.
    vel_disk_array  = calc_vel_stellar(x_array.T,y_array.T,i_stellar, vel_eq, diff_rot_rate)


    flux_disk_array_masked = flux_disk_array*mask.T[:,:,None] #This crops out the circle that is the planet.
    vel_disk_array_masked = vel_disk_array*mask.T[:,:,None]

    return(flux_disk_array_masked,vel_disk_array_masked,dxR)













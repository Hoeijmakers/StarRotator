import jax
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.dynamics import doppler_shift
from starrotator.lib.operations import vert_int_q_ld
# CONTINUE HERE NEXT TIME
# THESE ARRAYS CAN PROBABLY BE FLATTENED, 1D. BUT SHAPE MAYBE DOESNT MATTER?
# IN CASE OF NO DRR, THE VEL GRID SHOULD BE MUCH(!) SIMPLIFIED (IN QUARTERS).
# AND ONLY ONE INTERPOLATION PER VELOCITY (DUPLICATING INTERPOLATION IN THE VERTICAL DIRECTION IS UNNECESSARY, 
# ONLY COUNTING HOW MANY TIMES A PARTICULAR VELOCITY IS IN THE GRID, AND THAT GOES AS SIN(x) OR SMTH.
# NOTICE THAT IN THIS CASE, HAVING x=0 AND y=0 IS A PROBLEM (FOR MIRRORING THE QUARTERS).
@jit
def sum_stellar_spectrum_v1(wl,fx,x,vel_eq,i_stellar,a1,a2):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a SIMPLE VELOCITY GRID (that means: without drr)
        and a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation is 
        modelled by broad-band limb darkening. This version of the integration is fastest because
        it uses the axial symmetry of the disk, so it only computes the velocity of the disk along
        the projected equator.

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

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
    """
    v_axis = x*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.
    weights = vert_int_q_ld(x,a1,a2)
    fx_shifted = doppler_shift(wl,fx,v_axis).T*weights/jnp.nansum(weights)
    F_out = jnp.nansum(fx_shifted,axis=1)
    return(F_out)





@partial(jax.jit,static_argnames=['batched'])
def sum_stellar_spectrum_v2(wl,fx,vel_grid,flux_grid,batched=True):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a FREE VELOCITY GRID (that means: with drr)
        but a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation that is 
        available, is modelled by a broad-band, quadratic limb darkening model.
    
        In this scenario we cannot assume that any velocity in the grid repeats itself when drr is present, and there is
        probably no symmetry. This means that we need to interpolate the spectrum for all velocities and
        sum over the whole grid. We can do this directly, or batched (row-by-row in the grid) to 
        not load this huge array into memory. That may be a tiny bit faster but may only really be useful on
        big servers, so batching=True is the default.

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









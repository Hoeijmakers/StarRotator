from jax import jit, lax
from jax import numpy as jnp
from jaxoplanet.orbits import keplerian
from starrotator.lib.constants import rad_in_deg, R_sun, d_in_seconds, c as c_light

# @jit
# def pos_eccentric(phase, M = 1.0, m = 0.0, P = 365.0, e = 0.0, omega = 0.0,
#                 i=90.0):
#     """
#     This function calculates the radial velocity in km/s for a planet in 
#     an elliptical orbit using jaxoplanet. 
    
#     Input is provided in terms of the orbital phases at which the radial
#     velocity is required, and the system parameters, including the 
#     stellar mass. Stellar mass is assumed to be known to better precision 
#     than the semi-major axis a. If this is not the case, you need to 
#     proceed by calculating M from a, using Kepler III.

#     If the mass of the companion is non-negligible, then its mass can
#     also be set.

#     The coordinate system follows that of the exoplanet package:
#     https://docs.exoplanet.codes/en/latest/tutorials/data-and-models/


#     Parameters
#     ----------
#     phase : float, array-like
#         Orbital phase, typically between 0 and 1.0. 0.0 is mid-transit.
#         Equivalent to the Mean Anomaly divided by 2 pi.

#     M : float
#         Stellar mass in solar masses.

#     m : float
#         Planet mass in solar masses.

#     P : float
#         Orbital period in days.

#     e : float
#         eccentricity.

#     omega : float
#         argument of peri-apsis, following the exoplanet package coordinate system. 

#     i : float
#         Orbital inclination in degrees. 90 is transiting.

#     Returns
#     -------
#     rv : float, array-like, same as phase
#         The planet's radial velocity in km/s.

#     """


#     star = keplerian.Central(mass=M, radius=1.0)

#     system = keplerian.System(star).add_body(mass = m, period = P, 
#                         eccentricity = e, omega_peri = omega/rad_in_deg,
#                         inclination = i/rad_in_deg,time_transit=0.0
#                         )

#     planet = system.bodies[0]

#     # NEED TO CHANGE THIS TO PULL OUT POSITIONS IN X AND Y INSTEAD.
#     vz = planet.velocity(phase * P)[2] # along the z axis in Rsol per day.

#     # Convert to km/s, note that the positive z axis is the negative RV axis.
#     rv = vz * R_sun / 1e5 / d_in_seconds * -1 
   
#     return(xp,yp)




@jit
def doppler_factor(v):
    """
    This calculates the relativistic doppler factor given a radial velocity 
    in km/s.


    Parameters
    ----------
    v : float, array-like
        Radial velocity in km/s.

    Returns
    -------
    f : float, array-like, same as v
        The doppler factor to be multiplying the wavelength with.

    """

    beta = v * 1e5 / c_light # c is in cgs so convert v to km/s.
    f = jnp.sqrt((1 + beta)/(1 - beta))
    return(f)







@jit
def doppler_shift(wl,fx,v):
    """Doppler-shift a spectrum via linear interpolation. 
    
    This is the simplest,
    most generic and default doppler shifting function.

    Parameters
    ----------
    wl : array-like
        The wavelength axis of the spectrum, 1D array.
    fx : array-like
        The corresponding flux axis, 1D array.
    v : float, array-like
        The radial velocity (positive = redshift) in km/s.
        Multiple velocities can be provided in the same array.


    Returns
    -------
    fx_shifted: array-like
        The shifted spectrum evaluated on the original wavelength axis.
        If multiple velocities are given, this is a 2D array matching the 
        number of RVs to be shifted to. Note that edges are not extrapolated
        or handled: Edge values are repeated.
    """
    g = doppler_factor(v)
    fx_shifted = jnp.interp(wl/g[None,:].T, wl, fx)
    return(fx_shifted)


@jit
def doppler_shift_dlogl(dlogl,fx,v):
    """
    Doppler-shift with explicit linear interpolation on a constant-dlog(lambda) grid.

    Doppler-shifts a spectrum via index shifting. The flux spectrum needs to 
    be defined on a wavelength array with a constant log-lambda step. In this case,
    the velocity shift results in an element-wise shift that is independent of wavelength,
    meaning that the sub-grid shift is the same in each element, removing the need for a 
    generic interpolation step (the linear interpolation weights can be calculated once).

    This function is therefore faster than the doppler-shift function above
    (20x speedup with equivalent input measured).

    Parameters
    ----------
    delta : float
        The step-size in loglambda
    fx : array-like
        The corresponding flux axis, 1D array.
    v : float, array-like
        The radial velocity (positive = redshift) in km/s.
        Multiple velocities can be provided in the same array.

    Returns
    -------
    fx_int: array-like
        The shifted spectrum evaluated on the original wavelength axis.
        If multiple velocities are given, this is a 2D array matching the 
        number of RVs to be shifted to. Note that edges are not extrapolated
        or handled: Edge values are repeated.
    """
    nfx = fx.shape[0]
    I0 = jnp.arange(nfx+1,dtype=jnp.int32) #An index array for fx.
    v = jnp.atleast_1d(jnp.asarray(v))
    g = doppler_factor(v)
    D_i = -jnp.log(g)/dlogl #Number of indices to shift. Positive delta_i -> shift to the right in index space
    #This is a decimal number. We are going to decompose the shift as the sum of an integer part plus a decimal part.
    D_int = jnp.floor(D_i).astype(jnp.int32)
    D_rest = D_i % 1
    i_shifted_int = I0[None, :] + D_int[:, None] #These are the indices after shifting by an integer amount. Minus sign to shift to the right.
    i_shifted_int_legal = jnp.clip(i_shifted_int,0,nfx-2)
    fx_shifted = fx[i_shifted_int_legal].T
    fx_int = fx_shifted[0:-1]*(1-D_rest) + fx_shifted[1:]*D_rest#The interpolation step.
    return(fx_int.T)





@jit
def orbit_euclidian(phase, a = 5.0, m = 0.0, P = 4.0, e = 0.0, omega = 0.0,
                i=90.0):
    """
    This function calculates the x, y and z position of a planet in 
    an elliptical orbit using jaxoplanet. 
    
    Input is provided in terms of the orbital phases at which the position 
    is required, and the system parameters.

    If the mass of the companion is non-negligible, then its mass can
    also be set.

    The coordinate system follows that of the exoplanet package:
    https://docs.exoplanet.codes/en/latest/tutorials/data-and-models/


    Parameters
    ----------
    phase : float, array-like
        Orbital phase, typically between 0 and 1.0. 0.0 is mid-transit.
        Equivalent to the Mean Anomaly divided by 2 pi.

    a : float
        Orbital semi-major axis in solar radii.

    m : float
        Planet mass in solar masses.

    P : float
        Orbital period in days.

    e : float
        eccentricity.

    omega : float
        argument of peri-apsis, following the exoplanet package coordinate system. 

    i : float
        Orbital inclination in degrees. 90 is transiting.

    Returns
    -------
    x : float, array-like, same as phase
        The planet's x, y and z-positions in units of solar radii.

    """
    from jaxoplanet.orbits import keplerian
    from starrotator.lib.util import rad_in_deg, R_sun, d_in_seconds, G, M_sun
    import numpy as np


    a_cgs = a*R_sun
    P_cgs = P*d_in_seconds

    M = a_cgs**3 / P_cgs**2 * 4 * np.pi**2 / G - m


    star = keplerian.Central(mass=M/M_sun, radius=1.0)

    system = keplerian.System(star).add_body(mass = m, period = P, 
                        eccentricity = e, omega_peri = omega/rad_in_deg,
                        inclination = i/rad_in_deg,time_transit=0.0
                        )

    planet = system.bodies[0]

    x,y,z = planet.position(phase * P) # in Rsol.
    return(x,y,z*-1) #Invert z so that negative directions are towards the observer.



from jax import jit
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



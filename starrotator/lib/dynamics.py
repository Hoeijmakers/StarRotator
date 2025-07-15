from jax import jit
from jax import numpy as jnp
from jaxoplanet.orbits import keplerian
from starrotator.lib.constants import rad_in_deg, R_sun, d_in_seconds, c_light

@jit
def pos_eccentric(phase, M = 1.0, m = 0.0, P = 365.0, e = 0.0, omega = 0.0,
                i=90.0):
    """
    This function calculates the radial velocity in km/s for a planet in 
    an elliptical orbit using jaxoplanet. 
    
    Input is provided in terms of the orbital phases at which the radial
    velocity is required, and the system parameters, including the 
    stellar mass. Stellar mass is assumed to be known to better precision 
    than the semi-major axis a. If this is not the case, you need to 
    proceed by calculating M from a, using Kepler III.

    If the mass of the companion is non-negligible, then its mass can
    also be set.

    The coordinate system follows that of the exoplanet package:
    https://docs.exoplanet.codes/en/latest/tutorials/data-and-models/


    Parameters
    ----------
    phase : float, array-like
        Orbital phase, typically between 0 and 1.0. 0.0 is mid-transit.
        Equivalent to the Mean Anomaly divided by 2 pi.

    M : float
        Stellar mass in solar masses.

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
    rv : float, array-like, same as phase
        The planet's radial velocity in km/s.

    """


    star = keplerian.Central(mass=M, radius=1.0)

    system = keplerian.System(star).add_body(mass = m, period = P, 
                        eccentricity = e, omega_peri = omega/rad_in_deg,
                        inclination = i/rad_in_deg,time_transit=0.0
                        )

    planet = system.bodies[0]

    # NEED TO CHANGE THIS TO PULL OUT POSITIONS IN X AND Y INSTEAD.
    vz = planet.velocity(phase * P)[2] # along the z axis in Rsol per day.

    # Convert to km/s, note that the positive z axis is the negative RV axis.
    rv = vz * R_sun / 1e5 / d_in_seconds * -1 
   
    return(xp,yp)




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
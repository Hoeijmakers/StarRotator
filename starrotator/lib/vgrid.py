import jax
import jax.numpy as jnp
from jax import jit
import numpy as np
import starrotator.lib.operations as ops
from functools import partial

# The routines in this file calculate the velocity and flux grids of a 2D pixellated stellar disk, 
# meant for summing to create the broadened stellar spectrum.

@jit
def calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate):
    """
    This calculates the velocity grid of a rotating star.
    based on Cegla+2016. See Fig. 3 and equations 2 - 8 of [1]_.

    It takes the stellar parameters and then calculates the stellar velocities on a pixellated stellar disk.

    Parameters
    ----------
        x : array-like
            1D numpy arrays to create the stellar grid in units of stellar radius in the x-direction.

        y : array-like
            1D numpy arrays to create the stellar grid in units of stellar radius in the y-direction.

        i_stellar : float
            Inclination of the stellar spin axis in degrees.

        vel_eq : float
            Equatorial velocity of the star in km/s.
            
        diff_rot_rate : float
            Differential rotation parameter alpha as defined in Cegla et al. 2016.
            If you wish to use a non-differentially rotating star, you'll be better off calculating 
            the stellar spectrum without using this fully pixellated model (integration v1).

    Returns
    -------
        vel_stellar_grid : array-like
            A 2D, pixellated model of the stellar disk mapping the radial velocity as a function of 
            projected x and y position, in km/s.

    References
    ----------
    .. [1] Cegla, Lovis, et al. (2016), A&A, 588, A127
        (https://arxiv.org/pdf/1602.00322.pdf)        
    """
    z,x_full,y_full = calc_z(x,y)#This is jitted. x_full and y_full are simple coordinate arrays.

    # Doing equation 8 from Cegla 2016 would go like this:
    # vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*...
    # vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*...
    # (x_full*np.sin(alpha)-y_full*np.cos(alpha))))

     # However, we align the star to the coordinate system, and so the projected obliquity is not incorporated here. 
     # The tilt of the star compared to the planet has to be incorporated in the planet's path. So we do one rotation
     # fewer than Cegla 2016 and we don't align the x-y grid to the orbital plane.
    i_rad = jnp.radians(i_stellar)
    beta = jnp.pi/2 - i_rad # In the paragraph preceding eq 6 in C2016.
    y_per_prime = (z*jnp.sin(beta)+jnp.cos(beta)*y_full)# Equation 7 in C2016.
    vel_stellar_grid = x_full*vel_eq*jnp.sin(i_rad) *  (1.0 - diff_rot_rate * y_per_prime**2) # Equation 8 in C2016.
    return(vel_stellar_grid)



def calc_vel_stellar_old(x,y,i_stellar, vel_eq, diff_rot_rate):
    """
    Depricated version of the above. Note that there was a missing square in the line that activated drr!!
    """
    # Convert angles to rad
    # alpha = np.radians(proj_obliquity)
    beta = np.radians(i_stellar)
    #pre calculate matrices
    z,x_full,y_full = calc_z(x,y)#This is jitted.

    # Doing equation 8 from Cegla+2016 would go like this:
    # vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*(x_full*np.sin(alpha)-y_full*np.cos(alpha))))

     # However, we align the star to the coordinate system, and so the projected obliquity is not incorporated here. 
     # The tilt of the star compared to the planet has to be incorporated in the planet's path. So we do one rotation
     # fewer than Cegla 2016 and we don't align the x-y grid to the orbital plane.
    vel_stellar_grid = x_full*vel_eq*np.sin(beta) * (1.0 - diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*y_full))
    return(vel_stellar_grid)




@jit
def calc_z(x,y):
    """
    This calculates the z-coordinate of the star as used in [1]_. See Fig. 3 and equations 2 - 8.
    It takes the 1D x and y grids as input.

    Parameters
    ----------
        x : array
            1D numpy arrays to create the stellar grid in units of stellar radius in the x-direction.

        y : array
            1D numpy arrays to create the stellar grid in units of stellar radius in the y-direction.

        i_stellar : float
            Inclination of the stellar spin axis in degrees.

        vel_eq : float
            Equatorial velocity of the star in km/s.
            
        diff_rot_rate : float
            Differential rotation parameter as defined in Cegla et al. 2016.
            If you wish to use a non-differentially rotating star, you'll be better off calculating 
            the stellar spectrum without using this fully pixellated model.

    Returns
    -------
        d : array-like
            A 2D, pixellated model of the stellar disk mapping the z coordinate as used in [1]_, in
            units of x and y.

        x_full : array-like
            The x-coordinate broadcasted along the entire 2D grid.

        y_full : array-like
            The y-coordinate broadcasted along the entire 2D grid.

    References
    ----------
    .. [1] Cegla, Lovis, et al. (2016), A&A, 588, A127
        (https://arxiv.org/pdf/1602.00322.pdf)  
    """
    x_full,y_full = x[None, :],y[:, None]
    d = 1 - x_full**2 - y_full**2
    d_clipped = jnp.where(d < 0,jnp.nan,d)
    return(jnp.sqrt(d_clipped),x_full,y_full)


@partial(jax.jit,static_argnames=["norm"])
def calc_flux_stellar(x,y,u1,u2,norm=True):
    """
    This calculates a flux map of the stellar disk assuming quadratic limb darkening.

    Parameters
    ----------
        x : array
            1D numpy arrays to create the stellar grid in units of stellar radius in the x-direction.

        y : array
            1D numpy arrays to create the stellar grid in units of stellar radius in the y-direction.

        u1 : float
            Linear limb darkening term.

        u2 : float
            Quadratic limb darkening term.

        norm : bool
            If set to true, the flux map is divided by its sum, normalising the sum to 1.0
            This is meant to be physically relevant in case the full disk is considered; but may not
            always be desired (e.g. when comparing absolute flux values of parts of the disk). 

    Returns
    -------
        flux_grid : array-like
            A 2D, pixellated model of the stellar disk mapping its relative flux. If norm = True,
            this sums to 1.0.
    """
    # I timed this to be 1000 times faster than the older for-loop way. Thanks jax.
    # I also confirmed that this mathematically is identical to the old one.


    # The limb darkening function below is jitted. And x,y are passed such that they broadcast
    # into a 2D array.
    d = jnp.sqrt(x[None, :]**2 + y[:,None]**2)
    d_clipped = jnp.where(d > 1,jnp.nan,d)
    flux_grid = ops.limb_darkening(d_clipped,u1,u2)
    if norm:
        return(flux_grid/jnp.nansum(flux_grid))
    else:
        return(flux_grid)


def calc_flux_stellar_old(x,y,u1,u2):
    import numpy as np
    import starrotator.lib.operations as ops


    z,x_full,y_full = calc_z(x,y)
    z = np.array(z)
    flux_grid = z*0.0
    for i in range(len(x)):
        for j in range(len(y)):
            if np.sqrt(x[i]**2+y[j]**2) <= 1.0:
                flux_grid[j,i]=ops.limb_darkening((np.sqrt(x[i]**2+y[j]**2)),u1,u2)*(z[j,i]*0.0+1.0)#Last multiplication is to put NaNs back into place.
    flux_grid /= np.nansum(flux_grid)
    return(flux_grid)




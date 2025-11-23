import jax
import numpy as np
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.dynamics import doppler_shift, doppler_shift_dlogl
from starrotator.lib.operations import vert_int_q_ld, circ_int_q_ld, vert_int_q_ld_bounded
from starrotator.lib.vgrid import calc_vel_stellar
from starrotator.lib.vgrid import calc_flux_stellar





@partial(jax.jit,static_argnames=["N","constant_dlogl"])
def sum_stellar_spectrum_v1(wl,fx,vel_eq,i_stellar,a1,a2,N=400,constant_dlogl=False):
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

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees. 0 is pole-on.

        a1 : float
            Linear limb darkening coefficient.

        a2 : float
            Quadratic limb darkening coefficient.

        N : int
            The size of the grid over which the integration is carried out. Larger values produce more accurate results.

        constant_dlogl: bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, not an array. 

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
            Note that the flux is multiplied by the differential dx, so that the integral is 
            insensitive to N.
    """
    x = jnp.linspace(-1,1,N)
    dx = x[1]-x[0]
    v_axis = x*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.
    weights = vert_int_q_ld(x,a1,a2)

    if constant_dlogl:
        fx_shifted = doppler_shift_dlogl(wl,fx,v_axis).T*weights
    else:
        fx_shifted = doppler_shift(wl,fx,v_axis).T*weights

    total = circ_int_q_ld(a1,a2) #The analytical integral
    F_out = jnp.nansum(fx_shifted,axis=1)*2*dx#Factor of 2 because this is calculating a semi-circle.
    return(F_out/total)




@partial(jax.jit,static_argnames=["N"])
def sum_stellar_spectrum_v1_mu(wl,fx_array,vel_eq,i_stellar,mu,N=400):
    """This builds the stellar spectrum by doppler-shifting and interpolating an array of input
        spectra, corresponding to an array of mu angles.
        This version of the integration assumes a SIMPLE VELOCITY GRID (that means: without drr)
        but a MU-VARYING stellar spectrum. The computation therefore implicitly includes 
        center-to-limb variation, but also limb darkening: The continuum of the mu-dependent
        spectra is responsible for this. This version of the integration is fastest because
        it uses the axial symmetry of the disk, so it only computes the velocity of the disk along
        the projected equator, and weighs the mu-dependent spectra by their number of occurrences
        in each column of the disk. 

        This is the fastest possible approximation.

        Parameters
        ----------
        wl : array-like
            One-dimensional array of the stellar model wavelengths, corresponding to the flux values in fx_array.

        fx_array : array-like
            Two-dimensional array of the stellar model fluxes corresponding to wl, at different mu angles.
            If you intend to pass only one mu angle for some reason, then it should be passed explicitly with a leading dimension.

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees. 0 is pole-on.

        mu : array-like
            The array of mu-values corresponding to the array of spectra.


        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
            Note that the flux is multiplied by the differential dx, so that the integral is 
            insensitive to N.

        Notes
        -----
        In this function, the input spectra are broadcast to a shape of n_wl times N times n_mu_angles.
        If you have 1e5 wavelength points, 400 gridpoints and 10 mu angles, that's 3.6 GB.
        So do not choose too many mu angles or many wavelength points at the same time, or you will 
        overflow your RAM, and there is no explicit guard against that.
        If a real-life application occurs where you're overflowing your RAM, you could always batch the
        computation in wavelength manually.

    """
    x = jnp.linspace(-1,1,N)
    dx = x[1]-x[0]
    v_axis = x*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis.


    doppler_shift_batch = jax.vmap(doppler_shift, in_axes=(None, 0, None)) # Here is some of the magic.
    #This creates a function doppler_shift_batch that is a vectorised execution of the doppler_shift
    #function, looping over axis 0 of the second input, which is fx_array. So each row of fx_array
    #gets doppler shifted by v_axis.

    fx_array_shifted = doppler_shift_batch(wl,fx_array,v_axis)
    #So the output will have dimensions len(v) * len(fx_array) * len(wl) whereas the shape of the output
    # of a single call to doppler_shift would have dimensions len(v)*len(wl)).

    #Note that if any jitted function is called repeatedly with different shapes for input, this will trigger jax recompilations. 
    #Also note that the output can be big: For 100,000 wavelength points, 200 velocities and 10 mu angles, this could be 1.6GB in memory.




    #Now follows the main logic of this computation.
    #We create a map of bins of mu across the stellar disk.
    #This is constructed for each column (in x) along the y-axis of the grid.
    #So for each location at x, we determine which mu angles exist from the equator
    #up to the disk edge. Each mu-bin that exists in a column has a value mu_i.
    #Each of these mu-bins is thus associated with a spectrum fx_i.
    #The weight at which each of these spectra is added, is the area of that piece of surface.
    #That is dx * ymax-ymin.
    #So that's all we need to do. Calculate for each fx_i in which bins they occur on the disk,
    #and then add them. The only reason this looks a bit convolved is because we're summing over a disk
    #using square arrays, well here we go.

    #We first need to bin up the mu array to find bin edges in mu-space.
    mu_centers = (mu[0:-1]+mu[1:])/2
    mumin = jnp.concatenate([jnp.array([0]), mu_centers])
    mumax = jnp.concatenate([mu_centers,     jnp.array([1])])


    #We then compute the arrays of y-edges associated with these mu edges.
    #This comes from mu^2 = x^2+y^2, leading to y = sqrt(mu^2 - x^2).
    ymin = jnp.sqrt(jnp.clip(mumin[None,:]**2 - x[:, None]**2,0,None))
    ymax = jnp.sqrt(jnp.clip(mumax[None,:]**2 - x[:, None]**2,0,None))
    #This is where the magic happens, because the output of the below is a 
    #2D array for the minimum and maximum edges: N_mu times N_x.
    #For each x, it tells you what the y-values are of the edges of the mu-bins.
    #This is easiest to understand for x=0 (center of the disk up):
    #As you travel up in from y=0 to y=1, you start at the first mu-angle, until you
    #reach the y-coordinate of the first mu angle. Then you transition into the second mu-value, etc.
    #until y=1 and you have reached the top of the disk.
    #This is also true for x>0, but now, the starting value is non-zero mu. Meaning mu-bins disappear
    #as you travel outward in the x-direction. At x close to 1, most mu values don't occur anymore in the column
    #above you. This is automatically handled by the clip statement: For a mu bin that is not in the column at x=/=0,
    #ymin and ymax are both zero. So the area of this piece will be (ymax-ymin)*dx = 0. 
    # This also works for partial bins: If ymin gets clipped but not ymax. Lovely!

    W = (ymax-ymin) *dx
    
    total = circ_int_q_ld(0.0,0.0) #The analytical integral of a non-limb-darkened star.
    # print(np.shape(fx_array_shifted))
    # print(np.shape(W.T))
    fx_weighted = W.T * np.transpose(fx_array_shifted,(2,0,1))
    fx_sum = jnp.sum(fx_weighted,axis=(1,2))

    return(fx_sum*2 / total)#Factor of 2 because this has been assuming a semi-circle.




@partial(jax.jit,static_argnames=['N','batched'])
def sum_stellar_spectrum_v2(wl,fx,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=200,batched=True):
# def sum_stellar_spectrum_v2(wl,fx,vel_grid,flux_grid,batched=True):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a FREE VELOCITY GRID (that means: with drr)
        but a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation that is 
        available, is modelled by a broad-band, quadratic limb darkening model.
    
        In this scenario we cannot assume that any velocity in the grid repeats itself when drr is present, and there is
        probably no symmetry. This means that we need to interpolate the spectrum for all velocities and
        sum over the whole grid. We can do this fully vectorised, or batched (row-by-row in the grid) to 
        not load this huge array into memory. Full vectorization may be a tiny bit faster but may only really be useful on
        big servers with sufficient RAM, so batched=True is the default.

        This is one step short of full brute-force integration (freedom in the velocity field and in the flux field), so for this
        function there exists no mu-version (no v2_mu). That would be v3.

        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths. One-dimensional.

        fx : array-like
            The stellar model fluxes corresponding to wl. One-dimensional.

        vel_eq : float
            The equatorial velocity of the star in km/s.

        i_stellar : float
            The inclination of the stellar spin axis in degrees. 0 is pole-on.

        a1 : float
            Linear limb darkening term.

        a2 : float
            Quadratic limb darkening term.

        diff_rot_rate : float
            Differential rotation parameter as defined in Cegla et al. 2016.
            If you wish to use a non-differentially rotating star, you'll be better off calculating 
            the stellar spectrum without using this fully pixellated model.

        batched : bool
            Compute the sum row-by-row or entirely at once (warning: very memory hungry).

        Returns
        -------
        fx : array
            The flux axis of the summed spectrum corresponding to the input wavelength points.
    """
    x = jnp.linspace(-1,1,N)
    dx = x[1]-x[0]
    flux_grid = calc_flux_stellar(x,x,a1,a2,norm=False)
    vel_grid  = calc_vel_stellar(x,x,i_stellar, vel_eq, diff_rot_rate)
    
    total = circ_int_q_ld(a1,a2) #The analytical integral
    F_out = sum_by_vel_and_flux_row(wl,fx,vel_grid,flux_grid,batched=True)
    return(F_out*dx**2 / total)

    

@partial(jax.jit,static_argnames=["wl"])
def sum_stellar_spectrum_v3(wl,fx_array,fx_map,v_map,weights=None):
    """
    Brute-force integration of a custom stellar disk.
    
    This takes a custom velocity map (could be generated with rigid rotation, drr or anything else) as well as an array of
    spectra to be tiled onto the disk. The gridpoint locations of each spectrum in that list are mapped using the fx_map parameter.
    The gridsize is specified by the sizes of v_map and fx_map.

    """






def sum_by_vel_and_flux_row(wl,fx,vel_grid,flux_grid,batched=True):
    """This chops up the problem of summing up to the integrated spectrum by 
    looping over rows of the disk's grid. This can be used with an arbitrary
    grid of fluxes (that scale the spectrum fx), and velocity grids (that doppler
    shift the spectrum fx). The only limitation is that wl,fx is constant over
    the stellar disk -- so this is one step short of a completely free 
    brute-force integration, and NOT compatible with mu-dependent spectra."""
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
    else: # If a LOT of memory is available, then we may not need to batch and do it directly.
        # This is very memory-hungry because the Nx n Ny x Nwl array gets done all at once.
        # You very likely want to set batched=True instead....
        fx_shifted = doppler_shift(wl,fx,jnp.ravel(vel_grid)).T*jnp.ravel(flux_grid)
        F_out = jnp.nansum(fx_shifted,axis=1)
    return(F_out)
    






# Note that this approximation may be very good even for the case of mu-dependent spectra without drr 
# because the range of mu-angles covered by the planet is small. Mu dependence then still comes in via
# the change in xp,yp but that can be addressed by either:
# -writing the function such that accepts an fx-array and a mu-array.
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
            you have many planet positions and high N will result in high memory load. In case of overflowing RAM, you can always
            chunck the computation in wavelength yourself.

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is a 2D array as xp,yp are 1D). This integrates to an arbitrary number, but this scaling corresponds
            to the scaling of the full disk. Normalisation is therefore achieved by dividing by the integral of the full
            disk, that depends only on a1 and a2. Note that as part of this process, the flux is multiplied by the
            differential dx*R, so that the integral is insensitive to N.
    """
    x = jnp.linspace(-1,1,N) #This defines the grid over which the integration takes place.
    #The strategy here is to calculate (integrate) the amount of surface of the stellar disk
    #in the vertical direction, for each x-gridpoint. This is done analytically for a
    #limb darkened disk, but with y-limits that are different from the simple 1=x*2+y*2 if you are calculating
    #this integral on the full disk (sum_stellar_spectrum_v1). Instead, we do the same thing but just
    #evaluating the definite integral.
    #Everything below is there to calculate the y-limits of the PLANET disk, or the part of the planet disk that
    #is in front of the stellar disk (ignoring the part that is outside of the disk at in/egress).
    dxR = (x[1]-x[0])*Rp #Our grid is not from -1,1, but from -Rp to Rp. So much finer grid points than the star, for the same N.
    x_array = x*Rp + xp[:,None]#Shifting to the location of the planet.
    y_max = yp[:,None] + jnp.sqrt(1-x**2)*Rp #This is the edge of the planet on the top.
    y_min = yp[:,None] - jnp.sqrt(1-x**2)*Rp #The lower edge of the planet.
    y_star_max = jnp.sqrt(1-x_array**2)*0.999999#The maximum y you can have is the edge of the stellar disk.
    y_star_min = -jnp.sqrt(1-x_array**2)*0.999999
    y_final_max = jnp.clip(y_max,y_star_min,y_star_max)#Which is smaller? The y-coord of the planet or the y-coord of the stellar disk?
    y_final_min = jnp.clip(y_min,y_star_min,y_star_max)#jnp.clip is so awesome.
    v_axes = x_array*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis. No drr so this is nicely uniform along y.
    weights = vert_int_q_ld_bounded(x_array,y_final_min,y_final_max,a1,a2) #This here is all the magic, the analytical integration.
    #Weights now has the same shape as x_array. x_array could be a 2D array, N_planet_phases x N_gridpoints.
    #For each planet phase, weights expresses the amount of surface area (from y_min to y_max) at each x-position (meaning: at each RV).
    #Thats why a dxR appears below, that is the differential.
    fx_shifted = doppler_shift(wl,fx,v_axes).T*weights # This is a very big array isn't it? For large N this is large too?
    F_out = jnp.nansum(fx_shifted,axis=2).T #Multiply with the differential for normalization.
    total = circ_int_q_ld(a1,a2)
    return(F_out*dxR/total)




@jit
def find_mu_interval(mu_array, mu_value):
    # Find the index of the right boundary
    l2 = jnp.searchsorted(mu_array, mu_value, side="right")
    # Left boundary is one index before
    l1 = jnp.maximum(l2 - 1, 0)
    # Clip right_idx to avoid going out of bounds
    l2 = jnp.minimum(l2, mu_array.shape[0] - 1)
    return(l1, l2)

@jit
def mu_interp_weight(mu_array,l1,l2,mu_value):
    # Get interval endpoints
    mu_left = mu_array[l1]
    mu_right = mu_array[l2]
    # Compute fractional distance, safe against divide-by-zero
    t = (mu_value - mu_left) / (mu_right - mu_left + 1e-12)
    t = jnp.clip(t, 0.0, 1.0)  # avoid going outside [0,1]  
    return(t)

@jit
def interp_fx_array(mu_array,fx_array,mu_value):
    """
    This linearly interpolates the disk-resolved spectrum
    given a value of mu.
    It handles edge cases (so if the target value is outside the range
    of the my values provided, it edge clips as it should).
    It also allows for broadcasting, so you can get many
    interpolated spectra returned. Note that changing the number of mu angles
    will cause a jax recompile (shape change), but this is not a likely
    scenario inside a retrieval because the number of mu angles is an
    initial user choice. I imagine that this function can be hot-rodded
    to interpolate over other types of things (e.g. temperature_array instead
    of mu). It doesn't matter for recompilation as long as the number of 
    interpolates is fixed throughout a loop/retrieval.
    """
    l1,l2 = find_mu_interval(mu_array,mu_value)
    t = mu_interp_weight(mu_array,l1,l2,mu_value)
    fx_interpolated = (fx_array[l1].T*(1-t) + t*fx_array[l2].T).T
    return(fx_interpolated)
    



@partial(jax.jit,static_argnames=['N','small_planet'])
def sum_hidden_spectrum_v1_mu(wl,fx_array,xp,yp,Rp,vel_eq,i_stellar,mu,N=100,small_planet=True):
    """This builds the stellar spectrum that is hidden behind the planet by doppler-shifting 
        and interpolating the input spectrum, summing only over the range of coordinates that is
        obscured by the planet disk.

        This version of the integration assumes a SIMPLE VELOCITY GRID (that means: without drr)
        and a MU-DEPENDENT (mu-resolved) stellar spectrum. Center-to-limb variation is thus
        modelled by the variation in the limb-resolved spectra. 

        This function integrates in a manner that is equivalent to sum_stellar_spectrum_v1_mu,
        unless the small-planet approximation is made, in which case it integrates like sum_hidden_spectrum_v1.


        Parameters
        ----------
        wl : array-like
            The stellar model wavelengths.

        fx_array : array-like
            Two-dimensional array of the stellar model fluxes corresponding to wl, at different mu angles.
            If you intend to pass only one mu angle for some reason, then it should be passed explicitly with a leading dimension.

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

        mu : array-like
            The array of mu-values corresponding to the array of spectra.
            

        N : int
            Number of grid-points in the x direction onto which the planet is simulated. Note that the resolution of the planet is independent
            of the resolution of the grid of the star. Note that as this calculation is batched over planet positions, cases where
            you have many planet positions and high N will result in high memory load if batched is set to false.

        small_planet : bool
            Make a small-planet approximation (true) or not (false). In the small-planet approximation, the mu-variation across the planetary disk is neglected, and a single mu
            angle (at planet center) is assumed. The spectrum is interpolated from the two nearest mu angles for which spectra are provided (so any planet motion does result in a 
            change in the spectrum).

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is a 2D array as xp,yp are 1D). Note that as part of this process, the flux is multiplied by the
            differential dx*R, so that the integral is insensitive to N.
    """

    if small_planet:
        x = jnp.linspace(-1,1,N)#This whole strategy is the same as sum_hidden_spectrum_v1 above.
        #We can even compute the weights in the same way. The only thing we need to do differently
        #is determine which fx_spectrum is associated with which which mu angle.
        dxR = (x[1]-x[0])*Rp
        x_array = x*Rp + xp[:,None]
        y_max = yp[:,None] + jnp.sqrt(1-x**2)*Rp
        y_min = yp[:,None] - jnp.sqrt(1-x**2)*Rp
        y_star_max = jnp.sqrt(1-x_array**2)*0.999999
        y_star_min = -jnp.sqrt(1-x_array**2)*0.999999
        y_final_max = jnp.clip(y_max,y_star_min,y_star_max)
        y_final_min = jnp.clip(y_min,y_star_min,y_star_max)
        v_axes = x_array*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis. Same shape as x_array.
        weights = vert_int_q_ld_bounded(x_array,y_final_min,y_final_max,0.0,0.0) #No limb darkening since we are 
        #dealing with mu-dependent spectra.

        # Continue here by calculating the mu's corresponding to the xp,yp positions.
        mup = jnp.sqrt(xp**2+yp**2)
        # Create an array of fx's that match that (by interpolating what the stellar spectrum)
        # is at the mu of each planet location.
        fx_p = interp_fx_array(mu,fx_array,mup) # These are the interpolated planet spectra belonging to each xp,yp.


        doppler_shift_batch = jax.vmap(doppler_shift, in_axes=(None, 0, 0)) # Here is some of the magic.
        # #This creates a function doppler_shift_batch that is a vectorised execution of the doppler_shift
        # #function, looping over axis 0 of the second input, which is fx_array. So each row of fx_array
        # #gets doppler shifted by v_axis.


        total = circ_int_q_ld(0.0,0.0)
        #Actual normalization depends on the normalization of the mu-spectra, though.
        #This is just designed to normalize in case that the input spectrum is the same for all mu (no limb darkening).
        fx_shifted = doppler_shift_batch(wl,fx_p,v_axes).T * weights.T
        F_out = jnp.nansum(fx_shifted,axis=1).T * dxR
        return(F_out/total)
    else:
        raise Exception("NOT-SMALL PLANET IS NOT IMPLEMENTED YET.")







#sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100)
# def sum_hidden_spectrum_v2(wl,fx,vel_grid_array,flux_grid_array,N=100,batched=True):
@partial(jax.jit,static_argnames=['N','batched'])
def sum_hidden_spectrum_v2(wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=100,batched=True):
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
            you have many planet positions and high N will result in high memory load. In case of overflowing RAM, you can always
            chunck the computation in wavelength yourself.

        batched : bool
            Compute the integration phase-by-phase (True) or all at once (warning: very memory hungry).


        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is generally a 2D array). Note that the result is unnormalised and therefore depends on the
            normalization of the flux grid array and of the input spectrum.
    """


    x = jnp.linspace(-1,1,N)
    total = circ_int_q_ld(a1,a2) #The analytical integral
    flux_grid_array,vel_grid_array,mu_array,dxR = create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=N)


    def scan_fn(carry, inputs):
        vel_i, flux_i = inputs  # Each is a planetary flux grid.
        sum_i = sum_by_vel_and_flux_row(wl,fx,vel_i,flux_i,batched=batched)
        carry += 1
        return carry, sum_i  # No output needs to be collected.

    _, F_out = lax.scan(scan_fn, 0, (vel_grid_array.T, flux_grid_array.T))

    return(F_out*dxR**2 / total)





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
    mu = jnp.sqrt(x[None,:]**2 + x[:,None]**2) #Yes, this is the same mu as pySME mu jsyk (cos(theta)).
    mask = jnp.where(mu > 1,jnp.nan,mu)*0.0+1.0

    flux_disk_array  = calc_flux_stellar(x_array.T,y_array.T,a1,a2,norm=False) #This is really, really awesome.
    vel_disk_array  = calc_vel_stellar(x_array.T,y_array.T,i_stellar, vel_eq, diff_rot_rate)


    flux_disk_array_masked = flux_disk_array*mask.T[:,:,None] #This crops out the circle that is the planet.
    vel_disk_array_masked = vel_disk_array*mask.T[:,:,None]
    mu_array_masked = mu*1.0
    return(flux_disk_array_masked,vel_disk_array_masked,mu_array_masked,dxR)













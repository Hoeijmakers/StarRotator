import jax
import numpy as np
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.dynamics import doppler_shift, doppler_shift_dlogl
from starrotator.lib.operations import vert_int_q_ld, circ_int_q_ld, vert_int_q_ld_bounded
from starrotator.lib.vgrid import calc_vel_stellar
from starrotator.lib.vgrid import calc_flux_stellar
import matplotlib.pyplot as plt




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

        constant_dlogl : bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 

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




@partial(jax.jit,static_argnames=["N",'constant_dlogl'])
def sum_stellar_spectrum_v1_mu(wl,fx_array,vel_eq,i_stellar,mu_array,N=400,constant_dlogl=False):
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

        mu_array : array-like
            The array of mu-values corresponding to the array of spectra. Mu is defined as the cosine of the angle between the
            viewing angle and the normal from the surface, meaning that it varies between 1 to 0 from the center to the edge.
            Following mu = cos(theta), mu**2 + r**2 = R = 1. So mu = sqrt(1-r**2).

        N : int
            The size of the grid over which the integration is carried out. Larger values produce more accurate results.

        constant_dlogl : bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 


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


    if constant_dlogl:
        doppler_shift_batch = jax.vmap(doppler_shift_dlogl, in_axes=(None, 0, None))
    else:
        doppler_shift_batch = jax.vmap(doppler_shift, in_axes=(None, 0, None))
    #This creates a function doppler_shift_batch that is a vectorised execution of the doppler_shift
    #function, looping over axis 0 of the second input, which is fx_array. So each row of fx_array
    #gets doppler shifted by v_axis.

    #doppler_shift_batch = jax.vmap(doppler_shift, in_axes=(None, 0, None)) # Here is some of the magic.

    fx_array_shifted = doppler_shift_batch(wl,jnp.flip(fx_array,axis=0),v_axis) #Flip this because mu-array is sorted the wrong way.

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

    #We first need to convert the mu array to values of r, and then bin it up to find bin edges in r-space:
    r_array = jnp.sqrt(1 - jnp.flip(mu_array)**2)# Array in radial coordinate corresponding to mu angles. 
    #Also flip mu array so this is ordered ascending, otherwise the below doesn't work.
    r_centers = (r_array[0:-1]+r_array[1:])/2
    rmin = jnp.concatenate([jnp.array([0]), r_centers])
    rmax = jnp.concatenate([r_centers,     jnp.array([1])])


    #We then compute the arrays of y-edges associated with these mu edges.
    #This comes from mu^2 = x^2+y^2, leading to y = sqrt(mu^2 - x^2).
    ymin = jnp.sqrt(jnp.clip(rmin[None,:]**2 - x[:, None]**2,0,None))
    ymax = jnp.sqrt(jnp.clip(rmax[None,:]**2 - x[:, None]**2,0,None))
    #This is where the magic happens, because the output of the below is a 
    #2D array for the minimum and maximum edges: N_mu times N_x.
    #For each x, it tells you what the y-values are of the edges of the r-bins.
    #This is easiest to understand for x=0 (center of the disk up):
    #As you travel up in from y=0 to y=1, you start in first r-bin, until you
    #reach the y-coordinate of the second r bin, etc.
    #until y=1 and you have reached the top of the disk.
    #This is also true for x>0, but now, the starting value is non-zero r. Meaning r-bins disappear
    #as you travel outward in the x-direction. At x close to 1, most r values don't occur anymore in the column
    #above you. This is automatically handled by the clip statement: For an r bin that is not in the column at x=/=0,
    #ymin and ymax are both zero. So the area of this piece will be (ymax-ymin)*dx = 0. 
    # This also works for partial bins: If ymin gets clipped but not ymax. Lovely!

    W = (ymax-ymin) *dx

    
    total = circ_int_q_ld(0.0,0.0) #The analytical integral of a non-limb-darkened star.
    # print(np.shape(fx_array_shifted))
    # print(np.shape(W.T))
    fx_weighted = W.T * np.transpose(fx_array_shifted,(2,0,1))
    fx_sum = jnp.sum(fx_weighted,axis=(1,2))

    return(fx_sum*2 / total)#Factor of 2 because this has been assuming a semi-circle.








@partial(jax.jit,static_argnames=['N','batched','constant_dlogl'])
def sum_stellar_spectrum_v2(wl,fx,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=200,batched=True,constant_dlogl=False):
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
            the stellar spectrum without using this fully pixellated model, meaing v1 or v1_mu.

        batched : bool
            Compute the sum row-by-row or entirely at once (warning: very memory hungry).

        constant_dlogl : bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 

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
    F_out = sum_by_vel_and_flux_row(wl,fx,vel_grid,flux_grid,batched=True,constant_dlogl=constant_dlogl)
    return(F_out*dx**2 / total)

    

@jit
def sum_stellar_spectrum_v3(wl,fx,vel_eq,i_stellar,a1,a2,N=200,batched=True,constant_dlogl=False):
    """
    Brute-force integration of a stellar disk with a single input spectrum, as well as spots.
    Spots are defined using a sequence of x-y positions, creating a spot-mask.
    
    This takes a custom velocity map (could be generated with rigid rotation, drr or anything else) as well as an array of
    spectra to be tiled onto the disk. The gridpoint locations of each spectrum in that list are mapped using the fx_map parameter.
    The gridsize is specified by the sizes of v_map and fx_map.

    """


@jit
def sum_stellar_spectrum_v4(wl,fx_array,fx_map,v_map,weights=None):
    """
    Brute-force integration of a custom stellar disk.
    
    This takes a custom velocity map (could be generated with rigid rotation, drr or anything else) as well as an array of
    spectra to be tiled onto the disk. The gridpoint locations of each spectrum in that list are mapped using the fx_map parameter.
    The gridsize is specified by the sizes of v_map and fx_map.

    """


#Why is this not done column-by-column? One RV position at a time?
#Well, I guess because this is meant for the DRR case. So RV is not constant in columns.
#But you can imagine situations where summing in columns is better.
@partial(jax.jit,static_argnames=["batched","constant_dlogl"])
def sum_by_vel_and_flux_row(wl,fx,vel_grid,flux_grid,batched=True,constant_dlogl=False):
    """
    Batched lax.scan implementation of full-grid summing. Under normal circumstances,
    you will want batched=True, unless you have got tons of memory available.
    
    This chops up the problem of summing up to the integrated spectrum by 
    looping over rows of the disk's grid. This can be used with an arbitrary
    grid of fluxes (that scale the spectrum fx), and velocity grids (that doppler
    shift the spectrum fx). The only limitation is that wl,fx is constant over
    the stellar disk -- so this is one step short of a completely free 
    brute-force integration, and NOT compatible with mu-dependent spectra.
    
    This choice of batching has one drawback, and that is that the user cannot control
    batch-size. But arrays of size Nwl x N should in realistic cases not be too large
    for memory anyway; unless a bizarre number of wavelength points is used, or unless
    you are running this on a raspberry pi.
    
    """
    if batched: #We use lax.scan to sum over the grid row-by-row to reduce memory load.
        def scan_fn(carry, inputs):
            vel_row, flux_row = inputs  # Each is a row of len(y)
            if constant_dlogl:
                fx_shifted = doppler_shift_dlogl(wl,fx,vel_row).T * flux_row
            else:
                fx_shifted = doppler_shift(wl,fx,vel_row).T * flux_row
            # fx_shifted = doppler_shift(wl, fx, vel_row).T * flux_row  # shape (len(wl), len(y))
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
        # You very likely want to set batched=True instead, unless Nwl is small.
        # fx_shifted = doppler_shift(wl,fx,jnp.ravel(vel_grid)).T*jnp.ravel(flux_grid)
        if constant_dlogl:
            fx_shifted = doppler_shift_dlogl(wl,fx,jnp.ravel(vel_grid)).T*jnp.ravel(flux_grid)
        else:
            fx_shifted = doppler_shift(wl,fx,jnp.ravel(vel_grid)).T*jnp.ravel(flux_grid)

        F_out = jnp.nansum(fx_shifted,axis=1)
    return(F_out)
    






# Note that this approximation may be very good even for the case of mu-dependent spectra without drr 
# because the range of mu-angles covered by the planet is small. Mu dependence then still comes in via
# the change in xp,yp but that can be addressed by either:
# -writing the function such that accepts an fx-array and a mu-array.
# -calling this function in a jitted loop or vmap over different fx-es, corresponding to the mu's the planet traverses.
@partial(jax.jit,static_argnames=['N',"constant_dlogl"])
def sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100,constant_dlogl=False):
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

        constant_dlogl : bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 

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
    if constant_dlogl:
        fx_shifted = doppler_shift_dlogl(wl,fx,v_axes).T*weights
    else:
        fx_shifted = doppler_shift(wl,fx,v_axes).T*weights

    #fx_shifted = doppler_shift(wl,fx,v_axes).T*weights # This is a very big array isn't it? For large N this is large too?
    F_out = jnp.nansum(fx_shifted,axis=2).T 
    total = circ_int_q_ld(a1,a2)
    return(F_out*dxR/total)#Multiply with the differential for normalization.




@jit
def find_X_interval(X_array, X_value):
    # Find the index of the right boundary
    l2 = jnp.searchsorted(X_array, X_value, side="right")
    # Left boundary is one index before
    l1 = jnp.maximum(l2 - 1, 0)
    # Clip right_idx to avoid going out of bounds
    l2 = jnp.minimum(l2, X_array.shape[0] - 1)
    return(l1, l2)

@jit
def X_interp_weight(X_array,l1,l2,X_value):
    # Get interval endpoints
    X_left = X_array[l1]
    X_right = X_array[l2]
    # Compute fractional distance, safe against divide-by-zero
    t = (X_value - X_left) / (X_right - X_left + 1e-12)
    t = jnp.clip(t, 0.0, 1.0)  # avoid going outside [0,1]  
    return(t)

@jit
def interp_fx_array(r_array,fx_array,r_value):
    """
    This linearly interpolates the disk-resolved spectrum
    given a value of r.
    It handles edge cases (so if the target value is outside the range
    of the r values provided, it edge clips as it should).
    It also allows for broadcasting, so you can get many
    interpolated spectra returned. Note that changing the number of r values
    will cause a jax recompile (shape change), but this is not a likely
    scenario inside a retrieval because the number of mu angles is an
    initial user choice. I imagine that this function can be generalised
    to interpolate over other types of things (e.g. temperature_array instead
    of r). It doesn't matter for recompilation as long as the number of 
    interpolates is fixed throughout a loop/retrieval.
    
    Note that r_array needs to be sorted in ascending order!
    """
    l1,l2 = find_X_interval(r_array,r_value)
    t = X_interp_weight(r_array,l1,l2,r_value)
    fx_interpolated = (fx_array[l1].T*(1-t) + t*fx_array[l2].T).T
    return(fx_interpolated)
    



@partial(jax.jit,static_argnames=['N','small_planet','constant_dlogl'])
def sum_hidden_spectrum_v1_mu(wl,fx_array,xp,yp,Rp,vel_eq,i_stellar,mu_array,N=100,small_planet=True,constant_dlogl=False):
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
            angle (at planet center) is assumed. The spectrum is interpolated from the two nearest mu angles for which spectra are provided (so any amount of planet motion does result in a 
            change in the spectrum, and there is no chopping artefact that would result from a nearest-neighbour approach). This approximation effectively neglects the change in limb darkening
            across the planetary disk, since this is assumed to be included in the mu-dependent description of fx_array. Note however that the small planet approximation does not
            neglect the RV-variation across the planet disk.

        constant_dlogl : bool
            Whether or not to use fast doppler shifting by virtue of a constant log-lambda wavelength array. Setting this to True also requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 

        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet position
            (so this is a 2D array as xp,yp are 1D, as N_xp times N_wl). Note that as part of this process, the flux is multiplied by the
            differential dx*R, so that the value of the integral is insensitive to N.
    """

    if small_planet:
        x = jnp.linspace(-1,1,N)#This whole strategy is the same as sum_hidden_spectrum_v1 above.
        #We can even compute the weights in the same way. The only thing we need to do differently
        #is determine which fx_spectrum is associated with which which mu angle, which is not done above
        # because each fx_spectrum is the same in sum_hidden_spectrum_v1.
        # In the small planet approximation we adopt the approximation that at each planet position, the 
        # planet covers effectively one mu-value. The spectrum associated with this mu value is determined
        # by interpolating the fx_array of mu-dependent spectra, so there is no nearest-neighbour chopping.
        # Effectively, this assumes that the mu-dependent variation of the spectrum is negligible over
        # the small area of the planet. However, we do not ignore the fact that the planet covers a range
        # of radial velocities, even if small. Otherwise, a broadening effect would be missing.
        dxR = (x[1]-x[0])*Rp
        x_array = x*Rp + xp[:,None]
        y_max = yp[:,None] + jnp.sqrt(1-x**2)*Rp
        y_min = yp[:,None] - jnp.sqrt(1-x**2)*Rp
        y_star_max = jnp.sqrt(1-x_array**2)*0.999999
        y_star_min = -jnp.sqrt(1-x_array**2)*0.999999
        y_final_max = jnp.clip(y_max,y_star_min,y_star_max)
        y_final_min = jnp.clip(y_min,y_star_min,y_star_max)
        v_axes = x_array*vel_eq*jnp.sin(jnp.radians(i_stellar)) # velocities along the y=0 axis. Same shape as x_array.
        weights = vert_int_q_ld_bounded(x_array,y_final_min,y_final_max,0.0,0.0) #Flat weights (no limb darkening) 
        # since we are dealing with mu-dependent spectra.

        # Continue here by calculating the mu's corresponding to the xp,yp positions.
        # r_p = jnp.sqrt(xp**2+yp**2)
        mu_p = jnp.sqrt(1 - xp**2-yp**2)
        # r_array = jnp.sqrt(1 - mu_array**2)# Array in radial coordinate corresponding to mu angles.
        # Create an array of fx's that match that (by interpolating what the stellar spectrum is
        # at that position). Note the coordinate transformation between r and mu, as r is a more
        # logical coordinate in this x-y space we're working in.

        fx_p = interp_fx_array(mu_array,fx_array,mu_p) # These are the interpolated planet spectra belonging to each xp,yp.


        # Here follows a whole lot of debugging plots to check that interpolation works.
        # plt.figure(figsize=(10,6))
        # plt.plot(mu_array,r_array)
        # plt.show()
        # print('Inferred radial values:',r_p)

        # print('Mu and r pairs:')
        # plt.figure(figsize=(10,6))
        # for i in range(len(fx_array)):
        #     print(mu_array[i],r_array[i])
        #     muprint = np.round(mu_array[i],3)
        #     rprint = np.round(r_array[i],3)
        #     plt.plot(wl,fx_array[i],label=f'mu = {muprint:.3f}, r = {rprint:.3f}')
        # plt.legend(frameon=False)
        # plt.title('Mu-dependent spectra (input)')
        # plt.show()

        # plt.figure(figsize=(10,6))
        # for i in range(len(r_p)):
        #     plt.plot(wl,fx_p[i],label=f'r = {r_p[i]:.3f} mu = {mu_p[i]:.3f}')
        # plt.legend(frameon=False)
        # plt.title('Interpolated spectra')
        # plt.show()



        if constant_dlogl:
            doppler_shift_batch = jax.vmap(doppler_shift_dlogl, in_axes=(None, 0, 0))
        else:
            doppler_shift_batch = jax.vmap(doppler_shift, in_axes=(None, 0, 0))
        # This creates a function doppler_shift_batch that is a vectorised execution of the doppler_shift
        # function. Looping is done over axis 0 of the second and third inputs, which are fx_array and v_axes.
        # fx_array is an array of spectra corresponding to each x-y position on the disk.
        # v_axes is an array of N v values corresponding to each x-y position on the disk.
        # So this formulation of vmap with in_axes=(None, 0, 0) loops over each x-y position, shifting
        # each spectrum in fx_array to its corresponding array of v values for that position x,y.


        total = circ_int_q_ld(0.0,0.0)
        #Actual normalization depends on the normalization of the mu-spectra, though.
        #This is just designed to normalize in case that the input spectrum is the same for all mu (no limb darkening).

        # Again some debugging plots to check that batched doppler shifting works.
        # fx_shifted = doppler_shift_batch(wl,fx_p,v_axes)
        # i=1
        # fig,ax = plt.subplots(1,2,figsize=(13,5))
        # ax[0].plot(x_array[i],v_axes[i])
        # ax[1].plot(wl,fx_shifted[i,0])
        # ax[1].plot(wl,fx_shifted[i,99])
        # ax[1].plot(wl,fx_shifted[i,199])
        # ax[1].set_xlim(502.2,502.8)
        # plt.show()

        # i=3
        # fig,ax = plt.subplots(1,2,figsize=(13,5))
        # ax[0].plot(x_array[i],v_axes[i])
        # ax[1].plot(wl,fx_shifted[i,0])
        # ax[1].plot(wl,fx_shifted[i,99])
        # ax[1].plot(wl,fx_shifted[i,199])
        # ax[1].set_xlim(502.2,502.8)
        # plt.show()

        # i=5
        # fig,ax = plt.subplots(1,2,figsize=(13,5))
        # ax[0].plot(x_array[i],v_axes[i])
        # ax[1].plot(wl,fx_shifted[i,0])
        # ax[1].plot(wl,fx_shifted[i,99])
        # ax[1].plot(wl,fx_shifted[i,199])
        # ax[1].set_xlim(502.2,502.8)
        # plt.show()

        # fig,ax = plt.subplots(1,4,figsize=(13,5))
        # ax[0].plot(x_array[0],weights[0])
        # ax[1].plot(x_array[1],weights[1])
        # ax[2].plot(x_array[4],weights[4])
        # ax[3].plot(x_array[5],weights[5])
        # plt.plot()

        # print('fx_p shape',fx_p.shape)
        # print('v_axes shape',v_axes.shape)
        # print('fx_shifted shape',fx_shifted.shape)
        # print('Weights',np.shape(weights))
        fx_shifted_transposed = np.transpose(doppler_shift_batch(wl,fx_p,v_axes),(2,0,1))
        # print('fx_shifted transposed',fx_shifted_transposed.shape)
        fx_shifted_weighted = fx_shifted_transposed * weights
        # print('fx_shifted_weighted',fx_shifted_weighted.shape)


        F_out = jnp.nansum(fx_shifted_weighted,axis=2).T * dxR
        # print('fx_array_out',F_out.shape)
        return(F_out/total)
    else:
        raise Exception("NOT-SMALL PLANET IS NOT IMPLEMENTED YET.")





#sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=100)
# def sum_hidden_spectrum_v2(wl,fx,vel_grid_array,flux_grid_array,N=100,batched=True):
# @partial(jax.jit,static_argnames=['N','batched','constant_dlogl','spot'])
def sum_hidden_spectrum_v2(wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=100,batched=True,constant_dlogl=False,spot=False):
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

        constant_dlogl : bool
            Whether or not to use fast doppler shifting. This requires that the wavelength axis is a float, set to the
            constant value of dloglambda (and not an explicit wavelength array). 

        spot : bool
            Under normal circumstances, this will calculate the spectrum behind a circularly projected area on the 
            stellar disk. With spot = True, it uses instead projected circular areas on the star, i.e. spots.


        Returns
        -------
        F : array
            The flux axis of the summed spectrum corresponding to the input wavelength points, for each planet/spot position
            (so this is generally a 2D array). Note that the result is unnormalised and therefore depends on the
            normalization of the flux grid array and of the input spectrum.
    """


    x = jnp.linspace(-1,1,N)
    total = circ_int_q_ld(a1,a2) #The analytical integral

    if spot:
        flux_grid_array,vel_grid_array,mu_array,dxR = create_spot_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=N)
    else:
        flux_grid_array,vel_grid_array,mu_array,dxR = create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=N)


    def scan_fn(carry, inputs):
        vel_i, flux_i = inputs  # Each is a planetary flux grid.
        sum_i = sum_by_vel_and_flux_row(wl,fx,vel_i,flux_i,batched=batched,constant_dlogl=constant_dlogl)
        carry += 1
        return carry, sum_i  # No output needs to be collected.

    _, F_out = lax.scan(scan_fn, 0, (vel_grid_array.T, flux_grid_array.T))

    return((F_out.T*dxR**2 / total).T) #Transposing-detransposing because of the possibility that dxR is a vector (if spot True).





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

        Rp : float
            The radius of the planet, in stellar radii.

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
    r_squared = x[None,:]**2 + x[:,None]**2 
    mu_squared = 1 - r_squared #Yes, this is the same mu as pySME mu jsyk (cos(theta)).
    mask = jnp.where(r_squared > 1,jnp.nan,r_squared)*0.0+1.0

    flux_disk_array  = calc_flux_stellar(x_array.T,y_array.T,a1,a2,norm=False) #This is really, really awesome (the fact that it's all vectorised).
    vel_disk_array  = calc_vel_stellar(x_array.T,y_array.T,i_stellar, vel_eq, diff_rot_rate)


    flux_disk_array_masked = flux_disk_array*mask.T[:,:,None] #This crops out the circle that is the planet. 
    #Note that even though x and x_array do not measure the same thing (x is the star, x_array the planet),
    #this can be applied to the small region that is the planet because x (-1,1) has the same shape as x_array (-1,1)*Rp. 
    vel_disk_array_masked = vel_disk_array*mask.T[:,:,None]  
    #Bit hacky but it is robust.

    mu_array_masked = jnp.sqrt(mu_squared*mask)#Calculation of the mu array is inaccurate. Needs to do x_array**2+y_array**2 or something.
    #Implement this when needed. Until then, a jnp.nan is returned.
    return(flux_disk_array_masked,vel_disk_array_masked,mu_array_masked*jnp.nan,dxR)





@partial(jax.jit,static_argnames=['N','margin'])
def create_spot_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=50,margin=0.15):
    """This builds arrays of the flux grid and velocity grid of the stellar disk that are
       part of a spot. This is modeled after create_hidden_grid_array above.

       This is the simplest implementation of a spotty stellar disk, whereby the
       spectrum of the star gets replaced by that of a spot -- conceptually nearly
       identical to replacing the spectrum of the star with 0, in the case of when the
       hidden spectrum (see above function) is calculated.

        This version creates a FREE VELOCITY GRID (that means: with possible drr), and the 
        possibility to add center-to-limb variation is modelled by broad-band limb darkening. 

        The output are 3D arrays of flux and velocity grids, one for each spot,
        with areas that are outside the spot or stellar disk set to NaN. These are designed
        to be passed into sum_hidden_spectrum_v2.

    

        Note that the coordinates should be chosen such that spots do not overlap.
        If they overlap, over-subtraction will occur. This code does not check for this.

        Also note that there is no way to account for the planet overlapping the spots.
        So this is a model for unocculted spots only.

        Parameters
        ----------
        xp : array-like
            The horizontal location of the spots, either as a float or as a 1D array.

        yp : array-like
            The vertical location of the spots, either as a float or as a 1D array.

        Rp : array-like
            The radii of the spots, in stellar radii.

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

        margin : float
            For projected spots, the radius of the window needs to be slightly bigger than the radius of the spot.


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
    x = jnp.linspace(-1,1,N)*(1+margin)
    dx = x[1]-x[0]
    dxR = dx*Rp # Rescaling the size of the differential. This is an array of length S = N_spots

    zp = jnp.sqrt(1-xp**2-yp**2)

    X, Y = jnp.meshgrid(x, x)
    
    X_array = X*Rp[:,None,None]+xp[:,None,None]
    Y_array = Y*Rp[:,None,None]+yp[:,None,None]
    r2 = X_array**2 + Y_array**2
    Z_array = jnp.sqrt(jnp.clip(1 - r2, 0, 1))

    n_s = jnp.array([xp, yp, zp]).T #Unit vector in the direction normal to the spot.
    # stack grid normals
    grid = jnp.stack([X_array, Y_array, Z_array], axis=-1)   # (S, Ny, Nx, 3)
    # align spot normals to grid shape
    n_s = n_s[:, None, None, :]              # (S, 1, 1, 3)
    # dot product. This here is the magic, realising that in a spot, the unit vector normal to the surface should
    # be within a range of angles of the spot center, and the x-y plane maps onto surface normals (that is mu),
    # and then the fact that this "within range of angles" is effected by demanding that a dot product is greater
    # than a limit cos(alpha).
    dot_product = jnp.sum(grid * n_s, axis=-1)      # (S, Ny, Nx)
    # Note that with R equal to the spot radius, if we assume that to be the great-circle distance, then R = alpha.
    # Otherwise, alpha = jnp.asin(Rp). I use the latter because I think its more observationally meaningful.
    #Since we want to know cos(alpha), and alpha=arcsin(R), we get cos(alpha) = sqrt(1-R**2)
    cos_alpha = jnp.sqrt(1-Rp**2) #This is length S.
    cos_alpha = cos_alpha[:, None, None]


    mask = (r2 <= 1.0) & (dot_product >= cos_alpha) #To be inside a spot and inside the star at the same time.
    x_array,y_array = X_array[:,0,:],Y_array[:,:,0]
    flux_disk_array  = calc_flux_stellar(x_array.T,y_array.T,a1,a2,norm=False) #This is really, really awesome (the fact that it's all vectorised).
    vel_disk_array  = calc_vel_stellar(x_array.T,y_array.T,i_stellar, vel_eq, diff_rot_rate)
    

    flux_disk_array_masked =flux_disk_array * (mask.transpose(1,2,0)) #This crops out the circle that is the planet. 
    #Note that even though x and x_array do not measure the same thing (x is the star, x_array the planet),
    #this can be applied to the small region that is the planet because x (-1,1) has the same shape as x_array (-1,1)*Rp. 
    vel_disk_array_masked = vel_disk_array * (mask.transpose(1,2,0)) 
    #Bit hacky but it is robust.

    mu_array_masked = (1-r2) * mask  #Calculation of the mu array is inaccurate. Needs to do x_array**2+y_array**2 or something.
    #Implement this when needed. Until then, a jnp.nan is returned.
    return(flux_disk_array_masked,vel_disk_array_masked,mu_array_masked*jnp.nan,dxR)










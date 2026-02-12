#This file contains a collection of functions that operate on spectra.
def convolve(array,kernel,edge_degree=1,fit_width=2):
    """np.convolve with accounting for edge effects via extrapolation of the edge of the array using a polynomial fit
    to the edge elements. By default, I fit over a range that is twice the length of the kernel; but
    this value can be modified using the fit_width parameter.

    Parameters
    ----------
    array : list, np.ndarray
        The horizontal axis.

    kernel : list, np.ndarray
        The convolution kernel. It is required to have a length that is less than 25% of the size of the array.

    edge_degree : int
        The polynomial degree by which the array is extrapolated in order to

    fit_width : int
        The length of the area at the edges of array used to fit the polynomial, in units of the length of the kernel.
        Increase this number for small kernels or noisy arrays.
    Returns
    -------
    array_convolved : np.array
        The input array convolved with the kernel

    Example
    -------
    >>> import numpy as np
    >>> a=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    >>> b=[-0.5,0,0.5]
    >>> c=convolve(a,b,edge_degree=1)
    """

    import numpy as np

    array = np.array(array)
    kernel= np.array(kernel)

    if len(kernel) >= len(array)/4.0:
        raise Exception(f"Error in ops.convolve(): Kernel length is larger than a quarter of the array ({len(kernel)}, {len(array)}). Can't extrapolate over that length. And you probably don't want to be doing a convolution like that, anyway.")

    if len(kernel) % 2 != 1:
        raise Exception('Error in ops.convolve(): Kernel needs to have an odd number of elements.')

    #Perform polynomial fits at the edges.
    x=findgen(len(array))
    fit_left=np.polyfit(x[0:len(kernel)*2],array[0:len(kernel)*2],edge_degree)
    fit_right=np.polyfit(x[-2*len(kernel)-1:-1],array[-2*len(kernel)-1:-1],edge_degree)

    #Pad both the x-grid (onto which the polynomial is defined)
    #and the data array.
    pad=findgen(int((len(kernel)-1)/2))
    left_pad=pad-(len(kernel)-1)/2
    right_pad=np.max(x)+pad+1
    left_array_pad=np.polyval(fit_left,left_pad)
    right_array_pad=np.polyval(fit_right,right_pad)

    #Perform the padding.
    x_padded = np.append(left_pad , x)
    x_padded = np.append(x_padded , right_pad) #Pad the array with the missing elements of the kernel at the edge.
    array_padded = np.append(left_array_pad,array)
    array_padded = np.append(array_padded,right_array_pad)

    #Reverse the kernel because np.convol does that automatically.
    #(Imagine doing a derivative with a kernel [-1,0,1] and it gets reversed...)
    kr = kernel[::-1]
    #The valid keyword effectively undoes the padding, leaving only those values for which the kernel was entirely in the padded array.
    #This thus again has length equal to len(array).
    return np.convolve(array_padded,kr,'valid')



def airtovac(wlnm):
    """Convert air to vaccuum wavelengths.

        Parameters
        ----------
        wlnm : float, np.array()
            The wavelength(s) in nm.

        Returns
        -------
        wl: float, np.array()
            The wavelengths.
    """
    import lib.test as test
    import numpy
    #First run standard tests on the input
    test.typetest(wlnm,[int,float,numpy.ndarray],varname='wlnm in airtovac')
    test.notnegativetest(wlnm,varname='wlnm in airtovac')
    test.nantest(wlnm,varname='wlnm in airtovac')
    #Would still need to test that wl is in a physical range.
    wlA=wlnm*10.0
    s = 1e4 / wlA
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    return(wlA*n/10.0)



def clip_spectrum(wl,fx,wlmin,wlmax,pad=0):
    """This crops a spectrum wl,fx to wlmin and wlmax.
    If pad !=0, it should be set to a positive velocity
    by which the output spectrum will be padded around wlmin,wlmax,
    which is returned as well as the narrower cropped spectrum.

        Parameters
        ----------
        wl : np.ndarray()
            The stellar model wavelength(s) in nm.

        fx : np.ndarray()
            The stellar model flux.

        wlmin: float
            The minimum cropping wavelength.

        wlmax: float
            The maximum cropping wavelength.

        pad: float, m/s, optional
            A velocity corresponding to a shift around which the cropped array will
            be padded.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum.
            """
    import lib.test as test
    import numpy as np
    #First run standard tests on input
    test.typetest(wl,np.ndarray,varname='wl in clip_spectrum')
    test.typetest(fx,np.ndarray,varname='fx in clip_spectrum')
    test.typetest(wlmin,[int,float],varname='wlmin in clip_spectrum')
    test.typetest(wlmax,[int,float],varname='wlmax in clip_spectrum')
    test.typetest(pad,[int,float],varname='pad in clip_spectrum')
    test.notnegativetest(pad,varname='pad in clip_spectrum')
    test.nantest(wlmin,varname='wlmin in clip_spectrum')
    test.nantest(wlmax,varname='wlmax in clip_spectrum')
    test.nantest(pad,varname='pad in clip_spectrum')
    if wlmin >= wlmax:
        raise ValueError('wlmin in clip_spectrum should be smaller than wlmax.')
    wlc = wl[(wl >= wlmin) & (wl <= wlmax)]#This is the wavelength grid onto which we will interpolate the final result.
    fxc = fx[(wl >= wlmin) & (wl <= wlmax)]

    if pad > 0:
        wlmin_wide = wlmin/doppler(pad)
        wlmax_wide = wlmax*doppler(pad)
        wlc_wide = wl[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
        fxc_wide = fx[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
        return(wlc,fxc,wlc_wide,fxc_wide)
    else:
        return(wlc,fxc)

def crop_spectrum(wl,fx,pad):
    """This crops a spectrum wl,fx to wlmin and wlmax.
    If pad !=0, it should be set to a positive velocity
    by which the output spectrum will be padded within wlmin,wlmax, to allow
    for velocity shifts. The difference with clip_spectrum above is that this
    routine pads towards the inside (only returning the narrow spectrum), while
    clip_spectrum pads towards the outside, returning both the cropped and the
    padded cropped spectra.

        Parameters
        ----------
        wl : np.ndarray()
            The stellar model wavelength(s) in nm.

        fx : np.ndarray()
            The stellar model flux.

        pad: float, m/s
            A velocity corresponding to a shift around which the cropped array will
            be padded.

        Returns
        -------
        wl,fx: np.array(), np.array()
            The wavelength and flux of the integrated spectrum.
            """
    import lib.test as test
    import numpy as np
    wlmin = min(wl)
    wlmax = max(wl)
    #First run standard tests on input
    test.typetest(wl,np.ndarray,varname=       'wl in crop_spectrum')
    test.typetest(fx,np.ndarray,varname=       'fx in crop_spectrum')
    test.typetest(pad,[int,float],varname=    'pad in crop_spectrum')
    test.notnegativetest(pad,varname=         'pad in crop_spectrum')
    test.nantest(pad,varname=                 'pad in crop_spectrum')
    if wlmin >= wlmax:
        raise ValueError('wlmin in crop_spectrum should be smaller than wlmax.')
    wlmin_wide = wlmin*doppler(pad)#Crop towards the inside
    wlmax_wide = wlmax/doppler(pad)#Crop towards the inside
    wlc_narrow = wl[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
    fxc_narrow = fx[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
    return(wlc_narrow,fxc_narrow)





def vactoair(wlnm):
    """Convert vaccuum to air wavelengths.

        Parameters
        ----------
        wlnm : float, np.array()
            The wavelength(s) in nm.

        Returns
        -------
        wl: float, np.array()
            The wavelengths.
    """
    import lib.test as test
    import numpy
    #First run standard tests on the input
    test.typetest(wlnm,[int,float,numpy.ndarray],varname='wlnm in vactoair')
    test.notnegativetest(wlnm,varname='wlnm in vactoair')
    test.nantest(wlnm,varname='wlnm in vactoair')
    wlA = wlnm*10.0
    s = 1e4/wlA
    f = 1.0 + 5.792105e-2/(238.0185e0 - s**2) + 1.67917e-3/( 57.362e0 - s**2)
    return(wlA/f/10.0)







import jax
from jax import jit
import jax.numpy as jnp
from functools import partial
import numpy as np
from jax import lax



@jit
def limb_darkening(r,a1,a2):
    """Evaluate quadratic limb darkening, taken straight from wikipedia, eq 1.
    (https://en.wikipedia.org/wiki/Limb_darkening). Provide the z-coordinate of grid cell i,j
    and multiply the resulting weight against spectrum i,j.

    WARNING: WHAT FORMALISM DOES THIS USE? WIKIPEDIA IS NOT AUTHORITATIVE. WE NEED TO REFERENCE COMMONLY USED
    LIMB DARKENING LAWS.

    z = sin(psi), with psi the angle of incidence.
    psi = arcsin(z)

    I(z) = I(0) * (a0 + a1*cos(psi) + a2*cos(psi)**2 + ....)
    with a0 = 1 - (a1+a2+...).

        Parameters
        ----------
        r : float, np.array()
            The projected distance from the stellar center in units of stellar radii,
            bounded between 0 and 1.
        u1,u2: float,
            Linear and quadratic limb-darkening coefficients.
        Returns
        -------
        w: float
            The weight, <= 1.
    """

    # Note that this is equal as the route with the psi's:
    # psi = jnp.arcsin(r)
    # a0 = 1-a1-a2
    # I = a0+a1*np.cos(psi)+a2*np.cos(psi)**2
    # Because cos(arcsin(r)) = sqrt(1-r**2)

    # I suspect that this is trivally integrable.
    z = 1-r**2
    
    a0 = 1-a1-a2
    I = a0+a1*jnp.sqrt(z)+a2*z #The last z is the sqrt squared.
    return(I)


#Agruably this should be in integrate.py.
@jit
def vert_int_q_ld(x,a1,a2):
    """Evaluate the vertical integral of the quadratic limb darkening law, integrating in the +y
    direction from the x-axis. Evaluated at some value of x. This is analytical and obtaining
    this is not very pretty.

        Parameters
        ----------
        x : float, array-like
            The projected distance from the stellar center in units of stellar radii,
            bounded between 0 and 1, along the x-axis.
        a1,a2: float,
            Linear and quadratic limb-darkening coefficients.
        Returns
        -------
        II : float, array-like
            The integral at location x, of the same dimension as x.
    """ 
    # I have established that this returns the same answer as numerically
    # integrating the limb darkened flux disk in the y direction.
    U = 1-x**2
    a0 = 1-a1-a2
    return( jnp.sqrt(U) * (a0+a2-a2*x**2 - a2*U/3) + jnp.pi*a1/4 * U)


#Agruably this should be in integrate.py.
@jit
def vert_int_q_ld_bounded(x,ymin,ymax,a1,a2):
    """Evaluate the vertical integral of the quadratic limb darkening law, integrating in the +y
    direction from the x-axis between ymin and ymax. Evaluated at some value of x. This is analytical 
    and obtaining this is not very pretty. Equivalent to vert_int_q_ld but with boundaries specified.

        Parameters
        ----------
        x : float, array-like
            The projected distance from the stellar center in units of stellar radii,
            bounded between 0 and 1, along the x-axis.
        a1,a2: float,
            Linear and quadratic limb-darkening coefficients.
        Returns
        -------
        II : float, array-like
            The integral at location x, of the same dimension as x.
    """ 


    def indef_int(x,y_l,a1,a2):
        a0 = 1-a1-a2
        U1 = jnp.sqrt(1-x**2-y_l**2)
        U2 = 1-x**2
        I = y_l * (a0 + a2 - a2*x**2 -a2/3*y_l**2)   +   a1/2*(y_l*U1 + U2*jnp.atan(y_l / U1))
        return(I)
    

    return(indef_int(x,ymax,a1,a2)-indef_int(x,ymin,a1,a2))





#Agruably this should be in integrate.py.
@jit
def circ_int_q_ld(a1,a2):
    """This provides the disk-integrated flux of a limb-darkened disk. This acts as a 
    normalization constant. Deriving it is not pretty but the result is.

        Parameters
        ----------
        a1,a2: float,
            Linear and quadratic limb-darkening coefficients.

        Returns
        -------
        II : float
            The disk-integral.
    """ 
    a0 = 1-a1-a2
    return(2*jnp.pi * (a0/2 + a1/3 + a2/4))



def constant_velocity_grid(lam_min, lam_max, N):
    """
    Generate a log-uniform wavelength grid.

    Parameters
    ----------
    lam_min : float
        Starting wavelength (must be > 0).
    lam_max : float
        Ending wavelength (must be > 0 and > lam_min).
    N : int
        Number of wavelength points to generate.

    Returns
    -------
    lam : (N,) jnp.ndarray
        Wavelength array uniformly spaced in log(lambda).
    dloglam : float
        Constant spacing in log(lambda).
    """
    import numpy as np
    lam_min = float(lam_min)
    lam_max = float(lam_max)

    if lam_min <= 0 or lam_max <= 0:
        raise ValueError("Wavelengths must be positive for log spacing.")
    if lam_max <= lam_min:
        raise ValueError("lam_max must be greater than lam_min.")
    if N < 2:
        raise ValueError("N must be at least 2.")

    # Uniform grid in log lambda
    log_min = np.log(lam_min)
    log_max = np.log(lam_max)

    dloglam = (log_max - log_min) / (N - 1)
    loglam = log_min + dloglam * jnp.arange(N)

    lam = jnp.exp(loglam)
    return lam, dloglam




def prepare_gaussian_convolver(wl, R, nsigma=4,percentile=1):
    """
    Precompute kernel window half-width required by the jitted convolution function below.

    wl : 1D numpy array (not jax array!) of wavelength grid
    R  : resolving power
    nsigma : minimal truncation radius of Gaussian
    """

    wl = np.asarray(wl)


    if min(np.diff(wl)) <=0:
        raise Exception("Wavelength axis should be strictly ascending.")
    # Estimate an outlier-robust estimate of the smallest wavelength spacing in the grid
    # taking the 1% percentile by default (to avoid rare cases with insanely small wavelength shifts)
    # likely due to accidental repetition of wavelength values or sub-pixel shifted-overlaps.
    # This is quite pathological but at least now we control for it.
    dlam_min = np.percentile(np.diff(wl),percentile)

    # Largest wavelength -> largest sigma
    sigma_max = wl.max() / (R * 2.0 * jnp.sqrt(2.0 * jnp.log(2.0)))

    # Convert sigma from wavelength units to pixels
    radius = int(np.ceil(nsigma * sigma_max / dlam_min))

    return int(radius)



@partial(jax.jit,static_argnames=['radius'])
def convolve_gaussian_explicit(wl, flux, R, radius):
    """
    wl   : (nλ,)
    flux : (nλ,) or (nt, nλ)
    R    : resolving power
    radius: The half-width of the window that defines the mask onto which the gaussian kernel is defined.
        This needs to be large enough to cover even the widest FWHM (at the longest wl). Use the
        function prepare_gaussian_convolver to get a suitable estimate.
    """

    # Always operate on 2D arrays
    flux = jnp.atleast_2d(flux)     # (nt, nλ)

    nlam = wl.size
    nt   = flux.shape[0]
    window_size = 2*radius+1

    def convolve_one_pixel(i):
        lam_i = wl[i]
        # Gaussian sigma at this wavelength
        sigma_i = lam_i / (R * 2.0 * jnp.sqrt(2.0 * jnp.log(2.0)))
        # ---- fixed-size slice around pixel i ----
        start = jnp.clip(i - radius, 0, nlam - window_size)
        wl_win   = lax.dynamic_slice(wl,   (start,), (window_size,))
        flux_win = lax.dynamic_slice(flux, (0, start), (nt, window_size))
        # ---- Gaussian weights ----
        dlam = wl_win - lam_i
        weights = jnp.exp(-0.5 * (dlam / sigma_i) ** 2)
        # Normalize weights
        weights = weights / jnp.sum(weights)
        # Weighted sum over wavelength window
        return jnp.sum(flux_win * weights[None, :], axis=1)
    # Vectorize over wavelength pixels
    result = jax.vmap(convolve_one_pixel)(jnp.arange(nlam))
    # result shape: (nλ, nt) → transpose
    return result.T




@partial(jax.jit, static_argnames=("nsig"))
def convolve_gaussian_constant_dlogl(dlogl, flux, R, nsig=4):
    """
    Gaussian spectrograph convolution on constant log-lambda grid.
    Full-length kernel with Gaussian tails clipped beyond nsig sigma.

    Parameters
    ----------
    wl : (Nλ,)
        Log-spaced wavelength grid
    flux : (Nλ,) or (Nt, Nλ)
        Spectrum or time-series spectra
    R : float
        Resolving power
    nsig : float (static)
        Gaussian truncation radius in sigma

    Returns
    -------
    convolved : same shape as flux
    """

    # ---- ensure flux is 2D (Nt, Nλ) ----
    flux = jnp.atleast_2d(flux)
    Nt, Nlam = flux.shape

    # ---- Gaussian sigma in log λ ----
    c = 2.0 * jnp.sqrt(2.0 * jnp.log(2.0))   # FWHM→σ conversion
    sigma_loglam = 1.0 / (R * c)
    sigma_pix = sigma_loglam / dlogl

    # -------------------------------------------------
    # Build FULL-LENGTH circular Gaussian kernel
    # -------------------------------------------------

    # pixel indices
    i = jnp.arange(Nlam)

    # circular distance from pixel 0
    dist = jnp.minimum(i, Nlam - i)

    # Gaussian kernel
    kernel = jnp.exp(-0.5 * (dist / sigma_pix) ** 2)

    # ---- clip far Gaussian tails ----
    kernel = jnp.where(dist <= nsig * sigma_pix, kernel, 0.0)

    # normalize
    kernel = kernel / jnp.sum(kernel)

    # -------------------------------------------------
    # FFT convolution
    # -------------------------------------------------
    kernel_fft = jnp.fft.rfft(kernel)
    flux_fft   = jnp.fft.rfft(flux, axis=1)

    convolved_fft = flux_fft * kernel_fft[None, :]
    convolved = jnp.fft.irfft(convolved_fft, n=Nlam, axis=1)

    # return original dimensionality
    return convolved













def blur_spec(wl,spec,dv,truncsize = 20.0):
    """This function takes a spectrum, and blurs it using either a
    Gaussian kernel or a box kernel, which have a FWHM width of dv km/s everywhere.
    Meaning that the width changes dynamically on a constant d-lambda grid.
    Because the kernel needs to be recomputed on each element of the wavelength axis
    individually, this operation is much slower than convolution with
    a constant kernel, in which a simple shifting of the array, rather than a recomputation
    of the kernel is sufficient."""

    # print('I MAY NOT WANT TO USE BLUR-SPEC BECAUSE IT IS SLOW, AT LEAST IN BOX MODE.')
    # print('AND I HAVE NOT THOROUGHLY BENCHMARKED IT.')
    import numpy as np
    from matplotlib import pyplot as plt
    import astropy.constants as const
    import lib.test as test
    # Do not perform tests because this thing is in a double forloop.
    # test.typetest(dv,float,varname='dv in doppler(dv)')
    c = const.c.value#Comes out in m/s.
    test.typetest(dv,float,varname='dv in blur_spec')
    test.typetest(wl,np.ndarray,varname='w in blur_specl')
    test.typetest(spec,np.ndarray,varname='spec in blur_spec')
    test.typetest(truncsize,float,varname='truncsize in blur_spec')

    #truncsize=8.0#The gaussian is truncated at 8 sigma.

    sig_dv = dv / 2*np.sqrt(2.0*np.log(2)) #Transform FWHM to Gaussian sigma. In km/s.

    d_kernel=np.array([-1,0,1])/2.0
    deriv = convolve(wl,d_kernel)
    #l*dv/c=dl
    dwl=wl*dv/const.c.value*1000.0
    sig_wl=wl*sig_dv/const.c.value*1000.0#in nm
    sig_px=sig_wl/deriv
    trunc_dist=np.round(sig_px*truncsize).astype(int)

    spec_blurred=spec*0.0
    summm = []
    # pdb.set_trace()
    for i in range(0,len(wl)):
        #Im going to select wl in a bin so that I dont need to evaluate a gaussian over millions of points that are all zero
        binstart=max([0,i-trunc_dist[i]])
        binend=i+trunc_dist[i]
        k = gaussian(wl[binstart:binend],1.0,wl[i],sig_wl[i])
        k_n=k/np.sum(k)
        summm.append(np.sum(k))
        # try:
        # print(len(k_n),len(spec[binstart:binend]),len(wl[binstart:binend]))
        try:
            spec_blurred[i]=np.sum(k_n*spec[binstart:binend])
        except:
            print(len(wl),len(spec))
            pdb.set_trace()
        # except:
            # pdb.set_trace()
        #To speed up, need to select wl and then append with zeroes. <= what does that mean? Jens 03 mar 18
    # plt.plot(wl,summm)
    # plt.show()
    return(spec_blurred)


def smooth(fx,w,mode='gaussian',edge_degree=1):
    """This function takes a spectrum, and blurs it using either a
    Gaussian kernel or a box kernel, which have a FWHM width of w px everywhere.
    Meaning that the width changes dynamically on a constant d-lambda grid.
    """
    import numpy as np
    import lib.test as test

    test.typetest(w,float,'w')
    test.typetest(fx,np.ndarray,'fx')
    test.typetest(mode,str,'mode')
    test.typetest(edge_degree,int,'edge_degree')

    truncsize=8.0#The gaussian is truncated at 8 sigma.
    shape=np.shape(fx)

    sig_w = w / 2*np.sqrt(2.0*np.log(2)) #Transform FWHM to Gaussian sigma. In km/s.
    trunc_dist=np.round(sig_w*truncsize).astype(int)

    #First define the kernel.
    kw=int(np.round(truncsize*sig_w*2.0))
    if kw % 2.0 != 1.0:#This is to make sure that the kernel has an odd number of
    #elements, and that it is symmetric around zero.
        kw+=1

    kx=findgen(kw)
    kx-=np.mean(kx)#This must be centered around zero. Doing a hardcoded check:
    if (-1.0)*kx[-1] != kx[0]:
        print(kx)
        raise Exception("ERROR in box_smooth: Kernel could not be made symmetric somehow. Attempted kernel grid is printed above. Kernel width is %s pixels." % kw)

    if mode == 'gaussian':
        k=gaussian(kx,1.0,0.0,sig_w)

    k/=np.sum(k)
    return(convolve(fx,k,edge_degree))

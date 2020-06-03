#This file contains a collection of functions that operate on spectra.


def doppler(dv):
    """This computes the relativistic doppler parameter.
        Parameters
        ----------
        dv : float
            The radial velocity in m/s

        Returns
        -------
        x:  float
            Returns the doppler parameter.
    """
    import astropy.constants as const
    import numpy as np
    import lib.test as test
    # Do not perform tests because this thing is in a double forloop.
    # test.typetest(dv,float,varname='dv in doppler(dv)')
    c = const.c.value#Comes out in m/s.
    beta = dv/c
    return(np.sqrt((1+beta)/(1-beta)))#Relativistic Doppler effect.


def shift(wl,wl_wide,fx,dv):
    """Doppler-shift a spectrum.

        Parameters
        ----------
        wl : np.array()
            The wavelength axis of the spectrum.
        fx : np.array()
            The corresponding flux axis.
        dv : int, float
            The radial velocity (positive = redshift) in m/s.


        Returns
        -------
        fx: np.array()
            The shifted spectrum evaluated on the original wavelength axis.
            The missing edge is set to NaN. The function calling this should better
            be able to handle NaNs.
    """
    from scipy.interpolate import interp1d
    import astropy.constants as const
    import numpy as np
    import pdb

    #Do not run tests on the input here because this thing is in a double forloop somewhere.
    wl_shifted =  wl_wide*doppler(dv)#Relativistic Doppler effect.
    fx_i=interp1d(wl_shifted,fx,bounds_error=False)#Put NaNs at the edges.
    fx_shifted=fx_i(wl)#Interpolated onto the narrower wl.
    return(fx_shifted)


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



def limb_darkening_old(mu,u1,u2):
    """Evaluate quadratic limb darkening, taken from de Mooij 2017,
    (https://arxiv.org/pdf/1709.00680.pdf). Provide the mu-angle of grid cell i,j
    and multiply the resulting weight against spectrum i,j.

    This formula has been updated for a
    potential typo from (1-mu)**2 to (1-mu**2)

        Parameters
        ----------
        mu : float, np.array()
            the cosine of the angle between the line of sight and the stellar
            surface vertical.
        u1,u2: float,
            Linear and quadratic limb-darkening coefficients.
        Returns
        -------
        w: float
            The weight, <= 1.
    """
    w=(1-u1*(1-mu)-u2*(1-mu**2))

    return(w)

def limb_darkening(z,a1,a2):
    """Evaluate quadratic limb darkening, taken straight from wikipedia, eq 1.
    (https://en.wikipedia.org/wiki/Limb_darkening). Provide the z-coordinate of grid cell i,j
    and multiply the resulting weight against spectrum i,j.

    z = sin(psi), with psi the angle of incidence.
    psi = arcsin(z)

    I(z) = I(0) * (a0 + a1*cos(psi) + a2*cos(psi)**2 + ....)
    with a0 = 1 - (a1+a2+...).

        Parameters
        ----------
        z : float, np.array()
            The projected distance from the stellar center in units of stellar radii,
            bounded between 0 and 1.
        u1,u2: float,
            Linear and quadratic limb-darkening coefficients.
        Returns
        -------
        w: float
            The weight, <= 1.
    """
    import lib.test as test
    #First run standard tests on the input
    test.typetest(z,[int,float],varname='z in limb_darkening')
    test.typetest(a1,[int,float],varname='a1 in limb_darkening')
    test.typetest(a2,[int,float],varname='a2 in limb_darkening')
    test.nantest(z,varname='z in limb_darkening')
    test.nantest(a1,varname='a1 in limb_darkening')
    test.nantest(a2,varname='a2 in limb_darkening')

    import numpy as np
    if z > 1 or z < 0:
        raise ValueError("z coordinate should be in [0,1].")


    psi = np.arcsin(z)
    a0 = 1-a1-a2
    I = a0+a1*np.cos(psi)+a2*np.cos(psi)**2
    return(I)

def findgen(n,int=False):
    """This is basically IDL's findgen function.
    a = findgen(5) will return an array with 5 elements from 0 to 4:
    [0,1,2,3,4]
    """
    import numpy as np
    if int:
        return np.linspace(0,n-1,n).astype(int)
    else:
        return np.linspace(0,n-1,n)



def convolve(array,kernel,edge_degree=1):
    """It's unbelievable, but I could not find the python equivalent of IDL's
    /edge_truncate keyword, which truncates the kernel at the edge of the convolution.
    Therefore, unfortunately, I need to code a convolution operation myself.
    Stand by to be slowed down by an order of magnitude #thankspython.

    Nope! Because I can just use np.convolve for most of the array, just not the edge...

    So the strategy is to extrapolate the edge of the array using a polynomial fit
    to the edge elements. I fit over a range that is twice the length of the kernel.

    Example: y_blurred = convolve(x,y,edge_degree = 2)
    For a convolution where the edge is extrapolated with a second degree polynomial.
    """
    import numpy as np
    import lib.test as test
    test.typetest(edge_degree,int,'edge_degree in convolve')
    test.typetest(array,np.ndarray,'array in convolve')
    test.typetest(kernel,np.ndarray,'kernel in convolve')

    if len(kernel) >= len(array)/2:
        raise Exception("Error in convolution: Kernel length is larger than half of the array. Can't extrapolate over that length. And you probably don't want to be doing a convolution like that, anyway.")

    if len(kernel) % 2 != 1:
        raise Exception('Error in convolution: Kernel needs to have an odd number of elements.')

    #Perform polynomial fits at the edges.
    x=findgen(len(array))
    fit_left=np.polyfit(x[0:len(kernel)*2],array[0:len(kernel)*2],edge_degree)
    fit_right=np.polyfit(x[-2*len(kernel)-1:-1],array[-2*len(kernel)-1:-1],edge_degree)

    #Pad both the x-grid (onto which the polynomial is defined)
    #and the data array.
    pad=findgen( (len(kernel)-1)/2)
    left_pad=pad-(len(kernel)-1)/2
    right_pad=np.max(x)+pad+1
    left_array_pad=np.polyval(fit_left,left_pad)
    right_array_pad=np.polyval(fit_right,right_pad)

    #Perform the padding.
    x_padded = np.append(left_pad , x)
    x_padded = np.append(x_padded , right_pad) #Pad the array with the missing elements of the kernel at the edge.
    array_padded = np.append(left_array_pad,array)
    array_padded = np.append(array_padded,right_array_pad)

    #Reverse the kernel because np.convol does that automatically and I don't want that.
    #(Imagine doing a derivative with a kernel [-1,0,1] and it gets reversed...)
    kr = kernel[::-1]
    #The valid keyword effectively undoes the padding, leaving only those values for which the kernel was entirely in the padded array.
    #This thus again has length equal to len(array).
    return np.convolve(array_padded,kr,'valid')


def gaussian(x,A,mu,sig,cont=0.0):
    import numpy as np
    """This produces a gaussian function on the grid x with amplitude A, mean mu
    and standard deviation sig. Will need to expand it with a version that has
    a polynomial continuum in the same way that IDL does it."""
    return A * np.exp(-0.5*(x - mu)/sig*(x - mu)/sig)+cont

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
    import pdb
    from matplotlib import pyplot as plt
    import time
    import astropy.constants as const
    import numpy as np
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
    plt.plot(wl,summm)
    plt.show()
    return(spec_blurred)


def smooth(fx,w,mode='gaussian',edge_degree=1):
    """This function takes a spectrum, and blurs it using either a
    Gaussian kernel or a box kernel, which have a FWHM width of w px everywhere.
    Meaning that the width changes dynamically on a constant d-lambda grid.
    """

    import numpy as np
    import lib.test as test
    import pdb

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

#This file contains a collection of functions that operate on spectra.


def doppler(dv):
    """This computes the relativistic doppler parameter."""
    import astropy.constants as const
    import numpy as np
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

    wl_shifted =  wl_wide*doppler(dv)#Relativistic Doppler effect.
    fx_i=interp1d(wl_shifted,fx,bounds_error=False)#Put NaNs at the edges.
    fx_shifted=fx_i(wl)#Interpolated onto the narrower wl.
    # pdb.set_trace()
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
        wl : np.array()
            The stellar model wavelength(s) in nm.

        fx : np.array()
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



    wlc = wl[(wl >= wlmin) & (wl <= wlmax)]#This is the wavelength grid onto which we will interpolate the final result.
    fxc = fx[(wl >= wlmin) & (wl <= wlmax)]

    if pad > 0:
        wlmin_wide = wlmin*doppler(pad)#x2 to have some margin.
        wlmax_wide = wlmax*doppler(pad)
        wlc_wide = wl[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
        fxc_wide = fx[(wl >= wlmin_wide) & (wl <= wlmax_wide)]
        return(wlc,fxc,wlc_wide,fxc_wide)
    else:
        return(wlc,fxc)


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
    wlA = wlnm*10.0
    s = 1e4/wlA
    f = 1.0 + 5.792105e-2/(238.0185e0 - s**2) + 1.67917e-3/( 57.362e0 - s**2)
    return(wlA/f/10.0)



def limb_darkening(mu,u1,u2):
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

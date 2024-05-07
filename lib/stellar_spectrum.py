def get_spectrum(T,logg,Z,a):
    """Querying a PHOENIX photosphere model, either from disk or from the
        PHOENIX website if the spectrum hasn't been downloaded before, and
        returning the names of the files in the data subfolder. These may then
        be read by the wrapper function read_spectrum() below.

        Parameters
        ----------
        T : int, float
            The model photosphere temperature. Acceptable values are:
            2300 - 7000 in steps of 100, and 7200 - 12000 in steps of 200.
        logg : int, float
            The model log(g) value. Acceptable values are:
            0.0 - 6.0 in steps of 0.5.
        Z : int, float
            The model metallicity [Fe/H] value. Acceptable values are:
            -4.0, -3.0, -2.0 and -1.5 to +1.0 in steps of 0.5.
        a : int, float
            The model alpha element enhancement [alpha/M]. Acceptable values are:
            -0.2 to 1.2 in steps of 0.2, but only for Fe/H of -3.0 to 0.0.

        Returns
        -------
        wlname, specname: str, str
            The names of the files containing the wavelength and flux axes of
            the requested spectrum.
        """
    import requests
    import shutil
    import urllib.request as request
    from contextlib import closing
    import sys
    import os.path
    import lib.test as test

    #First run standard tests on the input
    test.typetest(T,[int,float],varname='T in get_spectrum')
    test.typetest(logg,[int,float],varname='logg in get_spectrum')
    test.typetest(Z,[int,float],varname='metallicity in get_spectrum')
    test.typetest(a,[int,float],varname='alpha in get_spectrum')
    test.nantest(T,varname='T in get_spectrum')
    test.nantest(logg,varname='logg in get_spectrum')
    test.nantest(Z,varname='Z in get_spectrum')
    test.nantest(a,varname='a in get_spectrum')
    test.notnegativetest(T,varname='T in get_spectrum')


    #This is where phoenix spectra are located.
    root = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/'

    #We assemble a combination of strings to parse the user input into the URL,
    z_string = '{:.1f}'.format(float(Z))
    if Z > 0:
        z_string = '+'+z_string
    else:
        z_string = '-'+z_string
    a_string=''
    if a > 0:
        a_string ='.Alpha=+'+'{:.2f}'.format(float(a))
    if a < 0:
        a_string ='.Alpha='+'{:.2f}'.format(float(a))
    t_string = str(int(T))
    if T < 10000:
        t_string = '0'+t_string
    g_string = '-'+'{:.2f}'.format(float(logg))


    #These are URLS for the input files.
    waveurl = root+'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
    specurl = root+'PHOENIX-ACES-AGSS-COND-2011/Z'+z_string+a_string+'/lte'+t_string+g_string+z_string+a_string+'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

    #These are the output filenames, they will also be returned so that the wrapper
    #of this function can take them in.
    wavename = 'data/PHOENIX/WAVE.fits'

    if os.path.isdir('data/PHOENIX/') == False:
        os.makedirs('data/PHOENIX/')


    specname = 'data/PHOENIX/lte'+t_string+g_string+z_string+a_string+'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

    #If it doesn't already exists, we download them, otherwise, we just pass them on:
    if os.path.exists(wavename) == False:
        with closing(request.urlopen(waveurl)) as r:
            with open(wavename, 'wb') as f:
                shutil.copyfileobj(r, f)
    if os.path.exists(specname) == False:
        print(specurl)
        with closing(request.urlopen(specurl)) as r:
            with open(specname, 'wb') as f:
                shutil.copyfileobj(r, f)
    return(wavename,specname)


def read_spectrum(T,logg,metallicity=0.0,alpha=0.0):
    """Wrapper for the function get_spectrum() above, that checks that the input
        T, log(g), metallicity and alpha are in accordance with what is provided by
        PHOENIX (as of November 1st, 2019), and passes them on to get_spectrum().


        Parameters
        ----------
        T : int,float
            The model photosphere temperature. Acceptable values are:
            2300 - 7000 in steps of 100, and 7200 - 12000 in steps of 200.
        logg : int,float
            The model log(g) value. Acceptable values are:
            0.0 - 6.0 in steps of 0.5.
        metallicity : int,float (optional, default = 0)
            The model metallicity [Fe/H] value. Acceptable values are:
            -4.0, -3.0, -2.0 and -1.5 to +1.0 in steps of 0.5.
            If no location is given, the ``location`` attribute of the Time
            object is used.
        alpha : int,float (optional, default = 0)
            The model alpha element enhancement [alpha/M]. Acceptable values are:
            -0.2 to 1.2 in steps of 0.2, but only for Fe/H of -3.0 to 0.0.

        Returns
        -------
        wl,f : np.array(),np.array()
            The wavelength (nm) and flux (erg/s/cm^2/cm) axes of the requested
            spectrum.
        """
    import numpy as np
    import sys
    import astropy.io.fits as fits
    import lib.test as test
    phoenix_factor = np.pi#This is a factor to correct the PHOENIX spectra to produce the correctly calibrated output flux.
    #If you compare PHOENIX to other models, a factor of 3 seems to be missing. Brett Morris encountered this problem when trying
    #to use PHOENIX to predict the received flux of real sources, to construct an ETC. He also needed a factor of 3, which he
    #interpreted as a missing factor of pi. So pi it is, until something else is deemed better (or until PHOENIX updates
    #their grid).


    #First run standard tests on the input:
    test.typetest(T,[int,float],varname='T in read_spectrum')
    test.typetest(logg,[int,float],varname='logg in read_spectrum')
    test.typetest(metallicity,[int,float],varname='metallicity in read_spectrum')
    test.typetest(alpha,[int,float],varname='alpha in read_spectrum')
    test.nantest(T,varname='T in get_spectrum')
    test.nantest(logg,varname='logg in get_spectrum')
    test.nantest(metallicity,varname='Z in get_spectrum')
    test.nantest(alpha,varname='a in get_spectrum')
    test.notnegativetest(T,varname='T in get_spectrum')


    #These contain the acceptable values.
    T_a = np.concatenate((np.arange(2300,7100,100),np.arange(7200,12200,200)))
    logg_a = np.arange(0,6.5,0.5)
    FeH_a = np.concatenate((np.arange(-4,-1,1),np.arange(-1.5,1.5,0.5)))
    alpha_a = np.arange(0,1.6,0.2)-0.2

    #Check that the input is contained in them:
    if T in T_a and metallicity in FeH_a and alpha in alpha_a and logg in logg_a:
        if (metallicity < -3 or metallicity > 0.0) and alpha !=0:
            print('Error: Alpha element enhancement != 0 only available for -3.0<[Fe/H]<0.0')
            sys.exit()
        #If so, retrieve the spectra:
        wavename,specname=get_spectrum(T,logg,metallicity,alpha)
        f = fits.getdata(specname)
        w = fits.getdata(wavename)
        return(w/10.0,f/phoenix_factor)#Implicit unit conversion here.
    else:
        print('Error: Provided combination of T, log(g), Fe/H and a/M is out of bounds.')
        print('T = %s, log(g) = %s, Fe/H = %s, a/M = %s' % (T,logg,metallicity,alpha))
        print('The following values are accepted:')
        print('T:')
        print(T_a)
        print('Log(g):')
        print(logg_a)
        print('Fe/H:')
        print(FeH_a)
        print('a/M:')
        print(alpha_a)
        print('Alpha element enhancement != 0 only available for -3.0<[Fe/H]<0.0')
        sys.exit()


def get_spectrum_pysme(wave_start, wave_end, T, logg, Z, linelist = '', mu=[], abund = {}, grid = ''):
    """Constructing a spectrum from pySME. pySME uses the MARCS 2014 grid by
    default. It is also possible to use the ATLAS12 model by setting grid='ATLAS12'.
    Individual elemental abundances can also be specified and updated. If an array of
    mu angles is specified, the spectrum for each mu angle is calculated and
    stored in a list.

        Parameters
        ----------
        wave_start : float
            Start of modelled wavelength range in nm in air
        wave_end :  float
            Ending wavelength range in nm in air
        T : int, float
            The model effective temperature of the star in Kelvin.
        logg : int, float
            The model log(g) value. The surface gravity of the star in log(cgs).
        Z : int, float
            The model metallicity of the star in log_10 relative to the indiviual
            abundances.
        linelist : str
            Path to VALD linelist used to generate spectrum
        mu : np.array()
            Array of mu angles to calculate the spectrum over. If an array is
            specified, a list of fluxes for each mu angle is returned.
        abund : dict
            Specific elemental abudances in dict form, eg. '{'X': 6.4}'. Default
            abundances are solar.
        grid : str
            Name of model atmosphere to be used. The default used by pySME is MARCS
            2014, which spans a temperature range of 2500 - 8000 K, logg between -0.5
            - 5.5 and metallicity -5 - 1. The other option is the ATLAS12 model, which
            can be used by specifying "atlas12.sav". This spans a temperature of
            3500 - 50 000 K and logg between 0 - 5.

        Returns
        -------
        wl,f : np.array(),np.array()
            The wavelength (nm) and flux (erg/s/cm^2/cm) axes of the requested
            spectrum.
        """
    from pysme.sme import SME_Structure as SME_Struct
    from pysme.abund import Abund
    from pysme.synthesize import synthesize_spectrum
    from pysme.solve import solve
    from pysme.linelist.vald import ValdFile
    import numpy as np
    import ast

    #initialise the pySME Object with values provided by the user
    sme = SME_Struct()
    sme.abund = Abund.solar()
    sme.teff, sme.logg, sme.monh = T, logg, Z
    sme.vsini = 0
    sme.normalize_by_continuum = False

    # Change from default grid
    if grid:
        sme.atmo.method = 'grid'
        sme.atmo.source = grid

    # Convert from nm to Angstrom
    wave_start *= 10
    wave_end *= 10

    # Account for velocity shift in wavelengths for StarRotator
    wave_start -= 5
    wave_end += 5

    if abund:
        for i in range(len(abund)):
            sme.abund.update_pattern(updates=ast.literal_eval(abund[i]))

    sme.wran = [[wave_start, wave_end]]
    vald = ValdFile(linelist)
    sme.linelist = vald

    if len(mu) > 0:
        fx = []
        for m in mu:#Loop over all mu angles
            sme.mu = m
            sme = synthesize_spectrum(sme)

            w, f = *sme.wave/10, *sme.synth
            wl = w
            fx.append(f.copy())
            if len(w) < 5:
                print("get_spectrum_pysme: It seems that your returned wavelength grid is shorter than expected. Please make sure that your line list has spectral lines in your specified wavelength range.")
        return(wl, fx)

    else:
        sme = synthesize_spectrum(sme)
        if len(*sme.wave) < 5:
                print("get_spectrum_pysme: It seems that your returned wavelength grid is shorter than expected. Please make sure that your line list has spectral lines in your specified wavelength range.")
        return(*sme.wave/10, *sme.synth)

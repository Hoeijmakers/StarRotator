def get_spectrum(T,logg,Z,a):
    """Querying a PHOENIX photosphere model, either from disk or from the
        PHOENIX website if the spectrum hasn't been downloaded before, and
        returning the names of the files in the data subfolder. These may then
        be read by the wrapper function read_spectrum() below.
        However, if limbs is set to a positive value, the code is switched to a
        mode in which the CLV effect is taken into account. In this event,
        PHOENIX spectra cannot be used, and the SPECTRUM package will be used to
        calculate spectra. For this functionality to work, the user needs to have
        built the SPECTRUM source code (included in this distribution). See the
        documentation for instructions.

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
    #If you compare PHOENIX to SPECTRUM, a factor of 3 seems to be missing. Brett Morris encountered this problem when trying
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

def download_kurucz(outpath):
    from contextlib import closing
    import shutil
    import urllib.request as request
    """ This is a convenience function that downloads all of the stellar
    atmosphere "supermodel" grids of Robert Kurucz to the data/KURUCZ folder.
    Just the standard ones. For use with SPECTRUM."""
    root = 'http://kurucz.harvard.edu/grids/grid'
    suffixes = ['M01','M02','M03','M05','M10','M15','M20','M25','M30','M35','M40','M45','M50','P00','P01','P02','P03','P05','P10']

    for s in suffixes:
        url = root+s.lower()+'/a'+s.lower()+'k2.dat'
        name = outpath + 'a' +s.lower()+'k2.dat'
        print(name,url)
        with closing(request.urlopen(url)) as r:
            with open(name, 'wb') as f:
                shutil.copyfileobj(r, f)

def available_kurucz_models():
    """These hard-code the available Kurucz models, as present on Nov 22, 2019 on
    'http://kurucz.harvard.edu/grids/grid'"""
    import numpy as np
    T_a = np.concatenate((np.arange(3500,13250,250),np.arange(13000,1000,51000)))
    logg_a = np.arange(0,5.5,0.5)
    Z_a = np.round(np.concatenate((np.arange(-5,0,0.5),np.arange(-0.3,0.4,0.1),np.array([0.5,1.0]))),decimals=1)#I need to do a rounding here because np.concatenate() makes a numerical error on the middle array... Weird!
    return(T_a,logg_a,Z_a)

def construct_kurucz_path(Z,root):
    """This constructs the path to the kurucz models after they are downloaded
    with download_kurucz()"""
    mstring = str(Z).replace('.','')
    if Z < 0:
        s = 'm'
        mstring=mstring.replace('-','')
    else:
        s = 'p'
    mpath = root+'a'+s+mstring+'k2.dat'
    return(mpath)






def compute_spectrum(T,logg,Z,mu,wlmin,wlmax,macroturbulence=0.0,mode='an',loud=False):
    """If the CLV effect needs to be taken into account, PHOENIX spectra cannot
        be used, and the SPECTRUM package will be used to calculate spectra.
        For this functionality to work, the user needs to have built the SPECTRUM
        source code (included in this distribution). See the documentation for
        instructions on how to do this.

        This code then acts as a wrapper for SPECTRUM.
        SPECTRUM needs to be made in the lib/SPECTRUM folder, in the same place where
        it is provided along with this distribution. This distribution also
        provides the necessary linelist files, as well as grids of stellar atmosphere
        files obtained from Kurucz. After spectrum is made, this should work out
        of the box.

        Provide a temperature, metallicity and logg, as well as a mu value or
        an np.array of mu values for the wrapper to loop over. The accepted values
        are different from those when using the disk-integrated PHOENIX models.

    T_a = np.concatenate((np.arange(3500,13250,250),np.arange(14000,1000,51000)))
    logg_a = np.arange(0,5.5,0.5)
    Z_a = np.round(np.concatenate((np.arange(-5,0,0.5),np.arange(-0.3,0.4,0.1),np.array([0.5,1.0]))),decimals=1)#I need to do a rounding here because np.concatenate() makes a numerical error on the middle array... Weird!


        Parameters
        ----------
        T : int, float
            The model photosphere temperature. Acceptable values are:
            3500 - 13,000 in steps of 250, and 14,000 - 50,000 in steps of 1,000.
        logg : int, float
            The model log(g) value. Acceptable values are:
            0.0 - 5.0 in steps of 0.5.
        Z : int, float
            The model metallicity [Fe/H] value. Acceptable values are:
            -5.0 - -0.5 in steps of 0.5; -0.3 - +0.3 in steps of 0.1; 0.5; 1.0.
        wlmin: float
            The minimum wavelength to be considered, in units of wl.
        wlmax: float
            The maximum wavelength to be considered, in units of wl.
        mu: float,np.array
            Limb angle(s) (cos(theta))
        macroturbulence: int,float
            Macroturbulent velocity (which broadens the lines) in km/s
        mode: str
            a: Prompt user for custom statoms.dat
            i: Isotope mode. For use with lukeall2.iso.lst linelist.
            n: Silent running.
            m/M: Specific intensities, normalised (m) or true (M).
            f: Gives output in specific intensities. Implied if m or M is set.
        Returns
        -------
        wl,fx
            The wavelengths (nm) and fluxes of the computed spectrum.
            If mu has more than 1 element, then fx is a list of flux axes, each
            of which matches to wl, and each of which matches to its respective
            mu value.
        """
    import subprocess
    import numpy as np
    import lib.test as test
    import os
    import pdb
    import astropy.io.ascii as ascii
    from lib.integrate import statusbar
    import lib.operations as ops
    from os import path
    import sys
    kurucz_root = 'data/KURUCZ/'#THIS IS THE SECOND TIME THIS IS DEFINED. SHOULD MATCH THE DEFINITION IN test_KURUCZ()!
    spath = './lib/SPECTRUM/selectmod'#This points to the supermod function.
    SPECTRUM = './lib/SPECTRUM/spectrum'#This points to spectrum.
    test.dir_exists('./lib/SPECTRUM/',varname = 'SPECTRUM root folder')
    test.dir_exists(kurucz_root,varname='Kurucz model folder')
    test.file_exists(spath,varname=spath+'in compute_spectrum()')
    test.file_exists(SPECTRUM,varname=spath+'in compute_spectrum()')#These four tests need to be moved into some dedicated test routine for SPECTRUM?


    #Hardcoded not to use isotope mode because itneeds the massive lukeall2.iso.lst that I dont want to upload to github.
    #mode = 'ainMf' #Mode of calling SPECTRUM.

    if 'i' in mode:
        isotope = True
    else:
        isotope = False
    if (('m' in mode) or ('M' in mode)) == False:
        mu = 0.0
    if 'a' in mode == False:
        print("WARNING: Compute_spectrum() mode does not include 'a'.")
        print("This will break SPECTRUM as it expects the stdatom file to be")
        print("in the root folder. 'a' points it to where it is supposed to be.")
        print("Adding 'a' to mode variable.")
        mode+='a'

    # if c_norm == True:
    #     mode = mode.lower()#Lower the case of the M.
    # if isotope == False:
    #     mode = mode.replace('i','')
    # if disk_integrated == True:
    #     mode = mode.replace('m','').replace('M','')
    #     mu = 0.0

    #Standard tests on input.
    test.typetest(T,[int,float],varname='T in compute_spectrum')
    test.typetest(logg,[int,float],varname='log(g) in compute_spectrum')
    test.typetest(Z,[int,float],varname='Metallicity in compute_spectrum')
    test.typetest(mu,[float,np.ndarray],varname='mu angle(s) in compute_spectrum')
    test.typetest(wlmin,[int,float],varname='wlmin in compute_spectrum')
    test.typetest(wlmax,[int,float],varname='wlmax in compute_spectrum')
    test.typetest(macroturbulence,[int,float],varname='macroturbulence in compute_spectrum')
    test.nantest(mu,'mu angle(s) in compute_spectrum')
    if wlmax <= wlmin:
        raise ValueError('Wlmax (%s) in compute_spectrum() should be smaller than wlmin (%s)' % (wlmax,wlmin))
    if wlmin < 90:
        raise ValueError('Wlmin (%s) in compute_spectrum() should be larger than 90 nm.' % wlmin)
    if wlmax > 4000:
        raise ValueError('Wlmax (%s) in compute_spectrum() should be less than than 4,000 nm.' % wlmin)

    #These hard-code the acceptable values.
    T_a,logg_a,Z_a = available_kurucz_models()

    if T in T_a and Z in Z_a and logg in logg_a:
        #This builds the Kurucz model atmosphere file, extracted from a 'supermodel'
        #grid. This follows the manual of SPECTRUM, section 4.10.
        #We first construct the filename of the supermodel, using the provided metallicity.
        mpath = construct_kurucz_path(Z,kurucz_root)
        test.file_exists(mpath,varname='Kurucz supermodel')
        #Check that this supermodel actually exists. It should, if test_KURUCZ() was run (which tests all allowed values of Z).

        opath = kurucz_root+'star.mod'#This will contain the temporary .mod file that is the model extracted from the Kurucz supermodel file.
        rpath = 'data/SPECTRUM.res'#This will contain the temporary response file used to < into the bash call of SPECTRUM.
        if isotope == True:
            lpath = 'lib/SPECTRUM/lukeall2.iso.lst'
        else:
            lpath = 'lib/SPECTRUM/merged_linelist.lst'#I made this linelist by manually mergin the other linelists. Will be part of our SPECTRUM fork.
        test.file_exists(lpath,varname=lpath+' (linelist file) in compute_spectrum()')
        sopath = 'data/SPECTRUM.spc'#Spectrum output path, temporary.

        try:
            res = subprocess.check_output([spath,mpath,opath, str(T),str(logg)])#This executes the command line as >>> selectmod supermodel.dat output.mod teff logg
            #So now we have made an atmosphere model file that SPECTRUM can use.
        except:
            print('Unexpected file input error.')
            print('Apparently the provided combination of T, log(g), Fe/H is')
            print('not included in this (%s) Kurucz model grid.' % mpath)
            print('T = %s, log(g) = %s, Fe/H = %s' % (T,logg,Z))

        if isinstance(mu,float):#Convert into a list so that we can loop.
            mu = [mu]
        fx = []#This will contain the flux axes, even if there is only one.

        bashcmd = SPECTRUM+' '+mode+' < '+rpath
        if loud == False:
            bashcmd+=' > temp'#Write to temp to suppress output.
        i=0#Just a counter for the statusbar.
        print('------ Executing bash cmd:  '+bashcmd)
        for m in mu:#Loop over all mu angles, keeping the atmosphere model constant.
            statusbar(i,mu)
            lines = []
            lines.append(opath)
            lines.append(lpath)
            lines.append('lib/SPECTRUM/stdatom.dat')
            if isotope == True:
                lines.append('lib/SPECTRUM/isotope.iso')
            lines.append(sopath)
            lines.append(str(macroturbulence))
            if ('m' in mode) or ('M' in mode):
                lines.append(str(m))
            lines.append(str(wlmin*10.0)+','+str(wlmax*10.0))#Convert to angstroms
            lines.append('0.01')#Hardcoded to 0.01A, i.e. R = 500,000 at 500nm.
            # line1 = opath#Input .mod file.
            # line3 = 'lib/SPECTRUM/stdatom.dat'
            # line5 = 'data/SPECTRUM.spc'#The output file.
            # line6 = str(macroturbulence)#Macroturbulence, default 2kms.
            # line7 = str(m)
            # line8 = str(wlmin*10.0)+','+str(wlmax*10.0)#Convert to angstroms
            # line9 = '0.01'
            # if isotope == True:#This is hardcoded to be off, at present.
            #     line2 = 'lib/SPECTRUM/lukeall2.iso.lst'#The linelist.
            #     line4 = 'lib/SPECTRUM/isotope.iso'#The isotope file.
            #     lines = [line1,line2,line3,line4,line5,line6,line7,line8,line9]
            # else:
            #     line2 = 'lib/SPECTRUM/merged_linelist.lst'
            #     lines = [line1,line2,line3,line5,line6,line7,line8,line9]#Load non-isotope linelist and ingore isotope file on line 4.

            with open(rpath,'w') as f:
                for line in lines:
                    f.write(line+'\n')#This works because the last line of the response
                    #file must end in a carriage return.
            void=os.system(bashcmd)#I use os.system() against the advice of teh internet because subprocesses.run() hangs forever.
            # print('Finished computing SPECTRUM')
            # with open('temp','w') as f:
            # process=subprocess.Popen(bashcmd.split(),stdout=subprocess.PIPE,shell=False)
            # stdout = process.communicate()
            test.file_exists(sopath,varname='SPECTRUM output ('+sopath+')')
            d=ascii.read(sopath)

            wl = d.columns[0].data#This is overwritten each time, but that doesn't matter because each time its the same.
            fx.append(d.columns[1].data)
            #Now we remove all the scratch files we made in this loop.
            # os.remove('temp')
            os.remove(sopath)
            os.remove(rpath)
            i+=1
        os.remove(opath)
        if len(mu) == 1:
            return(ops.airtovac(wl/10.0),fx[0])#Convert back to nm and do airtovac.
        else:
            return(ops.airtovac(wl/10.0),fx)#Convert back to nm and do airtovac.
    else:
        print('Error: Provided combination of T, log(g), Fe/H is out of bounds.')
        print('T = %s, log(g) = %s, Fe/H = %s' % (T,logg,Z))
        print('The following values are accepted:')
        print('T:')
        print(T_a)
        print('Log(g):')
        print(logg_a)
        print('Fe/H:')
        print(list(Z_a))
        sys.exit()

######################
#authors: Jens Hoeijmakers, Julia Seidel, Madeline Lam, Bibiana Prinoth
#Description: Calculates the stellar spectrum
#
#
#####################

#import statements
import sys
import numpy as np
import lib.test as test
import lib.vgrid as vgrid
import lib.plotting as pl
import lib.operations as ops
import lib.stellar_spectrum as spectrum
import lib.integrate as integrate
import pdb
import time
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Circle
import lib.planet_pos as ppos
import copy
#main body of code


class StarRotator(object):
    def __init__(self,wave_start,wave_end,grid_size,star_path='input/demo_star.txt',
    planet_path='input/demo_planet.txt',obs_path='input/demo_observations.txt', linelist_path='',
    input={}):
        """
            Welcome to StarRotator.
            ***********************
            The StarRotator object contains the main functionality of StarRotator.
            Upon initialization, the exoplanet system input is read from file, and the
            model is computed. The parameters needed to initialise are listed below:

            Parameters
            ----------
            wave_start : float
                Start of modelled wavelength range in nm in vacuum.
            wave_end : float
                Ending Wavelength range in nm in vacuum.
            grid_size: int
                Half-width of number of grid cells. Set to values greater than 200
                to limit numerical errors, or less if you are trying things out and just
                want speed.
            star_path : str
                Path to parameter file defining the star. This file should contain the following
                values on separate lines, and the default values are as follows:
                    50000.0     v_eq
                    90.0        stellar i
                    0.0         Differential rotation parameter (alpha)
                    5000.0      Stellar effective temperature (K)
                    0.0         Fe/H
                    4.5         logg
                    0.93        Limb-darkening coefficient u1
                    -0.23       Limb-darkening coefficient u2
                    0           Number of mu angles to consider. For values higher than zero, StarRotator switches to SPECTRUM rather than PHOENIX.
            planet_path: str
                Path to the parameter file defining the planet and its orbit. This file should contain
                the following values on separate lines, and the default values are as follows:
                    3.153       a/Rs
                    0.0         e
                    0.0         omega
                    86.79       Orbital inclination
                    85.0        Projected obliquity
                    0.08228     Rp/Rs
                    1.4811235   Orbital period
                    57095.68572 Transit center time - 2400000.
                    phases      mode, providing the interpretation of the timestamps of the observations:

            obs_path: str
                Path to the parameter file defining the timestamps of the observations.
                If mode (see previous) is set to 'phases', this file is assumed to contain
                a list of orbital phases. Otherwise it should be set to 'times', in which
                case the file should contain the times of the observations, in
                JD - 2400000. These times are assumed to mean the *start* times of
                each observation in the sequence. In addition, a second column should
                be provided giving the exposure time in seconds. This allows StarRotator
                to shift to the mid-times of the exposures. Note that for longer exposure times,
                the output of StarRotator will be less accurate because the signal of multiple
                exposures are convolved together. This effect is greater for aligned orbits.

                By default, the following phases are provided:
                    -0.06
                    -0.055
                    -0.05
                    -0.045
                    -0.04
                    -0.035
                    -0.03
                    -0.025
                    -0.02
                    -0.015
                    -0.01
                    -0.005
                    0.0
                    0.005
                    0.01
                    0.015
                    0.02
                    0.025
                    0.03
                    0.035
                    0.04
                    0.045
                    0.05
                    0.055
                    0.06
            linelist_path: str
                Path to the VALD linelist used for generating a spectrum using pySME.

            Class methods
            -------------
            After initializing the class like a=StarRotator(588.0,592.0,200.0),
            the following methods are available to manipulate the simulation after
            it has been calculated the first time, or to plot the simulation output.
            star.read_system() Reads in (new) exoplanet system files (see above).
            star.compute_spectrum() recomputes the simulation.
            star.plot_residuals() Produces a 2D plot of the residuals.
            star.animate() Produces an animation of the transiting system. One frame
            corresponding to each phase provided.

            Output attributes
            -----------------
            After initializing the class like star = StarRotator(588.0,592.0,200.0),
            the primary simulation output can be accessed as follows:
            star.wl (wavelength array)
            star.spectra (matrix of spectra, 2D np.array)
            star.lightcurve (list of fluxes)
            star.residual (residual after dividing out of transit spectra, 2D np.array)


        """
        self.wave_start=float(wave_start)
        self.wave_end=float(wave_end)
        self.grid_size=int(grid_size)
        self.read_system(star_path=star_path,planet_path=planet_path,obs_path=obs_path,input=input)
        self.linelist_path = linelist_path
        self.compute_spectrum()

        # return(self.wlF,F_out)#This is mad return statement. This whole function should be a class instead.

    def read_system(self,star_path='demo_star.txt',planet_path='demo_planet.txt',
    obs_path='demo_observations.txt',input={}):
        """Reads in the stellar, planet and observation parameters from file; performing
        tests on the input and lifting the read variables to the class object.

        Parameters
        ----------
            star_path : str
                Path to parameter file defining the star.
            planet_path: str
                Path to the parameter file defining the planet and its orbit.
            obs_path: str
                Path to the parameter file defining the timestamps of the observations.
            input: dict
                A dictionary with input parameters, instead of using textfiles as input.
                This allows programmatic control over all StarRotator inputs. The following keys
                need to be defined in any order:
                veq (equatorial velocity, float),
                stelinc (stellar inclination axis, float),
                drr (differential rotation, float),
                T (Teff, float)
                FeH (metallicity, float)
                logg (float)
                u1 (limb darkening parameter 1)
                u2 (limb darkening parameter 2)
                mus (number of mu angles, int)
                R (resolving power, float)
                model (model type, string, either PHOENIX or pySME)
                sma_Rs (a over Rs, float)
                e (eccentricity, float)
                omega (longitude of periastron, float)
                inclination (float)
                obliquity (float)
                RpRs (float)
                P (orbital period, float)
                mode (string, set either to phases or times)
                Tc (transit center time, float, mode set to times)
                phases (numpy array, set to the orbital phase values of the time series)



        """
        if len(input)==0:#If we read input from config files
            planetparams = open(planet_path,'r').read().splitlines()
            starparams = open(star_path,'r').read().splitlines()
            obsparams = open(obs_path,'r').read().splitlines()
            self.velStar = float(starparams[0].split()[0])
            self.stelinc = float(starparams[1].split()[0])
            self.drr = float(starparams[2].split()[0])
            self.T = float(starparams[3].split()[0])
            self.Z = float(starparams[4].split()[0])
            self.logg = float(starparams[5].split()[0])
            self.u1 = float(starparams[6].split()[0])
            self.u2 = float(starparams[7].split()[0])
            self.mus = int(starparams[8].split()[0])
            self.R = float(starparams[9].split()[0])
            self.model = str(starparams[10].split()[0])
            if self.model.lower() in ['pysme','sme']:
                self.grid_model = str(starparams[11].split()[0])
                self.abund = starparams[12:]
            self.sma_Rs = float(planetparams[0].split()[0])
            self.ecc = float(planetparams[1].split()[0])
            self.omega = float(planetparams[2].split()[0])
            self.orbinc = float(planetparams[3].split()[0])
            self.pob = float(planetparams[4].split()[0])#Obliquity.
            self.Rp_Rs = float(planetparams[5].split()[0])
            self.orb_p = float(planetparams[6].split()[0])
            self.transitC = float(planetparams[7].split()[0])#Is this used?
            self.mode = planetparams[8].split()[0]#This is a string.
            times = [] #These are in JD-24000000.0 or in orbital phase.
            for i in obsparams:
                times.append(float(i.split()[0]))
            self.times = np.array(times)
            self.exptimes = []
            if self.mode == 'times':
                for i in obsparams:
                    self.exptimes.append(float(i.split()[1]))
                self.exptimes = np.array(self.exptimes)
        else: #if the input dictionary is set.
            req_keys = ['veq','stelinc','drr','T','FeH','logg','u1','u2','R','mus','model','sma_Rs',
            'e','omega','inclination','obliquity','RpRs','P','Tc','mode']
            for k in req_keys:
                if k not in input:
                    raise Exception(f"Missing key ('{k}') in input dictionary.")
            self.velStar    = float(input['veq'])
            self.stelinc    = float(input['stelinc'])
            self.drr        = float(input['drr'])
            self.T          = float(input['T'])
            self.Z          = float(input['FeH'])
            self.logg       = float(input['logg'])
            self.u1         = float(input['u1'])
            self.u2         = float(input['u2'])
            self.R          = float(input['R'])
            self.mus        = int(input['mus'])
            self.model      = str(input['model'])
            if self.model.lower() in ['pysme','sme']:
                also_req_keys = ['grid_model','abund']
                for k in also_req_keys:
                    if k not in input:
                        raise Exception(f"Missing key ('{k}') in input dictionary.")
                self.grid_model = str(input['grid_model'])
                self.abund      = np.array(input['abund'])
            self.sma_Rs     = float(input['sma_Rs'])
            self.ecc        = float(input['e'])
            self.omega      = float(input['omega'])
            self.orbinc     = float(input['inclination'])
            self.pob        = float(input['obliquity'])#Obliquity.
            self.Rp_Rs      = float(input['RpRs'])
            self.orb_p      = float(input['P'])
            self.transitC   = float(input['Tc'])#Is this used?
            self.mode       = input['mode']#This is a string.
            if self.mode == 'times':
                also_req_keys = ['times','exptimes']
                for k in also_req_keys:
                    if k not in input:
                        raise Exception(f"Missing key ('{k}') in input dictionary.")
                self.times = np.array(input['times'])
                self.exptimes = np.array(input['exptimes'])
                if len(self.times) != len(self.exptimes):
                    raise Exception("Time and exptime arrays should have the same lengths.")
            elif self.mode == 'phases':
                also_req_keys = ['phases']
                for k in also_req_keys:
                    if k not in input:
                        raise Exception(f"Missing key ('{k}') in input dictionary.")
                self.times = np.array(input['phases'])
                self.exptimes = []



        self.Nexp = len(self.times)#Number of exposures.
        self.residual = None
        self.blurred = 0
        try:
            test.typetest(self.wave_start,float,varname='wave_start in input')
            test.nantest(self.wave_start,varname='wave_start in input')
            test.notnegativetest(self.wave_start,varname='wave_start in input')
            test.notnegativetest(self.velStar,varname='velStar in input')
            test.notnegativetest(self.stelinc,varname='stelinc in input')
            test.notnegativetest(self.T,varname='Teff in input')
            test.notnegativetest(self.logg,varname='logg in input')
            test.notnegativetest(self.R,varname='Resolution in input')
            test.notnegativetest(self.orb_p,varname='period in input')
            test.notnegativetest(self.Rp_Rs,varname='RpRs in input')
            test.notnegativetest(self.ecc,varname='e in input')
            test.notnegativetest(self.sma_Rs,varname='e in input')
        #add all the other input parameters
        except ValueError as err:
            print("Parser: ",err.args)

        if self.mus != 0:
            self.mus = np.linspace(0.1,1.0,self.mus)#Uncomment this to run in CLV mode with SPECTRUM.


    def compute_spectrum(self):
        """This is where the main computation takes place. The simulation output
        and other variables are raised to class-wide attributes.

        Parameters
        ----------
            None
        """
        import math
        #Two arrays for the x and y axes
        self.x = np.linspace(-1,1,num=2*self.grid_size) #in units of stellar radius
        self.y = np.linspace(-1,1,num=2*self.grid_size) #in units of stellar radius
        #Calculate the velocity and flux grids
        print('--- Computing velocity/flux grids')
        self.vel_grid = vgrid.calc_vel_stellar(self.x,self.y, self.stelinc, self.velStar, self.drr,  self.pob)
        self.flux_grid = vgrid.calc_flux_stellar(self.x,self.y,self.u1,self.u2)
        if isinstance(self.mus,np.ndarray) != True:#SWITCH BETWEEN PHOENIX (mus=0) AND SPECTRUM
            if self.model.lower() == 'phoenix':
                print('--- Reading spectrum from PHOENIX')
                print('-----T=%sK, log(g)=%s, Z=%s.' % (self.T,self.logg,self.Z))
                wl,fx = spectrum.read_spectrum(self.T,self.logg,metallicity=self.Z)
            elif self.model.lower() in ['pysme','sme']:
                print('--- Creating spectrum from pySME')
                print('-----T=%sK, log(g)=%s, Z=%s.' % (self.T,self.logg,self.Z))
                wl, fx= spectrum.get_spectrum_pysme(self.wave_start, self.wave_end, self.T, self.logg, self.Z, self.linelist_path, grid = self.grid_model)
            else:
                raise Exception('Invalid model spectrum chosen. Input either PHOENIX or pySME in '
                'star.txt')
            self.wl_in = wl*1.0
            self.fx_in = fx*1.0
            print('--- Integrating disk')
            if  self.drr == 0:
                print('------ Fast integration')
                wlF,F = integrate.build_spectrum_fast(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            else:
                print('------ Slow integration')
                wlF,F = integrate.build_spectrum_slow(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
        else:
            if self.model == 'pySME':
                print('--- Computing limb-resolved spectra with pySME')
                print('-----T=%sK, log(g)=%s, Z=%s.'% (self.T,self.logg,self.Z))
                wl, fx_list = spectrum.get_spectrum_pysme(self.wave_start, self.wave_end, self.T, self.logg, self.Z, self.linelist_path, self.mus, self.abund, grid=self.grid_model)
                print('--- Integrating limb-resolved disk')
                wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,self.mus, self.wave_start,self.wave_end,self.x,self.y,self.vel_grid, self.flux_grid)

            else:
                raise Exception('Invalid model spectrum chosen. Input pySME in star.txt')


        self.xp,self.yp,self.zp = ppos.calc_planet_pos(self.sma_Rs, self.ecc, self.omega, self.orbinc, self.pob, self.Rp_Rs, self.orb_p, self.transitC, self.mode, self.times, self.exptimes)


        F_out = np.zeros((self.Nexp,len(F)))
        F_planet = np.zeros((self.Nexp,len(F)))
        flux_out = []
        mask_out = []
        print('--- Building local spectrum')
        for i in range(self.Nexp):
            if isinstance(self.mus,np.ndarray) == True:
                wlp,Fp,flux,mask = integrate.build_local_spectrum_limb_resolved(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx_list,self.mus,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid, self.flux_grid)

            else:
                wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            integrate.statusbar(i,self.Nexp)

            F_out[i,:]=F-Fp
            F_planet[i,:] = Fp
            flux_out.append(flux)
            mask_out.append(mask)
        #This defines the output.
        self.wl = wlF
        self.stellar_spectrum = F
        self.Fp = copy.deepcopy(F_planet)
        self.spectra = copy.deepcopy(F_out)
        self.lightcurve = flux_out
        self.masks = mask_out
        self.residual = self.spectra/self.stellar_spectrum




    def plot_residuals(self):
        """Plot the residuals if available.

        Parameters
        ----------
            None
        Returns
        -------
            None
        """
        import matplotlib.pyplot as plt
        if self.residual is not None:
            plt.pcolormesh(self.wl,self.times,self.residual)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Phase')
            plt.show()


    def convolve_spectral_resolution(self):
        """Convolves the residual to conform to a spectral resolving power.
        An approximation is made is that the convolution of the ratio of the in/out of transit
        spectra is the same as the ratio of the convolutions of the in/out transit spectra. This
        is done because taking the ratio of the convolution leads to numerical errors.

        Parameters
        ----------
            None
        Returns
        -------
            None
        """
        import lib.operations as ops
        import astropy.constants as const
        import copy
        import scipy.ndimage.filters as SNF
        import scipy.interpolate as interp
        from lib.integrate import statusbar as statusbar
        dv = const.c.value / self.R / 1000.0 #in km/s

        if self.blurred != 0:
            raise Exception("Blurring has already been performed. Re-run compute-residuals if you "
            "wish to repeat this operation with a different value of the resolving power. ")

        if self.residual is not None:
            print('--- Blurring')
            for i in range(len(self.spectra)):
                statusbar(i,len(self.spectra))
                res = self.residual[i]
                wl_constant_v,res_constant_v,a = ops.constant_velocity_wl_grid(self.wl,res,3.0)
                res_smooth = ops.smooth(res_constant_v,dv/a)
                self.residual[i] = interp.interp1d(wl_constant_v,res_smooth,
                fill_value='extrapolate',bounds_error=False)(self.wl)
            self.blurred = 1
        else:
            print('--- Skipping blurring because no residuals have been computed.')


    def animate(self):
        """Plots an animation of the transit event, the stellar flux and velocity
        fields, and the resulting transit and line-shape variations. The animation
        is located in the subfolder `anim`. Anything located in this folder prior
        to running this function will be removed. The creation of the animation gif file
        requires that the imagemagick package is installed on your system. If it is not, an error
        will be raised. The PNG files will have been saved in the anim directory nonetheless.

        Parameters
        ----------
            None
        Returns
        -------
            None
        """
        import matplotlib.pyplot as plt
        import lib.integrate as integrate
        import numpy as np
        from matplotlib.patches import Circle
        import shutil
        import os
        import os.path
        if os.path.isdir('anim/') == True:
            shutil.rmtree('anim/')#First delete the contents of the anim folder.
        os.mkdir('anim/')
        minflux = min(self.lightcurve)
        F=self.stellar_spectrum
        for i in range(self.Nexp):
            mask = self.masks[i]
            fig,ax = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
            ax[0][0].pcolormesh(self.x,self.y,self.flux_grid*mask,vmin=0,
            vmax=1.0*np.nanmax(self.flux_grid),cmap='autumn')
            ax[1][0].pcolormesh(self.x,self.y,self.vel_grid*mask,cmap='bwr')
            ax[0][0].axes.set_aspect('equal')
            ax[1][0].axes.set_aspect('equal')
            if self.zp[i] > 0.0:
                planet1 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black',
                edgecolor='black', lw=1)
                planet2 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black',
                edgecolor='black', lw=1)
                ax[0][0].add_patch(planet1)
                ax[1][0].add_patch(planet2)
            ax[0][0].set_ylim((min(self.y),max(self.y)))
            ax[1][0].set_ylim((min(self.y),max(self.y)))
            ax[0][1].plot(self.times[0:i],self.lightcurve[0:i],'.',color='black')
            ax[0][1].set_xlim((min(self.times),max(self.times)))
            ax[0][1].set_ylim((minflux-0.1*self.Rp_Rs**2.0),1.0+0.1*self.Rp_Rs**2)
            ax[1][1].plot(self.wl,F/np.nanmax(F),color='black',alpha = 0.5)
            ymin = np.nanmin(F/np.nanmax(F))
            ymax = np.nanmax(F/np.nanmax(F))
            linedepth = ymax - ymin
            ax[1][1].plot(self.wl,self.spectra[i]/np.nanmax(self.spectra[i]),color='black')
            # ax[1][1].set_xlim((588.5,590.2))
            yl = (ymin-0.1*linedepth,ymax+0.3*linedepth)
            ax[1][1].set_ylim(yl)
            ax2 = ax[1][1].twinx()
            ax2.plot(self.wl,(self.spectra[i])*np.nanmax(F)/F/np.nanmax(self.spectra[i]),
            color='skyblue')
            sf = 30.0
            ax2.set_ylim((1.0-(1-yl[0])/sf,1.0+(yl[1]-1)/sf))
            ax2.set_ylabel('Ratio in transit / out of transit',fontsize = 7)
            ax2.tick_params(axis='both', which='major', labelsize=6)
            ax2.tick_params(axis='both', which='minor', labelsize=5)
            ax[0][0].set_ylabel('Y (Rs)',fontsize=7)
            ax[0][0].tick_params(axis='both', which='major', labelsize=6)
            ax[0][0].tick_params(axis='both', which='minor', labelsize=5)
            ax[1][0].set_ylabel('Y (Rs)',fontsize=7)
            ax[1][0].set_xlabel('X (Rs)',fontsize=7)
            ax[1][0].tick_params(axis='both', which='major', labelsize=6)
            ax[1][0].tick_params(axis='both', which='minor', labelsize=5)
            ax[0][1].set_ylabel('Normalised flux',fontsize=7)
            ax[0][1].set_xlabel('Orbital Phase',fontsize='small')
            ax[0][1].tick_params(axis='both', which='major', labelsize=6)
            ax[0][1].tick_params(axis='both', which='minor', labelsize=5)
            ax[1][1].set_ylabel('Normalised flux',fontsize=7)
            ax[1][1].set_xlabel('Wavelength (nm)',fontsize=7)
            ax[1][1].tick_params(axis='both', which='major', labelsize=6)
            ax[1][1].tick_params(axis='both', which='minor', labelsize=5)
            if len(str(i)) == 1:
                out = '000'+str(i)
            if len(str(i)) == 2:
                out = '00'+str(i)
            if len(str(i)) == 3:
                out = '0'+str(i)
            if len(str(i)) == 4:
                out = str(i)
            fig.savefig('anim/'+out+'.png', dpi=fig.dpi)
            integrate.statusbar(i,self.Nexp)
            plt.close()
        print('',end="\r")
        print('--- Saving to animation.gif')

        status = os.system('convert -delay 8 anim/*.png animation.gif')
        if status != 0:
            raise Exception('The conversion of the animation frames into a gif has failed; '
            'probably because the Imagemagick convert command was not found. The animation frames '
            'have been created in the anim/ folder. If you want to convert these into a .gif '
            'please do it in some alternative way, or install Imagemagick, see '
            'https://imagemagick.org')



def test_StarRotator():
    from StarRotator import StarRotator
    import matplotlib.pyplot as plt
    error_trigger = 0
    KELT9 = StarRotator(586.0,592.0,50.0)
    wl = KELT9.wl
    F_out = KELT9.stellar_spectrum
    spectra = KELT9.spectra
    residuals = KELT9.residual
    KELT9.convolve_spectral_resolution()



    #Test that multiple-convolution is detected and blocked:
    try:
        KELT9.convolve_spectral_resolution()
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Trying convolution twice in a row should be caught.")
    error_trigger=0

    #Testing another grid and another vartype.
    KELT9 = StarRotator(586,592.0,13)
    KELT9.convolve_spectral_resolution()


    #Now test the whole thing with a dictionary as input. First test that the input is well
    #tested:
    in_dict = {'lala':1.0}
    try:
        KELT9 = StarRotator(586,592.0,13,input=in_dict)
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Wrong input dictionary not caught.")
    error_trigger=0


    in_dict = {'veq':114000.0,'stelinc':90.0,'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
    'u1':0.93,'u2':-0.23,'R':115000,'mus':0,'model':'PHOENIX','sma_Rs':3.153,
            'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
            'Tc':57095.68572,'mode':'phases','phases':[-0.02,-0.01,0.0,0.01,0.02]}

    KELT9 = StarRotator(586,592.0,13,input=in_dict)


    in_dict = {'veq':114000.0,'stelinc':90.0,'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
    'u1':0.93,'u2':-0.23,'R':115000,'mus':0,'model':'PHOENIX','sma_Rs':3.153,
            'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
            'Tc':57095.68572,'mode':'times','times':[5000.0,5000.1,5000.2,5000.3],
            'exptimes':[10,10,10,10]}
    KELT9 = StarRotator(586,592.0,13,input=in_dict)

    print('Tests complete.')

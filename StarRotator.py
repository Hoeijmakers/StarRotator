######################
#authors: Jens Hoeijmakers and Julia Seidel
#Description: Calculates the stellar spectrum
#
#
#####################
#Call like: python3 main.py 586.0 592.0 110000.0 90.0 90.0 0.0 0.0 500
#import statements
import sys
import numpy as np
import argparse
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
    def __init__(self,wave_start,wave_end,grid_size,star_path='demo_star.txt',planet_path='demo_planet.txt',
                 obs_path='demo_observations.txt',linelist_path=''):
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
        self.read_system(star_path=star_path,planet_path=planet_path,obs_path=obs_path)
        self.linelist_path = linelist_path
        self.compute_spectrum()

        # return(self.wlF,F_out)#This is mad return statement. This whole function should be a class instead.

    def read_system(self,star_path='demo_star.txt',planet_path='demo_planet.txt',obs_path='demo_observations.txt'):
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
        """
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
        self.grid_model = str(starparams[11].split()[0])
        self.abund = starparams[12:]
        self.sma_Rs = float(planetparams[0].split()[0])
        self.ecc = float(planetparams[1].split()[0])
        self.omega = float(planetparams[2].split()[0])
        self.orbinc = float(planetparams[3].split()[0])
        self.pob = float(planetparams[4].split()[0])#Obliquity.
        self.Rp_Rs = float(planetparams[5].split()[0])
        self.orb_p = float(planetparams[6].split()[0])
        self.transitC = float(planetparams[7].split()[0])
        self.mode = planetparams[8].split()[0]#This is a string.
        times = [] #These are in JD-24000000.0 or in orbital phase.
        self.exptimes = []
        for i in obsparams:
            times.append(float(i.split()[0]))
        self.times = np.array(times)
        self.Nexp = len(self.times)#Number of exposures.
        if self.mode == 'times':
            for i in obsparams:
                self.exptimes.append(float(i.split()[1]))
            self.exptimes = np.array(self.exptimes)
        try:
            test.typetest(self.wave_start,float,varname='wave_start in input')
            test.nantest(self.wave_start,varname='wave_start in input')
            test.notnegativetest(self.wave_start,varname='wave_start in input')
            test.notnegativetest(self.velStar,varname='velStar in input')
            test.notnegativetest(self.stelinc,varname='stelinc in input')
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
            if self.model == 'PHOENIX':
                print('--- Reading spectrum from PHOENIX')
                print('-----T=%sK, log(g)=%s, Z=%s.' % (self.T,self.logg,self.Z))
                wl,fx = spectrum.read_spectrum(self.T,self.logg,metallicity=self.Z)
            elif self.model == 'pySME': 
                print('--- Creating spectrum from pySME')
                print('-----T=%sK, log(g)=%s, Z=%s.' % (self.T,self.logg,self.Z))
                wl, fx= spectrum.get_spectrum_pysme(self.wave_start, self.wave_end, self.T, self.logg, self.Z, self.linelist_path, grid = self.grid_model)
            else:
                print('ERROR: Invalid model spectrum chosen. Input either PHOENIX or pySME')
                sys.exit()
            print('--- Integrating disk')
            if  self.drr == 0:
                print('------ Fast integration')
                wlF,F = integrate.build_spectrum_fast(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            else:
                print('------ Slow integration')
                wlF,F = integrate.build_spectrum_slow(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
        else:#Meaning, if we have no mu's do:
            if self.model == 'SPECTRUM':
                test.test_KURUCZ()
                print('--- Computing limb-resolved spectra with SPECTRUM')
                print('-----T=%sK, log(g)=%s, Z=%s.'% (self.T,self.logg,self.Z))
                # print(self.wave_start,math.floor(ops.vactoair(self.wave_start)*10.0)/10.0)
                # print(self.wave_end,math.ceil(ops.vactoair(self.wave_end)*10.0)/10.0)
                # sys.exit()
                maxvel = math.ceil(np.nanmax(np.abs(self.vel_grid)))
                wl,fx_list = spectrum.compute_spectrum(self.T,self.logg,self.Z,self.mus,math.floor(ops.vactoair(self.wave_start*ops.doppler((-1.0)*maxvel))*10.0)/10.0,math.ceil(ops.vactoair(self.wave_end*ops.doppler(maxvel))*10.0)/10.0,mode='anM')
                # print(wl,self.wave_start)
                # sys.exit()
                wlF,F = integrate.build_spectrum_fast(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)

            elif self.model == 'pySME':
                print('--- Computing limb-resolved spectra with pySME')
                print('-----T=%sK, log(g)=%s, Z=%s.'% (self.T,self.logg,self.Z))
                wl, fx_list = spectrum.get_spectrum_pysme(self.wave_start, self.wave_end, self.T, self.logg, self.Z, self.linelist_path, self.mus, self.abund, grid=self.grid_model)
                print('--- Integrating limb-resolved disk')
                wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,self.mus, self.wave_start,self.wave_end,self.x,self.y,self.vel_grid, self.flux_grid)

            else:
                print('ERROR: Invalid model spectrum chosen. Input either SPECTRUM or pySME')
                sys.exit()


        self.xp,self.yp,self.zp = ppos.calc_planet_pos(self.sma_Rs, self.ecc, self.omega, self.orbinc, self.pob, self.Rp_Rs, self.orb_p, self.transitC, self.mode, self.times, self.exptimes)



        # for i in range(len(self.xp)):
        #     print(self.xp[i],self.yp[i],self.zp[i])
        # pdb.set_trace()
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
        self.compute_residuals()

    def compute_residuals(self):
        """Compute and return the residuals of the time_series of spectra.
        This may be the primary output of StarRotator.

        Parameters
        ----------
            None

        Returns
        -------
            residual : np.array
                2D matrix corresponding in dimensions to the length of the
                time-series times the number of wavelength points, containing
                the residuals of the spectra after dividing out the out-of-
                transit flux.
        """
        import pdb
        self.residual = self.spectra*0.0
        # self.residual_smooth = self.spectra*0.0
        for i in range(self.Nexp):
            self.residual[i,:]=self.spectra[i]/self.stellar_spectrum
        # pdb.set_trace()
        # self.apply_spectral_resolution(self.R)

    def plot_residuals(self):
        """Compute and plot the residuals.

        Parameters
        ----------
            None
        Returns
        -------
            None
        """
        import matplotlib.pyplot as plt
        # res = self.residuals()
        plt.pcolormesh(self.wl,self.times,self.residual)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Phase')
        plt.show()


    def convolve_spectral_resolution(self):
        import lib.operations as ops
        import astropy.constants as const
        import copy
        import scipy.ndimage.filters as SNF
        import scipy.interpolate as interp
        from lib.integrate import statusbar as statusbar
        dv = const.c.value / self.R / 1000.0 #in km/s
        # wl_high = np.linspace(np.min(self.wl),np.max(self.wl),num=len(self.wl)*30.0)
        print('---Blurring')

        for i in range(len(self.spectra)):
            statusbar(i,len(self.spectra))
            # spec_i = interp.interp1d(self.wl,self.spectra[i])
            # self.residual_smooth[i] = ops.blur_spec(self.wl,copy.deepcopy(self.residual[i]),dv)
            # self.spectra_smooth[i] = ops.blur_spec(self.wl,copy.deepcopy(self.spectra[i]),dv)
            # self.spectra_smooth[i] = SNF.gaussian_filter(copy.deepcopy(self.spectra[i]),20.0,truncate=6.0)
            # self.spectra_smooth[i] = ops.smooth(copy.deepcopy(self.spectra[i]),40.0)
            # self.spectra_smooth[i] = interp.interp1d(wl_high,ops.blur_spec(wl_high,spec_i(wl_high),dv))(self.wl)

            res = self.residual[i]
            wl_constant_v,res_constant_v,a = ops.constant_velocity_wl_grid(self.wl,res,3.0)
            res_smooth = ops.smooth(res_constant_v,dv/a)
            self.residual[i] = interp.interp1d(wl_constant_v,res_smooth,fill_value='extrapolate',bounds_error=False)(self.wl)
            # self.residual[i] = ops.blur_spec(self.wl,copy.deepcopy(self.residual[i]),dv)
        # self.residual_save = copy.deepcopy(self.residual)
        # self.residual = copy.deepcopy(self.residual_smooth)
        # self.spectra_save = copy.deepcopy(self.spectra)
        # self.spectra = copy.deepcopy(self.spectra_smooth)

    # def undo_spectral_resolution(self):
    #     self.spectra = copy.deepcopy(self.spectra_save)
        # self.residual = copy.deepcopy(self.residual_save)

    def animate(self):
        """Plots an animation of the transit event, the stellar flux and velocity
        fields, and the resulting transit and line-shape variations. The animation
        is located in the subfolder `anim`. Anything located in this folder prior
        to running this function will be removed.

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
            ax[0][0].pcolormesh(self.x,self.y,self.flux_grid*mask,vmin=0,vmax=1.0*np.nanmax(self.flux_grid),cmap='autumn')
            ax[1][0].pcolormesh(self.x,self.y,self.vel_grid*mask,cmap='bwr')
            ax[0][0].axes.set_aspect('equal')
            ax[1][0].axes.set_aspect('equal')
            if self.zp[i] > 0.0:
                planet1 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black', edgecolor='black', lw=1)
                planet2 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black', edgecolor='black', lw=1)
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
            ax2.plot(self.wl,(self.spectra[i])*np.nanmax(F)/F/np.nanmax(self.spectra[i]),color='skyblue')
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
            print('The conversion of the animation frames into a gif has')
            print('failed; probably because the Imagemagick convert command')
            print('was not found. The animation frames have been created in')
            print('the anim/ folder. If you want to convert these into a .gif,')
            print('please do it manually, or install Imagemagick, see')
            print('https://imagemagick.org')


        #This is for plotting random comparisons.
        # wl2,fx2 = spectrum.read_spectrum(T,logg)
        # wlF2,F2 = integrate.build_spectrum_fast(wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)
        # wlp,Fp,flux,mask = integrate.build_local_spectrum_limb_resolved(-0.3,0.0,0.1,wl,fx_list,mus, wave_start, wave_end,x,y,vel_grid)
        # wlp2,Fp2,flux2,mask2 = integrate.build_local_spectrum_fast(-0.3,0.0,0.1,wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)
        ##This overplots non-rotating SPECTRUM and PHOENIX spectra, normalised.
        ## plt.plot(wl,fx_list[-1]/max(fx_list[-1]),label='SPECTRUM')
        ## plt.plot(wl2,fx2/5e15,label='PHOENIX')
        ## plt.xlabel('Wavelength (nm)')
        ## plt.ylabel('Max-normalised flux')
        ## plt.title('T = %s K, log(g) = %s' % (T,logg))
        ## plt.legend()
        ## plt.show()

        # plt.plot(wlF,F/max(F),color='skyblue',alpha=0.5)
        # plt.plot(wlF,(F-Fp)/np.nanmax(F-Fp),color='skyblue',label='SPECTRUM')
        # plt.plot(wlF2,F2/max(F2),color='red',alpha=0.5)
        # plt.plot(wlF2,(F2-Fp2)/np.nanmax(F2-Fp2),color='red',label='PHOENIX')
        # plt.legend()
        # plt.xlabel('Wavelength (nm)')
        # plt.ylabel('Max-normalised flux')
        # plt.title('T = %s K, log(g) = %s, vsini = 110km/s' % (T,logg))
        # plt.show()
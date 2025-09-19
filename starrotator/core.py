######################
#authors: Jens Hoeijmakers, Julia Seidel, Madeline Lam, Bibiana Prinoth
#Description: Calculates the stellar spectrum
#
#
#####################

#import statements
import numpy as np
import starrotator.lib.vgrid as vgrid
import starrotator.lib.plotting as pl
import starrotator.lib.operations as ops
import starrotator.lib.stellar_spectrum as spectrum
import starrotator.lib.integrate_depr as integrate_depr
import starrotator.lib.dynamics as dynamics
import starrotator.lib.util as util
from importlib.resources import files
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import copy
import os.path

import jax
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.operations import vert_int_q_ld, circ_int_q_ld, vert_int_q_ld_bounded
from starrotator.lib.util import gaussian
import starrotator.lib.util as ut
from starrotator.lib.vgrid import calc_vel_stellar, calc_flux_stellar
from starrotator.lib.integrate import sum_hidden_spectrum_v1, sum_stellar_spectrum_v1
from starrotator.lib.integrate import sum_hidden_spectrum_v2, sum_stellar_spectrum_v2, create_hidden_grid_array
#main body of code


class StarRotator(object):
    def __init__(self,wave_start,wave_end,grid_size,grid_planet_size=None,system_path=None,
    obs_path=None,input={},linelist_path=''):
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
                Number of grid cells to make up the stellar disk. Set to values greater than 400
                to limit numerical errors, or less if you are trying things out and just
                want speed.
            grid_planet_size: None, int
                Number of grid cells to make up the planet disk. If set to None this is set to be
                equal to one quarter of the grid size of the star.
            system_path : str
                Path to parameter file defining the system. This file should contain the following keywords
                and values on separate lines, and the default values are as follows:
                    50000.0     v_eq
                    90.0        stellar i
                    0.0         Differential rotation parameter (alpha)
                    5000.0      Stellar effective temperature (K)
                    0.0         Fe/H
                    4.5         logg
                    0.93        Limb-darkening coefficient u1
                    -0.23       Limb-darkening coefficient u2
                    0           Number of mu angles to consider. For values higher than zero, StarRotator switches to SPECTRUM rather than PHOENIX.
                    3.153       a/Rs
                    0.0         e
                    0.0         omega
                    86.79       Orbital inclination
                    85.0        Projected obliquity
                    0.08228     Rp/Rs
                    1.4811235   Orbital period
                    57095.68572 Transit center time - 2400000.
                    phases      mode, providing the interpretation of the timestamps of the observations:

                Note that if mode is set to pysme, the limb darkening parameters that the user provides are overridden.
                If no system path is set, the input defaults to the demo data packaged along with the package.
            
            obs_path: str
                Path to the parameter file defining the timestamps of the observations.
                Tthis file is assumed to contain a list of orbital phases. Note that if you are 
                modelling observations of long exposure times, the output of StarRotator will be less accurate 
                because the planet moves significantly during an exposures. This effect is greater for orbits 
                aligned to the stellar equator (obliquity 0 degrees).

                By default, the following phases are provided in the demo file:
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
                This can also be provided in the input dictionary (see below). 
                
            input : dict
                A dictionary of the input parameters can be provided instead of input files.
                If a dictionary and file paths are provided simultaneously, the dictionary overrides.    

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
        import pdb
        self.status = 'initialised'
        self.wave_start=float(wave_start)
        self.wave_end=float(wave_end)
        self.grid_size=int(grid_size)
        if grid_planet_size is not None:
            self.grid_planet_size = int(grid_planet_size) # Making this accessible to the user if needed.
        else:
            self.grid_planet_size = int(grid_size/4)
        self.linelist_path = linelist_path
        self.input = input
        self.system_path = system_path
        self.obs_path = obs_path
        if (self.system_path is None or self.obs_path is None) and len(self.input) == 0:
            #Meaning, if no input is provided, default to the demo data.
            self.system_path = files("starrotator.data").joinpath("demo_system.txt")
            self.obs_path = files("starrotator.data").joinpath("demo_observations.txt")


        self.read_system(system_path=self.system_path,obs_path=self.obs_path,input=self.input)


        # Test if a confifile exists in the default location:
        if ut.CONFIG_FILE.exists():
            pass #If so, anything that uses the cache directory will read from it.
            # And we won't overwrite it.
        else: # If no configfile exists, we will write a default one:
            ut.save_default_config()#This contains a path to a default
            #location where app or cache data is expected.

        self.get_stellar_spectrum()
        self.compute_orbit()
        self.compute_spectrum()

    def read_system(self,system_path=None,obs_path=None,input={}):
        """Reads in the stellar, planet and observation parameters from file; performing
        tests on the input and raising the read variables to the class.

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
                veq (equatorial velocity, km/s, float),
                stelinc (inclination of stellar rotation axis, degrees, float),
                drr (differential rotation, float),
                T (Teff, float)
                FeH (metallicity, float)
                logg (float)
                u1 (limb darkening parameter 1)
                u2 (limb darkening parameter 2)
                mus (number of mu angles, int)
                R (resolving power, float)
                model (model type, string, either PHOENIX or pySME)
                sma_Rs (semi-major axis in stellar radii, NOT SOLAR radii!, float)
                Rstar (stellar radius in solar radii)
                e (eccentricity, float)
                omega (longitude of periastron, degrees, float)
                inclination (degrees, float)
                obliquity (degrees, float)
                RpRs (planet-star radius ratio, float)
                P (orbital period in days, float)
                phases (numpy array, set to the orbital phase values of the time series)


            Setting the input dictionary overrules the input parameter files.
            If the model is set to pySME, then the following parameters also need to be set:
                grid_model (str, either atlas12.sav or marcs2012.sav are provided by default),
                abund (list, empty by default) The elements are strings of definitions of
                single-key dictionaries of the form ["{X:6.4}","{Y:6.3}"] etc.
                linelist_path (str, path to VALD-style line-list for pysme to use)
        """
        self.status = 'start reading input'
        if len(input)==0:#If we read input from config files
            input = util.read_into_dictionary(system_path)
            util.check_integrity_input(input)

            phases = [] #These are in orbital phase.
            obsparams = open(obs_path,'r').read().splitlines()
            for i in obsparams:
                phases.append(float(i.split()[0]))
            self.phases = np.array(phases)
        else:
            util.check_integrity_input(input,['phases'])
            self.phases = np.array(input['phases'])


            
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
        self.sma_Rs     = float(input['sma_Rs'])
        self.Rstar      = float(input['Rstar'])
        self.ecc        = float(input['e'])
        self.omega      = float(input['omega'])
        self.orbinc     = float(input['inclination'])
        self.pob        = float(input['obliquity'])#Obliquity.
        self.Rp_Rs      = float(input['RpRs'])
        self.orb_p      = float(input['P'])
        if 'mp' not in input:
            self.mp = 0.0
        else:
            self.mp     = float(input['mp'])

        if self.model.lower() in ['pysme','sme']:
            also_req_keys = ['grid_model','abund','linelist_path']
            util.check_integrity_input(input,also_req_keys)


            self.grid_model = str(input['grid_model'])
            self.abund      = input['abund']
            self.linelist_path = input['linelist_path']
            if not os.path.isfile(self.linelist_path):
                raise Exception("pySME linelist_path does not point to an existing file.")




        self.Nexp = len(self.phases)#Number of exposures.
        self.residual = None
        self.blurred = 0
        try:
            util.vartest(self.wave_start,varname='wave_start in input',nonans=True,pos=True)
            util.vartest(self.wave_end,varname='wave_end in input',nonans=True,pos=True)
            util.vartest(self.velStar,varname='veq in input',nonans=True,notnegative=True)
            util.vartest(self.stelinc,varname='stelinc in input',nonans=True)
            util.vartest(self.T,varname='Teff in input',nonans=True,notnegative=True)
            util.vartest(self.logg,varname='logg in input',nonans=True)
            util.vartest(self.R,varname='Resolution in input',nonans=True,pos=True)
            util.vartest(self.orb_p,varname='period in input',nonans=True,notnegative=True)
            util.vartest(self.Rp_Rs,varname='RpRs in input',nonans=True,notnegative=True)
            util.vartest(self.ecc,varname='e in input',nonans=True,notnegative=True)
            util.vartest(self.sma_Rs,varname='sma_Rs in input',nonans=True,notnegative=True)
            util.vartest(self.Rstar,varname='Rstar in input',nonans=True,pos=True)
            util.vartest(self.u1,varname='u1 in input',nonans=True)
            util.vartest(self.u2,varname='u1 in input',nonans=True)
            util.vartest(self.model,varname='model in input')
            util.vartest(self.Z,varname='Metallicity in input',nonans=True)
            util.vartest(self.pob,varname='Obliquity in input',nonans=True)
            util.vartest(self.orbinc,varname='Inclination in input',nonans=True)
            util.vartest(self.omega,varname='Omega in input',nonans=True)
            util.vartest(self.drr,varname='Alpha (drr) in input',nonans=True)
            util.vartest(self.mp,varname='mp in input',nonans=True,notnegative=True)
        except ValueError as err:
            print("Parser: ",err.args)

        if self.mus != 0:
            self.mus = np.sqrt(0.5 * (2 * np.arange(self.mus) + 1) / self.mus)
        self.status = 'success reading input'


    def compute_orbit(self):
        """
        A wrapper for calling jaxoplanet Keplerian solver.
        Output is maintained as class attributes.
        """
        xp, yp, zp = dynamics.orbit_euclidian(self.phases, 
                                                a = self.sma_Rs * self.Rstar, 
                                                m = self.mp, 
                                                P = self.orb_p, 
                                                e = self.ecc, 
                                                omega = self.omega,
                                                i=self.orbinc
                                                )
        # Convert the output from solar radii to stellar radii:
        self.xp,self.yp,self.zp = xp/self.Rstar, yp/self.Rstar, zp/self.Rstar



    def get_stellar_spectrum(self):
        """
        Obtaining the unbraodened spectrum using one of StarRotator's 
        #default methos: PHOENIX or pySME.
        The input stored in the class attributes control a logic to 
        switch between PHOENIX (mus=0) and pySME.
        """

        if isinstance(self.mus,np.ndarray) != True:
            if self.model.lower() == 'phoenix':
                print('--- Reading spectrum from PHOENIX')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl_wide,fx_wide = spectrum.load_PHOENIX(self.T,
                                               self.logg,
                                               metallicity=self.Z)
                wl = wl_wide[(wl_wide > self.wave_start) & (wl_wide < self.wave_end)]
                fx = fx_wide[(wl_wide > self.wave_start) & (wl_wide < self.wave_end)]
            elif self.model.lower() in ['pysme','sme']:
                print('--- Generating spectrum using pySME')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl, fx= spectrum.get_spectrum_pysme(self.wave_start, 
                                                    self.wave_end, 
                                                    self.T, 
                                                    self.logg, 
                                                    self.Z, 
                                                    self.linelist_path, 
                                                    grid = self.grid_model)
            else:
                raise Exception('Invalid model spectrum chosen. Input either PHOENIX or pySME in '
                'star.txt')
            self.wl_in = wl*1.0
            self.fx_in = fx*1.0

    #Define a set of wrappers to avoid having to refer to self all the time in compute_spectrum().
    def compute_grids(self,x,y,i_stellar,vel_eq,diff_rot_rate,a1,a2):
        vel_grid  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)
        flux_grid  = calc_flux_stellar(x,y,a1,a2,norm=False)
        return(flux_grid,vel_grid)
    def calc_v1(self,wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N1=400,N2=200):
        F_out_v1 =  sum_stellar_spectrum_v1(wl,fx,vel_eq,i_stellar,a1,a2,N=N1,norm=False)
        F_in_v1 = sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=N2)
        return(F_out_v1,F_in_v1)
    def calc_v2(self,wl,fx,x,y,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,vel_grid,flux_grid,batched=True,N2=200):
        dx = x[1]-x[0]
        dy = y[1]-y[0]
        F_out_v2 = sum_stellar_spectrum_v2(wl,fx,vel_grid,flux_grid,batched=batched)*dx*dy
        flux_grid_array,vel_grid_array,mu_array,dxR = create_hidden_grid_array(xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=N2)
        F_in_v2 = sum_hidden_spectrum_v2(wl,fx,vel_grid_array,flux_grid_array,batched=batched) *dxR**2 
        F_out_v2.block_until_ready()
        F_in_v2.block_until_ready()
        return(F_out_v2,F_in_v2)
    #End of wrappers.

    def compute_spectrum(self):
        """This wraps the main computation, switching between modes of 
        integration depending on the input that is stored in the set of class 
        attributes. The simulation output and other variables are raised 
        to class-wide attributes. Note that jitting is not done at this 
        level, but all computations under the hood are jitted.

        Parameters
        ----------
        None
        """
        self.status = 'start computing spectra'
        #Two arrays for the x and y axes
        self.x = jnp.linspace(-1,1,self.grid_size) #in units of stellar radius
        self.y = jnp.linspace(-1,1,self.grid_size)


        if self.model == "pySME":
            self.flux_grid,self.vel_grid = self.compute_grids(
                self.x,self.y,self.stelinc,self.velStar,self.drr,0.0,0.0)#If pySME: override LD.
        else:
            self.flux_grid,self.vel_grid = self.compute_grids(
                self.x,self.y,self.stelinc,self.velStar,self.drr,self.u1,self.u2)


        # Now we do the integration, switching between modes as input requires.
        if self.drr == 0 and self.fx_in.ndim == 1:
            # print('Calculating v1')
            self.stellar_spectrum, self.Fp = self.calc_v1(self.wl_in,
                                            self.fx_in,
                                            self.xp,self.yp,
                                            self.Rp_Rs,
                                            self.velStar,
                                            self.stelinc,
                                            self.u1,self.u2,
                                            N1=self.grid_size,
                                            N2=self.grid_planet_size
                                            )
        elif self.drr != 0 and self.fx_in.ndim == 1:
            # print('Calculating v2')
            self.stellar_spectrum, self.Fp = self.calc_v2(self.wl_in,
                                            self.fx_in,
                                            self.x,self.y,
                                            self.xp,self.yp,
                                            self.Rp_Rs,
                                            self.velStar,
                                            self.stelinc,
                                            self.drr,
                                            self.u1,self.u2,
                                            self.vel_grid,self.flux_grid,
                                            batched=True,
                                            N2=self.grid_planet_size
                                            )
            
        else:
            raise Exception("Multi-dimensional stellar spectrum (mu dependence) is not supported yet.")
        self.status = 'success computing spectra'
        # F_out = np.zeros((self.Nexp,len(F)))
        # F_planet = np.zeros((self.Nexp,len(F)))        

        # CONTINUE HERE NEXT TIME:
        # compute_stellar_spectrum_1B in integrate.py for mu-dependence.

        # That concludes the computation.
        # Operation of pysme comes after that.
        # Also need to make a decision regarding control over the wavelength axis.
        # A start and a stop wavelength is a bit silly. 
        # Need a target (data) wavelength grid onto which the result is (to be) binned-interpolated. This can be done using shone.opacity.bin_opacity().
        # And an under-the-hood model wavelength grid. That can be the PHOENIX grid simply,
        # But for pySME it has to be set to something custom I suppose. And the choice very likely matters for computation time.
        # Then documentation.


        # flux_out = []
        # mask_out = []
        # print('--- Building local spectrum')
        # for i in range(self.Nexp):
        #     if isinstance(self.mus,np.ndarray) == True:
        #         wlp,Fp,flux,mask = integrate_depr.build_local_spectrum_limb_resolved(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx_list,self.mus,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid, self.flux_grid)
        #         self.fx_list = copy.deepcopy(fx_list)
        #     else:
        #         wlp,Fp,flux,mask = integrate_depr.build_local_spectrum_fast(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
        #     integrate_depr.statusbar(i,self.Nexp)

        #     F_out[i,:]=F-Fp
        #     F_planet[i,:] = Fp
        #     flux_out.append(flux)
        #     mask_out.append(mask)
        # #This defines the output.
        # self.wl = wlF
        # self.stellar_spectrum = F
        # self.Fp = copy.deepcopy(F_planet)
        # self.spectra = copy.deepcopy(F_out)
        # self.lightcurve = np.mean(F_out, axis=1) / np.max(np.mean(F_out, axis=1))
        # self.masks = mask_out
        # self.residual = self.spectra/self.stellar_spectrum
            
    

        # MAKE IT AT HABIT TO RUN PYTEST BEFORE COMMITS!


    def compute_spectrum_depr(self):
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
        if self.model == "pySME":
            self.flux_grid = vgrid.calc_flux_stellar(self.x,self.y,0,0) #Let pysme do the limb darkening
        else:
            self.flux_grid = vgrid.calc_flux_stellar(self.x,self.y,self.u1,self.u2)
        if isinstance(self.mus,np.ndarray) != True:#SWITCH BETWEEN PHOENIX (mus=0) AND pySME
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
                print('------ Fast integration (drr not set)')
                wlF,F = integrate_depr.build_spectrum_fast(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            else:
                print('------ Slow integration (drr set)')
                wlF,F = integrate_depr.build_spectrum_slow(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
        else:
            if self.model == 'pySME':
                print('--- Computing limb-resolved spectra with pySME')
                print('-----T=%sK, log(g)=%s, Z=%s.'% (self.T,self.logg,self.Z))
                wl, fx_list = spectrum.get_spectrum_pysme(self.wave_start, self.wave_end, self.T, self.logg, self.Z, self.linelist_path, self.mus, self.abund, grid=self.grid_model)
                print('--- Integrating limb-resolved disk')
                wlF,F = integrate_depr.build_spectrum_limb_resolved(wl,fx_list,self.mus, self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            else:
                raise Exception('Invalid model spectrum chosen. Make pySME the input model.')


        self.xp,self.yp,self.zp = ppos.calc_planet_pos(self.sma_Rs, self.ecc, self.omega, self.orbinc, self.pob, self.Rp_Rs, self.orb_p, self.times, self.exptimes)


        F_out = np.zeros((self.Nexp,len(F)))
        F_planet = np.zeros((self.Nexp,len(F)))
        flux_out = []
        mask_out = []
        print('--- Building local spectrum')
        for i in range(self.Nexp):
            if isinstance(self.mus,np.ndarray) == True:
                wlp,Fp,flux,mask = integrate_depr.build_local_spectrum_limb_resolved(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx_list,self.mus,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid, self.flux_grid)
                self.fx_list = copy.deepcopy(fx_list)
            else:
                wlp,Fp,flux,mask = integrate_depr.build_local_spectrum_fast(self.xp[i],self.yp[i],self.zp[i],self.Rp_Rs,wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            integrate_depr.statusbar(i,self.Nexp)

            F_out[i,:]=F-Fp
            F_planet[i,:] = Fp
            flux_out.append(flux)
            mask_out.append(mask)
        #This defines the output.
        self.wl = wlF
        self.stellar_spectrum = F
        self.Fp = copy.deepcopy(F_planet)
        self.spectra = copy.deepcopy(F_out)
        self.lightcurve = np.mean(F_out, axis=1) / np.max(np.mean(F_out, axis=1))
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
        from starrotator.lib.integrate_depr import statusbar as statusbar
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


    def compute_flux_grid_pysme(self):
        '''Calculate the flux grid used for visualisation purposes when using pySME-generated
        spectra. Note that this flux grid is not used in any calculations.
        
        Parameters:
        -----------
        None
        
        Returns:
        --------
        flux_grid: 2d np.array
            Flux grid calculated from mean flux of pySME spectra at each calculated mu angle.
        '''
        import numpy as np
        import lib.operations as ops
        import lib.stellar_spectrum as spectrum
        import sys
        import lib.test as test

        F = 0#output

        # Calculate radii at the edge of the annuli
        rmu = np.sqrt(1 - self.mus**2)
        rlist = np.sqrt(0.5 * (rmu[:-1] ** 2 + rmu[1:] ** 2))  # area midpoints between rmu
        rlist = np.concatenate(([1], rlist))

        z,x_full,y_full = vgrid.calc_z(self.x,self.y)
        flux_grid = z*0.0

        for i in range(len(self.x)):
            for j in range(len(self.y)):
                if np.sqrt(self.x[i]**2+self.y[j]**2) <= 1.0:
                    r = np.sqrt(self.x[i]**2 + self.y[j]**2)
                    index = np.where(r < rlist)[-1][-1]
                    if self.mus[index] > 0:
                        flux_grid[j,i] = np.nanmean(self.fx_list[index])

        flux_grid /= np.nansum(flux_grid)
        return flux_grid
    
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
        import starrotator.lib.integrate_depr as integrate_depr
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
        if self.model.lower() in ['pysme','sme'] and isinstance(self.mus,np.ndarray) == True:
            flux_grid = self.compute_flux_grid_pysme()
        else:
            flux_grid = self.flux_grid
        for i in range(self.Nexp):
            mask = self.masks[i]
            fig,ax = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
            ax[0][0].pcolormesh(self.x,self.y,flux_grid*mask,vmin=0, vmax=1.0*np.nanmax(flux_grid),cmap='autumn')
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
            ax[0][0].set_xlim((min(self.x),max(self.x)))
            ax[1][0].set_xlim((min(self.x),max(self.x)))
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
            integrate_depr.statusbar(i,self.Nexp)
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





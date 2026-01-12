######################
#authors: Jens Hoeijmakers, Julia Seidel, Madeline Lam, Bibiana Prinoth
#Description: Calculates the stellar spectrum of a rotating star during a transit event.
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
from starrotator.lib.integrate import sum_hidden_spectrum_v1, sum_stellar_spectrum_v1, sum_hidden_spectrum_v1_mu, sum_stellar_spectrum_v1_mu
from starrotator.lib.integrate import sum_hidden_spectrum_v2, sum_stellar_spectrum_v2
#main body of code


class StarRotator(object):
    def __init__(self,wave_start,wave_end,grid_size,
                grid_planet_size=None,system_path=None,
                obs_path=None,input={},linelist_path=''):
        """
        Welcome to StarRotator.

        The StarRotator object contains the main functionality of StarRotator.
        Upon initialization, the exoplanet system input is read from file, and the
        model is computed. Calling this class with the minimum required input
        makes the model default to demo input parameters shipped together with this code.
        Normal use cases habe the user provide input parameters, either via parameter
        files (system_path and obs_path) or with a dictionary (input). In the case of
        pySME usage (disk-resolved spectra), a path to a pySME-compatible linelist file
        needs to be provided as well.

        Note that in PHOENIX mode, Starrotator will download PHOENIX spectra automatically
        to a cache directorly, located in the default cache directory (access this path using
        the print_intput method -- see below).

        The calculation outputs the out-of-transit stellar spectrum, in-transit time-series,
        the in-transit time-series normalised to the wavelength-mean


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
            equal to half of the grid size of the star.

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
            phases      mode, providing the interpretation of the timestamps of the observations.
            Note that if mode is set to pysme, the limb darkening parameters that the user provides are overridden.
            If no system path is set, the input defaults to the demo data packaged along with the package.
            
        obs_path: str
            Path to the parameter file defining the timestamps of the observations.
            Tthis file is assumed to contain a list of orbital phases. Note that if you are 
            modelling observations of long exposure times, the output of StarRotator will be less accurate 
            because the planet moves significantly during an exposures. This effect is greater for orbits 
            aligned to the stellar equator (obliquity 0 degrees).

        linelist_path: str
            Path to the VALD linelist used for generating a spectrum using pySME.
            This can also be provided in the input dictionary (see below). 
                
            
        input : dict
            A dictionary of the input parameters can be provided instead of input files.
            If a dictionary and file paths are provided simultaneously, the dictionary overrides.    

            
        Class methods
        -------------
        THIS IS OUTDATED
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
        THIS IS OUTDATED
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
            self.grid_planet_size = int(grid_size/2)
        self.linelist_path = linelist_path
        self.input = input
        self.system_path = system_path
        self.obs_path = obs_path
        if (self.system_path is None or self.obs_path is None) and len(self.input) == 0:
            #Meaning, if no input is provided, default to the demo data.
            self.system_path = files("starrotator.data").joinpath("demo_system.txt")
            self.obs_path = files("starrotator.data").joinpath("demo_observations.txt")
            self.linelist_path = files('starrotator.data').joinpath("demo_linelist.dat")


        self.read_system(system_path=self.system_path,obs_path=self.obs_path,input=self.input)


        # Test if a confifile exists in the default location:
        if ut.CONFIG_FILE.exists():
            pass #If so, anything that uses the cache directory will read from it.
            # And we won't overwrite it.
        else: # If no configfile exists, we will write a default one:
            ut.save_default_config()#This contains a path to a default
            #location where app or cache data is expected.


        #This is the primary cascade that gets run.
        self.get_stellar_spectrum()
        self.compute_orbit()
        self.compute_spectrum()

    def read_system(self,system_path=None,obs_path=None,linelist_path=None,input={}):
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


            Setting the input dictionary overrules all keyword-based input parameters.
            If the model is set to pySME, then the following parameters also need to be set:
                grid_model (str, either atlas12.sav or marcs2012.sav are provided by default),
                abund (list, empty by default) The elements are strings of definitions of
                single-key dictionaries of the form ["{X:6.4}","{Y:6.3}"] etc.
                linelist_path (str, path to VALD-style line-list for pysme to use)
        """
        self.status = 'start reading input'
        if len(input)==0:#If we read input from config files, put them in a dictionary after all.
            input = util.read_into_dictionary(system_path)
            input['linelist_path'] = linelist_path
            util.check_integrity_input(input)

            phases = [] #These are in orbital phase.

            if obs_path is not None:
                obsparams = open(obs_path,'r').read().splitlines()
                for i in obsparams:
                    phases.append(float(i.split()[0]))
                self.phases = np.array(phases)
            else:
                raise Exception('Obs_path needs to be set if input dictionary is not.')
        else:
            util.check_integrity_input(input,['phases'])
            self.phases = np.array(input['phases'])

        self.input = input #Raise to class to make available.
            
        self.velStar    = float(input['veq'])
        self.stelinc    = float(input['stelinc'])
        self.drr        = float(input['drr'])
        self.T          = float(input['T'])
        self.Z          = float(input['FeH'])
        self.logg       = float(input['logg'])
        self.u1         = float(input['u1'])
        self.u2         = float(input['u2'])
        self.R          = float(input['R'])
        self.N_mu        = int(input['N_mu'])
        self.model      = str(input['model'])
        self.sma_Rs     = float(input['sma_Rs'])
        self.Rstar      = float(input['Rstar'])
        self.ecc        = float(input['e'])
        self.omega      = float(input['omega'])
        self.orbinc     = float(input['inclination'])
        self.pob        = float(input['obliquity'])#Obliquity.
        self.Rp_Rs      = float(input['RpRs'])
        self.orb_p      = float(input['P'])
        self.constant_dlogl = False
        self.small_planet = True

        # Optional planet mass:
        if 'mp' not in input:
            self.mp = 0.0
        else:
            self.mp     = float(input['mp'])

        # Set small planet approximation to True to default if we're in mu-resolved case:
        if self.N_mu > 0:
            if 'small_planet' not in input:
                self.small_planet = True
            else:
                self.small_planet = bool(input['small_planet'])

        # If pysme, check and set required pysme input:
        if self.model.lower() in ['pysme','sme']:
            also_req_keys = ['grid_model','abund','linelist_path']
            # Set linelist path to demo if not provided in input dictionary.
            if 'linelist_path' not in input:
                input['linelist_path'] = files('starrotator.data').joinpath("demo_linelist.dat")
            util.check_integrity_input(input,also_req_keys)
            self.grid_model = str(input['grid_model'])
            self.abund      = input['abund']
            self.linelist_path = input['linelist_path']
            if not os.path.isfile(self.linelist_path):
                raise Exception("pySME linelist_path does not point to an existing file.")
            






        self.Nexp = len(self.phases)#Number of exposures.
        self.residual = None
        self.wl_in = None
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
            util.vartest(self.N_mu,varname='N_mu in input',notnegative=True)
        except ValueError as err:
            print("Parser: ",err.args)


        # Deal with mu input: Equidistant spacing in sqrt-space.
        if self.N_mu != 0:
            self.mu_array = np.sqrt(np.linspace(0,1,int(self.N_mu)))
        self.status = 'success reading input'


    def print_input(self):
        """
        A function to print all StarRotator inputs / configurations that determine or
        debug how this model instance will or has behaved. 
        """

        print('Paths:')
        print(f'\tSystem path: \t {self.system_path}')
        print(f'\tObs path: \t {self.obs_path}')
        print(f'\tLinelist path: \t {self.linelist_path}')

        print('')
        print('Input dictionary:')
        for i in self.input:
            print('\t'+i,'\t',self.input[i])
            

        print('')
        print('Grid sizes:')
        print('\t',f'Star: {self.grid_size}')
        print('\t',f'Planet: {self.grid_planet_size}')

        print('')
        print('Wavelength:')
        print('\t',f'Start: {self.wave_start}')
        print('\t',f'End: {self.wave_end}')
        if self.wl_in is not None:
            print('\t',f'Min: {np.min(self.wl_in)}')
            print('\t',f'Max: {np.max(self.wl_in)}')
            print('\t',f'N: {len(self.wl_in)}')

        print('')
        print('Phases:')
        print(f'\tN_phases: \t {len(self.phases)}')
        for i in self.phases:
            print(f'\t{i}')

        print('')
        print('N_mu angles:')
        print('\t',self.N_mu)

        print('')


        print('')
        print('Cache dir:')
        print('\t'+str(ut.get_cache_dir()))



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
        Obtaining the unbroadened spectrum using one of StarRotator's 
        #default methods: PHOENIX or pySME.
        The input stored in the class attributes control a logic to 
        switch between PHOENIX (mus=0) and pySME.
        """

        if self.N_mu == 0: #If we have no mu-resolved case then...
            if self.model.lower() == 'phoenix': #Either use PHOENIX
                print('--- Reading spectrum from PHOENIX')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl_wide,fx_wide = spectrum.load_PHOENIX(self.T,
                                               self.logg,
                                               metallicity=self.Z)
                wl = wl_wide[(wl_wide > self.wave_start) & (wl_wide < self.wave_end)]
                fx = fx_wide[(wl_wide > self.wave_start) & (wl_wide < self.wave_end)]
            elif self.model.lower() in ['pysme','sme']: #Or use pySME
                print('--- Generating spectrum using pySME')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl, fx= spectrum.get_spectrum_pysme(self.wave_start, 
                                                    self.wave_end, 
                                                    self.T, 
                                                    self.logg, 
                                                    self.Z, 
                                                    self.linelist_path, 
                                                    grid = self.grid_model,
                                                    abund = self.abund
                                                    )
            else:
                raise Exception('Invalid model spectrum chosen. Input either PHOENIX or pySME in '
                'star.txt')
            self.wl_in = wl*1.0
            self.fx_in = fx*1.0
        else: #If we do have the mu-resolved case, then we must be in pySME:
            if self.model.lower() in ['pysme','sme']: #Or use pySME
                print('--- Generating spectrum using pySME with mu-dependence')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl, fx= spectrum.get_spectrum_pysme(self.wave_start, 
                                                    self.wave_end, 
                                                    self.T, 
                                                    self.logg, 
                                                    self.Z, 
                                                    linelist = self.linelist_path, 
                                                    mu = self.mu_array,
                                                    grid = self.grid_model,
                                                    abund = self.abund
                                                    )
                self.wl_in = wl*1.0
                self.fx_in = fx*1.0 # Note that this is now an array.
            else:
                raise Exception("With mu-dependence, model must be set to pysme.")


    #Define a set of wrappers to avoid having to refer to self all the time in compute_spectrum().
    #And to make the switching logic below look more clear.
    #Each passes through the input to calculate both the out-of-transit stellar spectrum as well as the in-transit time-series.
    #These wrapper functions wrap everything that needs to be done. If you're ever gong to try to dissect how to do 
    #bits and pieces of this calculation in some flavour, this is where to start looking.
    def compute_grids_v1_v2(self,x,y,i_stellar,vel_eq,diff_rot_rate,a1,a2):
        # This needs to be moved into the main logic, and the dependency on x,y be fixed...
        vel_grid  = calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate)
        flux_grid  = calc_flux_stellar(x,y,a1,a2,norm=False)
        return(flux_grid,vel_grid)
    def calc_v1(self,wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,Nstar=400,Nplanet=200,constant_dlogl=False):
        F_out_v1 =  sum_stellar_spectrum_v1(wl,fx,vel_eq,i_stellar,a1,a2,N=Nstar,constant_dlogl=constant_dlogl)
        F_in_v1 = sum_hidden_spectrum_v1(wl,fx,xp,yp,Rp,vel_eq,i_stellar,a1,a2,N=Nplanet,constant_dlogl=constant_dlogl)
        F_out_v1.block_until_ready()
        F_in_v1.block_until_ready()
        return(F_out_v1,F_in_v1)
    def calc_v2(self,wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,batched=True,Nstar=400,Nplanet=200,constant_dlogl=False):
        F_out_v2 = sum_stellar_spectrum_v2(wl,fx,vel_eq,i_stellar,a1,a2,diff_rot_rate,N=Nstar,batched=batched,constant_dlogl=constant_dlogl)
        F_in_v2 = sum_hidden_spectrum_v2(wl,fx,xp,yp,Rp,vel_eq,i_stellar,diff_rot_rate,a1,a2,N=Nplanet,batched=batched,constant_dlogl=constant_dlogl) 
        F_out_v2.block_until_ready()
        F_in_v2.block_until_ready()
        return(F_out_v2,F_in_v2)
    def calc_v1_mu(self,wl,fx_array,xp,yp,Rp,vel_eq,i_stellar,mu_array,Nstar=400,Nplanet=200,small_planet=True,constant_dlogl=False):
        F_out_v1 = sum_stellar_spectrum_v1_mu(wl,fx_array,vel_eq,i_stellar,mu_array,N=Nstar,constant_dlogl=constant_dlogl)
        F_in_v1 = sum_hidden_spectrum_v1_mu(wl,fx_array,xp,yp,Rp,vel_eq,i_stellar,mu_array,N=Nplanet,small_planet=small_planet,constant_dlogl=constant_dlogl) 
        F_out_v1.block_until_ready()
        F_in_v1.block_until_ready()
        return(F_out_v1,F_in_v1)
    #End of wrappers.




    def compute_spectrum(self):
        """This does the main logical switching between modes of 
        integration depending on the input that is stored in the set of class 
        attributes, and calls the calculation wrappers. 
        The simulation output and other variables are raised 
        to class-wide attributes. Note that jitting is not done at this 
        level, but all computations under the hood are jitted for speed.

        Parameters
        ----------
        None
        """
        self.status = 'start computing spectra'
        #Two arrays for the x and y axes
        self.x = jnp.linspace(-1,1,self.grid_size) #in units of stellar radius
        self.y = jnp.linspace(-1,1,self.grid_size)


        #These grids are not meant for computation (because under the hood, computation
        #might use analytical tricks, meaning that these grids are not exactly the same
        #as what is being calculated). Instead, these are for visualization purposes.
        #Note that this logic-switching here is now doubled, and this should be included in the main logic switches below.
        # if self.model == "pySME":
        #     self.flux_grid,self.vel_grid = self.compute_grids(
        #         self.x,self.y,self.stelinc,self.velStar,self.drr,0.0,0.0)#If pySME: override LD.
        # else:
        #     self.flux_grid,self.vel_grid = self.compute_grids(
        #         self.x,self.y,self.stelinc,self.velStar,self.drr,self.u1,self.u2)


        # Now we do the integration, switching between modes as input requires.
        if self.drr == 0 and self.fx_in.ndim == 1:
            print('Calculating v1 (no drr, no mu dependence)')
            self.stellar_spectrum, self.Fp = self.calc_v1(self.wl_in,
                                            self.fx_in,
                                            self.xp,self.yp,
                                            self.Rp_Rs,
                                            self.velStar,
                                            self.stelinc,
                                            self.u1,self.u2,
                                            Nstar=self.grid_size,
                                            Nplanet=self.grid_planet_size,
                                            constant_dlogl=self.constant_dlogl
                                            )
        elif self.drr != 0 and self.fx_in.ndim == 1:
            print('Calculating v2 (with drr, no mu dependence)')
            self.stellar_spectrum, self.Fp = self.calc_v2(self.wl_in,
                                            self.fx_in,
                                            self.xp,self.yp,
                                            self.Rp_Rs,
                                            self.velStar,
                                            self.stelinc,
                                            self.drr,
                                            self.u1,self.u2,
                                            batched=True,
                                            Nstar=self.grid_size,
                                            Nplanet=self.grid_planet_size,
                                            constant_dlogl=self.constant_dlogl
                                            )
        elif self.drr == 0. and self.fx_in.ndim > 1:
            print('Calculating v1_mu (no drr, with mu dependence)')
            self.stellar_spectrum, self.Fp = self.calc_v1_mu(self.wl_in,
                                            self.fx_in,
                                            self.xp,self.yp,
                                            self.Rp_Rs,
                                            self.velStar,
                                            self.stelinc,
                                            self.mu_array,
                                            Nstar = self.grid_size,
                                            Nplanet = self.grid_planet_size,
                                            small_planet=self.small_planet,
                                            constant_dlogl=self.constant_dlogl
                                            )
        else:
            print('Calculating v3 (with drr and with mu dependence)')
            raise Exception("Multi-dimensional stellar spectrum (mu dependence) plus drr is not supported yet.")
        



        # #This defines the output.
        self.wl = self.wl_in
        self.spectra = self.stellar_spectrum-self.Fp #Spectral time series
        self.lightcurve = np.mean(self.spectra, axis=1) / np.max(np.mean(self.spectra, axis=1))
        self.masks = None #This no longer exists as output but did previously
        self.spectra_norm = (self.spectra.T / np.nanmedian(self.spectra,axis = 1)).T
        self.residual = self.spectra/self.stellar_spectrum
        self.residual_norm = self.spectra_norm/self.stellar_spectrum * np.nanmedian(self.stellar_spectrum)




        #If this point is reached, at least the functions ran through, so we change the status message.
        self.status = 'success computing spectra'
        #That does not, of course, mean that the calculation is correct/accurate, just that it finished.
   


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
            plt.pcolormesh(self.wl,self.phases,self.residual)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Orbital phase')
            plt.show()




        # CONTINUE HERE NEXT TIME:
        # [DONE] Fix normalizations and uniformity of calling v1, v2, v1_mu.
        # [DONE] Implement optional switching to fast doppler shifting method (switched by wavelength input, either array or float). Done for v1.
        # [DONE] Complete the generation of output. 
        # [DONE] Fix definition of mu. Now mu is apparently assumed to be equal to r. But it is sqrt(1-r**2). Confirm this by testing new integrate.integrate_fluxdisk_mu function.
        # Finish logic switching.
        # Add dlogl input.
        # Implement spectral resolution.
        # Complete the plotting of basic output.
        # Complete / test the workflow with pySME.
        # Create a suite of working (KELT-9) examples.

        # Also need to make a decision regarding control over the wavelength axis:
        # A start and a stop wavelength is a bit silly. 
        # Need a target (data) wavelength grid onto which the result is (to be) binned-interpolated. This can be done using shone.opacity.bin_opacity().
        # And an under-the-hood model wavelength grid. This should be dloglambda-constant, so that doppler shifting can be sped up.
        # [DONE] Make fast doppler shifting
        # Add constant-dlogl switching in all functions. See if it can be made default.

        # TO DO LATER:
        # Fix installation instructions in Readme.
        # Fix documentation. [setup is complete. need to start populating with sections. docstring autosummaries are off because most functions have docstings that are poorly formatted]
        # Drop small-planet-approximation in integrate.py for mu-dependence in hidden-spectrum.
        # Add correction-polynomials for the integrals in the various functions v1,v2 to correct for the bias induced by the pixellation of the grid (approximating a circle with a grid of squares).
        # Can be measured empirically but can probably also be evaluated analytically (though perhaps not in the pseudo-analytical limb darkening model -- or, I don't want to bother). 
        # Add Retrieval workflow that makes use of the exposed functions (not the class object).
        # Add brute-force integration with no tricks, and free velocity and flux grids as input. That is v3. A spectrum index-map (that enables mu, or other variation in the spectrum).



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





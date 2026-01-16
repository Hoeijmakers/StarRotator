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
import os
from pathlib import Path
import jax
import json
from jax import jit, lax
from jax import numpy as jnp
from functools import partial
from starrotator.lib.operations import constant_velocity_grid
from starrotator.lib.util import gaussian
import starrotator.lib.util as ut
from starrotator.lib.vgrid import calc_vel_stellar, calc_flux_stellar
from starrotator.lib.integrate import sum_hidden_spectrum_v1, sum_stellar_spectrum_v1, sum_hidden_spectrum_v1_mu, sum_stellar_spectrum_v1_mu
from starrotator.lib.integrate import sum_hidden_spectrum_v2, sum_stellar_spectrum_v2
#main body of code


class StarRotator(object):
    def __init__(self,input=None,run_computation=True):
        """
        Welcome to StarRotator.

        The StarRotator object wraps the main functionality of StarRotator for easy use.
        This is not the only way to run the StarRotator computations, but the easiest one,
        and likely sufficient when you are running forward models and making plots.

        The StarRotator object is initialized with a dictionary or JSON file as input. 
        Calling this class with empty input makes the model run a set of demo input 
        parameters shipped together with this code.

        In normal use cases, the user can switch between modelling the stellar spectrum
        based on PHOENIX models (Husser et al. 2013) or pySME (Wehrhahn et al 2019). 
        pySME enables disk-resolved spectra (mu-dependence) as well as fine-tuned control
        over elemental abundances, but requires additional input, including a path to 
        a pySME-compatible linelist file as well as a choice of stellar model.

        Note that in PHOENIX mode, Starrotator will download PHOENIX spectra automatically
        to a cache directorly, located in the default cache directory (access this path using
        the print_intput method -- see below).

        Alternatively, a user has the freedom to supply their own stellar spectrum, or list of
        mu-dependent spectra if neither PHOENIX or pySME are sufficient.

        The calculation outputs the out-of-transit stellar spectrum, in-transit time-series,
        the in-transit time-series normalised to the wavelength-mean.



        Parameters
        ----------
        input : dict
            A dictionary of the input parameters can be provided instead of input files. The input dict has
            a large number of mandatory keywords as well as optional ones.   

            
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

        self.status = 'initialised'

        if input is None:
            input = files("starrotator.data").joinpath("demo_input.json")

        if isinstance(input, dict):
            pass
        elif isinstance(input, (str, os.PathLike)):
            infile = Path(input)
            if not infile.is_file():
                raise FileExistsError(f'Inpath {str(infile)} is not set to an existing file.')
                
            with infile.open("r") as f:
                input =  json.load(f)#This creates a dict from a JSON input file.
        else:
            raise TypeError("input must be a dict or a path-like object pointing to a valid JSON configuration file.")


        self.read_input(input)
        

        # Test if a confifile exists in the default location:
        if ut.CONFIG_FILE.exists():
            pass #If so, anything that uses the cache directory will read from it.
            # And we won't overwrite it.
        else: # If no configfile exists, we will write a default one:
            ut.save_default_config()#This contains a path to a default
            #location where app or cache data is expected.

        self.get_stellar_spectrum()
        self.compute_orbit()


        #This is the primary cascade that gets run.
        if run_computation:
            self.compute_spectrum()

        # That is the end of the simulation. Everything below here is the definition of class functions.
        # Including ones that expose output plots.




    def read_input(self,input={}):
        """Reads the input, performing tests on the input and raising the read variables to the class.

        Parameters
        ----------
            input: dict
                A dictionary with input parameters, instead of using textfiles as input.
                This allows programmatic control over all StarRotator inputs. 
                If the model is set to pySME, then the following parameters also need to be set:
                grid_model (str, either atlas12.sav or marcs2012.sav are provided by default),
                abund (list, empty by default) The elements are strings of definitions of
                single-key dictionaries of the form ["{X:6.4}","{Y:6.3}"] etc.
                linelist_path (str, path to VALD-style line-list for pysme to use)
        """
        self.status = 'start reading input'

        #All input needs to be set explicitly, even if sensible defaults could be defined.
        #The reason is that I want to expose the input to the user to create up-front understanding/
        #expectations of what is happening under the hood.
        mandatory_keywords = {'wavelength':[[list,'array-like'],['type']],
                              'wavelength_type':[[str],['type']],
                              'phases':[[list,'array-like'],['type']],
                              'grid_size_star':[[int],['type','pos','nonans']],
                              'grid_size_planet':[[int],['type','pos','nonans']],
                              'sma_Rs':[[float],['type','pos','nonans']],
                              'e':[[int,float],['type','notnegative','nonans']],
                              'inclination':[[int,float],['type','nonans']],
                              'obliquity':[[int,float],['type','nonans']],
                              'RpRs':[[int,float],['type','notnegative']],
                              'Rstar':[[int,float],['type','pos']],
                              'P':[[int,float],['type','pos','nonans']],
                              'mp':[[int,float],['type','notnegative','nonans']],
                              'veq':[[int,float],['type','notnegative','nonans']],
                              'stelinc':[[int,float],['type','nonans']],
                              'T':[[int,float],['type','pos','nonans']],
                              'FeH':[[int,float],['type','nonans']],
                              'logg':[[int,float],['type','nonans']],
                              'drr':[[int,float],['type','notnegative','nonans']],
                              'u1':[[int,float],['type','nonans']],
                              'u2':[[int,float],['type','nonans']],
                              'R':[[int,float],['type','nonans']],
                              'model':[[str],['type']]
                              }
        
        util.check_integrity_input(input,mandatory_keywords)
        #If this is met, then we can continue checking cases: e and pySME

        # Check that the wavelength axis is correctly specified.
        if input['wavelength_type'] not in ['explicit','linear','minmax','constant_dlogl']:
            raise Exception(f"wavelength_type in input needs to be set to 'explicit','linear','minmax' or 'constant_dlogl' ({input['wavelength_type']})")
        if input['wavelength_type'] == 'linear':
            if len(input['wavelength']) != 3:
                raise Exception(f"Length of wavelength input array is required to be equal to 3 (min, max, number of steps).")
        if input['wavelength_type'] == 'minmax':
            if len(input['wavelength']) != 2:
                raise Exception(f"Length of wavelength input array is required to be equal to 2 (min, max).")            
        if input['wavelength_type'] == 'constant_dlogl':
            if len(input['wavelength']) != 3:
                raise Exception(f"Length of wavelength input array is required to be equal to 3 (min, max, delta-logl).")          
        # The choice of wavelength axis is going to be passed to where the stellar spectrum is loaded/determined.
        # self.constant_dlogl is set in get_stellar_spectrum.


        #The following are only to be set if certain criteria are met.
        #If in pySME mode, then we require all the pySME input:
        if input['model'].lower() in ['pysme','sme']:
            additional_keywords = {'N_mu':[[int],['type','notnegative']],
                            'grid_model':[[str],['type']],
                            'abund':[[dict],['type']]
                            } #Note that the linelist path and 
                            # small-planet approximation are checked 
                            # for separately, below.
            util.check_integrity_input(input,additional_keywords)
            if int(input['N_mu']) > 0:
                if 'small_planet' not in input: #Set to default....
                    input['small_planet'] = True
                else:#... or check that it is OK.
                    util.check_integrity_input(input,{'small_planet':[[bool],['type']]})
            if 'linelist_path' not in input: #Set to default...
                demopath = files('starrotator.data').joinpath("demo_linelist.dat")
                util.file_exists(demopath)
                input['linelist_path'] = demopath
            else:#.... or check that it is OK.
                util.file_exists(input['linelist_path'])
        if input['e'] > 0:
            additional_keywords = {'omega':[[int,float],['type']]}
            util.check_integrity_input(input,additional_keywords)

        #Still good? Then we have finished our checks and we can raise everything
        #to the class-level.

        self.input = input #Raise to class to make available.
        self.wavelength_type = input['wavelength_type']
        self.wavelength = jnp.array(input['wavelength'])
        # self.wave_start = np.min(self.wavelength)
        # self.wave_end   = np.max(self.wavelength)      
        self.phases     = jnp.array(input['phases'])
        self.Nexp       = len(self.phases)#Number of exposures.
        self.grid_size  = int(input['grid_size_star'])
        self.grid_planet_size = int(input['grid_size_planet'])    
        self.sma_Rs     = float(input['sma_Rs'])    
        self.ecc        = float(input['e'])
        self.orbinc     = float(input['inclination'])
        self.pob        = float(input['obliquity'])#Obliquity.
        self.Rp_Rs      = float(input['RpRs'])
        self.Rstar      = float(input['Rstar'])
        self.orb_p      = float(input['P'])
        self.mp         = float(input['mp'])
        self.velStar    = float(input['veq'])
        self.stelinc    = float(input['stelinc'])
        self.T          = float(input['T'])
        self.Z          = float(input['FeH'])
        self.logg       = float(input['logg'])
        self.drr        = float(input['drr'])
        self.u1         = float(input['u1'])
        self.u2         = float(input['u2'])
        self.R          = float(input['R'])
        self.model      = str(input['model'])
        #These were all the mandatory ones.
        if 'omega' in input: 
            self.omega      = float(input['omega'])
        else:
            self.omega = 0.0


        # If pysme, check and set required pysme input:
        if self.model.lower() in ['pysme','sme']:
            # Set linelist path to demo if not provided in input dictionary.
            self.grid_model = str(input['grid_model'])
            self.abund      = input['abund']
            self.linelist_path = input['linelist_path']
            self.N_mu        = int(input['N_mu'])
            # Deal with mu input: Equidistant spacing in sqrt-space.
            self.mu_array = np.sqrt(np.linspace(0,1,int(self.N_mu)))
        else:
            self.N_mu = 0.0 #This variable needs to be defined because elsewhere I use it to 
            #distinguish the model case.

        self.residual = None
        self.wl_in = None
        self.blurred = 0
        self.status = 'success reading input'


    def print_input(self):
        """
        A function to print all StarRotator inputs / configurations that determine or
        debug how this model instance will or has behaved. 
        """


        print('')
        print('Input dictionary:')
        for i in self.input:
            if i not in ['wavelength','phases']:
                print('\t'+i,'\t',self.input[i]) # To keep it readable.
            

        print('')
        print('Grid sizes:')
        print('\t',f'Star: {self.grid_size}')
        print('\t',f'Planet: {self.grid_planet_size}')

        print('')
        print('Wavelength:')
        print('\t',f'Type: {self.wavelength_type}')
        print('\t',f'Coding: {self.wavelength}')
        if self.wl_in is not None:
            print('\t',f'Min: {np.min(self.wl_model)}')
            print('\t',f'Max: {np.max(self.wl_model)}')
            print('\t',f'N: {len(self.wl_model)}')
            if self.wavelength_type == 'constant_dlogl':
                print('\t',f'dlogl: {self.wl_in}')
            print(f'\t','To access wavelengths, print class.wl_model manually.')
        print('')
        print('Phases:')
        print(f'\t',f'N_phases: \t {len(self.phases)}')
        print(f'\t',f'To access phases, print class.phases manually.')

        if self.model.lower() in ['pysme','sme']:
            print('In pySME mode.')
            print(f'\tLinelist path: \t {self.linelist_path}')
            print(f'\t',f'N_mu angles: {self.N_mu}')

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
        switch between PHOENIX (mus=0) and pySME, and the definition
        of the wavelength axis. PHOENIX models are interpolated onto this
        axis when using explicit, linear or constant-dlogl grids (so be careful that
        sufficient resolution is available to avoid undersampling).
        Setting the wavelength type to min-max preserves the PHOENIX wavelength
        grid and uses a default pySME sampling.
        """


        if self.wavelength_type == 'explicit':
            wl_input = self.wavelength # In principle, the input to the spectrum functions is the wavelength axis directly.
            self.constant_dlogl = False
        if self.wavelength_type == 'minmax':
            # g = dynamics.doppler_factor(self.velStar) # Dont bother to compensate the edges. Up to the user (more predictable output).
            wl_input = self.wavelength #[np.min(self.wavelength)/g, np.max(self.wavelength)*g]
            self.constant_dlogl = False
        if self.wavelength_type == 'linear':
            # g = dynamics.doppler_factor(self.velStar)
            wl_input = np.linspace(*self.wavelength)# np.linspace(self.wavelength[0]/g, self.wavelength[1]*g,self.wavelength[2])
            #self.wavelength is coded as min,max,N, just as linspace.
            self.constant_dlogl = False
        if self.wavelength_type == 'constant_dlogl':
            wl_input,dlogl = constant_velocity_grid(*self.wavelength)
            self.constant_dlogl = True
            
        self.wl_requested_for_model = wl_input*1.0


        # To clear up what's happening with the wavelength input here and the various wavelength parameters:
        # The output of all of this is the wl_in and fx_in parameters, as they form the 
        # input to the computation later.

        # wl_in is either set to the wavelength axis that is coming out of the model, that is wl_model.
        # Or it is set to dlogl (scalar) in case constant_dlogl is requested. 
        # The wavelength output is also retrained in this case, in self.wl_model, but this is only
        # relevant for the case of constant_dlogl. In other cases, wl_in and wl_model are the same.

        # Now, wl_model is either equal to wl_input, unless minmax is requested, in which case it is 
        # equal to the PHOENIX model grid, cropped to min-max.
        # And wl_input is either an array with the wavelength samples, and then set to be equal to
        # wl_model, or it encoded minmax, in which case wl_model was set to the output of PHOENIX or pySME at minimax.

        if self.model.lower() == 'phoenix': #Either use PHOENIX
            print(f'--- Reading spectrum from PHOENIX in wavelength mode {self.wavelength_type}')
            print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
            wl_wide,fx_wide = spectrum.load_PHOENIX(self.T,
                                               self.logg,
                                               metallicity=self.Z)
            #This model is defined on a very wide wavelength range, from ~100nm to 5um.
            if self.wavelength_type == 'minmax':
                wl_model = wl_wide[(wl_wide > np.min(wl_input)) & (wl_wide < np.max(wl_input))]
                fx_model = fx_wide[(wl_wide > np.min(wl_input)) & (wl_wide < np.max(wl_input))]
            else:
                wl_model = wl_input
                fx_model = np.interp(wl_input,wl_wide,fx_wide)



        elif self.model.lower() in ['pysme','sme']: #Or use pySME
            if self.N_mu == 0:
                print(f'--- Generating spectrum using pySME without mu-dependence in wavelength mode {self.wavelength_type}.')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl, fx= spectrum.get_spectrum_pysme( 
                                                    wl_input, 
                                                    self.T, 
                                                    self.logg, 
                                                    self.Z, 
                                                    self.linelist_path, 
                                                    grid = self.grid_model,
                                                    abund = self.abund
                                                    )
            else:
                print(f'--- Generating spectrum using pySME with mu-dependence in wavelength mode {self.wavelength_type}.')
                print(f'-----T={self.T}K, log(g)={self.logg}, Z={self.Z}.')
                wl, fx= spectrum.get_spectrum_pysme(wl_input,
                                                    self.T, 
                                                    self.logg, 
                                                    self.Z, 
                                                    linelist = self.linelist_path, 
                                                    mu = self.mu_array,
                                                    grid = self.grid_model,
                                                    abund = self.abund
                                                    )
        else:
            raise Exception('Invalid model spectrum chosen. Choose either PHOENIX or pySME as spectral model.')
        
        
        self.wl_model = wl_model*1.0 #Retain the model output wavelength in case dlogl is set. 
        if self.wavelength_type == 'constant_dlogl':
            self.wl_in = dlogl
        else:
            self.wl_in = wl_model*1.0
        
        self.fx_in = fx_model*1.0




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
        # [DONE] Finish logic switching.
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





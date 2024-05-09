# This file contains the three fitting algorithms for StarRotator to fit the RM/DS.
#
# Algorithm 1: CCF-space
# Probably the most intuitive one is to fit the Doppler shadow trace in your CCF. 
# This algorithm assumes that the average stellar line of a given species is well approximated
# as a Gaussian function. 
#
# Algorithm 2: Single line in the rest frame of the star, NOT PLANET
# In general, you will not be able to see your absorption line in a time-resolved manner,
# at least not yet. So therefore we provide a fitting algorithm to fit the RM in the planetary restframe. 
#
# Algorithm 3: Time-resolved single line in the restframe of the star, NOT PLANET
# You belong to the lucky people with high SNR data? Well, look no further, this the generalisation of
# Algorithm 2 to the time-resolved component. 
#
#
#
import numpy as np
from pymultinest.solve import solve
import math
from scipy.special import erfcinv
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import astropy.constants as const

# The full magic of log-likelihoods with Gaussian probability distr. 
# Trust me, the funny way of accessing np functions and multiplications is 
# for speed. Julia tested this!


from numpy import pi as npi
from numpy import sqrt as npsqrt
from numpy import log as nplog
from numpy import fabs as npfabs
from numpy import sum as npsum

# These define the possible priors in the proper format for pymultinest
# You often find yourself to go for uniform priors :) 

def log_uniform_prior(rand,lower_upper):
    log10_lower = math.log10(lower_upper[0])
    log10_upper = math.log10(lower_upper[1])
    return 10.0**(log10_lower + rand*(log10_upper-log10_lower))

def uniform_prior(rand,lower_upper):
    return lower_upper[0] + rand * (lower_upper[1]-lower_upper[0])

def gaussian_prior(rand,mu_sigma):
    return mu_sigma[0] + mu_sigma[1] * math.sqrt(2) * erfcinv(2 * (1 - rand))

def var_prior(rand,output):
    return output

def compute_likelihood(data,errdata,model):
    sqrt_2Pi = npsqrt(2.*npi)
    diff_model_data = (model - data) * (model - data)
    result = npsum(-nplog(npfabs(errdata)*sqrt_2Pi)-diff_model_data/(2.*errdata*errdata))
    return result

# The following is all StarRotator related!

def remove_output(outpath):
    import os
    import glob
    for filename in glob.glob(outpath+'*'):
        print(f"--- [INFO] removing {filename}")
        os.remove(filename)
        
    return 0

def write_json(outpath, parameters):
    import os
    
    line = '['
    for i in range(len(parameters) - 1):
        line = line + '"' + parameters[i] +'",'
    
    line = line + '"' + parameters[-1] +'"]'
    
    if os.path.exists(outpath+'params.json'):
        print(f"--- [INFO] {outpath}params.json already exists. Writing a new one for you ")
        os.remove(outpath+'params.json')
    
    else:
        print(f"--- [INFO] Writing {outpath}params.json ")
    with open(outpath+'params.json', 'w') as f:
        f.write(line)  # python will convert \n to os.linesep

    f.close()
    
    return 0

def multinest_marginals_corner(prefix, save=False):
    import numpy
    from numpy import exp, log
    import matplotlib.pyplot as plt
    import sys, os
    import json
    import pymultinest
    import corner
    
    # #plt.rcParams.update({'font.size': 16})
    # plt.rc('font', size=12)          # controls default text sizes
    # #plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    # plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=15)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=16)    # legend fontsize
    # #plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    # plt.rcParams["xtick.direction"] = "in"
    # plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["ytick.major.width"] = 2
    # plt.rcParams["xtick.major.width"] = 2
    # plt.rcParams["xtick.major.size"] = 5
    # plt.rcParams["ytick.major.size"] = 5
    # plt.rcParams["xtick.minor.size"] = 3.5
    # plt.rcParams["ytick.minor.size"] = 3.5
    # #print(plt.rcParams.keys())
    # plt.rcParams["ytick.minor.width"] = 1
    # plt.rcParams["xtick.minor.width"] = 1
    # plt.rcParams['axes.linewidth'] = 2
    # plt.rcParams['legend.facecolor'] = 'white'
    # plt.rcParams['legend.edgecolor'] = 'white'
    # plt.rcParams['legend.framealpha'] = 1
    
    """
    Adapted from https://github.com/JohannesBuchner/PyMultiNest
    
    """

    print('model "%s"' % prefix)
    if not os.path.exists(prefix + 'params.json'):
        sys.stderr.write("""Expected the file %sparams.json with the parameter names.
    For example, for a three-dimensional problem:

    ["Redshift $z$", "my parameter 2", "A"]
    %s""" % (sys.argv[1], __doc__))
        sys.exit(2)
    parameters = json.load(open(prefix + 'params.json'))
    n_params = len(parameters)

    a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
    s = a.get_stats()

    json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)

    print('  marginal likelihood:')
    print('    ln Z = %.1f +- %.1f' % (s['global evidence'], s['global evidence error']))
    print('  parameters:')
    for p, m in zip(parameters, s['marginals']):
        lo, hi = m['1sigma']
        med = m['median']
        sigma = (hi - lo) / 2
        if sigma == 0:
            i = 3
        else:
            i = max(0, int(-numpy.floor(numpy.log10(sigma))) + 1)
        fmt = '%%.%df' % i
        fmts = '\t'.join(['    %-15s' + fmt + " +- " + fmt])
        print(fmts % (p, med, sigma))

    print('creating marginal plot ...')
    data = a.get_data()[:,2:]
    weights = a.get_data()[:,0]

    #mask = weights.cumsum() > 1e-5
    mask = weights > 1e-4

    corner.corner(data[mask,:], weights=weights[mask], 
        labels=parameters, show_titles=True)
    
    if save:
        plt.savefig(prefix + 'corner.pdf')
        plt.savefig(prefix + 'corner.png')
    
    plt.show()
    plt.close()


def generate_flux_grid(u1, u2, gridsize=25):
    # This is the function generating the flux grid for your retrieval in a minimal way.
    
    x = np.linspace(-1,1,2*gridsize) 
    y = np.linspace(-1,1,2*gridsize) # coordinate grid on the stellar surface, normalised to the stellar radius #one could move that out of the 

     #pre calculate matrices
    xx, yy = np.meshgrid(x, y)
    arg = 1-xx**2+yy**2
    mask = np.where(arg>0, 1, 0.)
    z = np.sqrt(mask*arg)
    
    mask = np.where((xx**2+yy**2) <= 1., 1., 0.)
    
    psi = np.arcsin(np.sqrt(mask*(xx**2+yy**2)))
    a0 = 1-u1-u2
    I = a0+u1*np.cos(psi)+u2*np.cos(psi)**2
    
    flux_grid = np.where((xx**2+yy**2) <= 1., I, 0.) # if not inside star, we just set it equal to one, because...                      
    flux_grid /= np.sum(flux_grid) # normalisation? limb darkening is constant, so we can do it outside!!
    
    return xx, yy, flux_grid
    
def doppler_shadow(parameters, radial_velocity, phases, xx, yy, flux_grid, orbinc, aRs, vsys, ecc=0.0, omega=0.0):
    #model function
    # ecc omega are currently 0. New version for non-circular orbits is coming soon.
    # parameters for model fit    
    pob, RpRs, vsini, amp, width = parameters
    #width = width*1000 # unit conversion!
    
    exp = len(phases) # number of exposures is equal to the number of phases
    alpha = np.radians(pob)
    
    # velocity grid
    gridsize = int(len(flux_grid)/2)
    vel_grid = xx*vsini*1000 #+ vsys*1000 #*(1.-drr*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*yy))
    v_avg = np.median(vel_grid,axis = 0) # this is fine, because it averages over the same projected velocities
    
    v = np.where(np.isnan(v_avg) == False, v_avg, 0) # we don't care about NaNs!
    
    # Now comes the shifting and adding part
    radial_velocity_shifted = v[:, np.newaxis] - radial_velocity*1000 + vsys*1000
    
    f = 1. - amp * np.exp(-0.5 * radial_velocity_shifted**2 / (width*1000)**2)
    
    ## next we are going to calculate the position of the planet on the orbit, for every phase. Wish me luck.
    # xp, yp, zp = calc_planet_pos(sma_Rs, ecc, omega, orbinc, pob, Rp_Rs, orb_p, times)
    
    inclin_bar = np.radians(orbinc)
    true_anom = 2.*np.pi*phases # CAREFUL, THIS ASSUMES A CIRCULAR ORBIT
    
    x_pl = aRs*np.sin(np.asarray(true_anom))
    y_pl = -aRs*np.cos(np.asarray(true_anom))*np.cos(inclin_bar)
    z_pl = aRs*np.cos(np.asarray(true_anom))*np.sin(inclin_bar)
    
    xp = x_pl*np.cos(np.radians(pob))-y_pl*np.sin(np.radians(pob))
    yp = x_pl*np.sin(np.radians(pob))+y_pl*np.cos(np.radians(pob))
    zp = z_pl
    
    
    ## Now we calculate the stellar flux behind the planet for each position on the orbit.
    
    xx_tile = np.reshape(np.tile(xx.flatten(), exp), (exp, 2*gridsize, 2*gridsize))
    yy_tile = np.reshape(np.tile(yy.flatten(), exp), (exp, 2*gridsize, 2*gridsize))

    x_full = (xx_tile.T - xp).T
    y_full = (yy_tile.T - yp).T
    
    
    planet_in_pixel = np.sqrt(x_full**2 + y_full**2) < RpRs # shape: (exp, 2*n, 2*n)
    planet_may_be_in_transit = np.broadcast_to(zp[:, None, None] > 0, x_full.shape)
    
    flux_weighting_with_occulting_planet = np.where(planet_in_pixel & planet_may_be_in_transit, flux_grid, 0)
    flux_weighting_no_planet = np.broadcast_to(flux_grid[None, ...], flux_weighting_with_occulting_planet.shape)
    
    #f = shift_and_coadd(radial_velocity_shifted, radial_velocity, stellar_flux) # this is the whole magic!
    
    F = (f.T @ np.sum(flux_weighting_no_planet, axis=1).T).T
    Fp = (f.T @ np.sum(flux_weighting_with_occulting_planet, axis=1).T).T
    
    F_out = F - Fp
    residuals = F_out / F
    
    residuals = residuals.T/np.mean(residuals, axis=1)

    return(residuals.T)


def rossiter_mclaughlin_2D(p, wl, phi, transit, wlS_wide, Fs_wide, xx, yy, flux_grid, orbinc, aRs, ecc=0.0, omega=0.0):
    pob, RpRs, vsini = p
    
    # planetary trace
    #x0 = xc+dxc*np.sin(2.0*np.pi*phi)*np.sin(np.radians(orbinc)) 
    #first_component = (amplitude / 1e3) * np.exp(-0.5 * (wl.T - x0)**2 / lw**2) 
    
    
    # Now comes the hard bit. We want to compute the stellar spectrum...
    gridsize = int(len(flux_grid)/2)
    exp = len(phi) # number of exposures is equal to the number of phases

    
    # VELOCITY GRID
    vel_grid = xx*vsini*1000 #+ vsys*1000 #*(1.-drr*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*yy))
        
    # OKAY, FLUX AND VELOCITY GRID ARE DONE. NEXT UP IS BUILDING THE SPECTRUM OF THE STAR
    
    flux_profile_star = np.sum(flux_grid, axis = 0) #This is the sum of the flux grid at the same velocities.
    v_avg = np.median(vel_grid,axis = 0) # this is fine, because it averages over the same projected velocities
    v = np.where(np.isnan(v_avg) == False, v_avg, 0)
    
    ## Now comes a linear interpolation for len(v) velocity steps Doppler-shifted
    
    c = const.c.value # Comes out in m/s.
    beta = v/c
    fac = 1 + beta #np.sqrt((1+beta)/(1-beta))
    wl_shifted = fac[:, np.newaxis] *  wlS_wide # doppler shifted wavelength each row is for one of the fluxes
    f = ((0*fac[:, np.newaxis]+1.) * Fs_wide) * 0.0 # pre-generating the size!
    
    # interpolation onto original wavelength grid
    for i in range(len(v)):
        f[i] = interp1d(wl_shifted[i], Fs_wide, bounds_error=None, fill_value='extrapolate')(wlS_wide)
    
    ## next we are going to calculate the position of the planet on the orbit, for every phase. Wish me luck.
    # xp, yp, zp = calc_planet_pos(sma_Rs, ecc, omega, orbinc, pob, Rp_Rs, orb_p, times)
    
    inclin_bar = np.radians(orbinc)
    
    true_anom = 2.*np.pi*phi # CAREFUL, THIS ASSUMES A CIRCULAR ORBIT
    
    x_pl = aRs*np.sin(np.asarray(true_anom))
    y_pl = -aRs*np.cos(np.asarray(true_anom))*np.cos(inclin_bar)
    z_pl = aRs*np.cos(np.asarray(true_anom))*np.sin(inclin_bar)
    
    xp = x_pl*np.cos(np.radians(pob))-y_pl*np.sin(np.radians(pob))
    yp = x_pl*np.sin(np.radians(pob))+y_pl*np.cos(np.radians(pob))
    zp = z_pl
    
    #import pdb
    #pdb.set_trace()
    
    ## Now we calculate the stellar flux behind the planet for each position on the orbit.
    
    xx_tile = np.reshape(np.tile(xx.flatten(), exp), (exp, 2*gridsize, 2*gridsize))
    yy_tile = np.reshape(np.tile(yy.flatten(), exp), (exp, 2*gridsize, 2*gridsize))

    x_full = (xx_tile.T - xp).T
    y_full = (yy_tile.T - yp).T
    
    
    planet_in_pixel = np.sqrt(x_full**2 + y_full**2) < RpRs # shape: (exp, 2*n, 2*n)
    planet_may_be_in_transit = np.broadcast_to(zp[:, None, None] > 0, x_full.shape)
    
    flux_weighting_with_occulting_planet = np.where(planet_in_pixel & planet_may_be_in_transit, flux_grid, 0)
    flux_weighting_no_planet = np.broadcast_to(flux_grid[None, ...], flux_weighting_with_occulting_planet.shape)
    
    
    #f = shift_and_coadd(wl_shifted, v, wlS_wide, Fs_wide)
    
    F = (f.T @ np.sum(flux_weighting_no_planet, axis=1).T).T
    Fp = (f.T @ np.sum(flux_weighting_with_occulting_planet, axis=1).T).T
    
    F_out = F - Fp
    residuals = F_out / F
    
    residuals = residuals.T/np.mean(residuals, axis=1)
    
    #transit_new = np.where(transit != 0., 1, 0.)
    
    #CCF_model = transit_new * first_component + residuals

    return(residuals.T)


def rossiter_mclaughlin_1D(p, wl, phi, transit, wlS_wide, Fs_wide, xx, yy, flux_grid, orbinc, aRs, ecc=0.0, omega=0.0):
    
    residual = rossiter_mclaughlin_2D()
    residual_time_averaged = np.nanmean(residual, axis=0)
    
    return residual_time_averaged


class DopplerShadow:
    
    def __init__(self,logL,prior,ndim,init_limits,data,errdata,radial_velocity,phases,xx,yy,flux_grid,orbinc,aRs,vsys):
        self.logL = logL
        self.prior = prior
        self.ndim = ndim
        self.data = data
        self.errdata = errdata
        self.init_limits = init_limits
        
        # these are the constant parameters
        self.radial_velocity = radial_velocity # unit conversion
        self.phases = phases
        self.xx = xx
        self.yy = yy
        self.flux_grid = flux_grid
        self.orbinc = orbinc
        self.aRs = aRs
        self.vsys = vsys #unit conversion    
        
    def cpte_model(self,param):
        
        model = doppler_shadow(parameters=param, # param are the fitting parameters!
                    radial_velocity=self.radial_velocity, 
                    phases=self.phases, 
                    xx=self.xx,
                    yy=self.yy,
                    flux_grid=self.flux_grid,
                    orbinc=self.orbinc, 
                    aRs=self.aRs,
                    vsys=self.vsys)
  
        
        return model

    def cpte_prior(self,unit_cube):
        pr = []
        limits = self.init_limits
        
        for k,i in enumerate(unit_cube):
            pr.append(self.prior(i,limits[k]))
        return pr
    
    def cpte_logL(self,param):
        if np.any(param == -np.inf) or np.any(param == np.inf):
            print("In SampleLine.py, MyModel(): Inf encountered.")
            return 0.
        
        model = self.cpte_model(param)
        like = self.logL(self.data,self.errdata,model)
        
        return like
    


class RM2D:
    
    def __init__(self,logL,prior,ndim,init_limits,data,errdata,wl,phases,xx,yy,flux_grid,orbinc,aRs,wlS,fxS, transit):
        self.logL = logL
        self.prior = prior
        self.ndim = ndim
        self.init_limits = init_limits
        self.data = data
        self.errdata = errdata
        
        # these are the constant parameters
        self.wavelength = wl # unit conversion
        self.phases = phases
        
        self.xx = xx
        self.yy = yy
        self.flux_grid = flux_grid
        self.orbinc = orbinc
        self.aRs = aRs
        self.wlS = wlS 
        self.fxS = fxS
        self.transit = transit
    
    def cpte_model(self,param):
        # p, wl, phi, transit, wlS_wide, Fs_wide, xx, yy, flux_grid, orbinc, aRs, ecc=0.0, omega=0.0
        
        model = rossiter_mclaughlin_2D(p=param, # param are the fitting parameters!
                             wl=self.wavelength, 
                             phi=self.phases, 
                             transit=self.transit,
                             wlS_wide=self.wlS,
                             Fs_wide=self.fxS, 
                             xx=self.xx,
                             yy=self.yy,
                             flux_grid=self.flux_grid,
                             orbinc=self.orbinc, 
                             aRs=self.aRs
                            )

        
        return model

    def cpte_prior(self,unit_cube):
        pr = []
        limits = self.init_limits
        
        for k,i in enumerate(unit_cube):
            pr.append(self.prior(i,limits[k]))
        return pr
    
    def cpte_logL(self,param):
        if np.any(param == -np.inf) or np.any(param == np.inf):
            print("In SampleLine.py, MyModel(): Inf encountered.")
            return 0.
        
        model = self.cpte_model(param)
        like = self.logL(self.data,self.errdata,model)
        
        return like
    
    
    
    
class RM1D:
    
    def __init__(self,logL,prior,ndim,init_limits,data,errdata,wl,phases,xx,yy,flux_grid,orbinc,aRs,wlS,fxS, transit):
        self.logL = logL
        self.prior = prior
        self.ndim = ndim
        self.init_limits = init_limits
        self.data = data
        self.errdata = errdata
        
        # these are the constant parameters
        self.wavelength = wl # unit conversion
        self.phases = phases
        
        self.xx = xx
        self.yy = yy
        self.flux_grid = flux_grid
        self.orbinc = orbinc
        self.aRs = aRs
        self.wlS = wlS 
        self.fxS = fxS
        self.transit = transit
    
    def cpte_model(self,param):
        # p, wl, phi, transit, wlS_wide, Fs_wide, xx, yy, flux_grid, orbinc, aRs, ecc=0.0, omega=0.0
        
        model = rossiter_mclaughlin_1D(p=param, # param are the fitting parameters!
                             wl=self.wavelength, 
                             phi=self.phases, 
                             transit=self.transit,
                             wlS_wide=self.wlS,
                             Fs_wide=self.fxS, 
                             xx=self.xx,
                             yy=self.yy,
                             flux_grid=self.flux_grid,
                             orbinc=self.orbinc, 
                             aRs=self.aRs
                            )

        
        return model

    def cpte_prior(self,unit_cube):
        pr = []
        limits = self.init_limits
        
        for k,i in enumerate(unit_cube):
            pr.append(self.prior(i,limits[k]))
        return pr
    
    def cpte_logL(self,param):
        if np.any(param == -np.inf) or np.any(param == np.inf):
            print("In SampleLine.py, MyModel(): Inf encountered.")
            return 0.
        
        model = self.cpte_model(param)
        like = self.logL(self.data,self.errdata,model)
        
        return like

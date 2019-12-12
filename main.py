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
#main body of code


class StarRotator(object):
    def __init__(self,wave_start,wave_end,grid_size,star_path='demo_star.txt',planet_path='demo_planet.txt',obs_path='demo_observations.txt'):
        """
            Welcome to StarRotator.
            ***********************
            Star Rotator needs the following input:
            wave_start: Wavelength range start in nm in vacuum. type:float
            wave_end: Wavelength range end in nm in vacuum. type:float
            grid_size: Number of grid cells, default 500. type:int
            star_path, planet_path and obs_path point to textfiles that contain
            the parameters of the star, the system and the time-series of the observations
            (either in phase or in time.)
            Good luck!

            After initializing the class like a=StarRotator(588.0,592.0,200.0),
            the simulated spectra can be accessed as:
            a.wl (wavelength array)
            a.spectra (matrix of spectra)
            a.lightcurve (list of fluxes)
        """
        self.wave_start=wave_start
        self.wave_end=wave_end
        self.grid_size=grid_size
        self.read_system(star_path=star_path,planet_path=planet_path,obs_path=obs_path)
        self.compute_spectrum()

        # return(self.wlF,F_out)#This is mad return statement. This whole function should be a class instead.

    def read_system(self,star_path='demo_star.txt',planet_path='demo_planet.txt',obs_path='demo_observations.txt'):
        """Reads in the stellar, planet and observation parameters from file; performing
        tests on the input and lifting the read variables to the class object.

        Input:
            star_path: str
            planet_path: str
            obs_path: str
        Returns:
            None
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
        self.mus = float(starparams[8].split()[0])
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
            self.mus = np.linspace(0.0,1.0,mus)#Uncomment this to run in CLV mode with SPECTRUM.


    def compute_spectrum(self):
        """This is where the main computation takes place. Output is passed to self."""
        #Two arrays for the x and y axes
        self.x = np.linspace(-1,1,num=2*self.grid_size) #in units of stellar radius
        self.y = np.linspace(-1,1,num=2*self.grid_size) #in units of stellar radius
        #Calculate the velocity and flux grids
        print('--- Computing velocity/flux grids')
        self.vel_grid = vgrid.calc_vel_stellar(self.x,self.y, self.stelinc, self.velStar, self.drr,  self.pob)
        self.flux_grid = vgrid.calc_flux_stellar(self.x,self.y,self.u1,self.u2)

        if isinstance(self.mus,np.ndarray) != True:#SWITCH BETWEEN PHOENIX (mus=0) AND SPECTRUM
            print('--- Reading spectrum from PHOENIX')
            wl,fx = spectrum.read_spectrum(self.T,self.logg)
            print('--- Integrating disk')
            if  self.drr == 0:
                print('------ Fast integration')
                wlF,F = integrate.build_spectrum_fast(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            else:
                print('------ Slow integration')
                wlF,F = integrate.build_spectrum_slow(wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
        else:#Meaning, if we have no mu's do:
            test.test_KURUCZ()
            print('--- Computing limb-resolved spectra with SPECTRUM')
            wl,fx_list = spectrum.compute_spectrum(self.T,self.logg,self.Z,self.mus,self.wave_start,self.wave_end,mode='anM')
            print('--- Integrating limb-resolved disk')
            wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,self.mus, self.wave_start,self.wave_end,self.x,self.y,self.vel_grid)


        self.xp,self.yp,self.zp = ppos.calc_planet_pos(self.sma_Rs, self.ecc, self.omega, self.orbinc, self.Rp_Rs, self.orb_p, self.transitC, self.mode, self.times, self.exptimes)

        F_out = np.zeros((self.Nexp,len(F)))
        flux_out = []
        mask_out = []
        for i in range(self.Nexp):
            wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(self.xp[i],self.yp[i],self.Rp_Rs,wl,fx,self.wave_start,self.wave_end,self.x,self.y,self.vel_grid,self.flux_grid)
            F_out[i,:]=F-Fp
            flux_out.append(flux)
            mask_out.append(mask)
        #This defines the output.
        self.wl = wlF
        self.stellar_spectrum = F
        self.spectra = F_out
        self.lightcurve = flux_out
        self.masks = mask_out

        def residuals(self):
            """Compute and return the residuals of the time_series of spectra.
            This may be the primary output of StarRotator."""
            self.residual = self.spectra*0.0
            for i in range(self.Nexp):
                self.residual[i,:]=self.spectra[i]/self.stellar_spectrum
            return(self.residual)

        def plot_residuals(self):
            """Compute and plot the residuals."""
            import matplotlib.pyplot as plt
            res = self.residuals()
            plt.pcolormesh(self.wl,self.times,res)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Phase')
            plt.show()

    def animate(self):
        import matplotlib.pyplot as plt
        import lib.integrate as integrate
        import numpy as np
        from matplotlib.patches import Circle
        import shutil
        import os
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
            planet1 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black', edgecolor='black', lw=1)
            planet2 = Circle((self.xp[i],self.yp[i]),self.Rp_Rs, facecolor='black', edgecolor='black', lw=1)
            ax[0][0].add_patch(planet1)
            ax[1][0].add_patch(planet2)
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
            sf = 20.0
            ax2.set_ylim((1.0-(1-yl[0])/sf,1.0+(yl[1]-1)/sf))
            ax2.set_ylabel('Ratio F_in/F_out',fontsize = 7)
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
            ax[0][1].set_xlabel('Timestep',fontsize='small')
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
        print('')
        os.system('convert -delay 6 anim/*.png animation.gif')

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




    # void1,void2,minflux,void3 = integrate.build_local_spectrum_fast(0,0,RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)

    #The following puts a large circular spot with a T 1000K less than the star in the center.
    # wl2,fx2 = spectrum.read_spectrum(T-1000,logg)
    # wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(0,0,0.2,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
    # wls,Fs,fluxs,masks = integrate.build_local_spectrum_fast(0,0,0.2,wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)
    # Ft = F-Fp+Fs
    # plt.plot(wlF,F)
    # plt.plot(wlF,Ft)
    # plt.show()
    # sys.exit()

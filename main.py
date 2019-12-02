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
#main body of code



def StarRotator(wave_start,wave_end,velStar,stelinc,orbinc,drr,pob,T,Z,logg,grid_size=500,flag=''):

    if flag == 'help' or flag == 'Help':
        print("Welcome to StarRotator.")
        print("***********************")
        print("Star Rotator needs the following input:")
        print("wave_start: Wavelength range start in nm in vacuum. type:float")
        print("wave_end: Wavelength range end in nm in vacuum. type:float")
        print("velStar: Equatiorial velocity of the star in m per s. type:float")
        print("stelinc: Inclination of the star. type:float")
        print("orbinc: Orbital inclination of the planet. type:float")
        print("drr: Differential rotation rate in units of stellar radii. type:float")
        print("pob: Projected obliquity of the star in degrees. type:float")
        print("T: Temperature of the star in K. type:float")
        print("Z: Metallicity of the star (dex). type:float")
        print("logg: Logarithmic gravity of the star in cgs. type:float") #4.5
        print("grid_size: Number of grid cells, default 500. type:int")
        print("Good luck!")
        sys.exit()

    try:
        test.typetest(wave_start,float,varname='wave_start in input')
        test.notnegativetest(wave_start,float,varname='wave_start in input')
        test.nantest(wave_start,float,varname='wave_start in input')
        #add all the other input parameters
    except ValueError as err:
        print("Parser: ",err.args)

    # test.test_integrators()
    u1 = 0.93
    u2 = -0.23#THESE ARE THE WINN2011 PAREMETERS, NOT KIPPING 2013.
    mus = 0#Normal operation without CLV.
    mus = np.linspace(0.0,1.0,3)#Uncomment this to run in CLV mode with SPECTRUM.

    #Two arrays for the x and y axes
    x = np.linspace(-1,1,num=2* grid_size) #in units of stellar radius
    y = np.linspace(-1,1,num=2* grid_size) #in units of stellar radius

    #Calculate the velocity and flux grids
    print('--- Computing velocity/flux grids')
    vel_grid = vgrid.calc_vel_stellar(x,y, stelinc, velStar, drr,  pob)
    flux_grid = vgrid.calc_flux_stellar(x,y,u1,u2)

    if isinstance(mus,np.ndarray) != True:#SWITCH BETWEEN PHOENIX (mus=0) AND SPECTRUM
        print('--- Reading spectrum from PHOENIX')
        wl,fx = spectrum.read_spectrum(T,logg)
        print('--- Integrating disk')
        if  drr == 0:
            print('------ Fast integration')
            wlF,F = integrate.build_spectrum_fast(wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
        else:
            print('------ Slow integration')
            wlF,F = integrate.build_spectrum_slow(wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
    else:#Meaning, if we have no mu's do:
        test.test_KURUCZ()
        print('--- Computing limb-resolved spectra with SPECTRUM')
        wl,fx_list = spectrum.compute_spectrum(T,logg,Z,mus, wave_start, wave_end,mode='anM')
        print('--- Integrating limb-resolved disk')
        wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,mus, wave_start, wave_end,x,y,vel_grid)
        wl2,fx2 = spectrum.read_spectrum(T,logg)
        wlF2,F2 = integrate.build_spectrum_fast(wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)

        wlp,Fp,flux,mask = integrate.build_local_spectrum_limb_resolved(-0.3,0.0,0.1,wl,fx_list,mus, wave_start, wave_end,x,y,vel_grid)
        wlp2,Fp2,flux2,mask2 = integrate.build_local_spectrum_fast(-0.3,0.0,0.1,wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)
        #This overplots non-rotating SPECTRUM and PHOENIX spectra, normalised.
        # plt.plot(wl,fx_list[-1]/max(fx_list[-1]),label='SPECTRUM')
        # plt.plot(wl2,fx2/5e15,label='PHOENIX')
        # plt.xlabel('Wavelength (nm)')
        # plt.ylabel('Max-normalised flux')
        # plt.title('T = %s K, log(g) = %s' % (T,logg))
        # plt.legend()
        # plt.show()

        plt.plot(wlF,F/max(F),color='skyblue',alpha=0.5)
        plt.plot(wlF,(F-Fp)/np.nanmax(F-Fp),color='skyblue',label='SPECTRUM')
        plt.plot(wlF2,F2/max(F2),color='red',alpha=0.5)
        plt.plot(wlF2,(F2-Fp2)/np.nanmax(F2-Fp2),color='red',label='PHOENIX')
        plt.legend()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Max-normalised flux')
        plt.title('T = %s K, log(g) = %s, vsini = 110km/s' % (T,logg))
        plt.show()

    sys.exit()

    #The following puts a large circular spot with a T 1000K less than the star in the center.
    # wl2,fx2 = spectrum.read_spectrum(T-1000,logg)
    # wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(0,0,0.2,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
    # wls,Fs,fluxs,masks = integrate.build_local_spectrum_fast(0,0,0.2,wl2,fx2, wave_start, wave_end,x,y,vel_grid,flux_grid)
    # Ft = F-Fp+Fs
    # plt.plot(wlF,F)
    # plt.plot(wlF,Ft)
    # plt.show()
    # sys.exit()


    #The following creates a transiting planet. We will need to offload it to some
    #tutorial-kind of script, but I put it here for now.

    RpRs = 0.12247
    nsteps = 100
    xp = np.linspace(-1.5,1.5,nsteps)
    yp = np.linspace(-0.55,-0.45,nsteps)

    lightcurve = []#flux points will be appended onto this.
    void1,void2,minflux,void3 = integrate.build_local_spectrum_fast(0,0,RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
    for i in range(len(xp)):
        wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(xp[i],yp[i],RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
        lightcurve.append(flux)
        fig,ax = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
        ax[0][0].pcolormesh(x,y,flux_grid*mask,vmin=0,vmax=1.0*np.nanmax(flux_grid),cmap='autumn')
        ax[1][0].pcolormesh(x,y,vel_grid*mask,cmap='bwr')
        ax[0][0].axes.set_aspect('equal')
        ax[1][0].axes.set_aspect('equal')
        planet1 = Circle((xp[i],yp[i]),RpRs, facecolor='black', edgecolor='black', lw=1)
        planet2 = Circle((xp[i],yp[i]),RpRs, facecolor='black', edgecolor='black', lw=1)
        ax[0][0].add_patch(planet1)
        ax[1][0].add_patch(planet2)


        ax[0][1].plot(lightcurve,'.',color='black')
        ax[0][1].set_xlim((0,nsteps))
        ax[0][1].set_ylim((minflux-0.1*RpRs**2.0),1.0+0.1*RpRs**2)
        ax[1][1].plot(wlF,F/np.nanmax(F),color='black',alpha = 0.5)
        ymin = np.nanmin(F/np.nanmax(F))
        ymax = np.nanmax(F/np.nanmax(F))
        linedepth = ymax - ymin
        ax[1][1].plot(wlF,(F-Fp)/np.nanmax(F-Fp),color='black')
        ax[1][1].set_xlim((588.5,590.2))
        yl = (ymin-0.1*linedepth,ymax+0.3*linedepth)
        ax[1][1].set_ylim(yl)
        ax2 = ax[1][1].twinx()
        ax2.plot(wlF,(F-Fp)*np.nanmax(F)/F/np.nanmax(F-Fp),color='skyblue')
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
        # plt.show()
        # sys.exit()
        if len(str(i)) == 1:
            out = '000'+str(i)
        if len(str(i)) == 2:
            out = '00'+str(i)
        if len(str(i)) == 3:
            out = '0'+str(i)
        if len(str(i)) == 4:
            out = str(i)

        fig.savefig('anim/'+out+'.png', dpi=fig.dpi)
        integrate.statusbar(i,xp)
        plt.close()
        #convert -delay 6 anim/*.png HD209458b.gif


if __name__ == '__main__':
   StarRotator('help')

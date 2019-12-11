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


def StarRotator(wave_start,wave_end,grid_size,star_path='demo_star.txt',planet_path='demo_planet.txt',obs_path='demo_observations.txt'):
#def StarRotator(wave_start,wave_end,velStar,stelinc,orbinc,drr,pob,T,Z,logg,grid_size=500,u1 = 0,u2=0):
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
    """

    planetparams = open(planet_path,'r').read().splitlines()
    starparams = open(star_path,'r').read().splitlines()
    obsparams = open(obs_path,'r').read().splitlines()

    velStar = float(starparams[0].split()[0])
    stelinc = float(starparams[1].split()[0])
    drr = float(starparams[2].split()[0])
    T = float(starparams[3].split()[0])
    Z = float(starparams[4].split()[0])
    logg = float(starparams[5].split()[0])
    u1 = float(starparams[6].split()[0])
    u2 = float(starparams[7].split()[0])
    mus = float(starparams[8].split()[0])

    sma_Rs = float(planetparams[0].split()[0])
    ecc = float(planetparams[1].split()[0])
    omega = float(planetparams[2].split()[0])
    orbinc = float(planetparams[3].split()[0])
    pob = float(planetparams[4].split()[0])
    Rp_Rs = float(planetparams[5].split()[0])
    orb_p = float(planetparams[6].split()[0])
    transitC = float(planetparams[7].split()[0])
    mode = planetparams[8].split()[0]#This is a string.

    times = [] #These are in JD-24000000.0 or in orbital phase.
    exptimes = []
    for i in obsparams:
        times.append(float(i.split()[0]))
    times = np.array(times)
    if mode == 'times':
        for i in obsparams:
            exptimes.append(float(i.split()[1]))
        exptimes = np.array(exptimes)



    try:
        test.typetest(wave_start,float,varname='wave_start in input')
        test.nantest(wave_start,varname='wave_start in input')
        test.notnegativetest(wave_start,varname='wave_start in input')
        test.notnegativetest(velStar,varname='velStar in input')
        test.notnegativetest(stelinc,varname='stelinc in input')

        #add all the other input parameters
    except ValueError as err:
        print("Parser: ",err.args)

    # test.test_integrators()
    # u1 = 0.93
    # u2 = -0.23#THESE ARE THE WINN2011 PAREMETERS, NOT KIPPING 2013.
    # mus = 0#Normal operation without CLV.

    if mus != 0:
        mus = np.linspace(0.0,1.0,mus)#Uncomment this to run in CLV mode with SPECTRUM.

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
            wlF,F = integrate.build_spectrum_fast(wl,fx,wave_start,wave_end,x,y,vel_grid,flux_grid)
        else:
            print('------ Slow integration')
            wlF,F = integrate.build_spectrum_slow(wl,fx,wave_start,wave_end,x,y,vel_grid,flux_grid)
    else:#Meaning, if we have no mu's do:
        test.test_KURUCZ()
        print('--- Computing limb-resolved spectra with SPECTRUM')
        wl,fx_list = spectrum.compute_spectrum(T,logg,Z,mus, wave_start, wave_end,mode='anM')
        print('--- Integrating limb-resolved disk')
        wlF,F = integrate.build_spectrum_limb_resolved(wl,fx_list,mus, wave_start, wave_end,x,y,vel_grid)




    xp,yp,zp = ppos.calc_planet_pos(sma_Rs, ecc, omega, orbinc, Rp_Rs, orb_p, transitC, mode, times, exptimes)
    plt.plot(xp,yp,'.')
    plt.xlim(-0.3,0.3)
    plt.ylim(-0.3,0.3)
    plt.show()
    sys.exit()
    for i in range(len(xp)):
        wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(xp[i],yp[i],RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)





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

    sys.exit()



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

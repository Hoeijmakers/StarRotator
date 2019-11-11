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
import lib.plotting as plt
import lib.operations as ops
import lib.stellar_spectrum as spectrum
import lib.integrate as integrate
import pdb
import time
import matplotlib.pyplot as pl
import math
#main body of code




if __name__ == '__main__':

    #call parser, gets the input parameters from the user
    try:
        parser = argparse.ArgumentParser(description='Input variables for StarRotator:')
        parser.add_argument('wave_start', metavar='waveS', type=float,
                            help='Wavelength range start in nm in vacuum. type:float')
        parser.add_argument('wave_end', metavar='waveE', type=float,
        help='Wavelength range end in nm in vacuum. type:float')
        parser.add_argument('velStar', metavar='velS', type=float,
        help='Equatiorial velocity of the star in m per s. type:float')
        parser.add_argument('stelinc', metavar='stelinc', type=float,
        help='Inclination of the star. type:float')
        parser.add_argument('orbinc', metavar='orbinc', type=float,
        help='Orbital inclination of the planet. type:float')
        parser.add_argument('drr', metavar='drr', type=float,
        help='Differential rotation rate in units of stellar radii. type:float')
        parser.add_argument('pob', metavar='pob', type=float,
        help='projected obliquity of the star in degrees. type:float')
        parser.add_argument('grid_size', metavar='grid', type=int, default=500,
        help='number of grid cells, default 500. type:int')
        args = parser.parse_args()
        test.test_parser(args)
    except ValueError as err:
        print("Parser: ",err.args)

    test.test_integrators()

    T = 10000.0
    logg = 4.5
    u1 = 0.387
<<<<<<< HEAD
    # u1 = 0.008
    u2 = 0.178
    # u2 = 0.1
=======
    #u1 = 0.0
    u2 = 0.178
#    u2 = 0.1
>>>>>>> c838650572c33aa34b147322916c94993c502f9e

    #two arrays for the x and y axis
    x = np.linspace(-1,1,num=2*args.grid_size) #in units of stellar radius
    y = np.linspace(-1,1,num=2*args.grid_size) #in units of stellar radius

    #calculate the velocity grid
    vel_grid = vgrid.calc_vel_stellar(x,y,args.stelinc,args.velStar,args.drr, args.pob)
    flux_grid = vgrid.calc_flux_stellar(x,y,u1,u2)
    wl,fx = spectrum.read_spectrum(T,logg)

    if args.drr == 0:
        wlF,F = integrate.build_spectrum_fast(wl,fx,args.wave_start,args.wave_end,x,y,vel_grid,flux_grid)
    else:
        wlF,F = integrate.build_spectrum_slow(wl,fx,args.wave_start,args.wave_end,x,y,vel_grid,flux_grid)


    #The following creates a transiting planet. We will need to offload it to some
    #tutorial-kind of script, but I put it here for now.
    from matplotlib.patches import Circle


    RpRs = 0.13
    nsteps = 100
    xp1 = np.linspace(-1.5,1.5,nsteps)
    yp1 = np.linspace(-0.5,0.2,nsteps)
    # xp2 = np.linspace(-0.1,-0.5,int(nsteps*1.0))
    # yp2 = np.linspace(2.5,-2.5,int(nsteps*1.0))
    # xp3 = np.linspace(1.3,0.0,nsteps)
    # yp3 = np.linspace(-0.5,1.3,nsteps)

    # xp = np.concatenate((xp1,xp2,xp3))
    # yp = np.concatenate((yp1,yp2,yp3))
    xp=xp1
    yp=yp1
    lightcurve = []
    void1,void2,minflux,void3 = integrate.build_local_spectrum_fast(0,0,0.16,wl,fx,args.wave_start,args.wave_end,x,y,vel_grid,flux_grid)
    for i in range(len(xp)):
        if i == nsteps-1:
            RpRs = 0.1
            lightcurve = []
        if i == 2*nsteps-1:
            lightcurve = []
            RpRs = 0.16
        # i=150
        wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(xp[i],yp[i],RpRs,wl,fx,args.wave_start,args.wave_end,x,y,vel_grid,flux_grid)
        lightcurve.append(flux)
        fig,ax = pl.subplots(nrows=2, ncols=2,figsize=(8,8))
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
        ax[1][1].set_ylim((ymin-0.1*linedepth,ymax+0.3*linedepth))
        ax2 = ax[1][1].twinx()
        ax2.plot(wlF,(F-Fp)*np.nanmax(F)/F/np.nanmax(F-Fp),color='skyblue')
        ax2.set_ylim((1.0-(1-ymin-0.1*linedepth)/2,1.0+(ymax+0.3*linedepth-1)/2))
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
        # pl.show()
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
        pl.close()

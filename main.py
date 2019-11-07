######################
#authors: Jens Hoeijmakers and Julia Seidel
#Description: Calculates the stellar spectrum
#
#
#####################

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

    T = 5000.0
    logg = 4.5
    u1 = 0.387
    u2 = 0.178

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

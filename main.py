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
import pdb
import time
import matplotlib.pyplot as pl
#main body of code


def statusbar(i,x):
    if type(x) == int:
        print('  '+f"{i/(float(x)-1)*100:.1f} %", end="\r")
    else:
        print('  '+f"{i/(len(x)-1)*100:.1f} %", end="\r")#Statusbar.

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


    T = 5000.0
    logg = 4.5


    #two arrays for the x and y axis
    x = np.linspace(-1,1,num=2*args.grid_size) #in units of stellar radius
    y = np.linspace(-1,1,num=2*args.grid_size) #in units of stellar radius

    #calculate the velocity grid
    vel_grid = vgrid.calc_vel_stellar(x,y,args.stelinc,args.velStar,args.drr, args.pob)

    plt.plot_star_2D(x,y,vel_grid,cmap="hot",quantities=['x-coordinate','y-coordinate','Radial velocity'],units=['R_*','R_*','(m/s)'],noshow=False)

    F = 0#output
    wl,fx = spectrum.read_spectrum(T,logg)

    wlc = wl[(wl >= args.wave_start) & (wl <= args.wave_end)]
    fxc = fx[(wl >= args.wave_start) & (wl <= args.wave_end)]

    #To do:
    #De-rotate the axes.
    #Pre-calculate mu grid along with v grid.
    #Deal with edge effects interpolation vis a vis wl clipping.
    #Debug shifting.
    start = time.time()
    for i in range(len(x)):
        for j in range(len(y)):
            F+=ops.shift(wlc,fxc,vel_grid[j,i])*ops.limb_darkening(np.sqrt(x[i]**2+y[j]**2),0.387,0.178)
        statusbar(i,len(x))
    print(time.time()-start)
    pl.plot(wlc,fxc)
    pl.plot(wlc,F)
    pl.show()

    F_behind_planet = =
    for i_in_planet:
        for j_in_planet:
    output = F -

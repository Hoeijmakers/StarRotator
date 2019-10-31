######################
#
#tbd
#
#
#####################

#import statements
import sys
import numpy as np
import argparse
import lib.test as test

#helper function to convert angles in radians
#comment Julia: maybe we can create a utils.py file?
def convert_angle_to_rad(angle):
    return 2.*np.pi/360.*angle
    
#import function files
def calc_vel_stellar(x,y,i_stellar, vel_eq, diff_rot_rate, proj_obliquity):
    """
    based on Cegla+2016. See Fig. 3 and equations 2 - 8
    https://arxiv.org/pdf/1602.00322.pdf
    It takes the stellar parameters and then calculates the stellar velocities in all bins for one quarter of the stellar disk
    input: x, y: vectors to create the stellar grid in units of stellar radius
    """
    #careful! all angles have to be in radiant!
    #Think carefully about which ones of these have to be transposed
    #x_sq = np.einsum('i,j->ij',x,x)
    #y_sq = np.einsum('i,j->ij',y,y)
    alpha = convert_angle_to_rad(proj_obliquity)
    beta = convert_angle_to_rad(i_stellar)
    xy = np.einsum('i,j->ij',x,y)
    x_full = np.tile(x,(len(x),1))
    y_full = np.tile(y,(len(y),1)).T
    x_sq = x_full*x_full
    y_sq=y_full*y_full
    z = np.sqrt(1.-x_sq-y_sq+2.*x_full*y_full*np.sin(2.*alpha))
    #I have to add the stellar radius back in here somewhere
    vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*vel_eq*np.sin(beta*(1-diff_rot_rate*z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*(x_full*np.sin(alpha)-y_full*np.cos(alpha))))
    
    return vel_stellar_grid
#maybe some global variables here


if __name__ == '__main__':

#main body of code

#call parser
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
        parser.add_argument('grid_size', metavar='grid', type=int, default=500,
        help='number of grid cells, default 500. type:int')
        args = parser.parse_args()
        test.test_parser(args)
    except ValueError as err:
        print("Parser: ",err.args)

    #comment Julia: drr differential rotation rate, where do we get this from
    #comment Julia: pob projected obliquity, where do we get this from
    #for testing
    drr = 1.
    pob = 80.
    x = np.linspace(0,1,num=args.grid_size) #in units of stellar radius
    y = np.linspace(0,1,num=args.grid_size) #in units of stellar radius

    vel_grid = calc_vel_stellar(x,y,args.stelinc,args.velStar,drr, pob)
    print(vel_grid)

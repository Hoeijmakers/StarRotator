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
    input: x, y: 1D numpy arrays to create the stellar grid in units of stellar radius
           i_stellar: inclination in degrees
           vel_eq: equatorial stellar velocity
           diff_rot_rate: differential rotation rate
           proj_obliquity: projected obliquity
    output: vel_stellar_grid: 2D numpy array of stellar velocities over one quarter of the stellar disk
    """
    #careful! all angles have to be in radiant!
    #Think carefully about which ones of these have to be transposed

    #convert angles to rad
    alpha = convert_angle_to_rad(proj_obliquity)
    beta = convert_angle_to_rad(i_stellar)
    #pre calculate matrices
    
    xy = np.einsum('i,j->ij',x,y)
    x_full = np.tile(x,(len(x),1))
    y_full = np.tile(y,(len(y),1)).T
    x_sq = x_full*x_full
    y_sq=y_full*y_full
    #this is the z coordinate in the tilted coordinate system from Cegla+2016
    # with some trigonometry magic
    z = np.sqrt(1.-x_sq-y_sq+4.*xy*np.sin(alpha)*np.cos(alpha))

    #equation 8 from Cegla+2016
    vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*(x_full*np.sin(alpha)-y_full*np.cos(alpha))))
    
    return vel_stellar_grid

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
    
    #two arrays for the x and y axis
    x = np.linspace(0,1,num=args.grid_size) #in units of stellar radius
    y = np.linspace(0,1,num=args.grid_size) #in units of stellar radius

    #calculate the velocity grid
    vel_grid = calc_vel_stellar(x,y,args.stelinc,args.velStar,args.drr, args.pob)
    print(vel_grid)

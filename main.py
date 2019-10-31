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

#import function files

#maybe some global variables here


if __name__ == '__main__':

#main body of code

#call parser
    parser = argparse.ArgumentParser(description='Input variables for StarRotator:')
    parser.add_argument('wave_start', metavar='waveS', type=float, nargs='+',
                        help='Wavelength range start in nm in vacuum. type:float')
    parser.add_argument('wave_end', metavar='waveE', type=float, nargs='+',
    help='Wavelength range end in nm in vacuum. type:float')
    args = parser.parse_args()
    print(args)

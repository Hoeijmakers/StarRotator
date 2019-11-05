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
    import numpy as np
    #careful! all angles have to be in radiant!
    #Think carefully about which ones of these have to be transposed

    #convert angles to rad
    alpha = np.radians(proj_obliquity)
    beta = np.radians(i_stellar)
    #pre calculate matrices

    xy = np.einsum('i,j->ij',x,y)
    x_full = np.tile(x,(len(x),1))
    y_full = np.tile(y,(len(y),1)).T
    x_sq = x_full*x_full
    y_sq=y_full*y_full
    #this is the z coordinate in the tilted coordinate system from Cegla+2016
    # with some trigonometry magic
    z = np.sqrt(1.-x_sq-y_sq)

    #equation 8 from Cegla+2016
    vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*(x_full*np.sin(alpha)-y_full*np.cos(alpha))))

    return vel_stellar_grid

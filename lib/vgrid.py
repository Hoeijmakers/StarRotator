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
    import pdb
    import lib.test as test

    #Standard tests on input
    test.typetest(x,np.ndarray,varname='x in calc_vel_stellar')
    test.typetest(y,np.ndarray,varname='x in calc_vel_stellar')
    test.typetest(i_stellar,float,varname='i_stellar in calc_vel_stellar')
    test.typetest(vel_eq,float,varname='vel_eq in calc_vel_stellar')
    test.typetest(diff_rot_rate,float,varname='diff_rot_rate in calc_vel_stellar')
    test.typetest(proj_obliquity,float,varname='proj_obliquity in calc_vel_stellar')
    test.nantest(x,varname='x in calc_vel_stellar')
    test.nantest(y,varname='y in calc_vel_stellar')
    test.nantest(i_stellar,varname='i_stellar in calc_vel_stellar')
    test.nantest(vel_eq,varname='vel_eq in calc_vel_stellar')
    test.nantest(diff_rot_rate,varname='diff_rot_rate in calc_vel_stellar')
    test.nantest(proj_obliquity,varname='proj_obliquity in calc_vel_stellar')
    test.notnegativetest(vel_eq,varname='vel_eq in calc_vel_stellar')

    #Careful! all angles have to be in radians!
    #Think carefully about which ones of these have to be transposed

    #convert angles to rad
    alpha = np.radians(proj_obliquity)
    beta = np.radians(i_stellar)
    #pre calculate matrices
    z,x_full,y_full = calc_z(x,y)
    #equation 8 from Cegla+2016
#    vel_stellar_grid = (x_full*np.cos(alpha)-y_full*np.sin(alpha))*vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*(x_full*np.sin(alpha)-y_full*np.cos(alpha))))

     #this has the star aligned to the coordinate system. Like this the projected obliquity is not incorporated here. The tilt of the star compared to the planet has to be incorporated in the planet's path.
    vel_stellar_grid = x_full*vel_eq*np.sin(beta)*(1.-diff_rot_rate*(z*np.sin(np.pi/2.-beta)+np.cos(np.pi/2.-beta)*y_full))

    return(vel_stellar_grid)

def calc_z(x,y):
    """This is a nice little wrapper to calculate the distance matrix of a circle."""
    import numpy as np
    xy = np.einsum('i,j->ij',x,y)
    x_full = np.tile(x,(len(y),1))
    y_full = np.tile(y,(len(x),1)).T
    x_sq = x_full*x_full
    y_sq=y_full*y_full
    #this is the z coordinate in the tilted coordinate system from Cegla+2016
    # with some trigonometry magic
    #it's also the z coordinate in the untilted system for similarity reasons
    d = 1.-x_sq-y_sq
    d[d<0]=np.nan#This supresses a runtime warning.
    return(np.sqrt(d),x_full,y_full)


def calc_flux_stellar(x,y,u1,u2):
    import numpy as np
    import lib.operations as ops
    import lib.test as test

    test.typetest(x,np.ndarray,varname='x in calc_flux_stellar')
    test.typetest(y,np.ndarray,varname='x in calc_flux_stellar')
    test.typetest(u1,[int,float],varname='i_stellar in calc_flux_stellar')
    test.typetest(u2,[int,float],varname='vel_eq in calc_flux_stellar')
    test.nantest(x,varname='x in calc_flux_stellar')
    test.nantest(y,varname='y in calc_flux_stellar')
    test.nantest(u1,varname='x in calc_flux_stellar')
    test.nantest(u2,varname='y in calc_flux_stellar')

    z,x_full,y_full = calc_z(x,y)
    flux_grid = z*0.0
    for i in range(len(x)):
        for j in range(len(y)):
            if np.sqrt(x[i]**2+y[j]**2) <= 1.0:
                flux_grid[j,i]=ops.limb_darkening((np.sqrt(x[i]**2+y[j]**2)),u1,u2)*(z[j,i]*0.0+1.0)#Last multiplication is to put NaNs back into place.

    # plotting.plot_star_2D(x,y,mu_grid,cmap="hot",quantities=['','',''],units=['','',''],noshow=False)
    flux_grid /= np.nansum(flux_grid)
    return(flux_grid)

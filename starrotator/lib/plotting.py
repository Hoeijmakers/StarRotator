
import numpy as np

def spot_masks(xp,yp,Rp,N=1000):
    """Plots spots on a stellar disk, where the spots have
    x-y locations in the projected image plane, and radii in units
    of stellar radius Rp. N is the number of pixels.
    
    Returns the mask of the stellar disk (m1) and the sum of all spots (m2, becomes greater than 1 for overlapping spots).
    
    Use as:
    m1,m2 = spot_masks(xp,yp,Rp)

    plt.imshow(m1+m2).
    """
    x = np.linspace(-1,1,N)
    dx = x[1]-x[0]
    zp = np.sqrt(1-xp**2-yp**2)

    X, Y = np.meshgrid(x, x)
    r2 = X**2 + Y**2
    Z = np.sqrt(np.clip(1 - r2, 0, 1))

    n_s = np.array([xp, yp, zp]).T #Unit vector in the direction normal to the spot.
    # stack grid normals


    grid = np.stack([X, Y, Z], axis=-1)   # (S, Ny, Nx, 3)
    grid = grid[None,:,:,:]
    # align spot normals
    n_s = n_s[:, None, None, :]              # (S, 1, 1, 3)


    # dot product
    dot_product = np.sum(grid * n_s, axis=-1)      # (S, Ny, Nx)


    # Note that with R equal to the spot radius, if we assume that to be the great-circle distance, then R = alpha.
    # Otherwise, alpha = jnp.asin(Rp). I use the latter because I think its more observationally meaningful.
    #Since we want to know cos(alpha), and alpha=arcsin(R), we get cos(alpha) = sqrt(1-R**2)
    cos_alpha = np.sqrt(1-Rp**2)
    # optional: broadcast cos_alpha
    cos_alpha = cos_alpha[:, None, None]


    mask1 = r2 <= 1.0
    mask2 = (r2 <= 1.0) & (dot_product >= cos_alpha)


    # img = r2*0.0+mask1+jnp.clip(np.sum(mask2,axis=0),0,1)

    return(mask1,np.sum(mask2,axis=0))


def latlon_to_xyz(lat, lon, radius=1.0):
    """
    Project latitude/longitude on a stellar surface to Cartesian x, y, z.

    Coordinate convention:
    - Stellar spin axis is the +y direction.
    - z is the line of sight.
    - x is the remaining sky-plane direction.
    - Inclination is 90 degrees, so the spin axis lies in the sky plane.
    - lat is stellar latitude, in radians.
    - lon is longitude, in radians, with lon=0 at the sub-observer meridian.

    Parameters
    ----------
    lat : float or array-like
        Latitude in degrees, from -90 to 90.
    lon : float or array-like
        Longitude in degrees.
    radius : float
        Stellar radius.

    Returns
    -------
    x, y, z : float or ndarray
        Cartesian coordinates with z in the line of sight.
    """
    lat = np.asarray(np.radians(lat))
    lon = np.asarray(np.radians(lon))

    x = radius * np.cos(lat) * np.sin(lon)
    y = radius * np.sin(lat)
    z = radius * np.cos(lat) * np.cos(lon)

    return x, y, z




def stellar_surface_xyz(lat, lon, radius=1.0, inclination=90.0):
    """
    Project latitude/longitude on a stellar surface to Cartesian x, y, z
    with an arbitrary stellar inclination.

    Coordinate convention:
    - Observer is along +z (line of sight).
    - x is horizontal in the sky plane.
    - y is vertical in the sky plane.
    - inclination = angle between spin axis and line of sight:
        i = 0   -> pole-on (spin axis along +z)
        i = pi/2 -> edge-on (spin axis along +y)

    Parameters
    ----------
    lat : float or array-like
        Latitude in radians (-pi/2 to +pi/2).
    lon : float or array-like
        Longitude in radians (lon=0 at sub-observer meridian).
    radius : float
        Stellar radius.
    inclination : float
        Inclination angle in radians.

    Returns
    -------
    x, y, z : float or ndarray
        Cartesian coordinates in observer frame.
    """
    lat = np.asarray(np.radians(lat))
    lon = np.asarray(np.radians(lon))

    # --- Step 1: coordinates in stellar frame (spin axis = +z) ---
    xs = radius * np.cos(lat) * np.sin(lon)
    ys = radius * np.sin(lat)
    zs = radius * np.cos(lat) * np.cos(lon)

    # --- Step 2: rotate around x-axis by inclination ---
    cosi = np.cos(np.radians(inclination-90))
    sini = np.sin(np.radians(inclination-90))

    x = xs
    y = ys * cosi - zs * sini
    z = ys * sini + zs * cosi

    return x, y, z










# The stuff below is old:

def plot_star_3D():
    #THIS IS NOT WORKING, ITS A WIP.

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.special import sph_harm     #import package to calculate spherical harmonics
    import pdb

    theta = np.linspace(0, 2*np.pi, 100)   #setting range for theta
    phi = np.linspace(0, np.pi, 100)       #setting range for phi
    phi, theta = np.meshgrid(phi, theta)   #setting the grid for phi and theta

    #Setting the cartesian coordinates of the unit sphere
    #Converting phi, theta, z to cartesian coordinates
    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(phi)

    #Setting the aspect ratio to 1 which makes the sphere look spherical and not elongated
    fig = plt.figure(figsize=plt.figaspect(1.))    #aspect ratio
    axes = fig.add_subplot(111, projection='3d')   #sets figure to 3d
    fig.suptitle('m=4   l=4', fontsize=18, x=0.52, y=.85)

    m, l = 4, 4   #m and l control the mode of pulsation and overall appearance of the figure
    #Calculating the spherical harmonic Y(l,m) and normalizing it
    axes.view_init(30, 45)
    plt.ion()
    plt.show()
    for idx,angle in enumerate(np.linspace(0,360,20)):
        figcolors = sph_harm(m, l, theta+angle, phi).real
        figmax, figmin = figcolors.max(), figcolors.min()
        figcolors = (figcolors-figmin)/(figmax-figmin)

        #Sets the plot surface and colors of the figure where seismic is the color scheme
        axes.plot_surface(x, y, z,  rstride=1, cstride=1,  facecolors=cm.autumn(figcolors))
        fig.canvas.draw_idle()
        pdb.set_trace()


def plot_star_2D(x,y,z,cmap="hot",quantities=['','',''],units=['','',''],noshow=False):
    """Plots the projected stellar disk.

        Parameters
        ----------
        x : np.array()
            The x coordinate of the map.
        y : np.array()
            The y coordinate of the map
        z : np.array()
            Two-dimensional image corresponding to the x and y axes.
        cmap: str (optional)
            The color map identifier corresponding to a matplotlib colormap.
            Defaults to "hot".
        quantities : list(str,str,str) (optional)
            A list of three strings corresponding to the axis labels (quantities).
        units : list(str,str,str) (optional)
            A list of three strings with the corresponding units.


        Returns
        -------
        An open matplotlib plot window.
    """
    import matplotlib.pyplot as plt
    if len(units) != 3:
        raise ValueError("For passing units, please provide a list containing three strings.")

    xlabel = quantities[0]
    ylabel = quantities[1]
    zlabel = quantities[2]

    if units[0] != '': xlabel+=' (%s)' % units[0]
    if units[1] != '': ylabel+=' (%s)' % units[1]
    if units[2] != '': zlabel+=' (%s)' % units[2]

    plt.imshow(z,cmap=cmap,extent=[min(x),max(x),min(y),max(y)])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cbar = plt.colorbar()
    cbar.set_label(zlabel, labelpad=0, rotation=270)
    if noshow == False:
        plt.show()
    return

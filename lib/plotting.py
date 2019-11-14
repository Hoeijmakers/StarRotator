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

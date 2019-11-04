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

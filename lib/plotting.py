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


#START HERE.
def plot_residuals(wlF,F,F_out,times):
    import matplotlib.pyplot as plt
    Nexp = len(F_out)
    residual = F_out*0.0
    for i in range(Nexp):
        residual[i,:]=F_out[i]/F
    plt.pcolormesh(wlF,times,residual)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Phase')
    plt.show()

# def animate_StarRotator(wave_start,wave_end,grid_size,star_path='demo_star.txt',planet_path='demo_planet.txt',obs_path='demo_observations.txt',out_path='anim.gif'):
def animate_StarRotator(wlF,F,F_out,flux_out,mask_out,times,Rp_Rs,xp,yp,zp,x,y,vel_grid,flux_grid):

    import matplotlib.pyplot as plt
    import lib.integrate as integrate
    import numpy as np
    from matplotlib.patches import Circle
    import shutil
    import os
    shutil.rmtree('anim/')#First delete the contents of the anim folder.
    os.mkdir('anim/')
    #The following creates a transiting planet. We will need to offload it to some
    #tutorial-kind of script, but I put it here for now.
    # RpRs = 0.12247
    # nsteps = 100
    # xp = np.linspace(-1.5,1.5,nsteps)
    # yp = np.linspace(-0.55,-0.45,nsteps)

    # lightcurve = []#flux points will be appended onto this.
    # void1,void2,minflux,void3 = integrate.build_local_spectrum_fast(0,0,RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
    lightcurve = np.array(flux_out)
    minflux = min(lightcurve)
    for i in range(len(xp)):
        # wlp,Fp,flux,mask = integrate.build_local_spectrum_fast(xp[i],yp[i],RpRs,wl,fx, wave_start, wave_end,x,y,vel_grid,flux_grid)
        # lightcurve.append(flux)
        mask = mask_out[i]

        fig,ax = plt.subplots(nrows=2, ncols=2,figsize=(8,8))
        ax[0][0].pcolormesh(x,y,flux_grid*mask,vmin=0,vmax=1.0*np.nanmax(flux_grid),cmap='autumn')
        ax[1][0].pcolormesh(x,y,vel_grid*mask,cmap='bwr')
        ax[0][0].axes.set_aspect('equal')
        ax[1][0].axes.set_aspect('equal')
        planet1 = Circle((xp[i],yp[i]),Rp_Rs, facecolor='black', edgecolor='black', lw=1)
        planet2 = Circle((xp[i],yp[i]),Rp_Rs, facecolor='black', edgecolor='black', lw=1)
        ax[0][0].add_patch(planet1)
        ax[1][0].add_patch(planet2)

        ax[0][1].plot(times[0:i],lightcurve[0:i],'.',color='black')
        ax[0][1].set_xlim((min(times),max(times)))
        ax[0][1].set_ylim((minflux-0.1*Rp_Rs**2.0),1.0+0.1*Rp_Rs**2)
        ax[1][1].plot(wlF,F/np.nanmax(F),color='black',alpha = 0.5)
        ymin = np.nanmin(F/np.nanmax(F))
        ymax = np.nanmax(F/np.nanmax(F))
        linedepth = ymax - ymin
        ax[1][1].plot(wlF,F_out[i]/np.nanmax(F_out[i]),color='black')
        # ax[1][1].set_xlim((588.5,590.2))
        yl = (ymin-0.1*linedepth,ymax+0.3*linedepth)
        ax[1][1].set_ylim(yl)
        ax2 = ax[1][1].twinx()
        ax2.plot(wlF,(F_out[i])*np.nanmax(F)/F/np.nanmax(F_out[i]),color='skyblue')
        sf = 20.0
        ax2.set_ylim((1.0-(1-yl[0])/sf,1.0+(yl[1]-1)/sf))
        ax2.set_ylabel('Ratio F_in/F_out',fontsize = 7)
        ax2.tick_params(axis='both', which='major', labelsize=6)
        ax2.tick_params(axis='both', which='minor', labelsize=5)

        ax[0][0].set_ylabel('Y (Rs)',fontsize=7)
        ax[0][0].tick_params(axis='both', which='major', labelsize=6)
        ax[0][0].tick_params(axis='both', which='minor', labelsize=5)

        ax[1][0].set_ylabel('Y (Rs)',fontsize=7)
        ax[1][0].set_xlabel('X (Rs)',fontsize=7)
        ax[1][0].tick_params(axis='both', which='major', labelsize=6)
        ax[1][0].tick_params(axis='both', which='minor', labelsize=5)

        ax[0][1].set_ylabel('Normalised flux',fontsize=7)
        ax[0][1].set_xlabel('Timestep',fontsize='small')
        ax[0][1].tick_params(axis='both', which='major', labelsize=6)
        ax[0][1].tick_params(axis='both', which='minor', labelsize=5)

        ax[1][1].set_ylabel('Normalised flux',fontsize=7)
        ax[1][1].set_xlabel('Wavelength (nm)',fontsize=7)
        ax[1][1].tick_params(axis='both', which='major', labelsize=6)
        ax[1][1].tick_params(axis='both', which='minor', labelsize=5)
        # plt.show()
        # sys.exit()
        if len(str(i)) == 1:
            out = '000'+str(i)
        if len(str(i)) == 2:
            out = '00'+str(i)
        if len(str(i)) == 3:
            out = '0'+str(i)
        if len(str(i)) == 4:
            out = str(i)

        fig.savefig('anim/'+out+'.png', dpi=fig.dpi)
        integrate.statusbar(i,xp)
        plt.close()
    print('')
    os.system('convert -delay 6 anim/*.png animation.gif')

def calc_orbit_times(time_stamp, transitC, exposure_time, orb_p):
    """
    Calculates the start, end and center of the transit from the provided times
    input:
        time_stamp: type: float, time
        transitC: type:float, transit center - 2400000.
        exposure_times: type: float, exposure time in seconds
        orb_p: type: float, orbital period in days
    output:
        orbit_start (, orbit_end, orbit_center): times of current exposure relative to transit start, end, center in days
    """
    #this allows to really think about start and end of the exposures. It might be interesting for long exposures, but I'm leaving this for later
    orbit_center = (time_stamp - transitC + 2400000.)% orb_p
    #orbit_start = (time_stamp-exposure_time/(2.*24.*3600.)- transitC + 2400000.) % orb_p
    #orbit_end = (time_stamp+exposure_time/(2.*24.*3600.)- transitC + 2400000.) % orb_p

    return  orbit_center# orbit_start, orbit_end,

def func_kepler(ecc_anom,mean_anom,ecc):
    """
    Find the eccentric anomaly using the mean anomaly and eccentricity
    via M = E - e sin(E)
    input:
        ecc_anom: type: list, eccentric anomaly
        mean_anom: type: list, mean anomaly
        ecc: type: float, eccentricity
    output:
        anome: type: list, eccentric anomaly
    """
    anom=ecc_anom-ecc*np.sin(ecc_anom)-mean_anom
    return anom

def calc_true_anom(ecc,phases,omega_bar):
    """
    Calculates the true anomaly and the mean anomaly of the system
    input:
        ecc: type: float, eccentricity
        phases: type: np.array, phases, either given by the user or calculated from time stamps
        omega_bar: type: float, angle
    output:

    """
    import numpy as np
    #Circular orbit
    if np.isclose(ecc,0.,1e-4):
        #True anomaly
        true_anom=2.*np.pi*phases
        #Eccentric anomaly
        ecc_anom=None
    #Eccentric orbit
    else:
        #True anomaly of the planet at mid-transit (in rad):
        #    - angle counted from 0 at the perisastron, to the star/Earth LOS
        #    - >0 counterclockwise, possibly modulo 2pi
        true_anom_mt=(np.pi*0.5)-omega_bar

        #True anomaly at the time of the transit
        #    - corresponds to 'dt_transit' (in years), time from periapsis to transit center
        #    - atan(X) is in -PI/2 ; PI/2
        ecc_anom_mt=2.*np.arctan(np.tan(true_anom_mt*0.5)*np.sqrt((1.-ecc)/(1.+ecc)))
        mean_anom_mt=ecc_anom_mt-ecc*np.sin(ecc_anom_mt)
        if (mean_anom_mt<0.):
            mean_anom_mt=mean_anom_mt+2.*np.pi

        #Mean anomaly
        #  - time origin of t_mean at the periapsis (t_mean=0 <-> M=0 <-> E=0)
        #  - M(t_mean)=M(dt_transit)+M(t_simu)
        mean_anom=2.*np.pi*phases+mean_anom_mt

        #Eccentric anomaly :
        #  - M = E - e sin(E)
        #    - >0 counterclockwise
        #  - angle, with origin at the ellipse center, between the major axis toward the periapsis and the line crossing the circle with radius 'a_Rs' at its intersection with the perpendicular to the major axis through the planet position
        ecc_anom=newton(func_kepler,mean_anom,args=(mean_anom,ecc,))

        #True anomaly of the planet at current time
        true_anom=2.*np.arctan(np.sqrt((1.+ecc)/(1.-ecc))*np.tan(ecc_anom/2.))
    return true_anom,ecc_anom

def  calc_planet_pos(sma_Rs, ecc, omega, inclin, l_spinorbit, Rp_Rs, orb_p, transitC, flag, step_grid, exposure_times=0.):
    """
    Takes the stellar and planet parameters as input and calulates the path of the planet
    in front of the cartesian stellar coordiante system
    input:
        sma_Rs: type: float, scaled semi major axis in solar radii
        ecc: type: float, eccentricity
        omega: type: float, Angle between the ascending node and the periastron, in the orbital plane (>0 counterclockwise)
        inclin: type: float, Inclination from the line of sight toward the normal to the orbital plane
        l_spinorbit: Orbit obliquity in degrees.
        Rp_Rs: type: float, ratio planet to star radii
        orb_p: type: float, orbital period in days
        transitC: type:float, transit center - 2400000.
        flag: type: string, can either be "times" for input as a real observation time array, or "phases" for an array of phases
        step_grid: type: np.array, if flag=="times", grid of times in jdb, if flag=="phases" grid of phases
        exposure_times: type: np.array, all exposure times, only needed if flag=="times", default 0.

    output: x_pl, y_pl, z_pl: type: np.arrays of floats, containing the position of the planet in units of stellar radii

    """
    import numpy as np


    inclin_bar = inclin*np.pi/180.
    omega_bar = omega*np.pi/180.
    obs_n = len(step_grid) #number of steps

    positions= np.empty([3,obs_n], dtype=float)

    if flag == "phases":
        phase = step_grid
    elif flag=="times":
        phase = np.zeros(obs_n,float)
        for i in range(obs_n):
            mid_time = calc_orbit_times(step_grid[i], transitC, exposure_times[i], orb_p)
            if mid_time<orb_p/2.:
                phase[i] = mid_time/orb_p
            else:
                phase[i] = (mid_time-orb_p)/orb_p
    else:
        raise Exception("Wrong flag option in calc_planet_pos. "
        "Possible values are 'phases' or 'times'.")


    #calc anomalies
    true_anom,ecc_anom = calc_true_anom(ecc,phase,omega_bar)

    #circular orbit
    if np.isclose(ecc,0.,1e-4):
        x_pl = sma_Rs*np.sin(np.asarray(true_anom))
        y_pl = -sma_Rs*np.cos(np.asarray(true_anom))*np.cos(inclin_bar)
        z_pl = sma_Rs*np.cos(np.asarray(true_anom))*np.sin(inclin_bar)
    #eccentric orbit
    else:
        #planet position in the orbital plane in Rstar
        X0_p = sma_Rs*(np.cos(ecc_anom)-ecc)
        Y0_p = sma_Rs*np.sqrt(1.-ecc*ecc)*np.sin(ecc_anom)
        #turn plane towards observer
        X1_p = X0_p*np.sin(omega_bar) + Y0_p*np.cos(omega_bar)
        Y1_p = -X0_p*np.cos(omega_bar) + Y0_p*np.sin(omega_bar)
        #translate to planet pos
        x_pl = Y1_p
        y_pl = -X1_p*np.cos(inclin_bar)
        z_pl = X1_p*np.sin(inclin_bar)
    return(x_pl*np.cos(np.radians(l_spinorbit))-y_pl*np.sin(np.radians(l_spinorbit)),
    x_pl*np.sin(np.radians(l_spinorbit))+y_pl*np.cos(np.radians(l_spinorbit)), z_pl)
    #CONVERT x_p and y_p to perpendicular x,y wrt stellar spin axis, see equation 4,5 of Cegla+ 2016.

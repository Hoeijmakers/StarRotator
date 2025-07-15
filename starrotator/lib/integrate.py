
# CONTINUE HERE NEXT TIME

def sum_stellar_spectrum_v1(wl,fx,vel_grid,flux_grid):
    """This builds the stellar spectrum by doppler-shifting and interpolating the input spectrum.
        This version of the integration assumes a FREE VELOCITY GRID (that means: with drr)
        but a MU-INDEPENDENT (static) stellar spectrum. The only center-to-limb variation is 
        modelled by broad-band limb darkening.

        Parameters
        ----------
        wl : np.array()
            The stellar model wavelengths.

        fx : np.array()
            The stellar model fluxes corresponding to wl.

        vel_grid : array-like
            The 2D velocity grid of the stellar disk. Generally this is a circular map in a square matrix.
            Values outside of the disk (in the corners of the square) are assumed to be set to NaN.

        flux_grid : array-like
            The 2D broad-band flux map of the stellar disk. Should have the same
            dimensions and axes as vel_grid.

        Returns
        -------
        wl : array
            The wavelength of the integrated spectrum of the same unit as the input.

        fx : array
            The flux axis of the summed spectrum corresponding to the wavelength points.
    """
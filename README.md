# StarRotator

StarRotator is a package that simulates a rotation-broadened stellar spectrum during an exoplanet transit event. The simulation is done via numerical integration of the stellar disk, with model photosphere spectra from either PHOENIX or SPECTRUM.

#### This is a minimal guide for for getting StarRotator to run.

StarRotator will run out of the box using a default exoplanet system defined in prepackaged configuration files (the demo .txt files included in the `./input/` folder).
In order to run StarRotator, perform the following steps:
1) Clone the repo / pull it to your machine.
2) To begin, open python in the root folder and hit `from StarRotator import StarRotator`.
3) Now StarRotator can be called as `KELT9 = StarRotator(586.0,592.0,200.0)` for a model spectrum of the Na-D lines of KELT-9, a fast-rotating 10,000K A-star orbited by KELT-9b, computed on a grid of (2x200)x(2x200) i.e. 400x400 square pixels.
4) Access the simulation output using attributes defined on the `KELT9` object: The wavelength axis of the model is accessed as `wl = KELT9.wl`, the out-of-transit stellar spectrum as `F_out = KELT9.stellar_spectrum`, the time-series of the modelled spectrum as `spectra = KELT9.spectra` and the residuals as obtained by dividing the out-of-transit spectrum out of the time-series: `residuals = KELT9.residual`. These can be blurred to some spectral resolution defined in `KELT9.R` (defined by default as 115,000 in the `demo_star.txt` parameter file), using `KELT9.convolve_spectral_resolution()`.
5) The StarRotator object contains methods to plot the simulation output. You can instantly plot the residuals of the time series like `KELT9.plot_residuals()`, and an animation of the entire time-series as `KELT9.animate()`, the result of which could look like the animation below:
6) To use the simulation output, in a typical use case the user would write the `KELT9.residual` numpy array to file (a fits file, numpy binary file, pickle, ascii table, etc.), and load it when needed in an external workflow.
7) To change the parameters of your exoplanet system you can create new input files. The parameters in the input files need to be given in the exact same order as in the demo files. Alternatively, you can call the StarRotator object with a dictionary that contains all these input values. `Planet = StarRotator(586.0,592.0,200.0,star_path='input/demo_star.txt',planet_path='input/demo_planet.txt',obs_path='input/demo_observations.txt')` if calling from input parameter files, or `Planet = StarRotator(586.0,592.0,200.0,input=dict)`, with `dict` set to a dictionary containing all the following named parameter values: `veq` (equatorial velocity, float), `stelinc` (stellar inclination axis, float), `drr` (differential rotation, float), `T` (Teff, float), `FeH` (metallicity, float), `logg` (float), `u1` (limb darkening parameter 1), `u2` (limb darkening parameter 2), `mus` (number of mu angles, int), `R` (resolving power, float), `model` (model type, string, either PHOENIX or pySME), `sma_Rs` (a over Rs, float), `e` (eccentricity, float), `omega` (longitude of periastron, float), `inclination` (float), `obliquity` (float), `RpRs` (float), `P` (orbital period, float), `mode` (string, set either to phases or times), `Tc` (transit center time, float, mode set to times), `phases` (numpy array, set to the orbital phase values of the time series).


![](demo.gif)

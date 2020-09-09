# StarRotator

StarRotator is a package that simulates a rotation-broadened stellar spectrum during an exoplanet transit event. The simulation is done via numerical integration of the stellar disk, with model photosphere spectra from either PHOENIX or SPECTRUM.

#### This is a minimal guide for for getting StarRotator to run.

In order to run StarRotator, perform the following steps:
1) Clone the repo / pull it to your machine. This will include the files to run SPECTRUM as well. If you are happy running StarRotator without SPECTRUM (i.e. without specific intensity-model spectra provided by the SPECTRUM code), you can skip ahead to step 14.
2) If not, navigate a terminal to the lib/SPECTRUM folder.
3) Type `make all` to compile SPECTRUM. It should compile with a few warnings.
4) Type `./spectrum`
5) Spectrum will open and prompt you for answers. Answer the following:
6) `vega.mod`
7) `merged_linelist.lst`
8) `../../temp`
9) `2.0`
10) `5000,5100`
11) `0.02`
12) Now, SPECTRUM should be running, while printing the output wavelength points to screen. You can see if the file called `temp` was created two folders up. This file should contain the spectrum.
13) If this succeeded, you can delete the `temp` file and are now ready to run StarRotator, with or without CLV enabled.


StarRotator will run out of the box using a not-so-fictional exoplanet system defined in prepackaged configuration files (the demo .txt files included in the root folder).

14) To begin, open python in the root folder and hit `from StarRotator import StarRotator`.
15) Now StarRotator can be called as `KELT9 = StarRotator(586.0,592.0,400.0)` for a model spectrum of the Na-D lines of KELT-9, a fast-rotating 10,000K A-star orbited by KELT-9b, computed on a grid of (2x400)x(2x400) i.e. 800x800 square pixels, and no CLV.
16) Access the simulation output using attributes defined on the `KELT9` object: The wavelength axis of the model is accessed as `wl = KELT9.wl`, the out-of-transit stellar spectrum as `F_out = KELT9.stellar_spectrum`, the time-series of the modelled spectrum as `spectra = KELT9.spectra` and the residuals as obtained by dividing the out-of-transit spectrum out of the time-series: `residuals = KELT9.residual`. These can be blurred to a spectral resolution defined in KELT9.R (defined by default as 115,000 in the `demo_star.txt` parameter file), using `KELT9.convolve_spectral_resolution()`.
17) The StarRotator object contains methods to plot the simulation output. You can instantly plot the residuals of the time series like `KELT9.plot_residuals()`, and an animation of the entire time-series as `KELT9.animate()`, the result of which could look like the animation below:

![](demo.gif)

# StarRotator

StarRotator is a package that simulates a rotation-broadened stellar spectrum during an exoplanet transit event. The simulation is done via numerical integration of the stellar disk, with model photosphere spectra from either PHOENIX or pySME.

#### This is a minimal guide for getting StarRotator to run.

StarRotator will run out of the box using a default exoplanet system defined in prepackaged configuration files (the demo `.txt` files included in the `./input/` folder).
In order to run StarRotator, perform the following steps:
1) Clone the repo / pull it to your machine.
2) To begin, open python in the root folder and hit `from StarRotator import StarRotator, testStarRotator`.
3) Run `testStarRotator()` to test that StarRotator functions normally. Functionality with PHOENIX models is possible even if you do not have pySME installed.
4) StarRotator can be called as `KELT9 = StarRotator(586.0,592.0,200.0)` for a model spectrum of the Na-D lines of KELT-9, a fast-rotating 10,000K A-star orbited by KELT-9b, computed on a grid of (2x200)x(2x200) i.e. 400x400 square pixels.
5) Access the simulation output using attributes defined on the `KELT9` object: The wavelength axis of the model is accessed as `wl = KELT9.wl`, the out-of-transit stellar spectrum as `F_out = KELT9.stellar_spectrum`, the modelled time-series as `spectra = KELT9.spectra` and the residuals as obtained by dividing the out-of-transit spectrum out of the time-series: `residuals = KELT9.residual`. These can be blurred to some spectral resolution defined in `KELT9.R` (defined by default as 115,000 in the `demo_star.txt` parameter file), using `KELT9.convolve_spectral_resolution()`. Convolution is applied to the residuals and not to the spectra that are divided by each other, because this introduces a numerical error.
6) The StarRotator object contains methods to plot the simulation output. You can instantly plot the residuals of the time series like `KELT9.plot_residuals()`, and an animation of the entire time-series as `KELT9.animate()`, the result of which could look like the animation below.
7) To use the simulation output, the user has access to the `KELT9.residual` numpy array, which can be written to file and loaded it when needed in an external workflow, or used as is.
8) To change the parameters of your exoplanet system you can create new input files. The parameters in the input files need to be given in the exact same order as in the demo files. Alternatively, you can call the StarRotator object with a dictionary that contains all these input values. See the example code block below':

```python
#When using the input files:
Planet = StarRotator(586.0,592.0,200.0,star_path='input/demo_star.txt',planet_path='input/demo_planet.txt',obs_path='input/demo_observations.txt')


#When using the dictionary:
dict = {
  'veq':114000.0, #
  'stelinc':90.0, #stellar inclination axis, degrees, float
  'drr':0.0, #Differential rotation parameter, float
  'T':10000.0, #Stellar effective temperature, K, float
  'FeH':0.0, #metallicity, float
  'logg':4.0, #log(g), cgs, float
  'u1':0.93, #limb darkening parameter 1
  'u2':-0.23, #Spectral resolution
  'R':115000, #limb darkening parameter 2
  'mus':0, #number of mu angles, int. Ignored if using PHOENIX.
  'model':'PHOENIX', #model type, string, either PHOENIX or pySME
  'sma_Rs':3.153, #a over Rs, float
  'e':0.0, #eccentricity, float
  'omega':0.0, #longitude of periastron, degrees, float
  'inclination':86.79, #degrees, float
  'obliquity':-84.8, #degrees, float
  'RpRs':0.08228, #Planet star radius ratio, float
  'P':1.4811235, #Orbital period, days, float
  'phases':[-0.02,-0.01,0.0,0.01,0.02] #numpy array, set to the orbital phases of the time series
}

Planet = StarRotator(586.0,592.0,200.0,input=dict)
```

![](demo.gif)

<br><br><br>

### A typical use case: running StarRotator with pySME as a forward model.

[PySME](https://github.com/AWehrhahn/SME) is used in StarRotator to allow a user to generate more precise forward-models of the Doppler-Shadow residuals for real exoplanet systems. Using known stellar parameters (T, log(g), Z, vsin i) and individual elemental abundances, in theory it should be possible to precisely forward model the residual and recalibrate observed spectra. To the broadened spectrum in the Na lines in a sodium-rich variation on the KELT-9 system with pySME, you can use the following example.


```python

from StarRotator import StarRotator
import matplotlib.pyplot as plt
import numpy as np

dict = {
  'veq':114000.0,'stelinc':90.0,'drr':0.0,
  'T':10500.0,'FeH':0.23,'logg':3.9,'u1':0.93,'u2':-0.23,
  'R':115000,'mus':5,'model':'pySME',
  'sma_Rs':3.153,'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,
  'P':1.4811235,'phases':np.arange(-0.03,0.03,40),'grid_model':'atlas12.sav','abund':{},
  'linelist_path':'input/demo_linelist.dat'
}

KELT9_nominal = StarRotator(586.0,592.0,200.0,input=dict)

dict['abund']=["{'Na':6.6}"]#Nominally it is 6.4
KELT9_rich = StarRotator(586.0,592.0,200.0,input=dict)

plt.plot(KELT9_nominal.wl,KELT9_nominal.stellar_spectrum,label='Nominal sodium')
plt.plot(KELT9_rich.wl,KELT9_rich.stellar_spectrum,label='Sodium enhanced')
plt.legend(frameon=False)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Relative flux')
plt.show()

```

This creates the following figure of the out-of-transit broadened stellar sodium lines. Adding extra broadening to account for the instrumental resolving power and  accessing the residuals can subsequently be done via `KELT9_rich.convolve_spectral_resolution()` and `KELT9_rich.residuals()`. In addition, in reality you will wish to use more complete line-lists (in VALD format) and data-driven stellar and system parameters. All this is left for the user to explore.
![](demo_spectrum.png)

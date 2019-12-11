# StarRotator



This is a minimal guide for for getting StarRotator to run.

#### Necessary steps to run StarRotator with CLV (i.e. <img src="https://latex.codecogs.com/gif.latex?\mu" /> -resolved photospheric spectra) enabled

In order to run StarRotator with specific intensity-model spectra provided by the SPECTRUM code, perform the following steps:
1) Clone the repo / pull it to your machine. This will include the files to run SPECTRUM as well. If you are happy running StarRotator without SPECTRUM, you can skip ahead past step 13.
2) If not, navigate a terminal to the lib/SPECTRUM folder
3) Type `make all` to compile SPECTRUM. I should compile with a few warnings.
4) Type `./spectrum`
5) Spectrum will open and prompt you for answers. Answer the following:
6) `vega.mod`
7) `merged_linelist.lst`
8) `../../temp`
9) `2.0`
10) `5000,5100`
11) `0.02`
12) Now, spectrum should be running, while printing the output wavelength points to screen. You can see if the file called temp was created two folders up. This file contains the spectrum. 
13) If this succeeded, you are ready to run StarRotator with CLV enabled. 


StarRotator will run out of the box using a fictional exoplanet system defined in prepackaged configuration files (the demo .txt files included in the root folder). To begin, open python in the root folder and `import main`.

14) Now StarRotator can be called as `out = main.StarRotator(586.0,592.0,400.0)` for a spectrum of the Na D lines, a fast-rotating star and a grid of (2x400)x(2x400) ie 600x600 pixels, and no CLV.

`out` is a massive tuple containing tons of output. Need to put this into a class in order to operate the plotting routines that are locked and loaded for business.





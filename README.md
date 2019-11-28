# StarRotator



Here comes a minimal Readme for getting main.py to run.

1) Clone the repo / pull it to your machine. This will include the files to run SPECTRUM as well.
2) Navigate a terminal to the lib/SPECTRUM folder
3) Type "make all". I should compile with a few warnings.
4) Type ./spectrum
5) Spectrum will open and prompt you for answers. Answer the following:
6) vega.mod
7) merged_linelist.lst
8) ../../temp
9) 2.0
10) 5000,5100
11) 0.02
12) Now, spectrum should run, printint the output wavelength points to screen. You can see if the file called temp was created
two folders up. This file contains the spectrum. You are now ready to run StarRotator, as follows:
13) python3 main.py 586.0 592.0 110000.0 90.0 90.0 0.0 0.0 300



Julia, right now some settings that determine the behaviour are hardcoded in main.py.
For example, there is a sys.exit() at the end of the decision tree, before the animation is made. You may want to remove it.
Also, you may want to set the parameter 'mus' to either zero, or mus = np.linspace(0.0,1.0,20) (for 20 mu angles). You can play with the stellar parameters in general. The first time you run it, the code should ask you to download Kurucz models. Its only a few dozen MB I believe, no big stuff.



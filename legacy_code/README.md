#RamPy
=====

This is a set of python scripts that are used to treat spectroscopic x y data

The core is in the spectratool.py module, with various functions to help organising the data, subtracting baseline (polynomial, linear, cubic splines), with access to linear, polynomial or gaussian/lorentzian function that can be used with peak fitting routine for instance. 

the /examples/ folder contain various scripts that uses the tools developped to do various spectral treatment, including subtraction of the 2nd order diamond signal from Raman spectra, as described for instance in Dalou C., Le Losq C., Mysen B. O. (2015) In situ study of the fractionation of hydrogen isotopes between aluminosilicate melts and coexisting aqueous fluids at high pressure and high temperature - Implications for the dD in magmatic processes. Earth and Planetary Science Letters 426, 158-166.

Developped under Mac with Anaconda Python.

C. Le Losq.
charles.lelosq@anu.edu.au

# DISCLAIMER

RAMPY IS CURRENTLY UNDER A NEW DEVELOPMENT PHASE AND MAY NOT BE USABLE FOR A FEW DAYS.

# REQUIREMENTS

For Python:
Scipy
Numpy
Matplotlib
lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
f2py

Others:
gfortran

# INSTALLATION INSTRUCTIONS

1) Be sure to setup correctly you $PATH variable in your .bash_profile file in your home directory. It must point to the python distribution you want to use. Be carefull, usually several python distributions can be installed in your system... Anaconda Python set up by itself the PATH during its installation. If you have doubt, simply google "setting up the .bash_profile file for python"

2) You need to install lmfit and aeval, follow the instructions from this webpage: http://cars9.uchicago.edu/software/python/lmfit/
I forked also the lmfit project in my Github...

4) You need also gfortran, install it if you don't have it;

5) Put the folder with the script where you want; Some files have to be moved and a particular folder construction have to be respected in order to work... Please read the following to learn about that...

6) You need to compile the gcvspl.f algorithm in the ./RamPy/src/gcvspl/ folder using the following command:

f2py -c -m gcvspl gcvspl.f

7) read the following instructions and descriptions to know what is useful for what!


# EXAMPLES

##### Diamond_Pressure.py

Can be placed in the folder you want... 

This script ask you the files for calculating the exact position of the 13C in your sample at room and high temperatures. It fits the diamond and neon peaks, then calculates the exact frequency of the diamond peaks and uses the equation provided in 
[Mysen and Yamashita, 2010. Speciation of reduced C-O-H volatiles in coexisting fluids and silicate melts determined in-situ to similar to 1.4 GPa and 800 degrees C. GCA 74, 4577-4588] to calculate the pressure.

Type P in the console to read the output.

##### HT_Treatment.py

Can be placed in the folder you want... 

This script is for subtracting the 2nd order of diamond from the Raman spectra of samples from HT/HP DAC experiments.

First, the script asks you the sample and diamond spectrum (acquired at the same P/T conditions)

Then it apply a spline calculated by the gcvspl.f algorithm. You can also use the univariate spline from scipy by changing the lines 69
and 73. Here the details:

Lines 69 and 73: select "unispline" or "gcvspline" in the linbaseline() function if you want the univariate or gcvspl splines. I advise the gcvspl, but univariate gives also good results.

Same lines: set up the spline factor at the end of the linbaseline() function. WARNING: i don't know why but univariate ask for quite high values of the spline factor to work, typically close to 10000000....

Then the software subtracts the baseline. It corrects any shift in frequency with using that of the ~2450 diamond peak, and then we assume the following:

Diamond_in_sample = K * Diamond.

We calculate K by adjusting by least-square Diamond_in_sample and Diamond in the 2140 - 2340 cm-1 region, where we assume that NO SAMPLE SIGNAL IS PRESENT. You can change those values if your experiments are different in the lines 112 and 113.

After line 129: the script is set up to calculate the D-O and H-O stretching areas in our case.

Last remark: Filters are also present lines 139 and 140 as well as line 75 go 77 to erase some of the noise in the data. You can comment them if not needed of course. You can also adjust the cutoff frequency if you need them.

Last step: Ask where it output a file containing the x and corrected y.

The figure is not saved.

##### LF_Treatment.py

Needs to follow a particular repertory structure: you have to put that in a ./root/ repertory (name root can be changed of course but it is for the example)

then put the spectra in ./root/raw for instance.

Write a list in a text file containing their name as:

./raw/name1.txt
./raw/name2.txt
./raw/name3.txt

create a ./root/treated/ directory (see line 132 if you want to change this name).

The script will ask you to open the list of name, it will read it, take the name, gives them to the function removebas_melt or removebas_fluid that you select line 128. Then, it outputs the a figure in .pdf and the new spectra in the ./root/treated/ directory.

This script fit the baseline of spectra in the 100-1500 cm-1 region. Two functions here for our local need. Select the function line 128. Adjust the baseline constrains in the functions line 58 or 98 and the spline coefficient line 59 or 99.

##### LF_dec.py

Deconvolution of Raman spectra with Gaussian and/or Lorentzian bands
Use the lmfit framework developed by Matt Newville for the curve_fit algorithm of Scipy 

##### IR_dec_comb.py and HT_dec 
Deprecated, need to be finished...






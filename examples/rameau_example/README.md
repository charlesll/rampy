# RAMEAU

Rameau is a set of functions and a class from the rampy library that provides a way to calculate the water content of a glass from its Raman spectrum.

At this stage, it provides access to external Raman calibrations, from the Le Losq et al. (2012) and Di Genova et al. (2017) articles.

# INSTALLATION

## On your machine

Works well on Mac OS or Linux (also very well in a virtual machine)

Problems with gcvspline in Windows can occur, due to the requirement of a FORTRAN compiler, but this should be solved as gcvspline has been precompiled for Python 2.7, 3.4, 3.5, 3.6, 3.7 and 3.8 for 32 and 64 bits Windows systems.

- Install a Python stack (see Anaconda Python for instance) with scipy, numpy, matplotlib, pandas, xlrd :

```
$ pip install --user numpy pandas scipy matplotlib xlrd
```

- Install rampy and gcvspline using the command in a terminal (in your working virtual environment) :

```
$ pip install rampy gcvspline
```

## Using a Vagrant box (easy)

This is the easy path that allows one to run rameau in a virtual environment, without having to tune a Python install.

- Install VirtualBox on your computer, see https://www.virtualbox.org/

- Install Vagrant, see https://www.vagrantup.com/downloads

- Point a terminal in your rameau working directory

- run the following commands in this terminal :

```
$ vagrant init charlesll/rameau
$ vagrant up
$ vagrant ssh
```

The first command will download the box, the second will start it, and the third will create a tunnel into the virtual machine.

The virtual machine will see everything in the working directory, so after you just need to run the below rameau commands to start working!

# USAGE

## Jupyter notebook or rameau.py script ?

- the Python_treatment.ipynb Jupyter notebook provides some interactive example with figures showcasing the different calibrations.

- however, the easiest usage is to work with the rameau.py script provided in this folder. This allows calling rameau directly from the terminal, providing a simple set of arguments.

## Folder prepration

- Put rameau.py or the notebook in your working directory (e.g. you can download this example directory and work in it)

- Spectra should be in a subfolder called "raw".

- Use the provided dataset.xlsx spreadsheet as a reference : it contains the name of the spectra, if you are calibrating your system, the glass standards' water contents, then iron contents (in all cases), and the ROIs for the Le Losq et al. (2012) method (abreviated LL2012). The Di Genova et al. method (DG2017) requires the iron contents (in wt%, total FeO). You can add sheets to the spreadsheet for your samples.

- then run your notebook or use the python rameau command as described in the next section.

# COMMANDS

```
$ python rameau.py -i <excelfile> -s <sheetname> -o <outputfile> -l <laser_wavelength> -c1 <LL2012_coefficients> -c2 <DG2017_coefficients>
```
-i or --ifile :

Input excel file, in the same directory as the rameau.py script. See the starting example for construction. See also Le Losq et al. (2012) for the choice of the ROIs.

-s or --sheetname :

Sheet that needs to be considered in the Excel file.

-o or --outputfile :

Output file name, csv format, where spectra names and predicted water contents will be stored.
Example : -o results.csv

Figures will be outputed automatically in folders ./figures/LL2012 and ./figures/DG2017/

Please note that any change can overwrite previous figures in those folders.

-l or --laser :

Laser wavelength that you are using.

-c1 or --coefs_LL2012 :

Coefficient of the Le Losq et al. 2012 calibration for your setup (to be determined using glass standards on a given Raman system).

-c2 or --coefs_DG2017 :

Coefficients of the Di Genova et al. 2017 calibration for your setup (to be determined using glass standards on a given Raman system).

NOTE : YOU NEED TO PROVIDE IRON CONTENTS IN THE INPUT EXCEL FILE FOR THE DG 2017 METHOD.

-m or --mode :

Either "calibration" or "prediction".

If "calibration", the software assumes you know the water contents of the glasses (stored in the input file) and determines the coefficients for the LL2012 and DG2017 methods. Those coefficients will be provided in the output file.

If "prediction", the software will use the coefficients provided (see c1 and c2 options) and will return the water contents in the output file.

# EXAMPLE

1) To predict things for the 2012 Le Losq et al. dataset, you can use the command :

```
$ python rameau.py -i dataset.xlsx -s 2012_AmMin -o results.csv -l 514.532 -m calibration
```

Please note that your terminal path should be set in the working folder.

2) same example but with a dataset from the Research School of Earth Sciences :

```
python rameau.py -i dataset.xlsx -s 2018_RSES -o results.csv -l 532.0 -m calibration
```

3) Prediction example, we will use the 2018_RSES for this example and provide previously calculated coefficients

```
$ python rameau.py -i dataset.xlsx -s 2018_RSES -o results.csv -l 532.0 -m prediction -c1 [0.009208] -c2 [-0.00326133  1.32225612]
```

# References

C. Le Losq, D. R. Neuville, R. Moretti, J. Roux, Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist. 97, 779–790 (2012).

D. Di Genova et al., Effect of iron and nanolites on Raman spectra of volcanic glasses: A reassessment of existing strategies to estimate the water content. Chemical Geology. 475, 76–86 (2017).

########################################
############ INSTALLATION ##############
########################################

Works best on Mac OS or Linux 

Problems with gcvspline in Windows are frequent, due to the requirement of a FORTRAN compiler.

Install a Python stack (see Anaconda Python for instance) with scipy, numpy, matplotlib, pandas

Install rampy and gcvspline using the command in a terminal (in your working virtual environment) :

> pip install rampy
> pip install gcvspline

########################################
################ USAGE #################
########################################

python rameau.py -i <excelfile> -s <sheetname> -o <outputfile> -l <laser_wavelength> -c1 <LL2012_coefficients> -c2 <DG2017_coefficients>

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

########################################
################ EXAMPLE ###############
########################################

1) To predict things for the 2012 Le Losq et al. dataset, you can use the command :

>>> python rameau.py -i dataset.xlsx -s 2012_AmMin -o results.csv -l 514.532 -m calibration

Please note that your terminal path should be set in the working folder.

2) same example but with a dataset from the Research School of Earth Sciences :

>>> python rameau.py -i dataset.xlsx -s 2018_RSES -o results.csv -l 532.0 -m calibration

3) Prediction example, we will use the 2018_RSES for this example and provide previously calculated coefficients

>>> python rameau.py -i dataset.xlsx -s 2018_RSES -o results.csv -l 532.0 -m prediction -c1 [0.009208] -c2 [-0.00326133  1.32225612]

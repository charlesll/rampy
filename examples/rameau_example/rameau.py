#!/usr/bin/env python
# coding: utf-8

########## Calling relevant libraries ##########
import numpy as np
import scipy, os, matplotlib, sys, getopt

from matplotlib import pyplot as plt

import pandas as pd
import rampy as rp

def main(argv):

    inputfile = None
    sheetname = None
    outputfile = None
    my_laser = None
    coefficient_calibration_LL2012 = [0.009296]
    coefficient_calibration_DG2017 = [-0.003261, 1.322256]
    mode = None
    
    try:
        opts, args = getopt.getopt(argv,"hi:s:o:l:c1:c2:m:",
                                 ["ifile=","sheetname=","ofile=","laser=","coefs_LL2012=","coefs_DG2017=","mode="])
    except getopt.GetoptError:
        print('rameau.py -i <excelfile> -s <sheetname> -o <outputfile> -l <laser_wavelength> -c1 <LL2012_coefficients> -c2 <DG2017_coefficients> -m <mode>')
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print('rameau.py -i <excelfile> -s <sheetname> -o <outputfile> -l <laser_wavelength> -c1 <LL2012_coefficients> -c2 <DG2017_coefficients> -m <mode>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-s", "--sheetname"):
            sheetname = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-l", "--laser"):
            laser = arg
        elif opt in ("-c1", "--coefs_LL2012"):
            coefficient_calibration_LL2012 = arg
        elif opt in ("-c2", "--coefs_DG2017"):
            coefficient_calibration_DG2017 = arg
        elif opt in ("-m", "--mode"):
            mode = arg

    print('')                              
    print('Input excel file is ', inputfile)
    print('Sheet selected is ', sheetname)
    print('Output file is ', outputfile)
    print('Selected laser wavelength is', laser)
    print('Selected mode is ', mode)
    print('')
    
    #    
    # DATA TREATMENT
    #
    liste_standards = pd.read_excel(inputfile, sheet_name=sheetname)

    # treatment of the dataset with both methods
    # ll_2012 and df_2012 are two objects with associated function. see documentation
    ll_2012 = rp.rameau(liste_standards)
    ll_2012.data_reduction(method="LL2012",laser=laser)

    dg_2017 = rp.rameau(liste_standards)
    dg_2017.data_reduction(method="DG2017",laser=laser)
    
    #
    # PLOTTING THE BASELINES
    #

    dirName = 'figures/LL2012/' # Create directory for storing figures. Change dirName accordingly
    try:
        os.makedirs(dirName) # Create target Directory
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

    dirName = 'figures/DG2017/'
    try:
        os.makedirs(dirName) # Create target Directory
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

    print('')
    # using the plot_spectra function...
    plot_spectra(ll_2012,path_out="./figures/LL2012/")
    plot_spectra(dg_2017,path_out="./figures/DG2017/")

    # 
    # MAKING PREDICTIONS OR SETTING THE CALIBRATION
    #                      
    if mode == "prediction":
        # if we want to make predictions, we use the provided coefficients
        ll_2012.p = np.array(coefficient_calibration_LL2012)
        dg_2017.p = np.array(coefficient_calibration_DG2017)

        ### WE MAKE PREDICTIONS USING THE PROVIDED COEFFICIENTS
        ll_2012.predict(method="LL2012") 
        dg_2017.predict(method="DG2017")

        out = pd.DataFrame()
        out["spectrum"] = liste_standards.Name
        out["LL2012_waterpredicted"] = ll_2012.water_predicted
        out["DG2017_waterpredicted"] = dg_2017.water_predicted
        out.to_csv(outputfile)

    elif mode == "calibration":
        ll_2012.calibrate(method="LL2012")
        ll_2012.predict(method="LL2012")

        dg_2017.calibrate(method="DG2017")
        dg_2017.predict(method="DG2017")

        # we print the calculated coefficients in a file named coefficients.txt
        file = open("coefficients.txt","w")
        file.write(f'Coefficient for the Le Losq et al. 2012 method is {ll_2012.p}\n')
        file.write(f'Coefficients for the Di Genova et al. 2017 method are {dg_2017.p}\n')
        file.close

        out = pd.DataFrame()
        out["spectrum"] = liste_standards.Name
        out["standard_water_content"] = ll_2012.water
        out["LL2012_waterpredicted"] = ll_2012.water_predicted
        out["DG2017_waterpredicted"] = dg_2017.water_predicted
        out.to_csv(outputfile)

    else:
        raise ValueError("mode should be set to calibration or prediction. Exiting...")
    
# Function definitions
def plot_spectra(rameau_object,**kwargs):

    idx = kwargs.get("idx",False)
    path_out = kwargs.get("path_out",[])

    if idx is not False:
        plt.figure(figsize=(6.5,3.25))
        ax1 = plt.subplot(1,2,1)
        ax1.plot(rameau_object.x,rameau_object.y[:,idx],'k.',ms=0.5)
        ax1.plot(rameau_object.x,rameau_object.y_base[:,idx],'r-')
        ax1.set_xlim(0,1500)
        ax1.set_ylim(0,np.max(rameau_object.y[np.where(rameau_object.x<1500),idx])+0.1*np.max(rameau_object.y[np.where(rameau_object.x<1500),idx]))

        ax1.set_xlabel(r'Raman shift, cm$^{-1}$')
        ax1.set_ylabel(r'Intensity, a. u.')

        ax2 = plt.subplot(1,2,2)
        ax2.plot(rameau_object.x,rameau_object.y[:,idx],'k.',ms=0.5)
        ax2.plot(rameau_object.x,rameau_object.y_base[:,idx],'r-')
        ax2.set_xlim(2700,4000)
        ax2.set_ylim(0,np.max(rameau_object.y[np.where(rameau_object.x>2700),idx])+0.1*np.max(rameau_object.y[np.where(rameau_object.x<4000),idx]))

        ax2.set_xlabel(r'Raman shift, cm$^{-1}$')
        ax2.set_ylabel(r'Intensity, a. u.')

        plt.tight_layout()

    else:
        print("Save spectra in provided path: %s"%path_out)
        for i in range(rameau_object.y.shape[1]):
            plt.figure(figsize=(6.5,3.25))
            ax1 = plt.subplot(1,2,1)
            ax1.plot(rameau_object.x,rameau_object.y[:,i],'k.',ms=0.5)
            ax1.plot(rameau_object.x,rameau_object.y_base[:,i],'r-')
            ax1.set_xlim(0,1500)
            ax1.set_ylim(0,np.max(rameau_object.y[np.where(rameau_object.x<1500),i])+0.1*np.max(rameau_object.y[np.where(rameau_object.x<1500),i]))

            ax1.set_xlabel(r'Raman shift, cm$^{-1}$')
            ax1.set_ylabel(r'Intensity, a. u.')

            ax2 = plt.subplot(1,2,2)
            ax2.plot(rameau_object.x,rameau_object.y[:,i],'k.',ms=0.5)
            ax2.plot(rameau_object.x,rameau_object.y_base[:,i],'r-')
            ax2.set_xlim(2700,4000)
            ax2.set_ylim(0,np.max(rameau_object.y[np.where(rameau_object.x>2700),i])+0.1*np.max(rameau_object.y[np.where(rameau_object.x<4000),i]))

            ax2.set_xlabel(r'Raman shift, cm$^{-1}$')
            ax2.set_ylabel(r'Intensity, a. u.')

            plt.tight_layout()
            plt.savefig(path_out+rameau_object.names[i]+".pdf")
            plt.close()
            
if __name__ == "__main__":
    main(sys.argv[1:])





# coding: utf-8
# (c) Copyright Le Losq et al. 2018
# Supplementary material of Le Losq, Berry, Kendrick, Neuville, O'Neill, 
# Determination of the oxidation state of iron in basaltic glasses by Raman spectroscopy

import numpy as np
import rampy as rp
from numpy.linalg import inv

# functions
def preparing_data(dataliste,**kwargs):
    """prepare the spectra before processing by the regression techniques
    
    Parameters
    ==========
    dataliste : Pandas dataframe
        A liste containing the name of the spectra, located in a folder indicated in pathin
        
    Options
    =======
    pathin : string
        the path of the spectra. Default = './raw/'
    cutoff : ndarray
        frequencies delimiting the region of interest for the regression. Default = np.array([850.,1040.])
    scale : float
        scaling coefficient for the intensity. Default = 1000
        
    Returns
    =======
    x : ndarray 
        the x axis as np.arange(300,1290,1.0)
    record : ndarray
        the y signal corrected from temperature and excitation line effects (23 °C, 532 nm)
    record_bas2 : ndarray
        the baseline fitted to record
    x_cut : ndarray
        the x axis of the region of interest
    record_hf_no_smo : ndarray
        the y signal in teh region of interest, scaled between 0 and 1 (no smoothing)
    record_hf : ndarray
        the y signal in the region of interest, scaled between 0 and 1 and smoothed with a whittaker algorithm.
    nb_exp : int
        number of experiments (= length of dataliste)
        
    Note
    ====
    Input spectra are assumed to have decreasing frequencies. If not, comment the line `data = rp.flipsp(data)`
    """
    
    #
    # Kwargs
    #
    
    cutoff = kwargs.get("cutoff",np.array([850.,1140.])) # roi of the linear baseline
    scale = kwargs.get("scale",1000) # scaling coefficient
    pathin = kwargs.get('pathin',"./raw/")
    
    # a new x axis for interpolation (all spectra may not have been sampled with the same x)
    x = np.arange(300,1300,1.0)
    
    # for the baseline fit, we grabe two points
    roi_cutoff = np.array([[cutoff[0]-0.4,cutoff[0]+0.4],[cutoff[1]-0.4,cutoff[1]+0.4]])
    
    # number of spectra
    nb_exp = len(dataliste)
    
    # array to record the treated spectra
    record = np.ones((x.shape[0],nb_exp))
    record_bas2 = np.ones((x.shape[0],nb_exp))
    
    # loop to read the spectra
    for i in range(nb_exp):
        data = np.genfromtxt(pathin+dataliste["spectra"].iloc[i],skip_header=1)
        
        # we need an increasing x axis for the interpolators so we check this point
        data = rp.flipsp(data)
        
        # finding the minimum between 1200 and 1300 to fit a constant baseline (bas1)
        idx_roi = np.where(data[:,1] == np.min(data[(data[:,0]>1200)&(data[:,0]<1300),1]))[0][0]  
        roi_bas1 = np.array([[data[idx_roi,0] - 15.,data[idx_roi,0] + 15.]])
        y_bas1, bas1 = rp.baseline(data[:,0],data[:,1],roi_bas1,"poly",polynomial_order=0)

        # resampling
        y_norm = rp.resample(data[:,0],y_bas1[:,0],x)
        
        # correcting from temperature and excitation line effect; the tlcorrection function automatically normalize to the area.
        trash, y_long, trash = rp.tlcorrection(x,y_norm,23.0,532.0)
    
        record[:,i] = y_long[:]*scale #with a scale factor to bring values closer to 1 for representation
        
        # now grabbing the signal above the cutting baseline (bas2) in the roi_cutoff portion of spectra
        y_corr, bas2 = rp.baseline(x,y_long[:],roi_cutoff,"poly",polynomial_order=1.0)
        
        x_cut = x[(roi_cutoff[0,0]<=x)&(x<=roi_cutoff[1,1])].reshape(-1,1)
        y_cut = y_corr[(roi_cutoff[0,0]<=x)&(x<=roi_cutoff[1,1])].reshape(-1,1)
        
        # initialisation of output arrays for signal of interest
        if i == 0: 
            record_hf = np.ones((y_cut.shape[0],nb_exp))
            record_hf_no_smo = np.ones((y_cut.shape[0],nb_exp))
        
        # Getting the good signal at HF (above the cut-off baseline) + Min-Max scaling
        record_hf_no_smo[:,i]= (y_cut[:,0]-np.min(y_cut[:,0]))/(np.max(y_cut[:,0])-np.min(y_cut[:,0]))
         
        # smoothing the signal with a Whittaker smoother = improves results
        record_hf[:,i] = rp.whittaker(record_hf_no_smo[:,i],Lambda = 10.0**3)
        
        # wew take care of correcting any deviation from 0 after smoothing
        y_r_2, _ = rp.baseline(x_cut,record_hf[:,i],roi_cutoff,"poly",p=1.0)
        record_hf[:,i] = ((y_r_2-np.min(y_r_2))/(np.max(y_r_2)-np.min(y_r_2))).reshape(-1)
        
        # for the baseline
        record_bas2[:,i] = bas2[:,0]*scale
        
    return x, record, record_bas2, x_cut, record_hf_no_smo, record_hf, nb_exp

#
# Function to perform the intensity method
#

def intensity(x,y,idx1=929.,idx2=931.):
    """calculate the intensity of the y signal as the mean value of y between idx1 and idx2
    
    Parameters
    ----------
    x : ndarray, shape m
        the x values
    y : ndarray, shape m * n
        the y values, with n samples
    idx1 : float
        the low x boundary for mean calculation. Default = 929.0
    idx2 : float
        the high x boundary for mean calculation. Default = 931.0
        
    Returns
    -------
    out : ndarray, shape n
        the mean intensity between idx1 and idx2 for the n samples
    """
    inten = np.ones(y.shape[1])
    x = x.reshape(-1)# safety
    for i in range(y.shape[1]):
        y_1 = np.mean(y[(idx1< x)&(x<idx2),i])
        inten[i] = y_1
    return inten

def chimie_control(datalist):
    """check that all needed oxides are there and setup correctly the Pandas datalist.

    Parameters
    ----------
    datalist : Pandas dataframe
        the user input list.

    Returns
    -------
    out : Pandas dataframe
        the output list with all required oxides.

    """

    list_oxides = ["sio2","al2o3","tio2","fe2o3","h2o","li2o","na2o","k2o","mgo","cao","bao","feo","nio","mno","p2o5"]

    for i in list_oxides:
        try:
            oxd = datalist[i]
        except:
            datalist[i] = 0.

    if (datalist["sio2"] > 1).any(): # if values were in percents
        datalist["sio2"] = datalist["sio2"]/100.0
        datalist["al2o3"] = datalist["al2o3"]/100.0
        datalist["tio2"] = datalist["tio2"]/100.0
        datalist["fe2o3"] = datalist["fe2o3"]/100.0
        datalist["h2o"] = datalist["h2o"]/100.0
        datalist["li2o"] = datalist["li2o"]/100.0
        datalist["na2o"] = datalist["na2o"]/100.0
        datalist["k2o"] = datalist["k2o"]/100.0
        datalist["mgo"] = datalist["mgo"]/100.0
        datalist["cao"] = datalist["cao"]/100.0
        datalist["bao"] = datalist["bao"]/100.0
        datalist["feo"] = datalist["feo"]/100.0
        datalist["nio"] = datalist["nio"]/100.0
        datalist["mno"] = datalist["mno"]/100.0
        datalist["p2o5"] = datalist["p2o5"]/100.0

    return datalist

def molarweights():
	"""returns a partial table of molecular weights for elements and oxides that can be used in other functions

    Returns
    =======
    w : dictionary 
        containing the molar weights of elements and oxides:

        - si, ti, al, fe, li, na, k, mg, ca, ba, o (no upper case, symbol calling)

        - sio2, tio2, al2o3, fe2o3, feo, li2o, na2o, k2o, mgo, cao, bao (no upper case, symbol calling)

    """
	w = {}

	# From IUPAC Periodic Table 2016, in g/mol
	w["si"] = 28.085
	w["ti"] = 47.867
	w["al"] = 26.982
	w["fe"] = 55.845
	w["h"] = 1.00794
	w["li"] = 6.94
	w["na"] = 22.990
	w["k"] = 39.098
	w["mg"] = 24.305
	w["ca"] = 40.078
	w["ba"] = 137.33
	w["o"] = 15.9994

	w["ni"] = 58.6934
	w["mn"] = 54.938045
	w["p"] = 30.973762

	# oxides
	w["sio2"] = w["si"] + 2* w["o"]
	w["tio2"] = w["ti"] + 2* w["o"]
	w["al2o3"] = 2*w["al"] + 3* w["o"]
	w["fe2o3"] = 2*w["fe"] + 3* w["o"]
	w["feo"] = w["fe"] + w["o"]
	w["h2o"] = 2*w["h"] + w["o"]
	w["li2o"] = 2*w["li"] +w["o"]
	w["na2o"] = 2*w["na"] + w["o"]
	w["k2o"] = 2*w["k"] + w["o"]
	w["mgo"] = w["mg"] + w["o"]
	w["cao"] = w["ca"] + w["o"]
	w["bao"] = w["ba"] + w["o"]

	w["nio"] = w["ni"] + w["o"]
	w["mno"] = w["mn"] + w["o"]
	w["p2o5"] = w["p"]*2 + w["o"]*5
	return w # explicit return

def wt_mol(chemtable):
	"""to convert weights in mol fraction

	Parameters
	==========
	chemtable: Pandas DataFrame 
		containing the fields sio2,tio2,al2o3,fe2o3,li2o,na2o,k2o,mgo,cao,feo

	Returns
	=======
	chemtable: Pandas DataFrame
		contains the fields sio2,tio2,al2o3,fe2o3,li2o,na2o,k2o,mgo,cao,feo in mol%
	"""

	w = molarweights()
	# conversion to mol in 100 grammes
	sio2 = chemtable["sio2"]/w["sio2"]
	tio2 = chemtable["tio2"]/w["tio2"]
	al2o3 = chemtable["al2o3"]/w["al2o3"]
	fe2o3 = chemtable["fe2o3"]/w["fe2o3"]
	h2o = chemtable["h2o"]/w["h2o"]
	li2o = chemtable["li2o"]/w["li2o"]
	na2o = chemtable["na2o"]/w["na2o"]
	k2o = chemtable["k2o"]/w["k2o"]
	mgo = chemtable["mgo"]/w["mgo"]
	cao = chemtable["cao"]/w["cao"]
	feo = chemtable["feo"]/w["feo"]
	bao = chemtable["bao"]/w["bao"]

	mno = chemtable["mno"]/w["mno"]
	nio = chemtable["nio"]/w["nio"]
	p2o5 = chemtable["p2o5"]/w["p2o5"]
	# renormalisation

	tot = sio2+tio2+al2o3+fe2o3+h2o+li2o+na2o+k2o+mgo+cao+feo+bao + mno + nio + p2o5

	chemtable["sio2"]=sio2/tot
	chemtable["tio2"]=tio2/tot
	chemtable["al2o3"]=al2o3/tot
	chemtable["fe2o3"]=fe2o3/tot
	chemtable["h2o"]=h2o/tot
	chemtable["li2o"]=li2o/tot
	chemtable["na2o"]=na2o/tot
	chemtable["k2o"]=k2o/tot
	chemtable["mgo"]=mgo/tot
	chemtable["cao"]=cao/tot
	chemtable["feo"]=feo/tot
	chemtable["bao"]=bao/tot

	chemtable["mno"]=mno/tot
	chemtable["nio"]=nio/tot
	chemtable["p2o5"]=p2o5/tot

	return chemtable


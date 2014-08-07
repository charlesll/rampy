# -*- coding: utf-8 -*-
"""
Created on Tue May 27 14:59:13 2014

@author: charleslelosq
Python script for removing interference fringe from FTIR or Raman spectra
It applies a bandstop butterfilter (from the scipy toolbox) implemented in the spectratools module
"""


import sys
sys.path.append("/Users/charleslelosq/anaconda/lib/python2.7/lib-charles")

import numpy as np
import scipy
import matplotlib
from pylab import *
from spectratools import spectrafilter
from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfile
#
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
rawspectre = np.genfromtxt(filename)#,skip_header=20, skip_footer=43) # Skipping lines from JASCO files

cutfq = np.array([0.021,0.025]) # PLEASE SELECT THE FREQUENCY DOMAIN TO CUT HERE (=1/P, with P the period of oscillation in cm-1)
treatsp = spectrafilter(rawspectre,'bandstop',cutfq,1,np.array([1]))

figure(1,[10,7])
plot(treatsp[:,0],treatsp[:,1],'r-')
plot(rawspectre[:,0],rawspectre[:,1],'k-')

# data are out in this directory
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
savefilename = asksaveasfile() # show an "Open" dialog box and return the path to the selected file
np.savetxt(savefilename,treatsp)
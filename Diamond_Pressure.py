# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:34:21 2014

@author: charleslelosq
Carnegie Institution for Science

This script calculates pressure (MPa) in DAC experiments from Raman shift of 13C diamonds
Put it anywhere, you need however the lib-charles library as well as numpy, scipy, matplotlib, and Tkinter
that usually come with any python distribution
"""
import sys
sys.path.append("/Users/charleslelosq/anaconda/lib/python2.7/lib-charles")

import numpy
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate
from scipy.optimize import curve_fit

from spectratools import linbaseline
from spectratools import pseudovoigt
from spectratools import spectrafilter

from Tkinter import *
from tkFileDialog import askopenfilename
import os


##############FOR THE 25° DIAMOND #################
tkMessageBox.showinfo(
            "Open file",
            "Please open the 25C Diamond spectrum")
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
print(filename)

rawspectre = np.genfromtxt(filename,skip_header=20, skip_footer=43) # We skip lines for the JASCO spectrometer, to be changed in some other case

figure(1)
plot(rawspectre[:,0],rawspectre[:,1],'k-')

# Individualize the 13C peak
pic13C = rawspectre[np.where((rawspectre[:,0]> 1220) & (rawspectre[:,0]<1305))]
##### SETUP THE FOLLOWING TO CUT HIGH FREQUENCY NOISE ON THE 13C PEAK #### COMMENT IF NOT NEEDED    
#cutfq = np.array([0.1]) #
#pic13C = spectrafilter(pic13C,'low',cutfq,1,np.array([1]))
#plot(pic13C[:,0],pic13C[:,1],'g-')

# SUbtract a linear baseline below it
bir13C = np.array([(1250,1270),(1298,1302)])
pic13C, baseline13C, coeffs13C = linbaseline(pic13C,bir13C,'linear',1)

#Fitting the peak with the pseudovoigt
para13C, cov13C = curve_fit(pseudovoigt, pic13C[:,0], pic13C[:,1],[1500,1283,2,0.5])

plot(pic13C[:,0],pic13C[:,1],'k-')
plot(pic13C[:,0],pseudovoigt(pic13C[:,0],para13C[0],para13C[1],para13C[2],para13C[3]),'r-')

# Individualize the Ne peak
picNe = rawspectre[np.where((rawspectre[:,0]> 1700) & (rawspectre[:,0]<1740))]

# SUbtract a linear baseline below it
birNe = np.array([(1700,1710),(1725,1730)])
picNe, baselineNe, coeffsNe = linbaseline(picNe,birNe,'linear',1)

#Fitting the peak with the pseudovoigt
paraNe, covNe = curve_fit(pseudovoigt, picNe[:,0], picNe[:,1],[4000,1718,2,1])

plot(picNe[:,0],picNe[:,1],'k-')
plot(picNe[:,0],pseudovoigt(picNe[:,0],paraNe[0],paraNe[1],paraNe[2],paraNe[3]),'r-')

# Ne line is usually at 584.72 nm so 1710.22 cm-1
correction = paraNe[1]-1710.22 # Exact shift of spectrometer
Freq13C_RT = para13C[1]-correction # Corrected frequency of the diamond at RT

##############FOR THE HT° DIAMOND #################
tkMessageBox.showinfo(
            "Open file",
            "Please open the HT Diamond spectrum")
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
print(filename)

rawspectre = np.genfromtxt(filename,skip_header=20, skip_footer=43) # Same as above
figure(2)
plot(rawspectre[:,0],rawspectre[:,1],'k-')

# Individualize the 13C peak
pic13C = rawspectre[np.where((rawspectre[:,0]> 1220) & (rawspectre[:,0]<1305))]
##### SETUP THE FOLLOWING TO CUT HIGH FREQUENCY NOISE ON THE 13C PEAK #### COMMENT IF NOT NEEDED    
#cutfq = np.array([0.1]) #
#pic13C = spectrafilter(pic13C,'low',cutfq,1,np.array([1]))
#plot(pic13C[:,0],pic13C[:,1],'g-')

# SUbtract a linear baseline below it
bir13C = np.array([(1240,1245),(1298,1302)])
pic13C, baseline13C, coeffs13C = linbaseline(pic13C,bir13C,'linear',1)

#Fitting the peak with the pseudovoigt
para13C_HT, cov_HT = curve_fit(pseudovoigt, pic13C[:,0], pic13C[:,1],[1000,1273,5,0.1])

plot(pic13C[:,0],pic13C[:,1],'k-')
plot(pic13C[:,0],pseudovoigt(pic13C[:,0],para13C_HT[0],para13C_HT[1],para13C_HT[2],para13C_HT[3]),'r-')

# Individualize the Ne peak
picNe = rawspectre[np.where((rawspectre[:,0]> 1700) & (rawspectre[:,0]<1740))]

# SUbtract a linear baseline below it
birNe = np.array([(1700,1710),(1725,1730)])
picNe, baselineNe, coeffsNe = linbaseline(picNe,birNe,'linear',1)

#Fitting the peak with the pseudovoigt
paraNe_HT, covNe_HT = curve_fit(pseudovoigt, picNe[:,0], picNe[:,1],[4000,1719,1,0.5])

plot(picNe[:,0],picNe[:,1],'k-')
plot(picNe[:,0],pseudovoigt(picNe[:,0],paraNe_HT[0],paraNe_HT[1],paraNe_HT[2],paraNe_HT[3]),'r-')

# Ne line is usually at 584.72 nm so 1710.22 cm-1
correction = paraNe_HT[1]-1710.22 # Exact shift of spectrometer
Freq13C_HT = para13C_HT[1]-correction

################ FLLOWING LINES CREATE A GUI FOR ASKING FOR TEMPERATURE #################

class Dialog(Toplevel):

    def __init__(self, parent, title = None):

        Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.pack(padx=5, pady=5)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden

        pass

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = Frame(self)

        w = Button(box, text="OK", width=10, command=self.ok, default=ACTIVE)
        w.pack(side=LEFT, padx=5, pady=5)
        w = Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):

        return 1 # override

    def apply(self):

        pass # override
        
class MyDialog(Dialog):

    def body(self, master):

        Label(master, text="Temperature:").grid(row=0)
        #Label(master, text="Second:").grid(row=1)

        self.e1 = Entry(master)
        #self.e2 = Entry(master)

        self.e1.grid(row=0, column=1)
        #self.e2.grid(row=1, column=1)
        return self.e1 # initial focus

    def apply(self):
        T = int(self.e1.get())
        self.result = T
        #second = int(self.e2.get())
        #print first, second # or something
        
root=Tk()
d=MyDialog(root)
T = d.result
root.destroy()

############ WE HAVE TEMPERATURE? THEN WE CALCULATE PRESSURE WITH EQ. FROM MYSEN YAMASHITA 2010

P = (Freq13C_HT-Freq13C_RT+1.065*10**-2*T + 1.769*10**-5*(T**2)) / 0.002707 #MPa




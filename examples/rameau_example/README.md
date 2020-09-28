# RAMEAU

Rameau is a set of functions and a class from the rampy library that provides a way to calculate the water content of a glass from its Raman spectrum.

At this stage, it provides access to external Raman calibrations, from the Le Losq et al. (2012) and Di Genova et al. (2017) articles.

# INSTALLATION

Works best on Mac OS or Linux 

Problems with gcvspline in Windows are frequent, due to the requirement of a FORTRAN compiler.

Install a Python stack (see Anaconda Python for instance) with scipy, numpy, matplotlib, pandas

Install rampy and gcvspline using the command in a terminal (in your working virtual environment) :

```
$ pip install rampy
$ pip install gcvspline
```

# Using rampy

- the notebook provides some interactive example with figures showcasing the different calibrations.

- however, the easiest usage is to work with the rameau.py script provided in this folder. This allows calling rameau directly from the terminal, providing a simple set of arguments. In this case, **read the README_RAMEAU.txt file.**

# References

C. Le Losq, D. R. Neuville, R. Moretti, J. Roux, Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist. 97, 779–790 (2012).

DG 2017 D. Di Genova et al., Effect of iron and nanolites on Raman spectra of volcanic glasses: A reassessment of existing strategies to estimate the water content. Chemical Geology. 475, 76–86 (2017).
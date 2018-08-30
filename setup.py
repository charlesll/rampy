#from numpy.distutils.core import setup
#from numpy.distutils.extension import Extension
from setuptools import setup, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='rampy',
      version='0.4.1',
      description='A Python module containing functions to treat spectroscopic (XANES, Raman, IR...) data',
      url='https://github.com/charlesll/rampy',
      author='Charles Le Losq',
      author_email='charles.lelosq@anu.edu.au',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='GNU-GPLv2',
      packages=['rampy', 'rampy.tests'],
      install_requires=['numpy>=1.12','scipy','sklearn','pandas','cvxpy>=1.0'],
      zip_safe=False)

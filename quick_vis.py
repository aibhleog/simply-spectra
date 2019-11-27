'''
Code used to plot up a quick extraction of raw data around the area of interest.
Requires the id name of the galaxy, the wavelength range of interest, and the 
center row to extract.

From there, the code will run a gaussian extraction with an aperture size of 7
and a default FWHM reminiscent of "decent" MOSFIRE seeing (~0.8"), then plotting
the signal and error from the multi-extension reduced data.
'''

import numpy as np
import matplotlib.pyplot as plt

def showme(id_gal,ycen, ...

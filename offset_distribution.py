'''
Measuring the distribution of offsets for all continuum objects.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.stats import mad_std
from scipy.optimize import curve_fit
import pandas as pd
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


# -- prompting user for mask name & date of observation -- #
# -- Mask name
mask = input('Type mask name & then press Enter:  ')
# making sure something was entered
assert len(mask) > 0, f"Need to enter a mask name. "\
	f"Currently, the mask is: \n \t\t  mask = {mask}"
print(f'Mask entered: "{mask}"',end='\n\n')

# -- Date of Observation
date = input('Type date & press Enter:  ')
assert len(date) > 0, f"Need to enter the date the data were taken. "\
	f"Currently, the date is: \n \t\t  date = {date}"
print(f'Date entered: "{date}"',end='\n\n')
# ------------------------------------------------------------- #

# reading in offsets measured from "continuum_offsets.py"
offsets = pd.read_csv(f'spectral-secrets/continuum_offsets_{mask}.txt',sep='\s+')
print(offsets,end='\n\n')

# making figure
plt.figure(figsize=(9,5))
hist = plt.hist(offsets.offset,bins=5,edgecolor='k')
plt.xlabel('continuum offset from DRP location (spatially) ')
plt.title(f'Offset Distribution for {mask} / {date}')
plt.gca().set_yticks([0,1,2,3,4,]) # because we only care about integer numbers here

plt.tight_layout()
plt.savefig(f'spectral-secrets/vis/{date}/offset_distribution_{mask}.pdf')
plt.close('all')

'''
Computing the distribution of offsets for all continuum objects.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.stats import mad_std
from scipy.optimize import curve_fit
from astropy.modeling.functional_models import Moffat1D
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

# for fitting the continuum shape
def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian function, to be used in the quick optimized extraction
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset


# defining paths -- this will be deleted once database of sources
# has been created in the private repository
home = '/run/media/aibhleog/eve_lily/observing-keck/'
path2d = 'rebecca_reductions/ContaminationSubtracted/'


# -- prompting user for galaxy ID name & date of observation -- #
# -- ID of galaxy
id_gal = input('Type ID name & press Enter:  ')
# making sure something was entered
assert len(id_gal) > 0, f"Need to enter an ID name. "\
	f"Currently, the ID name is: \n \t\t  id_gal = {id_gal}"
print(f'ID entered: "{id_gal}"',end='\n\n')

# -- Date of Observation
date = input('Type date & press Enter:  ')
assert len(date) > 0, f"Need to enter the date the data were taken. "\
	f"Currently, the date is: \n \t\t  date = {date}"
print(f'Date entered: "{date}"',end='\n\n')
# ------------------------------------------------------------- #


# reading in data
print(f'Reading in data for {id_gal}')
path = home + date + '/' + path2d
header = fits.getheader(path + id_gal + '_eps.fits')
signal = fits.getdata(path + id_gal + '_eps.fits')
error = fits.getdata(path + id_gal + '_sig.fits')
print(f'Dimensions: \t signal spectrum {signal.shape}\n' +
		f'\t\t error spectrum {error.shape}',end='\n\n')

# y center of the object (via the DRP)
ycen = abs(header['CRVAL2']) # pixels

# collapsing spectrum over the spectral direction
profile = np.nansum(signal,axis=1) # in case there are NaNs
med,mad = np.nanmedian(profile), mad_std(profile,ignore_nan=True)

# fitting the profile to find the continuum, but only focusing on the
# positive flux (will hopefully help eliminate noise counting as signal)
positive_profile = profile.copy()
positive_profile[profile < 0] = 0

# initial guesses for gaussian()
p0 = [ycen, np.nanmax(positive_profile[int(ycen)-10:int(ycen)+11]), 2, 2] 
popt,pcov = curve_fit(gaussian, range(len(profile)), positive_profile, p0=p0)

# plotting profile
plt.figure(figsize=(9,6))
plt.plot(profile)
plt.axvline(ycen,color='g',lw=3,label='sugg. center (via DRP)')
plt.plot(gaussian(range(len(profile)),*popt),label='Gaussian fit')

plt.legend()
#plt.ylim(med-mad*3, med+mad*3)
plt.xlabel('spatial axis [pixels]')
plt.ylabel('flux')

plt.tight_layout()
if os.path.exists(f'spectral-secrets/vis/{date}') == False:
	print(f'Creating spectral-secrets/vis/{date} directory.')
	os.system(f'mkdir spectral-secrets/vis/{date}')
plt.savefig(f'spectral-secrets/vis/{date}/profile_{id_gal}.pdf')
print('Figure saved.')
	
#plt.show()
plt.close('all')




















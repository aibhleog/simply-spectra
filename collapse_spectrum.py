'''
Takes a given spectrum and collapses it spatially to inspect for any
possible contamination (via continuum).

Currently just using this code for one galaxy in particular.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

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

# pixel scale
pixelscale = header['PSCALE'] # arcseconds/pixel

# collapsing spectrum over the spectral direction
profile = np.nansum(signal,axis=1) # in case there are NaNs

plt.figure(figsize=(9,6))
plt.plot(profile)

med = np.nanmedian(profile)
plt.ylim(med-med*0.2, med+med*0.3)
plt.xlabel('spatial axis [pixels]')
plt.ylabel('flux')

plt.tight_layout()
if os.path.exists(f'spectral-secrets/vis/{date}') == False:
	print(f'Creating spectral-secrets/vis/{date} directory.')
	os.system(f'mkdir spectral-secrets/vis/{date}')
plt.savefig(f'spectral-secrets/vis/{date}/profile_{id_gal}.pdf')
print('Figure saved.')
	
plt.show()
plt.close('all')











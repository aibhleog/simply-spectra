'''
Computing the distribution of offsets for all continuum objects.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.stats import mad_std
from scipy.optimize import curve_fit
from astropy.modeling.functional_models import Moffat1D
import pandas as pd
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


# -- prompting user for mask name, observation cycle, & date of observation -- #
# -- Mask name & cycle
mask_cycle = input('Type mask name & cycle with no spaces (ex. special_mask,2020B)\n & then press Enter:  ')
# making sure something was entered
assert len(mask_cycle) > 0, f"Need to enter an mask & cycle. "\
	f"Currently, the mask_cycle is: \n \t\t  mask_cycle = {mask_cycle}"
mask,cycle = mask_cycle.split(',')
print(f'\nMask entered: "{mask}"')
print(f'Cycle entered: "{cycle}"',end='\n\n')

# -- Date of Observation
date = input('Type date & press Enter:  ')
assert len(date) > 0, f"Need to enter the date the data were taken. "\
	f"Currently, the date is: \n \t\t  date = {date}"
print(f'Date entered: "{date}"',end='\n\n')
# ------------------------------------------------------------- #

targets = pd.read_csv(f'spectral-secrets/masks/{cycle}/{mask}.coords',sep='\s+',\
	names=['id','pri','mag','hh','hm','hs','dd','dm','ds','j2','j22','z','zz'])
targets.drop(index=targets.index.values[len(targets)-4:],inplace=True) # alignment stars

offsets = pd.DataFrame({'id_gal':[],'offset':[],'ycen':[],'mean':[],'A':[],'sig':[],'C':[]})
for id_gal in targets.id.values:
	# reading in data
	print(f'Reading in data for {id_gal}')
	path = home + date + '/' + path2d
	header = fits.getheader(path + id_gal + '_eps.fits')
	signal = fits.getdata(path + id_gal + '_eps.fits')
	error = fits.getdata(path + id_gal + '_sig.fits')
	print(f'Dimensions: \t signal spectrum {signal.shape}\n' +
			f'\t\t error spectrum {error.shape}',end='\n\n')

	# y center of the object (via the DRP)
	ycen = abs(header['CRVAL2']) - 7 # pixels, -7 because clipping on line 70
	bad = False # start out assuming fit will work

	# collapsing spectrum over the spectral direction
	profile = np.nansum(signal,axis=1) # in case there are NaNs
	profile = profile[7:-7].copy() # clipping off the dumb edges
	med,mad = np.nanmedian(profile), mad_std(profile,ignore_nan=True)

	# fitting the profile to find the continuum, but only focusing on the
	# positive flux (will hopefully help eliminate noise counting as signal)
	positive_profile = profile.copy()
	positive_profile[profile < 0] = 0

	# initial guesses for gaussian()
	p0 = [ycen, np.nanmax(positive_profile[int(ycen)-7:int(ycen)+7]), 3, 2] 
	try: 
		popt,pcov = curve_fit(gaussian, range(len(profile)), positive_profile, p0=p0,
			bounds=([ycen-7,1,2,-np.inf],[ycen+7,np.inf,5,med+mad*4]))
		print(f'Cen: {popt[0]}, A: {popt[1]}, sig: {popt[2]}, offset: {popt[3]}')
	except ValueError: bad = True; pass
	if any(0 in sublist for sublist in pcov) == True: print('Bad fit'); bad = True
	elif round(popt[2],5) == 2.: print('Bad fit'); bad = True

	# plotting profile
	plt.figure(figsize=(9,6))
	plt.plot(profile)
	plt.axvline(ycen,color='g',lw=3,label='sugg. center (via DRP)')
	if bad == False:
		plt.plot(gaussian(range(len(profile)),*popt),label='Gaussian fit')
		plt.text(0.03,0.05,f'Offset from sugg. center: {round(popt[0]-ycen,2)} pixels',\
			transform=plt.gca().transAxes)
	
	plt.legend()
	if id_gal[:4] != 'star':
		plt.ylim(med-mad*3, med+mad*5)
	plt.xlabel('spatial axis [pixels]')
	plt.ylabel('flux')
	plt.title(f'identifying continuum for {id_gal}')

	plt.tight_layout()
	if os.path.exists(f'spectral-secrets/vis/{date}') == False:
		print(f'Creating spectral-secrets/vis/{date} directory.')
		os.system(f'mkdir spectral-secrets/vis/{date}')
	plt.savefig(f'spectral-secrets/vis/{date}/profile_{id_gal}.pdf')
	print('Figure saved.',end='\n\n')
	
	plt.show()
	plt.close('all')

	# adding offset measurement to dataframe
	if bad == False:
		saveit = input('Accept fit? (y/n)  ')
		if saveit == 'y':	
			print('Okay, saving fit!')
			# making filler dataframe
			# note that I'm adding back in the 7 pixels that I had clipped out for 'profile'
			filler_df = pd.DataFrame({'id_gal':[id_gal],'offset':[round(popt[0]-ycen,3)],\
				'ycen':[ycen+7],'mean':[popt[0]+7],'A':[popt[1]],'sig':[popt[2]],'C':[popt[3]]})
			offsets = offsets.append(filler_df,ignore_index=True).copy()
	print('-------------',end='\n\n')

print('Offsets:\n',offsets,end='\n\n')
offsets.to_csv(f'spectral-secrets/continuum_offsets_{mask}.txt',sep='\t',index=False)
print('Offsets saved')
















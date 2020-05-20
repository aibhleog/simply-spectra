'''
Code used to plot up a quick extraction of raw data around the area of interest.
Requires the id name of the galaxy, date of observation, the wavelength range of
interest, and the center row to extract.
--> NOTE that this is for the reductions from RLL, therefore the files
    we care about are: [id]_eps.fits, [id]_sig.fits, [id]_snr.fits

From there, the code will run a gaussian extraction with an aperture size of 7
and a default FWHM reminiscent of "decent" MOSFIRE seeing (~0.8"), then plotting
the signal and error from the multi-extension reduced data.

To be implemented:
	Note that the the code uses the ID name of the galaxy to read in the file
	from the proper path.
'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

# for the optimized extraction shape
def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian function, to be used in the quick optimized extraction
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset
	
def showme(id_gal,date,ycen,xlims,aper=7,savefig=False,see=True):
	'''
	The main tool in this script, this function makes a quick image of a given reduced FITS
	image. The function takes the following input:
	
	INPUTS ---- id_gal:   	str, name of the reduced 2D MOSFIRE file to be read in
		        date:		str, date the data were taken
		        ycen:		int, row of the center of the emission (or thing of interest)
				xlims:		list, two integers describing the wavelength range of interest
				aper:		int, number of rows to extract the 1D spectrum over
				savefig:	bool, save the figure in the private repo?
				see:		bool, allow the plt.show() line?
				
	RETURNS ---	quick look figure			
	'''

	# making sure the aperture is an odd number
	assert aper%2 == 1, "Aperture size needs to be an odd number of pixels. "\
		f"Currently, the aperture size is: \n \t\t  aper = {aper} pixels (default is 7 pixels)."
	# making sure that xlims is a list with 2 items (for the plt.xlim)
	assert len(xlims) == 2, "The xlims kwarg needs to be a list with 2 values. "\
		f"Currently, xlims is: \n \t\t  xlims = {xlims}."

	# defining paths -- this will be deleted once database of sources
	# has been created in the private repository
	home = '/run/media/aibhleog/eve_lily/observing-keck/'
	path2d = 'rebecca_reductions/ContaminationSubtracted/'

	# reading in data
	print(f'\nReading in data for {id_gal}')
	path = home + date + '/' + path2d
	header = fits.getheader(path + id_gal + '_eps.fits')
	signal = fits.getdata(path + id_gal + '_eps.fits')
	error = fits.getdata(path + id_gal + '_sig.fits')
	print(f'Dimensions: \t signal spectrum {signal.shape}\n' +
			f'\t\t error spectrum {error.shape}')
	
	# defining the wavelength array
	stepsize = header['CDELT1']
	start_wave = header['CRVAL1']
	end_wave = start_wave + (len(signal[0])*stepsize)
	wave = np.arange(start_wave,end_wave,stepsize) # wavelength array [Angstroms]

	# defining optimized extraction gaussian
	fwhm = 0.8 / 0.18 # arcsec / [arcsec/pixel]
	gauss = gaussian(np.arange(aper),mean=3.,A=1.,sig=fwhm/2.35,offset=0.)
	gauss /= sum(gauss) # to make it sum to 1 to use as weights
	gauss_2D = np.zeros((len(gauss),len(wave))) # making 2D array of weights
	for i in range(aper):
		gauss_2D[i] = gauss[i]

	# optimally-extracting 1D spectra
	half = int(aper/2) # to make the cut out of the 2D image
	spec = np.nansum(signal[ycen-half:ycen+half+1].copy()*gauss_2D,axis=0)
	err = np.nansum(error[ycen-half:ycen+half+1].copy()*gauss_2D,axis=0)

	# plotting quick vis
	plt.figure(figsize=(8,5))
	plt.step(wave/1e4,spec,where='mid',label=f'aper: {aper} pixels')
	plt.fill_between(wave/1e4,err,-1*err,alpha=0.2,color='C0')
	plt.xlim(xlims[0]/1e4,xlims[1]/1e4)

	here = np.where((wave > xlims[0])&(wave < xlims[1]))[0]
	yvals = spec[here]
	plt.ylim(min(yvals)-min(yvals)*0.2,max(yvals)+max(yvals)*0.5)

	plt.legend()
	plt.ylabel('flux')
	plt.xlabel('wavelength [microns]')
	plt.title(f'Quick Vis of {date} / {id_gal}')

	# saving figure -- first checking if directory in the 
	# 'spectral-secrets' private repo exists; if not, making one
	if savefig == True:
		if os.path.exists(f'spectral-secrets/quick_vis/{date}') == False:
			print(f'Creating spectral-secrets/quick_vis/{date} directory.')
			os.system(f'mkdir spectral-secrets/quick_vis/{date}')
		plt.savefig(f'spectral-secrets/quick_vis/{date}/{id_gal}.pdf')
		print('Figure saved.')
	
	if see == True: plt.show()
	plt.close('all')


































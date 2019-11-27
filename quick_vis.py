'''
Code used to plot up a quick extraction of raw data around the area of interest.
Requires the id name of the galaxy, date of observation, the wavelength range of
interest, and the center row to extract.

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

# for the optimized extraction shape
def gaussian(xaxis, mean, A, sig, offset): 
    return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset
    
# quick conversion from wavelength to pixel
def wav2pix(lam,lam_0,step):
	return (lam - lam_0) / step

def pix2wav(pix,lam_0,step):
	return lam_0 + pix * step

def showme(id_gal,date,ycen,xlims,aper=7):
	# defining paths -- this will be deleted once database of sources
	# has been created in the private repository
	home = '/run/media/aibhleog/eve_lily/observing-keck/'
	path2d = 'mosdef_drp/Reduced/2D/'

	# reading in data
	path = home + date + '/' + path2d
	hdul = fits.open(path + id_gal)
	print('Opened {id_gal}, here are the extensions:',end='\n\n')
	print(hdul.info(),end='\n\n')

	print('Extracting extensions 1 & 3 for the signal & error spectra, respectively.',end='\n\n')
	
	# extracting the 2D signal & error spectra
	header = hdul[1].header # header informatin
	signal = hdul[1].data # 2D spectra
	error = hdul[3].data # 2D spectra

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
	spec = np.nansum(signal[ycen-3:ycen+4].copy()*gauss_2D,axis=0)
	err = np.nansum(error[ycen-3:ycen+4].copy()*gauss_2D,axis=0)
	
	plt.figure()
	plt.plot(wave,spec)
	plt.fill_between(wave,err,-1*err,alpha=0.2,color='C0')
	plt.xlim(xlims[0],xlims[1])

	here = np.where((wave > xlims[0])&(wave < xlims[1]))[0]
	yvals = spec[here]
	plt.ylim(min(yvals)-min(yvals)*0.2,max(yvals)+max(yvals)*0.5)


	plt.show()
	#print(list(spec))
	
	hdul.close() # make sure to close it at the end


































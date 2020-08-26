'''
This module allows an input *intrinsic* spectrum and E(B-V) and 
will output the estimated *observed* SED that has been reddened.
'''

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import numpy as np

# -- E(B-V)_star -- #
def ebv_star(ebv):
	return 0.44*ebv
	
# -- effective attenuation curve -- #
def k(lam):
	# must be in microns
	indx = np.arange(len(lam))
	low = indx[lam < 0.63]
	high = indx[lam >= 0.63]
	
	curve_low,curve_high = [],[]
	for i in range(len(low)):
		if lam[low[i]] >= 0.12:
			curve_low.append(2.659*(-2.156 + 1.509*lam[low[i]] \
				 - 0.198*(lam[low[i]])**2 + 0.011*(lam[low[i]])**3) + 4.05)
	for i in range(len(high)):	
		if lam[high[i]] <= 2.2:
			curve_high.append(2.659*(-1.857 + 1.04/lam[high[i]]) + 4.05)
	
	curve = np.concatenate((curve_low,curve_high))
	return curve
	
# -- calzetti law -- #
def calzetti(lam,spec,ebv,which='obs',vfv=True):
	# the which key allows you to choose whether you
	# are wanting the observed or intrinsic spectrum
	indx = np.arange(len(lam))
	below = indx[lam <= 2.2]
	here = indx[np.where(lam[below] >= 0.12)]
	fix_spec = spec[here]
	
	if which == 'obs':
		# when you input the *intrinsic* SED
		sed = fix_spec / np.power(10,0.4*k(lam)*ebv_star(ebv))
	else:
		# when you input the *observed* SED
		sed = fix_spec * np.power(10,0.4*k(lam)*ebv_star(ebv))
	
	if vfv == True:
		nu = 2.998e14 / lam[here]
		sed *= nu
	
	return sed


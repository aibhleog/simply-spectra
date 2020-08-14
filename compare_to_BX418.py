'''
Comparing the spectrum of BX418 (from Erb et al. 2010) to a 
particular galaxy's spectrum.  Will also be measuring the
velocity offsets from Lya, CIII], & other lines.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


# for the emission line profiles
def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian functionc
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset

# velocity offset plot using wavelengths instead of redshift
# where CIII] is the systemic (so that z_sys=0)
def velocity_offset(lam_obs,lam_rest):
    return 2.998e5 * ((lam_obs/lam_rest) - 1) # km/s


# reading in BX418 data
# wave [A], flam [cgs]
bx418 = pd.read_csv('BX418rest_flam.dat',names=['wave','flam'],sep='\s+')

# looking at spectrum
plt.figure(figsize=(11,5))
plt.step(bx418.wave,bx418.flam/1e-16,where='mid')

plt.xlabel('rest wavelength [$\AA$]')
plt.ylabel('flux [10$^{-16}$ erg/s/cm$^2$/$\AA$]')
plt.xlim(1190,bx418.loc[len(bx418)-1,'wave'])
plt.ylim(-0.05,0.5)

plt.tight_layout()
plt.show()
plt.close('all')


print('\nLooking at individual lines in the spectrum...',end='\n\n')
# zooming in on places to measure their redshifts
lines = {'c.iii]':[1906.68,1908.73],'lya':[1215.67],'c.iv':[1548.19,1550.76]}
fits = {}

for l,w in lines.items():
    print(f'Looking at line {l}')
    mean = np.mean(w) # for blended doublets
    sub_spec = bx418.query(f'wave > {mean-mean*0.05} and wave < {mean+mean*0.05}').copy()
    sub_spec.reset_index(inplace=True,drop=True)
    
    # fitting line
    # gaussian(xaxis, mean, A, sig, offset)
    popt,pcov = curve_fit(gaussian, sub_spec.wave, sub_spec.flam, \
        p0=[mean,1e-16,2,0.5e-16])
    
    if l == 'c.iv': # because it'll fit the absorption line first
        # marking the absorption line section and flattening it
        no_abs = sub_spec.copy()
        indx = no_abs.query('1535 < wave < 1551').index.values
        no_abs.loc[indx,'flam'] = np.ones(len(no_abs.loc[indx]))*np.median(no_abs.flam.values)
        popt,pcov = curve_fit(gaussian, no_abs.wave, no_abs.flam, \
            p0=[lines['c.iv'][1],1e-16,2,0.5e-16])
    
    print('-- best fit parameters:',popt)
    fits.update({f'{l}':popt}) # adding fits to new dictionary
    
    # calculating velocity offset from CIII]
    if l == 'c.iii]': ciii = velocity_offset(popt[0],mean) # slight offset, used to correct
    veloff = velocity_offset(popt[0],mean)-ciii # correcting the slight offset for CIII]
    print(f'-- velocity offset from systemic: {round(veloff,3)} km/s',end='\n\n')
    
    
    # zooming in on the line of interest
    plt.figure(figsize=(8,5.5))
    plt.step(sub_spec.wave,sub_spec.flam/1e-16,where='mid',label=l)
    plt.plot(sub_spec.wave,gaussian(sub_spec.wave,*popt)/1e-16,
            label=f'$\Delta$v: {round(veloff,3)}')
    
    for i in range(len(w)): plt.axvline(w[i],ls=':',color='g')
    
    if l == 'c.iii]':         
        plt.axvline(1882.71,ls=':',color='C1',label='si.iii]')
        plt.axvline(1892.03,ls=':',color='C1')
    
    plt.xlabel('rest wavelength [$\AA$]')
    plt.ylabel('flux [10$^{-16}$ erg/s/cm$^2$/$\AA$]')
    plt.xlim(sub_spec.wave.values[0],sub_spec.wave.values[-1])
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    plt.close('all')
    print()
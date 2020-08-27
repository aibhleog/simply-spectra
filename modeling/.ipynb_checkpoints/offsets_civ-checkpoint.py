'''
Investigating the offset of CIV emission in the Cloudy models
as a function of ionization, nebular metallicity, stellar metallicity,
stellar population type, age, etc.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from cloudy_func import * # written by TAH
import warnings
from scipy.optimize import OptimizeWarning

warnings.simplefilter("error", OptimizeWarning) # for when no CIV emission seen


__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


# for the emission line profiles
def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian functionc
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset

# velocity offset plot using wavelengths instead of redshift
# the spectrum is already at systemic (so that z_sys=0)
def velocity_offset(lam_obs,lam_rest):
    return 2.998e5 * ((lam_obs/lam_rest) - 1) # km/s


# model parameters
u = np.arange(-3.5,-1.4,0.2)    # full range, 11 ionization points
zneb = [0.1,0.2,0.3,0.5]        # full range, some have 0.4 as well
mass = 300                      # 300 or 100
stars = 'binary'                # binary or single
age = 7                         # 7, 7.477, or 8

neb = 3 # for the nebular metallicity
civ = 1548.19 # angstroms
for ion in u[:2]:
    # pulling model spectrum
    spec = get_cloudy_spec(f'{stars}_cont_{mass}',mass,age,zneb[neb],ioni=ion)
    spec['wavelength'] *= 1e4
    spec['spectrum'] /= (2.998e18/spec.wavelength.values) # nu*Fnu --> Fnu

    # zooming in around CIV
    spec = spec.query('1490 < wavelength < 1610').copy()
    spec['spectrum'] /= np.median(spec.spectrum.values) # normalizing it
    spec = reorder_spectrum(spec) # Cloudy orders it backwards


    # plotting CIV area of spectrum
    plt.figure(figsize=(9,6))
    plt.plot(spec.wavelength,spec.spectrum)
    text_kwargs = {'transform':plt.gca().transAxes,'fontsize':15}
    plt.text(0.025,0.94,f'logU: {round(ion,1)}',**text_kwargs)

    # fitting the CIV emission
    # gaussian(xaxis, mean, A, sig, offset)
    try:
        popt,pcov = curve_fit(gaussian,spec.wavelength,spec.spectrum,p0=[1548,1,2,1])
        plt.plot(spec.wavelength,gaussian(spec.wavelength,*popt))
        plt.axvline(popt[0])
        
        # calculating offset
        offset = velocity_offset(popt[0],civ)
        plt.text(0.025,0.05,f'offset: {round(offset,2)} km/s',**text_kwargs)
        
    except OptimizeWarning:
        print('\nNo emission detected, finding max value.',end='\n\n')
        zoomin = spec.query('1540 < wavelength < 1565').copy()
        # wavelength of peak value
        peak = zoomin.loc[zoomin['spectrum'].idxmax(),'wavelength']
        plt.axvline(peak,ls=':')
        
        # calculating offset using the max value (will likely stay the same)
        offset = velocity_offset(peak,civ)
        plt.text(0.025,0.05,f'offset: {round(offset,2)} km/s',**text_kwargs)
    
    plt.yscale('log')
    plt.ylim(0.25,6)
    plt.gca().set_yticks([0.3,1.,3,])
    plt.gca().set_yticklabels(['0.3','1.0','3.0',])
    plt.xlabel('rest wavelength [$\AA$]')
    plt.ylabel('normalized flux')

    plt.tight_layout()
    plt.show()
    plt.close()

    
print()
# ---------------------------------------------------------- #
# -- running through all models to build table of offsets -- #
# ---------------------------------------------------------- #
zstellar = [0.1,0.2,0.3,0.5]
zneb = [0.1,0.3,0.5]

offsets = pd.DataFrame({'z':[],'zneb':[],'u':[],'offset':[],'age':[],'mass':[],'stars':[]})
for stars in ['binary','single']:
    print('For stars:',stars)
    
    for mass in [300,100]:
        print('For mass:',mass)
        
        for met in zstellar: # stellar metallicity
            print('For Z_stellar:',met)
            
            for neb in zneb: # nebular metallicity
                # checking for when stellar == nebular when it's 0.3 or 0.5
                if neb == 0.1 and met == 0.3: pass # no need to run this model twice
                elif neb == 0.1 and met == 0.5: pass # no need to run this model twice
                else:
                    # need to check if matches stellar
                    if neb == 0.1: neb = met # fix nebular to stellar metallicity
                    print('For Z_neb:',neb)

                    for ion in u:
                        print('For logU:',round(ion,1),end=',\t')
                        # pulling model spectrum
                        spec = get_cloudy_spec(f'{stars}_cont_{mass}',mass,age,met,zneb=neb,ioni=ion)
                        spec['wavelength'] *= 1e4
                        spec['spectrum'] /= (2.998e18/spec.wavelength.values) # nu*Fnu --> Fnu

                        # zooming in around CIV
                        spec = spec.query('1490 < wavelength < 1610').copy()
                        spec['spectrum'] /= np.median(spec.spectrum.values) # normalizing it
                        spec = reorder_spectrum(spec) # Cloudy orders it backwards

                        # fitting the CIV emission
                        # gaussian(xaxis, mean, A, sig, offset)
                        try:
                            popt,pcov = curve_fit(gaussian,spec.wavelength,spec.spectrum,p0=[1548,1,2,1])

                            # calculating offset
                            offset = velocity_offset(popt[0],civ)
                            print(f'offset: {round(offset,2)} km/s')

                        except OptimizeWarning:
                            print('Bad fit/no emission detected.')
                            offset = np.nan

                        filldf = pd.DataFrame({'z':[met],'zneb':[neb],'u':[ion],'offset':[round(offset,3)],\
                                                  'age':[int(age)],'mass':[int(mass)],'stars':[stars]})
                        offsets = offsets.append(filldf,ignore_index=True)

                print()
            print(end='\n\n')
        print(end='\n\n')
    print(end='\n\n\n')

    
# ---------------------------------------- #
# ---------------------------------------- #
print('Saving table to file...',end='\n\n')

df_dtypes = {'zneb':float,'u':float,'offset':float,'age':int,'mass':int,'stars':str}
offsets = offsets.astype(df_dtypes) # to make sure column dtypes don't change
offsets.to_csv('plots-data/offsets_civ.txt',sep='\t',index=False)

print(offsets.head())
'''
Given a 1D spectrum of a galaxy (including an error spectrum), this script will plot
the S/N distribution.  It can also mark the location of particular places in the
spectrum (and where they fall in the distribution).
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import os

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


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

markit = [1132,1148]

# reading in data
print(f'Reading in data for {id_gal}')
path = 'spectral-secrets/vis/' + date + '/' + '1D/'
spec = pd.read_csv(path + id_gal + '_spectrum.txt',sep='\s+')
print(f'File contents: \n{spec.head()}',end='\n\n')

# plotting distribution
spec['snr'] = spec.flux/spec.ferr
plt.figure(figsize=(8,5))

# adding unit gaussian
n, bins, patches = plt.hist(spec.snr, 40, density=1)
mu, sigma = norm.fit(spec.snr)

y = norm.pdf(bins, mu, sigma)
plt.plot(bins, y, lw=2.5, label='$\mu$: %s \n$\sigma$: %s' %(round(mu,2), round(sigma,2)))

unit_gaussian = norm.pdf(bins, 0, 1)
plt.plot(bins, unit_gaussian, label='$\mu$: 0 \n$\sigma$: 1',lw=2, ls='--')
# marking 3sigma lines
plt.axvspan(-3*sigma,3*sigma,color='C0',alpha=0.2,zorder=0, label='$\leq|3\sigma|$')

# marking the peak S/N for each index
plt.axvline(spec.loc[markit[0],'snr'], ls=':', lw=3,label='CIV')
plt.axvline(spec.loc[markit[1],'snr'], ls=':', lw=3)
# plt.text()

plt.legend()
plt.gca().set_yticklabels([])
plt.xlabel('S/N')

plt.tight_layout()
plt.show()
plt.close()


'''
look = spec.query('snr > 4').index.values

for l in look:
    plt.figure()
    s = spec.query(f"{spec.loc[l-20,'wave']} < wave < {spec.loc[l+20,'wave']}").copy()
    plt.plot(s.wave,s.flux,label='flux')
    plt.plot(s.wave,s.smoothed,color='grey',label='smooth')
    plt.fill_between(s.wave,s.ferr,s.ferr*-1,alpha=0.2,zorder=0)
    
    plt.legend()
    plt.ylim(-6e-18,6e-18)

    plt.tight_layout()
    plt.show()
    plt.close()'''
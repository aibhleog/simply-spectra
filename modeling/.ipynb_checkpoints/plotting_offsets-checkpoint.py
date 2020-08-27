'''
Plotting the offset of CIV emission in the Cloudy models
as a function of ionization, nebular metallicity, stellar metallicity,
stellar population type, age, etc.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cloudy_func import * # written by TAH

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

# reading in table of offset
# the bad fits were marked with NaNs
df = pd.read_csv('plots-data/offsets_civ.txt',sep='\t')

try: print(df.loc[0,'wciii']) # if the column exists, it will work
except KeyError: # hopefully will only have to run this once
    cols = ['z', 'zneb', 'u', 'offset', 'age', 'mass', 'stars']
    df = df[cols].copy() # making sure the columns are how I want them for line 28

    # calculating CIII] EW
    w_ciii = []
    for i in df.index.values:
        met,neb,ion,offset,age,mass,stars = df.loc[i].values
        spec = get_cloudy_spec(f'{stars}_cont_{mass}',mass,age,met,zneb=neb,ioni=ion)
        spec['wavelength'] *= 1e4
        spec['spectrum'] /= (2.998e18/spec.wavelength.values) # nu*Fnu --> Fnu
        spec = reorder_spectrum(spec) # Cloudy orders it backwards

        ew = get_cloudy_ew('ciii',spec,manual=True)
        w_ciii.append(ew[0])

    df['wciii'] = w_ciii
    df.to_csv('plots-data/offsets_civ.txt',sep='\t',index=False)

    
# looking at figure
df.loc[df.offset.values > 800, 'offset'] = np.nan

plt.figure(figsize=(8,5.5))
plt.scatter(df.offset,df.wciii,s=100,edgecolor='k',alpha=0.8)


plt.xlabel('CIV offset from systemic [km/s]')
plt.ylabel('CIII] equivalent width [$\AA$]')

plt.tight_layout()
plt.savefig('../spectral-secrets/vis/ciii_vs_civ_offsets.pdf')
plt.show()
plt.close()



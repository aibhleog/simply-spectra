'''
Looking at the P Cygni profiles of the Cloudy models as a function of ionization.
Can also look at nebular metallicity, stellar population type, IMF upper mass limit.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cloudy_func import * # written by TAH

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


# ----------------------------------- #
# -- looking at the general UV ara -- #
# ----------------------------------- #
u,zneb,age,mass = -1.5,0.2,7,300

fig = plt.figure(figsize=(10,5.5))
bin_spec = get_cloudy_spec('binary_cont_300',mass,age,zneb)
sin_spec = get_cloudy_spec('single_cont_300',mass,age,zneb)
agn_spec = get_cloudy_spec_agn('agn',zneb)
name = ['binary-%sM$_{\odot}$ \nU=%s %sMyr \n%sZ$_{\odot}$'%(mass,u,age,zneb),\
        'single-%sM$_{\odot}$ \nU=%s %sMyr \n%sZ$_{\odot}$'%(mass,u,age,zneb),\
        'agn U=%s \n%sZ$_{\odot}$'%(u,zneb)]

count = 0
for spec in [bin_spec,sin_spec,agn_spec]:
    df = spec.loc[spec.u.values == round(u,1)].copy()        
    df['wavelength'] *= 1e4
    df['spectrum'] /= (2.998e18/df.wavelength.values) # nu*Fnu --> Fnu

    # zooming in around CIV
    df = df.query('1150 < wavelength < 1600').copy()
    df['spectrum'] /= np.median(df.spectrum.values) # normalizing it
    df = reorder_spectrum(df) # Cloudy orders it backwards

    plt.plot(df.wavelength,df.spectrum,label=name[count],zorder=5-count)
    count += 1

for l in [1548.19, 1550.76]:
    plt.axvline(l,ls=':',color='k')

plt.legend(frameon=False,fontsize=15,ncol=2)
plt.xlim(1150,1600)
plt.yscale('log')
plt.ylim(3e-1,110)

plt.ylabel(r'normalized F$_{\nu}$',fontsize=17)
plt.xlabel('rest wavelength [$\AA$]',fontsize=16,labelpad=5)

plt.tight_layout()
plt.show()
plt.close()
# ----------------------------------- #


# ------------------------------------------ #
# -- running through series of the models -- #
# ------------------------------------------ #
u = np.arange(-3.5,-1.4,0.2)    # full range, 11 ionization points
zneb = [0.1,0.2,0.3,0.5]        # full range, some have 0.4 as well
mass = 300                      # 300 or 100
stars = 'binary'                # binary or single
age = 7                         # 7, 7.477, or 8

for indx in range(len(zneb)):
    spec = get_cloudy_spec(f'{stars}_cont_{mass}',mass,age,zneb[indx])

    fig = plt.figure(figsize=(10,5.5))
    cmap = plt.get_cmap('viridis')
    colors = [cmap(i) for i in np.linspace(1,0,len(u))]

    for ion in u:
        i = u.tolist().index(ion)
        df = spec.loc[spec.u.values == round(ion,1)].copy()        
        df['wavelength'] *= 1e4
        df['spectrum'] /= (2.998e18/df.wavelength.values) # nu*Fnu --> Fnu

        # zooming in around CIV
        df = df.query('1490 < wavelength < 1610').copy()
        df['spectrum'] /= np.median(df.spectrum.values) # normalizing it
        df = reorder_spectrum(df) # Cloudy orders it backwards

        plt.plot(df.wavelength,df.spectrum,color=colors[i],label=round(ion,1))

    text_kwargs = {'transform':plt.gca().transAxes,'fontsize':16}
    plt.text(0.025,0.85,'P Cygni profile \naround CIV',**text_kwargs)
    plt.text(0.5,0.8,'log(U)',rotation=90,**text_kwargs)
    plt.text(0.025,0.05,'%s Myr BPASS %s stars, IMF $\leq$ %s M$\odot$, %s$ Z\odot$'\
             %(age,stars,mass,zneb[indx]),**text_kwargs)

    #for j in [1548.19, 1550.76]: plt.axvline(j,ls=':',color='k') # just looking at CIV lines
    plt.legend(ncol=3)
    plt.xlim(1490,1610)
    plt.yscale('log')
    plt.ylim(0.25,6)
    plt.gca().set_yticks([0.3,1.,3,])
    plt.gca().set_yticklabels(['0.3','1.0','3.0',])

    plt.ylabel(r'normalized F$_{\nu}$',fontsize=17)
    plt.xlabel('rest wavelength [$\AA$]',fontsize=16,labelpad=5)

    plt.tight_layout()
    #plt.savefig('plots-data/pcygni-profile_%s_age%sz%s_%s.pdf'%(stars,age,zneb[indx],mass))
    plt.show()
    plt.close()
# ------------------------------------------ #
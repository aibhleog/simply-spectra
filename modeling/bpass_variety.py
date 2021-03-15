'''
Looking into the line strengths of lines like CIV from the Cloudy models
to see how they compare with the output 1D spectra.

For example, in the line strength tables I generate, do they include 
stellar absorption in the line fluxes they list?  Like, the CIV values,
are they the absolute values or have they had the absorption subtracted out?
'''

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import smoothing as sm # written by TAH


# equation to get column for age in the BPASS datasets
def get_n(logage):
    return int((logage - 6)/0.1 + 2)


# let's start with something basic
bZs = ['001','002','004','006','010'] # 10% solar, Z=0.002
Zs = [0.05,0.1,0.2,0.3,0.5] # 10% solar
age = 6.5 # log(age) 
m_max = 300 # Msol

# reading in BPASS model
bhome = '/home/aibhleog/Desktop/catalogs/bpassv2/'
    
plt.figure(figsize=(9,5))
colors = [plt.get_cmap('Blues')(x) for x in np.linspace(0.4,1,len(Zs))]

# reading in BPASS model
bhome = '/home/aibhleog/Desktop/catalogs/bpassv2/'

for i in range(len(Zs)):
    bpass = pd.read_csv(bhome + \
                    f'bpassv2_imf135_{m_max}/OUTPUT_CONT/spectra-bin.z{bZs[i]}.dat',\
             usecols=[0,get_n(age)],memory_map=True,\
             delimiter='\s+',names=['wavelength','spectrum'])

    # zoom spec to be around CIV lines
    bpass = bpass.query('1500 < wavelength < 1600').copy()
    bpass['spectrum'] /= bpass['spectrum'].mean()

    plt.plot(bpass.wavelength,bpass.spectrum,label=f'{Zs[i]}',lw=1.5+0.15*i,color=colors[i])

plt.text(0.8,0.1,f'log(age)={age}',fontsize=15,transform=plt.gca().transAxes)

plt.legend(loc=3)
plt.xlabel('wavelength [$\AA$]')
plt.ylabel('normalized')

plt.tight_layout()
plt.show()
plt.close()



# the WR carbon-dominated contribution is most prominent at higher Zs
# the center of this "red bump" is at CIV 5808A

# running through ages
for age in np.arange(6,8.1,0.2):
    plt.figure(figsize=(9,8.5))
    gs = gridspec.GridSpec(2,1,height_ratios=[1,1])#,hspace=0.3)
    ax1 = plt.subplot(gs[0]) # CIV1548
    ax2 = plt.subplot(gs[1]) # CIV5808

    colors = [plt.get_cmap('Blues')(x) for x in np.linspace(0.4,1,len(Zs))]

    for i in range(len(Zs)):
        bpass = pd.read_csv(bhome + \
                        f'bpassv2_imf135_{m_max}/OUTPUT_CONT/spectra-bin.z{bZs[i]}.dat',\
                 usecols=[0,get_n(age)],memory_map=True,\
                 delimiter='\s+',names=['wavelength','spectrum'])

        # zoom spec to be around CIV1548 lines
        bpass_zoom = bpass.query('1500 < wavelength < 1600').copy()
        bpass_zoom['spectrum'] /= bpass_zoom['spectrum'].mean()
        ax1.plot(bpass_zoom.wavelength,bpass_zoom.spectrum,label=f'{Zs[i]}'+' $Z_{\odot}$',\
                 lw=1.5+0.15*i,color=colors[i])


        # zoom spec to be around CIV 5808 lines
        bpass_zoom = bpass.query('5700 < wavelength < 5900').copy()
        bpass_zoom['spectrum'] /= bpass_zoom['spectrum'].mean()
        ax2.plot(bpass_zoom.wavelength,bpass_zoom.spectrum,label=f'{Zs[i]}'+' $Z_{\odot}$',\
                 lw=1.5+0.15*i,color=colors[i])


    ax1.text(0.04,0.87,'BPASS models, $\leq%s$M$_{\odot}$'%m_max,fontsize=15, transform=ax1.transAxes)
    ax1.text(0.8,0.1,f'log(age)={"{:.1f}".format(age)}',fontsize=15,\
             transform=ax1.transAxes)
    ax2.text(0.8,0.87,f'log(age)={"{:.1f}".format(age)}',fontsize=15,\
             transform=ax2.transAxes)

    ax1.set_xlabel('wavelength [$\AA$]')
    ax1.set_ylabel('normalized')
    ax1.set_ylim(0.27,1.77)

    ax2.legend(loc=3,fontsize=13)
    ax2.set_xlabel('wavelength [$\AA$]')
    ax2.set_ylabel('normalized')
    ax2.set_ylim(0.88,1.1)

    plt.tight_layout()
    plt.savefig(f'plots-data/bpass-WR-carbon/redbump-age{"{:.1f}".format(age)}-{m_max}Msol.pdf')
    plt.show()
    plt.close()
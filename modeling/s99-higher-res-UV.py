'''
Playing around with reading in the S99 models, but the higher-res
modeling, too. 

'''


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cloudy_func import *


# reading in data
shome = '/home/aibhleog/Desktop/catalogs/starburst99/'

colors = [plt.get_cmap('Reds')(x) for x in np.linspace(0.3,1,4)]
ssp = ['Kroupa','IMF2','IMF17'] # single or binary populations
imf_name = ['IMF: Kroupa','IMF: -2.0','IMF: -1.7']
Z = ['001','002','008']  # metallicities used
spec = ssp[0]


# plotting
plt.figure(figsize=(9,5))
ax1 = plt.gca()

# LOWER RESOLUTION MODELS -- the standard spectrum
# -------------------------------------------------------------------------
for m in range(len(Z)):
    cols = ['time','wave','logtot','logstellar','lognebular']
    df = pd.read_csv(shome+'s99_cont_%s_v00_z%s.stb99'%
                     (spec,Z[m]),names=cols,skiprows=7,delimiter='\s+')
    s99 = df.loc[df['time'] == 1e7].copy()
    s99['lum'] = np.power(10,s99.logstellar.values)

    s99 = s99.query('1500 < wave < 1600').copy()
    s99.reset_index(inplace=True,drop=True)
    
    s99['lum'] /= s99.loc[3,'lum'] 
    ax1.plot(s99.wave,s99.lum,color=colors[m],lw=2.5-0.15*m) # ,label=lab
# -------------------------------------------------------------------------
    
    
# HIGHER RES MODELS (uvline, 1200-1800A)
# -------------------------------------------------------------------------
for m in range(len(Z)):
    cols = ['time','wave','loglum','normspec']
    df = pd.read_csv(shome+'cont_kroupa_100/s99_cont_%s_v00_z%s.uvline1'%
                     (spec,Z[m]),names=cols,skiprows=7,delimiter='\s+')
    s99 = df.loc[df['time'] == 1e7].copy()
    s99['lum'] = np.power(10,s99.loglum.values)

    s99 = s99.query('1500 < wave < 1600').copy()
    s99.reset_index(inplace=True,drop=True)

    s99['lum'] /= s99.loc[25,'lum'] # index matches about where low-R normalized
    ax1.plot(s99.wave,s99.lum,color=colors[m],lw=2.5-0.25*m)

# -------------------------------------------------------------------------
    
    
ax1.set_ylabel('normalized [erg/s/$\AA$]',fontsize=16)
ax1.set_xlabel('wavelength [$\AA$]',fontsize=15.5)

xlims = ax1.get_xlim() # getting these dimesions now becuase of legend stuff below
ylims = ax1.get_ylim() # cause I wanna change the scaling slightly

ax1.set_yticks([0.5,1,1.5,])
ax1.set_yticklabels(['0.5','1','1.5',])
ax1.set_xlim(xlims[0],xlims[1])
ax1.set_ylim(ylims[0]-0.3,ylims[1]+0.3)


plt.tight_layout()
# plt.savefig('plots-data/civ-stars-nebular.pdf')
# plt.savefig('/home/aibhleog/Documents/papers/mine/figures/civ-stars-nebular.pdf')
plt.show()
plt.close()

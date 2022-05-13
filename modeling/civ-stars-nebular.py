'''
making a three-panel plot zooming in around CIV for 
   1) just stellar populations
   2) cloudy+SF
   3) cloudy+AGN
'''

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import smoothing as sm # written by TAH
from cloudy_func import *


legend_kwargs = {'fontsize':13, 'handlelength':0.8,
                 'handletextpad':0.65, 'labelspacing':0.4}

# equation to get column for age in the BPASS datasets
def get_n(logage):
    return int((logage - 6)/0.1 + 2)


# making figure
plt.figure(figsize=(7,11.5))
gs = gridspec.GridSpec(3,1,height_ratios=[1,1,1],hspace=0)
ax1 = plt.subplot(gs[0]) # stellar pops
ax2 = plt.subplot(gs[1]) # cloudy+SF
ax3 = plt.subplot(gs[2]) # cloudy+AGN




# -------------------------------- #
# --------- STELLAR POPS --------- #
# -------------------------------- #


# reading in BPASS model
bhome = '/home/aibhleog/Desktop/catalogs/bpassv2/'

colors = [plt.get_cmap('Blues')(x) for x in np.linspace(0.3,1,len(Zs))]
bZs = ['001','002','004','006','010'] # Z
Zs = [0.05,0.1,0.4] # % solar
age = 7 # log(age) 
m_max = 300 # Msol


for i in range(len(Zs)):
    lab = 'BPASS v2.0, '+f'{Zs[i]}$\,Z_\odot$'
    bpass = pd.read_csv(bhome + \
                    f'bpassv2_imf135_{m_max}/OUTPUT_CONT/spectra-bin.z{bZs[i]}.dat',\
             usecols=[0,get_n(age)],memory_map=True,\
             delimiter='\s+',names=['wavelength','spectrum'])

    # zoom spec to be around CIV lines
    bpass = bpass.query('1500 < wavelength < 1600').copy()
    bpass.reset_index(inplace=True,drop=True)
    
    bpass['spectrum'] /= bpass.loc[24,'spectrum'] # just to normalize it
    ax1.plot(bpass.wavelength,bpass.spectrum,label=lab,
             lw=2.5-0.15*i,color=colors[i])

    
# adding S99 to the plot
shome = '/home/aibhleog/Desktop/catalogs/starburst99/cont_kroupa_100/'

colors = [plt.get_cmap('Reds')(x) for x in np.linspace(0.3,1,len(Zs))]
ssp = ['Kroupa','IMF2','IMF17'] # IMF slopes available
imf_name = ['IMF: Kroupa','IMF: -2.0','IMF: -1.7']
Z = ['001','002','008']  # metallicities used
spec = ssp[0] # just using Kroupa for now

for m in range(len(Z)):
    cols = ['time','wave','loglum','normspec']
    df = pd.read_csv(shome+'s99_cont_%s_v00_z%s.uvline1'%
                 (spec,Z[m]),names=cols,skiprows=7,delimiter='\s+')
    s99 = df.loc[df['time'] == 1e7].copy()
    s99['lum'] = np.power(10,s99.loglum.values)

	# zoom spec around CIV lines
    s99 = s99.query('1500 < wave < 1600').copy()
    s99.reset_index(inplace=True,drop=True)
    
    s99['lum'] /= s99.loc[25,'lum'] # normalizing
    ax1.plot(s99.wave,s99.lum,color=colors[m],lw=2.5-0.15*m)

    
ax1.text(0.97,0.08,f'log(age)$\,$=$\,${age}',ha='right',fontsize=15,transform=ax1.transAxes)

leg = ax1.legend(loc=3,**legend_kwargs,bbox_to_anchor=(0.01,0))
# set the linewidth of each legend object
for legobj in leg.legendHandles:
    legobj.set_linewidth(7.0)

ax1.set_xticklabels([])
ax1.set_ylabel('normalized [erg/s/$\AA$]',fontsize=16)

xlims = ax1.get_xlim() # getting these dimesions now becuase of legend stuff below
ylims = ax1.get_ylim() # cause I wanna change the scaling slightly

ax1.set_yticks([0.5,1,1.5,])
ax1.set_yticklabels(['0.5','1','1.5',])
ax1.set_xlim(xlims[0],xlims[1])
ax1.set_ylim(ylims[0]-0.3,ylims[1]+0.3)


# creating a twin axis so I can make another legend
ax1_twin = ax1.twiny()
for m in range(len(Z)):
    lab = 'S99, Kroupa'+f', {Zs[m]}$\,Z_\odot$'
    ax1_twin.plot(s99.wave-2000,s99.lum,label=lab,color=colors[m])

leg = ax1_twin.legend(loc=2,**legend_kwargs,bbox_to_anchor=(0.01,1))
# set the linewidth of each legend object
for legobj in leg.legendHandles:
    legobj.set_linewidth(7.0)
    
ax1_twin.set_xticklabels([])
# ax1_twin.set_yticklabels([])
ax1_twin.set_xlim(xlims[0],xlims[1])


# ------------------------------- #
# --------- Cloudy + SF --------- #
# ------------------------------- #

U = -1.5
Z = [0.1,0.2,0.3,0.5]
age = 7

colors = [[plt.get_cmap('Blues')(x) for x in np.linspace(0.2,1,len(Z))],
        [plt.get_cmap('Reds')(x) for x in np.linspace(0.2,1,len(Z))]]



for i in range(len(Z)):
    lab = f'{Z[i]}$\,Z_\odot$'
    spec = get_cloudy_spec(f'r1000',300,age,Z[i],mhome='testing',ioni=U)
    spec = get_fnu(spec) # to go from nu*Fnu to Fnu
    
    spec = spec.query('.1500 < wavelength < .1600').copy()
    spec.reset_index(inplace=True,drop=True)
    
    spec['spectrum'] /= spec.loc[52,'spectrum'] # just to normalize it
    ax2.plot(spec.wavelength.values*1e4,spec.spectrum,
             lw=2.5-0.15*i,color=colors[0][i],label=lab)


ax2.text(0.03,0.91,'Cloudy+BPASS',ha='left',fontsize=14,transform=ax2.transAxes)
ax2.text(0.97,0.08,f'log$\,$U$\,$=$\,-${abs(U)}\nlog(age)$\,$=$\,${age}',
         ha='right',fontsize=15,transform=ax2.transAxes)

ax2.set_yscale('log')
leg = ax2.legend(loc=2,**legend_kwargs,bbox_to_anchor=(0.01,0.92))
# set the linewidth of each legend object
for legobj in leg.legendHandles:
    legobj.set_linewidth(7.0)

ax2.set_xticklabels([])
ax2.set_yticks([0.3,1,3,])
ax2.set_yticklabels(['0.3','1','3',])

ax2.set_xlim(xlims[0],xlims[1])
ylims = ax2.get_ylim()
ax2.set_ylim(ylims[0]-0.05,ylims[1]+0.35)
ax2.set_ylabel('normalized [erg/s/cm$^2$/Hz]',fontsize=16)


# -------------------------------- #
# --------- Cloudy + AGN --------- #
# -------------------------------- #

U = [-1.5,-2.3,-2.5,-2.7,-3.5]
Z = [0.1,0.2,0.3,0.5]
i = 1

colors = [plt.get_cmap('Greys')(x) for x in np.linspace(0.3,1,len(U))]

for u in range(len(U)):
# for i in range(len(Z)):
    # lab = 'Cloudy+AGN, '+f'{Z[i]}$\,Z_\odot$'
    lab = f'log$\,$U =$\,-${abs(U[u])}'
    spec = get_cloudy_spec_agn(f'r1000/manual',Z[i],manual=True,mhome='testing',ioni=[U[u]])
    spec = get_fnu(spec) # to go from nu*Fnu to Fnu
    
    spec = spec.query('.1500 < wavelength < .1600').copy()
    spec.reset_index(inplace=True,drop=True)
    
    spec['spectrum'] /= spec.loc[52,'spectrum'] # just to normalize it
    ax3.plot(spec.wavelength.values*1e4,spec.spectrum,
             lw=2.5-0.15*i,color=colors[u],label=lab)


ax3.text(0.03,0.91,'Cloudy+AGN',ha='left',fontsize=14,transform=ax3.transAxes)
ax3.text(0.97,0.13,f'$Z\,$=$\,${Z[i]} $Z_\odot$',
         ha='right',fontsize=15,transform=ax3.transAxes)

ax3.set_yscale('log')
leg = ax3.legend(loc=2,**legend_kwargs,bbox_to_anchor=(0.01,0.92))
# set the linewidth of each legend object
for legobj in leg.legendHandles:
    legobj.set_linewidth(7.0)

ax3.set_yticks([1,3,10,30,])
ax3.set_yticklabels(['1','3','10','30',])

ylims = ax3.get_ylim()
ax3.set_ylim(ylims[0]-0.05,ylims[1]+0.35)
ax3.set_ylabel('normalized [erg/s/cm$^2$/Hz]',fontsize=16)


ax3.set_xlabel('wavelength [$\AA$]',fontsize=16)
ax3.set_xlim(xlims[0],xlims[1])


plt.tight_layout()
plt.savefig('plots-data/civ-stars-nebular.pdf')
plt.savefig('/home/aibhleog/Documents/papers/mine/figures/civ-stars-nebular.pdf')
# plt.show()
plt.close()




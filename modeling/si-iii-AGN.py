'''
Comparing how SiIII] is affected by AGN

'''


__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import matplotlib.pyplot as plt
from cloudy_func import * 


# agn_spec = get_cloudy_spec_agn('agn',agn_zneb[i])
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0,1,7)]

Z = [0.05,0.1,0.15,0.2,0.3,0.4,0.5]
lines = ['[SiIII]1883', 'SiIII]1892']

plt.figure(figsize=(8.5,9))
gs = gridspec.GridSpec(2,1,hspace=0,height_ratios=[1,1])

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

for zneb in Z:
    tab = get_cloudy_table_agn('agn',zneb)
    
    total =  tab[lines[0]].values + tab[lines[1]].values
    ax.plot(tab.u, total, lw=2.5, color=colors[Z.index(zneb)])
    
    # ratio = tab[lines[0]].values / tab['[CIII]1907'].values
    ratio = tab[lines[0]].values / (tab['[CIII]1907'].values
                                + tab['CIII]1909'].values)
    ax2.plot(tab.u, ratio, lw=2.5, label=f'{zneb} $Z_\odot$',
                color=colors[Z.index(zneb)])
    

ylims = ax2.set_ylim()
xlims = ax2.set_xlim()

# adding shaded region
measured1 = 0.92/4.4 # for SiIII] / tot(CIII]) 
measured2 = 1.08/6.6 # for SiIII] / tot(CIII]) 
# measured1 = 0.92/2.63 # for SiIII] / [CIII] 1907
# measured2 = 1.08/2.63 # for SiIII] / [CIII] 1907

ax2.axhspan(0,measured2,alpha=0.4,facecolor='#962A13',zorder=0,edgecolor='none')
ax2.axhspan(measured2,measured1,alpha=0.6,facecolor='#962A13',zorder=0,edgecolor='none')

plt.fill_between([-3.55,-2.96],[0.15,0.15],[0.17,0.17],\
		facecolor='w',edgecolor='k',zorder=4)
plt.fill_between([-3.55,-2.945],[0.182,0.182],[0.202,0.202],\
		facecolor='w',edgecolor='k',zorder=4)
ax2.text(0.025,0.144,'$<2\sigma$ if [CIII] 1907',fontsize=15,
         transform=ax2.transAxes,zorder=7)
ax2.text(0.025,0.042,'$<2\sigma$ if CIII] 1909',fontsize=15,
         transform=ax2.transAxes,zorder=7)

    
# axes    
ax.set_ylabel('Si$\,$III] $\lambda\lambda$1883,1892')
ax.set_xticklabels([])
ax.set_xticks([-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,])

ax2.legend(fontsize=14,labelspacing=0.2,handletextpad=0.3)#,ncol=2)
# ax2.set_ylabel('Si$\,$III] $\lambda$1883 / [CIII] $\lambda$1907')
ax2.set_ylabel('Si$\,$III] $\lambda$1883 / CIII] $\lambda\lambda$1907,1909')
ax2.set_xlabel('log$\,$U')
ax2.set_ylim(ylims)
ax2.set_xlim(xlims)
ax2.set_xticks([-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,])


plt.tight_layout()
plt.savefig('plots-data/si-iii-ionization-AGN.pdf')
plt.show()
plt.close('all')
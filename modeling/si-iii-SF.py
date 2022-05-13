'''
Comparing how SiIII] is affected by SF

'''


__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import matplotlib.pyplot as plt
from cloudy_func import * 


# agn_spec = get_cloudy_spec_agn('agn',agn_zneb[i])
colors = ['#707B7C','#5DADE2','#1F618D','#6C3483']

Z = [0.1,0.2,0.3,0.5]
lines = ['[SiIII]1883', 'SiIII]1892']
fold = ['binary_cont_300','single_cont_100']
ls = ['-','--']
a = 7 # log(Myr)

plt.figure(figsize=(8.5,9))
gs = gridspec.GridSpec(2,1,hspace=0,height_ratios=[1,1])

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

for s in range(2):
    for zneb in Z:
        tab = get_cloudy_table(fold[s],a,zneb)

        total =  tab[lines[0]].values + tab[lines[1]].values
        ax.plot(tab.u, total, color=colors[Z.index(zneb)],lw=2.5,ls=ls[s])

        # ratio = tab[lines[0]].values / tab['[CIII]1907'].values
        ratio = tab[lines[0]].values / (tab['[CIII]1907'].values
                                    + tab['CIII]1909'].values)
        ax2.plot(tab.u, ratio, color=colors[Z.index(zneb)],lw=2.5,ls=ls[s])


ylims = ax.set_ylim()
xlims = ax.set_xlim()
ylims2 = ax2.set_ylim()
xlims2 = ax2.set_xlim()

# setting up legends
for zneb in Z: ax.plot([0,0],[1,2],label=f'{zneb} $Z_\odot$',
                        color=colors[Z.index(zneb)],lw=2.5)
for s in range(2): ax2.plot([0,0],[1,2],ls=ls[s],color='k',lw=2.5,
                               label=fold[s])


# adding shaded region
measured1 = 0.92/4.4 # for SiIII] / tot(CIII]) 
measured2 = 1.08/6.6 # for SiIII] / tot(CIII]) 
# measured1 = 0.92/2.63 # for SiIII] / [CIII] 1907
# measured2 = 1.08/2.63 # for SiIII] / [CIII] 1907

y = 0.097
plt.fill_between([-3.55,-2.74],[0.15+y,0.15+y],[0.175+y,0.175+y],\
		facecolor='#962A13',edgecolor='k',zorder=4,alpha=0.2)
plt.fill_between([-3.55,-2.7],[0.19+y,0.19+y],[0.216+y,0.216+y],\
		facecolor='#962A13',edgecolor='k',zorder=4,alpha=0.4)
plt.fill_between([-3.55,-2.74],[0.15+y,0.15+y],[0.175+y,0.175+y],\
		facecolor='none',edgecolor='k',zorder=4,lw=1)
plt.fill_between([-3.55,-2.7],[0.19+y,0.19+y],[0.216+y,0.216+y],\
		facecolor='none',edgecolor='k',zorder=4,lw=1)

ax2.text(0.025,0.144,f'$<2\sigma$ = {round(measured1,2)} if [CIII] 1907',fontsize=15,
         transform=ax2.transAxes,zorder=7)
ax2.text(0.025,0.042,f'$<2\sigma$ = {round(measured2,2)} if CIII] 1909',fontsize=15,
         transform=ax2.transAxes,zorder=7)

    
# axes    
ax.set_ylabel('Si$\,$III] $\lambda\lambda$1883,1892')
ax.set_xticklabels([])
ax.set_ylim(ylims)
ax.set_xlim(xlims)
ax.set_xticks([-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,])
ax.legend(fontsize=14,labelspacing=0.2,handlelength=1,handletextpad=0.3)#,ncol=2)

ax2.legend(fontsize=14,labelspacing=0.3,handlelength=1.5)
ax2.set_ylabel('Si$\,$III] $\lambda$1883 / CIII] $\lambda\lambda$1907,1909')
ax2.set_xlabel('log$\,$U')
ax2.set_ylim(ylims2)#[0],0.547)
ax2.set_xlim(xlims2)
ax2.set_xticks([-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,])


plt.tight_layout()
plt.savefig('plots-data/si-iii-ionization-SF.pdf')
plt.show()
plt.close('all')
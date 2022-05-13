'''
Comparing how SiIII] is affected by AGN

'''


__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import matplotlib.pyplot as plt
from cloudy_func import * 



plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1,2,wspace=0,width_ratios=[1,1])

ax = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
       


# AGN
# ------
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0,1,7)]

Z = [0.05,0.1,0.15,0.2,0.3,0.4,0.5]

for zneb in Z:
    tab = get_cloudy_table_agn('agn',zneb)

    si3c3 = get_cloudy_line('si.iii]1',tab) / get_cloudy_line('c.iii]1',tab)
    # si3c3 = get_cloudy_line('si.iii]1',tab) / (get_cloudy_line('c.iii]1',tab)
    #                                 + get_cloudy_line('c.iii]2',tab))
    c43 = (get_cloudy_line('c.iv1',tab) + get_cloudy_line('c.iv2',tab)) / (get_cloudy_line('c.iii]1',tab) + get_cloudy_line('c.iii]2',tab))
    ax.plot(si3c3, c43, lw=2.5, color=colors[Z.index(zneb)], label=f'{zneb} $Z_\odot$')

ax.set_title('AGN')
ax.legend(fontsize=14,labelspacing=0.2,handlelength=1,handletextpad=0.3)#,ncol=2)
ax.set_yscale('log')
ax.set_xlabel('Si$\,$III] $\lambda$1883 / CIII] $\lambda$1907')
ax.set_ylabel('C$\,$IV $\lambda\lambda$1548,1551 / CIII] $\lambda\lambda$1907,1909')
    
ax.set_ylim(0.005,10)
ax.set_xlim(0.25,1)
    
ylims = ax.set_ylim()
xlims = ax.set_xlim()


# Star-forming galaxies
# ------------------
colors = ['#707B7C','#5DADE2','#1F618D','#6C3483']

Z = [0.1,0.2,0.3,0.5]
lines = ['[SiIII]1883', 'SiIII]1892']
fold = ['binary_cont_300','single_cont_100']
ls = ['-','--']
a = 7 # log(Myr)


for s in range(2):
    for zneb in Z:
        tab = get_cloudy_table(fold[s],a,zneb)

        si3c3 = get_cloudy_line('si.iii]1',tab) / get_cloudy_line('c.iii]1',tab)
        # si3c3 = get_cloudy_line('si.iii]1',tab) / (get_cloudy_line('c.iii]1',tab)
        #                                 + get_cloudy_line('c.iii]2',tab))
        c43 = (get_cloudy_line('c.iv1',tab) + get_cloudy_line('c.iv2',tab)) / (get_cloudy_line('c.iii]1',tab) + get_cloudy_line('c.iii]2',tab))
        ax2.plot(si3c3, c43, lw=2.5, color=colors[Z.index(zneb)],ls=ls[s])


ax3 = ax2.twiny()    

# setting up legends
for zneb in Z: ax2.plot([0,0],[1,2],label=f'{zneb} $Z_\odot$',
                        color=colors[Z.index(zneb)],lw=2.5)
for s in range(2): ax3.plot([0,0],[1,2],ls=ls[s],color='k',lw=2.5,
                               label=fold[s])

ax2.legend(fontsize=14,labelspacing=0.2,handlelength=1,handletextpad=0.3)#,ncol=2)
ax3.legend(fontsize=14,labelspacing=0.3,handlelength=1.5,loc=2)
        
ax2.set_title('SFG')
ax2.set_yscale('log')
ax2.set_xlabel('Si$\,$III] $\lambda$1883 / CIII] $\lambda$1907')
ax2.set_yticklabels([])
ax2.set_ylim(ylims)
ax2.set_xlim(xlims)
ax3.set_xlim(xlims)
ax3.set_xticklabels([])

           
           
           
           
plt.tight_layout()
plt.savefig('plots-data/si-iii-c43-AGN-SF.pdf')
plt.show()
plt.close('all')
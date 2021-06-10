'''
This script plots the velocity of Lya from systemic versus Muv.

Original script is from 2017, so I'm slowly converting it to 
better coding, etc.

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.gridspec as gridspec
from abs_uvmag import *

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


hutchison19_zlya = [7.5078,7.5078]
hutchison19_zciii = [7.5032,7.4941]
hutchison19_delv = [84.62+3.783, 405.89+3.783] 
hutchison19_muv = [-21.58,-21.58]


# z>6 LAEs
hiz = pd.read_csv('spectral-secrets/vis/highz-LAEs.txt',delimiter='\s+')


# z~3 LAEs
erb14 = pd.read_csv('/home/aibhleog/Desktop/keckData/information/erb14.txt',delimiter='\s+')
schen13 = pd.read_csv('/home/aibhleog/Desktop/keckData/information/schenker2013.txt',delimiter='\s+',skiprows=3)
sobral18 = pd.read_csv('/home/aibhleog/Desktop/keckData/information/sobral18.txt',delimiter='\s+',skiprows=2)
mclin = pd.read_csv('/home/aibhleog/Desktop/keckData/information/mclinden11-14.txt',delimiter='\s+')

    
# plot Muv vs delv
plt.figure(figsize=(9,7))
gs1 = gridspec.GridSpec(1,2,width_ratios=[4,1])
gs1.update(wspace=0.0,hspace=0.0)

ax = plt.subplot(gs1[0])


# for the love of gods... 2017 Taylor, what were you doing here lol
count,colors = 0, plt.cm.Blues(np.linspace(0,1,len(hiz)-3))
for t in range(len(hiz)-2):
    if hiz.loc[t,'tag'] == 1:
        ax.scatter(hiz.loc[t,'muv'],hiz.loc[t,'delv'],color=colors[t],marker=hiz.loc[t,'marker'],edgecolor='k',s=150,label=hiz.loc[t,'cite'])
        ax.scatter(hiz.loc[t,'muv'],hiz.loc[t,'delv'],color=colors[t],marker=hiz.loc[t,'marker'],edgecolor='k',s=300,zorder=3)
        count+=1
    else: 
        if hiz.loc[t,'tag'] == 0: ax.scatter(hiz.loc[t,'muv'],hiz.loc[t,'delv'],color=colors[count-1],marker=hiz.loc[t,'marker'],edgecolor='k',s=300,zorder=3)
        elif hiz.loc[t,'tag'] == 2: ax.scatter(hiz.loc[t,'muv'],hiz.loc[t,'delv'],color=colors[count],marker=hiz.loc[t,'marker'],edgecolor='k',s=300,zorder=3)
        elif hiz.loc[t,'tag'] == 3: ax.scatter(hiz.loc[t,'muv'],hiz.loc[t,'delv'],color=colors[count+1],marker=hiz.loc[t,'marker'],edgecolor='k',s=300,zorder=3)



# Hutchison et al 2019
ax.scatter(hutchison19_muv[0],hutchison19_delv[0],s=300,marker='*',color='#962A13',lw=1.5,edgecolors='k',label='Hutchison+19')
ax.scatter(hutchison19_muv[0],hutchison19_delv[0],s=900,marker='*',color='#962A13',edgecolors='k',label=None,zorder=5,lw=1.5)
ax.scatter(hutchison19_muv[1],hutchison19_delv[1],s=900,marker='*',facecolor='#F4DCD7',edgecolors='k',lw=1.15,label=None,zorder=3)
ax.plot(hutchison19_muv,hutchison19_delv,ls=':',color='k',zorder=4,lw=2.5)
        
    


# ---------------- literature z=2-3 LAEs -------------------- #
edge,mark = 'k','grey'#'#BF4018'#'#76AD79'
al = 0.6

ax.scatter(erb14['Muv'],erb14['delv'],color=mark,alpha=al,marker='o',s=50,label=None,zorder=1)
im = ax.scatter(erb14['Muv'],erb14['delv'],color='none',edgecolor=edge,alpha=0.73,marker='o',s=50,label=None,zorder=1)
ax.scatter(erb14['Muv'][0]-10,erb14['delv'][0],color=mark,marker='o',alpha=0.68,edgecolor=edge,s=100,label='Erb+14')

ax.scatter(schen13['Muv'][schen13.delv/schen13.derr > 0],schen13['delv'][schen13.delv/schen13.derr > 0],color=mark,alpha=al,marker='s',s=45,label=None,zorder=1)
ax.scatter(schen13['Muv'][schen13.delv/schen13.derr > 0],schen13['delv'][schen13.delv/schen13.derr > 0],color='none',edgecolor=edge,alpha=0.73,marker='s',s=45,label=None,zorder=1)
ax.scatter(schen13['Muv'][0]-10,schen13['delv'][0],color=mark,marker='s',edgecolor=edge,alpha=0.68,s=80,label='Schenker+13')

ax.scatter(sobral18['Muv'],sobral18['delv'],color=mark,alpha=al,marker='D',s=45,label=None,zorder=1)
ax.scatter(sobral18['Muv'],sobral18['delv'],color='none',edgecolor=edge,alpha=0.73,marker='D',s=45,label=None,zorder=1)
ax.scatter(sobral18['Muv'][0]-10,sobral18['delv'][0],color=mark,marker='D',edgecolor=edge,s=80,alpha=0.68,label='Sobral+18')


ax.scatter(mclin.muv,mclin.delv,color=mark,alpha=al,marker='>',s=50,label=None,zorder=1)
ax.scatter(mclin.muv,mclin.delv,color='none',alpha=0.73,marker='>',edgecolor=edge,s=50,label=None,zorder=1)
ax.scatter(np.asarray(mclin.muv)-10,mclin.delv,edgecolor=edge,color=mark,alpha=al,marker='>',s=82,label='McLinden+11,14',zorder=1)
    
from matplotlib.ticker import FormatStrFormatter
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.legend(fontsize=18,ncol=3,handletextpad=0.22,columnspacing=0.6,loc=(-0.07,1.02))
ax.set_xlim(-18,-23.5)
ax.set_ylim(-150,1050)
ax.set_xlabel(r'M$_{\rm UV}$ [AB]',fontsize=19.5)
ax.set_ylabel(r'$\Delta v_{Ly\alpha}$ [km/s]',labelpad=0,fontsize=20)
ax.tick_params(labelsize=18)

# ------------ histogram sideplot -------------------------------- #
ax2 = plt.subplot(gs1[1])
bins=np.histogram(np.hstack((hiz.delv.values[:-2],mclin.delv,erb14.delv,schen13.delv,sobral18.delv)), bins=20)[1]

highz = hiz.delv.values[:-2] # leaving out H19 for now
ax2.hist(highz,bins,orientation='horizontal',alpha=0.52,color=colors[6],label='high-z')
ax2.hist(highz,bins,orientation='horizontal',color='k',histtype='step',label='high-z')

lowz = np.concatenate((erb14.delv,schen13.delv,mclin.delv,sobral18.delv))
ax2.hist(lowz,bins,orientation='horizontal',alpha=0.4,color=mark,label='low-z')
ax2.hist(lowz,bins,orientation='horizontal',color='k',histtype='step',label='low-z')

ax.text(-18.0,1470,'------------ $z\geq6$ LAEs ------------',fontsize=20.5)
ax.text(-22.5,1470,'$z\simeq2\endash3$ LAEs',fontsize=20.5)

ax2.set_ylim(-150,1050)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
# ---------------------------------------------------------------- #


plt.show()
plt.close('all')
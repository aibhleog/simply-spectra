'''
Comparing Spitzer/IRAC [3.6]-[4.5] color to
various line ratios.

NOTE: for now, only my binary_cont_300 and single_cont_100 have the
lines that I want in their tables.
'''

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

import matplotlib.pyplot as plt
from cloudy_func import * 




# the plot

fig, ax1 = plt.subplots(1,1,figsize=(11,5.25))
gs1 = gridspec.GridSpec(1,2)
u = np.arange(-3.5,-1.4,0.2)
colors = ['#707B7C','#5DADE2','#1F618D','#6C3483']
metallicity = [0.1,0.2,0.3,0.5]
zmet = ['-','--',':']; zlw = [1.5,1.5,2.]
age = 7 #[7,7.477,8]
a = age
imf = [300,100]
names = [r'BPASS SPs with binaries, ${\rm M}^{up}_{\rm IMF}$: 300 M$_{\odot}$',\
         r'BPASS SPs no binaries, ${\rm M}^{up}_{\rm IMF}$: 100 M$_{\odot}$']
fold = ['binary_cont_300','single_cont_100']

new_u = [0,5,7,8,9,10] # the table covers a smaller range of u somehow

cmap = plt.get_cmap('viridis')
acount = 0
alp = 1.
for s in range(2):
    calz = []
    sindx = s
    print(sindx,imf[s],fold[s],names[s])
    ax = plt.subplot(gs1[s])
    for z in metallicity:
        ncount = 0
        indx = metallicity.index(z)
        if z == 0.1 or z == 0.2:
            for zn in [z,0.3,0.5]:
                spec = get_cloudy_spec(fold[s],imf[s],a,z,zneb=zn)
                df = get_cloudy_table(fold[s],a,z,zneb=zn)

                civ = df['CIV1548'].values.copy() + df['CIV1551'].values.copy()
                ciii = df['[CIII]1907'].values.copy() + df['CIII]1909'].values.copy() 
                ratio = civ / ciii
                
                m36,m45 = get_cloudy_norm_spec(spec)
                m36,m45 = m36[new_u], m45[new_u] # the table covers a smaller range of U
                
                ax.plot(m36-m45,ratio,color=colors[indx],lw=2.,ls=zmet[ncount],zorder=1)
                ax.scatter(m36[0]-m45[0],ratio[0],marker='o',lw=2.,facecolor='none',edgecolor=colors[indx],s=40)
                ax.scatter(m36[1]-m45[1],ratio[1],color=colors[indx],s=20)
                ax.scatter(m36[5]-m45[5],ratio[5],color=colors[indx],s=20)
                ncount += 1
        else:
            lwe = [2.5,3.]
            for zn in [0.3,0.5]:
                spec = get_cloudy_spec(fold[s],imf[s],a,z,zneb=zn)
                df = get_cloudy_table(fold[s],a,z,zneb=zn)

                civ = df['CIV1548'].values.copy() + df['CIV1551'].values.copy()
                ciii = df['[CIII]1907'].values.copy() + df['CIII]1909'].values.copy() 
                ratio = civ / ciii
                
                m36,m45 = get_cloudy_norm_spec(spec)
                m36,m45 = m36[new_u], m45[new_u] # the table covers a smaller range of U
                
                if ncount == 0: plt.plot(m36-m45,ratio,color=colors[indx],lw=0.6)
                ax.plot(m36-m45,ratio,color=colors[indx],lw=lwe[ncount],ls=zmet[ncount+1])
                ax.scatter(m36[0]-m45[0],ratio[0],marker='o',lw=2.,facecolor='none',edgecolor=colors[indx],s=40)
                ax.scatter(m36[1]-m45[1],ratio[1],color=colors[indx],s=20)
                ax.scatter(m36[5]-m45[5],ratio[5],color=colors[indx],s=20)
                ncount += 1
        alp -= 0.15
    acount += 1
    
#     txt = ax.text(0.95,0.9,'attenuation',ha='right',transform=ax.transAxes,fontsize=14)
#     ax.plot([1.03+0.02,1.03-0.01+0.09],[42,42],color='k',lw=4)
#     ax.scatter(1.03+0.048,42,marker=9,color='k',s=130)

#     ax.text(-0.375,53.5,'Stellar Z:',fontsize=14)
#     y = 40
#     for z in metallicity:
#         indx = metallicity.index(z)
#         ax.plot([-0.37,-0.21],[y,y],color=colors[indx],lw=2.)
#         ax.text(-0.14,y-0.078*y,'%s Z$_{\odot}$'%(z),fontsize=14)
#         y -= 0.265*y

#     count,y = 0,2.25
#     ax.text(0.762,3,'Nebular Z:',fontsize=14)
#     for zn in ['fix to Z$_{\star}$','0.3 Z$_{\odot}$','0.5 Z$_{\odot}$']:
#         ax.plot([0.76,0.93],[y,y],color='k',ls=zmet[count],lw=zlw[count])
#         ax.text(1.0,y-0.078*y,zn,fontsize=14)
#         y -= 0.265*y
#         count += 1
    
    ax.set_title(names[s],fontsize=15.5)
    ax.errorbar(irac[0],c32[0],xerr=iracerr,marker='*',color='#962A13',ms=25,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=8)
    ax.plot([irac[0],irac[0]],[c32[1],c32[0]],color='k')
    ax.errorbar(irac[0],c32[1],xerr=iracerr,marker='*',color='#F4DCD7',ms=25,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=6)
    ax.scatter(irac[0],c32[1],marker='*',s=700,color='none',edgecolor='k',zorder=4)
    
    # error bars to left
#     ax.errorbar(0.05,c32[0],yerr=c32err,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=8)   
#     ax.errorbar(0,c32[1],yerr=c32err,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=6)
    
    xmin = 0.535 # axis units, not data lol)
    xmax = 0.8082 # axis units, not data lol)
    ymin= c32[1]-c32err[0][0]
    ymax= c32[1]+c32err[1][0]
    ax.axhspan(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,facecolor='#E1A396',zorder=0,alpha=0.5)
    
    xmin = 0.535 # axis units, not data lol)
    xmax = 0.8082 # axis units, not data lol)
    ymin= c32[0]-c32err[0][0]
    ymax= c32[0]+c32err[1][0]
    ax.axhspan(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,facecolor='#962A13',zorder=0,alpha=0.5)
    
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.yscale('log')

    ax.set_xlim(-0.48,1.41)
    ax.set_ylim(0.0008,1.7)
    ax.set_ylabel(r'CIV$\,\lambda$1548,1551 / CIII]$\,\lambda$1907,1909',fontsize=14,labelpad=1)
    ax.set_xlabel('[3.6]$\endash$[4.5] [mag]',fontsize=14)

print()

plt.tight_layout()
plt.savefig(f'plots-data/measurements.pdf')
plt.savefig('/home/aibhleog/Documents/papers/mine/figures/irac-ratio.pdf')
plt.show()
plt.close()



# fig, ax1 = plt.subplots(1,1,figsize=(11,14))
# # gs1 = gridspec.GridSpec(3,2)
# gs1 = gridspec.GridSpec(2,2)
# u = np.arange(-3.5,-1.4,0.2)
# colors = ['#707B7C','#5DADE2','#1F618D','#6C3483']
# metallicity = [0.1,0.2,0.3,0.5]
# zmet = ['-','--',':']; zlw = [1.5,1.5,2.]
# age = 7 #[7,7.477,8]
# a = age
# sub = [0,2,1,3]
# imf = [300,300,100,100]
# names = [r'BPASS SPs with binaries, ${\rm M}^{up}_{\rm IMF}$: 300 M$_{\odot}$',\
#          r'BPASS SPs no binaries, ${\rm M}^{up}_{\rm IMF}$: 300 M$_{\odot}$',\
#          r'BPASS SPs with binaries, ${\rm M}^{up}_{\rm IMF}$: 100 M$_{\odot}$',\
#          r'BPASS SPs no binaries, ${\rm M}^{up}_{\rm IMF}$: 100 M$_{\odot}$',\
#          r'Starburst99 SPs, ${\rm M}^{up}_{\rm IMF}$: 100 M$_{\odot}$','']
# fold = ['binary_cont_300','single_cont_300','binary_cont_100','single_cont_100','','']

# new_u = [0,5,7,8,9,10] # the table covers a smaller range of u somehow

# cmap = plt.get_cmap('viridis')
# acount = 0
# alp = 1.
# for s in sub:
#     calz = []
#     sindx = sub.index(s)
#     print(sindx,imf[s],fold[s],names[s])
#     ax = plt.subplot(gs1[s])
#     for z in metallicity:
#         ncount = 0
#         indx = metallicity.index(z)
#         if z == 0.1 or z == 0.2:
#             for zn in [z,0.3,0.5]:
#                 spec = get_cloudy_spec(fold[s],imf[s],a,z,zneb=zn)
#                 df = get_cloudy_table(fold[s],a,z,zneb=zn)

#                 #ciii_ew = get_cloudy_ew('ciii',spec)
#                 civ = df['CIV1548'].values.copy() + df['CIV1551'].values.copy()
#                 ciii = df['[CIII]1907'].values.copy() + df['CIII]1909'].values.copy() 
#                 ratio = civ / ciii
                
#                 m36,m45 = get_cloudy_norm_spec(spec)
#                 m36,m45 = m36[new_u], m45[new_u] # the table covers a smaller range of U
                
#                 ax.plot(m36-m45,ratio,color=colors[indx],lw=2.,ls=zmet[ncount],zorder=1)
#                 ax.scatter(m36[0]-m45[0],ratio[0],marker='o',lw=2.,facecolor='none',edgecolor=colors[indx],s=40)
#                 ax.scatter(m36[1]-m45[1],ratio[1],color=colors[indx],s=20)
#                 ax.scatter(m36[5]-m45[5],ratio[5],color=colors[indx],s=20)
#                 ncount += 1
#         else:
#             lwe = [2.5,3.]
#             for zn in [0.3,0.5]:
#                 spec = get_cloudy_spec(fold[s],imf[s],a,z,zneb=zn)
#                 df = get_cloudy_table(fold[s],a,z,zneb=zn)

#                 #ciii_ew = get_cloudy_ew('ciii',spec)
#                 civ = df['CIV1548'].values.copy() + df['CIV1551'].values.copy()
#                 ciii = df['[CIII]1907'].values.copy() + df['CIII]1909'].values.copy() 
#                 ratio = civ / ciii
                
#                 m36,m45 = get_cloudy_norm_spec(spec)
#                 m36,m45 = m36[new_u], m45[new_u] # the table covers a smaller range of U
                
#                 if ncount == 0: plt.plot(m36-m45,ratio,color=colors[indx],lw=0.6)
#                 ax.plot(m36-m45,ratio,color=colors[indx],lw=lwe[ncount],ls=zmet[ncount+1])
#                 ax.scatter(m36[0]-m45[0],ratio[0],marker='o',lw=2.,facecolor='none',edgecolor=colors[indx],s=40)
#                 ax.scatter(m36[1]-m45[1],ratio[1],color=colors[indx],s=20)
#                 ax.scatter(m36[5]-m45[5],ratio[5],color=colors[indx],s=20)
#                 ncount += 1
#         alp -= 0.15
#     acount += 1
    
#     txt = ax.text(0.95,0.9,'attenuation',ha='right',transform=ax.transAxes,fontsize=14)
#     ax.plot([1.03+0.02,1.03-0.01+0.09],[42,42],color='k',lw=4)
#     ax.scatter(1.03+0.048,42,marker=9,color='k',s=130)

#     ax.text(-0.375,53.5,'Stellar Z:',fontsize=14)
#     y = 40
#     for z in metallicity:
#         indx = metallicity.index(z)
#         ax.plot([-0.37,-0.21],[y,y],color=colors[indx],lw=2.)
#         ax.text(-0.14,y-0.078*y,'%s Z$_{\odot}$'%(z),fontsize=14)
#         y -= 0.265*y

#     count,y = 0,2.25
#     ax.text(0.762,3,'Nebular Z:',fontsize=14)
#     for zn in ['fix to Z$_{\star}$','0.3 Z$_{\odot}$','0.5 Z$_{\odot}$']:
#         ax.plot([0.76,0.93],[y,y],color='k',ls=zmet[count],lw=zlw[count])
#         ax.text(1.0,y-0.078*y,zn,fontsize=14)
#         y -= 0.265*y
#         count += 1
    
#     ax.set_title(names[s],fontsize=15.5)
#     ax.errorbar(irac[0],rwc2[0],xerr=iracerr,marker='*',color='#F4DCD7',ms=25,ecolor='k',capsize=3,capthick=1.5,zorder=3)
#     ax.scatter(irac,rwc2,marker='*',s=700,color='none',edgecolor='k',zorder=4)
#     ax.plot([irac[0],irac[0]],[c32[0],rwc2[0]+0.3],ls=':',color='k',zorder=5)
#     ax.errorbar(irac[0],c32[0],xerr=iracerr,marker='*',color='#962A13',ms=25,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=6)
    
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)
#     plt.yscale('log')

#     ax.set_xlim(-0.48,1.41)
#     ax.set_ylim(0.79,85)
#     ax.set_ylabel('W$_{CIII],0}$ [$\AA$]',fontsize=14,labelpad=1)
#     ax.set_xlabel('[3.6]$\endash$[4.5] [mag]',fontsize=14)
    
# '''# -------------------------------------------------------------------------------------- #    
# s = 4
# ax = plt.subplot(gs1[s])
# acount = 0
# s99met = [0.05,0.1,0.4]
# vcolors = [cmap(i) for i in np.linspace(0,1,7)]
# cs = ['#707B7C','#5DADE2',vcolors[4],'#1F618D','#6C3483']
# alp = 1.

# print('Starburst99, 100')
# calz = []
# for z in s99met:
#     ncount = 0
#     indx = s99met.index(z)
#     for zn in [z,0.3,0.5]:
#         if zn == 0.05: pass
#         else:
#             spec = get_cloudy_spec('single_s99_kroupa',100,a,z,zneb=zn)
#             pre_spec = get_cloudy_spec('single_s99_kroupa',100,a,z,zneb=zn)

#             ciii_ew = get_cloudy_ew('ciii',spec)
#             m36,m45 = get_cloudy_norm_spec(spec)
            
#             ax.plot(m36-m45,ciii_ew,color=cs[indx],lw=2.,ls=zmet[ncount],zorder=1)
#             ax.scatter(m36[0]-m45[0],ciii_ew[0],marker='o',lw=2.,facecolor='none',edgecolor=cs[indx],s=40)
#             ax.scatter(m36[5]-m45[5],ciii_ew[5],color=cs[indx],s=20)
#             ax.scatter(m36[10]-m45[10],ciii_ew[10],color=cs[indx],s=20)

#         ncount += 1
#     alp -= 0.15

# ax.text(-0.375,53.5,'Stellar Z:',fontsize=14)
# y,zs = 40,[0.05,0.1,0.1,0.1,0.4]
# nameimfs = ['',r', $\alpha$: -1.7',r', $\alpha$: -2.0','','']
# for i in range(5):
#     ax.plot([-0.37,-0.21],[y,y],color=cs[i],lw=2.)
#     ax.text(-0.14,y-0.078*y,'%s Z$_{\odot}$'%(zs[i])+nameimfs[i],fontsize=14)
#     y = y - 0.265*y

# count,y = 0,2.25
# ax.text(0.762,3,'Nebular Z:',fontsize=14)
# for zn in ['fix to Z$_{\star}$','0.3 Z$_{\odot}$','0.5 Z$_{\odot}$']:
#     ax.plot([0.76,0.93],[y,y],color='k',ls=zmet[count],lw=zlw[count])
#     ax.text(1.0,y-0.078*y,zn,fontsize=14)
#     y -= 0.265*y
#     count += 1

# txt = ax.text(0.95,0.9,'attenuation',ha='right',transform=ax.transAxes,fontsize=14)
# ax.plot([1.03+0.02,1.03-0.01+0.09],[42,42],color='k',lw=4)
# ax.scatter(1.03+0.048,42,marker=9,color='k',s=130)

# imfs = [1.7,2.]
# for i in range(2):
#     ncount = 0 
#     for zn in [0.1,0.3,0.5]:
#         spec = get_cloudy_spec('single_s99_miscIMF',100,a,0.1,zneb=zn,imf=imfs[i])
#         ciii_ew = get_cloudy_ew('ciii',spec)
#         m36,m45 = get_cloudy_norm_spec(spec)
#         ax.plot(m36-m45,ciii_ew,color=cs[i+3],lw=2.,ls=zmet[ncount],zorder=1)
#         ax.scatter(m36[0]-m45[0],ciii_ew[0],marker='o',lw=2.,facecolor='none',edgecolor=cs[i+3],s=40)
#         ax.scatter(m36[5]-m45[5],ciii_ew[5],color=cs[i+3],s=20)
#         ax.scatter(m36[10]-m45[10],ciii_ew[10],color=cs[i+3],s=20)
        
#         ncount += 1
    
# ax.set_title(names[s],fontsize=16.5)
# ax.errorbar(irac[0],rwc2[0],xerr=iracerr,marker='*',color='#F4DCD7',ms=25,ecolor='k',capsize=3,capthick=1.5,zorder=3)
# ax.scatter(irac,rwc2,marker='*',s=700,color='none',edgecolor='k',zorder=4)
# ax.plot([irac[0],irac[0]],[c32[0],rwc2[0]+0.3],ls=':',color='k',zorder=5)
# ax.errorbar(irac[0],c32[0],xerr=iracerr,marker='*',color='#962A13',ms=25,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=6)

# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.yscale('log')

# ax.set_xlim(-0.48,1.41)
# ax.set_ylim(0.79,85)
# ax.set_ylabel('W$_{CIII],0}$ [$\AA$]',fontsize=14,labelpad=1)
# ax.set_xlabel('[3.6]$\endash$[4.5] [mag]',fontsize=14)
# # -------------------------------------------------------------------------------------- #
    
# s = 5
# ax = plt.subplot(gs1[s])
# cmap = plt.get_cmap('viridis')
# colors = [cmap(i) for i in np.linspace(0,1,7)]

# print('agn')
# agn_zneb = [0.05,0.1,0.15,0.2,0.3,0.4,0.5]
# calz = []
# for i in range(len(agn_zneb)):
#     agn_spec = get_cloudy_spec_agn('agn',agn_zneb[i])

    
    
#     agn_m36,agn_m45 = get_cloudy_norm_spec(agn_spec,norm=False)
#     agn_ciii = get_cloudy_ew('ciii',agn_spec)
#     ax.plot(agn_m36-agn_m45,agn_ciii,color=colors[i],label='%s Z$_{\odot}$'%(agn_zneb[i]))
#     ax.scatter(agn_m36[0]-agn_m45[0],agn_ciii[0],marker='o',lw=2.,facecolor='none',edgecolor=colors[i],s=40)
#     ax.scatter(agn_m36[5]-agn_m45[5],agn_ciii[5],color=colors[i],s=20)
#     ax.scatter(agn_m36[10]-agn_m45[10],agn_ciii[10],color=colors[i],s=20)
    
# txt = ax.text(0.65,0.5,'attenuation',ha='left',transform=ax.transAxes,fontsize=14)
# ax.plot([0.94+0.02,0.94+0.068],[6.5,6.5],color='k',lw=4)
# ax.scatter(0.94+0.048,6.5,marker=9,color='k',s=130)

# ax.set_title('AGN Models',fontsize=16.5)
# ax.errorbar(irac[0],rwc2[0],xerr=iracerr,marker='*',color='#F4DCD7',ms=25,ecolor='k',capsize=3,capthick=1.5,zorder=3)
# ax.scatter(irac,rwc2,marker='*',s=700,color='none',edgecolor='k',zorder=4)
# ax.plot([irac[0],irac[0]],[c32[0],rwc2[0]+0.3],ls=':',color='k',zorder=5)
# ax.errorbar(irac[0],c32[0],xerr=iracerr,marker='*',color='#962A13',ms=25,ecolor='k',markeredgecolor='k',capsize=3,capthick=1.5,zorder=6)
    
# plt.yscale('log')
# ax.set_xlim(-0.48,1.41)
# ax.set_ylim(0.79,85)
# ax.set_ylabel('W$_{CIII],0}$ [$\AA$]',fontsize=14,labelpad=1)
# ax.set_xlabel('[3.6]$\endash$[4.5] [mag]',fontsize=14)
# ax.text(-0.4,3.7,'Nebular Z:',ha='left',fontsize=14.5)
# ax.legend(fontsize=14,frameon=False,labelspacing=0.15,ncol=2)

# # -------------------------------------------------------------------------------------- #'''

# plt.tight_layout()
# plt.show()
# plt.close()
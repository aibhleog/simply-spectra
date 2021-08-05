'''
Reading in the Cloudy line flux tables I made to investigate 
the potential differences between the standard Cloudy resolution
(R=300) to a higher resolution (R=1000).

'''

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

from cloudy_func import * # written by TAH for her Cloudy models

# let's start with something basic
Z = 0.1 # 10% Solar
age = 7 # log(age)
U = -2.5

# reading in tables
r300 = get_cloudy_table(f'r300',age,Z,mhome='testing')
r1000 = get_cloudy_table(f'r1000',age,Z,mhome='testing')

# resetting indexing
r300.reset_index(inplace=True,drop=True)
r1000.reset_index(inplace=True,drop=True)

# getting list of lines
lines = r300.columns.values[5:] # just the emission lines

# printing the line fluxes (relative to Hbeta) for each line for both Rs
print('''Looking at the differences in flux for the two resolutions.

Showing\t\t R=300, \t R=1000''')
for line in lines:
    print(f'{line}: \t', r300.loc[2,line], '\t', r1000.loc[2,line])
    
# looking at the percent differences
print('\nPercent difference')
for line in lines:
    print(f'{line}: \t', 
          round(100*abs(r300.loc[2,line]-r1000.loc[2,line])/r300.loc[2,line],3))
    
    

# plotting comparison around CIV
lims = {-1.5:[5,50],-2.5:[0.7,2.5],-3.5:[0.05,0.27]}
r300_spec = get_cloudy_spec(f'r300',300,age,Z,mhome='testing',ioni=U)
r1000_spec = get_cloudy_spec(f'r1000',300,age,Z,mhome='testing',ioni=U)

plt.figure(figsize=(7,5))
plt.plot(r300_spec.wavelength.values*1e4,r300_spec.spectrum,label='R~300')
plt.plot(r1000_spec.wavelength.values*1e4,r1000_spec.spectrum,label='R~1000')

ymin,ymax = lims[U]

plt.xlim(1506,1589)
plt.ylim(ymin,ymax)
plt.xlabel('wavelength [$\AA$]')
plt.ylabel(r'f$_{\nu}$')
plt.text(0.04,0.83,f'logU = {U}, \n$Z$ = {Z}'+'$Z_{\odot}$',
        transform=plt.gca().transAxes,fontsize=16)

plt.tight_layout()
plt.show()
plt.close()


# FOR MAKING GIFS
# for Z in [0.1,0.2,0.3,0.5]:
#     for U in [-3.5,-2.5,-1.5]:
#         r300_spec = get_cloudy_spec(f'r300',300,age,Z,mhome='testing',ioni=U)
#         r1000_spec = get_cloudy_spec(f'r1000',300,age,Z,mhome='testing',ioni=U)

#         plt.figure(figsize=(7,5))
#         plt.plot(r300_spec.wavelength.values*1e4,r300_spec.spectrum,label='R~300')
#         plt.plot(r1000_spec.wavelength.values*1e4,r1000_spec.spectrum,label='R~1000')

#         ymin,ymax = lims[U]

#         plt.xlim(1506,1589)
#         plt.ylim(ymin,ymax)
#         plt.xlabel('wavelength [$\AA$]')
#         plt.ylabel(r'f$_{\nu}$')
#         plt.text(0.04,0.83,f'logU = {U}, \n$Z$ = {Z}'+'$Z_{\odot}$',
#                 transform=plt.gca().transAxes,fontsize=16)

#         plt.tight_layout()
#         plt.savefig(f'plots-data/resolution-comparisons/binary-age{age}z{Z}zneb{Z}u{U}_300.png',dpi=300,facecolor='w')
#         plt.show()
#         plt.close()



# saving line ratios to file
# --------------------------
myfile = open(f'plots-data/resolution-comparisons/line-ratios-z{Z}u{U}.txt','w')

# printing the line fluxes (relative to Hbeta) for each line for both Rs
print('''Looking at the differences in flux for the two resolutions.

Showing\t\t R=300, \t R=1000''',file=myfile)
for line in lines:
    print(f'{line}: \t', r300.loc[2,line], '\t', r1000.loc[2,line],file=myfile)
    
# looking at the percent differences
print('\nPercent difference',file=myfile)
for line in lines:
    print(f'{line}: \t', 
          round(100*abs(r300.loc[2,line]-r1000.loc[2,line])/r300.loc[2,line],3),file=myfile)

# closing file
myfile.close()
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
    

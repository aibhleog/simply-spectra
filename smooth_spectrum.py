'''
Script to play with smoothing a 1D spectrum using different filtering windows.
May also play with a double gaussian smoothing window.

'''

import numpy as np
import matplotlib.pyplot as plt
import smoothing as sm # smoothing module by TAH
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy.stats import norm

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'


def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian function, to be used in the quick optimized extraction
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset


id_gal = 'z7_GND_42912'
date = '2020feb27'

# -- prompting user for galaxy ID name & date of observation -- #
# -- ID of galaxy
# id_gal = input('Type ID name & press Enter:  ')
# # making sure something was entered
# assert len(id_gal) > 0, f"Need to enter an ID name. "\
# 	f"Currently, the ID name is: \n \t\t  id_gal = {id_gal}"
# print(f'ID entered: "{id_gal}"',end='\n\n')

# # -- Date of Observation
# date = input('Type date & press Enter:  ')
# assert len(date) > 0, f"Need to enter the date the data were taken. "\
# 	f"Currently, the date is: \n \t\t  date = {date}"
# print(f'Date entered: "{date}"',end='\n\n')
# # ------------------------------------------------------------- #


# reading in data
print(f'Reading in data for {id_gal}')
path = f'spectral-secrets/vis/{date}/1D/'
spec = pd.read_csv(path + id_gal + '_spectrum.txt',sep='\s+')
print(f'File contents: \n{spec.head()}',end='\n\n')


# PLAYING WITH SMOOTHING
# first, let's look at different smoothing windows
plt.figure(figsize=(13,4.5))
gs = gridspec.GridSpec(1,2,width_ratios=[1,0.4],wspace=0.05)

ax = plt.subplot(gs[0]) # the spectrum
ax2 = plt.subplot(gs[1]) # the histogram

x = np.linspace(-4,4,100)
ax.step(spec.wave,spec.flux/1e-18,where='mid',lw=1)
ax.fill_between(spec.wave,spec.ferr/1e-18,spec.ferr*-1/1e-18,alpha=0.2)
i = 0


# comment out the different paragraphs to see different smoothing effects
# -----------------------------------------------------------------------

# smoothing based upon window size
# for w in np.arange(5,15,2): # odd numbers only
#     smoothed = sm.smooth(spec.flux,window_len=w)
#     ax.step(spec.wave,smoothed/1e-18,where='mid',label=f'w={w}pix')
    
# #     smoothed_err = sm.smooth(spec.ferr,window_len=w)
#     smu,ssigma = norm.fit(smoothed/spec.ferr)
#     label = f'$\mu$={round(smu,2)}\n$\sigma$={round(ssigma,2)}'
#     ax2.plot(x,gaussian(x,smu,500,ssigma,0),label=label,color=f'C{i+1}')
#     i += 1
    
    
# smoothing based upon filtering window (not gaussian)
# for f in ['hann','hamming','bartlett','blackman']:
#     smoothed = sm.smooth(spec.flux,window_len=5,window=f)
#     ax.step(spec.wave,smoothed/1e-18,where='mid',label=f'f={f}')
    
# #     smoothed_err = sm.smooth(spec.ferr,window_len=5,window=f)
#     smu,ssigma = norm.fit(smoothed/spec.ferr)
#     label = f'$\mu$={round(smu,2)}\n$\sigma$={round(ssigma,2)}'
#     ax2.plot(x,gaussian(x,smu,500,ssigma,0),label=label,color=f'C{i+1}')
#     i += 1


# smoothing based upon sigma (gaussian)
for s in np.arange(0.6,1.2,0.2):
    smoothed = sm.smooth(spec.flux,window_len=3,window='gaussian',sigma=s)
    ax.step(spec.wave,smoothed/1e-18,where='mid',label=f'$\sigma_G$={s}')
    
#     smoothed_err = sm.smooth(spec.ferr,window_len=3,window='gaussian',sigma=s)
    smu,ssigma = norm.fit(smoothed/spec.ferr)
    label = f'$\mu$={round(smu,2)}\n$\sigma$={round(ssigma,2)}'
    ax2.plot(x,gaussian(x,smu,0.37,ssigma,0),label=label,color=f'C{i+1}')
    i += 1

# -----------------------------------------------------------------------

ax.legend(ncol=3,handlelength=1)

ax.set_ylabel('flux [10$^{-18}$ erg/s/cm$^2$/$\AA$]')
ax.set_xlabel('observed wavelength [microns]')
ax.set_ylim(-5,5.8)
ax.set_xlim(1.315,1.3215)

# ------------------------------ 
# adding a histogram for the S/N
snr = spec.flux.values/spec.ferr.values
mu,sigma = norm.fit(snr)
label = f'$\mu$={round(mu,2)}\n$\sigma$={round(sigma,2)}'

print(f'\nMeasured -- \tMean: {round(mu,3)}, \tSigma: {round(sigma,3)}')
print('Ideal ----- \tMean: 0.0, \tSigma: 1.0',end='\n\n')

ax2.hist(snr,bins=20,alpha=0.6,color='C0',density=1)
ax2.plot(x,gaussian(x,mu,0.37,sigma,0),label=label,color='C0')

ax2.legend(fontsize=13,handlelength=0.7)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.set_xlabel('S/N histograms')
ax2.set_ylim(0,0.45)
ax2.set_xlim(-4,5)

# plt.show()
plt.close()




# ========================= #
# DOUBLE GAUSSIAN SMOOTHING #
# ========================= #

plt.figure(figsize=(13,4.5))
gs = gridspec.GridSpec(1,2,width_ratios=[1,0.4],wspace=0.05)

ax = plt.subplot(gs[0]) # the spectrum
ax2 = plt.subplot(gs[1]) # the histogram

x = np.linspace(-4,4,100)
ax.step(spec.wave,spec.flux/1e-18,where='mid',lw=1)
ax.fill_between(spec.wave,spec.ferr/1e-18,spec.ferr*-1/1e-18,alpha=0.2)

width,window_len = 10,5
smoothed1 = sm.double_smooth(spec.flux,width=width,window_len=window_len,sigma=0.8)
ax.step(spec.wave,smoothed1/1e-18,where='mid',label=f'sep={width}pix \nwindow={window_len}pix')

width,window_len = 12,3
smoothed2 = sm.double_smooth(spec.flux,width=width,window_len=window_len,sigma=0.8)
ax.step(spec.wave,smoothed2/1e-18,where='mid',label=f'sep={width}pix \nwindow={window_len}pix')

window_len = 3
smoothed3 = sm.smooth(spec.flux,window_len=window_len,window='gaussian',sigma=0.8)
ax.step(spec.wave,smoothed3/1e-18,where='mid',label=f'single gaussian \nwindow={window_len}pix')

ax.set_ylabel('flux [10$^{-18}$ erg/s/cm$^2$/$\AA$]')
ax.set_xlabel('observed wavelength [microns]')
ax.set_ylim(-5,5.8)
ax.set_xlim(1.315,1.3215)
ax.legend(ncol=3,fontsize=13,handlelength=0.7,loc=8)


# histogram
mu,sigma = norm.fit(snr)
label = f'$\mu$={round(mu,2)}\n$\sigma$={round(sigma,2)}'

print(f'\nMeasured -- \tMean: {round(mu,3)}, \tSigma: {round(sigma,3)}')
print('Ideal ----- \tMean: 0.0, \tSigma: 1.0',end='\n\n')

ax2.hist(snr,bins=20,alpha=0.6,color='C0',density=1)
ax2.plot(x,gaussian(x,mu,0.37,sigma,0),label=label)

for smoothed in [smoothed1,smoothed2,smoothed3]:
    smu,ssigma = norm.fit(smoothed/spec.ferr)
    label = f'$\mu$={round(smu,2)}\n$\sigma$={round(ssigma,2)}'
    ax2.plot(x,gaussian(x,smu,0.37,ssigma,0),label=label)

ax2.legend(fontsize=13,handlelength=0.7)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
ax2.set_xlabel('S/N histograms')
ax2.set_ylim(0,0.45)
ax2.set_xlim(-4,5)

plt.savefig(f'spectral-secrets/vis/{date}/1D/{id_gal}-double-gaussian-smoothing.pdf')
plt.show()
plt.close()
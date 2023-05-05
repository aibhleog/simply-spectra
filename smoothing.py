'''
SMOOTHING A 1D SPECTRUM
-----------------------
The original code for this comes from 
https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

I've modified it to use scipy & allow for a gaussian shape.

'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal.windows as w


def smooth(x,window_len=11,window='hann',sigma=2):
    '''
    Window options:
    'flat','gaussian','hann','hamming','bartlett','blackman'
    
    To access them, use w.[window](window_len)
    '''
    if window_len < 3: # too small, doesn't smooth
        return x

    # np.r_ turns a list into an array, in this case it works like np.concatenate
    # s will be longer than x based on the window_len (will need to crop it)
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': # generates a moving average
        w = np.ones(window_len)
    elif window == 'gaussian': # requires sigma input
        w = eval('w.' + window + '(window_len,sigma)')
    else:
        w = eval('w.' + window + '(window_len)')

    # the actual smoothing line
    y = np.convolve(w/w.sum(), s, mode='valid')
    
    # y will be larger than x based upon the window_len size
    # for now, I've written in a return that slices out the original section
    h = int(window_len/2)
    y = y[h:-h] # removing the increase in array size due to window_len
    return y


def see_smoothing():
    '''
    Just a function that can be used to plot what this can look 
    like for the different filters.  Feel free to change it up.
    '''
    window_len = 7
    x = np.linspace(2,4,50)
    y = 3*np.cos(x)**3
    yn = y + np.random.normal(0,0.5,len(x))

    plt.figure(figsize=(8,4.5)); plt.axis('off')
    plt.plot(x,y,label='true',lw=2,ls=':')
    plt.plot(x,yn,label='w/noise')
    plt.text(0.4,0.8,f'window size: {window_len}',transform=plt.gca().transAxes,fontsize=16)
    
    for f in ['flat','gaussian','hann','hamming','bartlett','blackman']:
        s = smooth(yn,window_len,window=f)
        plt.plot(x,s,label=f,lw=1)

    plt.legend(ncol=3,loc=3,bbox_to_anchor=(0.1,-0.5))
    plt.tight_layout()
    plt.show()
    plt.close()
    
def see_windows():
    '''
    Just a function that can be used to visualize the different
    smoothing windows.
    '''
    window_len = 19
    filts = ['flat', 'gaussian', 'hann', 'hamming', 'bartlett', 'blackman']

    plt.figure(figsize=(7,4.5)); plt.axis('off')
    plt.text(0.7,0.82,f'window size: {window_len}',transform=plt.gca().transAxes,fontsize=13.5)

    i = 0
    for f in ['flat', 'gaussian', 'hann', 'hamming', 'bartlett', 'blackman']:
        label = f # assigning it so I can add the sigma for gaussian

        if f == 'flat': 
            wfilt = np.ones(window_len)
        elif f == 'gaussian': 
            wfilt = eval('w.'+f+'(window_len,2)')
            label = f'{f}\n$\sigma$=2' # adding note about sigma
        else: 
            wfilt = eval('w.'+f+'(window_len)')

        plt.plot(wfilt,label=label,color=f'C{i+2}') # so the color matches the other plot
        i += 1

    plt.legend(ncol=3,loc=3,bbox_to_anchor=(0.02,-0.5))
    plt.tight_layout()
    plt.show()
    plt.close()
    

def double_smooth(x,window_len=5,width=11,sigma=2):
    '''
    This function smooths a spectrum using a double gaussian with a
    specified window size & sigma, separated by a width.
    '''
    if window_len < 3: # too small, doesn't smooth
        return x

    # making double gaussian
    window = 'gaussian' # requires sigma input
    w_one = eval('w.' + window + '(window_len,sigma)')
    
    if width % 2 == 0: w = np.concatenate((w_one,np.zeros(width),w_one,[0]))
    else: w = np.concatenate((w_one,np.zeros(width),w_one))
    
    
    # np.r_ turns a list into an array, in this case it works like np.concatenate
    # s will be longer than x based on the window_len (will need to crop it)
    full_window_len = len(w)
    s = np.r_[x[full_window_len-1:0:-1],x,x[-2:-full_window_len-1:-1]]
    

    # the actual smoothing line
    y = np.convolve(w/w.sum(), s, mode='valid')
    
    # y will be larger than x based upon the window_len size
    # for now, I've written in a return that slices out the original section
    h0 = int(window_len/2)
    h = int(full_window_len/2)
    y = y[h0:-(h+(h-h0))] # removing the increase in array size due to window_len+width
    return y
    
    

'''
In order to add in the rebinning of the uncertainty -- first, I plan to use a variant of
the rebinning code used by Dr Jane Rigby:

https://github.com/janerigby/jrr/blob/master/spec.py#L74

Next, I plan to incorporate the uncertainty rebinning by taking inspiration from this MUSE
python package's code:

https://mpdaf.readthedocs.io/en/stable/spectrum.html
https://github.com/musevlt/mpdaf/blob/master/lib/mpdaf/obj/data.py#L1387

and getting to work with pandas (which Dr Rigby's code uses).

'''


def rebin_spec(wave, specin, wavnew):
	'''	
	Given wavelength, a spectrum, and new wavelength array, this
    function resamples the spectrum to match new array.
    
        NOTE: this works fine for the flux density, but uncertainty needs to
	be treated differently!
	'''
	import numpy as np
	from pysynphot import observation
	from pysynphot import spectrum

	spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin,keepneg=True)
	f = np.ones(len(wave))
	filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
	obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
	return obs.binflux


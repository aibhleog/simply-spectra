'''
SMOOTHING A 1D SPECTRUM
-----------------------
The original code for this comes from 
https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

I've modified it to use scipy & allow for a gaussian shape.

'''

import numpy as np
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


def see_windows():
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
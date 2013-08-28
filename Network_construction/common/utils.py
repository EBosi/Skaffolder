#!/usr/bin/env python
"""
utils

Common library

Spare parts
"""

# Borrowed from: www.garyrobinson.net
def slice_it(li, cols=10):
    start = 0
    for i in xrange(cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop
        
def get_span(li, span=4):
    i = 0
    while i < len(li):
        yield li[i : i+span]
        i += span

# Borrowed from stackoverflow
# python-most-idiomatic-way-to-convert-none-to-empty-string
def xstr(s):
    if s is None:
        return ''
    return str(s)

# Borrowed from http://www.scipy.org/Cookbook/SignalSmooth
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """
    from numpy.core.numeric import ones
    import numpy
    
    x = numpy.array(x)

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
    #if len(x) < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def compress(x, span=10):
    return [i[0] for i in get_span(x, span)]

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def rangeColors(minimum, maximum, colorrange):
    '''
    Given a min and a max value, returns a dict with stepwise colors
    '''
    hexs = {}
    prev = '#FFFFFF'
    i = minimum
    for color in slice_it(colorrange, cols=maximum-minimum+1):
        if len(color) == 0:
            hexs[i] = prev
        else:
            hexs[i] = rgb_to_hex(tuple([int(round(x*255))
                              for x in color[-1][:3]])).upper()
        prev = hexs[i]
        i += 1
    
    return hexs
        
import numpy as np
from scipy import signal

def smooth(x,window_len=11,window="hamming"):
    """smooth the data using a window with requested size.
    
    Taken from http://wiki.scipy.org/Cookbook/SignalSmooth
    
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
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)] #Not working pretty well....
    

def whittaker(x,**kwargs):
    """
    Whittaker smoother from Eilers 2003
    
    Inputs
    ------
    
    x: vector with the values to smooth (equally spaced)
    
    Options
    -------
    
    Lambda: the smoothing coefficient, the higher the smoother. Default = 10^5.
    
    Outputs
    -------
    
    A vector containing the smoothed values
    
    """
    # optional parameters
    lam = kwargs.get('Lambda',1.0*10**5)
    
    # starting the algorithm
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = sparse.linalg.spsolve(Z, w*y)
    
    return z

#### FILTERING OF DATA WITH A BUTTERWORTH FILTER 
   
def spectrafilter(spectre,filtertype,fq,numtaps,columns):
    """
    This function filters the HF of the spectra
    It use a low pass butterworth filter
    y is an array of spectra constructed with the previous spectraarray function
    filtertype contains the string defining which type of filter you want
    Choose between low, high, bandstop, bandpass
    fq is the frequency of the periodic signal you try to erase
    if it is a bandpass or band stop, f must be an array containing the cutoff frequencies
    columns is an array defining which columns you want to treat
    first column of y must contain the frequency
    """
    
    # we already say what is the output array
    out = np.zeros(spectre.shape)
    
    # Butterworth band stop filter caracteristics
    a = spectre[1,0] - spectre[0,0]
    samplerate = 1/a  #Hertz
    nyq_rate = samplerate/2 # frequence Nyquist
    cutf = fq # cutoff frequency
    #bandwidth = 0.005 # largeur filtre, for band pass/stop filters
    numtaps = 1 # ordre du filtre...
    
    for i in range(len(columns)):
        y = spectre[:,columns[i]]
        if (filtertype == 'low') or (filtertype == 'high'):
            b, a = signal.butter(numtaps, [(cutf/nyq_rate)], btype = filtertype)
            out[:,columns[i]] = signal.filtfilt(b, a, y) # filter with phase shift correction
        else:
            b, a = signal.butter(numtaps, [(cutf[0]/nyq_rate),(cutf[1]/nyq_rate)], btype = filtertype)
            out[:,columns[i]] = signal.filtfilt(b, a, y) # filter with phase shift correction

    # Note forgetting to register the x axis        
    out[:,0] = spectre[:,0]
    
    return out
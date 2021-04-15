import numpy as np
import matplotlib.pyplot as plt

from scipy import signal
from scipy.fft import fftshift

import ipywidgets as widgets

"""
Summary :
    Class implementing an adaptive notch filter

History :
    Apr 05/21 - File created by Daniele Borio
"""

class adaptivenotch :

    def __init__(self, ka = 0.9, mu = 0) :
        """
        Summary :
            object constructor.
    
        Arguments :
            ka - the pole contraction factor
            mu - the adaptation step
        """
        
        # Complex zero of the filter 
        self.z0 = 0.0 + 0.0j
        
        # Output of the AR part of the filter
        self.xi = 0.0
        
        # Pole contraction factor of the filter
        self.ka = ka 
        
        # Adaptation step of the adaptive algorithm
        if mu != 0 :
            self.mu = mu
        else :
            self.mu = 0.25 * (1 - ka)
        
    
    def filter( self, x ) :
        """
        Summary :
            Perform filtering on the input samples, x.
    
        Arguments :
            x - numpy array with the input complex samples
            
        Returns:
            y - the output samples.
            z0 - the zeros of the filter computed during the adaptation process.
        """    
        # Estimate the signal energy
        energy = np.var( x )
        mu = self.mu / energy
    
        # Allocate the output vectors 
        z0 = np.zeros( len( x ), dtype = np.complex128 )
        y = np.zeros( len( x ), dtype = np.complex128 )

        for ii in range( len(x) ) :    
            # output of the autoregressive part
            xi = x[ ii ] + self.ka * self.z0 * self.xi

            # ouput of the moving average part
            y[ ii ] = xi - self.z0 * self.xi

            # update of the filter zero
            self.z0 += mu * y[ ii ] * np.conj( self.xi )
            z0[ ii ] = self.z0

            # update the filter state
            self.xi = xi
            
        return y, z0
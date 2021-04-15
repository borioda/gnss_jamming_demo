#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 11:21:24 2020

@author: daniele
"""


import numpy as np
from scipy.stats import chi2

"""
Summary:
    Evaluate the search space for the code acquisition using the time 
    domain FFT technique   

Arguments:
   sig       : [vector] the Galileo/GPS input signal
   locC      : [vector] the local code replica

   N         : [integer] the input signal and local code length

   Nd        : [integer] number of Doppler bin used for the search space

   DopStep   : [Hz] Doppler bin width in Hz

   fs        : [Hz] sampling frequency

   fi        : [Hz] intermediate frequency 

Return:
     sspace       : [matrix] the search space

History:
   Nov 06/12 - File created by Daniele Borio
   Aug 16/20 - Reviewed version
"""

def DftParallelCodePhaseAcquisition(sig, locC, N, Nd, DopStep, fs, fi) :
    
    fif = fi/fs             # normalized intermediate frequency
    deltaf = DopStep/fs     # normalized Doppler step

    sspace = np.zeros((Nd, N)) + 1j*np.zeros((Nd, N))

    F_CA = np.conjugate( np.fft.fft(locC) )   # fft of the local code
    t = np.arange(0, N, 1)                    # time index

    for ff in range(0, Nd) :  
        fc = fif + (ff + 1 - np.ceil( Nd / 2 ) )*deltaf        
        IQ_comp = np.exp(-2*1j*np.pi*fc*t)*sig
        
        X = np.fft.fft(IQ_comp)
        
        sspace[ff,:] = np.fft.ifft( X *F_CA )
    # end for
    sspace = np.real( sspace*np.conjugate(sspace) )    
    return sspace

"""
Summary:
    ResampleCode - Resamples a ranging code.

Arguments:
   code - The code sequence (+/-1)
   numPoints - The number of points in the resampled code
   fs - The sampling rate to use for resampling
   tau0 - The phase of the first sample of the resampled code
   fc - The chipping rate of the code.

Returns:
    resampledCode - resampled code
"""
def ResampleCode(code, numPoints, fs, tau0, fc ) :
    
    #time instants
    tVec = np.arange( 0, numPoints ) / fs
    indeces = (tVec*fc + tau0).astype( int ) % code.size
    
    resampledCode = code[ indeces ]
    
    return resampledCode

"""
Summary:
    Determine the decision threshold given a false alarm probability
    
Arguments:from scipy.stats import chi2
    Pfa_sys - false alarm probability at the system level
    Nel - number of elements in the search space
    K - number of non-coherent integrations 
Returns:
    Th - normalized decision threshold
"""
def GetNormalizedDecisionThreshold( Pfa_sys, Nel, K ) :

    # False alarm probability at the cell level
    Pfa_cell = 1 - ( 1 - Pfa_sys )**( 1 / Nel )    
    
    # Normalized Threshold - use the inverse survival function of a chi 
    # squared random variable
    Th = chi2.isf( Pfa_cell, 2 * K )
    
    return Th

"""
Summary:
   Function that evaluates the noise floor for setting the acquisition 
   threshold.

Arguments:
   y :         [vector] contains the input samples
   fs:         [scalar] sampling frequency
   fc:         [scalar] code rate
   fif:        [scalar] the intermedinp.random.binomial(n, p, 1000)ate frequency

Returns:
   sigma2:     [scalar] the estimated noise variance
"""

def NoiseVarianceEstimator( y, fs, fc, fif ) :

    # First generate a fictitious code
    clen = int( np.round( y.size / fs * fc ) )
  
    # A bipolar random code usually has good correlation properties
    code = np.sign( 2 * np.random.binomial(1, 0.5, clen) - 1 )
                                           
    # Resample the code
    loc = ResampleCode( code, y.size, fs, 0, fc )
	
    # Now compute the correlators (for a single Doppler value is enough )
    correlators = DftParallelCodePhaseAcquisition( y, loc, y.size, 1, 0, fs, fif )

    #% Down-sample them to get uncorrelated values:
    step = int( np.round( fs / fc ) )
    correlators = correlators[::step]

    # Finally my noise variance estimate
    sigma2 = np.mean( correlators ) / 2
    
    return sigma2
       
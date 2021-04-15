import numpy as np
from matplotlib import cm

import ipyvolume as ipv
import ipywidgets as widgets

# Custom library for the generation of the C/A code and for signal acquisition
from cacode import cacode
import AcquiUtils as ac
import adaptivenotch as anf


class CAFInteractive :
    """
    Summary:
        Object constructor.
    
    Arguments:
        prn - the prn number to acquire
        data - vector containing the samples to be used for signal acquisition
        signal_char - dictionary describing the signal characteristics (sampling frequency and IF frequency)
        acq_char - parameters used for acquisition (Doppler step, number of Doppler bin, number of non-coherent integrations
    """
    
    def __init__(self, prn, data, signal_char, acq_char ):
        
        # signal parameters
        self.fs = signal_char["fs"]
        self.fi = signal_char["fi"]
        
        # Generate the local code
        locBipolar = cacode( prn )  # The code is bipolar format

        self.fc = 1.023e6    # Rate of the code
        Tc = 0.001      # Coherent integration time in ms
        self.Nc = int( Tc * self.fs )

        # Resample the code at the data sampling frequency
        self.locC = ac.ResampleCode( locBipolar, self.Nc, self.fs, 0, self.fc )    
        
        # Output search space
        self.K = acq_char["K"]    # Number of non-coherent integrations
        
        # Number of Doppler bins
        self.Nd = acq_char["Nd"]
        
        # Doppler step
        self.DopStep = acq_char["DopStep"]
        
        # Cross-Ambiguity Function (CAF)
        self.sspace = np.zeros( (self.Nd, self.Nc) )   # Search space were the results will be stored
        
        # Store the input data
        self.data = data
        
        # Compute the CAF
        self.evaluate_caf( data )
        
        ############################### Graphic and interactive elements #############################
        
        # Allocate an IPV figure
        self.figure = ipv.figure(width=800, height=500)
        
        self.codeDl = np.arange(0, self.Nc) / self.fs  # Code delays
        self.Freq = (np.arange(0, self.Nd) - np.floor( self.Nd / 2 ) ) * self.DopStep

        index = np.argmax( self.sspace )
        DopInd = int( index / self.sspace.shape[1] ) 
        codInd = np.mod( index, self.sspace.shape[1] ) 
        maxVal = self.sspace.max()
       
        # down-sampling factor (for the code)
        dsf = 40

        # initial index in order not to miss the correlation peak
        pstart = codInd % dsf

        code_val, freq_val = np.meshgrid(self.codeDl[pstart::dsf], self.Freq)
        
        colormap = cm.coolwarm
        Z = self.sspace[:, pstart::dsf] / maxVal
        color = colormap(Z)

        self.caf = ipv.plot_surface(freq_val, Z, code_val * 1000, color = color[...,:3])
        
        # ipv.ylabel("NCAF")
        # ipv.xlabel("Doppler [Hz]")
        # ipv.zlabel("Delay [ms]")
        
        ipv.show()
        
        # Add the widgets for interference mitigation
        Mitigation = ["None", "ANF08", "FDCS", "FDPB", "TDCS", "TDPB"]
        
        self.ui = widgets.interact(self.update, Th = widgets.FloatSlider(value=5, min=1, max=20.0, step = 0.1,  
                                   description='$T_h$', continuous_update=False), mtype = Mitigation, 
                                   K = widgets.IntSlider(value=1, min=1, max=10, step = 1, continuous_update=False))
        
    """
    Summary:
        Compute the Cross-Ambiguity function (CAF) using the data stored in the object.
    
    Arguments:
        data - the sample to use for the computation of the CAF
    """
    def evaluate_caf( self, data ) :
        
        self.sspace = np.zeros( (self.Nd, self.Nc) )
        
        for ii in range( 0, self.K ) :
            y = data[ (ii * self.Nc):((ii+1) * self.Nc) ]    # use just 1 period of code at the time

            # Compute the search space for a single coherent integration epoch
            Tsspace = ac.DftParallelCodePhaseAcquisition( y, self.locC, self.Nc, self.Nd, self.DopStep, self.fs, self.fi ) 
    
            self.sspace = self.sspace + Tsspace   # Non-coherently accumulate the results
    
    """
    Summary:
        Update the interactive plot of the CAF including interference mitigation
    
    Arguments:
        Th - Threshold to be for pulse blanking
        mtype - Type of interference mitigation
    """
    def update( self, Th, mtype, K ) :
        
        self.K = K
        
        # processing depending on the type of technique
        if mtype == "None" :
            # simply recompute the CAF without any interference mitigation
            self.evaluate_caf( self.data )
        
        elif mtype == "ANF08" :
            # allocate a notch filter and "clean" the data
            notch = anf.adaptivenotch(0.8)
            data, _ = notch.filter( self.data )
            self.evaluate_caf( data )
        
        elif mtype == "FDCS" :
            fdata = np.fft.fft( self.data )
            fdata = fdata / (abs(fdata) + 1e-6)
            data = np.fft.ifft( fdata )
            self.evaluate_caf( data )
        
        elif mtype == "FDPB" : 
            fdata = np.fft.fft( self.data )
            fdata = fdata * (abs(fdata) < 10000 * Th)
            data = np.fft.ifft( fdata )
            self.evaluate_caf( data )
            print("Th = %f" % (10000 * Th))
        
        elif mtype == "TDCS" :
            self.evaluate_caf( self.data / (abs(self.data) + 1e-6) )
            
        elif mtype == "TDPB" :
            self.evaluate_caf( self.data * ( abs(self.data) < Th ) )
                              
        # update the CAF
        index = np.argmax( self.sspace )
        DopInd = int( index / self.sspace.shape[1] ) 
        codInd = np.mod( index, self.sspace.shape[1] ) 
        maxVal = self.sspace.max()
       
        # down-sampling factor (for the code)
        dsf = 40

        # initial index in order not to miss the correlation peak
        pstart = codInd % dsf

        code_val, freq_val = np.meshgrid(self.codeDl[pstart::dsf], self.Freq)
        
        Z = self.sspace[:, pstart::dsf] / maxVal
        
        colormap = cm.coolwarm
        Z = self.sspace[:, pstart::dsf] / maxVal
        color = colormap(Z)
        
        self.caf.x = freq_val.flatten()
        self.caf.y = Z.flatten()
        self.caf.z = code_val.flatten() / 0.001
        
        self.caf.color = color[...,:3].reshape(-1, 3)
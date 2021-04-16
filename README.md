# GNSS Jamming Demo
Demonstration of the effects of interference mitigation techniques on GNSS signal acquisition.

In this repository, you will find the material used to demonstrate the impact of five **interference mitigation techniques** on GNSS signal acquisition.
The techniques considered are:

- Adaptive Notch Filter (ANF)
- Frequency Domain Complex Signum (FDCS)
- Frequency Domain Pulse Blanking (FDPB)
- Time Domain Complex Signum (TDCS)
- Time Domain Pulse Blanking (TDPB)

The material has been prepared as part of the Institute of Navigation (ION) webinar ["GNSS interference mitigation: A measurement and position domain assessment"](https://www.ion.org/publications/webinars.cfm) associated to the corresponding [paper](https://onlinelibrary.wiley.com/doi/full/10.1002/navi.391), which is available open access.

Two Jupyter notebooks are part of the demo:

- **AdaptiveNotchDemo:** demonstrate the working principles of the ANF. In particular, the impact of the *pole contraction factor*, k<sub>&alpha;</sub>, and of the *adaptation step*, &delta; are analysed
- **InteractiveJamDemo:** demonstrate the impact of five interference mitigation techniques on the Cross-Ambiguity Function (CAF). The impact of the blanking threshold is also investigated

Three short datasets (50 ms each) are provided in the *data* folder to experiment with interference mitigation and GNSS signal acquisition.

- **JammerData.bin**: GPS and Galileo signals affected by a sweept jamming signal with sweep range around 10 MHz. In this case, jamming can be effectively mitigated using frequency domain techniques
- **JamData400.bin** and **JamData500.bin**: snapshots extracted from the first test (Test 1) described in our [paper](https://onlinelibrary.wiley.com/doi/full/10.1002/navi.391). In this case, the jamming signal sweeps more than 35 MHz. So the jamming signal periodically enters and exists the receiver band resulting into a series of pulses. In this case, jamming can be effectively mitigated using time domain approaches. Suffixes "400" and "500" indicate the number of seconds from the start of Test 1.

All datasets have been collected with a 10 MHz sampling frequency, IQ 8-bit sampling. GNSS signals in JamData400.bin and JamData500.bin are affected by a clock drift of about 9 kHz.

## Dependencies:
In addition to standard python libraries such as numpy and scipy, you need to install the following libraries to run the notebooks:

- [**ipyvolume**](https://github.com/maartenbreddels/ipyvolume): for 3D rendering
- [**ipympl**](https://github.com/matplotlib/ipympl): for interactive matplotlib figures
- [**ipywidgets**](https://ipywidgets.readthedocs.io/en/stable/)  

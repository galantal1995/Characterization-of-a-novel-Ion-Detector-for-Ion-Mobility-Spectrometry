In Ion Mobility Spectrometry the detection of single ion pulses is in interest which might be difficult due to the presence of the detecor noise sources. 
Three different position detection algorithms were developed and evaluated (maximum, centeroid, gaussian). 
According to the symulations, the Gaussian method proves to be the most efficient

Process description of the gaussian_method.py:

- Uses a noise data obtained by the CMOS ion detector
- Adds artificially generated ion pulses
- Examines the frequency domain of the sum and applies a low pass filter to reduce the noise (theoritically calculated).
- Sigma clipping to determine a reference line.
- Peak detection algorithm with analyzing 5 neighbouring data points.
- Applies a Gaussian fit to determine the position of the detected peaks.

In another code it was repeated with 10.000 different signals to determine the accuracy of the detection method.

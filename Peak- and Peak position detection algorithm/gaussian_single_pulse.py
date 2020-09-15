from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, fftfreq, ifft
from scipy.signal import butter,filtfilt, lfilter
import scipy.signal
from scipy import signal as sg
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning

class implemented_functions():
    def gaussian_curve(self, amplitude):
        signal_gaussian = amplitude * scipy.signal.gaussian(points_on_gaussian, 20)

        return signal_gaussian

    def fast_fourier_transform(self, sample):
        # Frequency domain representation
        fourierTransform = np.fft.fft(sample) / len(sample)  # Normalize amplitude
        fourierTransform = fourierTransform[range(int(len(sample) / 2))]  # Exclude sampling frequency
        tpCount = len(sample)
        values = np.arange(int(tpCount / 2))
        timePeriod = tpCount / samplingFrequency
        frequencies = values / timePeriod

        return fourierTransform, frequencies

    def gaussian_fit(self, x, height_a, center_a, width_a): # Fitting function
        return (height_a * np.exp(-(x - center_a) ** 2 / (2 * width_a ** 2)))

    def butter_lowpass_filter(self, data, cutoff, samplingFrequency, order):
        normal_cutoff = cutoff / nyq
        # Get the filter coefficients
        b, a = butter(order, normal_cutoff, btype='lowpass')
        filtered = filtfilt(b, a, data)
        return filtered

func = implemented_functions()


npzfile = np.load('noise.npz')

# Input values, loaded from the zip file that contains the measured data
v_count = 19     # uV pnp.array per count
t_count = 6.25   # us per count for integration time

series_cds_rcf   = npzfile['series_cds_rcf'] * v_count
num_loops = npzfile['loops']
num_frames = npzfile['frames_to_load'] - 1

#n=num_frames*num_loops
n = 11000
noise = np.zeros(shape=(n))
for i in range(n):
    noise[i] = series_cds_rcf[0,11,i]

# At what intervals time points are sampled
samplingInterval = 6.25e-6
# How many time points are needed i,e., Sampling Frequency
samplingFrequency = 1/samplingInterval
# Begin time period of the signals
beginTime = 0
# End time period of the signals
endTime = n*samplingInterval
# Time points
time = np.arange(beginTime, endTime, samplingInterval)

# Creating the signal
position_of_signal = np.array([5000])
amplitude_of_signal = np.array([3000])
shift = 0
points_on_gaussian = 121  # do it with 101, 121, 141 and so on

original_position = np.zeros(shape=(len(position_of_signal)))
for i in range(len(position_of_signal)):
    original_position[i] = (position_of_signal[i]+5)*samplingInterval - (shift/10)*samplingInterval
print("Original position: ",original_position[0]*1e6)

signal_shortx = np.zeros(shape=(int(points_on_gaussian/10)-1, len(position_of_signal)))
signal_shorty = np.zeros(shape=(int(points_on_gaussian/10)-1, len(position_of_signal)))

signal = np.linspace(0, 0, n)
for h in range(len(position_of_signal)):
    signal_gaussian = func.gaussian_curve(amplitude_of_signal[h])
    for i in range(len(signal_shortx)):
        signal_shortx[i,h] = time[i + position_of_signal[h]]
        signal_shorty[i,h] = signal_gaussian[10 *(i+1) + shift]
        signal[position_of_signal[h] + i] = signal_shorty[i, h]

ions = np.zeros(shape=(len(position_of_signal)))
for i in range(len(position_of_signal)):
    ions[i] = int(np.sum(signal_shorty[:,i])/76)

# Adding the signal to the noise
signal_noise = noise + signal
signal_noise = signal_noise - np.mean(signal_noise)

# Filter requirements.
cutoff = 30e3     # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * samplingFrequency  # Nyquist Frequency
order = 5     # sin wave can be approx represented as quadratic

filtered_signal = func.butter_lowpass_filter(signal, cutoff, samplingFrequency, order)
filtered_noise = func.butter_lowpass_filter(noise, cutoff, samplingFrequency, order)
filtered_signal_noise = func.butter_lowpass_filter(signal_noise, cutoff, samplingFrequency, order)

# Fast Fourier Transform
(fft_signal, freq_signal) = func.fast_fourier_transform(signal)
(fft_noise, freq_noise) = func.fast_fourier_transform(noise)
(fft_filtered_noise, freq_filtered_noise) = func.fast_fourier_transform(filtered_noise)
(fft_signal_noise, freq_signal_noise) = func.fast_fourier_transform(signal_noise)
(fft_filtered_signal_noise, freq_filtered_signal_noise) = func.fast_fourier_transform(filtered_signal_noise)

clipped = sigma_clip(filtered_signal_noise, sigma=3, maxiters=None)
std = 3 * np.std(clipped)

detected_peak = np.zeros(shape=(n))
for i in range(2,n-3):
    if np.average(filtered_signal_noise[i-2 : i+3]) > std:
      detected_peak[i] = filtered_signal_noise[i]

detected_x = np.array([])
detected_y = np.array([])
for i in range(len(detected_peak)):
    if detected_peak[i]>0:
        detected_x = np.append(detected_x,time[i])
        detected_y = np.append(detected_y,detected_peak[i])

for i in range(n-1):
    if detected_peak[i] == 0 and detected_peak[i + 1] != 0:
        first = i
        break

# Considered data for the gaussian fit
x = time[first - 10:first + 5]                   #Time axis, between the first detected peak-10 and last detected peak+10
y = filtered_signal_noise[first - 10:first + 5]  #Signal axis, between the first detected peak-10 and last detected peak+10

height_a = np.max(y)
center_a = float(time[np.argwhere(detected_peak == np.max(detected_peak[first - 10:first + 5]))])
width_a = 15 * 1e-6
popt,pcov = curve_fit(func.gaussian_fit, x, y, p0=[height_a, center_a, width_a])
fitted_center = popt[1]


print("Found position: ", fitted_center * 1e6)
print('Time shift: ', (original_position[0] - fitted_center) * 1e6)

# --------------------------------------
#               PLOTTING
# --------------------------------------

# Plotting the signals
plt.figure(1)
ax1 = plt.subplot(211)
ax1.set_title('One ion pulse generated in %d frames with 6.25μs integration time' %(n), fontweight="bold", fontsize=14)
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Signal [μV]', fontsize=12)
ax1.plot(time, signal)
plt.grid(True)

ax2 = plt.subplot(223)
ax2.set_title('Ion pulse placed between the %dth and %dth frame' %(position_of_signal[-1], position_of_signal[-1]+len(signal_shortx)), fontweight="bold", fontsize=14)
ax2.set_xlabel('Time [s]', fontsize=12)
ax2.set_ylabel('Signal [μV]', fontsize=12)
ax2.plot(np.linspace(time[position_of_signal[-1] - 1]-shift/10*samplingInterval, time[position_of_signal[-1] + len(signal_shortx)]-shift/10*samplingInterval, points_on_gaussian), signal_gaussian, color='red',alpha=0.5)
ax2.stem(signal_shortx[:,-1], signal_shorty[:,-1], basefmt=' ', use_line_collection=True)
ax2.set_xticks(signal_shortx[0::5,-1])
ax2.axvline(original_position[-1], color='red', linestyle="--", label='Center of gravity')
plt.legend(fontsize=12)
plt.grid(True)

ax3 = plt.subplot(224)
ax3.set_title('Fast Fourier Transform', fontweight="bold", fontsize=14)
ax3.plot(freq_signal, abs(fft_signal))
ax3.set_xlabel('Frequency [Hz]', fontsize=12)
ax3.set_ylabel('PSD 'r'$[{μV^2}/Hz]$', fontsize=12)
plt.grid(True)

# Plotting noise
plt.figure(2)
ax1 = plt.subplot(211)
ax1.plot(time,noise, color='blue', alpha=0.5, label='Original noise')
ax1.plot(time, filtered_noise, color='orange', label='Noise with %dkHz lowpass filter' %(cutoff*1e-3))
ax1.set_title('Noise in %d frames with 6.25μs integration time' %(int(n)), fontweight="bold", fontsize=14)
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Signal [μV]', fontsize=12)
plt.grid(True)
ax1.legend()

ax2 = plt.subplot(223)
ax2.set_title('Fast Fourier Transform of noise', fontweight="bold", fontsize=14)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.plot(freq_noise, abs(fft_noise))
ax2.set_xlabel('Frequency [Hz]', fontsize=12)
ax2.set_ylabel('PSD 'r'$[{μV^2}/Hz]$', fontsize=12)
plt.grid(True)

ax3 = plt.subplot(224)
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.plot(freq_filtered_noise, abs(fft_filtered_noise))
ax3.set_title('Fast Fourier Transform of filtered noise', fontweight="bold", fontsize=14)
ax3.set_xlabel('Frequency [Hz]', fontsize=12)
ax3.set_ylabel('PSD 'r'$[{μV^2}/Hz]$', fontsize=12)
plt.grid(True)

# --------------------------
# Plotting signal plus noise
plt.figure(3)
ax1 = plt.subplot(311)
ax1.plot(time,signal, color='red')
ax1.set_title('Generated ion pulses with different amplitudes', fontweight="bold", fontsize=14)
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Signal [μV]', fontsize=12)
plt.grid(True)

ax2 = plt.subplot(312)
ax2.set_title('Signals plus noise with subtracting the mean value of the sum', fontweight="bold", fontsize=14)
ax2.plot(time,signal_noise, color='blue')
ax2.set_xlabel('Time [s]', fontsize=12)
ax2.set_ylabel('Signal [μV]', fontsize=12)
plt.grid(True)

ax3 = plt.subplot(313)
ax3.set_title('Signal plus Noise with the %dkHz lowpass filter' %(cutoff*1e-3), fontweight="bold", fontsize=14)
ax3.plot(time,filtered_signal_noise, color='orange')
ax3.set_xlabel('Time [s]', fontsize=12)
ax3.set_ylabel('Signal [μV]', fontsize=12)
plt.grid(True)

# --------------------------
# Plotting peak detection
plt.figure(4)
ax1 = plt.subplot(211)
ax1.plot(time, filtered_signal_noise, color='orange')
ax1.axhline(std, xmin=0, xmax=num_loops, color='red', linestyle="--", label='')
ax1.axhline(-std, xmin=0, xmax=num_loops, color='red', linestyle="--", label='3*$\sigma$ after \'sigma clip\'')
ax1.set_title('Signal plus Noise filtered with a %dkHz lowpass filter' %(cutoff*1e-3), fontweight="bold", fontsize=14)
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Signal [μV]', fontsize=12)
ax1.legend(fontsize=12)
plt.grid(True)

ax2 = plt.subplot(223)
ax2.set_title('Detected peak', fontweight="bold", fontsize=14)
ax2.plot(time,detected_peak, color='red')
ax2.set_xlabel('Time [s]', fontsize=12)
ax2.set_ylabel('Signal [μV]', fontsize=12)
plt.grid(True)

ax3 = plt.subplot(224)
ax3.set_title('Detected peak vs. Original pulse (last peak)', fontweight="bold", fontsize=14)
ax3.set_xlabel('Time [s]', fontsize=12)
ax3.set_ylabel('Signal [μV]', fontsize=12)
ax3.plot(np.linspace(time[position_of_signal[-1] - 1]-shift/10*samplingInterval, time[position_of_signal[-1] + len(signal_shortx)]-shift/10*samplingInterval, points_on_gaussian), signal_gaussian, color='red',alpha=0.5)
ax3.stem(signal_shortx[:,-1], signal_shorty[:,-1], basefmt=' ', use_line_collection=True)
ax3.plot(time[position_of_signal[-1] - 2:position_of_signal[-1] + len(signal_shortx)+2], detected_peak[position_of_signal[-1] - 2:position_of_signal[-1] + len(signal_shortx)+2],'o', color='red')
ax3.axvline(original_position[-1], color='red', linestyle="--", label='Center of gravity')
ax3.set_xticks(signal_shortx[0::5,-1])
plt.grid(True)

plt.figure(5)
ax3 = plt.subplot(111)
ax3.set_title('Detected vs. Original pulse', fontweight="bold", fontsize=14)
ax3.set_xlabel('Time [s]', fontsize=12)
ax3.set_ylabel('Signal [μV]', fontsize=12)
ax3.stem(signal_shortx[:,-1], signal_shorty[:,-1], basefmt=' ', use_line_collection=True, label='Original peaks')
#ax3.stem(signal_shortx[:,-2], signal_shorty[:,-2], basefmt=' ', use_line_collection=True)
ax3.plot(np.linspace(time[position_of_signal[-1] - 1]-shift/10*samplingInterval, time[position_of_signal[-1] + len(signal_shortx)]-shift/10*samplingInterval, points_on_gaussian), signal_gaussian, color='blue',alpha=0.3)
#ax3.plot(np.linspace(time[position_of_signal[-2] - 1]-shift/10*samplingInterval, time[position_of_signal[-2] + len(signal_shortx)]-shift/10*samplingInterval, points_on_gaussian), signal_gaussian, color='blue',alpha=0.3)

ax3.plot(time[first - 10:first + 20], filtered_signal_noise[first - 10:first + 20], 'o', color='black', label='Filtered Noise+Signal points')
ax3.plot(detected_x,detected_y,'o',color='red', label='Detected peaks')
x = np.linspace(time[first - 10], time[first + 20], 1000)
ax3.plot(x,func.gaussian_fit(x,*popt), '--', label='Fitted curve')
ax3.axvline(popt[1], color='red', linestyle="-.", label='Detected position', alpha=0.6)
ax3.axvline(original_position[-1], color='green', linestyle="-.", label='Original position', alpha=0.6)
plt.legend(fontsize=12)
plt.grid(True)

plt.tight_layout()
plt.subplots_adjust(left=0.13, bottom=0.13, right=0.91, top=0.92, wspace=0.34, hspace=0.56)
plt.show()
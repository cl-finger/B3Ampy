'''
Use scipy to design butterworth bandpass and lowpass filter and filter data
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
last modified: 2023-04-21
'''

from scipy.signal import butter,sosfiltfilt

def butter_bandpass(lowcut,highcut,fs,order):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	sos = butter(order, [low, high], analog=False, btype='band', output='sos')
	return sos
def butter_lowpass(highcut,fs,order):
	nyq = 0.5 * fs
	high = highcut / nyq
	sos = butter(order, high, analog=False, btype='lowpass', output='sos')
	return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order):
        sos = butter_bandpass(lowcut, highcut, fs, order)
        y = sosfiltfilt(sos, data)
        return y
def butter_lowpass_filter(data, highcut, fs, order):
        sos = butter_bandpass( highcut, fs, order)
        y = sosfiltfilt(sos, data)
        return y

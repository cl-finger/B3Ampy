'''
Functions for pre-processing data
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-21
'''
import numpy as np
from b3am import bp_filter
from scipy import signal

def downsampling(params,data):
	# pre-allocate new data array for downsampling
	step = int(params.sampling_rate/params.sampling_rate_new)
	print('Pre-processing all '+str(params.nstations)+' stations now.')
	if data.ndim == 2: #1C/DAS data: shape (nstations, nsamples), no component axis
		data_ds = np.zeros((data.shape[0],int(np.ceil(data.shape[1]/step))))
		for rr in range(params.nstations): #loop over all stations
			if np.sum(np.abs(data[rr,:]))>0: #exclude stations with no data at this day
				#demean, detrend, taper and bandpass filter
				data[rr,:] = data[rr,:]-np.mean(data[rr,:])
				data[rr,:] = data[rr,:]*signal.windows.tukey(data.shape[1],alpha=0.2)
				data[rr,:] = bp_filter.butter_lowpass_filter(data[rr,:], 0.9*params.sampling_rate_new, params.sampling_rate, 4)
				data[rr,:] = data[rr,:]*signal.windows.tukey(data.shape[1],alpha=0.2)
				data_ds[rr,:] = data[rr,::step]
		params.sampling_rate = params.sampling_rate_new
		return params,data_ds
	data_ds = np.zeros((data.shape[0],data.shape[1],int(np.ceil(data.shape[2]/step))))
	for rr in range(params.nstations): #loop over all stations
		if np.sum(np.abs(data[:,rr,:]))>0: #exclude stations with no data at this day
			#demean, detrend, taper and bandpass filter
			for ii in range(3): #loop over all 3 components
				data[ii,rr,:] = data[ii,rr,:]-np.mean(data[ii,rr,:])
				data[ii,rr,:] = data[ii,rr,:]*signal.windows.tukey(data.shape[2],alpha=0.2)
				data[ii,rr,:] = bp_filter.butter_lowpass_filter(data[ii,rr,:], 0.9*params.sampling_rate_new, params.sampling_rate, 4)
				data[ii,rr,:] = data[ii,rr,:]*signal.windows.tukey(data.shape[2],alpha=0.2)
				data_ds[ii,rr,:] = data[ii,rr,::step]
	params.sampling_rate = params.sampling_rate_new
	return params,data_ds

# def 1bit_norm
def onebit_norm(params, data):
        data_1bit = np.zeros_like(data)
        for rr in range(params.nstations): # Loop over all stations
                if np.sum(np.abs(data[:, rr, :])) > 0: # Exclude stations with no data at this day
                        data_1bit[:, rr, :] = np.sign(data[:, rr, :]) # Apply one-bit normalization for all components and at all time
        return params, data_1bit
# def time_norm

# def spectral_white

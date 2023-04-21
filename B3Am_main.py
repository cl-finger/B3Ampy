'''
B3Am_main.py
---
B3Ampy: Three-component beamformer for ambient noise analysis using python

This is the main script to be executed. Nothing should be changed in here.
It is recommended to run B3am_prep.py first to check all parameters and resolutions
B3Am_post.py will create output plots and pick dispersion curves

Required Input:
- set all parameters in params.py
- additional functions required (all in folder b3am/)

Output:
- (optional) Beam power plots for each time window and frequency
- text files for each frequency containing all beam power picks found 

-----------------------------------------------
Claudia Finger
claudia.finger@ieg.fraunhofer.de
May 2022
last modified 2023-04-21
-----------------------------------------------
Matlab Version created by Katrin Löer
https://
-----------------------------------------------
based on FK3C_FT.m and FK3C_FK.m by Nima Riahi
https://github.com/nimariahi/fk3c
'''
##################
### start time ###
##################
from datetime import datetime
import time
start_time = time.time()
print("--- starting now! It is %s ---" % (datetime.fromtimestamp(start_time)))

########################
### import libraries ###
########################
import numpy as np
from scipy import signal
import math
from multiprocessing import Pool

from b3am import compute_mode_vector as cmv
from b3am import calculate_and_save_beam_power_response as cbp
from b3am import pre_processing as pp
from b3am import FK_resolution as fkr

###############################################################################
### import parameters from params.py and station locations from stationfile ###
###############################################################################
import params
print('I found the following settings in params.py:')
for key, val in params.__dict__.items():
	if not key[0] == '_' and not key == 'os' and not key == 'date':
		print('\t'+str(key)+' = '+str(val))

# load station location file
stationfile = np.loadtxt(params.stationfile,skiprows=1,delimiter=',',dtype=str)
params.nstations = stationfile.shape[0]
print('Number of stations to analyse: '+str(params.nstations))

coords = np.zeros((params.nstations,2))
station_names = np.zeros(params.nstations,dtype='U4')
for ii in range(params.nstations):
	coords[ii,0] = stationfile[ii][1] 		# x, East
	coords[ii,1] = stationfile[ii][2] 		# y, North
	station_names[ii] = stationfile[ii][0]	# station names

########################################
### load ambient noise waveform data ###
########################################
data = np.load(params.indir+params.indata)
print('Loaded 3C waveform data with size '+str(data.shape))

######################################
### perform pre-processing of data ###
######################################
if params.want_down_sampling:
	params,data = pp.downsampling(params,data)
if params.want_onebit:
	#NOT IMPLEMENTED
	test = 0
if params.want_specwhite:
	#NOT IMPLEMENTED
	test = 0	

# rearrange data so that components and stations are one axis. Convention for component ordering: x,y,z!
data_x = np.zeros((params.nstations,data.shape[2]))
data_y = np.zeros((params.nstations,data.shape[2]))
data_z = np.zeros((params.nstations,data.shape[2]))
for ii in range(params.nstations):
	data_x[ii,:] = data[0,ii,:]
	data_y[ii,:] = data[1,ii,:]
	data_z[ii,:] = data[2,ii,:]
data_final = np.concatenate((data_x,data_y,data_z),axis=0)
print('Reorganised data for all '+str(params.nstations)+' stations. Data has new shape ',data_final.shape)

###########################################
### define resolution for FK parameters ###
###########################################
params = fkr.FK_resolution(params,coords)

###########################################################
### calculate STFFT (short time fast fourier transform) ###
###########################################################
longest_window = 4/params.fmin*params.sampling_rate
params.nwin = int(np.power(2,np.ceil(np.log(longest_window)/np.log(2)))) 
length_each_ft = params.sampling_rate/params.fstep
params.nfft = int(np.power(2,np.ceil(np.log(length_each_ft)/np.log(2))))
if params.nfft < params.nwin:
	params.nfft = params.nwin
print('Length of each time window is '+str(params.nwin)+' points = '+str(params.nwin/params.sampling_rate)+' s with 50% overlap')
print('Frequency resolution is '+str(params.sampling_rate/params.nfft)+' Hz')
for ii in range(params.nstations):
	f_,t,DFTE = signal.stft(np.squeeze(data_final[ii,:]),fs = params.sampling_rate,nperseg=params.nwin,nfft=params.nfft)   #default: noverlap = 0.5*nwin
	f_,t,DFTN = signal.stft(np.squeeze(data_final[ii+params.nstations,:]),fs = params.sampling_rate,nperseg=params.nwin,nfft=params.nfft)
	f_,t,DFTZ = signal.stft(np.squeeze(data_final[ii+2*params.nstations,:]),fs = params.sampling_rate,nperseg=params.nwin,nfft=params.nfft)
	if ii ==0:
		#create mask for desired frequencies
		f_mask1 = f_>=params.fmin
		f_mask2 = f_<=params.fmax
		f_mask = f_mask1*f_mask2
		params.f = f_[f_mask]
		DFTES = np.zeros((params.f.size,t.size,params.nstations),dtype=complex)
		DFTNS = np.zeros((params.f.size,t.size,params.nstations),dtype=complex)
		DFTZS = np.zeros((params.f.size,t.size,params.nstations),dtype=complex)
	# Store DFT data into spectrogram containers    
	DFTES[:,:,ii] = DFTE[f_mask,:]
	DFTNS[:,:,ii] = DFTN[f_mask,:]
	DFTZS[:,:,ii] = DFTZ[f_mask,:]
print('Computed spectrograms and stored them in container for each component with shape '+str(DFTES.shape))

###########################
### perform FK analysis ###
###########################
# compute mode vectors for all polarisation states #
ka,params = cmv.compute_mode_vector(params,coords)

# define maximum number of  maxima to be found in beam power
if params.want_multi_peak:
	params.no_maxima = params.no_maxima
else:
	params.no_maxima = 1

# calculate beam power response and save in list file
if params.want_parallel:
	f_pool = np.copy(params.f)
	def worker(f_pool):
    		results = cbp.calculate_results(DFTES,DFTNS,DFTZS,params,f_pool,ka)
    		return results
	pool = Pool(processes = len(f_pool))
	results = pool.map(worker, f_pool)    
	pool.close()
else:
	for ff in params.f:
		cbp.calculate_results(DFTES,DFTNS,DFTZS,params,ff,ka)

#####################################
### write used parameters to file ###
#####################################
fid = open(params.outdir+'params_out.txt','w') #will overwrite existing file
for key, val in params.__dict__.items():
	if not key[0] == '_' and not key == 'os' and not key == 'date':
		#print(key, val)
		tmp = str(key)+' = '+str(val)+' \n' 
		fid.write(tmp)
fid.close()

#################
### endy bits ###
#################
print("--- Finished after %s seconds ---" % (time.time()-start_time))

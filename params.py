'''
Define all parameters and settings for B3Am.
This script will be read by B3Am_main.py.
---------------------------------------
Claudia Finger
claudia.finger@ieg.fraunhofer.de
Mar 2023
'''
import os
###############
### General ###
###############
indir = './test_data/'				# path to 3C waveform data in .npy format 
indata = '3C_waveform_data_syn.npy'		# filename where 3C waveform data are stored (as .npy)
sampling_rate = 200				# Hz, sampling rate of data

stationfile = indir+'station_list_xy_cut.txt' 	# path and filename to stationfile (station coordinates in cartesian coordinates in m)

outdir = indir 					# path to folder where all output is saved

want_beampowerplots = True 			# if True: generate 1 figure per frequency and time window of beam power response

######################
### pre-processing ###
######################
want_specwhite = False				# use spectral whitening: NOT IMPLEMENTED
want_onebit = False				# use one bit time normalisation: NOT IMPLEMENTED
want_down_sampling = False			# if True: downsample data to sampling_rate in next line
if want_down_sampling:
	sampling_rate_new = 100			# Hz, desired sampling rate
else:
	sampling_rate_new = sampling_rate

#########################
### Fourier Transform ###
#########################
want_custom_frange = False 			# True: Use freuqency limits for Fourier Transform from next lines, False: calculate frequency range from station spacing
FT_fmin = 0.1					# minimum frequency (Hz) to save Fourier Transform
FT_fmax = 0.8					# maximum frequency (Hz) to save Fourier Transform
FT_fstep = 0.1 					# largest allowed frequency resolution (must be larger than fmin/4 or will be automatically corrected)

# Defaults for short-time fast-fourier transformation:
# time window length is calculated as next power of 2 of tw_factor/fmin
tw_factor = 10 # should be at least 4
 
# time windows have 50% overlap
######################
### FK calculation ###
######################
k_res = 201 					# number of points for wavenumber grid
custom_klimits = False				# False: use wavenumber limits calculated from station spacing, True: use wavenumber limits defined here (e.g. using array_response_function.py)
kmin = 0					# m
kmax = 0					# m
custom_azi_res = False				# False: use default resolution of 5° for backazimuth, True: use backazimuth resolution from next line
azi_res = 5					# degree
custom_dip_res = False				# False: use default resolution of 5° for dip, True: use dip resolution from next line
dip_res = 5					# in degree
custom_ell_res = False				# False: use default ellipticity resolution of 0.05, True: use ellipticity resolution from next line
ell_res = 0.05					# 

#########################################
### find peaks in beam power response ###
#########################################
min_beam_threshold = 0.7 			# range 0 to 1: threshold for extrema found in beam power response  
want_multi_peak = False 				# False: only look for one maximum in each f-t window, True: look for multiple peaks in each beam power response with limit in next line, 
no_maxima = 5					# maximum number of maxima to be found, e.g. try 5
min_dist_maxima = 10 				# minimum distance between multiple maxima per beam power plot, in points (multiply by resolution) 

##########################
### B3Ampy calculation ###
##########################
want_fast = True 				# True: use fast direct calculation, False: calculate SDM (requires more memory)

want_parallel = False				# True: use multiprocessing to run in parallel (max. one node)
# CAREFUL! Make sure that required memory is available (number of stations * number of time windows)

######################################################################
### copy this parameter file to output folder for future reference ###
######################################################################
from datetime import date
today=date.today() 		#attach current date to filename
os.system('cp ./params.py '+outdir+'/'+str(today)+'_params.py')


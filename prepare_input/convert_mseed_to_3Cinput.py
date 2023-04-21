'''
convert miniSEED data to .npy structure using obspy
- remove instrument response using NRL
- all other preprocessing is done in main B3Am script

This example uses data from the nodal network in Forge, Utah than can be downloaded from IRIS using FDSN code 8J (2016-2017)

This script uses obspy:
https://docs.obspy.org/
Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., and Wassermann, J. (May/June 2010), ObsPy: A Python Toolbox for Seismology, Seismological Research Letters, 81 (3), 530-533.
------------------------------------------------------
Claudia Finger
claudia.finger@ieg.fraunhofer.de
Mar 2023
'''
import obspy
from obspy import UTCDateTime
import numpy as np
from obspy.clients.nrl import NRL
nrl = NRL()
import matplotlib as mpl # for plotting without X server
mpl.use('Agg')
import matplotlib.pyplot as plt
font = {'size':22}
plt.rc('font',**font)
import os


indir = '' #path to mseed data 
outdir = '' #path where 3C waveform data should be stored
data_network_code = '8J'
data_year = '2016'
data_day = '20161220'           # YYYYMMDD, day to analyse
data_start_time = '00:00:00.000'        # HH:MM:SS.sss
data_end_time = '01:59:59.999'  # HH:MM:SS.sss
freq_min = 0.01 # for bandpass filtering
freq_max = 49 # Nyquist Frequency
sampling_rate = 100

path_to_stationfile = outdir+'station_list_xy.txt'

def julian_day(date):
    """
    Converts a date in format YYYYMMDD to the julian day.
    """
    date = str(date)
    if len(date) != 8:
        print("Date must be in the following format: YYYYMMDD")
        return -1
    else:
        days  = int(date[6:])
        month = int(date[4:6])
        year  = int(date[:4])
        if year % 4 != 0:
            month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        elif year % 4 == 0:
            month_days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        for i in range(month-1):
            days += month_days[i]
        return '{:03d}'.format(days)


# station names from stationfile 
stationfile = np.loadtxt(path_to_stationfile,skiprows=1,delimiter=',',dtype=str)
nstations = stationfile.shape[0]
print('Number of stations to analyse: '+str(nstations))
station = np.zeros(nstations,dtype='U4')
for ii in range(nstations):
	station[ii] = stationfile[ii][0]  # station names

starttime =  UTCDateTime(data_year+'-'+data_day[4:6]+'-'+data_day[6:]+'T'+data_start_time)
endtime =  UTCDateTime(data_year+'-'+data_day[4:6]+'-'+data_day[6:]+'T'+data_end_time)
print(starttime,endtime)
data_seis = np.zeros((3,nstations,int((endtime-starttime)*sampling_rate)+3)) #pre allocate for one whole day
for sta in range(nstations):
	print('Reading station ',station[sta])
	comp_order = ('E','N','Z')
	if not os.path.exists(outdir+'PICS_preproc'):
		os.mkdir(outdir+'PICS_preproc')
	fig,axs = plt.subplots(3,1,figsize=(20,15))	
	full_path = indir +data_network_code+'_' +str(station[sta])+'*'# Path for data
	print(full_path)
	try:
		st1 = obspy.read(full_path,fmt='MSEED',starttime=starttime,endtime=endtime)
	except:
		st1 = []
	for comp in range(3): # loop over all three components
		st = st1.select(channel='*'+comp_order[comp])
		print(st)
		if np.size(st)>0:
			print('Found data.')
			st.detrend('linear')
			st.detrend('demean')
			st.filter('bandpass',freqmin= freq_min, freqmax = freq_max, zerophase=True)

			#remove instrument response
			response = nrl.get_response( sensor_keys=['Magseis Fairfield','Generation 2','5 Hz'],datalogger_keys = ['Magseis Fairfield','Zland 1C or 3C','0 dB (1)','250','Linear Phase','Off']) #datalogger keys are unknown
			
			response_dict = response.get_paz().__dict__
			paz_dict = {"gain": response_dict["normalization_factor"],
				"poles": response_dict['_poles'],
				"sensitivity": response.instrument_sensitivity.value ,
				"zeros": response_dict['_zeros']}
			st.simulate(paz_remove=paz_dict)
			st.detrend('linear')
			st.detrend('demean')
			st.filter('bandpass',freqmin= freq_min, freqmax = freq_max, zerophase=True)
			st.merge(fill_value=0) #if multiple traces in mseed, this concatenates it to one file per day
			st.trim(starttime=starttime,endtime=endtime,pad=True,fill_value=0)
			st.taper(0.1)
			data_seis[comp,sta,:st[0].data.size] = st[0].data
		axs[comp].plot(data_seis[comp,sta,:],'k',lw=2)
	plt.tight_layout()
	fig.savefig(outdir+'PICS_preproc/'+'stations_'+str(station[sta])+'_'+str(data_day)+'_'+str(data_start_time)+'_'+str(data_end_time)+'.png')
	plt.close(fig)
	del st
print(np.max(np.max(data_seis,axis=0),axis=1))
#########################
### save data as .npy ###
#########################
np.save(outdir+'3C_waveform_data_'+str(data_day)+'_'+str(data_start_time)+'_'+str(data_end_time)+'.npy',data_seis,allow_pickle=True)

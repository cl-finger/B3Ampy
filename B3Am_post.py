'''
Analyse and plot results from B3Ampy created with B3Am_main.py

Input:
- list files (results_fXX.txt) created by B3Am_main.py
- parameter file params.py

Output:
- Figures
	- wavefield composition
	- Statistics for each wave type (for data quality check)
	- Histograms for each wave type (wavenumber, backazimuth, ellipticity and dip vs. frequency), including picked curves
- list files with picked paramters (wavenumber, backazimuth, ellipticity, dip) for each wave type 

---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-21
'''
import numpy as np
import matplotlib as mpl # for plotting without X server
mpl.use('Agg')
import matplotlib.pyplot as plt
font = {'size':24}
plt.rc('font',**font)
import math
import glob
from scipy import signal

import params 

filenames = glob.glob(params.outdir+'results_*')

########################
### read all results ###
########################
no_f = len(filenames)
no_tw = np.empty((no_f),dtype=int)
for ii in range(no_f):
	data = np.loadtxt(filenames[ii],skiprows=1,delimiter=',')
	no_tw[ii] = int(data.shape[0])
	if ii == 0:
		results = data
	else:
		results = np.append(results,data,axis=0)
print('total number of results',results.shape)

#results dim: 0: time windows, 1: parameters with:
	#0: wave type
	#1: frequency (Hz) 
	#2: backazimuth (deg)
	#3: dip (deg)
	#4: ellipticity
	#5: rotation y (deg)
	#6: wavenumber (/m)
	#7: velocity (m/s)
	#8: max absolute beam power
	#9: signal-to-noise ratio 

####################################################
### discard of unphysical velocities and low SNR ###
####################################################
#remove unphysical velocities
#vmax = 60000
#mask = results[:,7]<vmax
#results = results[mask,:] 

##############################
### divide into wave types ###
##############################
# P wave
mask_P = results[:,0]==0
results_P = results[mask_P,:]
# SV wave
mask_SV = results[:,0]==1
results_SV = results[mask_SV,:]
# retrograde Rayleigh wave
mask_retro = results[:,0]==3
results_retro = results[mask_retro,:]
# prograde Ralyeigh wave
mask_pro = results[:,0] == 4
results_pro = results[mask_pro,:]
# SH wave
mask_SH = results[:,0] == 2
results_SH = results[mask_SH,:] 

frequencies = np.sort(np.unique(results[:,1]))
no_f = len(frequencies)

#remove SNR>np.mean(SNR), unique for each wavetype
#SNR_mean = 0.5*np.mean(results_retro[:,9])
#mask = results_retro[:,9]>SNR_mean
#results_retro = results_retro[mask,:]
#SNR_mean = 0.5*np.mean(results_pro[:,9])
#mask = results_pro[:,9]>SNR_mean
#results_pro = results_pro[mask,:]

##############################
### wave composition plots ###
##############################
### plot wave type pie chart ###
labels = 'P','SV','SH/Love','retrograde Rayleigh','prograde Rayleigh'
sizes = [results_P.shape[0],results_SV.shape[0],results_SH.shape[0],results_retro.shape[0],results_pro.shape[0]]
fig,axs = plt.subplots(figsize=(20,15))
axs.set_title('Distribution of wave types for \n f = '+str(np.min(results[:,1]))+' Hz to f = '+str(np.max(results[:,1]))+' Hz')
axs.pie(sizes,labels=labels,autopct='%1.1f%%',labeldistance=1.1,pctdistance=0.6)
fig.savefig(params.outdir+'wavetype_pie.png')
 
### plot wave types against frequency:histogram ###
# count per frequeny for each wave type
count_f_rel = np.zeros((no_f,5))
count_f_abs = np.zeros((no_f,5))
for ii in range(no_f):
	idx = np.where(results_P[:,1]==frequencies[ii])[0]
	count_f_rel[ii,0] = len(idx)
	count_f_abs[ii,0] = np.sum(results_P[idx,8])
	idx = np.where(results_SV[:,1]==frequencies[ii])[0]
	count_f_rel[ii,1] = len(idx)
	count_f_abs[ii,1] = np.sum(results_SV[idx,8])
	idx = np.where(results_SH[:,1]==frequencies[ii])[0]
	count_f_rel[ii,2] = len(idx)
	count_f_abs[ii,2] = np.sum(results_SH[idx,8])
	idx = np.where(results_retro[:,1]==frequencies[ii])[0]
	count_f_rel[ii,3] = len(idx)
	count_f_abs[ii,3] = np.sum(results_retro[idx,8])
	idx = np.where(results_pro[:,1]==frequencies[ii])[0]
	count_f_rel[ii,4] = len(idx)
	count_f_abs[ii,4] = np.sum(results_pro[idx,8])

	count_f_rel[ii,:] = count_f_rel[ii,:]/np.sum(count_f_rel[ii,:])*100 #normalise to maximum number of time windows)
	

colors='blue','orange','green','red','magenta'

fig,axs = plt.subplots(figsize=(20,15))
axs.set_title('Wavefield composition: absolut')
y_offset = np.zeros(no_f)
for ii in range(5):
	axs.bar(x=frequencies,height=count_f_abs[:,ii],width=frequencies[1]-frequencies[0],bottom=y_offset,color=colors[ii],label=labels[ii])
	y_offset = y_offset + count_f_abs[:,ii]
axs.set_ylabel('absolute amplitude')
axs.set_xlabel('frequency (Hz)')
axs.legend()
fig.savefig(params.outdir+'wave_types_vs_f_hist.png')

fig,axs = plt.subplots(figsize=(20,15))
axs.set_title('Wave type distribution against frequency')
y_offset = np.zeros(no_f)
for ii in range(5):
	axs.bar(x=frequencies,height=count_f_rel[:,ii],width=frequencies[1]-frequencies[0],bottom=y_offset,color=colors[ii],label=labels[ii])
	y_offset = y_offset + count_f_rel[:,ii]
axs.set_ylabel('time window count (%)')
axs.set_xlabel('frequency (Hz)')
axs.legend()
fig.savefig(params.outdir+'wave_types_vs_f_hist_normalised.png')
 
### wave type againt frequency: line plot ###
fig,axs = plt.subplots(figsize=(20,15))
for ii in range(5):
	axs.plot(frequencies,count_f_rel[:,ii],color=colors[ii],label=labels[ii],lw=2,ls='-',marker='*',ms=20)
axs.legend()
axs.set_xlabel('frequency (Hz)')
axs.set_ylabel('time window count (%)')
fig.savefig(params.outdir+'wave_types_vs_f_line.png')

############################################################
### load resolutions that were used to calculate results ###
############################################################
ell = np.load(params.outdir+'ell.npy')
delta_ell = 0.5*(ell[1]-ell[0])
bins_ell = np.linspace(ell[0]-delta_ell,ell[-1]+delta_ell,ell.size+1,endpoint=True)

kgrid = np.load(params.outdir+'kgrid.npy')
delta_k = 0.5*(kgrid[1]-kgrid[0])
bins_k = np.linspace(kgrid[0]-delta_k,kgrid[-1]+delta_k,kgrid.size+1,endpoint=True)

dip = np.load(params.outdir+'dip.npy')*180/math.pi
delta_dip = 0.5*(dip[1]-dip[0])
bins_dip = np.linspace(dip[0]-delta_dip,dip[-1]+delta_dip,dip.size+1,endpoint=True)

kth = np.load(params.outdir+'kth.npy')*180/math.pi
delta_th = 0.5*(kth[1]-kth[0])
bins_th = np.linspace(kth[0]-delta_th,kth[-1]+delta_th,kth.size+1,endpoint=True)

delta_f = 0.5*(frequencies[1]-frequencies[0])
bins_frequencies = np.linspace(frequencies[0]-delta_f,frequencies[-1]+delta_f,frequencies.size+1,endpoint=True)

#######################
### plot statistics ###
#######################
def plot_statistics(results,wavetype):
	fig,axs=plt.subplots(8,1,figsize=(20,45))
	axs[0].hist(results[:,1],bins=bins_frequencies,color='k')
	axs[0].set_xlabel('frequency (Hz)')
	axs[1].hist(results[:,2],bins=bins_th,color='k')
	axs[1].set_xlabel('backazimuth (°)')
	axs[2].hist(results[:,3],bins=bins_dip,color='k')
	axs[2].set_xlabel('dip (°)')
	axs[3].hist(results[:,4],bins=bins_ell,color='k')
	axs[3].set_xlabel('ellipticity')
	axs[4].hist(results[:,5],bins=np.arange(-190,190,5),color='k')
	axs[4].set_xlabel('x rotation (°)')
	axs[5].hist(results[:,6],bins=bins_k,color='k')
	axs[5].set_xlabel('wavenumber (/m)')
	axs[6].hist(results[:,8],bins=np.linspace(0,np.max(results[:,8]),200,endpoint=True),color='k')
	axs[6].set_xlabel('absolute beam power')
	axs[7].hist(results[:,9],bins=np.linspace(0,np.max(results[:,9]),200,endpoint=True),color='k')
	axs[7].set_xlabel('SNR')
	plt.tight_layout()
	fig.savefig(params.outdir+str(wavetype)+'_statistics.png')
plot_statistics(results_P,'P')
plot_statistics(results_SV,'SV')
plot_statistics(results_SH,'SH')
plot_statistics(results_retro,'retro')
plot_statistics(results_pro,'pro')

############################
### calculate Histograms ###
############################

# P
H_k_P, yedges_k_P,xedges_k_P = np.histogram2d(results_P[:,1], results_P[:,6],bins=[bins_frequencies,bins_k],weights = results_P[:,9],density = True)
H_dip_P, yedges_dip_P,xedges_dip_P = np.histogram2d(results_P[:,1], results_P[:,3],bins=[bins_frequencies,bins_dip],weights = results_P[:,9],density=True)
H_th_P, yedges_th_P,xedges_th_P = np.histogram2d(results_P[:,1], results_P[:,2],bins=[bins_frequencies,bins_th],weights = results_P[:,9],density=True)
# SV
H_k_SV, yedges_k_SV,xedges_k_SV = np.histogram2d(results_SV[:,1], results_SV[:,6],bins=[bins_frequencies,bins_k],weights = results_SV[:,9],density = True)
H_dip_SV, yedges_dip_SV,xedges_dip_SV = np.histogram2d(results_SV[:,1], results_SV[:,3],bins=[bins_frequencies,bins_dip],weights = results_SV[:,9],density=True)
H_th_SV, yedges_th_SV,xedges_th_SV = np.histogram2d(results_SV[:,1], results_SV[:,2],bins=[bins_frequencies,bins_th],weights = results_SV[:,9],density=True)
# SH
H_k_SH, yedges_k_SH,xedges_k_SH = np.histogram2d(results_SH[:,1], results_SH[:,6],bins=[bins_frequencies,bins_k],weights = results_SH[:,9],density = True)
H_th_SH, yedges_th_SH,xedges_th_SH = np.histogram2d(results_SH[:,1], results_SH[:,2],bins=[bins_frequencies,bins_th],weights = results_SH[:,9],density=True)
# retro
H_k_retro, yedges_k_retro,xedges_k_retro = np.histogram2d(results_retro[:,1], results_retro[:,6],bins=[bins_frequencies,bins_k],weights = results_retro[:,9],density = True)
H_ell_retro, yedges_ell_retro,xedges_ell_retro = np.histogram2d(results_retro[:,1], results_retro[:,4],bins=[bins_frequencies,bins_ell],weights = results_retro[:,9],density=True)
H_th_retro, yedges_th_retro,xedges_th_retro = np.histogram2d(results_retro[:,1], results_retro[:,2],bins=[bins_frequencies,bins_th],weights = results_retro[:,9],density=True)
# pro
H_k_pro, yedges_k_pro,xedges_k_pro = np.histogram2d(results_pro[:,1], results_pro[:,6],bins=[bins_frequencies,bins_k],weights = results_pro[:,9],density = True)
H_ell_pro, yedges_ell_pro,xedges_ell_pro = np.histogram2d(results_pro[:,1], results_pro[:,4],bins=[bins_frequencies,bins_ell],weights = results_pro[:,9],density=True)
H_th_pro, yedges_th_pro,xedges_th_pro = np.histogram2d(results_pro[:,1], results_pro[:,2],bins=[bins_frequencies,bins_th],weights = results_pro[:,9],density=True)

##########################################
### pick local maxima from histograms  ###
##########################################
def get_peak_of_hist(bins,results_1dhist):
	#find local peaks
	peaks,properties = signal.find_peaks(results_1dhist,height=0.7*np.nanmax(results_1dhist),distance = 0.1*bins.size) #distance in indexes
	vel_max = bins[peaks] #in units
	widths = signal.peak_widths(results_1dhist, peaks, rel_height=0.5)
	vel_max_widths = (bins[1]-bins[0])*widths[0] #convert to real units, assumes linear sampling
	return vel_max,vel_max_widths

# P
k_pick_P = np.zeros((H_k_P.shape[0],10))
k_pick_P[:,:] = np.nan
k_width_P = np.zeros((H_k_P.shape[0],10))
k_width_P[:,:] = np.nan
for ii in range(H_k_P.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_k_P,H_k_P[ii,:])
	k_pick_P[ii,:peaks.size] = peaks+delta_k
	k_width_P[ii,:peaks.size] = widths
th_pick_P = np.zeros((H_th_P.shape[0],10))
th_pick_P[:,:] = np.nan
th_width_P = np.zeros((H_th_P.shape[0],10))
th_width_P[:,:] = np.nan
for ii in range(H_th_P.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_th_P,H_th_P[ii,:])
	th_pick_P[ii,:peaks.size] = peaks+delta_th
	th_width_P[ii,:peaks.size] = widths
dip_pick_P = np.zeros((H_dip_P.shape[0],10))
dip_pick_P[:,:] = np.nan
dip_width_P = np.zeros((H_dip_P.shape[0],10))
dip_width_P[:,:] = np.nan
for ii in range(H_dip_P.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_dip_P,H_dip_P[ii,:])
	dip_pick_P[ii,:peaks.size] = peaks+delta_dip
	dip_width_P[ii,:peaks.size] = widths
# SV
k_pick_SV = np.zeros((H_k_SV.shape[0],10))
k_pick_SV[:,:] = np.nan
k_width_SV = np.zeros((H_k_SV.shape[0],10))
k_width_SV[:,:] = np.nan
for ii in range(H_k_SV.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_k_SV,H_k_SV[ii,:])
	k_pick_SV[ii,:peaks.size] = peaks+delta_k
	k_width_SV[ii,:peaks.size] = widths
th_pick_SV = np.zeros((H_th_SV.shape[0],10))
th_pick_SV[:,:] = np.nan
th_width_SV = np.zeros((H_th_SV.shape[0],10))
th_width_SV[:,:] = np.nan
for ii in range(H_th_SV.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_th_SV,H_th_SV[ii,:])
	th_pick_SV[ii,:peaks.size] = peaks+delta_th
	th_width_SV[ii,:peaks.size] = widths
dip_pick_SV = np.zeros((H_th_SV.shape[0],10))
dip_pick_SV[:,:] = np.nan
dip_width_SV = np.zeros((H_th_SV.shape[0],10))
dip_width_SV[:,:] = np.nan
for ii in range(H_dip_SV.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_dip_SV,H_dip_SV[ii,:])
	dip_pick_SV[ii,:peaks.size] = peaks+delta_dip
	dip_width_SV[ii,:peaks.size] = widths
# SH
k_pick_SH = np.zeros((H_k_SH.shape[0],10))
k_pick_SH[:,:] = np.nan
k_width_SH = np.zeros((H_k_SH.shape[0],10))
k_width_SH[:,:] = np.nan
for ii in range(H_k_SH.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_k_SH,H_k_SH[ii,:])
	k_pick_SH[ii,:peaks.size] = peaks+delta_k
	k_width_SH[ii,:peaks.size] = widths
th_pick_SH = np.zeros((H_th_SH.shape[0],10))
th_pick_SH[:,:] = np.nan
th_width_SH = np.zeros((H_th_SH.shape[0],10))
th_width_SH[:,:] = np.nan
for ii in range(H_th_SH.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_th_SH,H_th_SH[ii,:])
	th_pick_SH[ii,:peaks.size] = peaks+delta_th
	th_width_SH[ii,:peaks.size] = widths
# retro
k_pick_retro = np.zeros((H_k_retro.shape[0],10))
k_pick_retro[:,:] = np.nan
k_width_retro = np.zeros((H_k_retro.shape[0],10))
k_width_retro[:,:,] = np.nan
for ii in range(H_k_retro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_k_retro,H_k_retro[ii,:])
	k_pick_retro[ii,:peaks.size] = peaks+delta_k
	k_width_retro[ii,:peaks.size] = widths
th_pick_retro = np.zeros((H_th_retro.shape[0],10))
th_pick_retro[:,:] = np.nan
th_width_retro = np.zeros((H_th_retro.shape[0],10))
th_width_retro[:,:,] = np.nan
for ii in range(H_th_retro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_th_retro,H_th_retro[ii,:])
	th_pick_retro[ii,:peaks.size] = peaks+delta_th
	th_width_retro[ii,:peaks.size] = widths
ell_pick_retro = np.zeros((H_ell_retro.shape[0],10))
ell_pick_retro[:,:] = np.nan
ell_width_retro = np.zeros((H_ell_retro.shape[0],10))
ell_width_retro[:,:,] = np.nan
for ii in range(H_ell_retro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_ell_retro,H_ell_retro[ii,:])
	ell_pick_retro[ii,:peaks.size] = peaks+delta_ell
	ell_width_retro[ii,:peaks.size] = widths
# pro
k_pick_pro = np.zeros((H_k_pro.shape[0],10))
k_pick_pro[:,:] = np.nan
k_width_pro = np.zeros((H_k_pro.shape[0],10))
k_width_pro[:,:] = np.nan
for ii in range(H_k_pro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_k_pro,H_k_pro[ii,:])
	k_pick_pro[ii,:peaks.size] = peaks+delta_k
	k_width_pro[ii,:peaks.size] = widths
th_pick_pro = np.zeros((H_th_pro.shape[0],10))
th_pick_pro[:,:] = np.nan
th_width_pro = np.zeros((H_th_pro.shape[0],10))
th_width_pro[:,:,] = np.nan
for ii in range(H_th_pro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_th_pro,H_th_pro[ii,:])
	th_pick_pro[ii,:peaks.size] = peaks+delta_th
	th_width_pro[ii,:peaks.size] = widths
ell_pick_pro = np.zeros((H_ell_pro.shape[0],10))
ell_pick_pro[:,:] = np.nan
ell_width_pro = np.zeros((H_ell_pro.shape[0],10))
ell_width_pro[:,:] = np.nan
for ii in range(H_ell_pro.shape[0]):
	peaks,widths = get_peak_of_hist(xedges_ell_pro,H_ell_pro[ii,:])
	ell_pick_pro[ii,:peaks.size] = peaks+delta_ell
	ell_width_pro[ii,:peaks.size] = widths
###################
### plot P wave ###
###################
fig,axs = plt.subplots(3,1,figsize=(15,20))
axs[0].pcolormesh(yedges_k_P, xedges_k_P, H_k_P.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_k_P))
for ii in range(10):
	axs[0].errorbar(frequencies,k_pick_P[:,ii],yerr=k_width_P[:,ii],c='r',lw=2)
axs[0].set_title('P wave')
axs[0].set_ylabel('wavenumber (/m)')
axs[0].set_xlabel('frequency (Hz)')
axs[1].pcolormesh(yedges_th_P, xedges_th_P, H_th_P.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_th_P))
for ii in range(10):
	axs[1].errorbar(frequencies,th_pick_P[:,ii],yerr=th_width_P[:,ii],c='r',lw=2)
axs[1].set_ylabel('azimuth (degree)')
axs[1].set_xlabel('frequency (Hz)')
axs[2].pcolormesh(yedges_dip_P, xedges_dip_P, H_dip_P.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_dip_P))
for ii in range(10):
	axs[2].errorbar(frequencies,dip_pick_P[:,ii],yerr=dip_width_P[:,ii],c='r',lw=2)
axs[2].set_ylabel('dip (degree)')
axs[2].set_xlabel('frequency (Hz)')
axs[2].invert_yaxis()
plt.tight_layout()
fig.savefig(params.outdir+'P_hist.png')
####################
### plot SV wave ###
####################
fig,axs = plt.subplots(3,1,figsize=(15,20))
axs[0].pcolormesh(yedges_k_SV, xedges_k_SV, H_k_SV.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_k_SV))
for ii in range(10):
	axs[0].errorbar(frequencies,k_pick_SV[:,ii],yerr=k_width_SV[:,ii],c='r',lw=2)
axs[0].set_title('SV wave')
axs[0].set_ylabel('wavenumber (/m)')
axs[0].set_xlabel('frequency (Hz)')
axs[1].pcolormesh(yedges_th_SV, xedges_th_SV, H_th_SV.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_th_SV))
for ii in range(10):
	axs[1].errorbar(frequencies,th_pick_SV[:,ii],yerr=th_width_SV[:,ii],c='r',lw=2)
axs[1].set_ylabel('azimuth (degree)')
axs[1].set_xlabel('frequency (Hz)')
axs[2].pcolormesh(yedges_dip_SV, xedges_dip_SV, H_dip_SV.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_dip_SV))
for ii in range(10):
	axs[2].errorbar(frequencies,dip_pick_SV[:,ii],yerr=dip_width_SV[:,ii],c='r',lw=2)
axs[2].set_ylabel('dip (degree)')
axs[2].set_xlabel('frequency (Hz)')
axs[2].invert_yaxis()
plt.tight_layout()
fig.savefig(params.outdir+'SV_hist.png')
####################
### plot SH wave ###
####################
fig,axs = plt.subplots(2,1,figsize=(15,20))
axs[0].pcolormesh(yedges_k_SH, xedges_k_SH, H_k_SH.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_k_SH))
for ii in range(10):
	axs[0].errorbar(frequencies,k_pick_SH[:,ii],yerr=k_width_SH[:,ii],c='r',lw=2)
axs[0].set_title('SH/ Love wave')
axs[0].set_ylabel('wavenumber (/m)')
axs[0].set_xlabel('frequency (Hz)')
axs[1].pcolormesh(yedges_th_SH, xedges_th_SH, H_th_SH.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_th_SH))
for ii in range(10):
	axs[1].errorbar(frequencies,th_pick_SH[:,ii],yerr=th_width_SH[:,ii],c='r',lw=2)
axs[1].set_ylabel('azimuth (degree)')
axs[1].set_xlabel('frequency (Hz)')
plt.tight_layout()
fig.savefig(params.outdir+'SH_hist.png')
#####################################
### plot retrograde Rayleigh wave ###
#####################################
fig,axs = plt.subplots(3,1,figsize=(15,20))
axs[0].pcolormesh(yedges_k_retro, xedges_k_retro, H_k_retro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_k_retro))
for ii in range(10):
	axs[0].errorbar(frequencies,k_pick_retro[:,ii],yerr=k_width_retro[:,ii],c='r',lw=2)
axs[0].set_title('retrograde Rayleigh wave')
axs[0].set_ylabel('wavenumber (/m)')
axs[0].set_xlabel('frequency (Hz)')
axs[1].pcolormesh(yedges_th_retro, xedges_th_retro, H_th_retro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_th_retro))
for ii in range(10):
	axs[1].errorbar(frequencies,th_pick_retro[:,ii],yerr=th_width_retro[:,ii],c='r',lw=2)
axs[1].set_ylabel('azimuth (degree)')
axs[1].set_xlabel('frequency (Hz)')
axs[2].pcolormesh(yedges_ell_retro, xedges_ell_retro, H_ell_retro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_ell_retro))
for ii in range(10):
	axs[2].errorbar(frequencies,ell_pick_retro[:,ii],yerr=ell_width_retro[:,ii],c='r',lw=2)
axs[2].set_ylabel('ellipticity')
axs[2].set_xlabel('frequency (Hz)')
axs[2].invert_yaxis()
plt.tight_layout()
fig.savefig(params.outdir+'retrograde_Rayleigh_hist.png')
##############################
### plot prograde Rayleigh ###
##############################
fig,axs = plt.subplots(3,1,figsize=(15,20))
axs[0].pcolormesh(yedges_k_pro, xedges_k_pro, H_k_pro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_k_pro))
for ii in range(10):
	axs[0].errorbar(frequencies,k_pick_pro[:,ii],yerr=k_width_pro[:,ii],c='r',lw=2)
axs[0].set_title('prograde Rayleigh wave')
axs[0].set_ylabel('wavenumber (/m)')
axs[0].set_xlabel('frequency (Hz)')
axs[1].pcolormesh(yedges_th_pro, xedges_th_pro, H_th_pro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_th_pro))
for ii in range(10):
	axs[1].errorbar(frequencies,th_pick_pro[:,ii],yerr=th_width_pro[:,ii],c='r',lw=2)
axs[1].set_ylabel('azimuth (degree)')
axs[1].set_xlabel('frequency (Hz)')
axs[2].pcolormesh(yedges_ell_pro, xedges_ell_pro, H_ell_pro.T, cmap='gray_r',vmin = 0,vmax = 0.7*np.nanmax(H_ell_pro))
for ii in range(10):
	axs[2].errorbar(frequencies,ell_pick_pro[:,ii],yerr=ell_width_pro[:,ii],c='r',lw=2)
axs[2].set_ylabel('ellipticity')
axs[2].set_xlabel('frequency (Hz)')
axs[2].invert_yaxis()
plt.tight_layout()
fig.savefig(params.outdir+'prograde_Rayleigh_hist.png')

###################################################
### plot polar plot f vs azimuth: all wavetypes ###
###################################################

fig=plt.figure(figsize=(20,15))
fig.suptitle('frequency vs. backazimuth')

ax=fig.add_subplot(231,polar=True)
ax.set_title('P wave')
ax.pcolormesh(np.radians(xedges_th_P),yedges_th_P,H_th_P,cmap='gray_r',shading='auto') #pcolormesh
ax.grid(True)

ax=fig.add_subplot(232,polar=True)
ax.grid(True)
ax.set_title('SV wave')
ax.pcolormesh(np.radians(xedges_th_SV),yedges_th_SV,H_th_SV,cmap='gray_r',shading='auto') #pcolormesh
ax.grid(True)

ax=fig.add_subplot(233,polar=True)
ax.grid(True)
ax.set_title('SH wave')
ax.pcolormesh(np.radians(xedges_th_SH),yedges_th_SH,H_th_SH,cmap='gray_r',shading='auto') #pcolormesh
ax.grid(True)

ax=fig.add_subplot(234,polar=True)
ax.grid(True)
ax.set_title('retrograde Rayleigh wave')
ax.pcolormesh(np.radians(xedges_th_retro),yedges_th_retro,H_th_retro,cmap='gray_r',shading='auto') #pcolormesh
ax.grid(True)

ax=fig.add_subplot(235,polar=True)
ax.grid(True)
ax.set_title('prograde Rayleigh wave')
ax.pcolormesh(np.radians(xedges_th_pro),yedges_th_pro,H_th_pro,cmap='gray_r',shading='auto') #pcolormesh
ax.grid(True)


plt.tight_layout()
fig.savefig(params.outdir+'all_polar.png')

##########################
### save picked things ###
##########################
#P
P_picks = np.zeros((no_f,7))
P_picks[:,0] = frequencies
P_picks[:,1] = k_pick_P[:,0]
P_picks[:,2] = k_width_P[:,0]
P_picks[:,3] = th_pick_P[:,0]
P_picks[:,4] = th_width_P[:,0]
P_picks[:,5] = dip_pick_P[:,0]
P_picks[:,6] = dip_width_P[:,0]
np.savetxt(params.outdir+'P_picks.txt',P_picks,
	header='freq (Hz), wavenumber (/m), wavenumber error (/m), azimuth(deg), azimuth error (deg), dip (deg), dip error (deg) ',
	fmt = '%.4f, %.6f, %.6f, %.1f, %.1f, %.1f, %.1f')
#SV
SV_picks = np.zeros((no_f,7))
SV_picks[:,0] = frequencies
SV_picks[:,1] = k_pick_SV[:,0]
SV_picks[:,2] = k_width_SV[:,0]
SV_picks[:,3] = th_pick_SV[:,0]
SV_picks[:,4] = th_width_SV[:,0]
SV_picks[:,5] = dip_pick_SV[:,0]
SV_picks[:,6] = dip_width_SV[:,0]
np.savetxt(params.outdir+'SV_picks.txt',SV_picks,
	header='freq (Hz), wavenumber (/m), wavenumber error (/m), azimuth(deg), azimuth error (deg), dip (deg), dip error (deg) ',
	fmt = '%.4f, %.6f, %.6f, %.1f, %.1f, %.1f, %.1f')

#SH
SH_picks = np.zeros((no_f,5))
SH_picks[:,0] = frequencies
SH_picks[:,1] = k_pick_SH[:,0]
SH_picks[:,2] = k_width_SH[:,0]
SH_picks[:,3] = th_pick_SH[:,0]
SH_picks[:,4] = th_width_SH[:,0]
np.savetxt(params.outdir+'SH_picks.txt',SH_picks,
	header='freq (Hz), wavenumber (/m), wavenumber error (/m), azimuth(deg), azimuth error (deg)',
	fmt = '%.4f, %.6f, %.6f, %.1f, %.1f')
# retro
retro_picks = np.zeros((no_f,7))
retro_picks[:,0] = frequencies
retro_picks[:,1] = k_pick_retro[:,0]
retro_picks[:,2] = k_width_retro[:,0]
retro_picks[:,3] = th_pick_retro[:,0]
retro_picks[:,4] = th_width_retro[:,0]
retro_picks[:,5] = ell_pick_retro[:,0]
retro_picks[:,6] = ell_width_retro[:,0]
np.savetxt(params.outdir+'retro_picks.txt',retro_picks,
	header='freq (Hz), wavenumber (/m), wavenumber error (/m), azimuth(deg), azimuth error (deg), ell, ell error ',
	fmt = '%.4f, %.6f, %.6f, %.1f, %.1f, %.2f, %.2f')
# pro
pro_picks = np.zeros((no_f,7))
pro_picks[:,0] = frequencies
pro_picks[:,1] = k_pick_pro[:,0]
pro_picks[:,2] = k_width_pro[:,0]
pro_picks[:,3] = th_pick_pro[:,0]
pro_picks[:,4] = th_width_pro[:,0]
pro_picks[:,5] = ell_pick_pro[:,0]
pro_picks[:,6] = ell_width_pro[:,0]
np.savetxt(params.outdir+'pro_picks.txt',pro_picks,
	header='freq (Hz), wavenumber (/m), wavenumber error (/m), azimuth(deg), azimuth error (deg), ell, ell error ',
	fmt = '%.4f, %.6f, %.6f, %.1f, %.1f, %.2f, %.2f')

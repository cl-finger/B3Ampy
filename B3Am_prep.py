'''
B3Am_prep.py
---
Compute and plot array response vector using parameters specified in params.py
---
Claudia Finger 
claudia.finger@ieg.fraunhofer.de
last modified 2023-04-21
'''
# load python libraries
import numpy as np
import math
from scipy import sparse
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
font = {'size':18}
plt.rc('font',**font)

# load B3Ampy functions
from b3am import FK_resolution as fkr
# load parameters
import params

print('Preparing array response vector for '+params.indir+params.indata)
################################
### load station coordinates ###
################################
data = np.loadtxt(params.stationfile,dtype=str,delimiter=',')
params.nstations = data.shape[0]
coords = np.zeros((params.nstations,2))
for ii in range(params.nstations):
	coords[ii,0] = data[ii][1]
	coords[ii,1] = data[ii][2]

###########################################
### calculate distance between stations ###
##########################################
dist = np.zeros((params.nstations,params.nstations))
for ii in range(params.nstations):
    for jj in range(params.nstations):
        dist[ii,jj] = np.sqrt((coords[ii,0]-coords[jj,0])**2+(coords[ii,1]-coords[jj,1])**2)
dist[dist==0] = np.nan
dmax = np.nanmax(dist)
dmin = np.nanmin(dist)

##############################
### get resolutions for FK ###
##############################
params = fkr.FK_resolution(params,coords)
kmax_plot = 2*params.kmax 		# max. wavenumber to plot ARF in
kgrid_plot = np.linspace(0,kmax_plot,2*params.k_res)
#Compute polar grid from angle and wavenumber grids
k = (sparse.kron(kgrid_plot, np.concatenate((np.cos(params.kth),np.sin(params.kth))).reshape((2,params.kth.size)) ) ).toarray()
km = 1/np.sqrt(params.nstations)*np.exp(1j*2*math.pi*np.dot(coords,k))
#sum over all stations
km_sum = np.sum(km,axis=0).reshape(kgrid_plot.size,params.kth.size)
#normalise
km_sum = np.real(np.abs(km_sum/np.nanmax(np.abs(km_sum))))

##################################
### plot array response vector ###
##################################
fig=plt.figure(figsize=(20,15))

ax2=fig.add_subplot(231)
ax2.set_title('station locations')
plt.scatter(coords[:,0],coords[:,1],s=100,c='k')
ax2.set_xlabel('Easting (m)')
ax2.set_ylabel('Northing (m)')
ax2.set_aspect('equal')

ax=fig.add_subplot(232,polar=True)
ax.grid(False)
ax.pcolormesh(params.kth,kgrid_plot,km_sum,cmap='gray_r',shading='auto') #pcolormesh
ax.set_rticks([params.kmin,params.kmax,kmax_plot])
circle = plt.Circle((0, 0), params.kmin, transform=ax.transData._b, edgecolor="red",facecolor='none',lw=2,alpha=0.5,label='k min = 1/(3 aperture)',ls='--')
ax.add_artist(circle)
circle = plt.Circle((0, 0), params.kmax, transform=ax.transData._b, edgecolor="red",facecolor='none',lw=2,alpha=0.5,label='k max = 1/(2 spacing)')
ax.add_artist(circle)


ax4 = fig.add_subplot(233)
ax4.plot((params.fmin,params.fmax),(params.fmin/params.kmax,params.fmax/params.kmax),lw=2,c='k')
ax4.plot((params.fmin,params.fmax),(params.fmin/params.kmin,params.fmax/params.kmin),lw=2,c='k')
ax4.plot((params.fmin,params.fmax),(params.fmin/params.kmin,params.fmin/params.kmin),lw=2,c='r')
ax4.plot((params.fmin,params.fmax),(params.fmax/params.kmax,params.fmax/params.kmax),lw=2,c='r')
ax4.plot((params.fmin,params.fmin),(params.fmax/params.kmax,params.fmin/params.kmin),lw=2,c='r')
ax4.plot((params.fmax,params.fmax),(params.fmax/params.kmax,params.fmin/params.kmin),lw=2,c='r',label='frequency and velocity limits')
ax4.set_ylabel('velocity (m/s)')
ax4.set_xlabel('frequency (Hz)')
ax4.legend(loc='lower left',bbox_to_anchor=(0,1))
ax4.set_aspect((params.fmax/(params.fmax/params.kmin)))

ax2=fig.add_subplot(235)
for ii in range(km_sum.shape[1]):
	plt.plot(kgrid_plot,km_sum[:,ii],'k')
plt.axvline(x=params.kmin,c='r',lw=2,label='wavenumber limits')
plt.axvline(x=params.kmax,c='r',lw=2)
plt.axhline(y=0.5,c='k')
ax2.legend(loc='lower left',bbox_to_anchor=(0,1))
ax2.set_xlabel('wavenumber ($m^{-1}$)')
ax2.set_ylabel('array response')
ax2.text(1.5*np.max(kgrid_plot),0,
	'k min = '+str(np.round(params.kmin,5))+' /m \n'+
	'k max = '+str(np.round(params.kmax,5))+' /m \n'+
	'f min = '+str(params.fmin)+' Hz \n'+
	'f max = '+str(params.fmax)+' Hz \n'+
	'v min = '+str(np.round(params.fmax/params.kmax,1))+' m/s \n'+
	'v max = '+str(np.round(params.fmin/params.kmin,1))+' m/s \n'
	,fontsize=26)
ax2.set_aspect(kmax_plot/1)

ax3=fig.add_subplot(234)
ax3.hist(dist.flatten(),color='k',bins=20)
ax3.set_xlabel('distance between receivers (m)')
ax3.set_ylabel('count')
ax3.set_aspect(4*np.nanmax(dist)/params.nstations**2)

#plt.tight_layout(w_pad=-9,h_pad=-2)
plt.tight_layout(w_pad=-12,h_pad=-2)
fig.savefig(params.outdir+'array_response_and_FK_resolution.png',dpi=300)

'''
For each discrete frequency, calculate and save beam power response
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
last modified: 2023-04-21
'''

import numpy as np
from . import get_wave_type as gwt
from . import calculate_beam_power as cbp
from . import compute_SDM as csdm
from . import extrema as ex 
import math
import os
import matplotlib as mpl # for plotting without X server
mpl.use('Agg')
import matplotlib.pyplot as plt
font = {'size':22}
plt.rc('font',**font)

def calculate_results(DFTES,DFTNS,DFTZS,params,f0,ka):
    #f0 = f[ff]
    ff = np.where(params.f == f0)
    print('Calculating beam power for all time windows for f = '+str(f0)+' Hz')
    data_E_f0 = np.squeeze(DFTES[ff,:,:])
    data_N_f0 = np.squeeze(DFTNS[ff,:,:])
    data_Z_f0 = np.squeeze(DFTZS[ff,:,:])
    data_allcomp = np.concatenate((data_E_f0,data_N_f0,data_Z_f0),axis=1).astype(np.csingle) #shape: 0=number of time windows, 1= 3*nstations, 19.07.24 changed dtype to np.csingle
    # calculate beam power and polarisation matrix
    if params.want_fast:
        P,Q,params = cbp.get_beampower_fast(data_allcomp,ka,params)
    else:
        S = csdm.compute_SDM(f0,data_allcomp)
        P,Q,params = cbp.get_beampower_SDM(S,ka,params)

    results_tmp = np.zeros((params.Nt,params.no_maxima,10)) #initialise array with results from maxima
    for ii in range(params.Nt): #loop over all time windows
        print('Getting wave parameters at maxima for time window '+str(ii+1)+' of '+str(params.Nt))
        #P_norm = ((P[:,:,ii])/np.nanmax(np.abs(P[:,:,ii])))**2 #normalise beam response
        P_norm = np.abs((P[:,:,ii])/np.nanmax(np.abs(P[:,:,ii]))) #normalise beam response
        if params.no_maxima==1:
            ind_max_x,ind_max_y = np.where(P_norm==np.nanmax(P_norm))
            ind_max_x = ind_max_x[0] #only first maximum if more points have maximal value, only failsafe
            ind_max_y = ind_max_y[0]
            if P_norm[int(ind_max_x),int(ind_max_y)] > params.min_beam_threshold:
                ind_kth = params.kth[ind_max_x]
                ind_kgrid = params.kgrid[ind_max_y]
            else:
                ind_kth = np.nan
                ind_kgrid = np.nan
            pID = int(Q[int(ind_max_x),int(ind_max_y),ii])
            #azimuth = polstates[pID,0]
            theta = params.polstates[pID,1]
            ellipticity = params.polstates[pID,2]
            rotation = params.polstates[pID,3]

            #get wave type number for these polarisation parameters: 0(P), 1(SV), 2(SH/Love), 3(retrograde Rayleigh), 4(prograde Rayleigh)
            wave_type = gwt.get_wave_type(ellipticity, rotation)
            jj = 0
            #put results together to write in output list
            results_tmp[ii,jj,:] = (wave_type,                              #wave type
                f0,                                     #frequency (Hz)
                float(ind_kth)*180/math.pi,     #backazimuth (degree)
                theta*180/math.pi,                      #dip (degree)
                ellipticity,                            #ellipticity (-)
                rotation*180/math.pi,                   #rotation around x (degree)
                float(ind_kgrid),               #wavenumber (/m)
                f0/float(ind_kgrid),            #velocity (m/s)
                P[int(ind_max_x),int(ind_max_y),ii],    #absolute beam power
                P_norm[int(ind_max_x),int(ind_max_y)]/np.mean(P_norm[:,:])) #signal-t-noise-ratio of beam power peak

            if ii == 0 and jj == 0:
                results = results_tmp[ii,jj,:].reshape(1,10)
            else:
                results = np.append(results,results_tmp[ii,jj,:].reshape(1,10),axis=0)
        else:
            ind_max_x,ind_max_y = ex.extrema(P_norm,params.min_beam_threshold,params.min_dist_maxima,params.no_maxima)
            ind_kth = np.zeros(ind_max_x.size)
            ind_kgrid = np.zeros(ind_max_x.size)
            for mm in range(ind_max_x.size):
                if P_norm[int(ind_max_x[mm]),int(ind_max_y[mm])] > params.min_beam_threshold:
                    ind_kth[mm] = params.kth[ind_max_x[mm]]
                    ind_kgrid[mm] = params.kgrid[ind_max_y[mm]]
                else:
                    ind_kth[mm] = np.nan
                    ind_kgrid[mm] = np.nan
            # get polaristaion parameters for this k-azimuth pair
            for jj in range(np.array(ind_max_x).size):
                pID = int(Q[int(ind_max_x[jj]),int(ind_max_y[jj]),ii])
                #azimuth = polstates[pID,0]
                theta = params.polstates[pID,1]	
                ellipticity = params.polstates[pID,2]	
                rotation = params.polstates[pID,3]

                #get wave type number for these polarisation parameters: 0(P), 1(SV), 2(SH/Love), 3(retrograde Rayleigh), 4(prograde Rayleigh)
                wave_type = gwt.get_wave_type(ellipticity, rotation)

                #put results together to write in output list
                results_tmp[ii,jj,:] = (wave_type,                              #wave type
                    f0,                                     #frequency (Hz)
                    float(ind_kth[jj])*180/math.pi, #backazimuth (degree)
                    theta*180/math.pi,                      #dip (degree)
                    ellipticity,                            #ellipticity (-)
                    rotation*180/math.pi,                   #rotation around x (degree)
                    float(ind_kgrid[jj]),           #wavenumber (/m)
                    f0/float(ind_kgrid[jj]),                #velocity (m/s)
                    P[int(ind_max_x[jj]),int(ind_max_y[jj]),ii],    #absolute beam power
                    P_norm[int(ind_max_x[jj]),int(ind_max_y[jj])]/np.mean(P_norm[:,:])) #signal-t-noise-ratio of beam power peak

                if ii == 0 and jj == 0:
                    results = results_tmp[ii,jj,:].reshape(1,10)
                else:
                    results = np.append(results,results_tmp[ii,jj,:].reshape(1,10),axis=0)

        if params.want_beampowerplots:
            #create folder for beam power pics	
            if not os.path.exists(params.outdir+'PICS_f'+str(f0)):
                os.mkdir(params.outdir+'PICS_f'+str(f0))
            wave_type_name = ['P','SV','L','retro R','pro R']
            fig,axs=plt.subplots(2,1,figsize=(20,15))
            #axs[0].set_title(wave_type[int(results_tmp[ii,jj,0])]+' wave with SNR = '+str(P_norm[ind_max_x,ind_max_y]/np.mean(P_norm[:,:]))    )
            im = axs[0].imshow(P_norm[:,:].T,extent=(np.min(params.kth*180/math.pi),np.max(params.kth*180/math.pi),np.max(params.kgrid),np.min(params.kgrid)),cmap='Greys')
            plt.colorbar(im,ax=axs[0])
            axs[0].plot(ind_kth*180/math.pi,ind_kgrid,'*r',ms=24)
            #axs[0].text(-50,ind_kgrid-0.00005,'azimuth = '+str(int(np.round(ind_kth*180/math.pi)))+', v = f/k = '+str(np.round(f0/ind_kgrid,2)),color='r')
            axs[0].set_xlabel('azimuth (°)')
            axs[0].set_ylabel('horizontal wavenumber ($m^{-1}$)')
            axs[0].set_aspect((0.5*360/np.abs(np.max(params.kgrid)-np.min(params.kgrid))))

            #custom colormap
            cmap = plt.cm.Greys  # define the colormap
            # extract all colors from the .jet map
            cmaplist = [cmap(i) for i in range(cmap.N)]	
            # create the new map
            cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)	
            # define the bins and normalize
            bounds = (0,params.dip.size,2*params.dip.size,2*params.dip.size+1,2*params.dip.size+1+params.ell.size,2*params.dip.size+1+2*params.ell.size)
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

            im = axs[1].imshow(Q[:,:,ii].T,extent=(np.min(params.kth*180/math.pi),np.max(params.kth*180/math.pi),np.max(params.kgrid),np.min(params.kgrid)),cmap=cmap,norm=norm)
            #cb = plt.colorbar(im,ax=axs[1],ticks = [0.5*dip.size,1.5*dip.size,2*dip.size+0.5,2*dip.size+1+0.5*ell.size,2*dip.size+1+1.5*ell.size])
            #cb.ax.set_yticklabels(['P','SV','L','retro R','pro R'])
            cb = plt.colorbar(im,ax=axs[1],ticks = [0.5*params.dip.size,1.5*params.dip.size,2*params.dip.size+0.5,2*params.dip.size+1+0.5*params.ell.size,2*params.dip.size+1+1.5*params.ell.size])
            cb.ax.set_yticklabels(['P','SV','SH/L','retro R','pro R'])
            axs[1].plot(ind_kth*180/math.pi,ind_kgrid,'*r',ms=24)
            #axs[1].text(-50,ind_kgrid-0.00005,wave_type_name[wave_type]+' wave, dip = '+str(np.round(int(theta*180/math.pi)))+'°, ell = '+str(np.round(int(ellipticity)))+', r    ot x = '+str(np.round(int(rotation*180/math.pi)))+'°',color='r')
            axs[1].set_xlabel('azimuth (°)')
            axs[1].set_ylabel('horizontal wavenumber ($m^{-1}$)')	
            axs[1].set_aspect((0.5*360/np.abs(np.max(params.kgrid)-np.min(params.kgrid))))
            plt.tight_layout()
            fig.savefig(params.outdir+'PICS_f'+str(f0)+'/'+'beampower_f'+str(f0)+'_t'+str(ii+1)+'.png')
            plt.close(fig)
    #save results in one file per frequency as list (each row is one time window or one maxima)
    np.savetxt(params.outdir+'results_f'+str(f0)+'.txt',results,header='wave type, freq (Hz), azimuth (deg), dip (deg), elliptcity, rotation x (deg), wavenumber (/m), velocity (m/s), max norm beam power, SNR',fmt = '%d, %.4f, %.2f, %.2f, %.3f, %.2f, %.4e, %.2f, %.4e, %.2f')

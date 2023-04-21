'''
Convert list of GPS coordinates to xy coordinates (in m UTM) and save in correct file format for B3Ampy
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-21
'''

from datetime import datetime
import time
start_time = time.time()
print("--- starting now! It is %s ---" % (datetime.fromtimestamp(start_time)))
import numpy as np
import utm

###############################################
### read station list in latitude,longitude ###
###############################################
path = './test_data/'
filename = path+'station_list_gps.txt'
stations_gps = np.loadtxt(filename,skiprows=1,dtype=str)

station_names = np.empty((stations_gps.shape[0]),dtype='U4')
station_lat = np.empty((stations_gps.shape[0]))
station_long = np.empty((stations_gps.shape[0]))
for ii in range(stations_gps.shape[0]):
    station_names[ii] = stations_gps[ii][0]
    station_lat[ii] = float(stations_gps[ii][1])
    station_long[ii] = float(stations_gps[ii][2])

##################################
### convert to UTM coordinates ###
##################################
utm_east = np.empty((stations_gps.shape[0]))
utm_north = np.empty((stations_gps.shape[0]))
utm_zone = np.empty((stations_gps.shape[0]))
utm_letter = np.empty((stations_gps.shape[0]),dtype='U1')
for ii in range(stations_gps.shape[0]):
    [utm_east[ii], utm_north[ii], utm_zone[ii], utm_letter[ii]] = utm.from_latlon(station_lat[ii], station_long[ii])

########################################################################
### convert to relative to reference station (use mean as reference) ###
########################################################################
station_x = np.empty((stations_gps.shape[0]))
station_y = np.empty((stations_gps.shape[0]))
for ii in range(stations_gps.shape[0]):
    station_x[ii] = utm_east[ii]-np.mean(utm_east)
    station_y[ii] = utm_north[ii]-np.mean(utm_north)

########################################
### save station coordinates in .txt ###
########################################
stations = np.empty((stations_gps.shape[0],2))
stations[:,0] = station_x
stations[:,1] = station_y
np.savetxt(path+'station_list_xy.txt',np.c_[station_names,stations],fmt='%s,%.8s,%.8s',delimiter=',',header='station name, x (m), y (m)')

### endy bits
print("--- Finished after %s seconds ---" % (time.time()-start_time))

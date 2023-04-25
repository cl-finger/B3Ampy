''''
Find local maxima in 2D array
inspired by Matlab function extrema.m:
Carlos Adrian Vargas Aguilera (2023). extrema.m, extrema2.m (https://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m-extrema2-m), MATLAB Central File Exchange. Retrieved April 21, 2023. 
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-21
'''

import numpy as np

def extrema(P_norm,threshold,dist_threshold,no_maxima):
	ind_y_max = np.zeros(P_norm.shape[0],dtype=int)
	for ii in range(P_norm.shape[0]):
		ind_y_max[ii] = int(np.argmax(P_norm[ii,:]))
	if P_norm[ii,int(ind_y_max[ii])]<threshold:
		ind_y_max[ii] = 0
	ind_x_max = np.zeros(P_norm.shape[1],dtype=int)
	for ii in range(P_norm.shape[1]):
		ind_x_max[ii] = int(np.argmax(P_norm[:,ii]))
	if P_norm[int(ind_x_max[ii]),ii]<threshold:
		ind_x_max[ii] = 0

	#find intersections
	xmax=[]
	ymax=[]
	for ii in range(P_norm.shape[0]):
		if ind_x_max[ind_y_max[ii]] == ii:
			xmax = np.append(xmax,ind_x_max[ind_y_max[ii]]).astype(int)
			ymax = np.append(ymax,ind_y_max[ii]).astype(int)
    
	#remove points smaller than threshold
	for ii in range(xmax.size):
		if P_norm[xmax[ii],ymax[ii]]<threshold:
			xmax[ii] = 0
			ymax[ii] = 0

	#remove points too close to each other
	for ii in range(xmax.size):
		for jj in range(ymax.size):
			if jj>ii: #upper triangle
				dist_max = np.sqrt((xmax[ii]-xmax[jj])**2+(ymax[ii]-ymax[jj])**2)
				if dist_max<dist_threshold:
					if P_norm[xmax[ii],ymax[ii]]>P_norm[xmax[jj],ymax[jj]]:
						xmax[jj] = 0
						ymax[jj] = 0
					elif P_norm[xmax[jj],ymax[jj]]>P_norm[xmax[ii],ymax[ii]]:
						xmax[ii] = 0
						ymax[ii] = 0
	mask1 = xmax == 0
	mask2 = ymax == 0
	mask = np.logical_not(mask1*mask2)
	xmax = xmax[mask]
	ymax = ymax[mask]

	#sort maxima by absolute value
	max_values = P_norm[xmax,ymax]
	idx_sort = np.argsort(max_values)[::-1] #returns indices of how to sort from highest to lowest
	xmax = xmax[idx_sort]	
	ymax = ymax[idx_sort]	

	#remove points that are over limit for max number of maxima
	if xmax.size>no_maxima:	
		xmax = xmax[:no_maxima]
		ymax = ymax[:no_maxima]
	return xmax,ymax

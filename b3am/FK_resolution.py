'''
Create arrays for parameters based on resolution preferences or default
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
last modified: 2023-04-21
'''
import numpy as np
import math
def FK_resolution(params,coords):
	# calculate station network parameters
	dist = np.zeros((params.nstations,params.nstations))
	for ii in range(params.nstations):
		for jj in range(params.nstations):
			dist[ii,jj] = np.sqrt((coords[ii,0]-coords[jj,0])**2+(coords[ii,1]-coords[jj,1])**2)
	dist[dist==0] = np.nan
	dmax = np.nanmax(dist)
	dmin = np.nanmin(dist)

	# Compute wavenumber range
	if params.custom_klimits:
		params.kmin = params.kmin
		params.kmax = params.kmax
	else:
		params.kmin = 1/(3*dmax)
		params.kmax = 1/(2*dmin)
	params.kgrid = np.linspace(params.kmin,params.kmax,params.k_res)
	np.save(params.outdir+'kgrid.npy',params.kgrid)
 
	# Compute backazimuth range
	if params.custom_azi_res:
		params.kth = np.arange(-180,180,params.azi_res)*math.pi/180
	else:
		params.kth = np.arange(-180,180,5)*math.pi/180 # Default angle resolution is 5 deg
	np.save(params.outdir+'kth.npy',params.kth)
	
	# Compute dip range
	if params.custom_dip_res:
		params.dip = np.arange(0,90+params.dip_res,params.dip_res)*math.pi/180
	else:
		params.dip = np.arange(0,90+5,5)*math.pi/180 #default resolution of 5 deg
	np.save(params.outdir+'dip.npy',params.dip)

	# Compute ellipticity range
	if params.custom_ell_res:
		params.ell = np.arange(params.ell_res,2,params.ell_res)
	else:
		params.ell = np.arange(0.05,2,0.05)
	np.save(params.outdir+'ell.npy',params.ell)

	# frequency range for FT
	if params.want_custom_frange:
		params.fmin = params.FT_fmin
		params.fmax = params.FT_fmax
		params.fstep = params.FT_fstep
	else:
		if 1/params.kmax < 500:
			a = 401
			b = 3.39
		elif 1/params.kmax > 500:
			a = 4058
			b = 3818
		params.fmax = np.round(a/(1/params.kmax + b),1)

		if 1/params.kmin < 500:
			a = 401
			b = 3.39
		elif 1/params.kmin > 500:
			a = 4058
			b = 3818
		params.fmin = np.round(a/(1/params.kmin + b),1)
		params.fstep = 0.1
	return params


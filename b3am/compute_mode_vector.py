'''
Calculate mode vectors from given resolution and station locations
------------------------------
Claudia Finger 
claudia.finger@ieg.fraunhofer.de
June 2021
''' 
import math
import numpy as np
from scipy import sparse
from . import polpar2cmplx as p2c


def compute_mode_vector(params,coords):
	#Compute polar grid from angle and wavenumber grids
	k = (sparse.kron(params.kgrid, np.concatenate((np.cos(params.kth),np.sin(params.kth))).reshape((2,params.kth.size)) ) ).toarray()
	#compute mode vectors for k grid
	km = 1/np.sqrt(params.nstations)*np.exp(1j*2*math.pi*np.dot(coords,k)) # in python j is imaginary number #(100,14472)
	
	#define polarisation states using dip and ell with defined resolutions
	#1: azimuth, 2: theta, 3: ellipticity, 4: xi (rotation about x-axis for azi pointing along x-axis)
	params.polstates = np.zeros((params.dip.size*2+1+params.ell.size*2,4))
	params.polstates[:,:] = np.squeeze(np.array([
 		[np.concatenate((np.zeros(params.dip.size*2),[0],np.zeros(params.ell.size*2)))],
 		[np.concatenate((params.dip,params.dip,[math.pi/2],np.ones(params.ell.size*2)*math.pi/2))],
		[np.concatenate((np.zeros(params.dip.size),2*np.ones(params.dip.size),[2],np.tile(params.ell,2)))],
		[np.concatenate((np.ones(params.dip.size*2)*math.pi,[math.pi/2],np.zeros(params.ell.size),np.ones(params.ell.size)*math.pi))]
		]).transpose())
	params.npolstates = params.polstates.shape[0]
	#compute complex polarity ellipse
	Z = p2c.polpar2cmplx(params.polstates)
         
	# Initialize 3C mode vector test matrix
	ka = np.zeros((3*km.shape[0],params.npolstates*km.shape[1]),dtype=complex)
	for ik in range(km.shape[1]): # for each wavenumber grid point
		#wave vector azimuth
		phi0=np.arctan2(k[1,ik],k[0,ik])
		#rotation matrix for this azimuth (rotation around z axis, counter clockwise when looking from above)
		Rz = np.array([[np.cos(phi0),-np.sin(phi0),0],[np.sin(phi0),np.cos(phi0),0],[0,0,1]])
		ka[:,ik*params.npolstates:(ik+1)*params.npolstates] = sparse.kron(np.dot(Rz,Z).T,km[:,ik]).toarray().T 
	return ka,params

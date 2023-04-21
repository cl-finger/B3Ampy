'''
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-02
'''

import math
import numpy as np
from scipy import sparse
 
def polpar2cmplx(ppar):
	# Pre-allocate mode vectors for polarization
	Z = np.zeros((3,ppar.shape[0]),dtype=complex);
	for ii in range(ppar.shape[0]):
		# Extract polarization parameteres and define rotation matrices
		re = ppar[ii,2]
		xi = ppar[ii,3]
		theta = ppar[ii,1]
		phi = ppar[ii,0]
		#phi is always 0 here, so we simplify the rotation as R = RzRyRx and put R directly in
		R = np.array([[-np.sin(theta),np.cos(theta)*np.sin(xi),np.cos(theta)*np.cos(xi)],[0,np.cos(xi),-np.sin(xi)],[np.cos(theta),np.sin(theta)*np.sin(xi),np.sin(theta)*np.cos(xi)]])
		if re<=1:
			la = 1
			lb = re
		else:
			la = 2-re
			lb = 1
		a = np.array([la,0,0])
		b = np.array([0,0,lb])
		a = np.dot(R,a)
		b = np.dot(R,b)
		Z[:,ii] = a - 1j*b
		Z[:,ii] = Z[:,ii]/np.sqrt(np.vdot(Z[:,ii],Z[:,ii])) #normalise Z
	return Z

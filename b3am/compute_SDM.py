'''
Compute cross-spectral density matrix
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-23
'''
import numpy as np

def compute_SDM(f,data):
	fnew = f
	#get M(number of freqs), N(number of time segs) and K(number of locations)
	N,K = np.shape(data)
	M = 1
	#compute SDM
	S = np.zeros((K,K,fnew.size,N),dtype=complex)
	#compute lower triangle SDM (and put conj(.) to upper simultaneously)
	for jj in range(K):
		for ii in range(jj,K): # only compute half of matrix and fill rest with conjugate at end of function
			#cross-correlation without normalization
			tmp = data[:,jj]*np.conj(data[:,ii])
			#reshape into spectrogram of (ii,jj) pair
			tmp = np.reshape(tmp,(M,N))
			#store that spectrogram and its complex conjugate at (i,j) and (j,i)
			S[ii,jj,0,:] = tmp
			S[jj,ii,0,:] = np.conj(tmp)
	return S

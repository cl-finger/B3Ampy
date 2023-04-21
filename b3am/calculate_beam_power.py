'''
Calculate beam power by comparing theoretical mode vectors to recorded data
---
Claudia Finger
claudia.finger@ieg.fraunhofer.de
2023-04-21
'''
import numpy as np

def get_beampower_fast(data,ka,params):
	### calculate beam power using direct multiplication (fast way)
	params.Nt = data.shape[0] #number of time windows
	K = 3*params.nstations
	P = np.zeros((params.kth.size,params.kgrid.size,params.Nt))
	Q = np.zeros((params.kth.size,params.kgrid.size,params.Nt))
	for ip in range(params.Nt): #loop over all time windows
		print('Calculating beam power for time window '+str(ip+1)+' of  '+str(params.Nt))
		Pmu = np.dot(ka.T,data[ip,:].T)
		Pmu = Pmu * np.conj(Pmu) * 1/K #* means elementwise computation, no matrix multiplication here
		Pmu = np.real(Pmu.reshape((params.npolstates,params.kth.size,params.kgrid.size),order='F'))
		kResp = np.max(Pmu,axis=0) 
		kind = np.argmax(Pmu,axis=0)
		P[:,:,ip] = kResp
		Q[:,:,ip] = kind
	return P,Q,params

def get_beampower_SDM(S,ka,params):
	params.Nt = S.shape[3]
	K = 3*params.nstations
	P = np.zeros((params.kth.size,params.kgrid.size,params.Nt))
	Q = np.zeros((params.kth.size,params.kgrid.size,params.Nt))
	S = S.reshape((K,K,params.Nt)) # only one frequency, this eliminates dimension for frequency
	for ip in range(S.shape[2]): #loop over all time windows
		print('Calculating beam power for time window '+str(ip)+' from '+str(S.shape[2]))
		SDM = np.copy(np.squeeze(S[:,:,ip]))
		ks = np.copy(ka)
		Pmu_temp = np.dot(np.conj(ks.T),SDM)
		Pmu2 =np.transpose(Pmu_temp).flatten()*ks.flatten()
		Pmu3 = (Pmu2.reshape((ks.shape[0],ks.shape[1]))) #(300,38*14472)
		Pmu4 = np.copy((1/(K**2))*np.real(np.sum(Pmu3,axis=0))) #sum over stations

		#store FK spectrum for this time window
		Pmu = Pmu4.reshape((params.npolstates,params.kth.size,params.kgrid.size),order='F') 
		kResp = np.max((Pmu),axis=0) 
		kind = np.argmax((Pmu),axis=0)
		P[:,:,ip] = kResp
		Q[:,:,ip] = kind
	return P,Q,params


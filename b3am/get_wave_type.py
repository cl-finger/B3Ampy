'''
Get wave type from polarisation parameters (ellipticity and rotation around x axis)
------------------
Claudia Finger
claudia.finger@ieg.fraunhofer.de
Mar 2023
'''

import math

def get_wave_type(ellipticity,rotation):
	if ellipticity==0:
		wave_type = 0 #P wave                 
	elif ellipticity==2:
		if rotation == math.pi:
			wave_type = 1 #SV wave
		elif rotation == math.pi/2:
			wave_type = 2 #SH/Love
	else:
		if rotation == 0:
			wave_type = 3 #retrograde Rayleigh wave
		elif rotation == math.pi:
			wave_type = 4 #prograde Rayleigh wave
	return wave_type

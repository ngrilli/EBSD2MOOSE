# Nicolò Grilli
# Università di Bristol
# 2 Settembre 2024

# Phase fields eta_i representing the EBSD data

import numpy as np

class PhaseFields:
	
	def __init__(self,file_name,EBSD_data,longitude_sections=1,latitude_sections=1):
		self.file_name = file_name
		self.EBSD_data = EBSD_data # an EBSD object
		# number of nodes along the x and y direction in a structured mesh
		self.Nx = max(self.EBSD_data.Nx_odd,self.EBSD_data.Nx_even) + 1
		self.Ny = self.EBSD_data.Ny + 1
		# the number of sections in which the unit sphere is split along the longitude and latitude
		self.longitude_sections = longitude_sections
		self.latitude_sections = latitude_sections
	
	# calculate section of the sphere in which the c crystal axis 
	# is rotated based on Euler angles
	def sphere_section(self,phi1,Phi,phi2):
		R = rotation_matrix(phi1,Phi,phi2)
		z = np.dot(R, np.array([0,0,1]))
		# calculate latitude and longitude
		if (abs(z[0]) > 1.0e-6):
			longitude = np.arctan(z[1]/z[0])
		elif (z[1] > 1.0e-6):
			longitude = np.pi / 2.0
		elif (z[1] < -1.0e-6):
			longitude = 3.0 * np.pi / 2.0
		else:
			longitude = 0.0
		return 1.0
		
	# calculate ZXZ rotation matrix, angles are in radians
	def rotation_matrix(phi1,Phi,phi2):
		R = np.zeros(shape=(3,3))
		R[0,0] = np.cos(phi1) * np.cos(phi2) - np.sin(phi1) * np.sin(phi2) * np.cos(Phi)
		R[0,1] = - np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(Phi)
		R[0,2] = np.sin(phi1) * np.sin(Phi)
		R[1,0] = np.sin(phi1) * np.cos(phi2) + np.cos(phi1) * np.sin(phi2) * np.cos(Phi)
		R[1,1] = - np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(Phi)
		R[1,2] = - np.cos(phi1) * np.sin(Phi)
		R[2,0] = np.sin(phi2) * np.sin(Phi)
		R[2,1] = np.cos(phi2) * np.sin(Phi)
		R[2,2] = np.cos(Phi)
		return R

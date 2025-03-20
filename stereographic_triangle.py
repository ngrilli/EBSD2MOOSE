# Nicolò Grilli
# Università di Bristol
# 20 Marzo 2025

# a stereographic triangle accounting for cubic symmetries

import numpy as np

class StereographicTriangle:
	
	def __init__(self,file_name,EBSD_data)
		self.file_name = file_name # file name to save image
		self.EBSD_data = EBSD_data # an EBSD object

	# generate rotation matrix (active rotation)
	# from Euler angles angles in radians
	def rotation_matrix(self,phi1,Phi,phi2):
		matrix = np.zeros(shape=(3,3))
		matrix[0,0] = np.cos(phi1) * np.cos(phi2) - np.sin(phi1) * np.sin(phi2) * np.cos(Phi)
		matrix[0,1] = - np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(Phi)
		matrix[0,2] = np.sin(phi1) * np.sin(Phi)
		matrix[1,0] = np.sin(phi1) * np.cos(phi2) + np.cos(phi1) * np.sin(phi2) * np.cos(Phi)
		matrix[1,1] = - np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(Phi)
		matrix[1,2] = - np.cos(phi1) * np.sin(Phi)
		matrix[2,0] = np.sin(phi2) * np.sin(Phi)
		matrix[2,1] = np.cos(phi2) * np.sin(Phi)
		matrix[2,2] = np.cos(Phi)
		return matrix

	# apply FCC symmetries to reduce a 3D vector to the standard triangle
	def apply_symmetry(self,v):
		# make all components positive by applying reflections x, y, z
		w = np.zeros(shape=(3))
		for i in range(3):
			if (v[i] < 0):
				w[i] = -v[i]
			else:
				w[i] = v[i]
		# if y component bigger than x component, apply reflection on (-110)
		if (w[1] > w[0]):
			temp = w[0]
			w[0] = w[1]
			w[1] = temp
		# if y component is bigger than z component, apply reflection on (01-1)
		if (w[1] > w[2]):
			temp = w[2]
			w[2] = w[1]
			w[1] = temp
		# if x component is bigger than z component, apply reflection on (10-1)
		if (w[0] > w[2]):
			temp = w[2]
			w[2] = w[0]
			w[0] = temp
		return w
		
	# calculate stereographic projection of a point (x,y,z) on the unit sphere
	def stereographic_projection(self,x,y,z):
		outstereo = np.zeros(shape=(3))
		outstereo[0] = x / (1.0 + z)
		outstereo[1] = y / (1.0 + z) 
		return outstereo

# example data structure to retrieve the Euler angles
# self.EBSD_data.phi1

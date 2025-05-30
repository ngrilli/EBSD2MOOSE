# Nicolò Grilli
# Università di Bristol
# 28 Maggio 2025

# orientation distribution function with random sampling

import numpy as np
import matplotlib.pyplot as plt
import random

class ODF:
	
	def __init__(self,EBSD_data,N_sampling):
		self.EBSD_data = EBSD_data # an EBSD object
		self.N_sampling = N_sampling # number of query points for random sampling
		self.kappa = 10.0 # concentration parameter of the probability function
		self.N_random_points = 10000 # number of random points on the unit sphere to generate the probability function
		# self.EBSD_points = np.zeros(shape=(len(self.EBSD_data.phi1),3)) # EBSD points on unit sphere
		# a data structure with probabilities? 
		# random choice according to probability
		# np.random.choice(5, 3, p=[0.1, 0, 0.3, 0.6, 0])

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
		
	def von_mises_fisher_kernel(self,query_point,sample_point,kappa):
		# query_point = 3D vector where function is evaluated
		# sample_point = 3D vector of the experimental point on the sphere
		# kappa = concentration parameter of the probability function at sample_point 
		dot_product = np.dot(query_point, sample_point)  # cosine of angle between query_point and sample_point
		C = kappa / (4.0 * np.pi * np.sinh(kappa)) # normalization constant
		return C * np.exp(kappa * dot_product)
		
	# calculate probability distribution function at query_point 
	def calculate_probability(self,query_point,kappa):
		N = len(self.EBSD_data.phi1)
		probability = 0.0 # initialize
		for i in range(N):
			phi1 = self.EBSD_data.phi1[i]
			Phi = self.EBSD_data.Phi[i]
			phi2 = self.EBSD_data.phi2[i]
			R = self.rotation_matrix((np.pi/180.0)*phi1,(np.pi/180.0)*Phi,(np.pi/180.0)*phi2)
			v = R.dot(np.array([1.0,0.0,0.0])) # [100] in the sample reference frame
			probability += self.von_mises_fisher_kernel(query_point,v,kappa)
		return (probability / N)
		
	# generate random points on the unit sphere and calculate probability function based on EBSD
	def generate_random_points(self):
		for i in range(self.N_random_points):
			phi1 = 2.0 * np.pi * random.random()
			Phi = np.arccos(2.0 * random.random() - 1.0)
			phi2 = 2.0 * np.pi * random.random()
			R = self.rotation_matrix(phi1,Phi,phi2)
			v = R.dot(np.array([1.0,0.0,0.0])) # [100] in the sample reference frame
			probability = calculate_probability(v,self.kappa)
		return 0.0 # return a vector with all probabilities?

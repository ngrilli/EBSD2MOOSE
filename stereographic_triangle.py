# Nicolò Grilli
# Zhuohao Song
# Università di Bristol
# 20 Marzo 2025

# a stereographic triangle accounting for cubic symmetries

import numpy as np
import matplotlib.pyplot as plt

class StereographicTriangle:
	
	def __init__(self,EBSD_data):
		self.EBSD_data = EBSD_data # an EBSD object
		# data structures to store coordinates of points in the stereographic projection
		self.x_stereographic = np.zeros(shape=(len(self.EBSD_data.phi1)))
		self.y_stereographic = np.zeros(shape=(len(self.EBSD_data.phi1)))
		# hard coded arbitrary orthonormal directions for stereographic projection
		# this will constitute a reference system for the stereographic triangle
		self.out_of_page_dir = np.array([0.0,0.0,1.0]) # crystallographic direction on the left
		self.towards_right_dir = np.array([1.0,0.0,0.0]) # crystallographic direction on the right
		self.towards_up_dir = np.cross(self.out_of_page_dir,self.towards_right_dir)

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
	# in the crystal reference frame
	def apply_symmetry(self,v):
		# make all components positive by applying reflections
		w = np.zeros(shape=(3))
		for i in range(3):
			w[i] = v[i]
		temp_dot = np.dot(w,self.towards_right_dir)
		if (temp_dot < 0):
			w = w - 2.0 * temp_dot * self.towards_right_dir
		temp_dot = np.dot(w,self.towards_up_dir)
		if (temp_dot < 0):
			w = w - 2.0 * temp_dot * self.towards_up_dir
		temp_dot = np.dot(w,self.out_of_page_dir)
		if (temp_dot < 0):
			w = w - 2.0 * temp_dot * self.out_of_page_dir
		# if towards_up component bigger than towards_right component, apply reflection with respect to n = (-110)
		if (np.dot(w,self.towards_up_dir) > np.dot(w,self.towards_right_dir)):
			n = self.towards_up_dir - self.towards_right_dir
			n = n / np.linalg.norm(n)
			w = w - 2.0 * np.dot(w,n) * n
		# if towards_up component bigger than out_of_page component, apply reflection with respect to n = (01-1)
		if (np.dot(w,self.towards_up_dir) > np.dot(w,self.out_of_page_dir)):
			n = self.towards_up_dir - self.out_of_page_dir
			n = n / np.linalg.norm(n)
			w = w - 2.0 * np.dot(w,n) * n
		# if towards_right component bigger than out_of_page component, apply reflection with respect to n = (10-1)
		if (np.dot(w,self.towards_right_dir) > np.dot(w,self.out_of_page_dir)):
			n = self.towards_right_dir - self.out_of_page_dir
			n = n / np.linalg.norm(n)
			w = w - 2.0 * np.dot(w,n) * n
		return w
		
	# calculate stereographic projection of a point on the unit sphere
	def stereographic_projection(self,x,y,z):
		outstereo = np.zeros(shape=(2))
		outstereo[0] = x / (1.0 + z)
		outstereo[1] = y / (1.0 + z) 
		return outstereo
		
	# fill data structure with stereographic projection points
	def generate_stereographic_projection(self):
		N = len(self.EBSD_data.phi1)
		for i in range(N):
			phi1 = self.EBSD_data.phi1[i]
			Phi = self.EBSD_data.Phi[i]
			phi2 = self.EBSD_data.phi2[i]
			R = self.rotation_matrix((np.pi/180.0)*phi1,(np.pi/180.0)*Phi,(np.pi/180.0)*phi2)
			transposeR = R.transpose() # transforms coordinates from sample reference frame to crystal reference frame
			v = transposeR.dot(np.array([1.0,0.0,0.0])) # stereographic projection with respect to X axis sample direction
			w = self.apply_symmetry(v)
			w0 = np.dot(w,self.towards_right_dir)
			w1 = np.dot(w,self.towards_up_dir)
			w2 = np.dot(w,self.out_of_page_dir)
			outstereo = self.stereographic_projection(w0,w1,w2)
			self.x_stereographic[i] = outstereo[0]
			self.y_stereographic[i] = outstereo[1]
			
	# Spherical linear interpolation (SLERP) between two unit vectors
	def slerp(self, v0, v1, t):
		dot = np.dot(v0, v1) # Dot product
		theta = np.arccos(np.clip(dot, -1.0, 1.0)) # Angle between vectors
		sin_theta = np.sin(theta)
		return (np.sin((1 - t) * theta) * v0 + np.sin(t * theta) * v1) / sin_theta

	# plot orientation points in the stereographic triangle
	def plot_stereographic_projection(self):
		fig_stereo, ax_stereo = plt.subplots()
		ax_stereo.scatter(self.x_stereographic,self.y_stereographic,s=10)
		ax_stereo.spines['top'].set_visible(False)
		ax_stereo.spines['bottom'].set_visible(False)
		ax_stereo.spines['left'].set_visible(False)
		ax_stereo.spines['right'].set_visible(False)
		PlotRadius = 0.515
		PlotAngle = 45.0
		ax_stereo.plot([0.0,PlotRadius*np.cos(PlotAngle*(np.pi/180.0))],[0.0,PlotRadius*np.sin(PlotAngle*(np.pi/180.0))], linewidth=2, color='black')
		ax_stereo.plot([0.0,0.414], [0.0,0.0], linewidth=2, color='black')
		# the stereographic projection of a circle connecting
		# the vectors [1,0,1] and [1,1,1] on the unit sphere
		v0 = np.array([1, 0, 1]) / np.sqrt(2)
		v1 = np.array([1, 1, 1]) / np.sqrt(3)
		angle = np.arccos(np.dot(v0, v1))
		t_values = np.linspace(0, 1, 100)
		arc_x = np.zeros(shape=(0))
		arc_y = np.zeros(shape=(0))
		for t in t_values:
			arc_points = np.array(self.slerp(v0, v1, t))
			temp_stereo = self.stereographic_projection(arc_points[0],arc_points[1],arc_points[2])
			arc_x = np.append(arc_x,temp_stereo[0])
			arc_y = np.append(arc_y,temp_stereo[1])
		ax_stereo.plot(arc_x,arc_y,linewidth=2,color='black')
		ax_stereo.axes.xaxis.set_ticklabels([])
		ax_stereo.axes.yaxis.set_ticklabels([])
		ax_stereo.set_xticks([])
		ax_stereo.set_yticks([])
		ax_stereo.set_xlim(-0.05,PlotRadius)
		ax_stereo.set_ylim(-0.03,0.4)
		fig_stereo.tight_layout()
		fig_stereo.savefig('stereographic_projection.png',dpi=200)

	# divide the stereographic triangle into different regions and classify directions based on their belonging
	def classify_vector(self):
		Nregions == 3
		v001 = np.array([0.0,0.0,1.0])
		v111 = np.array([1.0,1.0,1.0]) / np.sqrt(3.0)
		v101 = np.array([1.0,0.0,1.0]) / np.sqrt(2.0)
		v_from_001_to_101 = np.zeros(shape=(Nregions,3)) # vectors partitioning the stereographic triangle edges
		v_from_001_to_111 = np.zeros(shape=(Nregions,3))
		v_from_101_to_111 = np.zeros(shape=(Nregions,3))
		for i in range(Nregions):
			v_from_001_to_101[i,:] = self.slerp(v001, v101, (i+1.0)/(Nregions+1.0))
			v_from_001_to_111[i,:] = self.slerp(v001, v111, (i+1.0)/(Nregions+1.0))
			v_from_101_to_111[i,:] = self.slerp(v101, v111, (i+1.0)/(Nregions+1.0))
		# vectors partitioning the stereographic triangle areas
		n_horizontal_1 = np.cross(np.squeeze(v_from_001_to_111[0,:]),np.squeeze(v_from_001_to_101[0,:]))
		n_horizontal_2 = np.cross(np.squeeze(v_from_001_to_111[1,:]),np.squeeze(v_from_001_to_101[1,:]))
		n_horizontal_3 = np.cross(np.squeeze(v_from_001_to_111[2,:]),np.squeeze(v_from_001_to_101[2,:]))
		n_vertical_1 = np.cross(np.squeeze(v_from_001_to_111[0,:]),np.squeeze(v_from_101_to_111[0,:]))
		n_vertical_2 = np.cross(np.squeeze(v_from_001_to_111[1,:]),np.squeeze(v_from_101_to_111[1,:]))
		n_vertical_3 = np.cross(np.squeeze(v_from_001_to_111[2,:]),np.squeeze(v_from_101_to_111[2,:]))
		return 0

	# calculate latitude: the angle with respect to out_of_page_dir 
	def calculate_latitude(self,v):
		return np.arccos(np.dot(v,self.out_of_page_dir))

	# calculate longitude: the angle with respect to towards_right direction
	def calculate_longitude(self,v):
		return np.arccos(np.dot(v,self.towards_right_dir))

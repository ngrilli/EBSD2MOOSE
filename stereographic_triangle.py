# Nicolò Grilli
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
		
	# fill data structure with stereographic projection points
	def generate_stereographic_projection(self):
		N = len(self.EBSD_data.phi1)
		for i in range(N):
			phi1 = self.EBSD_data.phi1[i]
			Phi = self.EBSD_data.Phi[i]
			phi2 = self.EBSD_data.phi2[i]
			R = self.rotation_matrix((np.pi/180.0)*phi1,(np.pi/180.0)*Phi,(np.pi/180.0)*phi2)
			v = R.dot(np.array([0.0,0.0,1.0])) # stereographic projection with respect to [001] axis of cubic crystal
			w = self.apply_symmetry(v)
			outstereo = self.stereographic_projection(w[0],w[1],w[2])
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

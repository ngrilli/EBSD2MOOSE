# Nicolò Grilli
# Università di Bristol
# 28 Maggio 2025

# orientation distribution function with random sampling

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ive  # Modified Bessel function

class ODF:
	
	def __init__(self,EBSD_data):
		self.EBSD_data = EBSD_data # an EBSD object
		# data structures to store coordinates of points in the stereographic projection
		#self.x_stereographic = np.zeros(shape=(len(self.EBSD_data.phi1)))
		#self.y_stereographic = np.zeros(shape=(len(self.EBSD_data.phi1)))
		# hard coded arbitrary orthonormal directions for stereographic projection
		# this will constitute a reference system for the stereographic triangle
		#self.out_of_page_dir = np.array([0.0,0.0,1.0]) # crystallographic direction on the left
		#self.towards_right_dir = np.array([1.0,0.0,0.0]) # crystallographic direction on the right
		#self.towards_up_dir = np.cross(self.out_of_page_dir,self.towards_right_dir)
		
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
		
	def von_mises_fisher_kernel(query_point, sample_point, kappa):
		# query_point = 3D vector where function is evaluated
		# sample_point = 3D vector of the experimental point on the sphere
		# kappa = concetration parameter of the probability function at sample_point 
		dot_product = np.dot(query_point, sample_point)  # cosine of angle between query_point and sample_point
		C = kappa / (4.0 * np.pi * np.sinh(kappa)) # normalization constant
		return C * np.exp(kappa * dot_product)
		
	def spherical_odf():


	# ~ # plot orientation points in the stereographic triangle
	# ~ def plot_stereographic_projection(self):
		# ~ fig_stereo, ax_stereo = plt.subplots()
		# ~ ax_stereo.scatter(self.x_stereographic,self.y_stereographic,s=10)
		# ~ ax_stereo.spines['top'].set_visible(False)
		# ~ ax_stereo.spines['bottom'].set_visible(False)
		# ~ ax_stereo.spines['left'].set_visible(False)
		# ~ ax_stereo.spines['right'].set_visible(False)
		# ~ PlotRadius = 0.515
		# ~ PlotAngle = 45.0
		# ~ ax_stereo.plot([0.0,PlotRadius*np.cos(PlotAngle*(np.pi/180.0))],[0.0,PlotRadius*np.sin(PlotAngle*(np.pi/180.0))], linewidth=2, color='black')
		# ~ ax_stereo.plot([0.0,0.414], [0.0,0.0], linewidth=2, color='black')
		# ~ # the stereographic projection of a circle connecting
		# ~ # the vectors [1,0,1] and [1,1,1] on the unit sphere
		# ~ v0 = np.array([1, 0, 1]) / np.sqrt(2)
		# ~ v1 = np.array([1, 1, 1]) / np.sqrt(3)
		# ~ angle = np.arccos(np.dot(v0, v1))
		# ~ t_values = np.linspace(0, 1, 100)
		# ~ arc_x = np.zeros(shape=(0))
		# ~ arc_y = np.zeros(shape=(0))
		# ~ for t in t_values:
			# ~ arc_points = np.array(self.slerp(v0, v1, t))
			# ~ temp_stereo = self.stereographic_projection(arc_points[0],arc_points[1],arc_points[2])
			# ~ arc_x = np.append(arc_x,temp_stereo[0])
			# ~ arc_y = np.append(arc_y,temp_stereo[1])
		# ~ ax_stereo.plot(arc_x,arc_y,linewidth=2,color='black')
		# ~ ax_stereo.axes.xaxis.set_ticklabels([])
		# ~ ax_stereo.axes.yaxis.set_ticklabels([])
		# ~ ax_stereo.set_xticks([])
		# ~ ax_stereo.set_yticks([])
		# ~ ax_stereo.set_xlim(-0.05,PlotRadius)
		# ~ ax_stereo.set_ylim(-0.03,0.4)
		# ~ fig_stereo.tight_layout()
		# ~ fig_stereo.savefig('stereographic_projection.png',dpi=200)

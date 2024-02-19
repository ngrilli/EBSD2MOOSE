# Nicolò Grilli
# Università di Bristol
# 19 Febbraio 2024

# An Abaqus input file with structured mesh

import numpy as np

class AbaqusInputFile:
	
	def __init__(self, file_name, EBSD_data, thickness):
		self.file_name = file_name
		self.EBSD_data = EBSD_data # an EBSD object
		self.thickness = thickness
		self.nodes = np.zeros(shape=(0))
		self.elements = np.zeros(shape=(0))
		self.Nx = max(EBSD_data.Nx_odd,EBSD_data.Nx_even)
		self.Ny = EBSD_data.Ny
		self.Nz = self.thickness
		# list of nodes for the boundaries
		self.left = np.zeros(shape=(0))
		self.right = np.zeros(shape=(0))
		self.bottom = np.zeros(shape=(0))
		self.top = np.zeros(shape=(0))
		self.back = np.zeros(shape=(0))
		self.front = np.zeros(shape=(0))
		
	# create data structure with node coordinates
	def generate_nodes(self):
		# use unitary element size for the moment
		# frequency option missing
		self.nodes = np.zeros(shape=((self.Nx+1)*(self.Ny+1)*(self.Nz+1),3))
		for x in range(self.Nx+1):
			for y in range(self.Ny+1):
				for z in range(self.Nz+1):
					node_index = x + y*(self.Nx+1) + z*(self.Nx+1)*(self.Ny+1)
					self.nodes[node_index,0] = x # * element size
					self.nodes[node_index,1] = y
					self.nodes[node_index,2] = z
					if (x == 0):
						self.left = np.append(self.left, node_index+1)
					elif (x == self.Nx):
						self.right = np.append(self.right, node_index+1)
		
	# create data structure with node indices for each element
	def generate_elements(self):
		self.elements = np.zeros(shape=(self.Nx*self.Ny*self.Nz,8))
		for x in range(self.Nx):
			for y in range(self.Ny):
				for z in range(self.Nz):
					element_index = x + y*self.Nx + z*self.Nx*self.Ny
					node_indices = np.zeros(shape=(8))
					node_indices[0] = x + (self.Nx+1)*y + (self.Nx+1)*(self.Ny+1)*z + 1
					node_indices[1] = x + (self.Nx+1)*y + (self.Nx+1)*(self.Ny+1)*z + 2
					node_indices[2] = x + (self.Nx+1)*(y+1) + (self.Nx+1)*(self.Ny+1)*z + 2
					node_indices[3] = x + (self.Nx+1)*(y+1) + (self.Nx+1)*(self.Ny+1)*z + 1
					for lower_node in range(4):
						node_indices[4+lower_node] = node_indices[lower_node] + (self.Nx+1)*(self.Ny+1)
					for node in range(8):
						self.elements[element_index,node] = node_indices[node]
		
	def write_input_file(self):
		self.generate_nodes()
		self.generate_elements()
		fid = open(self.file_name,'w')
		fid.write("*Part, NAME=TEST" + "\n")
		fid.write("*NODE" + "\n")
		for node in range((self.Nx+1)*(self.Ny+1)*(self.Nz+1)):
			fid.write("{:d}".format(int(node+1)))
			fid.write(", ")
			fid.write("{:.1f}".format(self.nodes[node,0]))
			fid.write(", ")
			fid.write("{:.1f}".format(self.nodes[node,1]))
			fid.write(", ")
			fid.write("{:.1f}".format(self.nodes[node,2]))
			fid.write("\n")
		fid.write("*ELEMENT, TYPE=C3D8" + "\n")
		for elem in range(self.Nx*self.Ny*self.Nz):
			fid.write("{:d}".format(int(elem+1)))
			for node in range(8):
				fid.write(", ")
				fid.write("{:d}".format(int(self.elements[elem,node])))
			fid.write("\n")
		fid.write("*NSET, NSET=LEFT" + "\n")
		fid.write("*NSET, NSET=RIGHT")
		fid.write("*NSET, NSET=BOTTOM")
		fid.write("*NSET, NSET=TOP")
		fid.write("*NSET, NSET=BACK")
		fid.write("*NSET, NSET=FRONT")
		fid.close()

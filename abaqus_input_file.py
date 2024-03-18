# Nicolò Grilli
# Università di Bristol
# 19 Febbraio 2024

# An Abaqus input file with structured mesh

import numpy as np

class AbaqusInputFile:
	
	def __init__(self, file_name,EBSD_data,thickness=1):
		self.file_name = file_name
		self.EBSD_data = EBSD_data # an EBSD object
		self.thickness = thickness
		self.nodes = np.zeros(shape=(0))
		self.elements = np.zeros(shape=(0))
		self.Nelements = 0
		self.Nx = 0 #max(EBSD_data.Nx_odd,EBSD_data.Nx_even)
		self.Ny = 0 #EBSD_data.Ny
		self.Nz = self.thickness
		# list of nodes for the boundaries
		self.left = np.zeros(shape=(0))
		self.right = np.zeros(shape=(0))
		self.bottom = np.zeros(shape=(0))
		self.top = np.zeros(shape=(0))
		self.back = np.zeros(shape=(0))
		self.front = np.zeros(shape=(0))
		
	# calculate dimensions of this mesh
	def calculate_dimensions(self,frequency=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		if (nx_max < 0): # negative means no upper limit
			nx_max = max(self.EBSD_data.Nx_odd,self.EBSD_data.Nx_even)
		if (ny_max < 0):
			ny_max = self.EBSD_data.Ny
		self.Nx = int((nx_max-nx_min-2)/frequency)+1
		self.Ny = int((ny_max-ny_min-2)/frequency)+1
		
	# create data structure with node coordinates
	def generate_nodes(self,frequency=1):
		# use unitary pixel size for the moment
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
					if (y == 0):
						self.bottom = np.append(self.bottom, node_index+1)
					elif (y == self.Ny):
						self.top = np.append(self.top, node_index+1)
					if (z == 0):
						self.back = np.append(self.back, node_index+1)
					elif (z == self.Nz):
						self.front = np.append(self.front, node_index+1)
		
	# create data structure with node indices for each element
	def generate_elements(self,frequency=1):
		self.Nelements = self.Nx*self.Ny*self.Nz
		self.elements = np.zeros(shape=(self.Nelements,8))
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
		
	def write_input_file(self,frequency=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		self.calculate_dimensions(frequency,nx_min,nx_max,ny_min,ny_max)
		self.generate_nodes(frequency)
		self.generate_elements(frequency)
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
		for node in range(len(self.left)):
			fid.write("{:d}".format(int(self.left[node])))
			if (node == len(self.left)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")
		fid.write("*NSET, NSET=RIGHT" + "\n")
		for node in range(len(self.right)):
			fid.write("{:d}".format(int(self.right[node])))
			if (node == len(self.right)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")		
		fid.write("*NSET, NSET=BOTTOM" + "\n")
		for node in range(len(self.bottom)):
			fid.write("{:d}".format(int(self.bottom[node])))
			if (node == len(self.bottom)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")			
		fid.write("*NSET, NSET=TOP" + "\n")
		for node in range(len(self.top)):
			fid.write("{:d}".format(int(self.top[node])))
			if (node == len(self.top)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")	
		fid.write("*NSET, NSET=BACK" + "\n")
		for node in range(len(self.back)):
			fid.write("{:d}".format(int(self.back[node])))
			if (node == len(self.back)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")	
		fid.write("*NSET, NSET=FRONT" + "\n")
		for node in range(len(self.front)):
			fid.write("{:d}".format(int(self.front[node])))
			if (node == len(self.front)-1 or ((node+1)%6 == 0 and node > 0)):
				fid.write("\n")
			else:
				fid.write(", ")	
		# write an element set for each element
		for elem in range(self.Nelements):
			fid.write("*ELSET, ELSET=GRAIN" + str(elem+1) + "\n")
			fid.write("{:d}".format(int(elem+1)))
			fid.write("\n")
		fid.close()
		
	# convert inp file to med file
	def inp2med(self):
		import meshio
		inpmesh = meshio.read(self.file_name)
		medmesh = meshio.read(self.file_name) # initialize medmesh with fields from inpmesh
		medmesh.cell_tags = {} # add specific missing med tags dictionaries
		medmesh.point_tags = {}
		medmesh.cell_sets.clear() # clean unecessary info in medmesh
		medmesh.cell_sets_dict.clear()
		medmesh.point_sets.clear()
		del medmesh.cells[1:len(medmesh.cells)]
		list_cell_groups = [] # define cell groups
		for cellid in range(len(inpmesh.cells[0])):
			temp_cell_group = []
			for i, key in enumerate(inpmesh.cell_sets):
				if cellid in inpmesh.cell_sets[key][0]:
					temp_cell_group.append(key)       
			list_cell_groups.append(temp_cell_group)
		unique_cell_groups = [list(x) for x in set(tuple(x) for x in list_cell_groups) if x]
		for i, group in enumerate(unique_cell_groups):
			medmesh.cell_tags[-6-i] = group
		ctags = [0]*len(inpmesh.cells[0]) # define cell tags
		for i, item in enumerate(list_cell_groups):
			for j in medmesh.cell_tags:
				if item == medmesh.cell_tags[j]:
					ctags[i] = j
		medmesh.cell_data['cell_tags'] = [np.array(ctags)]
		medmesh.cell_data_dict['cell_tags'] = {list(inpmesh.cells_dict.keys())[0]:np.array(ctags)}
		list_point_groups = [] # define node groups
		for nodeid in range(len(inpmesh.points)):
			temp_point_group = []
			for i, key in enumerate(inpmesh.point_sets):
				if nodeid in inpmesh.point_sets[key]:
					temp_point_group.append(key)       
			list_point_groups.append(temp_point_group)
		unique_point_groups = [list(x) for x in set(tuple(x) for x in list_point_groups) if x]
		for i, group in enumerate(unique_point_groups): # define nodeset tags
			medmesh.point_tags[2+i] = group
		ptags = [0]*len(inpmesh.points)
		for i, item in enumerate(list_point_groups):
			for j in medmesh.point_tags:
				if item == medmesh.point_tags[j]:
					ptags[i] = j
		medmesh.point_data['point_tags'] = np.array(ptags)
		medmesh.write(self.file_name.split('.')[0] + '.med') # write mesh file

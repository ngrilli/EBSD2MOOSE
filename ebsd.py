# Nicolò Grilli
# Michael Salvini
# University of Bristol
# 7 Maggio 2023

# An EBSD measurement

import numpy as np
import matplotlib.pyplot as plt

class EBSD:
	
	def __init__(self, file_name):
		self.file_name = file_name
		self.file_type = self.file_name[-3:]
		self.Nx_odd = 0 # .ang file can have odd columns with a different size compared with even columns
		self.Nx_even = 0
		self.Ny = 0
		self.euler_start_line = 0
		self.phi1 = np.zeros(shape=(0)) # Euler angles in degrees
		self.Phi = np.zeros(shape=(0))
		self.phi2 = np.zeros(shape=(0))
		self.crystal_structure = np.zeros(shape=(0)) # phase
		self.phi1_map = np.zeros(shape=(0,0)) # 2D Euler angles maps
		self.Phi_map = np.zeros(shape=(0,0))
		self.phi2_map = np.zeros(shape=(0,0))
		self.crystal_structure_map = np.zeros(shape=(0,0)) # 2D phase map

	# Parse all information about the EBSD map
	def parse_ebsd_file(self):
		fid = open(self.file_name,'r')
		# ang file type
		if self.file_type == 'ang':
			for i, line in enumerate(fid):
				if line[:13] == '# NCOLS_ODD: ':
					self.Nx_odd = int(line.split('# NCOLS_ODD: ',1)[1])
				if line[:14] == '# NCOLS_EVEN: ':
					self.Nx_even = int(line.split('# NCOLS_EVEN: ',1)[1])
				if line[:9] == '# NROWS: ':
					self.Ny = int(line.split('# NROWS: ',1)[1])
				# start reading Euler angles
				if self.euler_start_line > 0:
					self.phi1[i-self.euler_start_line] = (180/np.pi) * float(line.split()[0])
					self.Phi[i-self.euler_start_line] = (180/np.pi) * float(line.split()[1])
					self.phi2[i-self.euler_start_line] = (180/np.pi) * float(line.split()[2])
					self.crystal_structure[i-self.euler_start_line] = int(line.split()[7])
				elif line[:2] == '  ': # start of Euler angles
					self.euler_start_line = i
					# size the data structures for storing Euler angles
					if (self.Ny % 2 == 0):
						size_of_EBSD_map = int((self.Nx_odd + self.Nx_even) * self.Ny / 2)
					else:
						size_of_EBSD_map = int((self.Nx_odd + self.Nx_even) * (self.Ny - 1) / 2 + self.Nx_odd)
					self.phi1 = np.zeros(shape=(size_of_EBSD_map))
					self.Phi = np.zeros(shape=(size_of_EBSD_map))
					self.phi2 = np.zeros(shape=(size_of_EBSD_map))
					self.crystal_structure = np.zeros(shape=(size_of_EBSD_map))
					# assign data in current line
					self.phi1[i-self.euler_start_line] = (180/np.pi) * float(line.split()[0])
					self.Phi[i-self.euler_start_line] = (180/np.pi) * float(line.split()[1])
					self.phi2[i-self.euler_start_line] = (180/np.pi) * float(line.split()[2])
					self.crystal_structure[i-self.euler_start_line] = int(line.split()[7])					
		# ctf file type: assuming odd columns size = even columns size
		if self.file_type == 'ctf':
			for i, line in enumerate(fid):
				if line[:7] == 'XCells	':
					self.Nx_odd = int(line.split('XCells	',1)[1])
					self.Nx_even = int(line.split('XCells	',1)[1])
				if line[:7] == 'YCells	':
					self.Ny = int(line.split('YCells	',1)[1])
				# start reading Euler angles
				if self.euler_start_line > 0:
					self.phi1[i-self.euler_start_line] = float(line.split()[5])
					self.Phi[i-self.euler_start_line] = float(line.split()[6])
					self.phi2[i-self.euler_start_line] = float(line.split()[7])
					self.crystal_structure[i-self.euler_start_line] = int(line.split()[0])
				elif line[:6] == 'Phase	':
					self.euler_start_line = i + 1
					# size the data structures for storing Euler angles
					self.phi1 = np.zeros(shape=(self.Nx_odd*self.Ny))
					self.Phi = np.zeros(shape=(self.Nx_odd*self.Ny))
					self.phi2 = np.zeros(shape=(self.Nx_odd*self.Ny))
					self.crystal_structure = np.zeros(shape=(self.Nx_odd*self.Ny))
		fid.close()
		
	# Generate an Euler angles file in MOOSE format
	# EBSD data are parsed with a specific frequency
	# and element number thickness can be selected
	def generate_MOOSE_Euler_angles_file(self,filename='EulerAngles.txt',frequency=1,thickness=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		euler_angles_file = open(filename,"w")
		max_Nx = max(self.Nx_odd,self.Nx_even)
		if (nx_max < 0): # negative means no upper limit
			nx_max = max_Nx
		if (ny_max < 0):
			ny_max = self.Ny
		for z in range(int(thickness)):
			for y in range(self.Ny):
				for x in range(max_Nx):
					# check frequency
					if (x%frequency == 0 and y%frequency == 0):
						# check boundaries
						if (x > nx_min and x < nx_max and y > ny_min and y < ny_max):
							euler_angles_file.write('{:0.2f}'.format(self.phi1_map[x,y]))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.2f}'.format(self.Phi_map[x,y]))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.2f}'.format(self.phi2_map[x,y]))
							euler_angles_file.write('\n')
		euler_angles_file.close()
		# print information about the size of the map
		print('Size of the Euler angles file')
		print('nx = ' + str(int((nx_max-nx_min-2)/frequency)+1))
		print('ny = ' + str(int((ny_max-ny_min-2)/frequency)+1))
		
	# Generate an Euler angles file in UMAT format
	def generate_UMAT_Euler_angles_file(self,filename='materials.dat',frequency=1,thickness=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		euler_angles_file = open(filename,"w")
		max_Nx = max(self.Nx_odd,self.Nx_even)
		if (nx_max < 0): # negative means no upper limit
			nx_max = max_Nx
		if (ny_max < 0):
			ny_max = self.Ny
		elem = 0 # incremental element counter
		for z in range(int(thickness)):
			for y in range(self.Ny):
				for x in range(max_Nx):
					# check frequency
					if (x%frequency == 0 and y%frequency == 0):
						# check boundaries
						if (x > nx_min and x < nx_max and y > ny_min and y < ny_max):
							elem += 1
							euler_angles_file.write('{:0.0f}'.format(elem))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.2f}'.format(self.phi1_map[x,y]))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.2f}'.format(self.Phi_map[x,y]))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.2f}'.format(self.phi2_map[x,y]))
							euler_angles_file.write(' ')
							euler_angles_file.write('{:0.0f}'.format(1))
							euler_angles_file.write(' ')
							# TO DO: multiphase in UMAT Euler angles file
							euler_angles_file.write('{:0.0f}'.format(1))
							euler_angles_file.write('\n')
		euler_angles_file.close()
		# print information about the size of the map
		print('Size of the Euler angles file')
		print('nx = ' + str(int((nx_max-nx_min-2)/frequency)+1))
		print('ny = ' + str(int((ny_max-ny_min-2)/frequency)+1))

	# plot the EBSD map
	def plot_EBSD_map(self,frequency=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		max_Nx = max(self.Nx_odd,self.Nx_even)
		if (nx_max < 0): # negative means no upper limit
			nx_max = max_Nx
		if (ny_max < 0):
			ny_max = self.Ny
		if (nx_min < 0):
			nx_min = 0
		if (ny_min < 0):
			ny_min = 0
		fig_phi1, ax_phi1 = plt.subplots()
		ax_phi1.contourf(np.transpose(self.phi1_map[nx_min:nx_max:frequency,ny_min:ny_max:frequency]),cmap='RdYlBu_r')
		ax_phi1.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=False,labelbottom=False,labelleft=False)
		fig_phi1.savefig('phi1.png',dpi=200)
		fig_Phi, ax_Phi = plt.subplots()
		ax_Phi.contourf(np.transpose(self.Phi_map[nx_min:nx_max:frequency,ny_min:ny_max:frequency]),cmap='RdYlBu_r')
		ax_Phi.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=False,labelbottom=False,labelleft=False)
		fig_Phi.savefig('Phi.png',dpi=200)
		fig_phi2, ax_phi2 = plt.subplots()
		ax_phi2.contourf(np.transpose(self.phi2_map[nx_min:nx_max:frequency,ny_min:ny_max:frequency]),cmap='RdYlBu_r')
		ax_phi2.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=False,labelbottom=False,labelleft=False)
		fig_phi2.savefig('phi2.png',dpi=200)
		
	# convert cell number to coordinates
	def cell2coords(self,elem):
		# check if elem belongs to an odd or even column
		if (elem % (self.Nx_odd + self.Nx_even) < self.Nx_odd): # an element in an odd column
			x = elem % (self.Nx_odd + self.Nx_even)
			y = int(elem / (self.Nx_odd + self.Nx_even)) * 2
		else: # an element in an even column
			x = elem % (self.Nx_odd + self.Nx_even) - self.Nx_odd
			y = int(elem / (self.Nx_odd + self.Nx_even)) * 2 + 1
		return [x,y]
		
	# convert coordinates to cell number
	def coords2cell(self,x,y):
		if (y % 2 == 0): # an element in an odd column
			elem = int(y / 2) * (self.Nx_odd + self.Nx_even) + x
		else: # an element in an even column
			elem = int(y / 2) * (self.Nx_odd + self.Nx_even) + self.Nx_odd + x 
		return elem

	# create a 2D data structure for Euler angles
	def generate_2D_Euler_angles_map(self):
		max_Nx = max(self.Nx_odd,self.Nx_even)
		self.phi1_map = np.zeros(shape=(max_Nx,self.Ny))
		self.Phi_map = np.zeros(shape=(max_Nx,self.Ny))
		self.phi2_map = np.zeros(shape=(max_Nx,self.Ny))
		self.crystal_structure_map = np.zeros(shape=(max_Nx,self.Ny))
		for elem in range(len(self.phi1)):
			[x,y] = self.cell2coords(elem)
			self.phi1_map[x,y] = self.phi1[elem]
			self.Phi_map[x,y] = self.Phi[elem]
			self.phi2_map[x,y] = self.phi2[elem]
			self.crystal_structure_map[x,y] = self.crystal_structure[elem]

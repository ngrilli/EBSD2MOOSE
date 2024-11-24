# Nicolò Grilli
# Università di Bristol
# 24 Novembre 2024

# multiphase material

import numpy as np

class Multiphase:
	
	def __init__(self,file_name,EBSD_data):
		self.file_name = file_name
		self.EBSD_data = EBSD_data # an EBSD object
		
	# Generate a phase file in MOOSE format
	def generate_MOOSE_phase_file(self,filename='Phase.txt',frequency=1,thickness=1,nx_min=-1,nx_max=-1,ny_min=-1,ny_max=-1):
		phase_file = open(filename,"w")
		max_Nx = max(self.EBSD_data.Nx_odd,self.EBSD_data.Nx_even)
		if (nx_max < 0): # negative means no upper limit
			nx_max = max_Nx
		if (ny_max < 0):
			ny_max = self.EBSD_data.Ny
		for z in range(int(thickness)):
			for y in range(self.EBSD_data.Ny):
				for x in range(max_Nx):
					# check frequency
					if (x%frequency == 0 and y%frequency == 0):
						# check boundaries
						if (x > nx_min and x < nx_max and y > ny_min and y < ny_max):
							phase_file.write(str(int(self.EBSD_data.crystal_structure_map[x,y])))
							phase_file.write('\n')
		phase_file.close()
		# print information about the size of the phase map
		print('Size of the phase file')
		print('nx = ' + str(int((nx_max-nx_min-2)/frequency)+1))
		print('ny = ' + str(int((ny_max-ny_min-2)/frequency)+1))

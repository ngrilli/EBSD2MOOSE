# Nicolò Grilli
# Arnav Peehal
# Università of Bristol
# 17 Agosto 2026

# A Neper generated polycrystal

import numpy as np
import meshio
from numpy import random as rd

class Neper:

    def __init__(self, file_name):
        self.file_name = file_name # Neper generated mesh file
        self.Nx = 0
        self.Ny = 0
        self.phi1_map = np.zeros(shape=(0,0)) # 2D Euler angles maps
        self.Phi_map = np.zeros(shape=(0,0))
        self.phi2_map = np.zeros(shape=(0,0))

    def parse_mesh_file(self):
         self.mesh = meshio.read(self.file_name)
         self.number_of_grains = len(self.mesh.cell_sets)
         points = self.mesh.points[:,:] # nodal coordinates
         element_type = list(self.mesh.cells_dict.keys())
         cells = self.mesh.cells_dict[element_type[0]] # elements and corresponding node lists
         centres = np.mean(points[cells], axis=1) # element center coordinates
         self.Nx = len(np.unique(centres[:, 0]))
         self.Ny = len(np.unique(centres[:, 1]))
         self.grain = np.zeros(shape=(self.Nx,self.Ny)) # grain index at each coordinate
         for grain_index, crystal in enumerate(self.mesh.cell_sets):
              element_set = self.mesh.cell_sets[crystal][0]
              for elem in element_set:
                [nx,ny] = self.cell2coords(elem)
                self.grain[nx,ny] = grain_index # zero based grain index
         print(self.grain)

	# convert cell number to coordinates
    def cell2coords(self,elem):
        x = elem % self.Nx
        y = int(elem / self.Nx) % self.Ny
        return [x,y]

    def generate_interface_file(self):
        # TO DO: grain boundary detection
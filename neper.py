# Nicolò Grilli
# Arnav Peehal
# Università of Bristol
# 17 Agosto 2026

# A Neper generated polycrystal

import numpy as np
import meshio
from numpy import random as rd
import matplotlib.pyplot as plt

class Neper:

    def __init__(self, file_name):
        self.file_name = file_name # Neper generated mesh file
        # example neper commands
        # neper -T -n 39 -domain 'cube(50.0,50.0,0.5)'
        # neper -M n39-id1.tess -elttype 'hex' -cl 0.5 -format 'inp' -order 1 -dim 3
        self.Nx = 0
        self.Ny = 0
        #self.phi1_map = np.zeros(shape=(0,0)) # 2D Euler angles maps
        #self.Phi_map = np.zeros(shape=(0,0))
        #self.phi2_map = np.zeros(shape=(0,0))

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

	# convert cell number to coordinates
    def cell2coords(self,elem):
        x = elem % self.Nx
        y = int(elem / self.Nx) % self.Ny
        return [x,y]

    def generate_interface_file(self):
        self.interface = np.zeros(shape=(self.Nx,self.Ny))
        delta_nx = 1
        delta_ny = 1
        interface_file = open("interface.txt","w")
        for nx in range(0,self.Nx): # search for interfaces
            for ny in range(0,self.Ny):
                for neigh_x in range(-delta_nx,delta_nx):
                    for neigh_y in range(-delta_ny,delta_ny):
                        if (nx+neigh_x >= 0 and nx+neigh_x < self.Nx):
                            if (ny+neigh_y >= 0 and ny+neigh_y < self.Ny):
                                if (self.grain[nx,ny] != self.grain[nx+neigh_x,ny+neigh_y]):
                                    self.interface[nx,ny] = 1
                interface_file.write('{:0.0f}'.format(self.interface[nx,ny]))
                interface_file.write('\n')
        interface_file.close()
        fig_GB, ax_GB = plt.subplots()
        ax_GB.contourf(np.transpose(np.squeeze(self.interface[:,:])),cmap='RdYlBu_r')
        ax_GB.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=False,labelbottom=False,labelleft=False)
        fig_GB.savefig('interface.png',dpi=200)

    def generate_euler_angles_file(self):
        self.phi1_map = np.zeros(shape=(self.Nx,self.Ny))
        self.Phi_map = np.zeros(shape=(self.Nx,self.Ny))
        self.phi2_map = np.zeros(shape=(self.Nx,self.Ny))
        self.phi1 = np.zeros(shape=(self.number_of_grains))
        self.Phi = np.zeros(shape=(self.number_of_grains))
        self.phi2 = np.zeros(shape=(self.number_of_grains))
        euler_angles_file = open("euler_angles.txt","w")
        # generate random Euler angles for each grain: uniform distribution on a sphere
        for grain_index in range(self.number_of_grains):
            self.phi1[grain_index] = 360.0 * rd.random()
            self.Phi[grain_index] = (180.0 / np.pi) * np.arccos(2.0 * rd.random() - 1.0)
            self.phi2[grain_index] = 360.0 * rd.random()
        for nx in range(0,self.Nx):
            for ny in range(0,self.Ny):
                self.phi1_map[nx,ny] = self.phi1[int(self.grain[nx,ny])]
                self.Phi_map[nx,ny] = self.Phi[int(self.grain[nx,ny])]
                self.phi2_map[nx,ny] = self.phi2[int(self.grain[nx,ny])]
                euler_angles_file.write('{:0.2f}'.format(self.phi1_map[nx,ny]))
                euler_angles_file.write(' ')
                euler_angles_file.write('{:0.2f}'.format(self.Phi_map[nx,ny]))
                euler_angles_file.write(' ')
                euler_angles_file.write('{:0.2f}'.format(self.phi2_map[nx,ny]))
                euler_angles_file.write('\n')
        euler_angles_file.close()
        fig_euler, ax_euler = plt.subplots()
        ax_euler.contourf(np.transpose(self.phi1_map[:,:]),cmap='RdYlBu_r')
        ax_euler.tick_params(axis='both',which='both',bottom=False,top=False,right=False,left=False,labelbottom=False,labelleft=False)
        fig_euler.savefig('phi1.png',dpi=200)
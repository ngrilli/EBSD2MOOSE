# Nicolò Grilli
# University of Bristol
# 1 Giugno 2023

# convert EBSD into Euler angles file
# compatible with MOOSE framework or UMAT

import numpy as np
import sys
import argparse
from ebsd import EBSD
from abaqus_input_file import AbaqusInputFile
from multiphase import Multiphase
from stereographic_triangle import StereographicTriangle

parser = argparse.ArgumentParser(prog='EBSD2MOOSE', \
                                 description='Convert EBSD ctf and ang files to MOOSE or UMAT Euler angles files', \
                                 epilog='Nicolò Grilli, University of Bristol')

parser.add_argument('filename') # positional argument
parser.add_argument('-f','--frequency',type=int,default=1) # option that takes a value
parser.add_argument('-t','--thickness',type=int,default=1)
parser.add_argument('-nx_min','--nx_min',type=int,default=-1)
parser.add_argument('-nx_max','--nx_max',type=int,default=-1)
parser.add_argument('-ny_min','--ny_min',type=int,default=-1)
parser.add_argument('-ny_max','--ny_max',type=int,default=-1)
parser.add_argument('-UMAT','--UMAT',action='store_true')  # on/off flag
parser.add_argument('-aster','--aster',action='store_true')
parser.add_argument('-multiphase','--multiphase',action='store_true')
parser.add_argument('-stereographic_projection','--stereographic_projection',action='store_true')

args = parser.parse_args()

data = EBSD(args.filename)

data.parse_ebsd_file()
data.generate_2D_Euler_angles_map()

data.plot_EBSD_map(args.frequency,args.nx_min,args.nx_max,args.ny_min,args.ny_max)

if (args.UMAT):
	data.generate_UMAT_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0]+'.txt',args.frequency,args.thickness,args.nx_min,args.nx_max,args.ny_min,args.ny_max)
	abaqus_file = AbaqusInputFile(args.filename.rsplit('.', maxsplit=1)[0]+'.inp',data,args.thickness)
	abaqus_file.write_input_file(args.frequency,args.nx_min,args.nx_max,args.ny_min,args.ny_max)
	if (args.aster):
		abaqus_file.inp2med()
else:
	data.generate_MOOSE_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0]+'.txt',args.frequency,args.thickness,args.nx_min,args.nx_max,args.ny_min,args.ny_max)
	if (args.multiphase):
		phase = Multiphase(args.filename,data)
		phase.generate_MOOSE_phase_file(args.filename.rsplit('.', maxsplit=1)[0]+'_phase.txt',args.frequency,args.thickness,args.nx_min,args.nx_max,args.ny_min,args.ny_max)
	if (args.stereographic_projection):
		stereo = StereographicTriangle(data)
		stereo.generate_stereographic_projection()
		stereo.plot_stereographic_projection()

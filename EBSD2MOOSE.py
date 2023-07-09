# Nicolò Grilli
# University of Bristol
# 1 Giugno 2023

# convert EBSD into Euler angles file
# compatible with MOOSE framework or UMAT

import numpy as np
from ebsd import EBSD
import sys
import argparse

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

args = parser.parse_args()

data = EBSD(1,args.filename)

data.parse_ebsd_file()
data.generate_2D_Euler_angles_map()

data.plot_EBSD_map(args.frequency,args.nx_min,args.nx_max,args.ny_min,args.ny_max)

if (args.UMAT):
	data.generate_UMAT_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0]+'.txt',args.frequency,args.thickness)
else:
	data.generate_MOOSE_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0]+'.txt',args.frequency,args.thickness,args.nx_min,args.nx_max,args.ny_min,args.ny_max)

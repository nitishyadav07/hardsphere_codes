#!/usr/bin/env python

from matplotlib.pyplot import imshow
import numpy as np
from scipy import *
from pylab import *
import math as math
import sys # to write data to file
import csv # for handling csv files, data appending
import os
from os import path		# for checking if data file already exists or not
from numpy import array, savetxt # for saving csv file
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d # for 3-d plots of coordinate data
from config_io_module import *

from func_hs import *		# function definitions
from constants import *		# physical constants


# USING ARGUMENT PARSER (argparse package) to use arguments (in case of more than 2 arguments passed from terminal to script):
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument('--n', nargs = 1)
parser.add_argument('--ds', nargs = 1)
parser.add_argument('--iter', nargs = 1)
#parser.add_argument('--box', nargs = 1)
args = parser.parse_args()
#print(args)
#n = int(args.n[0])
ds = float(args.ds[0])			# ds is the step size in real space
N_iter = float(args.iter[0])	# no. of MC iteration steps
#box = float(args.box[0])


################################################################################################
###		MAIN	###

[n, box, r] = read_cnf_atoms ( 'cluster.inp', with_v=False )
#[n, box, r] = read_cnf_atoms ( 'two_spheres.inp', with_v=False )

print('n = ', n)
print('box size = ', box)
print('epsilon = ', eps)
print('well depth = ', -U)
print('')
print('')
print('intial overlap status:', overlap(box, r))

Energy_ = n_energy(box, r, eps, U)

print('initial energy = ', Energy_)


for i in range(int(N_iter)):
	Move_( r, ds )
	Energy_ = n_energy(box, r, eps, U)
	savetxt('data'+str(i)+'.csv', r, delimiter=',')
	print('Iteration', i, '\t', 'Overlap', overlap( box , r ), '\t', 'Energy = ', Energy_)
	if i%10 == 0:
		fig = plt.figure()
		
		
		# syntax for 3-D plotting
		ax = plt.axes(projection ='3d')
		for j in range(int(len(r))):
			# syntax for plotting
			ax.scatter(r[j][0], r[j][1], r[j][2])
		
		fig.savefig(str(i)+'.png', format="png", dpi=1200)
		plt.close()
		
	i += 1



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
parser.add_argument('--atom_config', nargs = 1)
#parser.add_argument('--box', nargs = 1)
args = parser.parse_args()
#print(args)
#n = int(args.n[0])
atom_config_file = str(args.atom_config[0])
ds = float(args.ds[0])			# ds is the step size in real space
N_iter = int(args.iter[0])	# no. of MC iteration steps
#box = float(args.box[0])


################################################################################################
###		MAIN	###

[n, box, r] = read_cnf_atoms ( atom_config_file, with_v=False )

print('n = ', n)
print('box size = ', box)
print('epsilon = ', eps)
print('well depth = ', -U)
print('')
print('')
print('intial overlap status:', overlap(box, r))

Energy_initial = n_energy(box, r, eps, U)

print('initial energy = ', Energy_initial)

En_ = []
iteration_ = []
for i in range(N_iter):
	Energy_intial = n_energy(box, r, eps, U)
	r_new = Move_( r, ds, box )
	Energy_new = n_energy(box, r_new, eps, U)
	del_E = Energy_new-Energy_initial
	mp = move_prob(del_E, T)
	r_ = np.random.rand()
	
	print('')
	print('Iteration = ', i,'del_E = ', del_E, 'prob = ', mp, 'rand = ', r_)
	
	En_.append(Energy_new)
	iteration_.append(i)
	
	if Energy_new < Energy_initial:
		#En_.append(Energy_new)
		#iteration_.append(i)
		Energy_initial = Energy_new
		r = r_new
	elif Energy_new > Energy_initial:
		if mp < r_:
			#En_.append(Energy_new)
			#iteration_.append(i)
			Energy_initial = Energy_new
			r = r_new
		#elif mp > r_:
			#En_.append(Energy_new)
			#iteration_.append(i)		
		#else:
		#	r = r_old
	#elif mp > r_:
		#r = r_old
		
		#np.savetxt('data'+str(i)+'.csv', r, delimiter=',')
		#print('Iteration', i, '\t', 'Overlap', overlap( box , r ), '\t', 'Energy = ', Energy_)
		
	#if i%10 == 0:
	
	f = open("energy_data.csv", "w")
	f.write("{},{}\n".format("iteration", "Energy"))
	for x in zip(iteration_, En_):
		f.write("{},{}\n".format(x[0], x[1]))
	f.close()
	
	fig, ax = plt.subplots()
	# syntax for 3-D plotting
	ax = plt.axes(projection ='3d')
	for j in range(int(len(r))):
		# syntax for plotting
		ax.scatter(r[j][0], r[j][1], r[j][2])
		#ax.set_xlim(0, box);	ax.set_ylim(0, box);	ax.set_zlim(0, box)
		ax.set_xlim(-box, box);	ax.set_ylim(-box, box);	ax.set_zlim(-box, box)
		
	fig.savefig(str(i)+'.png', format="png", dpi=100)
	plt.close()
	print('Iteration', i, '\t', 'Overlap', overlap( box , r ), '\t', 'Energy = ', Energy_new)
	i += 1

fig = plt.figure()
plt.plot(iteration_, En_, 's'+'-', color='blue')
plt.show()
#fig.savefig('Energy.png', format="png", dpi=100)
plt.close()

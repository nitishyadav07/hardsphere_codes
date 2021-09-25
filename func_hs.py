import numpy as np
from scipy import *  # to define zeros((N), dtype=int) array
from math import *
import os

from constants import *	# Physical constants

######################### Function definitions: ###########################################################

def Initialize_(n):
    "Random initialization of 'n' hard spheres array r in 1x1x1 box"
    r=[]
    for i in range(n):
	    r.append([rand(), rand(), rand()])
    return r

def Move_(r, ds):
	"flipping spin chosen (indexed t) from edge, after doing energy calc"
	index = int(len(r)*rand())
	dx = rand()*ds
	dy = rand()*ds
	dz = rand()*ds
	if (r[index][0] + dx) < 1 and (r[index][0] + dx) > 0:
		r[index][0] = r[index][0] + dx
	if (r[index][1] + dy) < 1 and (r[index][1] + dy) > 0:
		r[index][1] = r[index][1] + dy
	if (r[index][2] + dz) < 1 and (r[index][2] + dz) > 0:
		r[index][2] = r[index][2] + dz
	return r

####################################################################################
# OVERLAP FUNCTIONS:

def overlap ( box, r ):
    """Takes in box and coordinate array, and signals any overlap."""

    # Actual calculation is performed by function overlap_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in overlap'

    for i in range(n-1):
        if overlap_1 ( r[i,:], box, r[i+1:,:] ):
            return True # Immediate return on detection of overlap

    return False

def overlap_1 ( ri, box, r ):
    """Takes in coordinates of an atom and signals any overlap.

    Values of box and partner coordinate array are supplied.
    """

    import numpy as np

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # It is assumed that positions are in units where box = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in overlap_1'
    assert ri.size==3, 'Dimension error for ri in overlap_1'
    
    inv_box_sq = 1.0 / box ** 2
    
    for rj in r:
    	rij = ri - rj            # Separation vector
    	rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
    	rij_sq = np.sum(rij**2)  # Squared separation
    	if rij_sq < inv_box_sq:  # Check within cutoff
    		return True # Immediate return on detection of overlap
    	return False
'''
    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        return np.any(rij_sq<inv_box_sq)
    
    else:
        for rj in r:
            rij = ri - rj            # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq < inv_box_sq:  # Check within cutoff
                return True # Immediate return on detection of overlap
	
	return False
'''

	

def n_overlap ( box, r ):
    """Takes in box and coordinate array, and counts overlaps."""

    # This routine is used in the calculation of pressure
    # Actual calculation is performed by function n_overlap_1

    n, d = r.shape
    assert d==3, 'Dimension error for r in n_overlap'

    n_ovr = 0
    for i in range(n-1):
        n_ovr = n_ovr + n_overlap_1 ( r[i,:], box, r[i+1:,:] )

    return n_ovr

def n_overlap_1 ( ri, box, r ):
    """Takes in coordinates of an atom and counts overlaps.

    Values of box and partner coordinate array are supplied.
    """

    import numpy as np

    # In general, r will be a subset of the complete set of simulation coordinates
    # and none of its rows should be identical to ri

    # It is assumed that positions are in units where box = 1

    nj, d = r.shape
    assert d==3, 'Dimension error for r in n_overlap_1'
    assert ri.size==3, 'Dimension error for ri in n_overlap_1'

    inv_box_sq = 1.0 / box ** 2

    if fast:
        rij = ri - r                    # Get all separation vectors from partners
        rij = rij - np.rint(rij)        # Periodic boundary conditions in box=1 units
        rij_sq = np.sum(rij**2,axis=1)  # Squared separations
        n_ovr = np.count_nonzero(rij_sq<inv_box_sq)

    else:
        n_ovr = 0
        for rj in r:
            rij = ri - rj            # Separation vector
            rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
            rij_sq = np.sum(rij**2)  # Squared separation
            if rij_sq < inv_box_sq:  # Check within cutoff
                n_ovr = n_ovr + 1

    return n_ovr


####################################################################################
#		ENERGY CALCULATIONS:

def energy ( ri, box, rj, eps, U ):
	rij = ri - rj            # Separation vector
    rij = rij - np.rint(rij) # Periodic boundary conditions in box=1 units
    rij_sq = np.sum(rij**2)  # Squared separation
	if rij_sq > 1 + eps		 # The sphere dia is taken to be 1
		E = 0
	elif rij_sq < 1 + eps and rij_sq > 1
		E = -U
	return E
	
def n_energy ( box, r, eps, U ):
	for ri in r:
		for rj in r and rj != ri:
			E = energy ( ri, box, rj, eps, U)
	return E


####################################################################################
# 		Probability of move


def move_prob ( del_E, T ): # del_E is change in system energy upon proposed movement of single particle
	P = np.exp ( del_E/(K_B*T))
	return P


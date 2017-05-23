
# coding: utf-8

# In[ ]:

# read crd and store lipids in upper or lower membrane
from htmd import *
import numpy as np
from copy import deepcopy
import sys
import getopt

def read_crd(fname):
	from htmd.molecule.readers import Topology
	t = Topology()

	mol = Molecule()
	coords = list()
	nAtoms = 0
	current_resnum = 1
	up = []
	down = []

	with open(fname, 'r') as f:
		for line in f:
			pieces = line.split()
			if (not pieces[0].startswith("*") and len(pieces) > 7):
				if pieces[3] == 'P':
					if float(pieces[6]) > 0:
						up.append(pieces[1])
					else:
						down.append(pieces[1])
				if (int(pieces[1]) > current_resnum):
					mol._parseTopology(t,fname)
					mol.coords = deepcopy(np.array(coords, dtype=np.float32))
					current_resnum = int(pieces[1])
				t.resid.append(pieces[1])
				t.resname.append(pieces[2])
				t.name.append(pieces[3])
				coords.append([[float(pieces[4])],[float(pieces[5])],[float(pieces[6])]])
				t.segid.append(pieces[7])
				nAtoms = nAtoms + 1

	mol._parseTopology(t, fname)
	mol.coords = np.array(coords, dtype=np.float32)
	return(mol, up, down)


# In[ ]:

# example molecules
m1 = read_crd('dlpe_mini/step5_assembly.crd')
m1[0].center
m2 = read_crd('sopa_1.crd') # has only 1 residue
m2[0].center


# In[ ]:

# store coord of residue to be rm of m1
coordp_m1 = m1[0].get('coords', 'resid 1 and name P')
# store coord of residue to be inserted of m1
coordp_m2 = m2[0].get('coords', 'resid 1 and name P')
# translate residue to be inserted to the final coord
m2[0].moveBy(coordp_m1-coordp_m2)
# remove residue
m1[0].remove('resid 1')
# insert residue
m1[0].append(m2[0])


# In[ ]:

# random sampling (for future)
import random
population = m1[1] # population can be a list or a string
k = 3
rand = random.sample(population, k)

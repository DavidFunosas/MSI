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
				t.resid.append(pieces[1])
				t.resname.append(pieces[2])
				t.name.append(pieces[3])
				coords.append([[float(pieces[4])],[float(pieces[5])],[float(pieces[6])]])
				t.segid.append(pieces[7])
				nAtoms = nAtoms + 1

	mol._parseTopology(t, fname)
	mol.coords = np.array(coords, dtype=np.float32)
	mol.filter('not water')
	return(mol, up, down)

if __name__ == "__main__":
	m1 = read_crd('dlpe_n128.crd')
	print (m1[0])
	print (m1[1])
	print (m1[2])

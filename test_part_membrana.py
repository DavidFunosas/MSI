from htmd import *
import numpy
from copy import deepcopy
import sys
import getopt

def read_crd(fname):
	from htmd.molecule.readers import Topology
	t = Topology()

	mol_list = list()
	mol = Molecule()
	coords = list()
	nAtoms = 0
	current_resnum = 1
	with open(fname, 'r') as f:
		for line in f:
			pieces = line.split()
			if (not pieces[0].startswith("*") and len(pieces) > 7):
				if (int(pieces[1]) > current_resnum):		
					print(current_resnum)
					mol._parseTopology(t,fname)
					mol.coords = deepcopy(coords)
					mol_list.append(mol.copy())
					mol = Molecule()
					coords = list()
					current_resnum = int(pieces[1])
				t.resname.append(pieces[2])
				t.name.append(pieces[3])
				coords.append([[float(pieces[4])],[float(pieces[5])],[float(pieces[6])]])
				t.segid.append(pieces[7])
				nAtoms = nAtoms + 1

	mol._parseTopology(t, fname)
	mol.coords = coords
	mol_list.append(mol)
	return mol_list


def main(argv):
	if (len(argv) < 1):
		print ('An input .crd file must be specified:')
		print ('\ttest.py <inputfile>')
		sys.exit(2)
	else:
		if (argv[0].endswith('.crd')):
      			#mol = read_crd('./chl1_1.crd')
			mols = read_crd(argv[0])
			mol = Molecule()
			for m in mols:
				mol.append(m)

			#If it is a membrane (how do you know?)
			if (mol.numResidues > 10):
				max_z = max(mol.z)
				min_z = min(mol.z)
				average_z = (max_z+min_z)/2.0

				upper_mols = list()
				lower_mols = list()

				for res in range(0, len(mols)):
					if (abs(max(mols[res].z)-average_z) > abs(min(mols[res].z)-average_z)):
						upper_mols.append(mols[res])
					else:
						lower_mols.append(mols[res])

				print ('Upper molecules:')
				print (len(upper_mols))
				print ('Lower molecules:')
				print (len(lower_mols))









			
		elif (argv[0].endswith('.pdb')):
			mol = Molecule(argv[0])
		else:
			print ('Not supported file format')
			sys.exit(2)

		#print(mol)
		#print(mol.coords)	
		
		#If it is a membrane (how do you know?)
		if (mol.numResidues > 10):
			max_z = max(mol.z)
			min_z = min(mol.z)
			average_z = (max_z+min_z)/2.0

			upper_mols = list()
			lower_mols = list()

			for res in range(1,mol.numResidues+1):
				print(str(res))
				mc = mol.copy()
				mc.filter('resid ' + str(res))

				if (abs(max(mc.z)-average_z) > abs(min(mc.z)-average_z)):
					upper_mols.append(mc)
				else:
					lower_mols.append(mc)

			print ('Upper molecules:')
			print (len(upper_mols))
			print ('Lower molecules:')
			print (len(lower_mols))
	
		#mol.view(viewer='vmd')


if __name__ == "__main__":
   main(sys.argv[1:])


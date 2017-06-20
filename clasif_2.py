#!/usr/bin/python

# read crd and store lipids in upper or lower membrane

import numpy as np
from copy import deepcopy
import sys, os
import getopt
import random
import argparse
from lipids_charm_list import lipid_parent, parent_atom, lipidSurfArea, sqrSize
import datetime
from htmd import *

#It creates a layer with user defined size instead of the original size of the template by copying the content of m11copy until the layer reaches the desired size
def layer_creation(t, us, m11copy):
	if (t != us):
		value = int(us/t)
		times = value*value - 1

		x = 0.0
		y = 0.0

		c = 1

		m1 = m11copy.copy()

		for i in range(0, times):
			m1copy = m11copy.copy()
			m1copy.center

			if c == value:
				y = y + t
				x = 0.0
				m1copy.moveBy([[ x,    y,  0.0]])
				m1.append(m1copy)
				c = 0


			else: #In this step I set the squares in first column
				x += t
				m1copy.moveBy([[ x,    y,  0.0]])
				m1.append(m1copy)

			c += 1

		return(m1)
	else:
		return(m11copy)

# get 2 molecules = up and low
def get_leaflets(fname1, fname2):
	crd1 = read_crd(fname1)
	crd2 = read_crd(fname2)
	upper = deepcopy(crd1[0])
	if (len(crd1[1])) > 0:
		upper.filter('resid ' + ' '.join(crd1[1]))
	if (len(crd2[2])) > 0:
		lower = deepcopy(crd2[0])
	lower.filter('resid ' + ' '.join(crd2[2]))
	return [upper, lower]
	
# Convert crd file into molecule object
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
                if pieces[3] == parent_atom[lipid_parent[fname.split('/')[-1].strip('.crd')]]:
                    if float(pieces[6]) > 0:
                        up.append(str(pieces[1]))
                    else:
                        down.append(str(pieces[1]))
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
	
def proportions (ratios, lipids, max_ratio):
    ratios_s = deepcopy(ratios) #If I don't do deepcopy(), when sorting ratios_s, ratios is also sorted (???)
    ratios_s.sort()
    repl_per_lip = {}

    if len(ratios) == 1:
        #Go straight to the final membrane size calculations
        next
    elif len(ratios) == 2:
        prop = ratios_s[0]/(sum(ratios_s))
        total_number_repl = 30*prop
        del lipids[max_ratio]
        repl_per_lip[lipids[0]] = total_number_repl
    else:
        numerator = sum(ratios_s[0:-1])
        denominator = sum(ratios_s)
        prop = numerator/denominator
        total_number_repl = 30*prop

        del ratios[max_ratio]

        ratios_s = ratios

        ratios_s.sort()



        i = 0

        while len(ratios_s) >= 1:
            if len(ratios_s) == 1:
                repl_per_lip[lipids[i]] = total_number_repl - pr_lip*total_number_repl
                break

            pr_lip = ratios_s[i]/sum(ratios_s)
            repl_per_lip[lipids[i]] = pr_lip*total_number_repl
            del ratios_s[i]
            i += 1
    
    
    return(prop, total_number_repl, repl_per_lip)



def propDifferencesAcceptable(target_props, current_props, acceptable_difference_threshold):
	acceptable = True
	i = 0
	while (acceptable and i < len(target_props)):
		acceptable = target_props[i] < current_props[i] + acceptable_difference_threshold and target_props[i] > current_props[i] - acceptable_difference_threshold
		i += 1
		
	return acceptable
	
def propDifferences(target_props, current_props):
	differences = 0.0
	for i in range(len(target_props)):
		differences += abs(target_props[i] - current_props[i])
		
	return differences

def getLipidToInsert(target_props, current_props, prop_names):
	index = 0
	for i in range(len(target_props)):
		if (target_props[i] - current_props[i] > target_props[index] - current_props[index]):
			index = i
	return prop_names[index]

def getLipidToDelete(target_props, current_props, prop_names):
	index = 0
	for i in range(len(target_props)):
		if (target_props[i] - current_props[i] < target_props[index] - current_props[index]):
			index = i
	return prop_names[index]

def get_proportion_changes(before, after, lips):
	proportion_changes = []
	for lip in lips:
		after_copy = after.copy()
		after_copy.filter('resname ' + lip)
		before_copy = before.copy()
		before_copy.filter('resname ' + lip)
		n_lip_after = len(list(set(after_copy.get('resid'))))
		n_total_after = len(list(set(after.get('resid'))))
		n_lip_before = len(list(set(before_copy.get('resid'))))
		n_total_before = len(list(set(before.get('resid'))))
		proportion_changes.append((n_lip_after/n_total_after - n_lip_before/n_total_before)*100)
	
	return proportion_changes


###################################################################################################################################
####### Argument parsing

parser = argparse.ArgumentParser(description="Build custom lipid bilayers using CHARMM-GUI's database http://www.charmm-gui.org/?doc=input/membrane_only&step=1.")
parser.add_argument('-u', '--upper',
					dest = "up_lip",
					action = "store",
					type=str,
					nargs='+',
					default = None,
					required = True,
					help = "Upper-layer lipids composition (separated by space)",)
parser.add_argument('-l', '--lower',
					dest = "low_lip",
					action = "store",
					type=str,
					nargs='+',
					default = None,
					required = True,
					help = "Lower-layer lipids composition (separated by space)",)
parser.add_argument('-ur', "--upratio",
					dest = "up_ratio",
					action = "store",
					type=int,
					nargs='+',
					default = None,
					required = True,
					help = "Upper-layer lipids ratio (separated by space)",)
parser.add_argument('-lr', '--lowratio',
					dest = "low_ratio",
					action = "store",
					type=int,
					nargs='+',
					default = None,
					required = True,
					help = "Lower-layer lipids ratio (separated by space)",)
parser.add_argument('-s', '--size',
					dest = "sys_size",
					action = "store",
					type=int,
					default = 30,
					required = False,
					help = "Average bilayer surface area (in Angstrom).")
parser.add_argument('-d', "--distance",
					dest = "bet_dist",
					action = "store",
					type=int,
					default = 0,
					required = False,
					help = "Distance between layers (in Angstrom).",)
parser.add_argument('-o', "--out_name",
					dest = "out_name",
					action = "store",
					type=str,
					default = '',
					required = False,
					help = "Output name of resulting PDB file.",)
parser.add_argument('-v', "--viewer",
					dest = "viewer",
					action = "store_true",
					required = False,
					help = "Visualization of output in VMD turned on.",)

options = parser.parse_args()

###################################################################################################################################
###################################################################################################################################
####### Main module.


if __name__ == "__main__":

###################################################################################################################################
####### Storing the parse arguments into local variables.
	up_lip = [x.upper() for x in options.up_lip] # convert lipids to uppercase
	low_lip = [x.upper() for x in options.low_lip]
	up_ratio = options.up_ratio
	low_ratio = options.low_ratio
	sys_size= options.sys_size
	bet_dist = options.bet_dist
	out_name = options.out_name
	viewer = options.viewer
		

###################################################################################################################################
####### Checking if the number of lipids coincide with their respective ratio.
	
	if len(up_ratio) != len(up_lip) or len(low_ratio) != len(low_lip):
		print("\nERROR:\n")
		if len(up_ratio) != len(up_lip):
			print("The number of lipids in the upper argument must coincide with the upper ratio argument size.\n")
		if len(low_ratio) != len(low_lip):
			print("The number of lipids in the lower argument must coincide with the lower ratio argument size.\n")
		print("Please, review your arguments and run again.\n")
		sys.exit()

###################################################################################################################################
####### Checking if user-provided lipids are included in our database.

	for a in up_lip:
		if a.upper() not in lipid_parent.keys():
			print("Sorry, %s is not in our database." %(a))
			sys.exit()

	for a in low_lip:
		if a.upper() not in lipid_parent.keys():
			print("Sorry, %s is not in our database." %(a))
			sys.exit()


## Look for major lipids


	if (len(up_lip) > 0 and len(low_lip) > 0):


		#This is an example of input lipids and ratios from the user
		ARGGG = [[(up_lip[i], up_ratio[i]) for i in range(len(up_lip))], [(low_lip[i], low_ratio[i]) for i in range(len(low_lip))]]

		ARGGG[0].sort(key=lambda x: x[1], reverse=True)
		ARGGG[1].sort(key=lambda x: x[1], reverse=True)


		my_membrane1 = './db/' + ARGGG[0][0][0] + '.crd'
		my_membrane2 = './db/' + ARGGG[1][0][0] + '.crd'
		leaflets = get_leaflets(my_membrane1, my_membrane2) # templates de arriba y de abajo

		
		# Calculate proportions (current --> target) 
		prop_numbers = [[i[-1] for i in ARGGG[0]], [i[-1] for i in ARGGG[1]]]
		prop_names = [[i[0] for i in ARGGG[0]], [i[0] for i in ARGGG[1]]]


		for leaflet in range(len(leaflets)): # work separately with each leaflet per loop-iteration

			#It determines the target proportions, in percentage, of each lipid for the current leaflet
			sum_prop = sum(prop_numbers[leaflet])
			mol = leaflets[leaflet]
			mol.center()
	
			target_props = []
			for i in range(len(prop_numbers[leaflet])):
				target_props.append(round(prop_numbers[leaflet][i]*100.0/sum_prop,0)) # convert absolute ratios into percentages

			#As we start with the template membrane formed only by the lipid with the highest proportion, it just has a 100% proportion currently, and the 
			#rest of lipids have a 0% proportion.
			current_props = [100]

	
			for i in range(1,len(prop_numbers[leaflet])):
				current_props.append(0)

			#It copies the data before the insertion/deletion(s). This way, in case they don't result in proportions closer to the desired ones, it 
			#can reverse the changes.
			previous_props = deepcopy(current_props)
			previous_differences = propDifferences(target_props, current_props)
			leaflet_copy = deepcopy(mol)

			previous_nLipids = mol.numResidues
			current_nLipids = previous_nLipids

			#If the attempts are unsuccessful 10 times in a row, it exits the loop
			nFailedAttemptsInARow = 0

			#Acceptable difference in percentage with respect to the target proportions.
			acceptable_difference_threshold = 5
	
			while (not propDifferencesAcceptable(target_props, current_props, acceptable_difference_threshold) and nFailedAttemptsInARow < 10):
				#Insertion/Deletion
				#It gets the lipid type that we should remove, that is, the one whose proportion is the highest with respect to the target one.
				#It gets all distinct residue identifiers that can be replaced.
				#It will probably only replace the type of lipid that conforms the whole membrane at the beginning.
				resids = list(set(mol.get('resid', 'resname ' + getLipidToDelete(target_props, current_props, prop_names[leaflet]))))
				#It randomly selects the residue to remove from among the ones selected (they have all the same type).
				res_to_delete = str(random.choice(resids))
				#For the moment, there are no changes in proportions as it hasn't done neither insertions nor deletions.
				proportion_changes = [0]*len(current_props)

				#It chooses the lipid to insert, that is, the one whose proportion is the lowest with respect to the target one.
				resname_to_insert = getLipidToInsert(target_props, current_props, prop_names[leaflet])
				#It gets the template from which it will take the residue that will be inserted, and then it just takes the first residue, as we 
				#just want to add one in each iteration
				crd_to_insert = read_crd('./db/' + resname_to_insert + '.crd')
				mol_to_insert = crd_to_insert[0]
				mol_to_insert.filter('resid ' + str(crd_to_insert[leaflet+1][0]))
				mol_to_insert.resid = np.array([int(res_to_delete)] * len(mol_to_insert.resid))
				mol_to_insert.center()


				# store coord of residue to be removed from mol
				coordp_m1 = np.array([mol.get('coords', 'resid ' + res_to_delete +  ' and name ' + parent_atom[lipid_parent[mol.get('resname', 'resid ' + res_to_delete)[0]]])[0]])
				# store coord of residue to be inserted in mol
				coordp_m2 = np.array([mol_to_insert.get('coords', 'name ' + parent_atom[lipid_parent[mol_to_insert.get('resname')[0]]])[0]])

				# translate residue to be inserted to the final coord
				mol_to_insert.moveBy(coordp_m1-coordp_m2)
								
				#It inserts the molecule into the membrane
				deleted_resids = mol_to_insert.append(mol, collisions=True)
				mol = mol_to_insert.copy()	
				
				#Applying the changes in proportions
				current_nLipids -= len(deleted_resids)
				current_nLipids += 1
				proportion_changes = get_proportion_changes(leaflet_copy, mol, prop_names[leaflet])
				current_props = [current_props[i]+proportion_changes[i] for i in range(len(current_props))]		

				#If proportions have improved, it preserves the changes
				current_differences = propDifferences(target_props, current_props)
				if (current_differences < previous_differences):
					previous_props = deepcopy(current_props)
					previous_differences = deepcopy(current_differences)
					leaflet_copy = deepcopy(mol)
					previous_nLipids = current_nLipids
					nFailedAttemptsInARow = 0
				#If proportions have not improved, it reverts the changes
				else:
					mol = leaflet_copy
					current_props = previous_props
					current_differences = previous_differences
					nFailedAttemptsInARow += 1
					current_nLipids = previous_nLipids

			leaflets[leaflet] = mol
	

		#The size we use in all our templates
		template_size_u = float(sqrSize[ARGGG[0][0][0]])
		template_size_l = float(sqrSize[ARGGG[1][0][0]])


		#It creates two layers with the size the user has indicated from the ones I already have
		uplayer = layer_creation(template_size_u, sys_size, leaflets[0])
		lowlayer = layer_creation(template_size_l, sys_size, leaflets[1])

		
		uplayer.center
		lowlayer.center
		
		#Calculating the max and min z coordinates of the lower leaflet to calculate afterwards the minimum distance (0)
		#between the layers
		zmax = max(lowlayer.z)
		zmin = min(lowlayer.z)

		#It moves the lower layer so as to have them separated by the distance (in the Z axis) defined by the user
		lowlayer.moveBy([[ 0.0,    0.0,   -(bet_dist+zmax-zmin)]])

		#It joins the two layers in a single molecule object and it writes it in the resulting .pdb file.
		uplayer.append(lowlayer)
	
		#uplayer.write(out_name)
		if (out_name is ''):
			out_name = 'charmmoutput_' + str(datetime.datetime.now()) + '.pdb'
		uplayer.write(out_name)
		print ('\nAwesome! The PDB file \"' + out_name + '\" has been created successfully!')


		if viewer:
			uplayer.view('VMD')
			input("\nPress enter to finish visualization and exit.")
		





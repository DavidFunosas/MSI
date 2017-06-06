
# coding: utf-8

# In[1]:

# read crd and store lipids in upper or lower membrane
from htmd import *
import numpy as np
from copy import deepcopy
import sys, os
import getopt
import random
import argparse
from lipids_charm_list import lipid_parent, parent_atom, lipidSurfArea, sqrSize
import datetime



def get_leaflets(fname1, fname2):
	crd1 = read_crd(fname1)
	crd2 = read_crd(fname2)
	upper = deepcopy(crd1[0])
	if (len(crd1[1])) > 0:
		upper.filter('resid ' + ' '.join(crd1[1]))
	if (len(crd2[2])) > 0:
		lower = deepcopy(crd2[0])
	lower.filter('resid ' + ' '.join(crd2[2]))
	return (upper, lower)
	

#Defining data retrieval function
def data_retrieval(ARG11):
	#Initialize some variables
	lipids = []
	#surfaces = []
	ratios = []

	for i in range(0, len(ARG11)): #Cambiar o de len
		lipids.append(ARG11[i][0])
		ratios.append(ARG11[i][1])
		#surfaces.append(float(d[ARG11[i][1]].rstrip()))

	#return(ratios, lipids, surfaces)
	return(ratios, lipids)
	
	
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
            print(ratios_s)
    
    
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



###################################################################################################################################
####### Argument parsing

parser = argparse.ArgumentParser(description="This program does a superimposition of two proteins based on Secondary Structure Elements (SEEs).")
parser.add_argument('-u', '--upper',
					dest = "up_lip",
					action = "store",
					type=str,
					nargs='+',
					default = None,
					required = True,
					help = "Upper-layer lipids. (Separated by spcae)",)
parser.add_argument('-l', '--lower',
					dest = "low_lip",
					action = "store",
					type=str,
					nargs='+',
					default = None,
					required = True,
					help = "Lower-layer lipids. (Separated by space)",)
parser.add_argument('-ur', "--upratio",
					dest = "up_ratio",
					action = "store",
					type=int,
					nargs='+',
					default = None,
					required = True,
					help = "Upper-layer lipid ratio. (Separated by space)",
					)
parser.add_argument('-lr', '--lowratio',
					dest = "low_ratio",
					action = "store",
					type=int,
					nargs='+',
					default = None,
					required = True,
					help = "Lower-layer lipid ratio. (Separated by space)",)
parser.add_argument('-s', '--size',
					dest = "sys_size",
					action = "store",
					type=int,
					nargs='+',
					default = 200,
					required = False,
					help = "Average bilayer surface area. (In Armstrong).",)
parser.add_argument('-d', "--distance",
					dest = "bet_dist",
					action = "store",
					type=int,
					default = 3,
					required = False,
					help = "Distance between layers. (In Armstrong).",
					)
parser.add_argument('-o', "--out_name",
					dest = "out_name",
					action = "store",
					type=str,
					default = '',
					required = False,
					help = "Output name of resulting PDB file.",
					)

options = parser.parse_args()

###################################################################################################################################
###################################################################################################################################
####### Main module.


if __name__ == "__main__":

###################################################################################################################################
####### Storing the parse arguments into local variables.
	up_lip = [x.upper() for x in options.up_lip]
	low_lip = [x.upper() for x in options.low_lip]
	up_ratio = options.up_ratio
	low_ratio = options.low_ratio
	sys_size= options.sys_size
	bet_dist = options.bet_dist
	out_name = options.out_name

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

	for a in low_lip:
		if a.upper() not in lipid_parent.keys():
			print("Sorry, %s is not in our database." %(a))


	if (len(up_lip) > 0 and len(low_lip) > 0):
		my_membrane1 = './db/' + up_lip[0] + '.crd'
		my_membrane2 = './db/' + low_lip[0] + '.crd'
		leaflets = get_leaflets(my_membrane1, my_membrane2)

		#Loading dictionary of surface areas per lipid
		d = {}
		with open("lipid_surface_areas.txt") as f:
			for line in f:
				(key,val) = line.split("\t")
				d[key] = val

		#This is an example of input lipids and ratios from the user
		#ARG11 =["DGPA", 3, "DYPA", 1]
		ARGGG = [[(up_lip[i], up_ratio[i]) for i in range(len(up_lip))], [(low_lip[i], low_ratio[i]) for i in range(len(low_lip))]]
		n = 30 #For dimensions of membrane template
		m = 100 #For dimensions from user

		#Retrieving information
		upper = data_retrieval(ARGGG[0])
		lower = data_retrieval(ARGGG[1])


		# In[80]:

		#Getting the index of the major lipid type and defining template membrane
		max_ratio_u = max(enumerate(upper[0]), key = lambda x: x[1])
		max_ratio_l = max(enumerate(lower[0]), key = lambda x: x[1])

		prop_numbers = [up_ratio, low_ratio]
		prop_names = [up_lip, low_lip]

		for leaflet in range(len(leaflets)):
			sum_prop = sum(prop_numbers[leaflet])
			mol = leaflets[leaflet]
			mol.center()
	
			target_props = []
			for i in range(len(prop_numbers[leaflet])):
				target_props.append(round(prop_numbers[leaflet][i]*100.0/sum_prop,0))

			#The first lipid is the most abundant one
			current_props = [100]
	
			for i in range(1,len(prop_numbers[leaflet])):
				current_props.append(0)
		
			#It copies the data before the insertion/deletion(s)
			previous_props = deepcopy(current_props)
			previous_differences = propDifferences(target_props, current_props)
			leaflet_copy = deepcopy(mol)

			#If the attempts are unsuccessful 10 times in a row, it exits the loop
			nFailedAttemptsInARow = 0

			acceptable_difference_threshold = 5
	
			while (not propDifferencesAcceptable(target_props, current_props, acceptable_difference_threshold) and nFailedAttemptsInARow < 10):
				#Insertion/Deletion
				#It gets all distinct residue identifiers that can be replaced
				#It will probably only replace the type of lipid that conforms the whole membrane at the beginning
				resids = list(set(mol.get('resid', 'resname ' + getLipidToDelete(target_props, current_props, prop_names[leaflet]))))
				#It randomly selects the residue to remove
				res_to_delete = str(random.choice(resids))
				proportion_changes = [0]*len(current_props)

				#It chooses the lipid to insert
				resname_to_insert = getLipidToInsert(target_props, current_props, prop_names[leaflet])
				mol_to_insert = read_crd('./db/' + resname_to_insert + '.crd')[0]
				mol_to_insert.filter('resid 1')
				mol_to_insert.center()


				# store coord of residue to be removed from mol
				coordp_m1 = mol.get('coords', 'resid ' + res_to_delete +  ' and name ' + parent_atom[lipid_parent[mol.get('resname', 'resid ' + res_to_delete)[0]]])
				# store coord of residue to be inserted in mol
				coordp_m2 = mol_to_insert.get('coords', 'name ' + parent_atom[lipid_parent[mol_to_insert.get('resname')[0]]])
				# translate residue to be inserted to the final coord
				mol_to_insert.moveBy(coordp_m1-coordp_m2)

				#It inserts the molecule into the membrane
				deleted_resids = mol.append(mol_to_insert)		

				#Applying the changes in proportions
				for res_id in deleted_resids:
					proportion_changes[prop_names[leaflet].index(leaflet_copy.get('resname', 'resid ' + res_id))] -= 1

				current_props = [current_props[i]+proportion_changes[i] for i in range(len(current_props))]		

				#If proportions have improved, it preserves the changes
				current_differences = propDifferences(target_props, current_props)
				if (current_differences < previous_differences):
					previous_props = deepcopy(current_props)
					previous_differences = deepcopy(current_differences)
					leaflet_copy = deepcopy(mol)
					nFailedAttemptsInARow = 0
				#If proportions have not improved, it reverts the changes
				else:
					mol = leaflet_copy
					current_props = previous_props
					current_differences = previous_differences
					nFailedAttemptsInARow += 1
	

		leaflets.view('VMD')
		if (out_name is None):
			leaflets.write('output - ' + str(datetime.datetime.now()) + '.pdb')
		else:
			leaflets.write(out_name)
		print ('The PDB file has been created successfully!')
		



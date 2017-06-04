
# coding: utf-8

# In[1]:

# read crd and store lipids in upper or lower membrane
from htmd import *
import numpy as np
from copy import deepcopy
import sys
import getopt
import random

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
	
def get_leaflets(fname):
	crd = read_crd(fname)
	upper = deepcopy(crd[0])
	upper.filter('resid ' + 'resid '.join(crd[1]))
	lower = deepcopy(crd[0])
	lower.filter('resid ' + 'resid '.join(crd[2]))
	return (upper, lower)
	
crd = 'dlpe_n128.crd'
leaflets = get_leaflets(crd)


# In[3]:

#Loading dictionary of surface areas per lipid
d = {}
with open("lipid_surface_areas.txt") as f:
    for line in f:
        (key,val) = line.split("\t")
        d[key] = val



# In[79]:


#This is an example of input lipids and ratios from the user
#ARG11 =["DGPA", 3, "DYPA", 1]
ARGGG = [["DGPA", 3, "DYPA", 1, "LAU", 1], ["DGPA", 7, "LAU", 1]] #First item always upper leaflet
n = 30 #For dimensions of membrane template
m = 100 #For dimensions from user


#Defining data retrieval function
def data_retrieval(ARG11):
    #Initialize some variables
    lipids = []
    surfaces = []
    ratios = []
    
    for i in range(0, len(ARG11)): #Cambiar o de len
        if (type(ARG11[i]) == int):
            ratios.append(ARG11[i])
        elif (type(ARG11[i]) != int):
            lipids.append(ARG11[i])
            surfaces.append(float(d[ARG11[i]].rstrip()))
            
    return(ratios, lipids, surfaces)

#Retrieving information
upper = data_retrieval(ARGGG[0])
lower = data_retrieval(ARGGG[1])



# In[80]:

#Getting the index of the major lipid type and defining template membrane
max_ratio_u = max(enumerate(upper[0]), key = lambda x: x[1])
max_ratio_l = max(enumerate(lower[0]), key = lambda x: x[1])

my_membrane1 = upper[1][max_ratio_u[0]] + '.crd'
my_membrane2 = lower[1][max_ratio_l[0]] + '.crd'

#print(my_membrane1)
#if membrane1 == my_membrane2:
#Put something here? and avoid the process of getting two different membrane templates, splitting them and ensambling them.
'''
if max_ratio_u[1] >= max_ratio_l[1]:
    max_ratio = max_ratio_u[0] #This is an index
    my_membrane = upper[1][max_ratio] + '.crd'
    main = lower[1][max_ratio]
    del lipid_number[0][max_ratio]
    del upper[1][max_ratio]
    
else:
    max_ratio = max_ratio_l[0] #This is an index)
    my_membrane = lower[1][max_ratio] + '.crd'
    main = upper[1][max_ratio]
    del lipid_number[1][max_ratio]
    del lower[1][max_ratio]
'''


# In[166]:
'''
#Number of replacements required per lipid type (for upper and lower leaflets separatedly)
def lipid_numbers(n, ratios, surfaces):
    from operator import mul

    smemb = n*n  #Surface of our template membrane
    smin = sum(map(mul, ratios, surfaces)) #Minimum area of the lipid set
    prop = int(smemb/smin)
    lipid_number = [i * prop for i in ratios] #List with number of lipids per type that needs to be replaced   
      
    return(lipid_number)


lipid_number = [lipid_numbers(n, upper[0], upper[2]), lipid_numbers(n, lower[0], lower[2])]
total_replacements = [sum(lipid_number[0]), sum(lipid_number[1])]
print(lipid_number[0])
'''

# In[59]:

# Every membrane template has 30 lipids
# When calculating number of lipids to replace, take this into account

#1:1  ->  1/2 lipids must be replaced
#1:2  ->  1/3 lipids must be replaced
#1:3:5 -> (3+1)/9 lipids must be replaced

#ARGGG = [["DGPA", 3, "DYPA", 1, "LAU", 1], ["DGPA", 7, "LAU", 1]]

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

up = proportions(upper[0], upper[1], max_ratio_u[0])
down = proportions(lower[0], lower[1], max_ratio_l[0])

#WHAT TO DO HERE WHEN THE NUMBER HAS DECIMALS?

print(up[0], up[2])
print(down[0], down[2])

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

#Number of lipids per leaflet: 30

for leaflet in range(len(leaflets)):
	sum_prop = sum(prop_numbers)
	mol = leaflets[leaflet]
	
	target_props = []
	for i in range(len(prop_numbers)):
		target_props.append(round(prop_numbers[i]*100.0/sum_prop,0))

	#The first lipid is the most abundant one
	current_props = [100]
	
	for i in range(1,len(prop_numbers)):
		current_props.append(0)
		
	#It copies the data before the insertion/deletion(s)
	previous_props = deepcopy(current_props)
	previous_differences = propDifferences(target_props, current_props)
	leaflet_copy = deepcopy(mol)

	#If the attempts are unsuccessful 10 times in a row, it exits the loop
	nFailedAttempsInARow = 0

	acceptable_difference_threshold = 5
	
	while (not propDifferencesAcceptable(target_props, current_props, acceptable_difference_threshold) and nFailedAttemptsInARow < 10):
		#Insertion/Deletion
		#It gets all distinct residue identifiers that can be replaced
		#It will probably only replace the type of lipid that conforms the whole membrane at the beginning
		resids = list(set(mol.get('resid', 'resname ' + getLipidToDelete(target_props, current_props, prop_names))))
		#It randomly selects the residue to remove
		res_to_delete = choice(resids)
		proportion_changes = [0]*len(current_props)

		#It chooses the lipid to insert
		resname_to_insert = getLipidToInsert(target_props, current_props, prop_names)
		mol_to_insert = read_crd(resname_to_insert + '.crd')[0]


		# store coord of residue to be removed from mol
		coordp_m1 = mol[0].get('coords', 'resid ' + res_to_delete)
		# store coord of residue to be inserted in mol
		coordp_m2 = mol_to_insert[0].get('coords')
		# translate residue to be inserted to the final coord
		mol_to_insert[0].moveBy(coordp_m1-coordp_m2)

		#It inserts the molecule into the membrane
		deleted_resids = mol.append(mol_to_insert)		

		#Applying the changes in proportions
		for res_id in deleted_resids:
			proportion_changes[prop_names.index(leaflet_copy.get('resname', 'resid ' + res_id))] -= 1

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



#k = sum(up[1], )

'''
# In[105]:



# In[73]:

#Removing values that correspond to the membrane type selected (it won't need replacement of this lipid type)
#del lipid_number[max_ratio]
#del lipids[max_ratio]

#Reading molecule membrane 1
m1 = read_crd(my_membrane)
m1[0].center

# random sampling (for future)
import random
population = m1[1] # population can be a list or a string

rand = []
for i in range(0, len(total_replacements)):
    rand1 = random.sample(m1[i+1], total_replacements[i])
    print(rand1)

rand_copy = rand


# In[ ]:


big_list = []
lip_repl = {}

if len(lipid_number) >= 1: #CHANGE THIS IF HERE AND ADAPT TO THE DIFFERENT SITUATIONS THAT WE WANT
    for i in range(0, len(lipid_number)):
        lip_repl[lipids[i]] = rand_copy[0:lipid_number[i]]
        del rand_copy[0:lipid_number[i]]

#for i in range(0, lipid_number-1):
#    small_list = rand[0:lipid_number[i]]
#    big_list.append([rand[]])
    
for key in lip_repl:
    print(key, lip_repl[key])


# In[113]:

#Reading molecules membrane (lipid donors)

m1_sa = d[my_membrane.split('.')[0]] #Surface area of main membrane lipid. Use this to compare before removing

for key in lip_repl:

    m2 = read_crd(key + '.crd')#Selecting membrane donor
    m2[0].filter('resid 1')
    
    m2_sa = d[key] #Surface area of m2 type lipid
        
    for value in lip_repl[key]:
        if m2_sa >= m1_sa:
            m2[0].center
            #m2[0].filter('resid 1')
            coordp_m2 = m2[0].get('coords', 'resid 1 and name P')#Selecting lipid
            # store coord of residue to be rm of m1
            residue = 'resid ' + str(value) + ' and name P'
            coordp_m1 = m1[0].get('coords', residue)

            res_short = 'resid ' + str(value)#Short string without atom specifications
            m1[0].remove(res_short)
            m2[0].moveBy(coordp_m1-coordp_m2)
            m1[0].append(m2[0])

        else:
            print("Not ready yet")
            
m1[0].view('VMD')


# In[ ]:
'''



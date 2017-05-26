
# coding: utf-8

# 

# In[23]:

#Load dictionary of surface areas per lipid
d = {}
with open("lipid_surface_areas.txt") as f:
	for line in f:
		(key,val) = line.split("\t")
		d[key] = val


#This is an example of input lipids and ratios from the user
ARG11 =["DGPA", 3, "DYPA", 1]
n = 30 #For dimensions


#Initialize some variables
lipids = []
ratios = []
surfaces = []

#Sort ARGV values into lipids or ratios and get the respective surface areas.
for i in range(0, len(ARG11)): #Cambiar o de len
	if (type(ARG11[i]) == int):
		ratios.append(ARG11[i])
	elif(type(ARG11) != int):
		lipids.append(ARG11[i])
		surfaces.append(float(d[ARG11[i]].rstrip()))

#Getting the main lipid type (higher in number) and select the membrane type from available ones
max_ratio = max(enumerate(ratios), key = lambda x: x[1])[0] #This is the index of the max value

for keys in d:
	if (keys == lipids[max_ratio]): #If more than one lipid with same ratio, it takes the first that appears
		my_membrane = keys + '.crd'
		break
        
#Calculation of number of lipids per type that will be necessary to replace
from operator import mul

sbig = n*n
ssmall = sum(map(mul, ratios, surfaces))
prop = int(sbig/ssmall)
lipid_number = [i * prop for i in ratios] #List with number of lipids per type that needs to be replace

#Removing values that correspond to the membrane type selected (it won't need replacement of this lipid tipe)
del lipid_number[max_ratio]
del lipids[max_ratio]

total_replacements = sum(lipid_number)   # THIS WILL BE THE K IN COTE'S .PY FILE

print(len(d))
print("ssmall =" + str(ssmall))
print("sbig =" + str(sbig))
print(prop)
print(lipid_number)
print(total_replacements)


# In[24]:

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


# In[25]:

# example molecules
m1 = read_crd('DGPA.crd')
m1[0].center
#m2 = read_crd('DYPA.crd') # Second
#m2[0].center


# In[8]:

#m1[0].view(viewer='VMD')


# In[28]:

# random sampling (for future)
import random
population = m1[1] # population can be a list or a string
#k = 3

rand = random.sample(population, total_replacements)
rand_copy = rand

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


# In[29]:

#Reading molecule membrane 1
m1 = read_crd(my_membrane)
m1[0].center

#Reading molecules membrane (lipid donors)
for key in lip_repl:
    #name = key + '.crd'
    m2 = read_crd(key + '.crd')#Selecting membrane donor
    m2[0].center
    coordp_m2 = m2[0].get('coords', 'resid 1 and name P')#Selecting lipid
    
    for value in lip_repl[key]:
        # store coord of residue to be rm of m1
        residue = 'resid ' + str(value) + ' and name P'
        residue_short = 'resid ' + str(value)#Short string without atom specifications
        coordp_m1 = m1[0].get('coords', residue)
        print(coordp_m1)
        
        #AQUÍ HABRÍA QUE INTEGRAR LA PARTE DE REMOVE AND INSERT
        #2[0].moveBy(coordp_m1-coordp_m2)
        # remove residue
        #1[0].remove(residue_short)
        # insert residue
        #1[0].append(m2[0])


# store coord of residue to be rm of m1
#coordp_m1 = m1[0].get('coords', residue)
#print(coordp_m1)
#coordp_m1 = m1[0].get('coords', 'resid 1 and name P')
#print(coordp_m1)
# store coord of residue to be inserted of m1
#coordp_m2 = m2[0].get('coords', 'resid 1 and name P')
# translate residue to be inserted to the final coord
#print(coordp_m2)





#m1[0].remove('resid 1 2 3')


# In[ ]:




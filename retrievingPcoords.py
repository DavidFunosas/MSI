from htmd import *

# Load one molecule, for instance 5JZQ. It has two P in resid 2.

mol = Molecule('5JZQ')

# Check which resid have P

mol.get('resid', sel='name P')

# Fetch the coordinates of both P in resid 2.

a = mol.get('coords', sel='name P and resid 2')

# Get the coordinates of only one of them.

lonelyP = a [0]

# Save the Z coordinate into a variable

pZcoord = lonelyP[2]


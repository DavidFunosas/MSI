#!/usr/bin/python

###################################################################################################################################
####### Importing the required modules.

import sys, os
import argparse
import numpy as np
from lipids_charm_list import lipid_parent, parent_atom, lipidSurfArea, sqrSize
import datetime

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
					type=int,,
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
		if a.upper() not in daug_parent.keys():
			print("Sorry, %s is not in our database." %(a))

	for a in low_lip:
		if a.upper() not in daug_parent.keys():
			print("Sorry, %s is not in our database." %(a))




###################################################################################################################################
####### Printing some things out.
		

	print(up_lip)
	uplist = up_lip
	print(len(uplist))
	print(low_lip)
	lowlist = low_lip	
	print(len(lowlist))
	print(up_ratio)	
	print(low_ratio)
	print(sys_size)
	print(bet_dist)
	parent_atom.keys()
	daug_parent.keys()

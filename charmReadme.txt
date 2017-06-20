#################################################################################################
#           				 CHARMM-ANDER            	                        #
#                                                                                               #
# Authors:  Maria José Falaguera, David Funosas, Lucía Rodríguez Vázquez & Edgar Sánchez Prados #                              
# Goal:     Build custom lipid bilayers using CHARMM-GUI's database                             #
#           (http://www.charmm-gui.org/?doc=input/membrane_only&step=1)                         # 
#                                                                                               #
#################################################################################################


0.REQUIREMENTS

The program has been tested on Linux (Ubuntu distro). A working python3 environment is required 
to run the scripts that are provided to automate the calculations and analysis.
The python3 environment must include the following modules: sys, os, numpy, copy, getopt, 
argparse and htmd.
User must first uncompress the db.zip file in the root directory of the script. 

1.SUMMARY

CHARMM-ANDER is a third-party-package able to build custom lipidic bilayers, either heterogeneous
or homogeneous. It is intended to reproduce the behaviour of CHARMM-GUI'S Membrane Builder by
command line.
The program needs a series of user-specified input arguments that include the lipidic 
composition of the membrane, the abundance ratio of each lipid type, the membrane size and 
distance between layers. 


2.USAGE

usage: clasif_2.py [-h] -u UP_LIP [UP_LIP ...] -l LOW_LIP [LOW_LIP ...] -ur
                   UP_RATIO [UP_RATIO ...] -lr LOW_RATIO [LOW_RATIO ...]
                   [-s SYS_SIZE] [-d BET_DIST] [-o OUT_NAME] [-v]

Build custom lipid bilayers using CHARMM-GUI's database http://www.charmm-
gui.org/?doc=input/membrane_only&step=1.

optional arguments:
  -h, --help            show this help message and exit
  -u UP_LIP [UP_LIP ...], --upper UP_LIP [UP_LIP ...]
                        Upper-layer lipids composition (separated by space)
  -l LOW_LIP [LOW_LIP ...], --lower LOW_LIP [LOW_LIP ...]
                        Lower-layer lipids composition (separated by space)
  -ur UP_RATIO [UP_RATIO ...], --upratio UP_RATIO [UP_RATIO ...]
                        Upper-layer lipids ratio (separated by space)
  -lr LOW_RATIO [LOW_RATIO ...], --lowratio LOW_RATIO [LOW_RATIO ...]
                        Lower-layer lipids ratio (separated by space)
  -s SYS_SIZE, --size SYS_SIZE
                        Average bilayer surface area (in Angstrom).
  -d BET_DIST, --distance BET_DIST
                        Distance between layers (in Angstrom).
  -o OUT_NAME, --out_name OUT_NAME
                        Output name of resulting PDB file.
  -v, --viewer          Visualization of output in VMD turned on.


3.EXAMPLE

The following example command will generate a lipid bilayer having SAPA and SAPC on the upper layer with a ratio 1 to 3 and SAPE and SAPG on the upper layer with a ratio 3 to 1. Distance between layers will be set to 

python charmm.py -u SAPA SAPC -l SAPE SAPG -ur 1 3 -lr 3 1 -d 10 -o example.pdb -v



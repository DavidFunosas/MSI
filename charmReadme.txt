			      _.--""`-..
			    ,'          `.
			  ,'          __  `.
			 /|          " __   \
			, |           / |.   .
			|,'          !_.'|   |
		      ,'             '   |   |
		     /              |`--'|   |
		    |                `---'   |
		     .   ,                   |                       ,".
		      ._     '           _'  |                    , ' \ `
		  `.. `.`-...___,...---""    |       __,.        ,`"   L,|
		  |, `- .`._        _,-,.'   .  __.-'-. /        .   ,    \
		-:..     `. `-..--_.,.<       `"      / `.        `-/ |   .
		  `,         """"'     `.              ,'         |   |  ',,
		    `.      '            '            /          '    |'. |/
		      `.   |              \       _,-'           |       ''
			`._'               \   '"\                .      |
			   |                '     \                `._  ,'
			   |                 '     \                 .'|
			   |                 .      \                | |
			   |                 |       L              ,' |
			   `                 |       |             /   '
#################################################################################################
#           				 CHARMM-ANDER            	                        #
#                                                                                               #
# Authors:  Maria José Falaguera, David Funosas, Lucía Rodríguez Vázquez & Edgar Sánchez Prados #                              
# Goal:     Build custom lipid bilayers using CHARMM-GUI's database                             #
#           http://www.charmm-gui.org/?doc=input/membrane_only&step=1.").                       # 
#                                                                                               #
#################################################################################################


0.REQUIREMENTS

The program has been tested on Linux (Ubuntu distro). A working python3 environment is required 
to run the scripts that are provided to automate the calculations and analysis.
The python3 environment must include the following modules: sys, os, numpy, copy, getopt, 
argparse and htmd.

1.SUMMARY

CHARMM-ANDER is a third-party-package able to build custom lipidic bilayers, either heterogeneous
or homogeneous. It is intended to reproduce the behaviour of CHARMM-GUI'S Membrane Builder by
command line.
The program needs a series of user-specified input arguments that include the lipidic 
composition of the membrane, the abundance ratio of each lipid type, the membrane size and 
distance between layers. 


2.USAGE

usage: clasif.py [-h] -u UP_LIP [UP_LIP ...] -l LOW_LIP [LOW_LIP ...] -ur
                 UP_RATIO [UP_RATIO ...] -lr LOW_RATIO [LOW_RATIO ...]
                 [-s SYS_SIZE [SYS_SIZE ...]] [-d BET_DIST] [-o OUT_NAME] [-v]

  -h, --help            show this help message and exit
  -u UP_LIP [UP_LIP ...], --upper UP_LIP [UP_LIP ...]
                        Upper-layer lipids. (Separated by spcae)
  -l LOW_LIP [LOW_LIP ...], --lower LOW_LIP [LOW_LIP ...]
                        Lower-layer lipids. (Separated by space)
  -ur UP_RATIO [UP_RATIO ...], --upratio UP_RATIO [UP_RATIO ...]
                        Upper-layer lipid ratio. (Separated by space)
  -lr LOW_RATIO [LOW_RATIO ...], --lowratio LOW_RATIO [LOW_RATIO ...]
                        Lower-layer lipid ratio. (Separated by space)
  -s SYS_SIZE [SYS_SIZE ...], --size SYS_SIZE [SYS_SIZE ...]
                        Average bilayer surface area. (In Armstrong).
  -d BET_DIST, --distance BET_DIST
                        Distance between layers. (In Armstrong).
  -o OUT_NAME, --out_name OUT_NAME
                        Output name of resulting PDB file.
  -v, --viewer          Visualization of output in VMD turned on.


3.EXAMPLE




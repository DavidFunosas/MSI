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
#           				 Charmm-Ander           	                        #
#                                                                                               #
# Authors:  Maria José Falaguera, David Funosas, Lucía Rodríguez Vázquez & Edgar Sánchez Prados #                              
# Goal:    The script allows the user to build custom lipid bilayers using CHARMM-GUI's         #
#          database.                                                                            #
#################################################################################################


0.REQUIREMENTS
The program has been tested on Linux (Ubuntu distro). A working python3 environment is required 
to run the scripts that are provided to automate the calculations and analysis.
The python3 environment must include the following modules: sys, os, numpy, copy, getopt, 
argparse and htmd.


1.SUMMARY
Charmm-Ander is a third-party-package able to build custom lipidic bilayers, either heterogeneous
or homogeneous. It is intended to reproduce the behaviour of CHARMM-GUI'S Membrane Builder by
command line.
The program needs a series of user-specified input arguments that include the lipidic 
composition of the membrane, the abundance ratio of each lipid type, the membrane size and 
distance between layers. 

2.INSTALLATION
#The package does not require any special installation procedure. The compressed
#folder must be downloaded and uncompressed and placed anywhere. In the main
#folder type:

#sudo python3 setup.py install

#With this command the databases needed to run the scripts will be stored in
#python containers and the modules needed for the analysis will be placed
#in the python3 library.


3.DATABASE FILES AND MODULES
There are four files extracted from CoDNaS, CATH and CSA databases from theirs
servers, parsed and stored as .csv and .tsv files (stored in "database/" directory).

-> codnas_headed.tsv	Contains one line per each conformer found in CoDNaS with its
			information. e.g.:

PDB Method Temperature PH Mutation Length Source Ligands Ligands_Biolip CATH Pool_repr
11AS_A XRAYDIFFRACTION 293 7.5 1 330 EscherichiacoliK12 ASN ASN 11asA00 11AS_A


-> cluster.csv	Contains in each line the group of conformers coming from the
		same protein structure e.g.:

11AS_A,11AS_B,12AS_A,12AS_B


-> csa.tsv		Contains in each line one active site residue of a CoDNaS
			conformer entry. e.g.:

11AS_A 46
11AS_A 100
11AS_A 116

-> cath_codes.csv	Contains in each line a CATH ID from CoDNaS and its CATH code
			e.g.:

11asA00,3.30.930.10

-> cath_nodes.tsv	Contains in each line a C.A.T.H. node and its description. e.g.:

3.30.930.10	Bira Bifunctional Protein; Domain 2

During installation these files will be stored in python containers that are going to be
persistently stored for future uses (stored in "database/" directory):

-> codnas_dic.p
-> clusters_list.p
-> csa_dic.p
-> cath_codes_dic.p
-> cath_nodes_dic.p


4.USAGE WITH ccas SCRIPT

Input files, ccas executable and database/ directory should be in the same folder.
If you want to execute the program in a folder that is not CCAS1.0, you must copy-paste
the executable and database/ directory to that folder. All the output files will be
created in that folder.

An example of the most basic use is the following:

./ccas query conformation

corrresponding:

query         query filename (without .pdb extension) to analyze. For NMR
              structures, introduce only one model.
conformation  PDB-ID_chain of a conformation belonging to query cluster.


The options supported by the script are the following:
-a --align        save aligned PDB file of the cluster conformers (e.g. "PDB-ID_aligned.pdb")
-p --print_table  print RMSD conformers hits table
-h --help         Show this help table

Outputs generated are three:
-> Conformation hits table.
-> Functoinal and Structural description of hits.
-> Aligned PDB files of the hits with respect to the query.

The outputs will be stored in its query respective output directory (e.g. "query_chain_output/").
If query file contains several chains, each chain will have its own directory.

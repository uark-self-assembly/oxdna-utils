This folder is made to contain different auxiliary files necessary to test the scripts. The following tests are:

generate-test-1:
This contains a simple test sequence for use with generate.py, generate-folded.py, or generate-RNA.py to make .top and .dat files out of the given sequence. Note that generate-folded.py and generate.py generate complement strands in reverse order from each other when the DOUBLE option is used.

hairpin18_MC_DOUBLE-trajectory.dat:
This is a .dat file from one of the tests in oxDNA/TEST. This file is placed hear for use with cullTraj.py. Do note that after running, cullTraj.py creates a file DO_NOT_DISCARD in the current directory that must be deleted before cullTraj.py will work again. This is to prevent accidentally culling the same file multiple times. 

##########################
#### Project 2 readme ####
##########################

Project:
This project is based on implemting global alignment with both a linear gap cost and affine gapcost
To use the files in this project you will need python.
Throughout the project multiple tests cases and questions were presented, and the majority of the 
results from running these tests can be found in the test_results folder

To answer the questions presented in the project 3 python scripts were made:
All the script require you to have the following files:
	- A fasta file with minimum 2 sequences in it that will be aligned
	- A txt file containing a phylip-like formatted matrix that is used as the substitution matrix
Run the scripts with the command-line and give the specifcied arguments. INPUt file with sequences will always be arg1.


##############################
### global_align_linear.py ### 
This script runs a global alignment with linear gapcost of 2 sequences given to in fasta format and return the optimal cost
and if unspecified also all possible alignments that lead to the optimal score through backtracking.
It both outputs the opt_score and alignments to the terminal, but also makes a fasta file containing
the results called aligned_sequences_linear in the current directory.

Run the script by:
python global_align_linear.py <SEQUENCE_FILE.FASTA> -g <GAPSCORE> -m <SUBSTITUTION_MATRIX_file> --hide-alignments
	
	note: 
	-m: the substition matrix file is a txt file containing a matrix in the phylip-like format
	--hide-alignments: If this argument is given the backtracking is skipped, increasing runtime

	Args:
	1. Fasta fille with sequences
	2. -g GAPSCORE
	3. -e GAP_EXTENSION
	4. -m SUBSTITUTION_MATRIX_file
	(5. --hide-alignments)

Example: 
python global_align_linear.py data/sequences_case1.fasta -g 5 -m substitution_matrix_phylip.txt



##############################
### global_align_affine.py ###
This script runs exactly like the linear version expect it takes an additional argument for the gap
extension cost.

Run the script by:
python global_align_affine.py <SEQUENCE_FILE.FASTA> -g <GAPSCORE> -e <GAP_EXTENSION> -m <SUBSTITUTION_MATRIX_file> --hide-alignments

	note: 
	-m: the substition matrix file is a txt file containing a matrix in the phylip-like format
	--hide-alignments: If this argument is given the backtracking is skipped, increasing runtime

	Args:
	1. Fasta fille with sequences
	2. -g GAPSCORE
	3. -e GAP_EXTENSION
	4. -m SUBSTITUTION_MATRIX_file
	5. -e GAP_EXTENSION
	(6. --hide-alignments)
	

example: 
python global_align_affine.py data/sequences_case1.fasta -g 5 -e 5 -m substitution_matrix_phylip.txt



#############################
### multiple_sequences.py ###
Imports the affine and linear function from the other scripts. Depending on given arg will do alignment
with either one or the other. 
Aligns all possible combinations of multiple sequences, but only returns the optimal cost in a table.
Only outputs the table to the terminal


Run the script by:
python global_align_linear.py <SEQUENCE_FILE.FASTA> -g <GAPSCORE> -e <GAP_EXTENSION> -m <SUBSTITUTION_MATRIX_file> -a
	
	note: 
	-e: not required if doing lienar cost
	-m: the substition matrix file is a txt file containing a matrix in the phylip-like format
	-a: if option is given will use affine gap cost, if not uses linear as default

	Args:
	1. Fasta fille with sequences
	2. -g GAPSCORE
	3. -e GAP_EXTENSION
	4. -m SUBSTITUTION_MATRIX_file
	5. -e GAP_EXTENSION
	(6. -a) 
	

example: 
# linear:
python multiple_sequences.py data/sequences_eval_all.fasta -g 5 -e 5 -m substitution_matrix_phylip.txt 
# affine:
python multiple_sequences.py data/sequences_eval_all.fasta -g 5 -e 5 -m substitution_matrix_phylip.txt -a
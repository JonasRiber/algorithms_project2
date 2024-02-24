#################################
#### Multiple seq alignments ####
#################################


# import functions from the other scripts 
from global_align_linear import global_linear
from global_align_affine import global_affine


import argparse

# parse args, added the -a to determine if to use affine or default (deafult = linear)
def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform global linear sequence alignment.")
    parser.add_argument("fasta_file", type=str, help="Path to input FASTA file containing sequences")
    parser.add_argument("-m", "--matrix", type=str, help="Path to scoring matrix file in Phylip-like format", required=True)
    parser.add_argument("-g", "--gap", type=int, help="Gap penalty", required=True)
    parser.add_argument("-e", "--extension", type=int, help="Gap extension penalty")
    parser.add_argument('-a', "--affine", action='store_true', help='use affine gap cost (if not specified use linear gap cost)')
    
    args = parser.parse_args()
    
    if args.matrix is None or args.gap is None:
        parser.error("-m <matrix_file> , -g <gapscore> is missing")

    return args


# Same read fasta function as before
def read_fasta_file(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:]
                sequence = ""
            else:
                sequence += line
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences

# same read phylip function as before
def read_phylip_like_matrix(matrix_file):
    matrix = {}
    with open(matrix_file, 'r') as file:
        num_chars = int(file.readline().strip())
        characters = ['A', 'C', 'G', 'T']
        for char in characters:
            matrix[char] = {}
        for line in file:
            line = line.strip().split()
            char = line[0]
            scores = list(map(int, line[1:]))
            for i, score in enumerate(scores):
                matrix[char][characters[i]] = score
    return matrix



def main():
    args = parse_arguments()
    fasta_file = args.fasta_file
    matrix_file = args.matrix
    gap_penalty = args.gap
    extension_penalty = args.extension
    gap_cost_type = args.affine

    # Load in substitution matrix - takes phylip makes nested dict
    scoring_matrix = read_phylip_like_matrix(matrix_file)

    # load sequences and define individual ones
    sequences = read_fasta_file(fasta_file)
    
    seq = [] # Set sequences in a list, seq1 corosponds to idx 0
    for s in sequences.values():
        seq.append(s)
        
    def compute_alignment_scores(sequences):
        alignment_scores = [[None for _ in range(len(sequences))] for _ in range(len(sequences))]
        
        if gap_cost_type == False:
            for i in range(len(sequences)):
                for j in range(len(sequences)):
                    global_align = global_linear(sequences[i], sequences[j], gap_penalty, scoring_matrix, hide_alignments = False) 
                    alignment_scores[i][j] = global_align[0]
        else:
            for i in range(len(sequences)):
                for j in range(len(sequences)):
                    global_align = global_affine(sequences[i], sequences[j], gap_penalty, extension_penalty, scoring_matrix, hide_alignments = False) 
                    if isinstance(global_align, int):
                        alignment_scores[i][j] = global_align
                    else: 
                        alignment_scores[i][j] = global_align[0]

        return alignment_scores

    # Compute the alignment scores for each pair of sequences
    alignment_scores_table = compute_alignment_scores(seq)
    
    # Print table with scores to terminal
    print("     ", end="")
    for label in sequences.keys():
        print("{:10}".format(label), end="")
    print()
    for i, row in enumerate(alignment_scores_table):
        print("{:5}".format(list(sequences.keys())[i]), end="")
        for score in row:
            print("{:10}".format(score), end="")
        print()


if __name__ == "__main__":
    main()



#python multiple_sequences.py <SEQUENCE_FILE.FASTE> -g <GAPSCORE> -e <GAP_EXTENSION> -m <SUBSTITUTION_MATIX> -a

# Example use
# use linear gapcost:    
# python multiple_sequences.py sequences_eval_all.fasta -g 5 -e 5 -m substitution_matrix_phylip.txt
    
# use affine gapcost:
# python multiple_sequences.py sequences_eval_all.fasta -g 5 -e 5 -m substitution_matrix_phylip.txt -a
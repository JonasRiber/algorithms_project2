
import argparse

# Define the args for the command-line call
def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform global linear sequence alignment.")
    parser.add_argument("fasta_file", type=str, help="Path to input FASTA file containing sequences")
    parser.add_argument("-m", "--matrix", type=str, help="Path to scoring matrix file in Phylip-like format", required=True)
    parser.add_argument("-g", "--gap", type=int, help="Gap penalty", required=True)
    parser.add_argument("-e", "--extension", type=int, help="Gap extension penalty", required=True)
    parser.add_argument('--hide-alignments', action='store_true', help='Hide aligned sequences')
    
    args = parser.parse_args()
    
    if args.matrix is None or args.gap is None or args.extension is None:
        parser.error("-m <matrix_file> , -g <gapscore> or -e <extension_penalty> is missing")

    return args


# read the fasta file containing the 2 sequences
# format looks like:
# >seq1
# acgtgtcaacgt 
# >seq2
# acgtcgtagcta 
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


# Read the .txt file containing the substitution matrix
# uses the phylip-like format
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


# Global alignment with linear gapcost
# takes 5 args, A and B are str representing sequences
# g is gap cost, sub_matrix the substitution_matrix
# and show_alignment defines if the alignment should be included

def global_affine(A, B, gapstart, gapextend, substitution_matrix, hide_alignments=True):
    A = A.upper()
    B = B.upper() 

    # Get cost out of substitution matrix
    def cost(x, y):
        return substitution_matrix[x][y]

    def g(k):
        return gapstart+gapextend*k 
    
    # Initialize matrix
    n = len(A)
    m = len(B)

    S = [[float('inf')] * (m + 1) for _ in range(n + 1)] # since minimization we use 'inf'
    D = [[float('inf')] * (m + 1) for _ in range(n + 1)]
    I = [[float('inf')] * (m + 1) for _ in range(n + 1)]
    
    S[0][0] = 0 # Nothing aligned

    # calcing edges
    for i in range(1, n + 1):
        D[i][0] = gapstart + gapextend*i
        S[i][0] = D[i][0]
    for j in range(1, m + 1):
        I[0][j] = gapstart + gapextend*j
        S[0][j] = I[0][j]
    
    
    # Loop through the matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            D[i][j] = min(S[i-1][j] + gapstart + gapextend, D[i-1][j] + gapextend) # when moving left (start new vs extension)
            I[i][j] = min(S[i][j-1] + gapstart + gapextend, I[i][j-1] + gapextend) # Moving down
            S[i][j] = min(S[i-1][j-1] + cost(A[i-1], B[j-1]), D[i][j], I[i][j]) # compare diagonal move vs del vs insert
            

    # print('\n'.join(' '.join(map(str, sub)) for sub in S)) # Print matrix

    opt_cost = S[-1][-1]
    if not hide_alignments:
        return S[-1][-1]
    
    
    # Backtracking - Needs to be recursive to get all the alignments 
    # The iterative version only gives one:
    def backtrack_recursive(A, B, S, i, j, alignmentA='', alignmentB='', alignments=None):
        if alignments is None: # initialize new lst
            alignments = []
        if i == 0 and j == 0: # base case 
            alignments.append([alignmentA, alignmentB])
        else: 
            if i > 0 and j > 0 and S[i][j] == S[i-1][j-1] + cost(A[i-1], B[j-1]): # diagonal move
                backtrack_recursive(A, B, S, i-1, j-1, A[i-1] + alignmentA, B[j-1] + alignmentB, alignments)
            k = 1 # if not diagonal move then gap, k need to determine gapchunk size
            while True: # each iteration represents finding a consecutive gap, increasing k+=1
                if k > i and k > j: # termination of while loop
                    break
                if i >= k and S[i][j] == S[i-k][j] + g(k): # gap in B
                    backtrack_recursive(A, B, S, i-k, j, A[i-k:i] + alignmentA, '-'*k + alignmentB, alignments)
                if j >= k and S[i][j] == S[i][j-k] + g(k): # gap in A
                    backtrack_recursive(A, B, S, i, j-k, '-'*k + alignmentA, B[j-k:j] + alignmentB, alignments)
                k += 1 
        return alignments
    
    alignments = backtrack_recursive(A, B, S, len(A), len(B))
    return opt_cost, alignments


# write the output to a fasta file 
def write_alignment_to_fasta(aligned_sequences, output_file):
    with open(output_file, 'w') as file:
        if isinstance(aligned_sequences, int):  # If aligned_sequences is an integer (optimal cost)
            file.write(">Optimal cost: {}\n".format(aligned_sequences))
        else:  # If aligned_sequences is a list of aligned sequences
            optimal_cost, sequence_pairs = aligned_sequences
            file.write(">Optimal cost: {}\n".format(optimal_cost))
            for i, sequence_pair in enumerate(sequence_pairs):
                file.write(">Alignment {}\n".format(i + 1))
                for j, sequence in enumerate(sequence_pair):
                    #file.write(">seq{}\n".format(j + 1))
                    file.write(sequence + "\n")


def main():
    args = parse_arguments()
    fasta_file = args.fasta_file
    matrix_file = args.matrix
    gap_penalty = args.gap
    extension_penalty = args.extension
    hide_alignments = args.hide_alignments

    sequences = read_fasta_file(fasta_file)
    scoring_matrix = read_phylip_like_matrix(matrix_file)

    sequence1 = sequences["seq1"]  # Assuming the sequence IDs are known
    sequence2 = sequences["seq2"]

    aligned_sequences = global_affine(sequence1, sequence2, gap_penalty, extension_penalty, scoring_matrix, hide_alignments = not hide_alignments)
    
    # Print the aligned sequences to terminal
    if isinstance(aligned_sequences, int):
        print("Optimal cost: {}".format(aligned_sequences))
    else:    
        for i in aligned_sequences:    
            if isinstance(i, list):
                for j in i:
                    print(j)
            else:
                print("Optimal cost: {}".format(i))
    
    # Use the write function to make putput file
    write_alignment_to_fasta(aligned_sequences, "affine_aligned_sequences.fasta")

if __name__ == "__main__":
    main()




# Command-line use cases:
# python global_align_linear.py <sequences_file> -m <substitution_matrix_file> -g <gap_score>     


# examples:
# python global_align_linear.py sequences.fasta -m substitution_matrix_phylip.txt -g 5 
# python global_align_linear.py sequences.fasta -m substitution_matrix_phylip.txt -g 5 --hide-alignments
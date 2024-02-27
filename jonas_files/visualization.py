import time
import matplotlib.pyplot as plt
from global_align_linear import global_linear
from global_align_affine import global_affine
from global_align_linear import read_phylip_like_matrix

scoring_matrix = read_phylip_like_matrix("./jonas_files/substitution_matrix_phylip.txt")

runtimes_linear = []
seq1 = "acgtgtcaacgt"
seq2 = "acgtcgtagcta"

for i in range(30):
    start_time = time.time()
    global_linear(seq1 * i, seq2 * i, 5, scoring_matrix, hide_alignments = False)
    end_time = time.time()
    runtime = end_time - start_time
    runtimes_linear.append(runtime)

runtimes_affine = []
seq1 = "acgtgtcaacgt"
seq2 = "acgtcgtagcta"

for i in range(30):
    start_time = time.time()
    global_affine(seq1 * i, seq2 * i, 5, 5, scoring_matrix, hide_alignments = False)
    end_time = time.time()
    runtime = end_time - start_time
    runtimes_affine.append(runtime)

# Plotting the runtimes
plt.rcParams['figure.figsize'] = [10, 6]
plt.plot(range(1, len(runtimes_linear) + 1), runtimes_linear, marker = 'o', label = "linear")
plt.plot(range(1, len(runtimes_affine) + 1), runtimes_affine, marker = 'o', label = "affine")
plt.xlabel('Length of Input Sequences')
plt.ylabel('Runtime [Sec]')
plt.title('Runtime of Function "global_linear" and "global_affine"')
plt.legend()
plt.grid(True)
plt.savefig("runtime_linear_affine")
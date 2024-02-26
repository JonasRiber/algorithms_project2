import time
import matplotlib.pyplot as plt
from jonas_files.global_align_linear import global_linear
from jonas_files.global_align_linear import read_phylip_like_matrix

scoring_matrix = read_phylip_like_matrix("./jonas_files/substitution_matrix_phylip.txt")

runtimes = []
seq1 = "acgtgtcaacgt"
seq2 = "acgtcgtagcta"

for i in range(30):
    start_time = time.time()
    global_linear(seq1 * i, seq2 * i, 5, scoring_matrix)
    end_time = time.time()
    runtime = end_time - start_time
    runtimes.append(runtime)

# Plotting the runtimes
plt.rcParams['figure.figsize'] = [10, 6]
plt.plot(range(1, len(runtimes) + 1), runtimes, marker = 'o')
plt.xlabel('Length of Input Sequences')
plt.ylabel('Runtime [Sec]')
plt.title('Runtime of Function "global_linear"')
plt.grid(True)
plt.show()
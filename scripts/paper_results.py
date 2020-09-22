#!/usr/bin/python3

"""
Input: A list of fasta/fastq to be sketched, a range (k_min, k_max) of k-mer lengths, epsilon, a working directory
Computations:
- for each dataset and k-mer value apply kmc
- get the L1 norm of the kmc databases
- get the skewness of the k-mer spectrum
- sketch the resulting kmc databases with fress sense
- run fress check to have (sum of errors, average error, max error)
- get the (theoretical) uncompressed size for each fress sketch
- get the compressed size of each fress sketch
Output:
- A big table in tsv format with the following columns:

dataset name | k-value | spectrum skew | threshold | L1 sum of deltas | average delta | max delta | uncompressed size | compressed size

- and relative plots:

k-value as x and (skew, size) as y
"""
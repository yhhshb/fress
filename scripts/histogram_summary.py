#!/usr/bin/python3

import sys
import math
 
nh = dict()
with open(sys.argv[1], "r") as histo:
    for line in histo:
        freq, kcnt = line.split('\t')
        freq = int(freq)
        kcnt = int(kcnt)
        if(kcnt != 0):
            nh[freq] = kcnt
             
histo = nh
L0 = len(histo)
L1 = sum(k*v for k,v in histo.items())
L0_kmers = sum(histo.values())
pd = [v/L0_kmers for k, v in histo.items()]
log_pd = [math.log(p, 2) for p in pd]
H0 = sum([-p*lp for p,lp in zip(pd, log_pd)])
sys.stderr.write("H0 = {}\n".format(H0))
sys.stderr.write("Total number of bits for an exact map = {}\n".format(L1*H0 + L1*math.log(H0 + 1, 2)))
sys.stderr.write("Number of different frequencies = {}\n".format(L0))
sys.stderr.write("Total number of elements = {}\n".format(L1))
sys.stderr.write("Total number of unique elements = {}\n".format(L0_kmers))
mcnt = min(histo.keys())
Mcnt = max(histo.keys())
sys.stderr.write("min count = {} with {} items having it\n".format(mcnt, histo[mcnt]))
sys.stderr.write("MAX count = {} with {} items having it\n".format(Mcnt, histo[Mcnt]))

"""
from math import factorial as fac

def binomial(x, y):
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom

def dumpMatrices(matrices: list):
    print("number of levels = ", lc)
    for m in matrices:
        for row in m: print(row)
        print("---------------------------------------------------")

def constructIndex(matrices: list):
    index = [dict() for _ in matrices]
    for k,mat in zip(range(len(matrices)), matrices):
        dim = len(mat)
        for i in range(dim):
            for j in range(dim):
                if (mat[i][j] != None):
                    index[k][mat[i][j]] = (i,j)
    return index

def getBinaryRepresentation(index: dict, counter: int):
    queue = [(counter, 0)]
    brepr = list()
    while(queue):
        cnt, level = queue.pop(0)
        if(level == len(index)):
            brepr.append(cnt)
        else:
            i, j = index[level][cnt]
            queue.append((i, level+1))
            queue.append((j, level+1))
    return brepr
 
def getBinaryWeight(brepr: list):
    w = 0
    for b in brepr:
        if(b): w += 1
    return w
    
------------------------------------------------------------------------------------------------

frequencies = list(histo.keys())
sorted_hist = sorted(histo.items(), key=lambda item: item[1], reverse=True)

sq = L0
lc = 0
matrices = list()
while(sq > 2):
    lc += 1
    sq = math.ceil(math.sqrt(sq))
    matrix = [[None for _ in range(sq)] for _ in range(sq)]
    nh = dict()
    counters = [k for k, v in sorted(histo.items(), key=lambda item: item[1], reverse=True)]
    index = 0
    slce = 0
    for slce in range(2*sq-1): #diagonals
        z = 0 if slce < sq else slce - sq + 1
        for j in range(z, slce-z+1):
            if index < len(counters):
                i = slce - j
                matrix[j][i] = counters[index]
                if i in nh: nh[i] += histo[counters[index]]
                else: nh[i] = histo[counters[index]]
                if j in nh: nh[j] += histo[counters[index]]
                else: nh[j] = histo[counters[index]]
            index += 1
    histo = nh
    matrices.append(matrix)

sys.stderr.write("{}\n".format(histo)) # {0: number of zeros, 1: number of ones}

index = constructIndex(matrices)
code_len = 2**len(index)
wcnts = [0 for _ in range(code_len)]

'''
for freq in frequencies:
    bit_representation = getBinaryRepresentation(index, freq)
    weight = getBinaryWeight(bit_representation)
    wcnts[weight] += 1
    sys.stdout.write("{:>{padding}} --> {} | w = {}\n".format(freq, bit_representation, weight, padding=len(str(frequencies[-1]))))
sys.stderr.write("{}\n".format(wcnts))
'''

npos1 = 0
nones = 0
start = 0
end = 0
num_codes = 0
while num_codes < L0:
    q = binomial(npos1, code_len)
    num_codes += q
    nones += npos1 * sum(v for _, v in sorted_hist[start:num_codes])
    start = num_codes
    npos1 += 1
    
sys.stderr.write("0: {}, 1: {}\n".format(code_len * L1 - nones, nones))
sys.stderr.write("Space for storing only the counters = {} B\n".format(code_len * L1 // 8))
"""





































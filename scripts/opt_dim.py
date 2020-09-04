#!/usr/bin/python3

import math

def coll_prob(c, b):
    return 1 - (math.e**(-c/b) if c < b else 0)

def pairwise_prob(f1, f2, c1, c2,r, b):
    return abs(f1 - f2) * (coll_prob(c1, b) * coll_prob(c2, b))**r

def row_prob(histo, idx, r, b):
    rs = 0
    for other, co in histo[idx:]:
        rs += pairwise_prob(histo[idx][0], other, histo[idx][1], co, r, b)
    return rs

def L1_error(histo, r, b):
    F = len(histo)
    error = 0
    for i in range(F):
        error += histo[i][1] * row_prob(histo, i, r, b)
    return error

def optimize(histo, e, r, b):
    L1 = sum([f * c for f, c in histogram])
    thr = e*L1
    constr = constb = False
    if b == 0 or b == None: 
        b = int(math.ceil(-histo[1][1] / math.log(0.5))) #TODO Compute B taking into account the absolute difference of the two heaviest elements (not necessarily 1 and 2)
    else: 
        constb = True
    if r == 0 or r == None: r = int(math.ceil(math.log(e) / math.log(0.5)))
    else: constr = True
    
    print("Starting (r, b) = ({}, {}) with threshold = {} for a L1 norm of {}".format(r, b, thr, L1))
    error = L1_error(histo, r, b)
    print("Starting error = " + str(error))
    if error < thr: return r, b
    elif constr and constb: raise RuntimeError("r = {} and b = {} do not allow to achieve the desired epsilon but only: {}".format(r, b, error/L1))
    
    rinc = 0
    binc = 0
    while(error > thr):
        if not constb: binc += 1
        if(not constr and binc >= b/r):
            binc = 0
            rinc += 1
        error = L1_error(histo, r + rinc, b + binc)
    return r, b

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--histo", help="histogram file", required=True)
    parser.add_argument("-e", "--epsilon", help="maximum error rate", type=float, required=True)
    parser.add_argument("-r", "--nrows", help="number of rows", type=int)
    parser.add_argument("-b", "--ncolumns", help="number of columns", type=int)
    
    args = parser.parse_args()
    histogram = list()
    with open(args.histo, "r") as hf:
        for line in hf:
            histogram.append(tuple(map(int, line.split('\t'))))
    histogram.sort(key=lambda tup: tup[1], reverse=True)
    optr, optb = tuple(map(math.ceil, optimize(histogram, args.epsilon, args.nrows, args.ncolumns)))
    print("Optimal (r, b) = ({}, {})".format(optr, optb))

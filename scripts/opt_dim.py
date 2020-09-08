#!/usr/bin/python3

import sys
import math
import functools

def L1_error(histo, r, b):
    #coll_prob = lambda c, b: 1 - (math.e**(-c/b) if c < b else 0)
    coll_probs = [1-math.e**(-c/b) if c < b else 1 for f, c in histo]
    no_coll_probs = [1-p for p in coll_probs]
    F = len(histo)
    
    #jmax_probs = [1 for _ in range(F)]
    #for j in range(F-2, -1, -1): jmax_probs[j] = jmax_probs[j+1] * no_coll_probs[j]
    #jmax_probs = [1 - cp*ncp for cp, ncp in zip(coll_probs, jmax_probs)]
    
    error = 0
    for i in range(F):
        rs = 0
        for j in range(i+1, F):
            diff = abs(histo[j][0] - histo[i][0])
            rs += (coll_probs[i] * coll_probs[j])**r * diff #* jmax_probs[j]**r
        error += histo[i][1] * rs
    return error

def optimize(histo, e, r, b):
    L1 = sum([f * c for f, c in histogram])
    thr = e*L1
    
    constr = constb = False
    if b == 0 or b == None:
        b = int(math.ceil(-histo[1][1] / math.log(0.5)))#TODO take into account the delta between the two heaviest elements
    else: 
        constb = True
    if r == 0 or r == None: r = int(math.ceil(math.log(e) / math.log(0.5)))
    else: constr = True
    
    error = L1_error(histo, r, b)
    print("Threshold = {} for a L1 norm of {}".format(thr, L1))
    print("(r, b) = ({}, {}) -> error = {}".format(r, b, round(error)))
    if constr and constb and error > thr: raise RuntimeError("r = {} and b = {} do not allow to achieve the desired epsilon but only: {}".format(r, b, error/L1))
    old_r = r
    while(error > thr):#increase r to the minimum value to achieve the desired threshold
        r += 1
        error = L1_error(histo, r, b)
    dim = r*b#total number of cells
    if constb:#if b was set by the user we are done
        return r, b
    elif constr:#reduce r to old_r
        r = old_r
    else:#reduce r and increase b as long as the error is below threshold
        while(error < thr):
            r -= 1
            b = math.ceil(dim/r)
            error = L1_error(histo, r, b)
            print("(r, b) = ({}, {}) -> error = {}".format(r, b, round(error)))
        r += 1
    b = math.ceil(dim/r)#set b to maintain constant memory
    #TODO optimize b 
    if b < len(histo):#FIXME this should be part of the optimization of b
        sys.stderr.write("Warning, b is smaller than the total number of labels, setting to that value\n")
        b = len(histo)
    return r, b
    
    """
    rinc = 0
    binc = 0
    while(error > thr):
        if not constb: binc += 1
        if(not constr and binc > b/r):
            binc = 0
            rinc += 1
        error = L1_error(histo, r + rinc, b + binc)
    r+=rinc
    b+=binc
    return r, b
    """
        

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
    histogram.sort(key=lambda tup: tup[1])#Changed from reverse = True to False because the min-counter item takes all
    optr, optb = optimize(histogram, args.epsilon, args.nrows, args.ncolumns)
    print("Optimal (r, b) = ({}, {})".format(optr, optb))

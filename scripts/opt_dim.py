#!/usr/bin/python3

import sys
import math
import functools

def L1_error(histo, r, b, explicit):
    coll_probs = [1-(1-1/b)**c for f, c in histo]
    sys.stderr.write(str(coll_probs) + '\n')
    F = len(histo)    
    error = 0
    for i in range(F):
        rs = 0
        for j in range(i+1, F):
            diff = abs(histo[j][0] - histo[i][0])
            rs += (coll_probs[j])**r * diff #* jmax_probs[j]**r
        error += histo[i][1] * rs
    return error

def optimize(histo, e, r, b, explicit):
    L1 = sum([f * c for f, c in histogram])
    thr = e*L1
    
    constr = constb = False
    if b == 0 or b == None:
        b = int(math.ceil(-histo[1][1] / math.log(0.5)))#TODO take into account the delta between the two heaviest elements
    else: 
        constb = True
    if r == 0 or r == None: r = int(math.ceil(math.log(e) / math.log(0.5)))
    else: constr = True
    
    error = L1_error(histo, r, b, explicit)
    print("Threshold = {} for a L1 norm of {}".format(thr, L1))
    print("(r, b) = ({}, {}) -> error = {}".format(r, b, round(error)))
    if constr and constb and error > thr: raise RuntimeError("r = {} and b = {} do not allow to achieve the desired epsilon but only: {}".format(r, b, error/L1))
    old_r = r
    while(error > thr):#increase r to the maximum value to achieve the desired threshold
        r += 1
        error = L1_error(histo, r, b, explicit)
    dim = r*b#total number of cells
    if constb:#if b was set by the user we are done
        return r, b
    elif constr:#reduce r to old_r
        r = old_r
    else:#reduce r and increase b as long as the error is below threshold
        while(error < thr):
            r -= 1
            b = math.ceil(dim/r)
            error = L1_error(histo, r, b, explicit)
            print("(r, b) = ({}, {}) -> error = {}".format(r, b, round(error)))
        r += 1
    b = math.ceil(dim/r)#set b to maintain constant memory
    #TODO optimize b 
    if b < len(histo):#FIXME this should be part of the optimization of b
        sys.stderr.write("Warning, b is smaller than the total number of labels, setting to that value\n")
        b = len(histo)
    return r, b

def str2bool(v):
    if isinstance(v, bool): return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'): return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'): return False
    else: raise argparse.ArgumentTypeError("Boolean value expected.")        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--histo", help="histogram file", required=True)
    parser.add_argument("-e", "--epsilon", help="maximum error rate", type=float, required=True)
    parser.add_argument("-r", "--nrows", help="number of rows", type=int)
    parser.add_argument("-b", "--ncolumns", help="number of columns", type=int)
    parser.add_argument("-x", "--explicit", help="The heaviest element is explicitly included in the sketch", type=str2bool, nargs='?', const=True, default=False)
    
    args = parser.parse_args()
    histogram = list()
    with open(args.histo, "r") as hf:
        for line in hf:
            histogram.append(tuple(map(int, line.split('\t'))))
    histogram.sort(key=lambda tup: tup[1], reverse=True)
    optr, optb = optimize(histogram, args.epsilon, args.nrows, args.ncolumns, args.explicit)
    print("Optimal (r, b) = ({}, {})".format(optr, optb))

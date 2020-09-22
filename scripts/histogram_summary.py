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
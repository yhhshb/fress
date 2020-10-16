import os
import sys
import math
import ntpath
import random
import numpy as np

import time

enc_nuc_table = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

"""
enc_nuc_table_extended deals with letters different from 
{A, C, G, T, N} by simply replacing them with one non-ambigous DNA
U -> T
R -> G | Y -> C | K -> G | M -> A | S -> C |W -> A
B -> A | D -> C | H -> G | V -> T
Note that B, D, H and V are mapped to the letters they
do not represent for sure. This breaks the definition
but it is only to have a deterministic choice for the
replacement. DO NOT USE THIS TABLE other than set comparisons.
"""
enc_nuc_table_extended = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4, #5 is the gap character '-' -> it should be filtered out
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 0, 1,  1, 4, 4, 2,  2, 4, 4, 2,  4, 0, 4, 4, 
    4, 4, 3, 2,  3, 3, 3, 0,  4, 2, 4, 4,  4, 4, 4, 4, 
    4, 0, 0, 1,  1, 4, 4, 2,  2, 4, 4, 2,  4, 0, 4, 4, 
    4, 4, 3, 2,  3, 3, 3, 4,  4, 2, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

dec_nuc_table = ['A', 'C', 'G', 'T', 'N']

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

error_code = lambda k: 4 ** k

def path_leaf(path):
    head, tail = ntpath.split(path)
    if tail: return head, tail
    else: return head, ntpath.basename(head)

def is_fastx(filename):
    if filename.endswith("fa"): return ("fasta", False)
    elif filename.endswith("fq"): return ("fastq", False)
    elif filename.endswith("fna"): return ("fasta", False)
    elif filename.endswith("fnq"): return ("fastq", False)
    elif filename.endswith("fasta"): return ("fasta", False)
    elif filename.endswith("fastq"): return ("fastq", False)
    elif filename.endswith("fa.gz"): return ("fasta", True)
    elif filename.endswith("fq.gz"): return ("fastq", True)
    elif filename.endswith("fna.gz"): return ("fasta", True)
    elif filename.endswith("fnq.gz"): return ("fasta", True)
    elif filename.endswith("fasta.gz"): return ("fasta", True)
    elif filename.endswith("fastq.gz"): return ("fasta", True)
    else: return False

def getFastx(path: str):
    fastx_files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and is_fastx(f)]
    fastx_files.sort(key = lambda v : v.upper())
    return fastx_files

def paths2names(paths: list):
    names = [str(os.path.basename(path)).split('.')[0] for path in paths]
    return names

def getTSV(path: str):
    tsv_files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and os.path.splitext(f)[1] == ".tsv"]
    return tsv_files

def make_qgram_dist(wj, l1, l2):
    """Transforms a weighted jaccard index to a Q-gram distance"""
    return (l1 + l2) * ((1.0 - wj)/(1.0 + wj))

def make_mash_dist(k, J):
    if (J < 0): raise ValueError("Negative Jaccard index")
    return -1.0/k * math.log(2*J / (1+J)) if J != 0.0 else 1.0

def encode(kmer, table) -> int:
    code = 0
    for ch in kmer:
        code *= 4
        base_code = table[ord(ch)]
        if (base_code > 3):
            return error_code(len(kmer)) #this is outside the max value
        code += base_code
    return code

def decode(enc: int) -> str:
    code, length = enc
    ret = ''
    for _ in range(length):
        index = code & 3
        code >>= 2
        ret = dec_nuc_table[index] + ret
    return ret

def find_first_good_kmer(seq: str, k: int, start: int, table):
    while(encode(seq[start:start+k], table) == error_code(k)):
        #sys.stderr.write("error kmer = " + kmer + '\n')
        start += 1
    return start

def kmercount(seq: str, k: int, table):
    kmc = dict()
    l = len(seq)
    if l < k: return
    for i in range(l - k + 1):
        kmer_for = seq[i:(i+k)]
        if 'N' in kmer_for: continue
        kmer_rev = kmer_for.translate(comp_tab)[::-1]
        if kmer_for < kmer_rev: kmer = kmer_for
        else: kmer = kmer_rev
        kmer = encode(kmer, table)
        if kmer in kmc: kmc[kmer] += 1
        else: kmc[kmer] = 1
    return kmc

'''
def kmercount(seq: str, k: int, table):
    mask = 4**k - 1
    kmc = dict()
    idxmax = len(seq) - k # The +1 (last kmer) is handled at the end
    index = 0
    interrupted = True
    while(index < idxmax):
        if(interrupted):
            interrupted = False
            index = find_first_good_kmer(seq, k, index, table)
            fkmer = encode(seq[index:index+k-1], table)
            rkmer = 
            base = table[ord(seq[index+k-1])]
        kmer = (kmer << 2 + base) & mask
        #check = encode(seq[index:index+k], table)
        #sys.stderr.write(str(kmer) + " == " + str(check) + '\n')
        if(kmer in kmc): kmc[kmer] += 1
        else: kmc[kmer] = 1
        base = table[ord(seq[index+k])]
        if (base >= 4): interrupted = True
        index += 1
    if(base < 4):
        kmer = (kmer * 4 + base) & mask
        if(kmer in kmc): kmc[kmer] += 1
        else: kmc[kmer] = 1
    return kmc
'''

def kmc2ffp(kmc: dict, k: int):
    profile = np.zeros(4 ** k, dtype=np.int32)
    for hval, count in kmc.items():
        profile[hval] = count
    return profile

def kmc2minHash(kmc: dict, s: int):
    from datasketch import MinHash
    mh = MinHash(num_perm=s)
    for hval in kmc.keys():
        mh.update(hval.to_bytes(3, "big"))
    return mh

def kmc2WmH(kmc: dict, s: int):
    from datasketch import MinHash
    mh = MinHash(num_perm=s)
    for hval, count in kmc.items():
        for i in range(count):
            mh.update(hval.to_bytes(3, "big") + i.to_bytes(3, "big"))
    return mh

def getICWSGenerator(k: int, s: int):
    from datasketch import WeightedMinHashGenerator
    return WeightedMinHashGenerator(4 ** k, sample_size=s)

def kmc2ICWS(generator, kmc: dict, k: int, s: int):
    profile = kmc2ffp(kmc, k)
    return generator.minhash(profile)

'''
def ffp(seq: str, k: int, s: int, table=enc_nuc_table):
    seq_set = seq.split('N')
    kmc = defaultdict(lambda: 0)
    for token in seq_set:
        tklen = len(token)
        for i in range(tklen - k + 1):
            kmc[token[i:i+k]] += 1
    profile = np.zeros(4 ** k, dtype=np.int32)
    msk = MinHash(num_perm=s)
    for kmer, count in kmc.items():
        msk.update(kmer.encode("utf8"))
        profile[encode(kmer, table)] += count
    return msk, profile
'''

def QWJ(kmc1: dict, kmc2: dict):
    union = set(kmc1.keys()).union(set(kmc2.keys()))
    intersection_size = 0
    minsum = 0
    maxsum = 0
    qgram = 0
    for elem in union:
        if((elem in kmc1) and (elem in kmc2)):
            intersection_size += 1
            minsum += min(kmc1[elem], kmc2[elem])
            maxsum += max(kmc1[elem], kmc2[elem])
            qgram += abs(kmc1[elem] - kmc2[elem])
        elif((elem in kmc1) and (elem not in kmc2)):
            maxsum += kmc1[elem]
            qgram += kmc1[elem]
        elif((elem not in kmc1) and (elem in kmc2)):
            maxsum += kmc2[elem]
            qgram += kmc2[elem]
        else:
            raise ValueError("The union contains something alien to both datasets -> The impossible happened")
    return qgram, minsum / maxsum, intersection_size, len(union)

'''
def QWJ(ffp1, ffp2):
    if ffp1.shape != ffp2.shape:
        raise RuntimeError("the two FFPs don't have the same shape")
    union_size = 0
    intersection_size = 0
    min_sum = 0
    max_sum = 0
    qgram = 0
    for i in range(ffp1.size):
        if (ffp1[i] != 0 and ffp2[i] != 0):
            union_size += 1
            intersection_size += 1
        elif (ffp1[i] != 0 or ffp2[i] != 0):
            union_size += 1
        min_sum += min(ffp1[i], ffp2[i])
        max_sum += max(ffp1[i], ffp2[i])
        qgram += abs(ffp1[i] - ffp2[i])
    return qgram, min_sum / max_sum, intersection_size, union_size
'''

def edit_from_cigar(CIGAR: str):
    matches = mismatches = insertions = deletions = opens = 0
    i = 0
    while (i != len(CIGAR)):
        j = i + 1
        while (CIGAR[i:j].isdigit()):
            j += 1
        j -= 1
        if (CIGAR[j] == '='):
            matches += int(CIGAR[i:j])
        elif (CIGAR[j] == 'X'):
            mismatches += int(CIGAR[i:j])
        elif (CIGAR[j] == 'I'):
            insertions += int(CIGAR[i:j])
            opens += 1
        elif (CIGAR[j] == 'D'):
            deletions += int(CIGAR[i:j])
            opens += 1
        else:
            raise RuntimeError("Unable to parse CIGAR string: " + CIGAR)
        i = j + 1
    return mismatches + insertions + deletions, opens + mismatches, mismatches

def replace_bad_nucleotides(seq: str) -> str:
    """Replace letters different from {A, C, G, T, N}

    Substitutions are:
    - -> 'N' | U -> T
    R -> G | Y -> C | K -> G | M -> A | S -> C |W -> A
    B -> A | D -> C | H -> G | V -> T

    Note that B, D, H and V are mapped to the letters they
    do not represent for sure. This breaks the definition
    but it is only to have a deterministic choice for the
    replacement. DO NOT USE THIS TABLE other than set comparisons.
    """
    seq = seq.upper()
    seq = seq.replace('-', 'N').replace('U', 'T')
    seq = seq.replace('R', 'G').replace('Y', 'C').replace('K', 'G').replace('M', 'A').replace('S', 'C').replace('W', 'A')
    seq = seq.replace('B', 'A').replace('D', 'C').replace('H', 'G').replace('V', 'T')
    return seq

class DenseMarkovChain:

    def __init__(self, k: int = 0):
        self.k = k
        self.mk = np.zeros((4 ** k) * 4 if k != 0 else 0, dtype=np.uint64)
        self.L1 = 0

    def __add__(self, other):
        if (self.k != other.k): raise RuntimeError("The two markov chains are incompatibles")
        result = DenseMarkovChain()
        result.k = self.k
        result.L1 = self.L1 + other.L1
        result.mk = self.mk + other.mk
        return result

    def add_observation(self, seq: str):
        seq_set = seq.split('N')
        for token in seq_set:
            for i in range(len(token) - self.k): # -2 because we want to look at the last base as a separate element
                kmer = token[i:i + self.k]
                self.mk[encode(kmer, enc_nuc_table_extended) + encode(token[i+self.k], enc_nuc_table_extended)] += 1 #Remember: all ambiguous bases to ACGT
                self.L1 += 1

    def simulate(self, length: int, start: str = None) -> str:
        if (start == None):
            start = random.randint(0, 4 ** self.k - 1)
        elif (len(start) != self.k or encode(start, enc_nuc_table) == 4 ** self.k):
            raise ValueError("Incompatible starting point")
        kmer = start[:]
        for _ in range(length):
            r = random.randint(0, self.L1)
            code = encode(kmer, enc_nuc_table)
            if (self.mk[code] > r):
                next = 'A'
            elif (self.mk[code] + self.mk[code + 1] > r):
                next = 'C'
            elif (self.mk[code] + self.mk[code + 1] + self.mk[code + 2] > r):
                next = 'G'
            else:
                next = 'T'
            start = start + next
            kmer = kmer[1:] + next
        return start

    def open4store(self, path: str):
        return open(path + ".bin", "wb")

    def open4load(self, path: str):
        return open(path, "rb")

    def store(self, fd):
        np.array(self.k, dtype=np.uint8).tofile(fd)
        self.mk.tofile(fd)

    def load(self, fd):
        self.k = int(np.fromfile(fd, dtype=np.uint8, count = 1))
        dim = (4 ** self.k) * 4 if self.k != 0 else 0
        self.mk = np.fromfile(fd, dtype=np.uint64, count=dim)
        self.L1 = np.sum(self.mk)

class SparseMarkovChain:

    def __init__(self, k: int = 0):
        self.k = k
        self.mk = dict()
        self.L1 = 0

    def __add__(self, other):
        if (self.k != other.k): raise RuntimeError("The two markov chains are incompatibles")
        result = SparseMarkovChain()
        result.k = self.k
        result.L1 = self.L1 + other.L1
        result = self.mk.copy()
        result.mk.setdefault(None)
        for key, value in other.mk:
            if (key in result.mk):
                result.mk[key] += value
            else:
                result.mk[key] = value
        return result

    def add_observation(self, seq: str):
        self.mk.setdefault(1)
        seq_set = seq.split('N')
        for token in seq_set:
            for i in range(len(token) - self.k + 2): # -2 because we want to look at the last base as a separate element
                key = (token[i:i + self.k], token[i+self.k])
                self.mk[key] += 1 #Remember: all ambiguous bases to ACGT
                self.L1 += 1
        self.mk.setdefault(None)

    def simulate(self, length: int, start: str = None):
        if (start == None):
            start = random.randint(0, 4 ** self.k - 1)
        elif (len(start) != self.k or encode(start, enc_nuc_table) == 4 ** self.k):
            raise ValueError("Incompatible starting point")
        kmer = start[:]
        for _ in range(length):
            r = random.randint(0, self.L1)
            p1 = self.mk[(kmer, 'A')] if (kmer, 'A') in self.mk else 0
            p2 = self.mk[(kmer, 'C')] if (kmer, 'C') in self.mk else 0
            p3 = self.mk[(kmer, 'G')] if (kmer, 'G') in self.mk else 0
            if (p1 > r):
                next = 'A'
            elif (p1 + p2 > r):
                next = 'C'
            elif (p1 + p2 + p3 > r):
                next = 'G'
            else:
                next = 'T'
            start = start + next
            kmer = kmer[1:] + next
        return start

    def open4store(self, path: str):
        return open(path + ".txt", "w")

    def open4load(self, path: str):
        return open(path, "r")

    def store(self, fd):
        fd.write(str(self.k) + '\n')
        for (kmer, b), value in self.mk:
            fd.write(kmer + ' ' + b + ' ' + str(value) + '\n')
        fd.flush()

    def load(self, fd):
        self.k = int(fd.readline())
        for line in fd.readline():
            kmer, b, svalue = line.split()
            value = int(svalue)
            self.mk[(kmer, b)] = value
            self.L1 += value

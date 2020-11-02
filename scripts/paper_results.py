#!/usr/bin/python3

"""Run experiments for the SMS paper

This script requires the kmc.py module of the wgram repository (https://github.com/yhhshb/wgram)
Remember to put the install kmc in the wgram folder by running download_tools.sh and install_tools.sh
Remember to add the wgram folder to the sys.path of this script.
"""

import os
import sys
import zipfile
import math
import subprocess
import pandas
from scipy.stats import skew

import kmc

fressdir = os.path.dirname(os.path.normpath(os.path.dirname(os.path.abspath(__file__))))
fress = os.path.join(fressdir, "fress")

def compress(path: str, files: list, compressed_path: str):
    # Select the compression mode ZIP_DEFLATED for compression
    # or zipfile.ZIP_STORED to just store the file
    compression = zipfile.ZIP_DEFLATED

    # create the zip file first parameter path/name, second mode
    zf = zipfile.ZipFile(compressed_path, mode="w")
    try:
        for fname in files:
            # Add file to the zip file
            # first parameter file to zip, second filename in zip
            zf.write(os.path.join(path, fname), fname, compress_type=compression)
    except FileNotFoundError:
        sys.stderr.write("An error occurred during compression\n")
    finally:
        zf.close()

def run_fress_histogram(kmc_name: str, histo_name: str):
    out = subprocess.run([fress, "histogram", "-i", kmc_name, "-o", histo_name])
    if(out.returncode != 0): raise Exception("Error while computing the histogram for {}".format(kmc_name))

def run_fress_sense(kmc_name: str, sketch_name: str, epsilon: float, r: int = None, b: int = None):
    out = subprocess.run([fress, "sense", "-i", kmc_name, "-o", sketch_name, "-e", str(epsilon)] + (["-r", str(r)] if r else []) + (["-b", str(b)] if r else []), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_check(kmc_name: str, sketch_name: str, merged: int = 0, new_freq: float = 0):
    out = subprocess.run([fress, "check", "-i", kmc_name, "-d", sketch_name] + (["-g", str(merged)] if merged else []) + (["-f", str(new_freq)] if new_freq else []), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while checking the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_info(sketch_name: str):
    out = subprocess.run([fress, "info", "-d", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while getting info about {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_mms(kmc_name: str, sketch_name: str, epsilon: float, ignored: int = None):
    out = subprocess.run([fress, "mms", "-i", kmc_name, "-o", sketch_name, "-e", str(epsilon)] + (["-g", str(ignored)] if ignored else []), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the MM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_cms(kmc_name: str, sketch_name: str, epsilon: float, ignored: int = None):
    out = subprocess.run([fress, "cms", "-i", kmc_name, "-o", sketch_name, "-e", str(epsilon)] + (["-g", str(ignored)] if ignored else []), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the CM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_cmschk(kmc_name: str, sketch_name: str, ignored: int = None):
    out = subprocess.run([fress, "cmschk", "-i", kmc_name, "-d", sketch_name] + (["-g", str(ignored)] if ignored else []), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_bbhash(kmc_name: str, sketch_name: str):
    out = subprocess.run([fress, "bbhash", "-i", kmc_name, "-o", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the BBHash mphf {}".format(sketch_name))
    return out.stdout.decode("utf-8").split('\n')[-1].split()

def run_sms_for(fastx: str, k: int, epsilon: float, args):
    """Build and check one single SM sketch

    Input: 
    - one fasta/fastq to be sketched
    - the k-mer length
    - the approximation factor epsilon
    - a working directory
    - the output directory for kmc databases
    - the output directory for the sketch
    - a temporary directory
    - maximum allowed memory

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
    dataset name | epsilon | k | number of collisions | threshold | L1 sum of deltas | average delta | max delta | uncompressed size | compressed size
    
    ATTENTION: average delta is not (L1 sum of deltas / number of collisions) but an average computed over using collisions computed 
    as intersections of size different than one (instead of the wrong frequency)
    """
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)
    pre_file = kmcdb + ".kmc_pre"
    suf_file = kmcdb + ".kmc_suf"
    if(not os.path.exists(pre_file) or not os.path.exists(suf_file)): kmc.count(k, fastx, args.c, args.w, args.m, True)

    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    histo_name = sketch_name + ".shist.txt"
    cmb_name = sketch_name + ".cmb.txt"
    bin_name = sketch_name + ".bin"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(args.f, sketch_name)
    histo_path = os.path.join(args.f, histo_name)
    cmb_path = os.path.join(args.f, cmb_name)
    #bin_path = os.path.join(args.f, bin_name)
    arch_path = os.path.join(args.f, arch_name)

    L1, dim = run_fress_sense(kmcdb, sketch_path, epsilon)
    L1 = int(L1)
    dim = int(dim)
    ncolls, ntrue_colls, sod, avgd, maxd = run_fress_check(kmcdb, sketch_path)
    ncolls = int(ncolls)
    ntrue_colls = int(ntrue_colls)
    avgd = float(avgd)
    maxd = int(maxd)

    #histo = pandas.read_csv(histo_path, sep='\t', header=None)
    #skewness = skew(histo.to_numpy()[:,1])
    ncombinations = 0
    with open(cmb_path, "r") as hc:
        for _ in hc: ncombinations += 1
    theoretical_udim = round(dim * math.ceil(math.log(ncombinations, 2)) / 8) + os.stat(histo_path).st_size + os.stat(cmb_path).st_size
    compress(args.f, [histo_name, cmb_name, bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, ntrue_colls, round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim)

def run_sms_check(fastx: str, k: int, epsilon: float, args):
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)

    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    sketch_path = os.path.join(args.f, sketch_name)

    ncolls, ntrue_colls, sod, avgd, maxd = run_fress_check(kmcdb, sketch_path)
    ncolls = int(ncolls)
    ntrue_colls = int(ntrue_colls)
    avgd = float(avgd)
    maxd = int(maxd)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, ncolls, ntrue_colls, sod, avgd, maxd)

def run_sms_info(fastx: str, k: int, epsilon: float, args):
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, _ = kmc.getKMCPaths(k, fastx, args.c)
    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    sketch_path = os.path.join(args.f, sketch_name)

    R, B, RB = run_fress_info(sketch_path)
    R = int(R)
    B = int(B)
    RB = int(RB)
    return "{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, R, B, RB)

def run_mms_for(fastx: str, k: int, epsilon: float, args):
    """Build and check one single MM sketch

    See run_sms_for documetantion.
    """
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)
    pre_file = kmcdb + ".kmc_pre"
    suf_file = kmcdb + ".kmc_suf"
    if(not os.path.exists(pre_file) or not os.path.exists(suf_file)): kmc.count(k, fastx, args.c, args.w, args.m, True)

    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    bin_name = sketch_name + ".cms"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(args.f, sketch_name)
    arch_path = os.path.join(args.f, arch_name)

    L1, dim, max_val = run_fress_cms(kmcdb, sketch_path, epsilon, args.g)
    L1 = int(L1)
    dim = int(dim)
    max_val = int(max_val)
    ncolls, ntrue_colls, sod, avgd, maxd = run_fress_cmschk(kmcdb, sketch_path, args.g)
    ncolls = int(ncolls)
    ntrue_colls = int(ntrue_colls)
    avgd = float(avgd)
    maxd = int(maxd)

    theoretical_udim = round(dim * math.ceil(math.log(max_val, 2)) / 8)
    compress(args.f, [bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, ntrue_colls, round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim)

def run_cms_for(fastx: str, k: int, epsilon: float, args):
    """Build and check one single CM sketch

    See run_sms_for documetantion.
    """
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)
    pre_file = kmcdb + ".kmc_pre"
    suf_file = kmcdb + ".kmc_suf"
    if(not os.path.exists(pre_file) or not os.path.exists(suf_file)): kmc.count(k, fastx, args.c, args.w, args.m, True)

    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    bin_name = sketch_name + ".cms"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(args.f, sketch_name)
    arch_path = os.path.join(args.f, arch_name)

    L1, dim, max_val = run_fress_cms(kmcdb, sketch_path, epsilon, args.g)
    L1 = int(L1)
    dim = int(dim)
    max_val = int(max_val)
    ncolls, ntrue_colls, sod, avgd, maxd = run_fress_cmschk(kmcdb, sketch_path, args.g)
    ncolls = int(ncolls)
    ntrue_colls = int(ntrue_colls)
    avgd = float(avgd)
    maxd = int(maxd)

    theoretical_udim = round(dim * math.ceil(math.log(max_val, 2)) / 8)
    compress(args.f, [bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, ntrue_colls, round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim)

def run_bbhash_for(fastx: str, k: int, _, args):
    """Build and check BBHash MPHF with = 1

    See run_sms_for documentation about input parameters.
    Output:
    - A big table in tsv format with the following columns:
    dataset name | k-value | mphf uncompressed size | total uncompressed size | total compressed size
    """
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)
    pre_file = kmcdb + ".kmc_pre"
    suf_file = kmcdb + ".kmc_suf"
    if(not os.path.exists(pre_file) or not os.path.exists(suf_file)): kmc.count(k, fastx, args.c, args.w, args.m, True)

    sketch_name = "{}k{}".format(filename, k)
    mphf_name = sketch_name + ".bbh"
    payload_name = sketch_name + ".pld"
    arch_name = sketch_name + "_BBH.gz"
    sketch_path = os.path.join(args.f, sketch_name)
    mphf_path = os.path.join(args.f, mphf_name)
    arch_path = os.path.join(args.f, arch_name)

    max_val, L0 = run_fress_bbhash(kmcdb, sketch_path)
    max_val = int(max_val)
    L0 = int(L0)

    mphf_size = os.stat(mphf_path).st_size
    theoretical_udim = round(L0 * math.ceil(math.log(max_val, 2)) / 8) + mphf_size
    compress(args.f, [mphf_name, payload_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}".format(filename, k, mphf_size, theoretical_udim, cdim)

def run_combination(args, command):
    """Run fress for multiple parameters

    The input file must contain one line for each dataset in the following format:
    <dataset name> [k1, ..., kn] [epsilon1, ..., epsilon<m>]
    one line for each combination of parameter will be generated.
    """
    #description_file = args.file
    #output_file = args.o
    #kmc_outdir = args.c
    #fress_outdir = args.f
    #tmpdir = args.w
    #max_mem = args.m
    with open(args.o, "w") as oh:
        with open(args.file, "r") as dh:
            for line in dh:
                blocks = [block.replace(' ', '').strip() for block in line.split("[")]
                if(len(blocks) != 3): raise ValueError("Found a non-properly formatted line: {}".format(line))
                if (blocks[1][-1] != ']'): raise ValueError("Found a non-properly formatted line: {}".format(line))
                else: blocks[1] = blocks[1].replace(']', '')
                if (blocks[2][-1] != ']'): raise ValueError("Found a non-properly formatted line: {}".format(line))
                else: blocks[2] = blocks[2].replace(']', '')
                dataset = blocks[0].strip()
                ks = [int(block.strip()) for block in blocks[1].split(',')]
                epsilons = [float(block.strip()) for block in blocks[2].split(',')] if command != run_bbhash_for else [None]
                for e in epsilons:
                    for k in ks:
                        sys.stderr.write("Now fressing {} with k = {} and epsilon = {}\n".format(dataset, k, e))
                        oh.write(command(dataset, k, e, args) + "\n")
                        oh.flush()

def read_histo(histo_path: str):
    histogram = list()
    with open(histo_path, "r") as hf:
        for line in hf:
            histogram.append(tuple(map(int, line.split('\t'))))
    return histogram

def L1_error(histo: list, r: int, b: int, coll_probs: list):
    F = len(histo)
    error = 0
    for i in range(F):
        rs = 0
        for j in range(i+1, F):
            diff = abs(histo[j][0] - histo[i][0])
            rs += (coll_probs[j])**r * diff
        error += histo[i][1] * rs
    return error

def optimize_by_heuristic(histo: list, e: float, r: int, b: int):
    L1 = sum([f * c for f, c in histo])
    thr = e*L1
    constr = constb = False
    if b == 0 or b == None: b = int(math.ceil(-histo[1][1] / math.log(0.5)))
    else: constb = True
    if r == 0 or r == None: r = int(math.ceil(math.log(e) / math.log(0.5)))
    else: constr = True

    coll_probs = [1-(1-1/b)**c for f, c in histo]

    error = L1_error(histo, r, b, coll_probs)
    sys.stderr.write("Threshold = {} for a L1 norm of {}\n".format(thr, L1))
    sys.stderr.write("(r, b) = ({}, {}) -> error = {}\n".format(r, b, round(error)))
    if constr and constb and error > thr: raise RuntimeError("r = {} and b = {} do not allow to achieve the desired epsilon but only: {}".format(r, b, error/L1))
    old_r = r
    while(error > thr):#increase r to the maximum value to achieve the desired threshold
        r += 1
        error = L1_error(histo, r, b, coll_probs)
    dim = r*b#total number of cells
    if constb:#if b was set by the user we are done
        return r, b
    elif constr:#reduce r to old_r
        r = old_r
    else:#reduce r and increase b as long as the error is below threshold
        while(error < thr):
            r -= 1
            b = math.ceil(dim/r)
            error = L1_error(histo, r, b, coll_probs)
            print("(r, b) = ({}, {}) -> error = {}".format(r, b, round(error)))
        r += 1
    b = math.ceil(dim/r)#set b to maintain constant memory
    return r, b

def optimize_by_theorem(histo: list, e: float, r: int, b: int):
    rb = 2**63
    oldrb = rb + 1
    r = 0
    b = 0
    while(oldrb > rb):
        oldrb = rb
        r += 1
        b = histo[1][1] * math.exp(math.log(1/e)/r)
        rb = r * b
    r -= 1
    return r, math.ceil(histo[1][1] * math.exp(math.log(1/e)/r))

def opt_dim_main(histo_name: str, epsilon: float, r: int, b: int, to_collapse: int):
    histo = read_histo(histo_name)
    histo.sort(key=lambda tup: tup[1], reverse=True)
    new_frequency = None
    unavoidable_error = 0
    if (to_collapse > 1):
        merged = histo[:to_collapse]
        new_count = sum([c for _, c in merged])
        new_frequency = round(sum([f * c for f, c in merged]) / new_count, 6)
        unavoidable_error = math.ceil(sum([abs(new_frequency - f) * c for f, c in merged]))
        histo = [(new_frequency, new_frequency)] + histo[to_collapse:]
    hoptr, hoptb = optimize_by_heuristic(histo, epsilon, r, b)
    toptr, toptb = optimize_by_theorem(histo, epsilon, r, b)
    coll_probs = [1-(1-1/toptb)**c for f, c in histo]
    sys.stderr.write("New frequency = {}. Unavoidable error = {}\n".format(new_frequency, unavoidable_error))
    sys.stderr.write("Optimal by my heuristic (r, b) = ({}, {})\n".format(hoptr, hoptb))
    sys.stderr.write("Optimal by theorem (r, b) = ({}, {}) with error = {}\n".format(toptr, toptb, round(L1_error(histo, toptr, toptb, coll_probs))))
    return hoptr, hoptb, toptr, toptb, new_frequency, unavoidable_error

def run_merged_sms_for(fastx: str, k: int, epsilon: float, args):
    """Build and check one single SM sketch with column merging

    Input: 
    - one fasta/fastq to be sketched
    - the k-mer length
    - the approximation factor epsilon
    - a working directory
    - the output directory for kmc databases
    - the output directory for the sketch
    - a temporary directory
    - maximum allowed memory
    - number of columns to merge

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
    dataset name | epsilon | k-value | spectrum skew | threshold | L1 sum of deltas | average delta | max delta | uncompressed size | compressed size
    """
    import time
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, args.c)
    pre_file = kmcdb + ".kmc_pre"
    suf_file = kmcdb + ".kmc_suf"
    if(not os.path.exists(pre_file) or not os.path.exists(suf_file)): kmc.count(k, fastx, args.c, args.w, args.m, True)

    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    histo_name = sketch_name + ".shist.txt"
    cmb_name = sketch_name + ".cmb.txt"
    bin_name = sketch_name + ".bin"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(args.f, sketch_name)
    histo_path = os.path.join(args.f, histo_name)
    cmb_path = os.path.join(args.f, cmb_name)
    #bin_path = os.path.join(args.f, bin_name)
    arch_path = os.path.join(args.f, arch_name)

    tmp_name = str(time.time()) + ".hist.txt"
    run_fress_histogram(kmcdb, tmp_name)
    r, b, _, _, freq, unerr = opt_dim_main(tmp_name, epsilon, None, None, args.g)
    os.remove(tmp_name)

    L1, dim = run_fress_sense(kmcdb, sketch_path, epsilon, r, b)
    L1 = int(L1)
    dim = int(dim)
    ncolls, ntrue_colls, sod, avgd, maxd = run_fress_check(kmcdb, sketch_path, args.g, freq)
    ncolls = int(ncolls)
    ntrue_colls = int(ntrue_colls)
    sod = float(sod)
    avgd = float(avgd)
    maxd = float(maxd)
    tavg = sod/ntrue_colls

    histo = pandas.read_csv(histo_path, sep='\t', header=None)
    skewness = skew(histo.to_numpy()[:,1])
    ncombinations = 0
    with open(cmb_path, "r") as hc:
        for _ in hc: ncombinations += 1
    theoretical_udim = round(dim * math.ceil(math.log(ncombinations, 2)) / 8) + os.stat(histo_path).st_size + os.stat(cmb_path).st_size
    compress(args.f, [histo_name, cmb_name, bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, r, b, skewness, ncolls, ntrue_colls, round(L1 * epsilon), sod, avgd, tavg, maxd, theoretical_udim, cdim, args.g, freq, unerr)

def histogram_info(histo_path: str):
    histo = read_histo(histo_path)
    sys.stderr.write("Number of different frequencies = {}\n".format(len(histo)))
    sys.stderr.write("min frequency = {}\n".format(min([f for f, _ in histo])))
    sys.stderr.write("MAX frequency = {}\n".format(max([f for f, _ in histo])))
    sys.stderr.write("L0 norm = {}\n".format(sum([c for _, c in histo])))
    sys.stderr.write("L1 norm = {}\n".format(sum([f * c for f, c in histo])))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_sms = subparsers.add_parser("sms", help="Run set-min sketch pipeline for fixed values. This generates only one line of the final result table")
    parser_sms.add_argument("file", help="fastx file to process", type=str)
    parser_sms.add_argument("-k", help="k-mer value", type=int)
    parser_sms.add_argument("-e", help="epsilon", type=float)
    parser_sms.add_argument("-c", help="kmc folder")
    parser_sms.add_argument("-f", help="fress output folder")
    parser_sms.add_argument("-w", help="tmp dir")
    parser_sms.add_argument("-m", help="max memory for kmc", type=int)

    parser_smsm = subparsers.add_parser("smsm", help="Run set-min sketch pipeline for multiple values.")
    parser_smsm.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_smsm.add_argument("-o", help="output file (tsv format)", required=True)
    parser_smsm.add_argument("-c", help="kmc folder", required=True)
    parser_smsm.add_argument("-f", help="fress output folder", required=True)
    parser_smsm.add_argument("-w", help="tmp dir", required=True)
    parser_smsm.add_argument("-m", help="max memory for kmc", type=int)

    parser_cmsm = subparsers.add_parser("cmsm", help="Run count-min sketch pipeline for multiple values.")
    parser_cmsm.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_cmsm.add_argument("-o", help="output file (tsv format)", required=True)
    parser_cmsm.add_argument("-c", help="kmc folder", required=True)
    parser_cmsm.add_argument("-f", help="fress output folder", required=True)
    parser_cmsm.add_argument("-w", help="tmp dir", required=True)
    parser_cmsm.add_argument("-m", help="max memory for kmc", type=int)
    parser_cmsm.add_argument("-g", help="ignore this element during insertion", type=int)

    parser_bbhm = subparsers.add_parser("bbhm", help="Run BBHash pipeline for multiple values.")
    parser_bbhm.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_bbhm.add_argument("-o", help="output file (tsv format)", required=True)
    parser_bbhm.add_argument("-c", help="kmc folder", required=True)
    parser_bbhm.add_argument("-f", help="fress output folder", required=True)
    parser_bbhm.add_argument("-w", help="tmp dir", required=True)
    parser_bbhm.add_argument("-m", help="max memory for kmc", type=int)

    parser_checkonly = subparsers.add_parser("check", help="Run check for set-min sketches produced by a pipeline")
    parser_checkonly.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_checkonly.add_argument("-o", help="output file (tsv format)", required=True)
    parser_checkonly.add_argument("-c", help="kmc folder", required=True)
    parser_checkonly.add_argument("-f", help="fress output folder", required=True)

    parser_info = subparsers.add_parser("info", help="Run info for set-min sketches produced by a pipeline")
    parser_info.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_info.add_argument("-o", help="output file (tsv format)", required=True)
    parser_info.add_argument("-c", help="kmc folder", required=True)
    parser_info.add_argument("-f", help="fress output folder", required=True)

    parser_optdim = subparsers.add_parser("optdim", help="Dimension the sketch for given a histogram")
    parser_optdim.add_argument("histo", help="Histogram")
    parser_optdim.add_argument("-e", "--epsilon", help="maximum error rate", type=float, required=True)
    parser_optdim.add_argument("-r", "--nrows", help="number of rows", type=int)
    parser_optdim.add_argument("-b", "--ncolumns", help="number of columns", type=int)
    parser_optdim.add_argument("-m", "--merge", help="Merge the first columns of the histogram to augment the skewness", type=int, default=1)

    parser_mergecols = subparsers.add_parser("mergecols", help="Run set-min sketch pipeline for merged columns of the histogram.")
    parser_mergecols.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_mergecols.add_argument("-o", help="output file (tsv format)", required=True)
    parser_mergecols.add_argument("-c", help="kmc folder", required=True)
    parser_mergecols.add_argument("-f", help="fress output folder", required=True)
    parser_mergecols.add_argument("-w", help="tmp dir", required=True)
    parser_mergecols.add_argument("-m", help="max memory for kmc", type=int)
    parser_mergecols.add_argument("-g", help="number of columns to merge", type=int)

    parser_histinfo = subparsers.add_parser("histinfo", help="Get summary information about histogram")
    parser_histinfo.add_argument("histo", help="Input histogram file")

    args = parser.parse_args(sys.argv)
    if (args.command == "sms"): print(run_sms_for(args.file, args.k, args.e, args))
    elif (args.command == "smsm"): run_combination(args, run_sms_for)
    elif (args.command == "cmsm"): run_combination(args, run_cms_for)
    elif (args.command == "bbhm"): run_combination(args, run_bbhash_for)
    elif (args.command == "check"): run_combination(args, run_sms_check)
    elif (args.command == "info"): run_combination(args, run_sms_info)
    elif (args.command == "optdim"): opt_dim_main(args.histo, args.epsilon, args.nrows, args.ncolumns, args.merge)
    elif (args.command == "mergecols"): run_combination(args, run_merged_sms_for)
    elif (args.command == "histinfo"): histogram_info(args.histo)
    else: parser.print_help(sys.stderr)
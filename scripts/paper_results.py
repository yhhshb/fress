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

----------------------------------------------------------------

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

sys.path.append("/media/shibuya/workspace/wgram") #directory to wgram repository
sys.path.append("/home/igm/workspace/wgram")
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

def run_fress_sense(kmc_name: str, sketch_name: str, epsilon: float):
    out = subprocess.run([fress, "sense", "-i", kmc_name, "-o", sketch_name, "-e", str(epsilon)], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_check(kmc_name: str, sketch_name):
    out = subprocess.run([fress, "check", "-i", kmc_name, "-d", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while checking the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_cms(kmc_name: str, sketch_name: str, epsilon: float):
    out = subprocess.run([fress, "cms", "-i", kmc_name, "-o", sketch_name, "-e", str(epsilon)], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the CM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_cmschk(kmc_name: str, sketch_name):
    out = subprocess.run([fress, "cmschk", "-i", kmc_name, "-d", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the SM sketch {}".format(sketch_name))
    return out.stdout.decode("utf-8").split()

def run_fress_bbhash(kmc_name: str, sketch_name: str):
    out = subprocess.run([fress, "bbhash", "-i", kmc_name, "-o", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the BBHash mphf {}".format(sketch_name))
    return out.stdout.decode("utf-8").split('\n')[-1].split()

def run_sms_for(fastx: str, k: int, epsilon: float, kmc_outdir:str, fress_outdir: str, tmpdir: str, max_mem: int):
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    kmc.count(k, fastx, kmc_outdir, tmpdir, max_mem, True)

    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, kmc_outdir)
    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    histo_name = sketch_name + ".shist.txt"
    cmb_name = sketch_name + ".cmb.txt"
    bin_name = sketch_name + ".bin"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(fress_outdir, sketch_name)
    histo_path = os.path.join(fress_outdir, histo_name)
    cmb_path = os.path.join(fress_outdir, cmb_name)
    #bin_path = os.path.join(fress_outdir, bin_name)
    arch_path = os.path.join(fress_outdir, arch_name)

    L1, dim = run_fress_sense(kmcdb, sketch_path, epsilon)
    L1 = int(L1)
    dim = int(dim)
    ncolls, sod, avgd, maxd = run_fress_check(kmcdb, sketch_path)
    ncolls = int(ncolls)
    avgd = float(avgd)
    maxd = int(maxd)

    histo = pandas.read_csv(histo_path, sep='\t', header=None)
    skewness = skew(histo.to_numpy()[:,1])
    ncombinations = 0
    with open(cmb_path, "r") as hc:
        for _ in hc: ncombinations += 1
    theoretical_udim = round(dim * math.ceil(math.log(ncombinations, 2)) / 8) + os.stat(histo_path).st_size + os.stat(cmb_path).st_size
    compress(fress_outdir, [histo_name, cmb_name, bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{:.2f}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, skewness, round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim, "NA")

def run_cms_for(fastx: str, k: int, epsilon: float, kmc_outdir:str, fress_outdir: str, tmpdir: str, max_mem: int):
    if (epsilon < 0 or epsilon > 1): raise ValueError("epsilon must be a number between 0 and 1")
    kmc.count(k, fastx, kmc_outdir, tmpdir, max_mem, True)

    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, kmc_outdir)
    sketch_name = "{}k{}e{}".format(filename, k, str(epsilon).split('.')[1])
    bin_name = sketch_name + ".cms"
    arch_name = sketch_name + ".gz"
    sketch_path = os.path.join(fress_outdir, sketch_name)
    arch_path = os.path.join(fress_outdir, arch_name)

    L1, dim, max_val = run_fress_cms(kmcdb, sketch_path, epsilon)
    L1 = int(L1)
    dim = int(dim)
    max_val = int(max_val)
    ncolls, sod, avgd, maxd = run_fress_cmschk(kmcdb, sketch_path)
    ncolls = int(ncolls)
    avgd = float(avgd)
    maxd = int(maxd)

    theoretical_udim = round(dim * math.ceil(math.log(max_val, 2)) / 8)
    compress(fress_outdir, [bin_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{}".format(filename, epsilon, k, "NA", round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim, "NA")

def run_bbhash_for(fastx: str, k: int, _, kmc_outdir:str, fress_outdir: str, tmpdir: str, max_mem: int):
    kmc.count(k, fastx, kmc_outdir, tmpdir, max_mem, True)

    filename, _, _, _, kmcdb = kmc.getKMCPaths(k, fastx, kmc_outdir)
    sketch_name = "{}k{}".format(filename, k)
    mphf_name = sketch_name + ".bbh"
    payload_name = sketch_name + ".pld"
    arch_name = sketch_name + "_BBH.gz"
    sketch_path = os.path.join(fress_outdir, sketch_name)
    mphf_path = os.path.join(fress_outdir, mphf_name)
    arch_path = os.path.join(fress_outdir, arch_name)

    max_val, L0 = run_fress_bbhash(kmcdb, sketch_path)
    max_val = int(max_val)
    L0 = int(L0)

    mphf_size = os.stat(mphf_path).st_size
    theoretical_udim = round(L0 * math.ceil(math.log(max_val, 2)) / 8) + mphf_size
    compress(fress_outdir, [mphf_name, payload_name], arch_path)
    cdim = os.stat(arch_path).st_size
    os.remove(arch_path)
    return "{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{}".format(filename, "NA", k, "NA", 0, 0, 0, 0, theoretical_udim, cdim, mphf_size)

def run_combination(description_file: str, output_file: str, kmc_outdir: str, fress_outdir: str, tmpdir: str, max_mem: int, command):
    """Run fress for multiple parameters

    The input file must contain one line for each dataset in the following format:
    <dataset name> [k1, ..., kn] [epsilon1, ..., epsilon<m>]
    one line for each combination of parameter will be generated.
    """
    with open(output_file, "w") as oh:
        with open(description_file, "r") as dh:
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
                        oh.write(command(dataset, k, e, kmc_outdir, fress_outdir, tmpdir, max_mem) + "\n")
                        oh.flush()

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

    parser_cmsm = subparsers.add_parser("cmsm", help="Run set-min sketch pipeline for multiple values.")
    parser_cmsm.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_cmsm.add_argument("-o", help="output file (tsv format)", required=True)
    parser_cmsm.add_argument("-c", help="kmc folder", required=True)
    parser_cmsm.add_argument("-f", help="fress output folder", required=True)
    parser_cmsm.add_argument("-w", help="tmp dir", required=True)
    parser_cmsm.add_argument("-m", help="max memory for kmc", type=int)

    parser_bbhm = subparsers.add_parser("bbhm", help="Run set-min sketch pipeline for multiple values.")
    parser_bbhm.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_bbhm.add_argument("-o", help="output file (tsv format)", required=True)
    parser_bbhm.add_argument("-c", help="kmc folder", required=True)
    parser_bbhm.add_argument("-f", help="fress output folder", required=True)
    parser_bbhm.add_argument("-w", help="tmp dir", required=True)
    parser_bbhm.add_argument("-m", help="max memory for kmc", type=int)

    args = parser.parse_args(sys.argv)
    if (args.command == "sms"): print(run_sms_for(args.file, args.k, args.e, args.c, args.f, args.w, args.m))
    elif (args.command == "smsm"): run_combination(args.file, args.o, args.c, args.f, args.w, args.m, run_sms_for)
    elif (args.command == "cmsm"): run_combination(args.file, args.o, args.c, args.f, args.w, args.m, run_cms_for)
    elif (args.command == "bbhm"): run_combination(args.file, args.o, args.c, args.f, args.w, args.m, run_bbhash_for)
    else: parser.print_help(sys.stderr)
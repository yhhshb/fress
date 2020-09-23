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
    if(out.returncode != 0): raise Exception("Error while building the sketch")
    return out.stdout.decode("utf-8").split()

def run_fress_check(kmc_name: str, sketch_name):
    out = subprocess.run([fress, "check", "-i", kmc_name, "-d", sketch_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("Error while building the sketch")
    return out.stdout.decode("utf-8").split()

def run_fress_for(fastx: str, k: int, epsilon: float, kmc_outdir:str, fress_outdir: str, tmpdir: str, max_mem: int):
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
    bin_path = os.path.join(fress_outdir, bin_name)
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
    return "{}\t{}\t{}\t{:.2f}\t{}\t{}\t{:.4f}\t{}\t{}\t{}".format(filename, epsilon, k, skewness, round(L1 * epsilon), sod, avgd, maxd, theoretical_udim, cdim)

def run_fress_combination(description_file: str, output_file: str, kmc_outdir: str, fress_outdir: str, tmpdir: str, max_mem: int):
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
                epsilons = [float(block.strip()) for block in blocks[2].split(',')]
                for e in epsilons:
                    for k in ks:
                        sys.stderr.write("Now fressing {} with k = {} and epsilon = {}\n".format(dataset, k, e))
                        oh.write(run_fress_for(dataset, k, e, kmc_outdir, fress_outdir, tmpdir, max_mem) + "\n")
                        oh.flush()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_single = subparsers.add_parser("single", help="Run pipeline for fixed values. This generates only one line of the final result table")
    parser_single.add_argument("file", help="fastx file to process", type=str)
    parser_single.add_argument("-k", help="k-mer value", type=int)
    parser_single.add_argument("-e", help="epsilon", type=float)
    parser_single.add_argument("-c", help="kmc folder")
    parser_single.add_argument("-f", help="fress output folder")
    parser_single.add_argument("-w", help="tmp dir")
    parser_single.add_argument("-m", help="max memory for kmc", type=int)

    parser_multiple = subparsers.add_parser("multiple", help="Run pipeline for multiple values.")
    parser_multiple.add_argument("file", help="File containing one line per dataset in the form of <dataset> [k1, ...,kn] [epsilon1, ..., epsilon<m>]")
    parser_multiple.add_argument("-o", help="output file (tsv format)")
    parser_multiple.add_argument("-c", help="kmc folder")
    parser_multiple.add_argument("-f", help="fress output folder")
    parser_multiple.add_argument("-w", help="tmp dir")
    parser_multiple.add_argument("-m", help="max memory for kmc", type=int)

    args = parser.parse_args(sys.argv)
    if (args.command == "single"): print(run_fress_for(args.file, args.k, args.e, args.c, args.f, args.w, args.m))
    elif (args.command == "multiple"): run_fress_combination(args.file, args.o, args.c, args.f, args.w, args.m)
    else: parser.print_help(sys.stderr)
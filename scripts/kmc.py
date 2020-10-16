#!/usr/bin/python3

import os
import sys
import subprocess
import kmerencdec

def getKMCDatabaseNames(path: str):
    kmc_files = [os.path.join(path, str(os.path.basename(f)).split('.')[0]) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.endswith("kmc_suf")]
    kmc_files.sort(key = lambda v : v.upper())
    return kmc_files

def getWMHSketches(path: str):
    wmh_files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.endswith(".wmh")]
    wmh_files.sort(key = lambda v : v.upper())
    return wmh_files

def getKMCPaths(k: int, input_file: str, output_path: str):
    _, filename = kmerencdec.path_leaf(input_file)
    filename = filename.split('.')[0]
    
    cwd = os.getcwd()
    kfolder = os.path.join(cwd, output_path)
    kfolder = os.path.join(kfolder, "k"+str(k))
    intermidiate_name = "dummy_" + filename
    intermidiate_file = os.path.join(kfolder, intermidiate_name)
    sorted_name = "S_" + filename
    sorted_file = os.path.join(kfolder, sorted_name)
    final_name = filename
    final_file = os.path.join(kfolder, final_name)
    return filename, kfolder, intermidiate_file, sorted_file, final_file

def count(k: int, input_file: str, output_path: str, kmc_dummy_folder: str, mmemory: int, unsorted: bool):
    """Call kmc to count all k-mers.
    
    The database in output is sorted by k-mer.
    """
    _, kfolder, intermidiate_file, sorted_file, final_file = getKMCPaths(k, input_file, output_path)
    fmt = kmerencdec.is_fastx(input_file)
    if fmt:
        fmt, _ = fmt
        if fmt == "fasta": fmt = "-fm"
        elif fmt == "fastq": fmt = "-fq"
        else: raise ValueError("Unrecognized file type")
    else:
        raise ValueError("Input file not in fastx format")
    try:
        os.makedirs(kfolder)
    except FileExistsError:
        pass
    this_path = os.path.dirname(os.path.abspath(__file__))
    kmc = os.path.join(this_path, "kmc")
    kmc_tools = os.path.join(this_path, "kmc_tools")
    subprocess.run([kmc, "-b", "-k"+str(k), "-ci0", "-cs4294967295", "-cx4294967295", "-m"+str(mmemory), fmt, input_file, intermidiate_file, kmc_dummy_folder], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if (not unsorted): subprocess.run([kmc_tools, "transform", intermidiate_file, "sort", sorted_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if os.path.exists(sorted_file+".kmc_pre"): os.remove(intermidiate_file+".kmc_pre")
    else: os.rename(intermidiate_file+".kmc_pre", sorted_file+".kmc_pre")
    if os.path.exists(sorted_file+".kmc_suf"): os.remove(intermidiate_file+".kmc_suf")
    else: os.rename(intermidiate_file+".kmc_suf", sorted_file+".kmc_suf")
    os.replace(sorted_file+".kmc_pre", final_file+".kmc_pre")
    os.replace(sorted_file+".kmc_suf", final_file+".kmc_suf")

def countFolder(infolder: str, outfolder: str, ks: list, kmc_dummy_folder: str, mmemory: int, unsorted: bool):
    """Build kmc databases for the input folder for all k values.
    
    The output folder will have a folder of the same name as the input with multiple folders named as k<value>.
    This function also creates the temporary folder for kmc hard-disk computations.
    """
    fastxs = kmerencdec.getFastx(infolder)
    _, basename = kmerencdec.path_leaf(infolder)
    working_folder = os.path.join(outfolder, basename)
    try: os.makedirs(working_folder)
    except FileExistsError: pass
    try: os.makedirs(kmc_dummy_folder)
    except FileExistsError: pass
    for k in ks:
        for fpath in fastxs:
            count(k, fpath, working_folder, kmc_dummy_folder, mmemory, unsorted)

def wmhSketch(kmc_name: str, s: int, out_name: str):
    this_path = os.path.dirname(os.path.abspath(__file__))
    i2cws = os.path.join(this_path, "cws")
    if(not out_name.endswith(".wmh")): out_name = out_name + ".wmh"
    ret = subprocess.run([i2cws, "sketch", "-i", kmc_name, "-o", out_name, "-s", str(s)], stdout = subprocess.PIPE)
    if(ret.returncode != 0): raise Exception("subprocess error")

def sketchFolder(kmcfolder: str, s: int, outfolder: str):
    """Build sketches for all kmc databases in the given folder.
    
    The sketches will be saved in the output folder with the same basename of each file and extension wmh
    """
    kmcs = getKMCDatabaseNames(kmcfolder)
    for db in kmcs:
        wmhSketch(db, s, os.path.join(outfolder, os.path.basename(db)) + ".wmh")

def histogram(kmc_name: str, output_file: str):
    filepath, filename = kmerencdec.path_leaf(kmc_name)
    filename = filename.split('.')[0]
    input_file = os.path.join(filepath, filename)
    this_path = os.path.dirname(os.path.abspath(__file__))
    kmc_tools = os.path.join(this_path, "kmc_tools")
    subprocess.run([kmc_tools, "transform", input_file, "histogram", output_file])

def dumpWindows(k: int, w: int):
    """Print count vectors for non-overlapping windows
    """
    from Bio import SeqIO
    fd = sys.stdin
    start = 0
    end = w
    for record in SeqIO.parse(fd, "fasta"):
        seq = str(record.seq)
        while(end <= len(seq)):
            kmc = kmerencdec.kmercount(seq[start:end], k, kmerencdec.enc_nuc_table_extended)
            ffp = kmerencdec.kmc2ffp(kmc, k)
            print(str(ffp).replace('\n', ''))
            start = end
            end += w
        if(start != len(seq) - w):
            kmc = kmerencdec.kmercount(seq[start:], k, kmerencdec.enc_nuc_table_extended)
            ffp = kmerencdec.kmc2ffp(kmc, k)
            print(str(ffp).replace('\n', ''))

def JWJ(path1: str, path2: str):
    """Compute Jaccard and Weighted Jaccard between two kmc outputs."""
    this_path = os.path.dirname(os.path.abspath(__file__))
    i2cws = os.path.join(this_path, "cws")
    ret = subprocess.run([i2cws, "full", path1, path2], stdout = subprocess.PIPE)
    if(ret.returncode == 0):
        return [float(x) for x in ret.stdout.decode("utf8").split('|')]
    else: 
        raise Exception("subprocess error")

def eWJ(path1:str , path2: str):
    """Comapre two weighted minHash sketches"""
    this_path = os.path.dirname(os.path.abspath(__file__))
    i2cws = os.path.join(this_path, "cws")
    ret = subprocess.run([i2cws, "compare", path1, path2], stdout = subprocess.PIPE)
    if(ret.returncode == 0):
        return float(ret.stdout.decode("utf8"))
    else: 
        raise Exception("subprocess error")

def minHash(ref: str, query: str, k: int, s: int):
    """Compute minHash similarity and MASH distance by calling MASH."""
    out = subprocess.run(["mash", "dist", "-k", str(k), "-s", str(s), ref, query], encoding = 'utf8', stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    if(out.returncode != 0): raise Exception("mash subprocess error")
    ol = out.stdout.split('\t')
    iu = ol[4].split('/')
    return int(iu[0])/int(iu[1]), ol[2]

def checkFastx(filepath: str, outfolder: str, k: int):
    """check if the given file hash been kmced in the given working directory"""
    path, filename = kmerencdec.path_leaf(filepath)
    _, dataset_name = kmerencdec.path_leaf(path)
    filename = filename.split('.')[0]
    #foldername = os.path.basename(path)
    working_folder = os.path.join(outfolder, dataset_name)
    k_folder = os.path.join(working_folder, "k"+str(k))
    kmc_name = os.path.join(k_folder, filename)
    kmc_pre = kmc_name + ".kmc_pre"
    kmc_suf = kmc_name + ".kmc_suf"
    if os.path.exists(kmc_pre) and os.path.exists(kmc_suf): return kmc_name
    else: return ""

def setup_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_count = subparsers.add_parser("count", help = "Count k-mers")
    parser_count.add_argument("file", help = "path to a fasta file", type=str)
    parser_count.add_argument("-k", help = "k-mer size", type=int, required = True)
    parser_count.add_argument("-m", help = "Maximum RAM in GB", type=int, default=6)
    parser_count.add_argument("-w", help = "working path if disk is used", type=str, required=True)
    parser_count.add_argument("-o", help = "output path", type=str, required=True)
    parser_count.add_argument("-y", help = "disable sorting of the output database", action="store_true")

    parser_all = subparsers.add_parser("all", help = "Count for all files in given folder and for different ks")
    parser_all.add_argument("folder", help = "Folder with FASTX files")
    parser_all.add_argument("-l", help = "min k", type=int, required=True)
    parser_all.add_argument("-u", help = "max k", type=int, required=True)
    parser_all.add_argument("-m", help = "Maximum RAM in GB", type=int, default=6)
    parser_all.add_argument("-w", help = "working path if disk is used", type=str, required=True)
    parser_all.add_argument("-o", help = "output path", type=str, required=True)
    parser_all.add_argument("-y", help = "disable sorting of the output databases", action="store_true")

    parser_cws = subparsers.add_parser("cws", help = "Sketch a KMC database using the CWS algorithm (by sketching pairs of (k-mer, id) to discriminate the repetitions)")
    parser_cws.add_argument("-i", help = "input file")
    parser_cws.add_argument("-o", help = "output file (.wmh extension added automatically if not specified)")
    parser_cws.add_argument("-s", help = "sketch size", type=int)
    
    parser_cwsall = subparsers.add_parser("cwsall", help = "sketch all KMC databases inside a folder")
    parser_cwsall.add_argument("-i", help = "folder containing the KMC databases")
    parser_cwsall.add_argument("-o", help = "output folder")
    parser_cwsall.add_argument("-s", help = "sketch size (common for all the output sketches)", type=int)

    parser_wdist = subparsers.add_parser("dumpw", help = "print the count vectors of non-overlapping windows")
    #parser_wdist.add_argument("file", help = "path to a fasta file", type=str)
    parser_wdist.add_argument("-k", help = "k-mer length", type=int, required=True)
    parser_wdist.add_argument("-w", help = "window size", type=int, required=True)

    return parser

def main(args):
    if   (args.command == "count"): count(args.k, args.file, args.o, args.w, args.m, args.y)
    elif (args.command == "all"): countFolder(args.folder, args.o, range(args.l, args.u), args.w, args.m, args.y)
    elif (args.command == "cws"): wmhSketch(args.i, args.s, args.o)
    elif (args.command == "cwsall"): sketchFolder(args.i, args.s, args.o)
    elif (args.command == "dumpw"): dumpWindows(args.k, args.w) 
    else: parser.print_help(sys.stderr)

if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args(sys.argv)
    main(args)

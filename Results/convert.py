#!/usr/bin/python3

"""
Post-processing script for the paper.
BE CAREFUL: the column format returned by fress is different now.
This script primary purpose is:
1) to merge together all the intermidiate tables (asms + check + info + kmc), (cms + kmc)
2) make two tables based on the epsilon values (0.01, 0.001)
3) save everything in csv format for LaTeX
"""

import functools
import pandas as pd

smsdf = pd.read_table("asms_pipeline.out.tsv", sep='\t', header=None)
checkdf = pd.read_table("check_pipeline.out.tsv", sep='\t', header=None)
infodf = pd.read_table("info_pipeline.out.tsv", sep='\t', header=None)
kmcdf = pd.read_table("kmc_pipeline.out.tsv", sep='\t', header=None)
cmsdf = pd.read_table("cms_pipeline.out.tsv", sep='\t', header=None)
bbhdf = pd.read_table("bbhash_pipeline.out.tsv", sep='\t', header=None)

smsdf.columns = ["name", "epsilon", "k", "skew", "threshold", "ssod", "savg", "smax", "ssize", "scsize"]
checkdf.columns = ["name", "epsilon", "k", "ncolls", "ntrue", "sod", "avg", "max"]
infodf.columns = ["name", "epsilon", "k", "R", "B", "RB"]
kmcdf.columns = ["name", "k", "kmc"]
cmsdf.columns = ["name", "epsilon", "k", "skew", "threshold", "csod", "cavg", "cmax", "csize", "ccsize"]
bbhdf.columns = ["name", "k", "mphfsize", "bsize", "bcsize"]

smsdf = smsdf.drop("skew", 1)
smsdf["name"] = smsdf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
smsdf["name"] = smsdf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

checkdf = checkdf.drop("sod", 1)
checkdf = checkdf.drop("avg", 1)
checkdf = checkdf.drop("max", 1)
checkdf["name"] = checkdf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
checkdf["name"] = checkdf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

infodf["name"] = infodf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
infodf["name"] = infodf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

kmcdf["name"] = kmcdf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
kmcdf["name"] = kmcdf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

cmsdf = cmsdf.drop("skew", 1)
cmsdf["name"] = cmsdf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
cmsdf["name"] = cmsdf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

bbhdf["name"] = bbhdf["name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
bbhdf["name"] = bbhdf["name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)

ldfs = [smsdf, checkdf, infodf]

joined = functools.reduce(lambda  left,right: pd.merge(left, right, on=["name", "epsilon", "k"], how="outer"), ldfs)
joined = pd.merge(joined, cmsdf, on=["name", "epsilon", "k", "threshold"], how="outer")
joined = pd.merge(joined, kmcdf, on=["name", "k"], how="outer")
joined = pd.merge(joined, bbhdf, on=["name", "k"], how="outer")
joined["tavg"] = joined["ssod"] / joined["ntrue"]
joined["savg"] = joined["savg"].map(lambda x: '{:.2f}'.format(x))
joined["cavg"] = joined["cavg"].map(lambda x: '{:.2f}'.format(x))
joined["tavg"] = joined["tavg"].map(lambda x: '{:.2f}'.format(x))
joined = joined[["name", "epsilon", "k", "kmc", "R", "B", "RB", "ssize", "scsize", "threshold", "ncolls", "ntrue", "ssod", "savg", "tavg", "smax", "csod", "cavg", "cmax", "csize", "ccsize"]]
joined.to_csv("fress_all.csv", header=True, index=False)
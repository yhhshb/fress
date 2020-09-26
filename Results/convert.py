#!/usr/bin/python3

import pandas as pd

df = pd.read_table("additive_fress_results.tsv", sep='\t', header=None)
df.columns = ["Name", "epsilon", "k", "skew", "threshold", "dsum", "avg", "max", "size", "csize"]
df = df.drop("skew", 1)
df["avg"] = df["avg"].map(lambda x: '{:.2f}'.format(x))
df["Name"] = df["Name"].map(lambda x: "USakai" if x == "U_Sakai" else x)
df["Name"] = df["Name"].map(lambda x: "SRR622461" if x == "SRR622461_1" else x)
df.to_csv("fress_results.csv", header=True, index=False)
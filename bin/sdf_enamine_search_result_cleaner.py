#! /home/siu/anaconda3/envs/rdkit/bin/python

import argparse
import os
import rdkit
import pandas as pd
from rdkit.Chem import PandasTools
import glob

parser = argparse.ArgumentParser(description = 'Read all result-search*.sdf and clean and merge them')
args = parser.parse_args()

sdf_list = glob.glob('result-search*.sdf')

for sdf in sdf_list:
    with open(sdf, 'r') as f:
        lines = f.readlines()

    with open(sdf, 'w') as f:
        for line in lines:
            if "> <" in line:
                line = line.replace("> <", "\n> <")
            elif ">  <" in line:
                line = line.replace(">  <", "\n> <")
            elif "$$$$" in line:
                line = line.replace("$$$$", "\n$$$$")
            f.write(line)

    f.close()

dfs = []
for sdf in sdf_list:
    df = PandasTools.LoadSDF(sdf)
    dfs.append(df)

cdf = pd.concat(dfs, axis=0, ignore_index=True)
rdkit.Chem.PandasTools.WriteSDF(cdf, 'result-search_merged.sdf', molColName='ROMol', idName='CatalogID', properties=list(cdf.columns))

print('result-search_merged.sdf')
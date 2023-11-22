#! /home/siu/anaconda/envs/chem/bin/python

import pandas as pd
import os
import argparse
import rdkit
from rdkit import Chem
from sys import exit
from myrdkit import neutralize_atoms

parser = argparse.ArgumentParser(description = 'Write columns names on a text file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument("--smiles", "-s", required=True, help='smiles column') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv
fn = os.path.splitext(csv)[0]
smiles = args.smiles

df = pd.read_csv(csv)

# mols 생성에 실패하면 none이 추가되는데 그걸 제외하고 길이 비교하여 검증
mols = [Chem.MolFromSmiles(x) for x in df[smiles]]
mols = list(filter(None, mols))

if len(df[smiles]) != len(mols):
    print('Invalid smiles are found')
    exit()

# 
n_mols = [neutralize_atoms(x) for x in mols]
n_smiles = [Chem.MolToSmiles(x, kekuleSmiles=True) for x in n_mols]

df[f"{smiles}_copy"] = df[smiles]
df[smiles] = n_smiles

df.to_csv(f"{fn}_nt.csv", index=False)

print("<filen name>_nt.csv is generated")



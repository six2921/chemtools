#! /home/siu/anaconda/bin/python

import pandas as pd
import argparse
import rdkit
from rdkit import Chem

parser = argparse.ArgumentParser(description = 'Write columns names on a text file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument("--smiles", "-s", required=True, help='smiles column') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv
col = args.smiles

df = pd.read_csv(csv)
df['valid_smi'] = [None if pd.isna(x) else Chem.MolFromSmiles(x) for x in df['smiles']]

df.to_csv(csv, index=False)         

count = df['valid_smi'].isnull().sum()
print('')
print(f"{count} smiles is not valid. You can find it in 'valid_smi' column.")

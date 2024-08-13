#! /home/siu/anaconda/bin/python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'extract single column and write it on a smi file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument('titles', help='title column') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv
smiles = args.titles

df = pd.read_csv(csv)
f = open("smiles.smi", 'w')
for i in df[smiles]:
    f.write(f'{i}\n')
f.close()

print("The smiles.smi is generated ")

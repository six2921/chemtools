#! /home/siu/anaconda3/envs/rdkit/bin/python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'extract all values for a column. ex) To make smi file from smiles column')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument('titles', help='title column') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv
titles = args.titles

df = pd.read_csv(csv)
f = open("titles.txt", 'w')
for i in df[titles]:
    f.write(f'{i}\n')
f.close()

print("The titles.txt is generated ")

#! /home/siu/anaconda3/envs/rdkit/bin/python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Write columns names on a text file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv

df = pd.read_csv(csv)
f = open("features.txt", 'w')
for i in df.columns:
    f.write(f'{i}\n')
f.close()

print("The features.txt is generated ")

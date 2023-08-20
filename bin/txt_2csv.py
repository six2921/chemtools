#! /home/siu/anaconda3/envs/rdkit/bin/python

import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description = 'Convert txt(from Data Worrior) into csv')
parser.add_argument('txt', help='txt file, seperator ; is used.') # 필요한 인수를 추가
args = parser.parse_args()

txt = args.txt
df = pd.read_csv(txt, sep='\t')

for col in df.columns:
    if 'Structure of' in col:
        del df[col]
    
    if 'Structure [idcode]' in col:
        del df[col]

    if 'Unnamed: ' in col:
        del df[col]

    if 'SmilesFragFp' in col:
        del df[col]

name = txt.split('.')[0]+'.csv'
df.to_csv(name, index=False)

print('Converted')

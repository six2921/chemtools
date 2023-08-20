#! /home/siu/anaconda3/envs/rdkit/bin/python

import argparse
import os
import pandas as pd
import rdkit
from rdkit.Chem import PandasTools


parser = argparse.ArgumentParser(description = 'convert maegz to sdf and csv')
parser.add_argument('sdf', metavar='sdf', help='input .sdf files')
parser.add_argument('csv', metavar='csv', help='input .csv files')
args = parser.parse_args()

sdf = args.sdf
csv = args.csv

df1 = PandasTools.LoadSDF(sdf, molColName='Molecule')
for col in df1.columns:
    try:
        print(col, '|', df1[col][0][:10])
    except:
        pass

key1 = 'dZ1%mO'
while key1 not in list(df1.columns):
    key1 = input('첫번째 파일에 키 컬럼: ')
print(' ')
print(' ')

df2 = pd.read_csv(csv)
for col in df2.columns:
    try:
        print(col, '|', df2[col][0][:10])
    except:
        pass

key2 = 'dZ1%mO'
while key2 not in list(df2.columns):
    key2 = input('두번째 파일에 키 컬럼: ')
print(' ')
print(' ')

# how = input('병합 방식 [outer, inner, left, right]: ') # 무조건 left (sdf)에 넣어야 함

df = pd.merge(df1, df2, left_on=key1, right_on=key2, how='left', suffixes=('', '_1'))

fn = os.path.splitext(sdf)[0]
rdkit.Chem.PandasTools.WriteSDF(df, fn+'_merged.sdf', molColName='Molecule', idName='compound_id', properties=df.columns, allNumeric=False)

print('_merged.sdf is generated')

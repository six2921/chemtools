#! /home/siu/anaconda3/envs/rdkit/bin/python

import argparse
import os
import rdkit
from rdkit.Chem import PandasTools

parser = argparse.ArgumentParser(description = 'convert maegz to sdf and csv')
parser.add_argument('sdf', help='input .sdf files')
args = parser.parse_args()

file = args.sdf

df = PandasTools.LoadSDF(sdf, molColName='Molecule')
for col in df.columns:
    try:
        print(col, '|', df[col][0][:10])
    except:
        pass

key = 'dZ1%mO'
while key not in list(df.columns):
    key = input('첫번째 파일에 키 컬럼: ')
print(' ')

fn = os.path.splitext(sdf)[0]
rdkit.Chem.PandasTools.WriteSDF(df, fn+'_new.sdf', molColName='Molecule', idName='compound_id', properties=df.columns, allNumeric=False)

print('_new.sdf is generated')

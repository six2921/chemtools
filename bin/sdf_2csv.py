#! /home/siu/anaconda/envs/chem/bin/python

import argparse
import os
from rdkit.Chem import PandasTools


parser = argparse.ArgumentParser(description = 'convert maegz to sdf and csv')
parser.add_argument('sdf', metavar='sdf', help='input .sdf files')
args = parser.parse_args()

sdf = args.sdf

df = PandasTools.LoadSDF(sdf, smilesName='smiles_rdkit', molColName='molecule_rdkit')

# smiles_rdkit으로 생성된 컬럼이름 변경
if 'smiles' in df.index:
    pass
else:
    df.rename(columns={'smiles_rdkit':'smiles'})

# 컬럼 프린트
for col in df.columns:
    try:
        print(col, '|', df[col][0][:10])
    except:
        pass

# 입력 스크립트
name = 'dZ1%mO'
while name not in list(df.columns):
    name = input('Write name column: ')
print(' ')
print(' ')

df.drop('molecule_rdkit', axis=1, inplace=True) # molecule 컬럼 삭제

fn = os.path.splitext(sdf)[0]
df.to_csv(fn+'.csv', index=False)

print('sdf is converted into csv')

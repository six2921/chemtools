#! /home/siu/anaconda/envs/chem/bin/python

import argparse
from rdkit.Chem import PandasTools
import pandas as pd

# argparse 설정
parser = argparse.ArgumentParser(description='You can generate subset sdf file by input compound code number (ex. 0153 1125 1535 -> Copy PF-0153 PF-1125 PF-1535)')
parser.add_argument('sdf', type=str, help='sdf file')
args = parser.parse_args()

sdf = args.sdf

# 파일 읽기
df = PandasTools.LoadSDF(sdf)
print("Total number of molecules: ", len(df))
pd.set_option('display.max_columns', None)
print(df.head(3))

# id 컬럼을 복사해서, 화합물 코드 앞에 prefix를 제거하고 정수로 변환
id_column = input("\n Enter the column name for cpd ID: ")
df['code_temp1127'] =  df[id_column].apply(lambda x: int(x.split('-')[1]))

sub_merge = pd.DataFrame()
while 1:
    # 화합물 번호 입력받고 정수로 변환
    codes = input("ex) 0153 1125 1535 / Write 'exit' to exit: ")
    if codes == 'exit':
        break
    else:
        codes = codes.split()
        try:
            codes = [int(code) for code in codes]
            sub = df.loc[df['code_temp1127'].isin(codes)]
            sub_merge = pd.concat([sub_merge, sub])
            print('TOTAL:', len(list(set(sub_merge))), sub_merge[id_column].unique(), )
        except ValueError:
            print("<Please enter numbers only>")
            continue

sub_merge = sub_merge.drop_duplicates(subset=[id_column])

if len(sub_merge) > 0:
    # 파일 쓰기
    new_filename = sdf.split('.')[0] + "_subset.sdf"
    PandasTools.WriteSDF(sub_merge, new_filename,  idName=id_column, properties=list(df.columns))
    print(f"{sub_merge[id_column].unique()} cpds are written to {new_filename}")
#! /home/siu/anaconda3/bin/python

import pandas as pd
import os
from pathos.multiprocessing import ProcessingPool as Pool
import argparse
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    return None

def calculate_molecular_weight_with_index(smiles):
    index, smiles = smiles
    mw = calculate_molecular_weight(smiles)
    return index, mw

def update_result(results, df):
    for index, mw in results:
        df.at[index, 'MW'] = mw

# Argument parsing
parser = argparse.ArgumentParser(description="Calculate molecular weight from SMILES in a CSV file.")
parser.add_argument("input_csv", help="Path to the input CSV file.")
args = parser.parse_args()

# CSV 파일 읽기
df = pd.read_csv(args.input_csv, low_memory=False)

# NaN 값을 빈 문자열로 대체
df['smiles'] = df['smiles'].fillna('')

# Progress bar 설정
pbar = tqdm(total=len(df['smiles']), desc="Processing SMILES", ncols=100)

def calculate_molecular_weight_with_progress(smiles):
    result = calculate_molecular_weight_with_index(smiles)
    pbar.update()
    return result

# 병렬 처리 설정
with Pool() as pool:
    results = pool.map(calculate_molecular_weight_with_progress, enumerate(df['smiles']))

# Progress bar 닫기
pbar.close()

# 결과를 DataFrame에 추가
update_result(results, df)

# 출력 파일 이름 생성
base, ext = os.path.splitext(args.input_csv)
output_csv = f"{base}-MW{ext}"

# 결과를 새로운 CSV 파일로 저장
df.to_csv(output_csv, index=False)
print(f"분자량이 '{output_csv}' 파일에 기록되었습니다.")


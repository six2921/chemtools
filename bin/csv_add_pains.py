#! /home/siu/anaconda/bin/python

import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog
import os
from pathos.multiprocessing import ProcessingPool as Pool
import argparse
from tqdm import tqdm
from rdkit import RDLogger

# RDKit 경고 억제
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# PAINS 필터 카탈로그 생성
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

def filter_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    entry = catalog.GetFirstMatch(mol)
    if entry is not None:
        return f"PAINS match found: {entry.GetDescription()}"
    return "No PAINS match"

def filter_molecule_with_index(smiles):
    index, smiles = smiles
    result = filter_molecule(smiles)
    return index, result

def update_result(results, df):
    for index, result in results:
        df.at[index, 'PAINS_Result'] = result

# Argument parsing
parser = argparse.ArgumentParser(description="Apply PAINS filter to SMILES in a CSV file.")
parser.add_argument("input_csv", help="Path to the input CSV file.")
args = parser.parse_args()

# CSV 파일 읽기
df = pd.read_csv(args.input_csv, low_memory=False)

# NaN 값을 빈 문자열로 대체
df['smiles'] = df['smiles'].fillna('')

# Progress bar 설정
pbar = tqdm(total=len(df['smiles']), desc="Processing SMILES", ncols=100)

def filter_molecule_with_progress(smiles):
    result = filter_molecule_with_index(smiles)
    pbar.update()
    return result

# 병렬 처리 설정
with Pool() as pool:
    results = pool.map(filter_molecule_with_progress, enumerate(df['smiles']))

# Progress bar 닫기
pbar.close()

# 결과를 DataFrame에 추가
update_result(results, df)

# 출력 파일 이름 생성
base, ext = os.path.splitext(args.input_csv)
output_csv = f"{base}-PAINS{ext}"

# 결과를 새로운 CSV 파일로 저장
df.to_csv(output_csv, index=False)
print(f"PAINS 필터 적용 결과가 '{output_csv}' 파일에 저장되었습니다.")

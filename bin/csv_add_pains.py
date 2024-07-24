#! /home/siu/anaconda/bin/python

import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog
import argparse
import os

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

def main(input_csv):
    # CSV 파일 읽기
    df = pd.read_csv(input_csv)

    # PAINS 필터 적용
    df['PAINS_Result'] = df['smiles'].apply(filter_molecule)

    # 출력 파일 이름 생성
    base, ext = os.path.splitext(input_csv)
    output_csv = f"{base}-PAINS{ext}"

    # 결과를 새로운 CSV 파일로 저장
    df.to_csv(output_csv, index=False)
    print(f"PAINS 필터 적용 결과가 '{output_csv}' 파일에 저장되었습니다.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Apply PAINS filter to SMILES in a CSV file.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    
    args = parser.parse_args()
    
    main(args.input_csv)


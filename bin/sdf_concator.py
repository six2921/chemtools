#! /home/siu/anaconda/bin/python

from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser(description = 'concont sdf files in folder')

def read_sdf_files(file_list):
    sdf_data = []
    for file in file_list:
        suppl = Chem.SDMolSupplier(file)
        sdf_data.extend([mol for mol in suppl if mol is not None])
    return sdf_data

def add_duplicate_column(mol_list):
    df = PandasTools.LoadSDF(mol_list, molColName='ROMol')
    smiles = df['ROMol'].apply(Chem.MolToSmiles)
    df['duplicated'] = smiles.duplicated().astype(int)
    unique_smiles = smiles.drop_duplicates().reset_index()
    unique_smiles['duplicated'] = unique_smiles.index + 1
    duplicate_dict = dict(zip(unique_smiles['ROMol'], unique_smiles['duplicated']))
    df['duplicated'] = df['ROMol'].apply(lambda mol: duplicate_dict[Chem.MolToSmiles(mol)])
    return df

def save_to_sdf(df, output_file):
    PandasTools.WriteSDF(df, output_file, properties=list(df.columns))

if __name__ == "__main__":
    # SDF 파일 목록 가져오기
    sdf_files = glob.glob("*.sdf")

    # SDF 파일 읽기
    mol_list = read_sdf_files(sdf_files)

    # 중복 구조에 대해 'duplicated' 컬럼 추가
    df = add_duplicate_column(mol_list)

    # 결과를 새로운 SDF 파일로 저장
    save_to_sdf(df, "output.sdf")


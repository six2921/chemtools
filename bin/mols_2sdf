#! /home/siu/anaconda/bin/python

import os
from rdkit import Chem
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Move to directory containing mol files and run this. Merge them into sdf")

os.chdir(os.getcwd())
folder_path = os.getcwd()
mollist = os.listdir(folder_path)

mol_dict = {}
for filename in mollist:
    if filename.endswith(".mol"):
        mol = Chem.MolFromMolFile(filename)
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            name = filename.split(".mol")[0]  # 확장자 제거
            mol_dict[filename] = {
                "name": name, 
                "mol": mol,
                "smiles": smiles
            }

df = pd.DataFrame.from_dict(mol_dict, orient="index")

# SDF 파일 작성
with Chem.SDWriter("output.sdf") as writer:
    for index, row in df.iterrows():
        mol = row["mol"]
        mol.SetProp("_Name", row["name"])  # 이름 추가
        mol.SetProp("SMILES", row["smiles"])  # SMILES 추가
        writer.write(mol)
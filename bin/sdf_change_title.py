#! /home/siu/anaconda/envs/chem/bin/python

import argparse
import os

from rdkit import Chem

# argparse 설정
parser = argparse.ArgumentParser(description='Rename molecules in an SDF file based on a property value')
parser.add_argument('filename', metavar='filename', type=str, help='the name of the SDF file')
args = parser.parse_args()

filename = args.filename

# 파일 읽기
suppl = Chem.SDMolSupplier(filename)

# 첫번째 분자의 프로퍼티와 값 출력
mol = suppl[0]
for prop in mol.GetPropNames():
    print(f"{prop}: {mol.GetProp(prop)}")

# 사용자 입력 받기
prop_name = input("Enter the property name to copy to the molecule name: ")

# 분자 이름 변경
new_suppl = []
suppl = Chem.SDMolSupplier(filename)
for i, mol in enumerate(suppl):
    print(f"Processing molecule", {i}, end="\r")
    if mol is None:
        print(f"Error: Failed to read molecule at index {i}")
        continue
    try:
        prop_value = mol.GetProp(prop_name)
        mol.SetProp("_Name", prop_value)
        new_suppl.append(mol)
    except KeyError:
        print(f"Property {prop_name} not found for molecule {i+1}")

# 파일 쓰기
new_filename = os.path.splitext(filename)[0] + "_newtitle.sdf"
writer = Chem.SDWriter(new_filename)
for mol in new_suppl:
    writer.write(mol)
writer.close()
print(f"New molecule names written to {new_filename}")
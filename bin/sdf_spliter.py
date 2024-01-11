#! /home/siu/anaconda/envs/chem/bin/python

from rdkit import Chem
import math
import os
import argparse

parser = argparse.ArgumentParser(description = 'Organize your sdf file')
parser.add_argument('sdf', help='your sdf file') # 필요한 인수를 추가
parser.add_argument('bunch', help='how many molecules per a file') # 필요한 인수를 추가
args = parser.parse_args()

sdf = args.sdf
bunch = args.bunch
bunch = int(bunch)

supplier = Chem.SDMolSupplier(sdf)
fn = os.path.splitext(sdf)[0]

for i in range(0, math.ceil(len(supplier)/bunch)):
    tf = [mol for mol in supplier[i*bunch:((i+1)*bunch)] if mol is not None]
    print(i*bunch, ((i+1)*bunch))
    writer = Chem.SDWriter(f'{fn}_{i}.sdf')
    for mol in tf:
        writer.write(mol)
    writer.close()

print("File_N.sdf are generated")
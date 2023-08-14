#! /home/siu/anaconda3/envs/rdkit/bin/python

import rdkit
from rdkit.Chem import rdmolfiles
import gzip
import argparse
import os

parser = argparse.ArgumentParser(description = 'count num of cpds in mae or maegz')
parser.add_argument('file', metavar='file', help='input .maegz/.mae/.sdf files')
args = parser.parse_args()

name = args.file
print(name)

if '.mae' in name:
    suppl = rdmolfiles.MaeMolSupplier(name)
    
    GetNumAtoms = []
    for mol in suppl:
        if mol is not None: GetNumAtoms.append(mol.GetNumAtoms())
    print(len(GetNumAtoms))

elif '.maegz' in name:
    import gzip
    suppl = rdmolfiles.MaeMolSupplier(gzip.open(name))

    GetNumAtoms = []
    for mol in suppl:
        if mol is not None: GetNumAtoms.append(mol.GetNumAtoms())
    print(len(GetNumAtoms))

elif '.sdf' in name:
    suppl = rdmolfiles.ForwardSDMolSupplier(name)
    
    GetNumAtoms = []
    for mol in suppl:
        if mol is not None: GetNumAtoms.append(mol.GetNumAtoms())
    print(len(GetNumAtoms))
    
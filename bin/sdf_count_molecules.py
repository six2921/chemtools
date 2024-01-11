#! /home/siu/anaconda/envs/chem/bin/python

from rdkit.Chem import rdmolfiles
import argparse

parser = argparse.ArgumentParser(description = 'count num of cpds in mae or maegz')
parser.add_argument('file', metavar='file', help='input .sdf files')
args = parser.parse_args()

name = args.file
print(name)

if '.sdf' in name:
    suppl = rdmolfiles.ForwardSDMolSupplier(name)
    
    GetNumAtoms = []
    for mol in suppl:
        if mol is not None: GetNumAtoms.append(mol.GetNumAtoms())
    print(len(GetNumAtoms))
else:
    print('not sdf file')
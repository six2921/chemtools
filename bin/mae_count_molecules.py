from rdkit.Chem import rdmolfiles
import argparse

parser = argparse.ArgumentParser(description = 'count num of cpds in mae or maegz')
parser.add_argument('file', metavar='file', help='input .mae files')
args = parser.parse_args()

name = args.file
print(name)

if '.mae' in name:
    suppl = rdmolfiles.MaeMolSupplier(name)
    
    GetNumAtoms = []
    for mol in suppl:
        if mol is not None: GetNumAtoms.append(mol.GetNumAtoms())
    print(len(GetNumAtoms))
else:
    print('not mae file')
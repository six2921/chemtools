import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel

# argparse 설정 및 인자 파싱
parser = argparse.ArgumentParser(description="Convert SMILES to MOL2 file with optional 3D structure optimization.")
parser.add_argument("smiles", type=str, help="SMILES string to convert")
parser.add_argument("--name", type=str, help="Output file name without extension (default is SMILES string).")
parser.add_argument("--opt", action="store_true", help="Generate 3D structure and optimize using UFF.")
args = parser.parse_args()

def smiles_to_mol2(smiles, file_name=None, optimize=False):
    # RDKit을 사용해 SMILES를 Mol 객체로 변환
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES string")
        return
    
    if optimize:
        # 3D 구조 생성 및 최적화
        mol = Chem.AddHs(mol)  # 수소 추가
        if AllChem.EmbedMolecule(mol) == -1:
            print("Failed to generate 3D structure")
            return
        # UFF로 최적화
        if mol.GetNumConformers():
            AllChem.UFFOptimizeMolecule(mol)
            print("3D structure generated and optimized using UFF.")
        else:
            print("Conformer generation failed; skipping optimization.")
    else:
        print("Generating 2D structure without optimization.")
    
    # RDKit Mol 객체를 MOL 형식으로 변환
    mol_block = Chem.MolToMolBlock(mol)
    
    # Open Babel로 MOL2 변환
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("mol", "mol2")
    
    ob_mol = openbabel.OBMol()
    ob_conversion.ReadString(ob_mol, mol_block)
    
    # 파일 이름 설정
    if file_name:
        mol2_file = f"{file_name}.mol2"
    else:
        mol2_file = f"{smiles}.mol2"
    
    ob_conversion.WriteFile(ob_mol, mol2_file)
    print(f"MOL2 file saved as {mol2_file}")

# 함수 호출
smiles_to_mol2(args.smiles, args.name, args.opt)


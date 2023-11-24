#! /home/siu/anaconda3/envs/rdkit/bin/python

import argparse, os, glob, shutil
from rdkit import Chem
from rdkit.Chem import AllChem
from natsort import natsorted
import pickle
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# argparse 객체 생성
parser = argparse.ArgumentParser(description="Read sdf and save as pkls")
parser.add_argument("filename", help="SDF file")
args = parser.parse_args()

filename = args.filename

# ------------- SPLIT FILE -------------
limit=10000  # 파일 쪼개는 기준

# 파일 읽기
with open(filename, "r") as f:
    text = f.read()

# "$$$$" 문자열 갯수 확인
num_delimiters = text.count("$$$$")
print(f"Number of Cpds: {num_delimiters}")

# limit개 단위로 파일 쪼개기
if num_delimiters > limit:
    # "$$$$" 문자열을 기준으로 파일 쪼개기
    parts = text.split("$$$$\n")
    for i in range(0, num_delimiters, limit):
        # 파일 이름 생성
        new_filename = f"{os.path.splitext(filename)[0]}_{i//limit+1}.sdf"

        # 새로운 파일에 "$$$$" 문자열과 그 아래 텍스트 쓰기
        with open(new_filename, "w") as f:
            for j, part in enumerate(parts[i:i+limit]):
                if j > 0:
                    f.write("$$$$\n")
                f.write(part)
            f.write("$$$$\n")
else:
    new_filename = f"{os.path.splitext(filename)[0]}_1.sdf"
    shutil.copy(filename, new_filename)  # 그냥 _1 붙여서 복사

# 파일 리스트 읽기
pattern = f"{os.path.splitext(filename)[0]}_*.sdf"
files = glob.glob(pattern)
files = natsorted(files)  # 파일 이름 순서대로 정렬
print(f'Number of splited files: {len(files)} | ({limit} cpds per file)')

# ------------- CHECK PROPS -------------
# 첫번째 분자의 프로퍼티와 값 출력
suppl = Chem.SDMolSupplier(files[0])
mol = suppl[0]
for prop in mol.GetPropNames():
    print(f"{prop}: {mol.GetProp(prop)}")

# 사용자 입력 받기
title = input("Enter the property name to copy to the molecule name: ")

if title not in mol.GetPropNames():
    print(f"Error: Property '{title}' not found in the file.")
    exit()

# ------------- READ AND MERGE -------------
# 파일 읽기 (병렬)
def read_file(file):
    suppl = Chem.SDMolSupplier(file)
    mols = [mol for mol in suppl if mol is not None and mol.HasProp(title)]
    names = [mol.GetProp(title) for mol in mols]
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
    return mols, names, smiles, fps

with Pool(processes=cpu_count()) as p:
    results = list(tqdm(p.imap(read_file, files), total=len(files)))

# 딕셔너리에 데이터 저장
data = {
    'mols': [],
    'names': [],
    'smiles': [],
    'fps': []
}

for result in results:
    data['mols'].append(result[0])
    data['names'].append(result[1])
    data['smiles'].append(result[2])
    data['fps'].append(result[3])

# 딕서녀리에 저장된 리스트를 평평하게
data = {
    'mols': [mol for file_mols in data['mols'] for mol in file_mols],
    'names': [name for file_names in data['names'] for name in file_names],
    'smiles': [smile for file_smiles in data['smiles'] for smile in file_smiles],
    'fps': [fp for file_fps in data['fps'] for fp in file_fps]
}

print("Number of mols, names, smiles, fps:")
print(len(data['mols']), len(data['names']), len(data['smiles']), len(data['fps']))
print("Unconverted molecules:", num_delimiters-len(data['mols']))

# 피클로 저장
print("Saving...")
fn = os.path.splitext(os.path.basename(filename))[0]

for key, value in data.items():
    with open(f'{fn}_{key}.pkl', 'wb') as f:
        pickle.dump(value, f)

print("Done!")

# 파일 삭제
for file in files:
    os.remove(file)

#! /home/siu/anaconda/bin/python

import argparse, os, time, glob, shutil
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
from natsort import natsorted
from multiprocessing import Pool, cpu_count
import sys

# argparse 설정
# 변수는 add_X 형태로 설정해야 함
parser = argparse.ArgumentParser(description="Add properties to an SDF file")
parser.add_argument("filename", help="SDF file")
parser.add_argument("--add_smiles", action="store_true", help="Add SMILES as a property")
parser.add_argument("--add_fp", action="store_true", help="Add FP as a property")
parser.add_argument("--add_mw", action="store_true", help="Add MW as a property")
args = parser.parse_args()

filename = args.filename
add_smiles = args.add_smiles
add_mw = args.add_mw

# 옵션 변수가 주어져 있지 않으면 에러 메시지 출력하고 종료
arg_list = [arg for arg in vars(args)]  # arg 변수 리스트

count = 0
for arg in arg_list:
    if getattr(args, arg) == True:
        count += 1
if count < 1:
    print("Error: No arguments given.")
    sys.exit()

# ------------- FUNCTIONS -------------
# 함수 이름은 반드시 argparse 변수 이름과 같아야 함
def add_smiles(filename):
    suppl = Chem.SDMolSupplier(filename)
    res = [Chem.MolToSmiles(x, isomericSmiles=True) for x in suppl]
    print(f"Generated SMILES: {len(res)}")
    return res

def add_mw(filename):
    suppl = Chem.SDMolSupplier(filename)
    res = [int(Chem.rdMolDescriptors.CalcExactMolWt(x)) for x in suppl]
    print(f"Generated MW: {len(res)}")
    return res

def add_fp(filename):
    suppl = Chem.SDMolSupplier(filename)
    res = [AllChem.RDKFingerprint(x) for x in suppl]
    print(f"Generated FP: {len(res)}")
    return res

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


# ------------- PARALLEL CALC -------------
dict_results = {}  # 결과를 저장할 딕셔너리

for arg in arg_list:  # action이 store_true인 arg 변수만 실행
    if getattr(args, arg) == True:
        with Pool(processes=cpu_count()//2) as p:
            start_time = time.time()
            results = p.map(eval(arg), files)
            
            # 결과 저장
            list_results = []
            for result in results:
                list_results.extend(result)
            dict_results[arg] = list_results
            
            # 결과 출력
            print(f"Generated {arg}: {len(list_results)}")

            end_time = time.time()
            print(f"Time: {end_time - start_time:.2f} seconds")

# suppl에 프로퍼티 추가 (병렬 계산 아님)
start_time = time.time()

suppl = Chem.SDMolSupplier(filename)
new_suppl = []
for i, mol in enumerate(suppl):
    for k in dict_results.keys():
        prop_name = k.replace("add_", "rdkit_")
        mol.SetProp(prop_name, str(dict_results[k][i]))
    new_suppl.append(mol)
    #print(f"Reading and adding new props", {i}, end="\r")

    suppl_len = len(suppl)
    print(f'Reading and adding new props: {i+1} / {suppl_len}, ({time.time()-start_time:.2f} sec)', end='\r')


# ------------- WRITE -------------
print('Saving as single file...')
start_time = time.time()

# 새로운 파일에 저장
writer = Chem.SDWriter("output.sdf")

suppl_len = len(new_suppl)
for i, mol in enumerate(new_suppl):
    writer.write(mol)
    print(f'Wrting mols as SDF: {i+1} / {suppl_len}, ({time.time()-start_time:.2f} sec)', end='\r')
writer.close()

# CSV 파일로 저장
with open("output.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    props = new_suppl[0].GetPropNames()  # 프로퍼티 리스트 가져오기
    props = [prop for prop in props] # 리스트로 변환
    writer.writerow(["Name"] + props)  # 헤더 쓰기
    for i, mol in enumerate(new_suppl):
        print(f"Writing csv", {i}, end="\r")
        name = mol.GetProp("_Name")
        row = [name] + [mol.GetProp(prop) for prop in props]  # 데이터 쓰기
        writer.writerow(row)

# SMILES 파일로 저장
with open("output.smi", "w") as smifile:
    writer = Chem.SmilesWriter(smifile, includeHeader=False)
    for i, mol in enumerate(new_suppl):
        print(f"Writing csv", {i}, end="\r")
        writer.write(mol)
    writer.close()

# 시간 측정
end_time = time.time()
print(f"Writing Time: {end_time - start_time:.2f} seconds")

# 생성된 파일 이름들 프린트
print('Generated files: output.sdf, output.csv, output.smi')

# 파일 삭제
for file in files:
    os.remove(file)
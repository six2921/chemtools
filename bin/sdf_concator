from rdkit import Chem
import argparse
import glob

parser = argparse.ArgumentParser(description="Concatenate SDF files from a directory into one SDF file.")

def read_sdf_files(file_list):
    sdf_data = []
    for file in file_list:
        suppl = Chem.SDMolSupplier(file)
        for mol in suppl:
            if mol is not None:
                sdf_data.append(mol)
            else:
                print(f"Error reading molecule from file: {file}")
    return sdf_data

def save_to_sdf(mol_list, output_file):
    writer = Chem.SDWriter(output_file)
    for mol in mol_list:
        writer.write(mol)
    writer.close()

# SDF 파일 목록 가져오기
sdf_files = glob.glob("*.sdf")

# 결과 파일 스킵
sdf_files = [f for f in sdf_files if f != "sdf_concator_output.sdf"]

# SDF 파일 읽기
mol_list = read_sdf_files(sdf_files)

# 결과를 새로운 SDF 파일로 저장
save_to_sdf(mol_list, "sdf_concator_output.sdf")

print("sdf_concator_output.sdf is generated")


import argparse
from rdkit import Chem
import pandas as pd
import os

def read_sdf_file(file):
    sdf_data = []
    suppl = Chem.SDMolSupplier(file)
    for mol in suppl:
        if mol is not None:
            sdf_data.append(mol)
        else:
            print(f"Error reading molecule from file: {file}")
    return sdf_data

def add_duplicate_column(mol_list):
    smiles_list = [Chem.MolToSmiles(mol) for mol in mol_list]
    df = pd.DataFrame(smiles_list, columns=["SMILES"])
    df["duplicate_id"] = 0
    duplicate_counter = 1
    for smiles in df["SMILES"].unique():
        duplicates = df[df["SMILES"] == smiles]
        if len(duplicates) > 1:
            df.loc[df["SMILES"] == smiles, "duplicate_id"] = duplicate_counter
            duplicate_counter += 1
    df["duplicate_id"] = df["duplicate_id"].astype(str)
    df.loc[df["duplicate_id"] == "0", "duplicate_id"] = ""
    return df

def save_to_sdf(df, mol_list, output_file):
    writer = Chem.SDWriter(output_file)
    for idx, row in df.iterrows():
        mol = mol_list[idx]
        if row["duplicate_id"]:
            mol.SetProp("duplicate_id", row["duplicate_id"])
        writer.write(mol)
    writer.close()

def save_unique_to_sdf(df, mol_list, output_file):
    unique_indices = df.drop_duplicates(subset="SMILES", keep="first").index
    writer = Chem.SDWriter(output_file)
    for idx in unique_indices:
        mol = mol_list[idx]
        writer.write(mol)
    writer.close()
    return len(df) - len(unique_indices)

def save_removed_to_sdf(df, mol_list, output_file):
    removed_indices = df[df.duplicated(subset="SMILES", keep="first")].index
    writer = Chem.SDWriter(output_file)
    for idx in removed_indices:
        mol = mol_list[idx]
        writer.write(mol)
    writer.close()

def get_property_values(mol, property_names):
    """Get properties from the molecule in a case-insensitive manner."""
    properties = {}
    lower_props = {prop.lower(): prop for prop in mol.GetPropNames()}
    for prop in property_names:
        prop_lower = prop.lower()
        if prop_lower in lower_props:
            properties[prop] = mol.GetProp(lower_props[prop_lower])
    return properties

def main(input_file):
    # 입력 파일 이름과 확장자를 분리
    file_root, file_ext = os.path.splitext(input_file)
    
    # SDF 파일 읽기
    mol_list = read_sdf_file(input_file)

    # 중복 구조에 대해 'duplicate_id' 컬럼 추가
    df = add_duplicate_column(mol_list)

    # 결과 파일 이름 생성
    output_file = f"{file_root}_marked{file_ext}"
    unique_output_file = f"{file_root}_unique{file_ext}"
    removed_output_file = f"{file_root}_removed{file_ext}"

    # 결과를 새로운 SDF 파일로 저장
    save_to_sdf(df, mol_list, output_file)
    print(f"{output_file} is generated")

    # 중복을 제거한 결과를 또 다른 SDF 파일로 저장하고, 중복으로 삭제된 항목 수 반환
    num_removed = save_unique_to_sdf(df, mol_list, unique_output_file)
    print(f"{unique_output_file} is generated")
    print(f"Number of duplicate molecules removed: {num_removed}")

    # 삭제된 항목들을 포함한 SDF 파일로 저장
    save_removed_to_sdf(df, mol_list, removed_output_file)
    print(f"{removed_output_file} is generated")

    # 삭제된 항목들의 프로퍼티 출력
    property_names = ["id", "name", "title"]
    duplicate_indices = df[df.duplicated(subset="SMILES", keep="first")].index
    for idx in duplicate_indices:
        mol = mol_list[idx]
        properties = get_property_values(mol, property_names)
        if properties:
            print(f"Removed molecule properties: {properties}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate SDF files into one SDF file with duplicate checks.")
    parser.add_argument("input_file", type=str, help="Input SDF file to process.")
    args = parser.parse_args()

    main(args.input_file)

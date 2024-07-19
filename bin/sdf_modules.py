from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from pathos.multiprocessing import ProcessingPool
from tqdm import tqdm
import os
import pandas as pd

def sdf_calc_props(file):
    mols = [mol for mol in Chem.SDMolSupplier(file, removeHs=False) if mol is not None]
    with Chem.SDWriter(file) as writer:
        for mol in mols:
            mol.SetProp("rdSMILES", Chem.MolToSmiles(mol))
            mol.SetProp("rdMW", "{:.2f}".format(Descriptors.MolWt(mol)))
            mol.SetProp("rdLogP", "{:.2f}".format(Descriptors.MolLogP(mol)))
            mol.SetProp("rdRotB", str(Descriptors.NumRotatableBonds(mol)))
            mol.SetProp("rdTPSA", "{:.2f}".format(Descriptors.TPSA(mol)))
            writer.write(mol)
    return file

def split_sdf_file(input_file, chunk_size=10000):
    """Split an SDF file into smaller chunks."""
    file_base_name, file_ext = os.path.splitext(input_file)
    generated_files = []
    with open(input_file, 'r') as file:
        content = file.read().strip().split('$$$$\n')
    
    for i in range(0, len(content), chunk_size):
        chunk = content[i:i + chunk_size]
        if not chunk:
            continue
        output_file = f"{file_base_name}_{i//chunk_size+1}{file_ext}"
        with open(output_file, 'w') as out_file:
            out_file.write('$$$$\n'.join(chunk) + '$$$$\n')
        generated_files.append(output_file)
    
    return generated_files

def merge_sdf_files(file_list):
    """Merge a list of SDF files into a single file."""
    output_file = file_list[0].replace("_1", "_merged")
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            with open(filename, 'r') as infile:
                outfile.write(infile.read())
    return output_file

def sdf_change_title(file, prop_name):
    suppl = Chem.SDMolSupplier(file)
    new_suppl = []
    for mol in suppl:
        if mol is None:
            continue
        try:
            prop_value = mol.GetProp(prop_name)
            mol.SetProp("_Name", prop_value)
            new_suppl.append(mol)
        except KeyError:
            continue

    with Chem.SDWriter(file) as writer:
        for mol in new_suppl:
            writer.write(mol)

    return file

def sdf_to_csv(file):
    mols = [mol for mol in Chem.SDMolSupplier(file, removeHs=False) if mol is not None]
    records = []
    for mol in mols:
        record = {}
        for prop in mol.GetPropNames():
            record[prop] = mol.GetProp(prop)
        record["SMILES"] = Chem.MolToSmiles(mol)
        records.append(record)
    
    df = pd.DataFrame(records)
    csv_file = file.replace('.sdf', '.csv')
    df.to_csv(csv_file, index=False)
    return csv_file

def concat_csv_files(file_list):
    """Concatenate a list of CSV files into a single file."""
    df_list = [pd.read_csv(file) for file in file_list]
    concatenated_df = pd.concat(df_list, ignore_index=True)
    output_file = file_list[0].replace("_1.csv", "_merged.csv")
    concatenated_df.to_csv(output_file, index=False)
    return output_file
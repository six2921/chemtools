import os
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem import PandasTools
from tqdm import tqdm

# ---------------------
# SDF Property Calculation Function
# ---------------------

def sdf_calc_props(file):
    """
    Calculate and add chemical properties to each molecule in an SDF file,
    then write the molecules back to the same file.
    """
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

# ---------------------
# SDF Splitting and Merging Functions
# ---------------------

def split_sdf_file_auto(input_file):
    """
    Automatically split an SDF file into chunks based on the total number of compounds
    and the system's CPU count.
    
    Let N be the total number of compounds and C be the CPU count.
      - If N/C <= 10, use number of chunks = C.
      - If 11 <= N/C <= 1000, use number of chunks = max(1, C//2).
      - If 1001 <= N/C <= 10000, use number of chunks = C.
      - If N/C > 10000, then each chunk will contain 10,000 compounds.
    The chunk size is calculated as ceil(N / (number of chunks)).
    Returns a list of generated SDF file names.
    """
    cpu_count = os.cpu_count() or 1

    with open(input_file, 'r') as file:
        content = file.read().strip().split('$$$$\n')
    compounds = [x for x in content if x.strip() != ""]
    total = len(compounds)
    ratio = total / cpu_count

    if ratio <= 10:
        num_chunks = cpu_count
        chunk_size = math.ceil(total / num_chunks)
    elif ratio <= 1000:
        num_chunks = max(1, cpu_count // 2)
        chunk_size = math.ceil(total / num_chunks)
    elif ratio <= 10000:
        num_chunks = cpu_count
        chunk_size = math.ceil(total / num_chunks)
    else:
        chunk_size = 10000

    file_base_name, file_ext = os.path.splitext(input_file)
    generated_files = []
    for i in range(0, total, chunk_size):
        chunk = compounds[i:i + chunk_size]
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

# ---------------------
# SDF Property Listing and Selection Functions
# ---------------------

def list_properties(input_file: str) -> list:
    """
    Print property names of the first valid molecule in a formatted table.

    Returns
    -------
    list
        Property names in order of appearance.
    """
    supplier = Chem.SDMolSupplier(input_file)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        print("No valid molecule found in the SDF file.")
        return []

    props = list(mol.GetPropNames())

    print(f"{'idx':<5} {'property':<25} {'type':<10} {'value(first)':<25}")
    print("-" * 70)
    for idx, prop in enumerate(props):
        value = mol.GetProp(prop)
        print(f"{idx:<5} {prop[:25]:<25} {'object':<10} {value[:25]:<25}")

    return props


def select_properties(input_file: str, *prop_types) -> tuple:
    """
    Ask the user to select property indices shown by `list_properties`.

    Returns
    -------
    tuple
        Selected property names (strings).
    """
    prop_list = list_properties(input_file)
    if not prop_list:
        exit(1)

    chosen = []
    for prop_type in prop_types:
        prompt = f"Enter the number for the {prop_type} property: "
        while True:
            try:
                idx = int(input(prompt).strip())
                if idx < 0 or idx >= len(prop_list):
                    raise ValueError
                chosen.append(prop_list[idx])
                break
            except ValueError:
                print("Invalid input. Please enter a valid number.")
    return tuple(chosen)


"""
File: csv_add_prop.py
Usage: python csv_add_prop.py input.csv [--mw] [--coo] [--tpsa] [--pains]

This script adds chemical properties to a CSV file based on a SMILES column.
It calculates properties for the SMILES values in the selected column.
Available flags:
    --mw     : Add molecular weight (MW) property.
    --coo    : Add COO group count property.
    --tpsa   : Add TPSA property.
    --pains  : Add PAINS filter property.
The SMILES column is chosen interactively.
"""

import argparse
import os
import sys
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import FilterCatalog
from rdkit import RDLogger
from helpers_csv import select_columns, is_valid_smiles

# Suppress RDKit warnings
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Setup PAINS filter catalog for the 'pains' property.
pains_params = FilterCatalog.FilterCatalogParams()
pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
pains_catalog = FilterCatalog.FilterCatalog(pains_params)

def calculate_mw(smiles: str):
    """Return molecular weight for the given SMILES, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.MolWt(mol) if mol else None

def count_coo(smiles: str):
    """Return the count of COO groups in the SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 0
    pattern_cooh = Chem.MolFromSmarts("C(=O)O")
    pattern_carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    matches_cooh = mol.GetSubstructMatches(pattern_cooh)
    matches_carboxylate = mol.GetSubstructMatches(pattern_carboxylate)
    matches = set(matches_cooh + matches_carboxylate)
    return len(matches)

def calculate_tpsa(smiles: str):
    """Return TPSA for the given SMILES, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcTPSA(mol) if mol else None

def check_pains(smiles: str):
    """Return a tuple (PAINS status, description) for the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "none", None
    entry = pains_catalog.GetFirstMatch(mol)
    if entry is not None:
        return "has", entry.GetDescription()
    return "none", None

# Parse arguments
parser = argparse.ArgumentParser(
    description="Add chemical properties to a CSV file based on a SMILES column.",
    epilog="Example: python csv_add_prop.py data.csv --mw --pains"
)
parser.add_argument("input_csv", help=".csv file (requires a column for classification)")
parser.add_argument("--mw", action="store_true", help="Add molecular weight (MW) property")
parser.add_argument("--coo", action="store_true", help="Add COO group count property")
parser.add_argument("--tpsa", action="store_true", help="Add TPSA property")
parser.add_argument("--pains", action="store_true", help="Add PAINS filter property")
args = parser.parse_args()

if not (args.mw or args.coo or args.tpsa or args.pains):
    parser.error("No property flag provided. Choose at least one of: --mw, --coo, --tpsa, --pains.")

input_csv = args.input_csv

# Read the CSV file
df = pd.read_csv(input_csv, low_memory=False)

# Use the helper function to select the SMILES column interactively.
# Since we only need the SMILES column, unpack as a one-element tuple.
(smiles_col,) = select_columns(df, "smiles")

# Ensure the selected SMILES column is non-null and convert values to strings.
df[smiles_col] = df[smiles_col].fillna('').astype(str)

# Validate the SMILES values.
print("\nValidating SMILES...")
invalid_smiles = [s for s in df[smiles_col] if not is_valid_smiles(s)]
if invalid_smiles:
    print("Error: The following invalid SMILES were found. Please fix your file:")
    for s in invalid_smiles:
        print(s)
    sys.exit(1)
else:
    print("All SMILES are valid.")

# Use ProcessPoolExecutor to calculate properties in parallel.
def parallel_map(func, data, desc):
    with ProcessPoolExecutor() as executor:
        results = list(tqdm(executor.map(func, data), total=len(data), desc=desc))
    return results

if args.mw:
    print("Calculating molecular weight (MW)...")
    mw_results = parallel_map(calculate_mw, df[smiles_col].tolist(), "Calculating MW")
    df['MW'] = mw_results

if args.coo:
    print("Counting COO groups...")
    coo_results = parallel_map(count_coo, df[smiles_col].tolist(), "Counting COO groups")
    df['COO_count'] = coo_results

if args.tpsa:
    print("Calculating TPSA...")
    tpsa_results = parallel_map(calculate_tpsa, df[smiles_col].tolist(), "Calculating TPSA")
    df['TPSA'] = tpsa_results

if args.pains:
    print("Checking PAINS filter...")
    pains_results = parallel_map(check_pains, df[smiles_col].tolist(), "Checking PAINS")
    df['PAINS'] = [r[0] for r in pains_results]
    df['PAINS_desc'] = [r[1] for r in pains_results]

# Save the updated CSV file (append '_with_props' before the extension)
base, ext = os.path.splitext(input_csv)
output_csv = f"{base}_with_props{ext}"
df.to_csv(output_csv, index=False)
print(f"Done. Results saved in {output_csv}")


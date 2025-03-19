import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import FilterCatalog
import argparse
import os
from rdkit import RDLogger
from tqdm import tqdm
from pathos.multiprocessing import ProcessingPool as Pool

# PAINS 필터 카탈로그 생성
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

# RDKit 경고 억제
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Carboxyl acid 서브구조 SMARTS
COOH_SMARTS = [
    'C(=O)O',
    'C(=O)[O-]'
]

def check_pains(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "none", None
    entry = catalog.GetFirstMatch(mol)
    if entry is not None:
        return 'has', entry.GetDescription()
    return 'none', None

def check_cooh(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    for smarts in COOH_SMARTS:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            return 'COOH'
    return 'No COOH'

def calculate_mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)

def apply_function_with_progress(args):
    index, smiles, func = args
    return index, func(smiles)

def apply_check_pains(args):
    index, smiles = args
    return index, check_pains(smiles)

def apply_check_cooh(args):
    index, smiles = args
    return index, check_cooh(smiles)

def apply_calculate_mw(args):
    index, smiles = args
    return index, calculate_mw(smiles)

def process_column(df, column_name, func, smiles_col, multi_output=False):
    tasks = [(i, smiles) for i, smiles in enumerate(df[smiles_col])]
    with Pool() as pool:
        if func == check_pains:
            results = list(tqdm(pool.imap(apply_check_pains, tasks), total=len(df), desc=f"Processing {column_name}"))
        elif func == check_cooh:
            results = list(tqdm(pool.imap(apply_check_cooh, tasks), total=len(df), desc=f"Processing {column_name}"))
        elif func == calculate_mw:
            results = list(tqdm(pool.imap(apply_calculate_mw, tasks), total=len(df), desc=f"Processing {column_name}"))

    if multi_output:
        return [result for index, (result, detail) in sorted(results)], [detail for index, (result, detail) in sorted(results)]
    return [result for index, result in sorted(results)]

parser = argparse.ArgumentParser(description='Add properties to CSV based on SMILES.')
parser.add_argument('input_file', type=str, help='Input CSV file with SMILES column')
parser.add_argument('--smiles', type=str, default='smiles', help='Name of the column containing SMILES strings')
parser.add_argument('--pains', action='store_true', help='Check for PAINS')
parser.add_argument('--cooh', action='store_true', help='Check for carboxyl acid (COOH)')
parser.add_argument('--mw', action='store_true', help='Calculate molecular weight')

args = parser.parse_args()

input_file = args.input_file
output_file = os.path.splitext(input_file)[0] + '-props.csv'
smiles_col = args.smiles

print("# Reading CSV")
df = pd.read_csv(input_file, low_memory=False)
print("# Done")

invalid_smiles_count = 0

if args.pains:
    pains_results, pains_details = process_column(df, 'PAINS', check_pains, smiles_col, multi_output=True)
    df['PAINS'] = pains_results
    df['PAINS_Type'] = pains_details
    invalid_smiles_count += df['PAINS'].isna().sum()

if args.cooh:
    df['COOH'] = process_column(df, 'COOH', check_cooh, smiles_col)
    invalid_smiles_count += df['COOH'].isna().sum()

if args.mw:
    df['MW'] = process_column(df, 'MW', calculate_mw, smiles_col)
    invalid_smiles_count += df['MW'].isna().sum()

df.to_csv(output_file, index=False)

# 카테고리별 개수 출력
if args.pains:
    pains_counts = df['PAINS'].value_counts(dropna=True)
    print("PAINS Categories:")
    print(pains_counts)

if args.cooh:
    cooh_counts = df['COOH'].value_counts(dropna=True)
    print("COOH Categories:")
    print(cooh_counts)

print(f"Number of invalid SMILES: {invalid_smiles_count}")

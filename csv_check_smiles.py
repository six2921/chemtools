"""
File: csv_check_smiles.py
Usage: python csv_check_smiles.py input.csv

This script verifies the SMILES in the selected column of a CSV file.
If any invalid SMILES are found, they are printed and the script exits.
If all SMILES are valid, the script immediately prompts you to choose which transformation(s)
you want to apply interactively to the SMILES values.
Available options:
   (t) Standardize (convert to canonical SMILES)
   (n) Neutralize (apply neutralization; placeholder implementation)
   (p) Strip salts (keep only the largest fragment)
   (s) Save and finish
Once a transformation is applied, that option will not be offered again.
The resulting CSV is saved as _smiles_checked.csv.
"""

import argparse
import sys
import os
import pandas as pd
from tqdm import tqdm
from helpers_csv import select_columns, is_valid_smiles, sanitize_smiles

# Parse arguments
parser = argparse.ArgumentParser(
    description="Verify SMILES in a CSV file and optionally transform them interactively.",
    epilog="Example: python csv_check_smiles.py data.csv"
)
parser.add_argument("input_csv", help="Path to the input CSV file")
args = parser.parse_args()

input_csv = args.input_csv
df = pd.read_csv(input_csv)

# Use the helper function to interactively select the SMILES column.
# The select_columns function lists column indices, names, and first value, and returns the selected column(s).
(smiles_col,) = select_columns(df, "smiles")

# Check that the selected column exists and fill NaN values.
if smiles_col not in df.columns:
    print(f"Error: Column '{smiles_col}' not found in CSV file.")
    sys.exit(1)
df[smiles_col] = df[smiles_col].fillna('')

# Validate SMILES values.
print("Checking SMILES validity...")
tqdm.pandas(desc="Validating SMILES")
df['valid'] = df[smiles_col].progress_apply(is_valid_smiles)

invalid = df.loc[~df['valid'], smiles_col]
if not invalid.empty:
    print("Invalid SMILES found:")
    for s in invalid:
        print(s)
    sys.exit(1)
else:
    print("All SMILES are valid. Proceeding with transformation...")

df.drop(columns=['valid'], inplace=True)

# Interactively sanitize the SMILES values.
df[smiles_col] = sanitize_smiles(df[smiles_col], smiles_col)

# Save the resulting CSV file with a new name.
output_csv = os.path.splitext(input_csv)[0] + "_smiles_checked.csv"
df.to_csv(output_csv, index=False)
print(f"Sanitized CSV saved as '{output_csv}'.")


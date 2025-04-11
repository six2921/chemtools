"""
File: csv_to_sdf.py
Usage: python csv_to_sdf.py input.csv

This script converts a CSV file to an SDF file using RDKit after interactively sanitizing SMILES.
It performs the following steps:
  1. Reads the CSV file.
  2. Displays the column names (with indices and an example value) and prompts you to select:
       - The molecule name column.
       - The SMILES column.
  3. Validates that all SMILES in the selected column are valid.
     If any invalid SMILES are found, they are printed and the script exits.
  4. Interactively applies transformations (standardize, neutralize, strip salts) to the SMILES.
     The available options are:
       (t) Standardize (convert to canonical SMILES)
       (n) Neutralize (apply neutralization; placeholder implementation)
       (p) Strip salts (retain only the largest fragment)
       (s) Save and finish
     Once an option is applied, it will not be offered again.
     The final transformed SMILES are stored in a new column (<SMILES column> + "_sanit").
  5. Converts the sanitized SMILES into RDKit molecule objects.
  6. Writes an SDF file (with the same base name as the CSV) with all DataFrame columns as properties.
"""

import argparse
import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from helpers_csv import select_columns, is_valid_smiles, sanitize_smiles

# Set up the argument parser
parser = argparse.ArgumentParser(
    description="Convert a CSV file to an SDF file using RDKit after interactively sanitizing SMILES."
)
parser.add_argument("input_csv", help="Path to the input CSV file")
args = parser.parse_args()

# Read the CSV file into a DataFrame
csv_file = args.input_csv
df = pd.read_csv(csv_file)

# Use the helper function to select the molecule name and SMILES columns.
# This function prints the column names with indices and a sample value, then prompts for input.
name_col, smiles_col = select_columns(df, "name", "smiles")

# Ensure the SMILES column is filled and converted to strings.
df[smiles_col] = df[smiles_col].fillna('').astype(str)

# Validate all SMILES values.
print("\nValidating SMILES...")
invalid_smiles = [s for s in df[smiles_col] if not is_valid_smiles(s)]
if invalid_smiles:
    print("Error: The following invalid SMILES were found. Please fix your file:")
    for s in invalid_smiles:
        print(s)
    sys.exit(1)
else:
    print("All SMILES are valid.")

# Interactively sanitize SMILES values.
print("\nProceeding with SMILES sanitization...")
sanitized_smiles = sanitize_smiles(df[smiles_col], smiles_col)

# Save the sanitized SMILES in a new column (original SMILES column name + "_sanit").
sanitized_col = smiles_col + "_sanit"
df[sanitized_col] = sanitized_smiles

# Convert the sanitized SMILES into RDKit molecule objects and store them in the "ROMol" column.
df['ROMol'] = df[sanitized_col].apply(lambda s: Chem.MolFromSmiles(s))

# Generate the output SDF filename by replacing the .csv extension with .sdf.
output_sdf = os.path.splitext(csv_file)[0] + ".sdf"

# Write the DataFrame to an SDF file with all columns as properties.
PandasTools.WriteSDF(df, output_sdf, molColName='ROMol', idName=name_col, properties=list(df.columns), allNumeric=False)

print(f"\nSDF file '{output_sdf}' has been generated.")


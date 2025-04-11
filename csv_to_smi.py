import argparse
import pandas as pd
import os
from helpers_csv import select_columns

# Set up the argument parser
parser = argparse.ArgumentParser(
    description="Extract selected columns from a CSV file and write them to a SMI file."
)
parser.add_argument("csv", help="Path to the input CSV file")
args = parser.parse_args()

# Read the CSV file into a DataFrame
df = pd.read_csv(args.csv)

# Use the helper function to prompt the user for selecting the molecule name and SMILES columns.
name_col, smiles_col = select_columns(df)

# Generate the output file name (replace .csv with .smi)
output_file = os.path.splitext(args.csv)[0] + ".smi"

# Write the SMI file with the format: SMILES <tab> Molecule Name
with open(output_file, 'w', encoding="utf-8") as f:
    for idx, row in df.iterrows():
        f.write(f"{row[smiles_col]}\t{row[name_col]}\n")

print(f"\nThe file '{output_file}' has been generated.")


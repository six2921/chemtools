import argparse
import os
import pandas as pd
from rdkit.Chem import PandasTools

# Set up argument parsing
parser = argparse.ArgumentParser(description='Convert SDF and CSV files and merge them')
parser.add_argument('sdf', metavar='sdf', help='Input .sdf file')
parser.add_argument('csv', metavar='csv', help='Input .csv file')
args = parser.parse_args()

# Load the SDF file into a DataFrame
sdf = args.sdf
csv = args.csv

df1 = PandasTools.LoadSDF(sdf, molColName='Molecule')
print("Columns in the SDF file:")
for col in df1.columns:
    try:
        print(f"{col} | {df1[col][0][:10]}")
    except Exception:
        pass

# Prompt the user for the key column in the SDF file
key1 = 'dZ1%mO'
while key1 not in df1.columns:
    key1 = input('Enter the key column for the first (SDF) file: ')
print()

# Load the CSV file into a DataFrame
df2 = pd.read_csv(csv)
print("Columns in the CSV file:")
for col in df2.columns:
    try:
        print(f"{col} | {df2[col][0][:10]}")
    except Exception:
        pass

# Prompt the user for the key column in the CSV file
key2 = 'dZ1%mO'
while key2 not in df2.columns:
    key2 = input('Enter the key column for the second (CSV) file: ')
print()

# Merge the DataFrames
df = pd.merge(df1, df2, left_on=key1, right_on=key2, how='left', suffixes=('', '_1'))

# Write the merged DataFrame to a new SDF file
output_sdf = os.path.splitext(sdf)[0] + '_merged.sdf'
PandasTools.WriteSDF(df, output_sdf, molColName='Molecule', idName=key2, properties=list(df.columns), allNumeric=False)

print(f"{output_sdf} is generated")

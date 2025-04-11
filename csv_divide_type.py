"""
File: csv_divide_type.py
Usage: python csv_divide_type.py input.csv

This script divides a CSV file into separate CSV files based on a classification column.
The column to use for division is chosen interactively.
For each unique value in the selected column, a new CSV file is created (saved in the same directory)
with the input filename appended by an underscore and the safe version of the unique value.
"""

import argparse
import os
import sys
import pandas as pd
from tqdm import tqdm
from helpers_csv import select_columns

# Set up the argument parser
parser = argparse.ArgumentParser(
    description="Divide a CSV file into multiple files based on a classification column in the CSV. For example, if the CSV file contains a 'type' column with values such as 'human' and 'rat', it will create new files named with suffixes '_human' and '_rat' respectively.",
    epilog="Example: python csv_divide_type.py data.csv"
)
parser.add_argument("input_csv", help=".csv file (the classification column will be selected interactively)")
args = parser.parse_args()

input_csv = args.input_csv
df = pd.read_csv(input_csv, low_memory=False)

# Use the helper function to interactively select the classification column.
# Since we want a column regardless of type, we pass "any" as a placeholder.
(class_col,) = select_columns(df, "type")

# Ensure the selected column is non-null.
df[class_col] = df[class_col].fillna('')

# Get unique values in the selected classification column.
unique_vals = df[class_col].unique()

# Divide the CSV file into separate files for each unique value.
for val in tqdm(unique_vals, desc="Dividing CSV"):
    # Create a subset of the DataFrame where the classification column equals the current unique value.
    subset = df[df[class_col] == val]
    # Replace spaces with underscores in the unique value (if it's a string) for safe file naming.
    safe_val = val.replace(" ", "_") if isinstance(val, str) else str(val)
    # Generate the output filename using the input file's directory and base name.
    output_csv = os.path.join(
        os.path.dirname(input_csv),
        os.path.basename(input_csv).replace(".csv", f"_{safe_val}.csv")
    )
    subset.to_csv(output_csv, index=False)
    print(f"Saved {output_csv}")

print("\nDivision complete. _type.csv files have been generated.")


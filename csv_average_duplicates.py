"""
File: csv_average_duplicates.py
Usage: python csv_average_duplicates.py input.csv --column1 name --column2 inhibition

This script checks for duplicate rows in a CSV file based on a specified text column and calculates
the average of another numeric column for those duplicate groups.
Note: Blank values (or values that cannot be converted to a number) 
      are ignored in the average calculation.

Example:
  Duplicate rows exist as follows:
    name, ..., inhibition, ...
    cpd_01, ..., 80, ...
    cpd_01, ..., 100, ...

  They will be processed as follows:
    name, ..., inhibition, ...
    cpd_01, ..., 90, ...

Usage:
  Provide the CSV file as input. The script will display the column information and then
  prompt you to select:
    - A text column (e.g., name)
    - A numeric column (e.g., inhibition)
  Duplicate groups are then displayed with the following format:
      dups   column1              column2
        4    CHEMBL4163523        4.85, none
  Then, you will be prompted:
    "Average values?
     (y) Yes
     (n) No (I need to check this file)"
  If 'y' is selected, the script computes averages for each duplicate group.
  Next, the script prompts:
    "Retain the other columns? (The value from the first duplicate row is used.)
     (y) Yes
     (n) No"
  In both cases, blank values are ignored in the average calculation.
  The result is saved as _avg.csv.
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from helpers_csv import select_columns, is_valid_smiles

# Parse arguments (only the CSV file is needed from command line)
parser = argparse.ArgumentParser(
    description="Average duplicate rows by averaging values from a specified numeric column based on duplicates in a text column.",
    epilog="""Example:
  Duplicate rows exist as follows:
    name, ..., inhibition, ...
    cpd_01, ..., 80, ...
    cpd_01, ..., 100, ...

  They will be processed as follows:
    name, ..., inhibition, ...
    cpd_01, ..., 90, ...

Usage:
  The script displays the column information and prompts you to select a text column and a numeric column.
  Then, you'll be asked:
    "Average values?
     (y) Yes
     (n) No (I need to check this file)"
  If 'y', averages are computed for each duplicate group.
  Next, you'll be asked:
    "Retain the other columns? (The value from the first duplicate row is used.)
     (y) Yes
     (n) No"
  The final result is saved as _avg.csv.
""",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("input_csv", help="Path to the input CSV file")
args = parser.parse_args()

input_csv = args.input_csv
df = pd.read_csv(input_csv, low_memory=False)

# Use the helper function to select two columns: first for the text (e.g., name) and second for numeric (e.g., inhibition)
# Usage: select_columns(df, "name", "numer")
name_col, value_col = select_columns(df, "name", "numer")

# Ensure the selected numeric column is non-null and convert its values to string (if needed conversion later).
df[value_col] = df[value_col].fillna('').astype(str)

# (Optionally, you could verify the numeric column's first value is convertible to float within select_columns.)

# Create duplicate mask based on text column (name_col)
dup_mask = df.duplicated(subset=[name_col], keep=False)
dup_df = df[dup_mask].copy()

if dup_df.empty:
    print(f"No duplicate rows found based on column '{name_col}'.")
    sys.exit(0)

# For average calculation: Group by name_col converting value_col to float (ignoring blanks or non-convertible values)
grouped_avg = dup_df.groupby(name_col)[value_col].apply(lambda x: pd.to_numeric(x, errors='coerce').dropna().tolist())
grouped_avg = grouped_avg.sort_index()

# For display: Group by name_col preserving original values, converting blanks to "none"
grouped_display = dup_df.groupby(name_col)[value_col].apply(
    lambda x: [str(v).strip() if str(v).strip() != "" else "none" for v in x.tolist()]
)
grouped_display = grouped_display.sort_index()

# Display duplicate groups with the specified formatting.
# Format: idx (width 5), column (width 15), type (width 10), value(first) (width 15)
print("\n{:<5} {:<15} {:<10} {:<15}".format("idx", "column", "type", "value(first)"))
print("-" * 50)
for idx, col in enumerate(df.columns):
    col_type = str(df[col].dtype)
    val_first = str(df[col].iloc[0]) if not df.empty else ""
    print("{:<5} {:<15} {:<10} {:<15}".format(idx, col[:15], col_type[:10], val_first[:15]))

print("\nDuplicate groups (for averaging):")
print("{:<5} {:<20} {:<15}".format("dups", "column1", "column2"))
for key in grouped_display.index:
    display_values = grouped_display.loc[key]
    disp_str = ", ".join(display_values)
    count = len(display_values)
    print(f"{count:<5d} {key:<20} {disp_str:<15}")

# Prompt user whether to average these duplicate groups.
user_input = input("\n  Average values?\n  (y) Yes\n  (n) No (I need to check this file)\n  : ").strip().lower()
if user_input != 'y':
    print("Exiting without processing.")
    sys.exit(0)

# Compute the average for each duplicate group (as float) ignoring blank/invalid values.
avg_results = grouped_avg.apply(lambda vals: np.mean(vals) if len(vals) > 0 else np.nan)

# Ask user for output mode: whether to retain the other columns.
output_mode = input("\n  Retain the other columns?\n  (The value from the first duplicate row is used.)\n  (y) Yes\n  (n) No\n  : ").strip().lower()

# Process duplicate groups based on the selected output mode.
if output_mode == 'y':
    # Replace duplicates in the original DataFrame: for each group, keep the first row and update value_col with the average.
    for key, avg_val in avg_results.items():
        indices = df.index[df[name_col] == key].tolist()
        if indices:
            first_idx = indices[0]
            df.at[first_idx, value_col] = avg_val
            for idx in indices[1:]:
                df.drop(index=idx, inplace=True)
    output_file = os.path.splitext(input_csv)[0] + "_avg.csv"
    df.to_csv(output_file, index=False)
    print(f"\nDuplicates replaced. File saved as '{output_file}'.")
else:
    # Create a summary file with only name_col and the computed average for each duplicate group.
    summary_df = avg_results.reset_index()
    summary_df.columns = [name_col, value_col + "_avg"]
    output_file = os.path.splitext(input_csv)[0] + "_avg.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\nSummary file created. File saved as '{output_file}'.")


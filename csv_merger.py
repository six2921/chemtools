import argparse
import os
import sys
import pandas as pd
from helpers_csv import select_columns

# Set up the argument parser
parser = argparse.ArgumentParser(
    description="Merge two CSV files interactively using column selection.",
    epilog="Example: python csv_merger.py left.csv right.csv"
)
parser.add_argument("file1", help="Left CSV file")
parser.add_argument("file2", help="Right CSV file")
args = parser.parse_args()

# Read the left and right CSV files
df1 = pd.read_csv(args.file1)
df2 = pd.read_csv(args.file2)

# Use the helper function to select the key column for each file.
# Since any type is acceptable, we pass "any" as a placeholder.
(key1,) = select_columns(df1, "key1")
(key2,) = select_columns(df2, "key2")

# Define a mapping for merge method based on single-letter input.
merge_options = {
    "o": "outer",
    "i": "inner",
    "l": "left",
    "r": "right"
}

# Prompt user to select the merge method (input: o, i, l, r)
choice = input("Enter merge method [o: outer, i: inner, l: left, r: right]: ").strip().lower()
while choice not in merge_options:
    print("Invalid input. Please enter one of: o, i, l, r.")
    choice = input("Enter merge method [o: outer, i: inner, l: left, r: right]: ").strip().lower()
how = merge_options[choice]

# Merge the two DataFrames using the selected key columns and merge method.
merged_df = pd.merge(df1, df2, left_on=key1, right_on=key2, how=how, suffixes=('', '_1'))
print("Rows in Merged: ", len(merged_df))

# Save the merged DataFrame to a CSV file.
output_file = "merged_out.csv"
merged_df.to_csv(output_file, index=False)
print(f"Merged CSV saved as '{output_file}'.")


import pandas as pd
import math
import os
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description="Split a CSV file into multiple files by row count.")
parser.add_argument("csv", help="Path to the input CSV file")
args = parser.parse_args()

csv_file = args.csv

# Read the CSV file into a DataFrame.
df = pd.read_csv(csv_file)
total_rows = len(df)
print(f"Total rows in the CSV file: {total_rows}")

# Prompt the user to enter the number of rows per file (bunch size)
bunch_input = input("Enter the number of rows per file: ").strip()
try:
    bunch = int(bunch_input)
except ValueError:
    print("Invalid input. Please enter a valid integer for the number of rows per file.")
    exit(1)

# Calculate the number of files required (round up)
num_files = math.ceil(total_rows / bunch)

# Extract the base filename without extension.
fn = os.path.splitext(csv_file)[0]

for i in range(num_files):
    start_index = i * bunch
    end_index = min((i + 1) * bunch, total_rows)
    chunk_df = df[start_index:end_index]
    output_file = os.path.join(os.path.dirname(csv_file), f"{fn}_{i}.csv")
    chunk_df.to_csv(output_file, index=False)
    print(f"{start_index:>6} - {end_index:>6} -> {output_file}")

print("Done")


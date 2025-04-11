import argparse
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description="Select columns from a CSV file and create a new CSV file.")
parser.add_argument('input_file', type=str, help="Path to the input CSV file")
args = parser.parse_args()

# Read the CSV file
df = pd.read_csv(args.input_file)
print("Columns in the CSV file:")
for idx, col in enumerate(df.columns):
    print(f"{idx}: {col}")

# Prompt the user to enter the column indices to extract (e.g., 1 3 10 14)
columns = input("Enter the column indices to extract (separated by space, e.g., 1 3 10 14): ").split()

# Create a new DataFrame with the selected columns
selected_columns = [df.columns[int(col)] for col in columns]
df_selected = df[selected_columns]

# Generate the output file name by appending '-subset' to the input file name
output_file_name = args.input_file.replace('.csv', '') + '-subset.csv'
df_selected.to_csv(output_file_name, index=False)

print(f"A new CSV file '{output_file_name}' with the selected columns has been created.")


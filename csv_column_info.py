import argparse
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description="Display column information for a CSV file.")
parser.add_argument("input_file", type=str, help="Path to the input CSV file")
args = parser.parse_args()

# Read the CSV file
df = pd.read_csv(args.input_file)

# Print the list of columns with their indices
print("Columns in the CSV file:")
for idx, col in enumerate(df.columns):
    print(f"{idx}: {col}")

# Build a DataFrame to hold column statistics
# The statistics include: Column name, Data type, Min, Max, Unique count, NA count, Duplicate count, and a Sample value (first row)
stat = pd.DataFrame(columns=["Column", "Dtype", "Min", "Max", "Unique", "NA Count", "Duplicate Count", "Sample"])
row_index = 0

for col in df.columns:
    # Convert the column to numeric values if possible; values that cannot be converted become NaN.
    converted = pd.to_numeric(df[col], errors='coerce')
    # Check if the conversion resulted in at least one non-NaN value.
    if converted.notna().sum() > 0:
        col_min = converted.min()
        col_max = converted.max()
    else:
        col_min = None
        col_max = None
    
    uniq = df[col].nunique()
    na_count = df[col].isna().sum()
    dupl = df[col].duplicated().sum()
    sample = df[col].iloc[0] if not df[col].empty else None
    
    stat.loc[row_index] = [col, df[col].dtype, col_min, col_max, uniq, na_count, dupl, sample]
    row_index += 1

print("\nColumn Information:")
print(stat)


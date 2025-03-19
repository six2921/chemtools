import argparse
from sdf_modules import split_sdf_file, sdf_to_csv, concat_csv_files, merge_sdf_files
from pathos.multiprocessing import ProcessingPool
from tqdm import tqdm
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert SDF file to CSV.")
parser.add_argument("input_sdf", help="The SDF file to process.")
args = parser.parse_args()

# Split the SDF file
print("Splitting the SDF file...")
generated_files = split_sdf_file(args.input_sdf, chunk_size=1000)

# Define a function to convert SDF to CSV
def process_file(file):
    print(f"Processing file: {file}")  # Debugging line to ensure correct file paths
    return sdf_to_csv(file)

# Process the split SDF files in parallel
print("Processing split SDF files in parallel...")
with ProcessingPool() as pool:
    list_with_progress = tqdm(generated_files, desc="Processing Files")
    processed_files = pool.map(process_file, list_with_progress)

# Merge the processed CSV files
print("Merging processed CSV files...")
merged_csv_file = concat_csv_files(processed_files)

print(f"Process completed. Merged CSV file: {merged_csv_file}")

# Delete the split SDF files and temporary CSV files after merging
for file in generated_files:
    os.remove(file)
for file in processed_files:
    os.remove(file)
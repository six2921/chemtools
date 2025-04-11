import argparse
from helpers_sdf import split_sdf_file_auto, merge_sdf_files, sdf_calc_props
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Set up argument parsing
parser = argparse.ArgumentParser(description="Split, process, and merge SDF files.")
parser.add_argument("sdf", help="The SDF file to process.")
args = parser.parse_args()

sdf_file = args.sdf

# Automatically split the SDF file using the new function, which determines chunk size based on CPU count
print("Splitting the SDF file automatically based on total compounds and CPU count...")
generated_files = split_sdf_file_auto(sdf_file)

# Define a function to process a single SDF file
def process_file(file):
    print(f"Processing file: {file}")  # Debug: show which file is being processed
    return sdf_calc_props(file)

# Process the split SDF files in parallel using ProcessPoolExecutor
print("Processing split SDF files in parallel...")
with ProcessPoolExecutor() as executor:
    processed_files = list(tqdm(executor.map(process_file, generated_files),
                                  total=len(generated_files),
                                  desc="Processing Files"))

# Merge the processed SDF files into a single SDF file
print("Merging processed SDF files...")
merged_file = merge_sdf_files(processed_files)

print(f"Process completed. Merged file: {merged_file}")

# Optional: Delete the temporary split files after merging
for file in generated_files:
    os.remove(file)


import argparse
from sdf_modules import split_sdf_file, merge_sdf_files, sdf_calc_props
from pathos.multiprocessing import ProcessingPool
from tqdm import tqdm
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description="Split, process, and merge SDF files.")
parser.add_argument("input_sdf", help="The SDF file to process.")
args = parser.parse_args()

# Split the SDF file
print("Splitting the SDF file...")
generated_files = split_sdf_file(args.input_sdf, chunk_size=1000)

# Define a function to process a single SDF file
def process_file(file):
    print(f"Processing file: {file}")  # Debugging line to ensure correct file paths
    return sdf_calc_props(file)

# Process the split SDF files in parallel
print("Processing split SDF files in parallel...")
with ProcessingPool() as pool:
    list_with_progress = tqdm(generated_files, desc="Processing Files")
    processed_files = pool.map(process_file, list_with_progress)

# Merge the processed SDF files
print("Merging processed SDF files...")
merged_file = merge_sdf_files(processed_files)

print(f"Process completed. Merged file: {merged_file}")

# Optional: Delete the split files after merging
for file in generated_files:
    os.remove(file)


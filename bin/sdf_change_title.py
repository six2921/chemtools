#!/home/siu/anaconda/bin/python

import argparse
from sdf_modules import split_sdf_file, merge_sdf_files, sdf_change_title
from pathos.multiprocessing import ProcessingPool
from tqdm import tqdm
import os
from rdkit import Chem

# Set up argument parsing
parser = argparse.ArgumentParser(description="Copy or Replace Title Column of SDF File")
parser.add_argument("input_sdf", help="The SDF file to process.")
args = parser.parse_args()

# Split the SDF file
print("Splitting the SDF file...")
generated_files = split_sdf_file(args.input_sdf, chunk_size=1000)

# Display the properties of the first molecule in the first file
suppl = Chem.SDMolSupplier(generated_files[0])
mol = None
for mol in suppl:
    if mol is not None:
        break

if mol is None:
    raise ValueError("No valid molecule found in the first file.")

print("Properties of the first molecule:")
props = mol.GetPropNames()
for prop in props:
    print(f"{prop}: {mol.GetProp(prop)}")

# Ask user for the property name
while True:
    prop_name = input("Enter the property name to copy to the molecule name: ")
    if prop_name in props:
        break
    else:
        print(f"Property '{prop_name}' not found. Please enter a valid property name.")

# Define a function to change the title of an SDF file
def process_file(file):
    print(f"Processing file: {file}")  # Debugging line to ensure correct file paths
    return sdf_change_title(file, prop_name)

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

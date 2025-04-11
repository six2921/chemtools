"""
File: sdf_change_title.py
Usage: python sdf_change_title.py sdf

This script changes molecule titles in an SDF file based on a user-selected property.
Steps:
  1. It reads the SDF file and displays the _Name (title) of the first two valid molecules.
  2. It asks whether you want to change these title values.
  3. If yes, it lists available properties from the first molecule (using helper_sdf's list_properties function) 
     and prompts you to select a property by number.
  4. The _Name property of every molecule is then updated with the value from the selected property.
  5. The updated molecules are written to a new SDF file (input file name with "_title_changed" appended).
"""

import argparse
import os
import sys
from rdkit import Chem
from helpers_sdf import list_properties   # from our helpers_sdf module

# Set up argument parsing with argument name "sdf"
parser = argparse.ArgumentParser(description="Change molecule titles in an SDF file based on a chosen property.")
parser.add_argument("sdf", help="The SDF file to process.")
args = parser.parse_args()

sdf_file = args.sdf

# Read the SDF file using RDKit's SDMolSupplier and collect the _Name properties of the first two valid molecules.
supplier = Chem.SDMolSupplier(sdf_file)
titles = []
for mol in supplier:
    if mol is not None:
        try:
            title_val = mol.GetProp("_Name")
        except KeyError:
            title_val = "<No _Name>"
        titles.append(title_val)
        if len(titles) == 2:
            break

if not titles:
    print("No valid molecules found in the SDF file.")
    sys.exit(1)


# List available properties from the first valid molecule using list_properties from helpers_sdf.
print("\nAvailable properties from the first molecule:")
list_properties(sdf_file)

print("\nCurrent title values for the first two molecules:")
for i, title in enumerate(titles, start=1):
    print(f"{i}: {title}")

# Retrieve property names from the first valid molecule.
supplier = Chem.SDMolSupplier(sdf_file)
first_mol = None
for mol in supplier:
    if mol is not None:
        first_mol = mol
        break

if first_mol is None:
    print("No valid molecule found in the SDF file.")
    sys.exit(1)

props = list(first_mol.GetPropNames())

# Prompt user to select a property by number.
while True:
    try:
        index = int(input("\nProperty to use for the new title: ").strip())
        if index < 0 or index >= len(props):
            raise ValueError
        new_title_prop = props[index]
        break
    except ValueError:
        print("Invalid input. Please enter a valid property number.")

print(f"Selected property: {new_title_prop}")

# Update the '_Name' property for each molecule based on the selected property.
supplier = Chem.SDMolSupplier(sdf_file)
updated_mols = []
for mol in supplier:
    if mol is None:
        continue
    try:
        new_title = mol.GetProp(new_title_prop)
        mol.SetProp("_Name", new_title)
    except KeyError:
        # If the selected property is missing, leave the title unchanged.
        pass
    updated_mols.append(mol)

# Write updated molecules to a new SDF file.
output_file = os.path.splitext(sdf_file)[0] + "_title_changed.sdf"
writer = Chem.SDWriter(output_file)
for mol in updated_mols:
    writer.write(mol)
writer.close()

print(f"\nTitle updated SDF file saved as '{output_file}'.")


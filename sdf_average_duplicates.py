"""
File: sdf_average_duplicates.py
Usage: python sdf_average_duplicates.py input.sdf

This script checks for duplicate molecules based on SMILES in an SDF file.
If duplicates are found, it displays the duplicated groups and prompts the user:
  (k) Keep the first occurrence
  (a) Average numeric properties and keep first text values
  (n) Skip and leave as is

Final result is saved as:
  - *_rd_first.sdf if first kept
  - *_rd_avg.sdf if averaging applied
"""

import argparse
from rdkit import Chem
import pandas as pd
import os

# -------------------- Parse Input --------------------
parser = argparse.ArgumentParser(description="Check and remove duplicate molecules in an SDF file based on SMILES.")
parser.add_argument("sdf", help="Input SDF file")
args = parser.parse_args()
sdf_file = args.sdf

# -------------------- Load Molecules --------------------
mols = [mol for mol in Chem.SDMolSupplier(sdf_file) if mol is not None]
records = []
for mol in mols:
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    smiles = Chem.MolToSmiles(mol)
    records.append((name, smiles, mol))
df = pd.DataFrame(records, columns=["name", "smiles", "mol"])

# -------------------- Find Duplicates --------------------
df_sorted = df.sort_values("smiles")
df_dups = df_sorted[df_sorted.duplicated("smiles", keep=False)]

if df_dups.empty:
    print("✅ No duplicate SMILES found.")
    exit()

print("\n[ Duplicate molecules based on SMILES ]")
print("name           smiles")
print("-" * 60)
for _, row in df_dups.iterrows():
    print(f"{row['name']:<15} {row['smiles']:<40}")

num_groups = df_dups["smiles"].nunique()
print(f"\n🔎 Found {len(df_dups)} duplicated molecules in {num_groups} structure group(s).")

# -------------------- Show first group's properties --------------------
first_smiles = df_dups["smiles"].iloc[0]
group = df_dups[df_dups["smiles"] == first_smiles]
first_mol = group["mol"].iloc[0]

print("\n🔬 Properties for the first duplicate group (first molecule):")
for prop in first_mol.GetPropNames():
    print(f"  {prop}: {first_mol.GetProp(prop)}")

# -------------------- Ask how to handle duplicates --------------------
print("\nHow would you like to handle duplicates?")
print("  (k) Keep first only")
print("  (a) Average numeric, keep first text")
print("  (n) Skip")
choice = input("Your choice: ").strip().lower()

if choice == "n":
    print("❌ No changes made. Exiting.")
    exit()

# -------------------- Deduplicate --------------------
output_file = ""

if choice == "k":
    df_unique = df.drop_duplicates("smiles", keep="first")
    output_file = os.path.splitext(sdf_file)[0] + "_rd_first.sdf"
    writer = Chem.SDWriter(output_file)
    for mol in df_unique["mol"]:
        writer.write(mol)
    writer.close()

    print(f"\n✅ Duplicates removed by keeping first.")
    print(f"Original molecules: {len(df)}")
    print(f"Final molecules: {len(df_unique)}")
    print(f"Saved as: {output_file}")

elif choice == "a":
    grouped = df.groupby("smiles")
    averaged_mols = []

    for smiles, group in grouped:
        mols_in_group = group["mol"].tolist()
        base = mols_in_group[0]
        new_mol = Chem.Mol(base)
        new_mol.ClearProp("_Name")

        for prop in base.GetPropNames():
            values = []
            for mol in mols_in_group:
                if mol.HasProp(prop):
                    val = mol.GetProp(prop).strip()
                    if val.lower() in ["", "nan", "none"]:
                        continue
                    try:
                        values.append(float(val))
                    except ValueError:
                        values = None
                        break
            if values is None or not values:
                new_mol.SetProp(prop, base.GetProp(prop))
            else:
                avg_val = str(round(sum(values) / len(values), 4))
                new_mol.SetProp(prop, avg_val)

        averaged_mols.append(new_mol)

    output_file = os.path.splitext(sdf_file)[0] + "_rd_avg.sdf"
    writer = Chem.SDWriter(output_file)
    for mol in averaged_mols:
        writer.write(mol)
    writer.close()

    print(f"\n✅ Duplicates removed with averaging.")
    print(f"Original molecules: {len(df)}")
    print(f"Final molecules: {len(averaged_mols)}")
    print("\nℹ️ Numeric values averaged. Text fields kept from first entry.")
    print(f"Saved as: {output_file}")


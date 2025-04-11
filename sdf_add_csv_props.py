"""
File: sdf_add_csv_props.py
Merge selected columns from a CSV file into an SDF file using a key column.

Usage
-----
    python sdf_add_csv_props.py input.sdf input.csv
"""

import argparse
import os
import sys
import pandas as pd
from rdkit.Chem import PandasTools
from helpers_csv import select_columns as select_columns_csv, list_columns as list_columns_csv
from helpers_sdf import select_properties  # 남겨두면 SDF‑측 선택이 필요할 때 활용 가능

# ---------- argument parsing ----------
parser = argparse.ArgumentParser(
    description="Merge properties from a CSV file into an SDF file using a key column."
)
parser.add_argument("sdf", help="Input .sdf file (first)")
parser.add_argument("csv", help="Input .csv file (second)")
args = parser.parse_args()

if not args.sdf.endswith(".sdf") or not args.csv.endswith(".csv"):
    print("❌  Usage:  python sdf_add_csv_props.py input.sdf input.csv")
    sys.exit(1)

sdf_file, csv_file = args.sdf, args.csv

# ---------- load SDF ----------
try:
    sdf_df = PandasTools.LoadSDF(sdf_file, molColName="Molecule", removeHs=False)
except Exception as e:
    print(f"❌  Failed to load SDF: {e}")
    sys.exit(1)

# Ensure _Name column exists
if "_Name" not in sdf_df.columns:
    sdf_df["_Name"] = [
        mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else ""
        for mol in sdf_df["Molecule"]
    ]

print("\n🔎  First 5 titles in SDF (_Name):")
print("\n".join(f"  {i+1}. {v}" for i, v in enumerate(sdf_df['_Name'].head(5))))

# ---------- load CSV ----------
try:
    csv_df = pd.read_csv(csv_file)
except Exception as e:
    print(f"❌  Failed to load CSV: {e}")
    sys.exit(1)

# ---------- choose key column in CSV ----------
print("\nSelect the key column in the CSV that matches the SDF _Name:")
(csv_key_col,) = select_columns_csv(csv_df, "key")

# ---------- matching / duplicate statistics ----------
sdf_keys = sdf_df["_Name"].astype(str)
csv_keys = csv_df[csv_key_col].astype(str)

matched     = sdf_keys.isin(csv_keys).sum()
sdf_dup_cnt = sdf_df.duplicated("_Name").sum()
csv_dup_cnt = csv_df.duplicated(csv_key_col).sum()

print(f"\n📊  SDF rows          : {len(sdf_df)}")
print(f"📊  CSV rows          : {len(csv_df)}")
print(f"📊  Matched titles    : {matched}")
print(f"📊  SDF duplicates    : {sdf_dup_cnt}")
print(f"📊  CSV duplicates    : {csv_dup_cnt}")

# CSV duplicates → abort
if csv_dup_cnt > 0:
    print("\n❌  CSV key column contains duplicates. Please deduplicate the CSV and rerun.")
    sys.exit(1)

# SDF duplicates → warn and confirm
if sdf_dup_cnt > 0:
    print(f"\n⚠️  {sdf_dup_cnt} duplicate title(s) found in the SDF file.")
    print("    The same CSV data will be copied to each duplicate molecule.")
    if input("Continue anyway? (y/n): ").strip().lower() != 'y':
        print("Aborted by user.")
        sys.exit(0)

if matched == 0:
    print("\n❌  No matching keys found between SDF and CSV.")
    sys.exit(1)

# ---------- choose CSV columns to import (multi‑select) ----------
print("\nSelect CSV columns to merge into SDF (you can input several indices).")
list_columns_csv(csv_df)                       # show table once

while True:
    raw = input("prop columns (e.g. 2 4 6): ").strip()
    parts = raw.split()                        # space‑separated only
    try:
        idxs = [int(p) for p in parts]
        if any(i < 0 or i >= len(csv_df.columns) for i in idxs):
            raise ValueError
        import_cols = [csv_df.columns[i] for i in idxs]
        break
    except ValueError:
        print("Invalid input. Please enter valid column numbers separated by spaces.")

# ---------- merge ----------
merged_df = pd.merge(
    sdf_df,
    csv_df[[csv_key_col] + import_cols],
    how="left",
    left_on="_Name",
    right_on=csv_key_col,
    suffixes=('', '_csv')
)

out_sdf = os.path.splitext(sdf_file)[0] + "_merged.sdf"
PandasTools.WriteSDF(
    merged_df,
    out_sdf,
    molColName="Molecule",
    idName="_Name",
    properties=list(merged_df.columns),
    allNumeric=False
)

print(f"\n✅  Merged SDF saved as: {out_sdf}")


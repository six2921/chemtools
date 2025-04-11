"""
File: sdf_extract_by_id.py
Goal
----
Pick compounds from an SDF by typing the **numeric part** of their IDs.

Quick rules
-----------
•  PF‑00153  →  153  
•  Type 153 at the prompt to pick that compound.  
•  You can enter many at once, e.g. 153 158 160
"""

import argparse
import sys
import pandas as pd
from rdkit.Chem import PandasTools

# ── CLI argument ──────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Interactively extract molecules from an SDF by numeric ID."
)
parser.add_argument("sdf", help="Input SDF file")
args = parser.parse_args()
sdf_path = args.sdf

# ── Load SDF ──────────────────────────────────────────────────────
try:
    df = PandasTools.LoadSDF(
        sdf_path,
        molColName="Molecule",
        idName="_Name",          # store title line in '_Name'
        removeHs=False,
    )
except Exception as e:
    sys.exit(f"❌  Failed to read SDF: {e}")

if "_Name" not in df.columns:  # back‑stop for empty titles
    df["_Name"] = [m.GetProp("_Name") if m.HasProp("_Name") else ""
                   for m in df["Molecule"]]

print(f"\nTotal molecules : {len(df)}")

# ── Preview first 5 titles + numeric parts ───────────────────────
print("\n🔎  First 5 titles and numeric IDs")
for i, row in df.head(5).iterrows():
    title = row["_Name"]
    numeric = int(str(title).split("-")[-1]) if "-" in str(title) else "n/a"
    print(f"  {i+1:>2}. {title:<25} → {numeric}")

print(
    "\n•  PF‑00153  →  153\n"
    "•  Type 153 to pick that compound.\n"
    "•  You can enter many at once, e.g. 153 158 160\n"
)

# ── Ask for ID column ─────────────────────────────────────────────
id_col = input("Column that contains the IDs: ").strip()
if id_col not in df.columns:
    sys.exit(f"❌  Column '{id_col}' not found.")

# helper numeric column (PREFIX‑0123 → 123)
df["_numeric_id"] = df[id_col].apply(
    lambda x: int(str(x).split("-")[-1]) if "-" in str(x) else None
)

# ── Interactive selection loop ───────────────────────────────────
subset = pd.DataFrame()
while True:
    codes = input("\nIDs to keep (or 'exit'): ").strip()
    if codes.lower() == "exit":
        break
    try:
        nums = [int(c) for c in codes.split()]
    except ValueError:
        print("⚠️  Numbers only.")
        continue
    picked = df[df["_numeric_id"].isin(nums)]
    subset = pd.concat([subset, picked]).drop_duplicates(subset=[id_col])
    print(f"Selected: {len(subset)}")

# ── Write output ─────────────────────────────────────────────────
if subset.empty:
    print("\nNo molecules selected – nothing written.")
else:
    out_file = f"{sdf_path.rsplit('.',1)[0]}_subset.sdf"
    PandasTools.WriteSDF(
        subset,
        out_file,
        idName=id_col,
        properties=list(df.columns),
    )
    print(f"\n🎉  {len(subset)} molecule(s) saved to {out_file}")


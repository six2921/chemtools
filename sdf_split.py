"""
File: sdf_split.py
Split one SDF file into N smaller SDF files.

Usage
-----
python sdf_split.py input.sdf
  → script asks “How many chunks?”
Output:  <input>_1.sdf, <input>_2.sdf, …
"""

import argparse
import math
import os
from rdkit import Chem

# ---------- args ----------
parser = argparse.ArgumentParser(description="Split a large SDF into N chunks.")
parser.add_argument("sdf", help="input SDF file")
args = parser.parse_args()

in_file = args.sdf

# ---------- read molecules ----------
mols = [m for m in Chem.SDMolSupplier(in_file) if m is not None]
total = len(mols)
if total == 0:
    print("❌  No valid molecules found.")
    exit(1)

# ---------- ask chunk count ----------
while True:
    try:
        n_chunks = int(input(f"Total {total} molecules.  How many files? "))
        if n_chunks <= 0:
            raise ValueError
        break
    except ValueError:
        print("Enter a positive integer.")

chunk_size = math.ceil(total / n_chunks)
base, ext = os.path.splitext(in_file)

print("\n🔨  Splitting …")
for i in range(n_chunks):
    start = i * chunk_size
    end   = min((i + 1) * chunk_size, total)
    out_name = f"{base}_{i+1}{ext}"

    w = Chem.SDWriter(out_name)
    for mol in mols[start:end]:
        w.write(mol)
    w.close()

    print(f"  {start:>6}-{end:<6}  →  {out_name}")

print(f"\n✅  Done.  {n_chunks} files written.")


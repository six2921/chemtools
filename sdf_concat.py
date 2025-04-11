"""
File: sdf_concat.py
Concatenate numbered SDF chunks back into a single SDF.

Usage
-----
python sdf_concat.py base_0.sdf
    • The argument must be the **first** chunk
      (e.g. base_0.sdf, base-0.sdf, …).

Flow
----
1. Detect every file that shares the same prefix and ends with “_N.sdf” or “-N.sdf”.
2. Show the list and ask the user to confirm merging.
3. Merge in numeric order into <base>_concatenated.sdf.
4. Ask whether to delete the original chunk files.
"""

import argparse
import os
import re
from glob import glob
from rdkit import Chem
from tqdm import tqdm

# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Concatenate numbered SDF chunks.")
parser.add_argument("first_chunk", help="first chunk file (e.g. base_0.sdf)")
args = parser.parse_args()
first = args.first_chunk

# Identify chunk family
m = re.search(r"(.+?)([_-])(\d+)\.sdf$", first, re.IGNORECASE)
if not m:
    print("❌  File name must end with _0.sdf / -0.sdf etc.")
    exit(1)

base_prefix, sep = m.group(1), m.group(2)
pattern = f"{base_prefix}{sep}[0-9]*.sdf"
files = sorted(
    glob(pattern),
    key=lambda x: int(re.search(r"(\d+)\.sdf$", x).group(1))
)

if not files:
    print("❌  No matching chunk files found.")
    exit(1)

print("\n🔎  Files detected (numeric order):")
for f in files:
    print(" ", f)

# Confirm merge
if input("\nMerge these files? (y/n): ").strip().lower() != 'y':
    print("Aborted.")
    exit(0)

# Merge chunks
merged_mols = []
for f in tqdm(files, desc="Merging"):
    merged_mols.extend([mol for mol in Chem.SDMolSupplier(f) if mol is not None])

out_name = f"{base_prefix}_concatenated.sdf"
writer = Chem.SDWriter(out_name)
for mol in merged_mols:
    writer.write(mol)
writer.close()

print(f"\n✅  {len(files)} files merged.")
print(f"   Total molecules : {len(merged_mols)}")
print(f"   Output          : {out_name}")


# Optionally delete original chunks
if input("\nDelete original chunk files? (y/n): ").strip().lower() == 'y':
    for f in files:
        try:
            os.remove(f)
        except OSError:
            print(f"  ⚠️  Could not delete {f}")
    print("🗑️  Original chunks deleted.")
else:
    print("Original chunk files kept.")


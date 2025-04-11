"""
File: sdf_extract_samples.py
Extract the first *N* molecules from an SDF file and save them as
<original>_sample.sdf.

Example
-------
python sdf_extract_samples.py big.sdf -s 20
    → big_sample.sdf   (contains the first 20 records)
"""

import argparse
import os

# ----------------------- CLI -----------------------
parser = argparse.ArgumentParser(
    description="Extract the first N molecules from an SDF file."
)
parser.add_argument("sdf", help="input SDF file")
parser.add_argument(
    "-s", "--size", type=int, default=10,
    help="number of molecules to keep (default: 10)"
)
args = parser.parse_args()

# ------------------- read & copy -------------------
lines = []
count = 0
with open(args.sdf, "r") as f:
    for line in f:
        if line.startswith("$$$$"):
            count += 1
        lines.append(line)
        if count >= args.size:
            break

# Write the sample SDF
base, ext = os.path.splitext(args.sdf)
out_name = f"{base}_sample{ext}"
with open(out_name, "w") as f:
    f.writelines(lines)

print(f"✅  Saved {count} molecules to {out_name}")


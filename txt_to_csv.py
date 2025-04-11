"""
txt2csv_dw.py
-------------
Convert a Data Warrior TXT export (TAB‑separated) to a clean CSV.

• Removes columns starting with:
    "Structure of", "Structure [idcode]", "Unnamed:", "SmilesFragFp"
• Saves result as <input>.csv in the same directory.
"""

import argparse, os, pandas as pd

# 1) parse CLI argument -------------------------------------------------
parser = argparse.ArgumentParser(
    description="Convert a Data Warrior TXT (TAB‑separated) to CSV and drop structure columns."
)
parser.add_argument("txt", help="Data Warrior TXT file exported with TAB delimiter")
txt_file = parser.parse_args().txt

# 2) read TXT -----------------------------------------------------------
df = pd.read_csv(txt_file, sep="\t")

# 3) drop unwanted columns ---------------------------------------------
drop_cols = [
    c for c in df.columns
    if c.startswith(("Structure of", "Structure [idcode]", "Unnamed:", "SmilesFragFp"))
]
df.drop(columns=drop_cols, inplace=True)

# 4) write CSV ----------------------------------------------------------
out_csv = os.path.splitext(txt_file)[0] + ".csv"
df.to_csv(out_csv, index=False)
print(f"✅  Converted: {out_csv}")


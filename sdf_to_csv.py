"""
sdf_to_csv.py
Convert a large SDF to one CSV:

1. split_sdf_file_auto (helpers_sdf) 로 자동 분할
2. 각 조각을 병렬로 CSV 변환
3. CSV 조각 합치기
4. 임시 파일 삭제
"""

import argparse, os, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from rdkit import Chem
from helpers_sdf import split_sdf_file_auto

# ── CLI ─────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("sdf", help="input SDF")
args = parser.parse_args()
sdf_path = args.sdf

# ── 1. split ───────────────────────────────
chunks = split_sdf_file_auto(sdf_path)
print(f"split → {len(chunks)} chunks")

# ── 2. chunk→CSV (parallel) ────────────────
def to_csv(chunk):
    mols = [m for m in Chem.SDMolSupplier(chunk) if m]
    rows = [{p: m.GetProp(p) for p in m.GetPropNames()} | {"SMILES": Chem.MolToSmiles(m)}
            for m in mols]
    df = pd.DataFrame(rows)
    out = chunk.replace(".sdf", ".csv")
    df.to_csv(out, index=False)
    return out

with ProcessPoolExecutor() as pool:
    csv_chunks = list(tqdm(pool.map(to_csv, chunks), total=len(chunks), desc="chunks→csv"))

# ── 3. concat CSV ──────────────────────────
merged = pd.concat((pd.read_csv(p) for p in csv_chunks), ignore_index=True)
out_csv = os.path.splitext(sdf_path)[0] + "_merged.csv"
merged.to_csv(out_csv, index=False)
print(f"merged csv → {out_csv}  (rows={len(merged)})")

# ── 4. cleanup ─────────────────────────────
for f in chunks + csv_chunks:
    try: os.remove(f)
    except OSError: pass
print("tmp files deleted")


"""
File: smiles_string_to_mol.py
Convert a SMILES string to a single‑molecule MOL file with optional H‑addition
and 3 D optimisation.

Workflow
--------
1. Ask for an output *base* name (blank → use the SMILES string).
2. Ask whether to add explicit hydrogens.                  [y / n]
3. Ask whether to generate a 3 D conformer and optimise with UFF.   [y / n]
   └─ 3 D optimisation **does not** add hydrogens automatically.
4. Save the result as `<base>.mol`.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
import textwrap
import sys

# ── 1. SMILES argument ─────────────────────────────────────────
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent(
        """
        Convert a SMILES string to a MOL file.
        All other options (file name, add hydrogens, 3 D optimisation)
        are asked interactively.
        """
    ),
)
parser.add_argument("smiles", help="SMILES string to convert")
args = parser.parse_args()

smiles = args.smiles.strip()
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    sys.exit("❌  Invalid SMILES string – aborting.")

# ── 2. Interactive options ────────────────────────────────────
base = input("\nOutput file name (blank → use the SMILES string): ").strip()
if not base:
    base = smiles.replace("/", "_").replace("\\", "_")

if input("Add explicit hydrogens? [y / n] : ").strip().lower() == "y":
    mol = Chem.AddHs(mol)
    print("✅  Hydrogens added.")

if input("Generate 3D conformer and optimise with UFF? [y / n] : ").strip().lower() == "y":
    print("⏳  Generating 3D conformer …")
    if AllChem.EmbedMolecule(mol, randomSeed=0xF00D) == -1:
        sys.exit("❌  3D conformer generation failed – aborting.")
    print("⏳  UFF optimisation …")
    AllChem.UFFOptimizeMolecule(mol)
    print("✅  3D structure optimised.")
else:
    print("ℹ️  3D optimisation skipped (2D structure).")

# ── 3. Write MOL file ─────────────────────────────────────────
mol_path = f"{base}.mol"
Chem.MolToMolFile(mol, mol_path)
print(f"\n🎉  MOL file saved as: {mol_path}")


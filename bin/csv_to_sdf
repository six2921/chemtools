import pandas as pd
import argparse
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools

parser = argparse.ArgumentParser(description = '')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument("--name", "-n", required=True, help='name column') # 필요한 인수를 추가
parser.add_argument("--smiles", "-s", required=True, help='smiles column') # 필요한 인수를 추가
parser.add_argument('--prep', action='store_true')
args = parser.parse_args()

csv = args.csv
name = args.name
smiles = args.smiles
prep = args.prep

df = pd.read_csv(csv)
df['ROMol'] = df[smiles].apply(lambda x: Chem.MolFromSmiles(x))

fn = os.path.splitext(csv)[0]
rdkit.Chem.PandasTools.WriteSDF(df, f'{fn}.sdf', molColName='ROMol', idName=name, properties=list(df.columns), allNumeric=False)

print(f"{fn}.sdf is generated")

if prep == True:
    print("ligprep is running")
    os.system(f"/opt/schrodinger/suites2022-2/ligprep -bff 16 -pht 0.0 -epik -HOST localhost:8 -NSTRUCT 500 -isd {fn}.sdf -osd {fn}_prep.sdf")

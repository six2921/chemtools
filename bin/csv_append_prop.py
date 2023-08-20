#! /home/siu/anaconda3/envs/rdkit/bin/python

import pandas as pd
import glob
import os
import argparse

parser = argparse.ArgumentParser(description = 'Organize your csv file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv

fn = os.path.splitext(csv)[0]
files = glob.glob(f'{fn}*')

df = pd.concat([pd.read_csv(x) for x in files], axis=1)
df = df.loc[:, ~df.T.duplicated()]

df.to_csv(f"{fn}_appended.csv", index=False)

print("_appended.csv is generated.")

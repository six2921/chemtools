#! /home/siu/anaconda3/envs/rdkit/bin/python

import pandas as pd
import math
import os
import argparse

parser = argparse.ArgumentParser(description = 'Organize your csv file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
parser.add_argument('bunch', help='how many row per a file') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv
bunch = args.bunch
bunch = int(bunch)

df = pd.read_csv(csv)
fn = os.path.splitext(csv)[0]

for i in range(0, math.ceil(len(df)/bunch)):
    tf = df[i*bunch:((i+1)*bunch)]
    print(i*bunch, ((i+1)*bunch))
    tf.to_csv(f'{fn}_{i}.csv', index = False)

print("File_N.csv are generated")


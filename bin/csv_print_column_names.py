#! /home/siu/anaconda/envs/chem/bin/python

import pandas as pd
from pandas.api.types import is_numeric_dtype
import os
import argparse

parser = argparse.ArgumentParser(description = 'Organize your csv file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv

def stat(df):
    stat = pd.DataFrame(columns = ["Column", "Type", "Min", "Max", "Unique", "Null", "Dupl", "Sample"])
    num = 0
    for col in df.columns:
        min = df[col].min() if df[col].dtype in ['int64', 'float64'] else None
        max = df[col].max() if df[col].dtype in ['int64', 'float64'] else None
        dupl = None if len(df[col].unique()) == 1 else df[col].duplicated().sum()
        uniq = len(df[col].unique())
        stat.loc[num]=[col, df[col].dtype, min, max, uniq, df[col].isna().sum(), dupl, df[col][df.index.min()]]
        num +=1
    print(stat)
    return stat


df = pd.read_csv(csv)

stat(df)
print('')
print(f"# Number of molecules and features: {df.shape}".center(100))


#! /home/siu/anaconda/envs/chem/bin/python

import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description = 'Merge Two csvs')
parser.add_argument('file1', help='left csv file')
parser.add_argument('file2', help='right csv file')

args = parser.parse_args()

file1 = args.file1
file2 = args.file2

df1 = pd.read_csv(file1)
for col in df1.columns:
    try:
        print(col, '|', str(df1[col][0])[:20])
    except:
        pass

key1 = 'dZ1%mO'
while key1 not in list(df1.columns):
    key1 = input('첫번째 파일에 키 컬럼: ')
print(' ')
print(' ')


df2 = pd.read_csv(file2)
for col in df2.columns:
    try:
        print(col, '|', str(df2[col][0])[:20])
    except:
        pass

key2 = 'dZ1%mO'
while key2 not in list(df2.columns):
    key2 = input('두번째 파일에 키 컬럼: ')
print(' ')
print(' ')

how = input('병합 방식 [outer, inner, left, right]: ')

df = pd.merge(df1, df2, left_on=key1, right_on=key2, how=how, suffixes=('', '_1'))
print('Rows in Merged: ', len(df))

df.to_csv('merged_out.csv', index=False)

#! /home/siu/anaconda/bin/python

import pandas as pd 
import itertools
import argparse
import os

parser = argparse.ArgumentParser(description = 'ex) Assay concentration column have two types of values, 10uM, 3uM. You can specify this column.')
parser.add_argument('file', help='csv file')

args = parser.parse_args()

file = args.file
fn = os.path.splitext(file)[0]

df = pd.read_csv(file)
print(*list(df.columns), sep='\n')

name = 'dZ1%mO'
while name not in list(df.columns):
    name = input('물질 이름 컬럼(name): ')
print(' ')

value = 'dZ1%mO'
while value not in list(df.columns):
    value = input('실험값 컬럼(inhibition): ')
print(' ')

sort = 'dZ1%mO'
while sort not in list(df.columns):
    sort = input('실험값 종류 컬럼(species): ')
print(' ')

# make df for merge reference 
x = list(df.columns)
y = [sort, value]

merged = df[sorted(set(x) - set(y), key=x.index)].drop_duplicates(subset=name)
print("Names in this file:", len(merged))

# print unique values in sort column
sort_values = df[sort].unique()
print("Unique values in type:", sort_values)

# separates columns and merge
dups = []
for sv in sort_values:
    tf = df.loc[df[sort] == sv, [name, value]]
    tf[value+'_'+str(sv)] = tf[value]
    tf = tf.drop([value], axis=1)
    dups.append(list(tf.loc[tf.duplicated(subset=name) == True, name]))
    merged = pd.merge(merged, tf, on=name, how='outer')

dups = list(itertools.chain(*dups))

merged.to_csv(f"{fn}_dv.csv" , index=False)

print("Rows in this file:", len(merged))
print("Duplicated Names in output:", dups)

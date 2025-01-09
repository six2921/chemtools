import pandas as pd 
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description = '특정 컬럼의 카테고리 값이 다른 것을 분리할 때')
parser.add_argument('file', help='csv file')

args = parser.parse_args()

file = args.file


df = pd.read_csv(file)
print(*list(df.columns), sep='\n')

col_id = 'dZ1%mO'
while col_id not in list(df.columns):
    col_id = input('물질 이름 컬럼(name): ')
print(' ')

col_value = 'dZ1%mO'
while col_value not in list(df.columns):
    col_value = input('실험값 컬럼(inhibition): ')
print(' ')

dups = df.loc[df[col_id].duplicated() == True, col_id].unique()
dups

for dup in dups:
    df.loc[df[col_id]==dup]
    values = np.array(df.loc[df[col_id]==dup, col_value].to_list())
    mean, diff, std = np.mean(values), np.max(values) - np.min(values), np.std(values)
    
    tf = df.drop(col_value, axis=1)
    after_drop = tf.loc[df[col_id] == dup].drop_duplicates(keep=False)
    if len(after_drop) != 0 :
        print(f"{dup} has different value in their duplicates")
    else:
        keep_idx = df.loc[df[col_id] == dup].index.to_list()[0]
        df.loc[keep_idx, col_value] = mean
        
        dup_idx_list = df.loc[df[col_id] == dup].index.to_list()[1:]
        for i in dup_idx_list:
            df.drop(i, inplace=True)
        print(f"{dup}, replaced")


fn = os.path.splitext(file)[0]
df.to_csv(f"{fn}_edited.csv", index=False)

print(" ")
print("_edited.csv file was generated.")
